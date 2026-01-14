//! Algebraic Multigrid (AMG) preconditioner implementation

use super::coarsening::{
    aggregation_coarsening, falgout_coarsening, hmis_coarsening, hybrid_coarsening,
    pmis_coarsening, ruge_stueben_coarsening,
};
use super::interpolation::{
    create_classical_interpolation, create_direct_interpolation, create_standard_interpolation,
};
use super::{
    AMGConfig, AMGHierarchy, AMGStatistics, CoarseningStrategy, GaussSeidelSmoother,
    InterpolationStrategy, JacobiSmoother, MultigridLevel, MultigridSmoother, SmootherType,
    SymmetricGaussSeidelSmoother,
};
use crate::error::{MathError, Result};
use crate::linear_solver::traits::Preconditioner;
use crate::sparse::{sparse_sparse_mul, SparseMatrix};
use nalgebra::DVector;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use std::time::Instant;

/// Algebraic Multigrid preconditioner
#[derive(Clone)]
pub struct AlgebraicMultigrid<T: RealField + Copy> {
    /// Multigrid hierarchy levels
    levels: Vec<MultigridLevel<T>>,
    /// AMG configuration
    config: AMGConfig,
    /// Performance statistics
    statistics: AMGStatistics,
    /// Setup completion flag
    is_setup: bool,
    /// Optional cached hierarchy operators
    hierarchy: Option<AMGHierarchy<T>>,
}

impl<T: RealField + Copy + FromPrimitive> AlgebraicMultigrid<T> {
    /// Create a new AMG preconditioner
    pub fn new(matrix: &SparseMatrix<T>, config: AMGConfig) -> Result<Self> {
        let mut amg = Self {
            levels: Vec::new(),
            config,
            statistics: AMGStatistics::default(),
            is_setup: false,
            hierarchy: None,
        };

        amg.setup(matrix)?;
        Ok(amg)
    }

    /// Create a new AMG preconditioner with a specified configuration
    pub fn with_config(matrix: &SparseMatrix<T>, config: AMGConfig) -> Result<Self> {
        Self::new(matrix, config)
    }

    /// Create a new AMG preconditioner with a pre-existing hierarchy
    pub fn with_hierarchy(
        matrix: &SparseMatrix<T>,
        config: AMGConfig,
        hierarchy: AMGHierarchy<T>,
    ) -> Result<Self> {
        let mut amg = Self {
            levels: Vec::new(),
            config,
            statistics: AMGStatistics::default(),
            is_setup: false,
            hierarchy: Some(hierarchy),
        };

        amg.setup(matrix)?;
        Ok(amg)
    }

    /// Get the current hierarchy for caching
    pub fn get_hierarchy(&self) -> AMGHierarchy<T> {
        AMGHierarchy::from_levels(&self.levels)
    }

    /// Recompute the hierarchy operators for a new matrix with same sparsity pattern
    pub fn recompute(&mut self, fine_matrix: &SparseMatrix<T>) -> Result<()> {
        if self.levels.is_empty() {
            return Err(
                MathError::InvalidInput("AMG hierarchy not initialized".to_string()).into(),
            );
        }

        // Update the finest level matrix
        self.levels[0].matrix = fine_matrix.clone();
        self.levels[0].smoother = self.create_smoother(fine_matrix);

        // Recompute coarse matrices using existing transfer operators
        for i in 0..self.levels.len() - 1 {
            let (restriction, interpolation) = {
                let current = &self.levels[i];
                match (&current.restriction, &current.interpolation) {
                    (Some(r), Some(p)) => (r, p),
                    _ => {
                        return Err(MathError::InvalidInput(
                            "Missing transfer operators".to_string(),
                        )
                        .into())
                    }
                }
            };

            // A_coarse = R * A_fine * P
            let temp = sparse_sparse_mul(restriction, &self.levels[i].matrix);
            let coarse_matrix = sparse_sparse_mul(&temp, interpolation);

            // Update next level
            self.levels[i + 1].matrix = coarse_matrix.clone();
            self.levels[i + 1].smoother = self.create_smoother(&coarse_matrix);
        }

        Ok(())
    }

    /// Setup the AMG hierarchy
    fn setup(&mut self, fine_matrix: &SparseMatrix<T>) -> Result<()> {
        let setup_start = Instant::now();

        // Create finest level
        let finest_level = MultigridLevel {
            matrix: fine_matrix.clone(),
            restriction: None,
            interpolation: None,
            smoother: self.create_smoother(fine_matrix),
        };

        self.levels.push(finest_level);
        self.statistics.level_sizes.push(fine_matrix.nrows());

        // Build hierarchy
        if let Some(ref hierarchy) = self.hierarchy {
            // Use cached operators to build coarse matrices
            for (restriction, interpolation) in &hierarchy.operators {
                let coarse_matrix = {
                    let current_level = self.levels.last_mut().unwrap();
                    let coarse_matrix = restriction * &current_level.matrix * interpolation;

                    // Store operators in the finer level
                    current_level.restriction = Some(restriction.clone());
                    current_level.interpolation = Some(interpolation.clone());

                    coarse_matrix
                };

                let smoother = self.create_smoother(&coarse_matrix);
                self.levels.push(MultigridLevel {
                    matrix: coarse_matrix,
                    restriction: None,
                    interpolation: None,
                    smoother,
                });
                self.statistics
                    .level_sizes
                    .push(self.levels.last().unwrap().matrix.nrows());
            }
        } else {
            // Standard build
            while self.levels.len() < self.config.max_levels
                && self.levels.last().unwrap().matrix.nrows() > self.config.min_coarse_size
            {
                self.build_next_level()?;
            }
        }

        // Record statistics
        self.statistics.num_levels = self.levels.len();
        self.statistics.setup_time = setup_start.elapsed().as_secs_f64();
        self.compute_complexities();

        self.is_setup = true;
        Ok(())
    }

    /// Build the next coarser level in the hierarchy
    fn build_next_level(&mut self) -> Result<()> {
        let (restriction, interpolation, coarse_matrix) = {
            let current_level = self.levels.last().unwrap();

            // Perform coarsening using the configured strategy
            let coarsening_result = match self.config.coarsening_strategy {
                CoarseningStrategy::RugeStueben => ruge_stueben_coarsening(
                    &current_level.matrix,
                    T::from_f64(self.config.strength_threshold).unwrap_or_else(T::zero),
                ),
                CoarseningStrategy::Aggregation => {
                    aggregation_coarsening(&current_level.matrix, 8) // Use default aggregate size
                }
                CoarseningStrategy::Hybrid => hybrid_coarsening(
                    &current_level.matrix,
                    T::from_f64(self.config.strength_threshold).unwrap_or_else(T::zero),
                    4,
                ),
                CoarseningStrategy::Falgout => falgout_coarsening(
                    &current_level.matrix,
                    T::from_f64(self.config.strength_threshold).unwrap_or_else(T::zero),
                ),
                CoarseningStrategy::PMIS => pmis_coarsening(
                    &current_level.matrix,
                    T::from_f64(self.config.strength_threshold).unwrap_or_else(T::zero),
                ),
                CoarseningStrategy::HMIS => hmis_coarsening(
                    &current_level.matrix,
                    T::from_f64(self.config.strength_threshold).unwrap_or_else(T::zero),
                    T::from_f64(0.5).unwrap_or_else(T::zero),
                ),
            }
            .map_err(|e| MathError::InvalidInput(format!("Coarsening failed: {e}")))?;

            // Create interpolation operator based on the coarsening result
            let interpolation = match self.config.interpolation_strategy {
                InterpolationStrategy::Classical => create_classical_interpolation(
                    &current_level.matrix,
                    &coarsening_result.coarse_points,
                    &coarsening_result.strength_matrix,
                    self.config.max_interpolation_points,
                ),
                InterpolationStrategy::Direct => Ok(create_direct_interpolation(
                    &coarsening_result.fine_to_coarse_map,
                    current_level.matrix.nrows(),
                    coarsening_result.coarse_points.len(),
                )),
                InterpolationStrategy::Standard => create_standard_interpolation(
                    &current_level.matrix,
                    &coarsening_result.coarse_points,
                    &coarsening_result.strength_matrix,
                ),
            }
            .map_err(|e| MathError::InvalidInput(format!("Interpolation failed: {e}")))?;

            // Create restriction operator (transpose of interpolation)
            // In nalgebra-sparse, transpose of CSR is CSC, then convert back to CSR
            let restriction = SparseMatrix::from(&interpolation.clone().transpose_as_csc());

            // Create coarse matrix: R * A_fine * P
            let temp_matrix = sparse_sparse_mul(&restriction, &current_level.matrix);
            let coarse_matrix = sparse_sparse_mul(&temp_matrix, &interpolation);

            (restriction, interpolation, coarse_matrix)
        };

        // Update the current level with its operators to the next level
        let last_idx = self.levels.len() - 1;
        self.levels[last_idx].restriction = Some(restriction);
        self.levels[last_idx].interpolation = Some(interpolation);

        // Create smoother for this level
        let smoother = self.create_smoother(&coarse_matrix);

        // Add new level (coarsest level starts with no restriction/interpolation)
        let level_size = coarse_matrix.nrows();
        let new_level = MultigridLevel {
            matrix: coarse_matrix,
            restriction: None,
            interpolation: None,
            smoother,
        };

        self.levels.push(new_level);
        self.statistics.level_sizes.push(level_size);

        Ok(())
    }

    /// Create smoother for a given matrix
    fn create_smoother(&self, _matrix: &SparseMatrix<T>) -> Box<dyn MultigridSmoother<T>> {
        let relaxation = T::from_f64(self.config.relaxation_factor).unwrap_or_else(T::one);
        match self.config.smoother_type {
            SmootherType::SymmetricGaussSeidel => {
                Box::new(SymmetricGaussSeidelSmoother::new(relaxation))
            }
            SmootherType::Jacobi => Box::new(JacobiSmoother::new(relaxation)),
            _ => Box::new(GaussSeidelSmoother::new(relaxation)),
        }
    }

    /// Compute operator and grid complexities
    fn compute_complexities(&mut self) {
        if self.levels.is_empty() {
            return;
        }

        let fine_nnz = self.levels[0]
            .matrix
            .values()
            .iter()
            .filter(|&&x| x.abs() > T::from_f64(1e-15).unwrap_or_else(T::zero))
            .count();
        let mut total_nnz = 0;
        let mut total_vars = 0;

        for level in &self.levels {
            total_nnz += level
                .matrix
                .values()
                .iter()
                .filter(|&&x| x.abs() > T::from_f64(1e-15).unwrap_or_else(T::zero))
                .count();
            total_vars += level.matrix.nrows();
        }

        self.statistics.operator_complexity = total_nnz as f64 / fine_nnz as f64;
        self.statistics.grid_complexity = total_vars as f64 / self.levels[0].matrix.nrows() as f64;
    }

    /// Perform a single V-cycle
    fn v_cycle(&self, level_idx: usize, b: &DVector<T>, x: &mut DVector<T>) {
        let current_level = &self.levels[level_idx];

        // 1. Pre-smoothing
        current_level.smoother.apply(
            &current_level.matrix,
            x,
            b,
            self.config.pre_smooth_iterations,
        );

        if level_idx < self.levels.len() - 1 {
            // 2. Compute residual: r = b - A*x
            let r = b - &current_level.matrix * &*x;

            // 3. Restriction: r_coarse = R * r
            let r_coarse = current_level.restriction.as_ref().unwrap() * r;

            // 4. Recursive call to coarse level
            let mut e_coarse = DVector::zeros(r_coarse.len());
            self.v_cycle(level_idx + 1, &r_coarse, &mut e_coarse);

            // 5. Interpolation: e = P * e_coarse
            let e = current_level.interpolation.as_ref().unwrap() * e_coarse;

            // 6. Correction: x = x + e
            *x += e;

            // 7. Post-smoothing
            current_level.smoother.apply(
                &current_level.matrix,
                x,
                b,
                self.config.post_smooth_iterations,
            );
        } else {
            // Coarsest level solve (can be more iterations or a direct solve)
            current_level
                .smoother
                .apply(&current_level.matrix, x, b, 10);
        }
    }
}

impl<T: RealField + Copy + FromPrimitive> Preconditioner<T> for AlgebraicMultigrid<T> {
    fn apply_to(&self, r: &DVector<T>, z: &mut DVector<T>) -> Result<()> {
        z.fill(T::zero());

        self.v_cycle(0, r, z);

        Ok(())
    }
}
