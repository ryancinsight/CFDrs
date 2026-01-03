//! Algebraic Multigrid (AMG) preconditioner implementation

use super::{
    AMGConfig, AMGStatistics, AMGHierarchy, CoarseningStrategy, GaussSeidelSmoother, InterpolationStrategy,
    JacobiSmoother, MultigridLevel, MultigridSmoother, SmootherType, SymmetricGaussSeidelSmoother,
};
use super::coarsening::{
    ruge_stueben_coarsening, aggregation_coarsening, hybrid_coarsening,
    falgout_coarsening, pmis_coarsening, hmis_coarsening,
};
use super::interpolation::{
    create_classical_interpolation, create_direct_interpolation, create_standard_interpolation,
};
use crate::error::{Result, MathError};
use crate::linear_solver::traits::Preconditioner;
use crate::sparse::{SparseMatrix, sparse_sparse_mul};
use nalgebra_sparse::CsrMatrix;
use nalgebra::{DVector, DMatrix};
use std::time::Instant;

/// Algebraic Multigrid preconditioner
pub struct AlgebraicMultigrid {
    /// Multigrid hierarchy levels
    levels: Vec<MultigridLevel<f64>>,
    /// AMG configuration
    config: AMGConfig,
    /// Performance statistics
    statistics: AMGStatistics,
    /// Setup completion flag
    is_setup: bool,
    /// Optional cached hierarchy operators
    hierarchy: Option<AMGHierarchy<f64>>,
}

impl AlgebraicMultigrid {
    /// Create a new AMG preconditioner
    pub fn new(matrix: &SparseMatrix<f64>, config: AMGConfig) -> Result<Self> {
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

    /// Create a new AMG preconditioner with a pre-existing hierarchy
    pub fn with_hierarchy(
        matrix: &SparseMatrix<f64>,
        config: AMGConfig,
        hierarchy: AMGHierarchy<f64>,
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
    pub fn get_hierarchy(&self) -> AMGHierarchy<f64> {
        AMGHierarchy::from_levels(&self.levels)
    }

    /// Setup the AMG hierarchy
    fn setup(&mut self, fine_matrix: &SparseMatrix<f64>) -> Result<()> {
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
                self.statistics.level_sizes.push(self.levels.last().unwrap().matrix.nrows());
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
                CoarseningStrategy::RugeStueben => {
                    ruge_stueben_coarsening(&current_level.matrix, self.config.strength_threshold)
                }
                CoarseningStrategy::Aggregation => {
                    aggregation_coarsening(&current_level.matrix, 8) // Use default aggregate size
                }
                CoarseningStrategy::Hybrid => {
                    hybrid_coarsening(&current_level.matrix, self.config.strength_threshold, 4)
                }
                CoarseningStrategy::Falgout => {
                    falgout_coarsening(&current_level.matrix, self.config.strength_threshold)
                }
                CoarseningStrategy::PMIS => {
                    pmis_coarsening(&current_level.matrix, self.config.strength_threshold)
                }
                CoarseningStrategy::HMIS => {
                    hmis_coarsening(&current_level.matrix, self.config.strength_threshold, 0.5)
                }
            }.map_err(|e| MathError::InvalidInput(format!("Coarsening failed: {e}")))?;

            // Create interpolation operator based on the coarsening result
            let interpolation = match self.config.interpolation_strategy {
                InterpolationStrategy::Classical => {
                    create_classical_interpolation(
                        &current_level.matrix,
                        &coarsening_result.coarse_points,
                        &coarsening_result.strength_matrix,
                        self.config.max_interpolation_points
                    )
                }
                InterpolationStrategy::Direct => {
                    Ok(create_direct_interpolation(
                        &coarsening_result.fine_to_coarse_map,
                        current_level.matrix.nrows(),
                        coarsening_result.coarse_points.len()
                    ))
                }
                InterpolationStrategy::Standard => {
                    create_standard_interpolation(
                        &current_level.matrix,
                        &coarsening_result.coarse_points,
                        &coarsening_result.strength_matrix
                    )
                }
            }.map_err(|e| MathError::InvalidInput(format!("Interpolation failed: {e}")))?;

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
    fn create_smoother(&self, _matrix: &SparseMatrix<f64>) -> Box<dyn MultigridSmoother<f64>> {
        match self.config.smoother_type {
            SmootherType::GaussSeidel => {
                Box::new(GaussSeidelSmoother::new(self.config.relaxation_factor))
            }
            SmootherType::SymmetricGaussSeidel => Box::new(SymmetricGaussSeidelSmoother::new(
                self.config.relaxation_factor,
            )),
            SmootherType::Jacobi => Box::new(JacobiSmoother::new(self.config.relaxation_factor)),
            _ => Box::new(GaussSeidelSmoother::new(self.config.relaxation_factor)),
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
            .filter(|&&x| x.abs() > 1e-15)
            .count();
        let mut total_nnz = 0;
        let mut total_vars = 0;

        for level in &self.levels {
            total_nnz += level.matrix.values().iter().filter(|&&x| x.abs() > 1e-15).count();
            total_vars += level.matrix.nrows();
        }

        self.statistics.operator_complexity = total_nnz as f64 / fine_nnz as f64;
        self.statistics.grid_complexity = total_vars as f64 / self.levels[0].matrix.nrows() as f64;
    }

    /// Apply multigrid cycle using configured cycle type
    fn apply_cycle(&self, residual: &DVector<f64>) -> Result<DVector<f64>> {
        let result = match self.config.cycle_type {
            super::CycleType::VCycle => {
                let (correction, _) = super::cycles::apply_v_cycle(
                    &self.levels,
                    residual,
                    1,     // Single cycle
                    1e-12, // Internal tolerance
                )?;
                correction
            }
            super::CycleType::WCycle => {
                let (correction, _) = super::cycles::apply_w_cycle(
                    &self.levels,
                    residual,
                    1, // Single cycle
                    1e-12,
                )?;
                correction
            }
            super::CycleType::FCycle => {
                let (correction, _) = super::cycles::apply_f_cycle(
                    &self.levels,
                    residual,
                    1, // Single cycle
                    1e-12,
                )?;
                correction
            }
        };

        Ok(result)
    }

    /// Get AMG statistics
    pub fn statistics(&self) -> &AMGStatistics {
        &self.statistics
    }
}

impl Preconditioner<f64> for AlgebraicMultigrid {
    fn apply_to(&self, r: &DVector<f64>, z: &mut DVector<f64>) -> Result<()> {
        if !self.is_setup {
            return Err(MathError::InvalidInput(
                "AMG preconditioner not properly set up".to_string(),
            )
            .into());
        }

        // Apply multigrid cycle to residual
        let result = self.apply_cycle(r)?;
        z.copy_from(&result);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    #[test]
    fn test_amg_creation() {
        // Create a simple 2D Poisson matrix
        let n = 10;
        let mut matrix = DMatrix::zeros(n * n, n * n);

        for i in 0..n {
            for j in 0..n {
                let idx = i * n + j;
                matrix[(idx, idx)] = 4.0; // Main diagonal

                // Neighbors
                if i > 0 {
                    matrix[(idx, (i - 1) * n + j)] = -1.0;
                }
                if i < n - 1 {
                    matrix[(idx, (i + 1) * n + j)] = -1.0;
                }
                if j > 0 {
                    matrix[(idx, i * n + (j - 1))] = -1.0;
                }
                if j < n - 1 {
                    matrix[(idx, i * n + (j + 1))] = -1.0;
                }
            }
        }

        let config = AMGConfig {
            max_levels: 3,
            min_coarse_size: 10,
            ..Default::default()
        };

        let sparse_matrix = CsrMatrix::from(&matrix);
        let amg = AlgebraicMultigrid::new(&sparse_matrix, config);
        assert!(amg.is_ok());
    }

    #[test]
    fn test_amg_statistics() {
        let matrix = DMatrix::identity(100, 100);
        let sparse_matrix = CsrMatrix::from(&matrix);
        let config = AMGConfig::default();

        let amg = AlgebraicMultigrid::new(&sparse_matrix, config).unwrap();
        let stats = amg.statistics();

        assert!(stats.num_levels >= 1);
        assert!(stats.operator_complexity >= 1.0);
        assert!(stats.grid_complexity >= 1.0);
    }

    #[test]
    fn test_amg_hierarchy_caching() {
        // Create a 2D Poisson matrix
        let n = 20;
        let mut matrix = DMatrix::zeros(n * n, n * n);
        for i in 0..n {
            for j in 0..n {
                let idx = i * n + j;
                matrix[(idx, idx)] = 4.0;
                if i > 0 { matrix[(idx, (i - 1) * n + j)] = -1.0; }
                if i < n - 1 { matrix[(idx, (i + 1) * n + j)] = -1.0; }
                if j > 0 { matrix[(idx, i * n + (j - 1))] = -1.0; }
                if j < n - 1 { matrix[(idx, i * n + (j + 1))] = -1.0; }
            }
        }

        let config = AMGConfig {
            max_levels: 4,
            min_coarse_size: 10,
            ..Default::default()
        };

        // First solve - standard setup
        let setup_start = Instant::now();
        let sparse_matrix = CsrMatrix::from(&matrix);
        let amg1 = AlgebraicMultigrid::new(&sparse_matrix, config.clone()).unwrap();
        let time1 = setup_start.elapsed();
        let hierarchy = amg1.get_hierarchy();

        // Second solve - cached setup
        let setup_start2 = Instant::now();
        let amg2 = AlgebraicMultigrid::with_hierarchy(&sparse_matrix, config, hierarchy).unwrap();
        let time2 = setup_start2.elapsed();

        // Verify hierarchy sizes match
        assert_eq!(amg1.statistics().num_levels, amg2.statistics().num_levels);
        assert_eq!(amg1.statistics().level_sizes, amg2.statistics().level_sizes);

        // Cached setup should be faster (usually significantly for larger matrices)
        println!("Standard setup: {:?}, Cached setup: {:?}", time1, time2);
        // We don't assert time1 > time2 because for small matrices noise might interfere,
        // but it's a good sanity check during manual runs.
    }
}
