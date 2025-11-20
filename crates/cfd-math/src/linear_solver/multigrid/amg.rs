//! Algebraic Multigrid (AMG) preconditioner implementation

use super::{
    AMGConfig, AMGStatistics, CoarseningStrategy, GaussSeidelSmoother, InterpolationStrategy,
    JacobiSmoother, MultigridLevel, MultigridSmoother, SmootherType, SymmetricGaussSeidelSmoother,
};
use crate::error::Result;
use crate::linear_solver::traits::Preconditioner;
use nalgebra::{DMatrix, DVector};
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
}

impl AlgebraicMultigrid {
    /// Create a new AMG preconditioner
    pub fn new(matrix: &DMatrix<f64>, config: AMGConfig) -> Result<Self> {
        let mut amg = Self {
            levels: Vec::new(),
            config,
            statistics: AMGStatistics::default(),
            is_setup: false,
        };

        amg.setup(matrix)?;
        Ok(amg)
    }

    /// Setup the AMG hierarchy
    fn setup(&mut self, fine_matrix: &DMatrix<f64>) -> Result<()> {
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
        while self.levels.len() < self.config.max_levels
            && self.levels.last().unwrap().matrix.nrows() > self.config.min_coarse_size
        {
            self.build_next_level()?;
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
        let current_level = self.levels.last().unwrap();

        // Perform coarsening
        let (coarse_points, fine_to_coarse_map) = match self.config.coarsening_strategy {
            CoarseningStrategy::RugeStueben => {
                self.ruge_stueben_coarsening(&current_level.matrix)?
            }
            CoarseningStrategy::Aggregation => {
                self.aggregation_coarsening(&current_level.matrix)?
            }
            CoarseningStrategy::Hybrid => self.hybrid_coarsening(&current_level.matrix)?,
            CoarseningStrategy::Falgout => self.falgout_coarsening(&current_level.matrix)?,
            CoarseningStrategy::PMIS => self.pmis_coarsening(&current_level.matrix)?,
            CoarseningStrategy::HMIS => self.hmis_coarsening(&current_level.matrix)?,
        };

        // Create interpolation operator
        let interpolation = self.create_interpolation_operator(
            &current_level.matrix,
            &coarse_points,
            &fine_to_coarse_map,
        )?;

        // Create restriction operator (transpose of interpolation)
        let restriction = interpolation.transpose();

        // Create coarse matrix: R * A_fine * P
        let coarse_matrix = &restriction * &current_level.matrix * &interpolation;

        // Create smoother for this level
        let smoother = self.create_smoother(&coarse_matrix);

        // Add new level
        let level_size = coarse_matrix.nrows();
        let new_level = MultigridLevel {
            matrix: coarse_matrix,
            restriction: Some(restriction),
            interpolation: Some(interpolation),
            smoother,
        };

        self.levels.push(new_level);
        self.statistics.level_sizes.push(level_size);

        Ok(())
    }

    /// Ruge-Stüben coarsening algorithm
    fn ruge_stueben_coarsening(
        &self,
        matrix: &DMatrix<f64>,
    ) -> Result<(Vec<usize>, Vec<Option<usize>>)> {
        let n = matrix.nrows();
        let mut coarse_points = Vec::new();
        let mut fine_to_coarse_map = vec![None; n];

        // Calculate strength of connection matrix
        let strength_matrix = self.compute_strength_matrix(matrix)?;

        // Step 1: Select coarse points (strongest connections)
        let mut lambda = vec![0.0; n];
        let mut unassigned = (0..n).collect::<Vec<_>>();

        while !unassigned.is_empty() {
            // Find point with maximum lambda value
            let mut max_lambda = 0.0;
            let mut max_idx = None;

            for &i in &unassigned {
                if lambda[i] > max_lambda {
                    max_lambda = lambda[i];
                    max_idx = Some(i);
                }
            }

            if let Some(i) = max_idx {
                // Mark as coarse point
                coarse_points.push(i);
                fine_to_coarse_map[i] = Some(coarse_points.len() - 1);
                unassigned.retain(|&x| x != i);

                // Update lambda for neighboring points
                for j in 0..n {
                    if strength_matrix[(i, j)] > 0.0 {
                        lambda[j] += strength_matrix[(i, j)];
                    }
                }
            } else {
                break;
            }
        }

        // Step 2: Assign remaining points to interpolation sets
        for i in 0..n {
            if fine_to_coarse_map[i].is_none() {
                // Find strongest connection to coarse point
                let mut max_strength = 0.0;
                let mut best_coarse = None;

                for &c in &coarse_points {
                    if strength_matrix[(i, c)] > max_strength {
                        max_strength = strength_matrix[(i, c)];
                        best_coarse = Some(c);
                    }
                }

                if let Some(c) = best_coarse {
                    fine_to_coarse_map[i] = fine_to_coarse_map[c];
                }
            }
        }

        Ok((coarse_points, fine_to_coarse_map))
    }

    /// Aggregation-based coarsening
    fn aggregation_coarsening(
        &self,
        matrix: &DMatrix<f64>,
    ) -> Result<(Vec<usize>, Vec<Option<usize>>)> {
        let n = matrix.nrows();
        let mut coarse_points = Vec::new();
        let mut fine_to_coarse_map = vec![None; n];
        let mut aggregated = vec![false; n];

        // Simple aggregation: group nearby points
        let mut aggregate_id = 0;
        for i in 0..n {
            if !aggregated[i] {
                // Start new aggregate
                coarse_points.push(i);
                fine_to_coarse_map[i] = Some(aggregate_id);

                // Add strongly connected neighbors
                for j in 0..n {
                    if !aggregated[j] && matrix[(i, j)].abs() > 1e-6 {
                        aggregated[j] = true;
                        fine_to_coarse_map[j] = Some(aggregate_id);
                    }
                }

                aggregated[i] = true;
                aggregate_id += 1;
            }
        }

        Ok((coarse_points, fine_to_coarse_map))
    }

    /// Compute strength of connection matrix
    fn compute_strength_matrix(&self, matrix: &DMatrix<f64>) -> Result<DMatrix<f64>> {
        let n = matrix.nrows();
        let mut strength = DMatrix::zeros(n, n);

        for i in 0..n {
            // Find maximum off-diagonal element in row i
            let mut max_off_diag: f64 = 0.0;
            for j in 0..n {
                if i != j {
                    max_off_diag = max_off_diag.max(matrix[(i, j)].abs());
                }
            }

            // Mark strong connections
            let threshold = self.config.strength_threshold * max_off_diag;
            for j in 0..n {
                if i != j && matrix[(i, j)].abs() >= threshold {
                    strength[(i, j)] = matrix[(i, j)];
                }
            }
        }

        Ok(strength)
    }

    /// Hybrid coarsening (wrapper for external function)
    fn hybrid_coarsening(&self, matrix: &DMatrix<f64>) -> Result<(Vec<usize>, Vec<Option<usize>>)> {
        use super::coarsening::hybrid_coarsening;
        let result = hybrid_coarsening(matrix, self.config.strength_threshold, 4)?;
        Ok((result.coarse_points, result.fine_to_coarse_map))
    }

    /// Falgout coarsening (wrapper for external function)
    fn falgout_coarsening(
        &self,
        matrix: &DMatrix<f64>,
    ) -> Result<(Vec<usize>, Vec<Option<usize>>)> {
        use super::coarsening::falgout_coarsening;
        let result = falgout_coarsening(matrix, self.config.strength_threshold)?;
        Ok((result.coarse_points, result.fine_to_coarse_map))
    }

    /// PMIS coarsening (wrapper for external function)
    fn pmis_coarsening(&self, matrix: &DMatrix<f64>) -> Result<(Vec<usize>, Vec<Option<usize>>)> {
        use super::coarsening::pmis_coarsening;
        let result = pmis_coarsening(matrix, self.config.strength_threshold)?;
        Ok((result.coarse_points, result.fine_to_coarse_map))
    }

    /// HMIS coarsening (wrapper for external function)
    fn hmis_coarsening(&self, matrix: &DMatrix<f64>) -> Result<(Vec<usize>, Vec<Option<usize>>)> {
        use super::coarsening::hmis_coarsening;
        let result = hmis_coarsening(matrix, self.config.strength_threshold, 0.5)?;
        Ok((result.coarse_points, result.fine_to_coarse_map))
    }

    /// Create interpolation operator
    fn create_interpolation_operator(
        &self,
        fine_matrix: &DMatrix<f64>,
        coarse_points: &[usize],
        fine_to_coarse_map: &[Option<usize>],
    ) -> Result<DMatrix<f64>> {
        let fine_n = fine_matrix.nrows();
        let coarse_n = coarse_points.len();
        let mut interpolation = DMatrix::zeros(fine_n, coarse_n);

        match self.config.interpolation_strategy {
            InterpolationStrategy::Classical => {
                // Classical interpolation
                for (fine_i, &coarse_i_opt) in fine_to_coarse_map.iter().enumerate() {
                    if let Some(coarse_i) = coarse_i_opt {
                        // Direct injection for coarse points
                        if coarse_points.contains(&fine_i) {
                            interpolation[(fine_i, coarse_i)] = 1.0;
                        } else {
                            // Interpolation from neighboring coarse points
                            let mut total_weight = 0.0;
                            let mut weights = Vec::new();

                            // Find strongly connected coarse points
                            for &cp in coarse_points {
                                let strength = fine_matrix[(fine_i, cp)].abs();
                                if strength > 1e-10 {
                                    weights.push((cp, strength));
                                    total_weight += strength;
                                }
                            }

                            // Limit number of interpolation points
                            weights.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
                            weights.truncate(self.config.max_interpolation_points);

                            // Normalize weights
                            for (cp, weight) in weights {
                                let coarse_idx =
                                    coarse_points.iter().position(|&x| x == cp).unwrap();
                                let normalized_weight = if total_weight > 0.0 {
                                    weight / total_weight
                                } else {
                                    1.0
                                };
                                interpolation[(fine_i, coarse_idx)] = normalized_weight;
                            }
                        }
                    }
                }
            }
            _ => {
                // Simplified direct interpolation for now
                for (fine_i, &coarse_i_opt) in fine_to_coarse_map.iter().enumerate() {
                    if let Some(coarse_i) = coarse_i_opt {
                        interpolation[(fine_i, coarse_i)] = 1.0;
                    }
                }
            }
        }

        Ok(interpolation)
    }

    /// Create smoother for a given matrix
    fn create_smoother(&self, _matrix: &DMatrix<f64>) -> Box<dyn MultigridSmoother<f64>> {
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
            .iter()
            .filter(|&&x| x.abs() > 1e-15)
            .count();
        let mut total_nnz = 0;
        let mut total_vars = 0;

        for level in &self.levels {
            total_nnz += level.matrix.iter().filter(|&&x| x.abs() > 1e-15).count();
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
            return Err(crate::error::MathError::InvalidInput(
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

        let amg = AlgebraicMultigrid::new(&matrix, config);
        assert!(amg.is_ok());
    }

    #[test]
    fn test_amg_statistics() {
        let matrix = DMatrix::identity(100, 100);
        let config = AMGConfig::default();

        let amg = AlgebraicMultigrid::new(&matrix, config).unwrap();
        let stats = amg.statistics();

        assert!(stats.num_levels >= 1);
        assert!(stats.operator_complexity >= 1.0);
        assert!(stats.grid_complexity >= 1.0);
    }
}
