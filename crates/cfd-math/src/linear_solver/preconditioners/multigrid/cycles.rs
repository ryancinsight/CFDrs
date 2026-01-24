//! Multigrid cycle algorithms for AMG

use super::MultigridLevel;
use crate::SparseMatrix;
use cfd_core::error::{Error, NumericalErrorKind, Result};
use nalgebra::{DMatrix, DVector};
use std::time::Instant;

/// Apply V-cycle multigrid algorithm
pub fn apply_v_cycle(
    levels: &[MultigridLevel<f64>],
    residual: &DVector<f64>,
    max_cycles: usize,
    tolerance: f64,
) -> Result<(DVector<f64>, CycleStatistics)> {
    let start_time = Instant::now();

    if levels.is_empty() {
        return Err(cfd_core::error::Error::InvalidConfiguration(
            "No multigrid levels available".to_string(),
        ));
    }

    let mut correction = DVector::zeros(residual.len());
    let mut cycle_count = 0;
    let mut residual_history: Vec<f64> = Vec::new();

    // Continue cycling until convergence or max cycles reached
    for cycle in 0..max_cycles {
        cycle_count = cycle + 1;

        // Apply one V-cycle
        apply_multigrid_cycle(1, levels, residual, &mut correction)?;

        let residual_norm = (residual - &levels[0].matrix * &correction).norm();
        residual_history.push(residual_norm);
        if residual_norm < tolerance {
            break;
        }
    }

    let cycle_time = start_time.elapsed().as_secs_f64();

    let convergence_factor = if residual_history.len() >= 2 {
        let mut ratio_product = 1.0;
        for i in 1..residual_history.len() {
            let previous = residual_history[i - 1];
            if previous == 0.0 {
                ratio_product = 0.0;
                break;
            }
            ratio_product *= residual_history[i] / previous;
        }
        ratio_product.powf(1.0 / (residual_history.len() as f64 - 1.0))
    } else {
        0.0
    };

    let stats = CycleStatistics {
        cycle_type: CycleType::VCycle,
        total_cycles: cycle_count,
        convergence_factor,
        total_time: cycle_time,
        time_per_cycle: cycle_time / cycle_count as f64,
    };

    Ok((correction, stats))
}

/// Apply a single multigrid cycle iteration
fn apply_multigrid_cycle(
    gamma: usize,
    levels: &[MultigridLevel<f64>],
    residual: &DVector<f64>,
    correction: &mut DVector<f64>,
) -> Result<()> {
    if levels.is_empty() {
        return Ok(());
    }

    if levels.len() == 1 {
        // Coarsest level - solve exactly
        return solve_coarsest_level(&levels[0].matrix, residual, correction);
    }

    let current_level = &levels[0];

    // 1. Pre-smoothing
    current_level.smoother.apply(
        &current_level.matrix,
        correction,
        residual,
        2, // pre_smooth_iterations
    );

    // 2. Compute residual after pre-smoothing
    // r = b - A*x
    let r_fine = residual - &current_level.matrix * &*correction;

    // 3. Restrict residual
    let r_coarse = if let Some(ref restriction) = current_level.restriction {
        restriction * r_fine
    } else {
        return Err(Error::InvalidConfiguration(
            "Missing restriction operator".to_string(),
        ));
    };

    // 4. Recursive call for coarser levels
    let mut coarse_correction = DVector::zeros(r_coarse.len());
    for _ in 0..gamma {
        apply_multigrid_cycle(gamma, &levels[1..], &r_coarse, &mut coarse_correction)?;
    }

    // 5. Interpolate correction back
    if let Some(ref interpolation) = current_level.interpolation {
        let fine_correction = interpolation * coarse_correction;
        *correction += fine_correction;
    } else {
        return Err(Error::InvalidConfiguration(
            "Missing interpolation operator".to_string(),
        ));
    }

    // 6. Post-smoothing
    current_level.smoother.apply(
        &current_level.matrix,
        correction,
        residual,
        2, // post_smooth_iterations
    );

    Ok(())
}

/// Apply W-cycle multigrid algorithm
pub fn apply_w_cycle(
    levels: &[MultigridLevel<f64>],
    residual: &DVector<f64>,
    max_cycles: usize,
    tolerance: f64,
) -> Result<(DVector<f64>, CycleStatistics)> {
    let start_time = Instant::now();

    if levels.is_empty() {
        return Err(cfd_core::error::Error::InvalidConfiguration(
            "No multigrid levels available".to_string(),
        ));
    }

    let mut correction = DVector::zeros(residual.len());
    let mut cycle_count = 0;
    let mut residual_history: Vec<f64> = Vec::new();

    for cycle in 0..max_cycles {
        cycle_count = cycle + 1;

        // Apply W-cycle (gamma = 2)
        apply_multigrid_cycle(2, levels, residual, &mut correction)?;

        // Check convergence
        let residual_norm = (residual - &levels[0].matrix * &correction).norm();
        residual_history.push(residual_norm);
        if residual_norm < tolerance {
            break;
        }
    }

    let cycle_time = start_time.elapsed().as_secs_f64();

    let convergence_factor = if residual_history.len() >= 2 {
        let mut ratio_product = 1.0;
        for i in 1..residual_history.len() {
            let previous = residual_history[i - 1];
            if previous == 0.0 {
                ratio_product = 0.0;
                break;
            }
            ratio_product *= residual_history[i] / previous;
        }
        ratio_product.powf(1.0 / (residual_history.len() as f64 - 1.0))
    } else {
        0.0
    };

    let stats = CycleStatistics {
        cycle_type: CycleType::WCycle,
        total_cycles: cycle_count,
        convergence_factor,
        total_time: cycle_time,
        time_per_cycle: cycle_time / cycle_count as f64,
    };

    Ok((correction, stats))
}

/// Apply Full Multigrid (F-cycle) algorithm
pub fn apply_f_cycle(
    levels: &[MultigridLevel<f64>],
    rhs: &DVector<f64>,
    max_cycles: usize,
    tolerance: f64,
) -> Result<(DVector<f64>, CycleStatistics)> {
    let start_time = Instant::now();

    if levels.is_empty() {
        return Err(cfd_core::error::Error::InvalidConfiguration(
            "No multigrid levels available".to_string(),
        ));
    }

    if levels.len() == 1 {
        let (correction, stats) = apply_v_cycle(levels, rhs, max_cycles, tolerance)?;
        let mut final_stats = stats;
        final_stats.cycle_type = CycleType::FCycle;
        final_stats.total_time = start_time.elapsed().as_secs_f64();
        return Ok((correction, final_stats));
    }

    let mut rhs_levels: Vec<DVector<f64>> = Vec::with_capacity(levels.len());
    rhs_levels.push(rhs.clone());

    for level_idx in 0..levels.len() - 1 {
        let restriction = levels[level_idx]
            .restriction
            .as_ref()
            .ok_or_else(|| Error::InvalidConfiguration("Missing restriction operator".to_string()))?;
        let coarse_rhs = restriction * &rhs_levels[level_idx];
        rhs_levels.push(coarse_rhs);
    }

    let coarsest_idx = levels.len() - 1;
    let mut correction = DVector::zeros(rhs_levels[coarsest_idx].len());
    solve_coarsest_level(
        &levels[coarsest_idx].matrix,
        &rhs_levels[coarsest_idx],
        &mut correction,
    )?;

    let mut cycle_count = 0;
    let mut residual_history: Vec<f64> = Vec::new();

    for level_idx in (0..coarsest_idx).rev() {
        let interpolation = levels[level_idx]
            .interpolation
            .as_ref()
            .ok_or_else(|| Error::InvalidConfiguration("Missing interpolation operator".to_string()))?;
        let mut fine_correction = interpolation * &correction;

        apply_multigrid_cycle(
            1,
            &levels[level_idx..],
            &rhs_levels[level_idx],
            &mut fine_correction,
        )?;

        if level_idx == 0 {
            cycle_count += 1;
            let mut residual_norm =
                (&rhs_levels[0] - &levels[0].matrix * &fine_correction).norm();
            residual_history.push(residual_norm);

            while cycle_count < max_cycles && residual_norm >= tolerance {
                apply_multigrid_cycle(
                    1,
                    &levels[level_idx..],
                    &rhs_levels[level_idx],
                    &mut fine_correction,
                )?;
                cycle_count += 1;
                residual_norm =
                    (&rhs_levels[0] - &levels[0].matrix * &fine_correction).norm();
                residual_history.push(residual_norm);
            }
        }

        correction = fine_correction;
    }

    let cycle_time = start_time.elapsed().as_secs_f64();

    let convergence_factor = if residual_history.len() >= 2 {
        let mut ratio_product = 1.0;
        for i in 1..residual_history.len() {
            let previous = residual_history[i - 1];
            if previous == 0.0 {
                ratio_product = 0.0;
                break;
            }
            ratio_product *= residual_history[i] / previous;
        }
        ratio_product.powf(1.0 / (residual_history.len() as f64 - 1.0))
    } else {
        0.0
    };

    let stats = CycleStatistics {
        cycle_type: CycleType::FCycle,
        total_cycles: cycle_count,
        convergence_factor,
        total_time: cycle_time,
        time_per_cycle: cycle_time / cycle_count as f64,
    };

    Ok((correction, stats))
}

/// Cycle type enumeration
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CycleType {
    /// V-cycle: efficient, most commonly used
    VCycle,
    /// W-cycle: more work per cycle, sometimes better convergence
    WCycle,
    /// F-cycle: full multigrid, theoretically optimal
    FCycle,
}

/// Statistics for multigrid cycle performance
#[derive(Debug, Clone)]
pub struct CycleStatistics {
    /// Type of cycle used
    pub cycle_type: CycleType,
    /// Total number of cycles performed
    pub total_cycles: usize,
    /// Average convergence factor per cycle
    pub convergence_factor: f64,
    /// Total time spent in cycles (seconds)
    pub total_time: f64,
    /// Average time per cycle (seconds)
    pub time_per_cycle: f64,
}

/// Solve system on coarsest level using direct or iterative method
fn solve_coarsest_level(
    matrix: &SparseMatrix<f64>,
    rhs: &DVector<f64>,
    solution: &mut DVector<f64>,
) -> Result<()> {
    let n = matrix.nrows();

    if n <= 100 {
        // For small systems, convert to dense and use Gaussian elimination
        let mut dense = DMatrix::zeros(n, n);
        for i in 0..n {
            let row = matrix.row(i);
            for (&j, &val) in row.col_indices().iter().zip(row.values().iter()) {
                dense[(i, j)] = val;
            }
        }
        gaussian_elimination_solve(&dense, rhs, solution)?;
    } else {
        // For larger systems, use iterative method
        let mut residual = rhs - matrix * &*solution;
        let tolerance = 1e-12;
        let max_iter = 100;

        for _ in 0..max_iter {
            // Simple Jacobi iteration for coarsest level
            let solution_old = solution.clone();

            for i in 0..n {
                let mut sum = 0.0;
                let row = matrix.row(i);
                let mut diag: f64 = 0.0;

                for (&j, &val) in row.col_indices().iter().zip(row.values().iter()) {
                    if i == j {
                        diag = val;
                    } else {
                        sum += val * solution_old[j];
                    }
                }

                if diag.abs() > 1e-15 {
                    solution[i] = (rhs[i] - sum) / diag;
                }
            }

            residual = rhs - matrix * &*solution;
            if residual.norm() < tolerance {
                break;
            }
        }
    }

    Ok(())
}

/// Simple Gaussian elimination solver for small matrices
fn gaussian_elimination_solve(
    matrix: &DMatrix<f64>,
    rhs: &DVector<f64>,
    solution: &mut DVector<f64>,
) -> Result<()> {
    let n = matrix.nrows();
    let mut augmented = DMatrix::zeros(n, n + 1);

    // Create augmented matrix [A | b]
    for i in 0..n {
        for j in 0..n {
            augmented[(i, j)] = matrix[(i, j)];
        }
        augmented[(i, n)] = rhs[i];
    }

    // Forward elimination
    for i in 0..n {
        // Find pivot
        let mut max_row = i;
        for k in i + 1..n {
            if augmented[(k, i)].abs() > augmented[(max_row, i)].abs() {
                max_row = k;
            }
        }

        // Swap rows
        augmented.swap_rows(i, max_row);

        // Check for singularity
        if augmented[(i, i)].abs() < 1e-15 {
            return Err(Error::Numerical(NumericalErrorKind::SingularMatrix));
        }

        // Eliminate
        for k in i + 1..n {
            let factor = augmented[(k, i)] / augmented[(i, i)];
            for j in i..=n {
                augmented[(k, j)] -= factor * augmented[(i, j)];
            }
        }
    }

    // Backward substitution
    for i in (0..n).rev() {
        let mut sum = 0.0;
        for j in i + 1..n {
            sum += augmented[(i, j)] * solution[j];
        }
        solution[i] = (augmented[(i, n)] - sum) / augmented[(i, i)];
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn compute_residual(
        matrix: &SparseMatrix<f64>,
        rhs: &DVector<f64>,
        solution: &DVector<f64>,
    ) -> DVector<f64> {
        rhs - matrix * solution
    }

    fn create_test_multigrid_level() -> MultigridLevel<f64> {
        // Create a simple 3x3 matrix
        let mut coo = nalgebra_sparse::CooMatrix::new(3, 3);
        coo.push(0, 0, 2.0);
        coo.push(0, 1, -1.0);
        coo.push(1, 0, -1.0);
        coo.push(1, 1, 2.0);
        coo.push(1, 2, -1.0);
        coo.push(2, 1, -1.0);
        coo.push(2, 2, 2.0);

        let matrix = nalgebra_sparse::CsrMatrix::from(&coo);

        // Simple Gauss-Seidel smoother
        let smoother = super::super::smoothers::GaussSeidelSmoother::new(1.0);

        MultigridLevel {
            matrix,
            restriction: None, // Single level test
            interpolation: None,
            smoother: Box::new(smoother),
        }
    }

    #[test]
    fn test_v_cycle_single_level() {
        let level = create_test_multigrid_level();
        let levels = vec![level];

        let rhs = nalgebra::DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let initial_solution = nalgebra::DVector::from_vec(vec![0.3, -0.2, 0.4]);
        let residual = compute_residual(&levels[0].matrix, &rhs, &initial_solution);

        let (correction, stats) = apply_v_cycle(&levels, &residual, 5, 1e-6).unwrap();

        // Check that we got a result
        assert_eq!(correction.len(), 3);
        assert_eq!(stats.cycle_type, CycleType::VCycle);
        assert!(stats.total_cycles > 0);
        assert!(stats.total_time > 0.0);
    }

    #[test]
    fn test_w_cycle_single_level() {
        let level = create_test_multigrid_level();
        let levels = vec![level];

        let rhs = nalgebra::DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let initial_solution = nalgebra::DVector::from_vec(vec![0.3, -0.2, 0.4]);
        let residual = compute_residual(&levels[0].matrix, &rhs, &initial_solution);

        let (correction, stats) = apply_w_cycle(&levels, &residual, 3, 1e-6).unwrap();

        assert_eq!(correction.len(), 3);
        assert_eq!(stats.cycle_type, CycleType::WCycle);
        assert!(stats.total_cycles > 0);
    }

    #[test]
    fn test_cycle_statistics() {
        let level = create_test_multigrid_level();
        let levels = vec![level];

        let rhs = nalgebra::DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let initial_solution = nalgebra::DVector::from_vec(vec![0.3, -0.2, 0.4]);
        let residual = compute_residual(&levels[0].matrix, &rhs, &initial_solution);

        // Use 0.0 tolerance to ensure it runs for all 2 cycles
        let (_, stats) = apply_v_cycle(&levels, &residual, 2, 0.0).unwrap();

        assert!(stats.time_per_cycle > 0.0);
        assert!(stats.total_time >= stats.time_per_cycle);
        assert_eq!(stats.total_cycles, 2);
    }

    #[test]
    fn test_gaussian_elimination_small_matrix() {
        let mut matrix = nalgebra::DMatrix::zeros(3, 3);
        matrix[(0, 0)] = 1.0;
        matrix[(1, 1)] = 1.0;
        matrix[(2, 2)] = 1.0;
        matrix[(0, 1)] = 1.0;
        matrix[(1, 2)] = 1.0;

        let rhs = nalgebra::DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let mut solution = nalgebra::DVector::zeros(3);

        gaussian_elimination_solve(&matrix, &rhs, &mut solution).unwrap();

        // Check that solution is reasonable (exact values depend on matrix)
        assert!(solution.iter().all(|&x| x.is_finite()));
        assert!(solution.iter().any(|&x| x.abs() > 1e-10)); // Should not be zero
    }

    #[test]
    fn test_empty_levels_error() {
        let levels = Vec::new();
        let residual = nalgebra::DVector::from_vec(vec![1.0, 2.0, 3.0]);

        let result = apply_v_cycle(&levels, &residual, 1, 1e-6);
        assert!(result.is_err());
    }
}
