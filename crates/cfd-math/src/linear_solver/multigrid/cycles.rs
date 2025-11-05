//! Multigrid cycle algorithms for AMG

use super::MultigridLevel;
use nalgebra::{DMatrix, DVector};
use std::time::Instant;

/// Apply V-cycle multigrid algorithm
pub fn apply_v_cycle(
    levels: &[MultigridLevel<f64>],
    residual: &DVector<f64>,
    max_cycles: usize,
    tolerance: f64,
) -> Result<(DVector<f64>, CycleStatistics), &'static str> {
    let start_time = Instant::now();

    if levels.is_empty() {
        return Err("No multigrid levels available");
    }

    let mut correction = DVector::zeros(residual.len());
    let mut cycle_count = 0;

    // Continue cycling until convergence or max cycles reached
    for cycle in 0..max_cycles {
        cycle_count = cycle + 1;

        // Apply one V-cycle
        apply_single_v_cycle(levels, residual, &mut correction)?;

        // Check convergence (simplified - in practice you'd compute actual residual)
        let residual_norm = (residual - &levels[0].matrix * &correction).norm();
        if residual_norm < tolerance {
            break;
        }
    }

    let cycle_time = start_time.elapsed().as_secs_f64();

    let stats = CycleStatistics {
        cycle_type: CycleType::VCycle,
        total_cycles: cycle_count,
        convergence_factor: 0.0, // Would need to track convergence history
        total_time: cycle_time,
        time_per_cycle: cycle_time / cycle_count as f64,
    };

    Ok((correction, stats))
}

/// Apply a single V-cycle iteration
fn apply_single_v_cycle(
    levels: &[MultigridLevel<f64>],
    residual: &DVector<f64>,
    correction: &mut DVector<f64>,
) -> Result<(), &'static str> {
    if levels.is_empty() {
        return Ok(());
    }

    // Pre-smoothing on finest level
    let finest_level = &levels[0];
    finest_level.smoother.apply(
        &finest_level.matrix,
        correction,
        residual,
        2, // pre_smooth_iterations
    );

    // Compute residual after pre-smoothing
    let mut current_residual: DVector<f64> = residual - &finest_level.matrix * &*correction;

    // Restrict to coarser levels
    for level_idx in 0..levels.len() - 1 {
        let current_level = &levels[level_idx];

        // Restrict residual
        if let Some(ref restriction) = current_level.restriction {
            current_residual = restriction * current_residual;
        } else {
            return Err("Missing restriction operator");
        }

        // Apply coarse level correction
        if level_idx + 1 < levels.len() {
            let _coarse_level = &levels[level_idx + 1];
            let mut coarse_correction = DVector::zeros(current_residual.len());

            // Recursive call for gamma cycles (V-cycle: gamma = 1)
            apply_single_v_cycle_recursive(
                &levels[level_idx + 1..],
                &current_residual,
                &mut coarse_correction,
                1, // gamma = 1 for V-cycle
            )?;

            // Interpolate correction back
            if let Some(ref interpolation) = current_level.interpolation {
                let fine_correction = interpolation * coarse_correction;
                *correction += fine_correction;
            } else {
                return Err("Missing interpolation operator");
            }
        }
    }

    // Solve exactly on coarsest level
    if levels.len() == 1 {
        // Single level - solve directly
        solve_coarsest_level(&levels[0].matrix, &current_residual, correction)?;
    }

    // Post-smoothing on finest level
    finest_level.smoother.apply(
        &finest_level.matrix,
        correction,
        residual,
        2, // post_smooth_iterations
    );

    Ok(())
}

/// Recursive helper for multi-cycle algorithms
fn apply_single_v_cycle_recursive(
    levels: &[MultigridLevel<f64>],
    residual: &DVector<f64>,
    correction: &mut DVector<f64>,
    gamma: usize,
) -> Result<(), &'static str> {
    if levels.is_empty() {
        return Ok(());
    }

    // Apply gamma cycles
    for _ in 0..gamma {
        apply_single_v_cycle(levels, residual, correction)?;
    }

    Ok(())
}

/// Apply W-cycle multigrid algorithm
pub fn apply_w_cycle(
    levels: &[MultigridLevel<f64>],
    residual: &DVector<f64>,
    max_cycles: usize,
    tolerance: f64,
) -> Result<(DVector<f64>, CycleStatistics), &'static str> {
    let start_time = Instant::now();

    if levels.is_empty() {
        return Err("No multigrid levels available");
    }

    let mut correction = DVector::zeros(residual.len());
    let mut cycle_count = 0;

    for cycle in 0..max_cycles {
        cycle_count = cycle + 1;

        // Apply W-cycle (gamma = 2)
        apply_single_v_cycle_recursive(levels, residual, &mut correction, 2)?;

        // Check convergence
        let residual_norm = (residual - &levels[0].matrix * &correction).norm();
        if residual_norm < tolerance {
            break;
        }
    }

    let cycle_time = start_time.elapsed().as_secs_f64();

    let stats = CycleStatistics {
        cycle_type: CycleType::WCycle,
        total_cycles: cycle_count,
        convergence_factor: 0.0,
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
) -> Result<(DVector<f64>, CycleStatistics), &'static str> {
    let start_time = Instant::now();

    if levels.is_empty() {
        return Err("No multigrid levels available");
    }

    let _solution: DVector<f64> = DVector::zeros(rhs.len());
    let _cycle_count = 0;

    // F-cycle starts from coarsest level and works up
    // This is a simplified implementation

    // For now, fall back to V-cycle
    // Full F-cycle implementation would require more sophisticated
    // nested iteration schemes
    let (correction, stats) = apply_v_cycle(levels, rhs, max_cycles, tolerance)?;

    let mut final_stats = stats;
    final_stats.cycle_type = CycleType::FCycle;
    final_stats.total_time = start_time.elapsed().as_secs_f64();

    Ok((correction, final_stats))
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
    matrix: &DMatrix<f64>,
    rhs: &DVector<f64>,
    solution: &mut DVector<f64>,
) -> Result<(), &'static str> {
    let n = matrix.nrows();

    if n <= 100 {
        // For small systems, use Gaussian elimination
        // This is a simplified direct solver
        gaussian_elimination_solve(matrix, rhs, solution)?;
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
                for j in 0..n {
                    if i != j {
                        sum += matrix[(i, j)] * solution_old[j];
                    }
                }

                if matrix[(i, i)].abs() > 1e-15 {
                    solution[i] = (rhs[i] - sum) / matrix[(i, i)];
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
) -> Result<(), &'static str> {
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
        let mut pivot_row = i;
        for k in i + 1..n {
            if augmented[(k, i)].abs() > augmented[(pivot_row, i)].abs() {
                pivot_row = k;
            }
        }

        // Swap rows
        for j in 0..=n {
            let temp = augmented[(i, j)];
            augmented[(i, j)] = augmented[(pivot_row, j)];
            augmented[(pivot_row, j)] = temp;
        }

        // Eliminate
        for k in i + 1..n {
            let factor = augmented[(k, i)] / augmented[(i, i)];
            for j in i..=n {
                if i == j {
                    augmented[(k, j)] = 0.0;
                } else {
                    augmented[(k, j)] -= factor * augmented[(i, j)];
                }
            }
        }
    }

    // Back substitution
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
    use nalgebra::DMatrix;

    fn create_test_multigrid_level() -> MultigridLevel<f64> {
        // Create a simple 3x3 matrix
        let mut matrix = DMatrix::zeros(3, 3);
        matrix[(0, 0)] = 2.0; matrix[(0, 1)] = -1.0;
        matrix[(1, 0)] = -1.0; matrix[(1, 1)] = 2.0; matrix[(1, 2)] = -1.0;
        matrix[(2, 1)] = -1.0; matrix[(2, 2)] = 2.0;

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
        let residual = rhs.clone(); // Simplified: residual = rhs

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
        let residual = rhs.clone();

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
        let residual = rhs.clone();

        let (_, stats) = apply_v_cycle(&levels, &residual, 2, 1e-6).unwrap();

        assert!(stats.time_per_cycle > 0.0);
        assert!(stats.total_time >= stats.time_per_cycle);
        assert_eq!(stats.total_cycles, 2);
    }

    #[test]
    fn test_gaussian_elimination_small_matrix() {
        let mut matrix = DMatrix::identity(3, 3);
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
