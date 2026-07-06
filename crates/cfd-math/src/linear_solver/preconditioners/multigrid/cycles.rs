//! Multigrid cycle algorithms for AMG

use super::{MultigridLevel, MultigridVector, SparseMatrix};
use crate::linear_solver::dense_bridge::solve_leto_csr_with_leto_dense_array;
use cfd_core::error::{Error, Result};
use leto_ops::spmv as leto_spmv;
use std::time::Instant;

fn l2_norm(vector: &MultigridVector<f64>) -> f64 {
    let mut sum = 0.0;
    for i in 0..vector.shape()[0] {
        sum += vector[i] * vector[i];
    }
    sum.sqrt()
}

fn vector_sub(lhs: &MultigridVector<f64>, rhs: &MultigridVector<f64>) -> MultigridVector<f64> {
    let mut output = MultigridVector::zeros([lhs.shape()[0]]);
    for i in 0..lhs.shape()[0] {
        output[i] = lhs[i] - rhs[i];
    }
    output
}

fn vector_add_assign(lhs: &mut MultigridVector<f64>, rhs: &MultigridVector<f64>) {
    for i in 0..lhs.shape()[0] {
        lhs[i] += rhs[i];
    }
}

fn residual(
    matrix: &SparseMatrix<f64>,
    rhs: &MultigridVector<f64>,
    solution: &MultigridVector<f64>,
) -> MultigridVector<f64> {
    let applied =
        leto_spmv(matrix, &solution.view()).expect("invariant: multigrid SpMV inputs are valid");
    vector_sub(rhs, &applied)
}

fn sparse_vector_mul(
    matrix: &SparseMatrix<f64>,
    vector: &MultigridVector<f64>,
) -> MultigridVector<f64> {
    leto_spmv(matrix, &vector.view()).expect("invariant: multigrid transfer SpMV inputs are valid")
}

/// Apply V-cycle multigrid algorithm
pub fn apply_v_cycle(
    levels: &[MultigridLevel<f64>],
    residual: &MultigridVector<f64>,
    max_cycles: usize,
    tolerance: f64,
) -> Result<(MultigridVector<f64>, CycleStatistics)> {
    let start_time = Instant::now();

    if levels.is_empty() {
        return Err(cfd_core::error::Error::InvalidConfiguration(
            "No multigrid levels available".to_string(),
        ));
    }

    let mut correction = MultigridVector::zeros([residual.shape()[0]]);
    let mut cycle_count = 0;
    let mut residual_history: Vec<f64> = Vec::new();

    // Continue cycling until convergence or max cycles reached
    for cycle in 0..max_cycles {
        cycle_count = cycle + 1;

        // Apply one V-cycle
        apply_multigrid_cycle(1, levels, residual, &mut correction)?;

        let residual_norm = l2_norm(&self::residual(&levels[0].matrix, residual, &correction));
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
    residual: &MultigridVector<f64>,
    correction: &mut MultigridVector<f64>,
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
    let r_fine = self::residual(&current_level.matrix, residual, correction);

    // 3. Restrict residual
    let r_coarse = if let Some(ref restriction) = current_level.restriction {
        sparse_vector_mul(restriction, &r_fine)
    } else {
        return Err(Error::InvalidConfiguration(
            "Missing restriction operator".to_string(),
        ));
    };

    // 4. Recursive call for coarser levels
    let mut coarse_correction = MultigridVector::zeros([r_coarse.shape()[0]]);
    for _ in 0..gamma {
        apply_multigrid_cycle(gamma, &levels[1..], &r_coarse, &mut coarse_correction)?;
    }

    // 5. Interpolate correction back
    if let Some(ref interpolation) = current_level.interpolation {
        let fine_correction = sparse_vector_mul(interpolation, &coarse_correction);
        vector_add_assign(correction, &fine_correction);
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
    residual: &MultigridVector<f64>,
    max_cycles: usize,
    tolerance: f64,
) -> Result<(MultigridVector<f64>, CycleStatistics)> {
    let start_time = Instant::now();

    if levels.is_empty() {
        return Err(cfd_core::error::Error::InvalidConfiguration(
            "No multigrid levels available".to_string(),
        ));
    }

    let mut correction = MultigridVector::zeros([residual.shape()[0]]);
    let mut cycle_count = 0;
    let mut residual_history: Vec<f64> = Vec::new();

    for cycle in 0..max_cycles {
        cycle_count = cycle + 1;

        // Apply W-cycle (gamma = 2)
        apply_multigrid_cycle(2, levels, residual, &mut correction)?;

        // Check convergence
        let residual_norm = l2_norm(&self::residual(&levels[0].matrix, residual, &correction));
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
    rhs: &MultigridVector<f64>,
    max_cycles: usize,
    tolerance: f64,
) -> Result<(MultigridVector<f64>, CycleStatistics)> {
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

    let mut rhs_levels: Vec<MultigridVector<f64>> = Vec::with_capacity(levels.len());
    rhs_levels.push(rhs.clone());

    for level_idx in 0..levels.len() - 1 {
        let restriction = levels[level_idx].restriction.as_ref().ok_or_else(|| {
            Error::InvalidConfiguration("Missing restriction operator".to_string())
        })?;
        let coarse_rhs = sparse_vector_mul(restriction, &rhs_levels[level_idx]);
        rhs_levels.push(coarse_rhs);
    }

    let coarsest_idx = levels.len() - 1;
    let mut correction = MultigridVector::zeros([rhs_levels[coarsest_idx].shape()[0]]);
    solve_coarsest_level(
        &levels[coarsest_idx].matrix,
        &rhs_levels[coarsest_idx],
        &mut correction,
    )?;

    let mut cycle_count = 0;
    let mut residual_history: Vec<f64> = Vec::new();

    for level_idx in (0..coarsest_idx).rev() {
        let interpolation = levels[level_idx].interpolation.as_ref().ok_or_else(|| {
            Error::InvalidConfiguration("Missing interpolation operator".to_string())
        })?;
        let mut fine_correction = sparse_vector_mul(interpolation, &correction);

        apply_multigrid_cycle(
            1,
            &levels[level_idx..],
            &rhs_levels[level_idx],
            &mut fine_correction,
        )?;

        if level_idx == 0 {
            cycle_count += 1;
            let mut residual_norm = l2_norm(&self::residual(
                &levels[0].matrix,
                &rhs_levels[0],
                &fine_correction,
            ));
            residual_history.push(residual_norm);

            while cycle_count < max_cycles && residual_norm >= tolerance {
                apply_multigrid_cycle(
                    1,
                    &levels[level_idx..],
                    &rhs_levels[level_idx],
                    &mut fine_correction,
                )?;
                cycle_count += 1;
                residual_norm = l2_norm(&self::residual(
                    &levels[0].matrix,
                    &rhs_levels[0],
                    &fine_correction,
                ));
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
    rhs: &MultigridVector<f64>,
    solution: &mut MultigridVector<f64>,
) -> Result<()> {
    let n = matrix.nrows();

    if n <= 100 {
        let dense_solution = solve_leto_csr_with_leto_dense_array(matrix, rhs)?;
        for i in 0..solution.shape()[0] {
            solution[i] = dense_solution[i];
        }
    } else {
        // For larger systems, use iterative method
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

            let residual = self::residual(matrix, rhs, solution);
            if l2_norm(&residual) < tolerance {
                break;
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::csr_from_parts;
    use super::*;

    fn create_test_multigrid_level() -> MultigridLevel<f64> {
        let matrix = csr_from_parts(
            3,
            3,
            vec![0, 2, 5, 7],
            vec![0, 1, 0, 1, 2, 1, 2],
            vec![2.0, -1.0, -1.0, 2.0, -1.0, -1.0, 2.0],
            "cycle test matrix",
        )
        .unwrap();

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

        let rhs = MultigridVector::from_shape_vec([3], vec![1.0, 2.0, 3.0]).unwrap();
        let initial_solution = MultigridVector::from_shape_vec([3], vec![0.3, -0.2, 0.4]).unwrap();
        let residual = self::residual(&levels[0].matrix, &rhs, &initial_solution);

        let (correction, stats) = apply_v_cycle(&levels, &residual, 5, 1e-6).unwrap();

        // Check that we got a result
        assert_eq!(correction.shape(), [3]);
        assert_eq!(stats.cycle_type, CycleType::VCycle);
        assert!(stats.total_cycles > 0);
        assert!(stats.total_time > 0.0);
    }

    #[test]
    fn test_w_cycle_single_level() {
        let level = create_test_multigrid_level();
        let levels = vec![level];

        let rhs = MultigridVector::from_shape_vec([3], vec![1.0, 2.0, 3.0]).unwrap();
        let initial_solution = MultigridVector::from_shape_vec([3], vec![0.3, -0.2, 0.4]).unwrap();
        let residual = self::residual(&levels[0].matrix, &rhs, &initial_solution);

        let (correction, stats) = apply_w_cycle(&levels, &residual, 3, 1e-6).unwrap();

        assert_eq!(correction.shape(), [3]);
        assert_eq!(stats.cycle_type, CycleType::WCycle);
        assert!(stats.total_cycles > 0);
    }

    #[test]
    fn test_cycle_statistics() {
        let level = create_test_multigrid_level();
        let levels = vec![level];

        let rhs = MultigridVector::from_shape_vec([3], vec![1.0, 2.0, 3.0]).unwrap();
        let initial_solution = MultigridVector::from_shape_vec([3], vec![0.3, -0.2, 0.4]).unwrap();
        let residual = self::residual(&levels[0].matrix, &rhs, &initial_solution);

        // Use 0.0 tolerance to ensure it runs for all 2 cycles
        let (_, stats) = apply_v_cycle(&levels, &residual, 2, 0.0).unwrap();

        assert!(stats.time_per_cycle > 0.0);
        assert!(stats.total_time >= stats.time_per_cycle);
        assert_eq!(stats.total_cycles, 2);
    }

    #[test]
    fn test_leto_coarsest_solve_small_matrix() {
        let matrix = csr_from_parts(
            3,
            3,
            vec![0, 2, 4, 5],
            vec![0, 1, 1, 2, 2],
            vec![1.0, 1.0, 1.0, 1.0, 1.0],
            "cycle coarsest test matrix",
        )
        .unwrap();

        let rhs = MultigridVector::from_shape_vec([3], vec![1.0, 2.0, 3.0]).unwrap();
        let mut solution = MultigridVector::zeros([3]);

        solve_coarsest_level(&matrix, &rhs, &mut solution).unwrap();

        assert!((solution[0] - 2.0_f64).abs() < 1e-12_f64);
        assert!((solution[1] + 1.0_f64).abs() < 1e-12_f64);
        assert!((solution[2] - 3.0_f64).abs() < 1e-12_f64);
    }

    #[test]
    fn test_empty_levels_error() {
        let levels = Vec::new();
        let residual = MultigridVector::from_shape_vec([3], vec![1.0, 2.0, 3.0]).unwrap();

        let result = apply_v_cycle(&levels, &residual, 1, 1e-6);
        assert!(result.is_err());
    }
}
