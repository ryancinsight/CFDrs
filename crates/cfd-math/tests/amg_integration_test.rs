//! Integration test for AMG preconditioner with BiCGSTAB and GMRES solvers

use cfd_math::linear_solver::preconditioners::multigrid::{
    AMGConfig, AlgebraicMultigrid, CoarseningStrategy, CycleType,
};
use cfd_math::linear_solver::{
    BiCGSTAB, IterativeLinearSolver, IterativeSolverConfig, Preconditioner, GMRES,
};
use cfd_math::sparse::spmv;
use nalgebra::{DVector, RealField};
use nalgebra_sparse::CsrMatrix;
use num_traits::FromPrimitive;

/// Create a test matrix representing a 2D Poisson equation
fn create_poisson_matrix<T: RealField + From<f64> + FromPrimitive>(n: usize) -> CsrMatrix<T> {
    let size = n * n;
    let mut values = Vec::new();
    let mut indices = Vec::new();
    let mut indptr = vec![0];

    for i in 0..size {
        let row = i / n;
        let col = i % n;

        // Top neighbor (i-n)
        if row > 0 {
            values.push(T::from_f64(-1.0).unwrap());
            indices.push(i - n);
        }

        // Left neighbor (i-1)
        if col > 0 {
            values.push(T::from_f64(-1.0).unwrap());
            indices.push(i - 1);
        }

        // Diagonal element (i)
        values.push(T::from_f64(4.0).unwrap());
        indices.push(i);

        // Right neighbor (i+1)
        if col < n - 1 {
            values.push(T::from_f64(-1.0).unwrap());
            indices.push(i + 1);
        }

        // Bottom neighbor (i+n)
        if row < n - 1 {
            values.push(T::from_f64(-1.0).unwrap());
            indices.push(i + n);
        }

        indptr.push(values.len());
    }

    CsrMatrix::try_from_csr_data(size, size, indptr, indices, values).unwrap()
}

fn create_exact_solution<T: RealField + From<f64>>(size: usize) -> DVector<T> {
    DVector::from_fn(size, |i, _| T::from_f64((i as f64).sin()).unwrap())
}

fn create_rhs<T: RealField + From<f64> + Copy>(
    matrix: &CsrMatrix<T>,
    solution: &DVector<T>,
) -> DVector<T> {
    let mut rhs = DVector::zeros(matrix.nrows());
    spmv(matrix, solution, &mut rhs);
    rhs
}

#[test]
fn test_amg_with_bicgstab() {
    let n = 8; // 8x8 grid = 64 unknowns
    let matrix = create_poisson_matrix::<f64>(n);
    let exact_solution = create_exact_solution(matrix.nrows());
    let rhs = create_rhs(&matrix, &exact_solution);

    // Create AMG preconditioner
    let amg_config = AMGConfig {
        cycle_type: CycleType::VCycle,
        pre_smooth_iterations: 2,
        post_smooth_iterations: 2,
        max_levels: 5,
        min_coarse_size: 10,
        coarsening_strategy: CoarseningStrategy::RugeStueben,
        ..Default::default()
    };
    let amg = AlgebraicMultigrid::new(&matrix, amg_config).unwrap();

    // Create BiCGSTAB solver and solve with AMG preconditioner
    let solver_config = IterativeSolverConfig {
        max_iterations: 1000,
        tolerance: 1e-8,
        ..Default::default()
    };
    let solver = BiCGSTAB::new(solver_config);

    // Solve the system with AMG preconditioning
    let mut solution = DVector::zeros(matrix.nrows());
    let result = solver.solve_preconditioned(&matrix, &rhs, &amg, &mut solution);

    if let Err(ref e) = result {
        println!("BiCGSTAB Error: {:?}", e);
    }
    assert!(result.is_ok(), "BiCGSTAB with AMG should converge");

    // Check solution accuracy
    let error = (&solution - &exact_solution).norm();
    assert!(error < 1e-6, "Solution error should be small: {}", error);
}

#[test]
fn test_amg_with_gmres() {
    let n = 6; // 6x6 grid = 36 unknowns
    let matrix = create_poisson_matrix::<f64>(n);
    let exact_solution = create_exact_solution(matrix.nrows());
    let rhs = create_rhs(&matrix, &exact_solution);

    // Create AMG preconditioner
    let amg = AlgebraicMultigrid::new(&matrix, AMGConfig::default()).unwrap();

    // Create GMRES solver with AMG preconditioner
    let solver_config = IterativeSolverConfig {
        max_iterations: 100,
        tolerance: 1e-8,
        ..Default::default()
    };
    let solver = GMRES::new(solver_config, 20); // GMRES needs restart dimension

    // Solve the system with AMG preconditioning
    let mut solution = DVector::zeros(matrix.nrows());
    let result = solver.solve_preconditioned(&matrix, &rhs, &amg, &mut solution);

    if let Err(ref e) = result {
        println!("GMRES Error: {:?}", e);
    }
    assert!(result.is_ok(), "GMRES with AMG should converge");

    // Check solution accuracy
    let error = (&solution - &exact_solution).norm();
    assert!(error < 1e-6, "Solution error should be small: {}", error);
}

#[test]
fn test_amg_different_cycles() {
    let n = 4; // 4x4 grid = 16 unknowns
    let matrix = create_poisson_matrix::<f64>(n);
    let exact_solution = create_exact_solution(matrix.nrows());
    let rhs = create_rhs(&matrix, &exact_solution);

    // Test V-cycle
    let amg_v = AlgebraicMultigrid::new(
        &matrix,
        AMGConfig {
            cycle_type: CycleType::VCycle,
            ..Default::default()
        },
    )
    .unwrap();

    // Test W-cycle
    let amg_w = AlgebraicMultigrid::new(
        &matrix,
        AMGConfig {
            cycle_type: CycleType::WCycle,
            ..Default::default()
        },
    )
    .unwrap();

    // Both should work
    let mut solution_v = DVector::zeros(matrix.nrows());
    let mut solution_w = DVector::zeros(matrix.nrows());

    amg_v.apply_to(&rhs, &mut solution_v).unwrap();
    amg_w.apply_to(&rhs, &mut solution_w).unwrap();

    // Solutions should be reasonable (not checking exactness since AMG is a preconditioner)
    assert!(
        solution_v.norm() > 0.0,
        "V-cycle solution should be non-zero"
    );
    assert!(
        solution_w.norm() > 0.0,
        "W-cycle solution should be non-zero"
    );
}

#[test]
fn test_amg_construction_edge_cases() {
    // Test with very small matrix
    let small_matrix = create_poisson_matrix::<f64>(2); // 4 unknowns
    let amg_small = AlgebraicMultigrid::new(&small_matrix, AMGConfig::default());
    assert!(amg_small.is_ok(), "AMG should handle small matrices");

    // Test with larger matrix
    let large_matrix = create_poisson_matrix::<f64>(16); // 256 unknowns
    let amg_large = AlgebraicMultigrid::new(&large_matrix, AMGConfig::default());
    assert!(amg_large.is_ok(), "AMG should handle larger matrices");

    // Test custom configuration
    let custom_config = AMGConfig {
        cycle_type: CycleType::FCycle,
        pre_smooth_iterations: 3,
        post_smooth_iterations: 3,
        max_levels: 8,
        min_coarse_size: 5,
        coarsening_strategy: CoarseningStrategy::RugeStueben, // Classical was renamed or doesn't exist
        ..Default::default()
    };
    let amg_custom = AlgebraicMultigrid::new(&large_matrix, custom_config);
    assert!(
        amg_custom.is_ok(),
        "AMG should accept custom configurations"
    );
}
