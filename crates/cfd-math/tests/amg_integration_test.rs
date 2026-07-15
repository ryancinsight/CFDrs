//! Integration test for AMG preconditioner with BiCGSTAB and GMRES solvers

use cfd_math::linear_solver::preconditioners::multigrid::{
    AMGConfig, AlgebraicMultigrid, CoarseningStrategy, CycleType,
};
use cfd_math::linear_solver::{BiCGSTAB, IterativeSolverConfig, Preconditioner, GMRES};
use cfd_math::sparse::{spmv, SparseMatrix, SparseMatrixBuilder};
use eunomia::{FloatElement, RealField};
use leto::Array1;
use leto_ops::{CsrMatrix as AtlasSparseMatrix, Scalar as LetoScalar};

struct PoissonSystem<T: RealField + Copy + LetoScalar> {
    solver_matrix: SparseMatrix<T>,
    amg_matrix: AtlasSparseMatrix<T>,
}

/// Create a test matrix representing a 2D Poisson equation
fn create_poisson_system<T: RealField + FloatElement + LetoScalar>(n: usize) -> PoissonSystem<T> {
    let size = n * n;
    let mut builder = SparseMatrixBuilder::new(size, size);
    let mut values = Vec::new();
    let mut indices = Vec::new();
    let mut indptr = vec![0];

    for i in 0..size {
        let row = i / n;
        let col = i % n;

        // Top neighbor (i-n)
        if row > 0 {
            let value = <T as FloatElement>::from_f64(-1.0);
            builder
                .add_entry(i, i - n, value)
                .expect("invariant: top-neighbor entry is in bounds");
            values.push(value);
            indices.push(i - n);
        }

        // Left neighbor (i-1)
        if col > 0 {
            let value = <T as FloatElement>::from_f64(-1.0);
            builder
                .add_entry(i, i - 1, value)
                .expect("invariant: left-neighbor entry is in bounds");
            values.push(value);
            indices.push(i - 1);
        }

        // Diagonal element (i)
        let diagonal = <T as FloatElement>::from_f64(4.0);
        builder
            .add_entry(i, i, diagonal)
            .expect("invariant: diagonal entry is in bounds");
        values.push(diagonal);
        indices.push(i);

        // Right neighbor (i+1)
        if col < n - 1 {
            let value = <T as FloatElement>::from_f64(-1.0);
            builder
                .add_entry(i, i + 1, value)
                .expect("invariant: right-neighbor entry is in bounds");
            values.push(value);
            indices.push(i + 1);
        }

        // Bottom neighbor (i+n)
        if row < n - 1 {
            let value = <T as FloatElement>::from_f64(-1.0);
            builder
                .add_entry(i, i + n, value)
                .expect("invariant: bottom-neighbor entry is in bounds");
            values.push(value);
            indices.push(i + n);
        }

        indptr.push(values.len());
    }

    let solver_matrix = builder
        .build()
        .expect("invariant: generated Poisson entries build a valid solver matrix");
    let amg_matrix = AtlasSparseMatrix::from_parts(values, indices, indptr, size, size)
        .expect("invariant: generated Poisson CSR is valid for AMG matrix");

    PoissonSystem {
        solver_matrix,
        amg_matrix,
    }
}

fn create_exact_solution<T: FloatElement>(size: usize) -> Array1<T> {
    Array1::from_shape_vec(
        [size],
        (0..size)
            .map(|i| <T as FloatElement>::from_f64((i as f64).sin()))
            .collect(),
    )
    .expect("invariant: exact-solution shape matches data")
}

fn create_rhs<T: RealField + Copy + LetoScalar>(
    matrix: &SparseMatrix<T>,
    solution: &Array1<T>,
) -> Array1<T> {
    let mut rhs = Array1::zeros([matrix.nrows()]);
    spmv(matrix, solution, &mut rhs);
    rhs
}

fn apply_preconditioner<P: Preconditioner<f64>>(
    preconditioner: &P,
    rhs: &Array1<f64>,
    out: &mut Array1<f64>,
) {
    preconditioner
        .apply_to(rhs, out)
        .expect("preconditioner application succeeds");
}

fn vector_norm(values: &Array1<f64>) -> f64 {
    (0..values.shape()[0])
        .map(|idx| values[idx] * values[idx])
        .sum::<f64>()
        .sqrt()
}

fn vector_distance(lhs: &Array1<f64>, rhs: &Array1<f64>) -> f64 {
    assert_eq!(
        lhs.shape(),
        rhs.shape(),
        "invariant: compared arrays share shape"
    );
    (0..lhs.shape()[0])
        .map(|idx| {
            let diff = lhs[idx] - rhs[idx];
            diff * diff
        })
        .sum::<f64>()
        .sqrt()
}

fn energy_norm(matrix: &SparseMatrix<f64>, values: &Array1<f64>) -> f64 {
    let mut applied = Array1::zeros([matrix.nrows()]);
    spmv(matrix, values, &mut applied);
    let dot = (0..values.shape()[0])
        .map(|idx| values[idx] * applied[idx])
        .sum::<f64>();
    dot.max(0.0).sqrt()
}

#[test]
fn test_amg_with_bicgstab() {
    let n = 8; // 8x8 grid = 64 unknowns
    let system = create_poisson_system::<f64>(n);
    let exact_solution = create_exact_solution(system.solver_matrix.nrows());
    let rhs = create_rhs(&system.solver_matrix, &exact_solution);

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
    let amg = AlgebraicMultigrid::new(&system.amg_matrix, amg_config).unwrap();

    // Create BiCGSTAB solver and solve with AMG preconditioner
    let solver_config = IterativeSolverConfig {
        max_iterations: 1000,
        tolerance: 1e-8,
        ..Default::default()
    };
    let solver = BiCGSTAB::new(solver_config);

    // Solve the system with AMG preconditioning
    let mut solution_array = Array1::zeros([system.solver_matrix.nrows()]);
    let result =
        solver.solve_preconditioned(&system.solver_matrix, &rhs, &amg, &mut solution_array);

    if let Err(ref e) = result {
        println!("BiCGSTAB Error: {:?}", e);
    }
    assert!(result.is_ok(), "BiCGSTAB with AMG should converge");

    // Check solution accuracy
    let error = vector_distance(&solution_array, &exact_solution);
    assert!(error < 1e-6, "Solution error should be small: {}", error);
}

#[test]
fn test_amg_with_gmres() {
    let n = 6; // 6x6 grid = 36 unknowns
    let system = create_poisson_system::<f64>(n);
    let exact_solution = create_exact_solution(system.solver_matrix.nrows());
    let rhs = create_rhs(&system.solver_matrix, &exact_solution);

    // Create AMG preconditioner
    let amg = AlgebraicMultigrid::new(&system.amg_matrix, AMGConfig::default()).unwrap();

    // Create GMRES solver with AMG preconditioner
    let solver_config = IterativeSolverConfig {
        max_iterations: 100,
        tolerance: 1e-8,
        ..Default::default()
    };
    let solver = GMRES::new(solver_config, 20); // GMRES needs restart dimension

    // Solve the system with AMG preconditioning
    let mut solution_array = Array1::zeros([system.solver_matrix.nrows()]);
    let result =
        solver.solve_preconditioned(&system.solver_matrix, &rhs, &amg, &mut solution_array);

    if let Err(ref e) = result {
        println!("GMRES Error: {:?}", e);
    }
    assert!(result.is_ok(), "GMRES with AMG should converge");

    // Check solution accuracy
    let error = vector_distance(&solution_array, &exact_solution);
    assert!(error < 1e-6, "Solution error should be small: {}", error);
}

#[test]
fn test_amg_different_cycles() {
    let n = 4; // 4x4 grid = 16 unknowns
    let system = create_poisson_system::<f64>(n);
    let exact_solution = create_exact_solution(system.solver_matrix.nrows());
    let rhs = create_rhs(&system.solver_matrix, &exact_solution);

    // Test V-cycle
    let amg_v = AlgebraicMultigrid::new(
        &system.amg_matrix,
        AMGConfig {
            cycle_type: CycleType::VCycle,
            ..Default::default()
        },
    )
    .unwrap();

    // Test W-cycle
    let amg_w = AlgebraicMultigrid::new(
        &system.amg_matrix,
        AMGConfig {
            cycle_type: CycleType::WCycle,
            ..Default::default()
        },
    )
    .unwrap();

    // Both should work
    let mut solution_v = Array1::zeros([system.solver_matrix.nrows()]);
    let mut solution_w = Array1::zeros([system.solver_matrix.nrows()]);

    apply_preconditioner(&amg_v, &rhs, &mut solution_v);
    apply_preconditioner(&amg_w, &rhs, &mut solution_w);

    // Solutions should be reasonable (not checking exactness since AMG is a preconditioner)
    assert!(
        vector_norm(&solution_v) > 0.0,
        "V-cycle solution should be non-zero"
    );
    assert!(
        vector_norm(&solution_w) > 0.0,
        "W-cycle solution should be non-zero"
    );
}

#[test]
fn test_amg_construction_edge_cases() {
    // Test with very small matrix
    let small_system = create_poisson_system::<f64>(2); // 4 unknowns
    let amg_small = AlgebraicMultigrid::new(&small_system.amg_matrix, AMGConfig::default());
    assert!(amg_small.is_ok(), "AMG should handle small matrices");

    // Test with larger matrix
    let large_system = create_poisson_system::<f64>(16); // 256 unknowns
    let amg_large = AlgebraicMultigrid::new(&large_system.amg_matrix, AMGConfig::default());
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
    let amg_custom = AlgebraicMultigrid::new(&large_system.amg_matrix, custom_config);
    assert!(
        amg_custom.is_ok(),
        "AMG should accept custom configurations"
    );
}

#[test]
fn test_amg_two_grid_convergence_factor() {
    let n = 6;
    let system = create_poisson_system::<f64>(n);
    let amg = AlgebraicMultigrid::new(
        &system.amg_matrix,
        AMGConfig {
            max_levels: 2,
            min_coarse_size: 4,
            cycle_type: CycleType::VCycle,
            pre_smooth_iterations: 2,
            post_smooth_iterations: 2,
            ..Default::default()
        },
    )
    .unwrap();
    let size = system.solver_matrix.nrows();
    let mut worst_observed_ratio = 0.0f64;
    for i in 0..size {
        let mut e = Array1::zeros([size]);
        e[i] = 1.0;
        let mut r_array = Array1::zeros([size]);
        spmv(&system.solver_matrix, &e, &mut r_array);
        let mut z = Array1::zeros([size]);
        apply_preconditioner(&amg, &r_array, &mut z);
        let mut e_next = Array1::zeros([size]);
        for row in 0..size {
            e_next[row] = e[row] - z[row];
        }
        let e_norm = energy_norm(&system.solver_matrix, &e);
        let e_next_norm = energy_norm(&system.solver_matrix, &e_next);
        if e_norm > 0.0 {
            worst_observed_ratio = worst_observed_ratio.max(e_next_norm / e_norm);
        }
    }
    assert!(
        worst_observed_ratio < 1.0,
        "Two-grid smoothing should reduce A-norm basis errors; worst ratio {}",
        worst_observed_ratio
    );
}
