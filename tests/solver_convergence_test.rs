//! Test numerical solver convergence

use cfd_math::{
    sparse::{SparseMatrixBuilder, SparseMatrixExt},
    linear_solver::{LinearSolver, ConjugateGradient},
};
use cfd_core::solver::LinearSolverConfig;
use nalgebra::DVector;

#[test]
fn test_cg_convergence_on_poisson() {
    // Create 2D Poisson problem: -∇²u = f
    // Discretized as Au = b
    let n = 20; // 20x20 grid
    let n_total = n * n;
    let h = 1.0 / (n as f64 + 1.0);
    let h2 = h * h;
    
    // Build 5-point stencil matrix
    let mut builder = SparseMatrixBuilder::new(n_total, n_total);
    
    for i in 0..n {
        for j in 0..n {
            let idx = i * n + j;
            
            // Diagonal entry
            builder.add_entry(idx, idx, 4.0 / h2).unwrap();
            
            // Off-diagonal entries
            if i > 0 {
                builder.add_entry(idx, idx - n, -1.0 / h2).unwrap();
            }
            if i < n - 1 {
                builder.add_entry(idx, idx + n, -1.0 / h2).unwrap();
            }
            if j > 0 {
                builder.add_entry(idx, idx - 1, -1.0 / h2).unwrap();
            }
            if j < n - 1 {
                builder.add_entry(idx, idx + 1, -1.0 / h2).unwrap();
            }
        }
    }
    
    let a = builder.build().unwrap();
    
    // Create RHS vector (unit source)
    let b = DVector::from_element(n_total, 1.0);
    
    // Configure solver with reasonable tolerance
    let config = LinearSolverConfig {
        tolerance: 1e-6,
        max_iterations: 1000,
        ..Default::default()
    };
    
    let solver = ConjugateGradient::new(config);
    
    // Solve the system
    let x = solver.solve(&a, &b, None).expect("Solver should converge");
    
    // Verify solution by checking residual
    let residual = &b - &a * &x;
    let residual_norm = residual.norm();
    
    assert!(residual_norm < 1e-5, "Residual too large: {}", residual_norm);
    
    // Check that solution is non-trivial
    let solution_norm = x.norm();
    assert!(solution_norm > 0.01, "Solution is too small: {}", solution_norm);
}

#[test]
fn test_solver_with_diagonal_preconditioner() {
    use cfd_math::linear_solver::{JacobiPreconditioner, Preconditioner};
    
    // Create simple SPD system
    let n = 100;
    let mut builder = SparseMatrixBuilder::new(n, n);
    
    // Tridiagonal matrix
    for i in 0..n {
        builder.add_entry(i, i, 4.0).unwrap();
        if i > 0 {
            builder.add_entry(i, i - 1, -1.0).unwrap();
        }
        if i < n - 1 {
            builder.add_entry(i, i + 1, -1.0).unwrap();
        }
    }
    
    let a = builder.build().unwrap();
    let b = DVector::from_element(n, 1.0);
    
    // Create Jacobi preconditioner
    let precond = JacobiPreconditioner::new(&a).expect("Should create preconditioner");
    
    let config = LinearSolverConfig {
        tolerance: 1e-8,
        max_iterations: 500,
        ..Default::default()
    };
    
    let solver = ConjugateGradient::new(config);
    
    // Solve with preconditioning
    let x = solver.solve_preconditioned(&a, &b, &precond, None)
        .expect("Preconditioned solver should converge");
    
    // Check residual
    let residual = &b - &a * &x;
    assert!(residual.norm() < 1e-7, "Preconditioned residual too large");
}

#[test]
fn test_solver_robustness() {
    // Test with ill-conditioned matrix
    let n = 50;
    let mut builder = SparseMatrixBuilder::new(n, n);
    
    // Create matrix with varying diagonal dominance
    for i in 0..n {
        let diag_value = 2.0 + (i as f64) * 0.1; // Varying diagonal
        builder.add_entry(i, i, diag_value).unwrap();
        
        if i > 0 {
            builder.add_entry(i, i - 1, -0.5).unwrap();
        }
        if i < n - 1 {
            builder.add_entry(i, i + 1, -0.5).unwrap();
        }
    }
    
    let a = builder.build().unwrap();
    
    // Create smooth RHS
    let b: DVector<f64> = DVector::from_fn(n, |i, _| {
        ((i as f64) * std::f64::consts::PI / (n as f64)).sin()
    });
    
    let config = LinearSolverConfig {
        tolerance: 1e-5,
        max_iterations: 2000, // More iterations for ill-conditioned system
        ..Default::default()
    };
    
    let solver = ConjugateGradient::new(config);
    
    // Should still converge, even if slowly
    let result = solver.solve(&a, &b, None);
    assert!(result.is_ok(), "Solver should handle ill-conditioned matrix");
    
    if let Ok(x) = result {
        let residual = &b - &a * &x;
        assert!(residual.norm() < 1e-4, "Residual acceptable for ill-conditioned system");
    }
}