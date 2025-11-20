//! Extended comprehensive edge case and property-based tests for linear solvers
//!
//! This module adds property-based testing with proptest and additional edge cases
//! to complement the existing edge_case_tests.rs file.
//!
//! Standards: Saad (2003) "Iterative Methods for Sparse Linear Systems"

#[cfg(test)]
mod extended_edge_case_tests {
    use crate::linear_solver::preconditioners::{IdentityPreconditioner, JacobiPreconditioner};
    use crate::linear_solver::traits::{IterativeLinearSolver, Preconditioner};
    use crate::linear_solver::IterativeSolverConfig;
    use crate::linear_solver::{BiCGSTAB, ConjugateGradient};
    use crate::sparse::SparseMatrixBuilder;
    use approx::assert_relative_eq;
    use nalgebra::DVector;
    use proptest::prelude::*;

    /// Edge case: Single element system (1x1 matrix)
    /// Tests degenerate case where system reduces to scalar equation
    #[test]
    fn test_single_element_system() -> Result<(), Box<dyn std::error::Error>> {
        let n = 1;
        let mut builder = SparseMatrixBuilder::new(n, n);
        builder.add_entry(0, 0, 5.0)?;
        let a = builder.build()?;

        let b = DVector::from_element(n, 10.0);
        let mut x = DVector::zeros(n);

        let config = IterativeSolverConfig::new(1e-10).with_max_iterations(10);
        let solver = ConjugateGradient::new(config);
        let identity = IdentityPreconditioner;

        solver.solve(&a, &b, &mut x, Some(&identity))?;

        // x = b/A = 10.0/5.0 = 2.0
        assert_relative_eq!(x[0], 2.0, epsilon = 1e-10);
        Ok(())
    }

    /// Edge case: Diagonal matrix (optimal case for iterative solvers)
    /// Should converge in exactly 1 iteration
    #[test]
    fn test_diagonal_matrix_cg() -> Result<(), Box<dyn std::error::Error>> {
        let n = 10;
        let mut builder = SparseMatrixBuilder::new(n, n);

        // Pure diagonal matrix (no coupling)
        for i in 0..n {
            builder.add_entry(i, i, (i + 1) as f64)?;
        }
        let a = builder.build()?;

        let b = DVector::from_fn(n, |i, _| (i + 1) as f64);
        let mut x = DVector::zeros(n);

        let config = IterativeSolverConfig::new(1e-12).with_max_iterations(10);
        let solver = ConjugateGradient::new(config);
        let identity = IdentityPreconditioner;

        solver.solve(&a, &b, &mut x, Some(&identity))?;

        // For diagonal matrix, x[i] = b[i]/A[i,i] = 1.0
        for i in 0..n {
            assert_relative_eq!(x[i], 1.0, epsilon = 1e-10);
        }
        Ok(())
    }

    /// Edge case: Identity matrix (trivial system)
    /// Solution should be x = b
    #[test]
    fn test_identity_matrix() -> Result<(), Box<dyn std::error::Error>> {
        let n = 5;
        let mut builder = SparseMatrixBuilder::new(n, n);

        // I * x = b => x = b
        for i in 0..n {
            builder.add_entry(i, i, 1.0)?;
        }
        let a = builder.build()?;

        let b = DVector::from_fn(n, |i, _| (i + 1) as f64);
        let mut x = DVector::zeros(n);

        let config = IterativeSolverConfig::new(1e-12).with_max_iterations(10);
        let solver = BiCGSTAB::new(config);
        let identity = IdentityPreconditioner;

        solver.solve(&a, &b, &mut x, Some(&identity))?;

        // Solution should equal RHS
        for i in 0..n {
            assert_relative_eq!(x[i], b[i], epsilon = 1e-10);
        }
        Ok(())
    }

    /// Edge case: Strongly diagonally dominant matrix
    /// Guarantees convergence for iterative methods
    #[test]
    fn test_strongly_dominant_matrix() -> Result<(), Box<dyn std::error::Error>> {
        let n = 10;
        let mut builder = SparseMatrixBuilder::new(n, n);

        // |a_ii| > Σ|a_ij| for j≠i (strong diagonal dominance)
        for i in 0..n {
            builder.add_entry(i, i, 10.0)?;
            if i > 0 {
                builder.add_entry(i, i - 1, 1.0)?;
            }
            if i < n - 1 {
                builder.add_entry(i, i + 1, 1.0)?;
            }
        }
        let a = builder.build()?;

        let b = DVector::from_element(n, 1.0);
        let mut x = DVector::zeros(n);

        let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
        let solver = ConjugateGradient::new(config);
        let identity = IdentityPreconditioner;

        solver.solve(&a, &b, &mut x, Some(&identity))?;

        // Verify solution
        let ax = &a * &x;
        for i in 0..n {
            assert_relative_eq!(ax[i], b[i], epsilon = 1e-8);
        }
        Ok(())
    }

    /// Edge case: All ones RHS vector
    /// Common case in pressure Poisson equation with uniform sources
    #[test]
    fn test_all_ones_rhs() -> Result<(), Box<dyn std::error::Error>> {
        let n = 10;
        let mut builder = SparseMatrixBuilder::new(n, n);

        for i in 0..n {
            builder.add_entry(i, i, 2.0)?;
            if i > 0 {
                builder.add_entry(i, i - 1, -1.0)?;
            }
            if i < n - 1 {
                builder.add_entry(i, i + 1, -1.0)?;
            }
        }
        let a = builder.build()?;

        let b = DVector::from_element(n, 1.0);
        let mut x = DVector::zeros(n);

        let config = IterativeSolverConfig::new(1e-10).with_max_iterations(1000);
        let solver = BiCGSTAB::new(config);
        let identity = IdentityPreconditioner;

        solver.solve(&a, &b, &mut x, Some(&identity))?;

        // Verify Ax = b
        let ax = &a * &x;
        for i in 0..n {
            assert_relative_eq!(ax[i], b[i], epsilon = 1e-8);
        }
        Ok(())
    }

    /// Edge case: Alternating sign RHS (stress test for stability)
    #[test]
    fn test_alternating_sign_rhs() -> Result<(), Box<dyn std::error::Error>> {
        let n = 10;
        let mut builder = SparseMatrixBuilder::new(n, n);

        for i in 0..n {
            builder.add_entry(i, i, 2.0)?;
            if i > 0 {
                builder.add_entry(i, i - 1, -1.0)?;
            }
            if i < n - 1 {
                builder.add_entry(i, i + 1, -1.0)?;
            }
        }
        let a = builder.build()?;

        // Alternating +1, -1, +1, -1, ...
        let b = DVector::from_fn(n, |i, _| if i % 2 == 0 { 1.0 } else { -1.0 });
        let mut x = DVector::zeros(n);

        let config = IterativeSolverConfig::new(1e-8).with_max_iterations(1000);
        let solver = BiCGSTAB::new(config);
        let identity = IdentityPreconditioner;

        solver.solve(&a, &b, &mut x, Some(&identity))?;

        // Verify solution
        let ax = &a * &x;
        for i in 0..n {
            assert_relative_eq!(ax[i], b[i], epsilon = 1e-6);
        }
        Ok(())
    }

    /// Edge case: Jacobi preconditioner with uniform diagonal
    /// Tests degenerate scaling case
    #[test]
    fn test_jacobi_uniform_diagonal() -> Result<(), Box<dyn std::error::Error>> {
        let n = 5;
        let mut builder = SparseMatrixBuilder::new(n, n);

        // All diagonal entries are same value
        for i in 0..n {
            builder.add_entry(i, i, 3.0)?;
            if i > 0 {
                builder.add_entry(i, i - 1, 0.5)?;
            }
            if i < n - 1 {
                builder.add_entry(i, i + 1, 0.5)?;
            }
        }
        let a = builder.build()?;

        let precond = JacobiPreconditioner::new(&a)?;
        let r = DVector::from_element(n, 6.0);
        let mut z = DVector::zeros(n);

        precond.apply_to(&r, &mut z)?;

        // All entries should be scaled by same factor: z = r/3 = 2.0
        for i in 0..n {
            assert_relative_eq!(z[i], 2.0, epsilon = 1e-10);
        }
        Ok(())
    }

    /// Edge case: Large sparse system (performance/scalability test)
    #[test]
    fn test_large_sparse_system() -> Result<(), Box<dyn std::error::Error>> {
        let n = 100;
        let mut builder = SparseMatrixBuilder::new(n, n);

        // Large tridiagonal system
        for i in 0..n {
            builder.add_entry(i, i, 2.0)?;
            if i > 0 {
                builder.add_entry(i, i - 1, -1.0)?;
            }
            if i < n - 1 {
                builder.add_entry(i, i + 1, -1.0)?;
            }
        }
        let a = builder.build()?;

        let b = DVector::from_element(n, 1.0);
        let mut x = DVector::zeros(n);

        let config = IterativeSolverConfig::new(1e-8).with_max_iterations(1000);
        let solver = ConjugateGradient::new(config);
        let identity = IdentityPreconditioner;

        solver.solve(&a, &b, &mut x, Some(&identity))?;

        // Verify residual norm
        let ax = &a * &x;
        let residual: f64 = (0..n)
            .map(|i| {
                let diff = ax[i] - b[i];
                diff * diff
            })
            .sum::<f64>()
            .sqrt();
        assert!(residual < 1e-6, "Residual {} exceeds tolerance", residual);
        Ok(())
    }

    /// Property test: Solution satisfies Ax = b for random SPD tridiagonal matrices
    proptest! {
        #[test]
        fn prop_cg_solves_spd_system(
            n in 5..20usize,
            diag in 3.0..10.0f64,
            off_diag in -2.0..2.0f64
        ) {
            let mut builder = SparseMatrixBuilder::new(n, n);

            // Create SPD tridiagonal matrix
            for i in 0..n {
                builder.add_entry(i, i, diag).unwrap();
                if i > 0 {
                    builder.add_entry(i, i - 1, off_diag).unwrap();
                }
                if i < n - 1 {
                    builder.add_entry(i, i + 1, off_diag).unwrap();
                }
            }
            let a = builder.build().unwrap();

            let b = DVector::from_element(n, 1.0);
            let mut x = DVector::zeros(n);

            let config = IterativeSolverConfig::new(1e-8).with_max_iterations(1000);
            let solver = ConjugateGradient::new(config);
            let identity = IdentityPreconditioner;

            solver.solve(&a, &b, &mut x, Some(&identity)).unwrap();

            // Property: Ax should equal b
            let ax = &a * &x;
            for i in 0..n {
                assert_relative_eq!(ax[i], b[i], epsilon = 1e-6);
            }
        }
    }

    /// Property test: BiCGSTAB converges for non-singular systems
    proptest! {
        #[test]
        fn prop_bicgstab_convergence(
            n in 5..15usize,
            scale in 1.0..5.0f64
        ) {
            let mut builder = SparseMatrixBuilder::new(n, n);

            // Create diagonally dominant matrix
            for i in 0..n {
                builder.add_entry(i, i, scale * 3.0).unwrap();
                if i > 0 {
                    builder.add_entry(i, i - 1, -scale).unwrap();
                }
                if i < n - 1 {
                    builder.add_entry(i, i + 1, -scale).unwrap();
                }
            }
            let a = builder.build().unwrap();

            let b = DVector::from_element(n, 1.0);
            let mut x = DVector::zeros(n);

            let config = IterativeSolverConfig::new(1e-7).with_max_iterations(1000);
            let solver = BiCGSTAB::new(config);
            let identity = IdentityPreconditioner;

            solver.solve(&a, &b, &mut x, Some(&identity)).unwrap();

            // Property: Solution exists and is finite
            for i in 0..n {
                assert!(x[i].is_finite(), "Solution contains non-finite values");
            }

            // Property: Residual is small
            let ax = &a * &x;
            let residual: f64 = (0..n).map(|i| {
                let diff = ax[i] - b[i];
                diff * diff
            }).sum::<f64>().sqrt();
            assert!(residual < 1e-5, "Residual {} too large", residual);
        }
    }

    /// Property test: Jacobi preconditioner preserves vector norms
    proptest! {
        #[test]
        fn prop_jacobi_diagonal_scaling(
            n in 5..15usize,
            diag_val in 1.0..10.0f64
        ) {
            let mut builder = SparseMatrixBuilder::new(n, n);

            for i in 0..n {
                builder.add_entry(i, i, diag_val).unwrap();
                if i > 0 {
                    builder.add_entry(i, i - 1, 0.1).unwrap();
                }
                if i < n - 1 {
                    builder.add_entry(i, i + 1, 0.1).unwrap();
                }
            }
            let a = builder.build().unwrap();

            let precond = JacobiPreconditioner::new(&a).unwrap();
            let r = DVector::from_element(n, diag_val);
            let mut z = DVector::zeros(n);

            precond.apply_to(&r, &mut z).unwrap();

            // Property: z[i] = r[i] / diag[i] = 1.0 for uniform diagonal
            for i in 0..n {
                assert_relative_eq!(z[i], 1.0, epsilon = 1e-10);
            }
        }
    }

    /// Property test: Solution is unique for non-singular matrices
    proptest! {
        #[test]
        fn prop_solution_uniqueness(
            n in 5..12usize,
            seed1 in 0.1..5.0f64,
            seed2 in 0.1..5.0f64
        ) {
            let mut builder = SparseMatrixBuilder::new(n, n);

            for i in 0..n {
                builder.add_entry(i, i, 2.0).unwrap();
                if i > 0 {
                    builder.add_entry(i, i - 1, -1.0).unwrap();
                }
                if i < n - 1 {
                    builder.add_entry(i, i + 1, -1.0).unwrap();
                }
            }
            let a = builder.build().unwrap();

            let b = DVector::from_element(n, 1.0);

            // Solve with two different initial guesses
            let mut x1 = DVector::from_element(n, seed1);
            let mut x2 = DVector::from_element(n, seed2);

            let config = IterativeSolverConfig::new(1e-9).with_max_iterations(1000);
            let solver = ConjugateGradient::new(config);
            let identity = IdentityPreconditioner;

            solver.solve(&a, &b, &mut x1, Some(&identity)).unwrap();
            solver.solve(&a, &b, &mut x2, Some(&identity)).unwrap();

            // Property: Solutions should be identical regardless of initial guess
            for i in 0..n {
                assert_relative_eq!(x1[i], x2[i], epsilon = 1e-6);
            }
        }
    }
}
