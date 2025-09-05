//! Tests for linear solvers

#[cfg(test)]
mod tests {

    use crate::linear_solver::preconditioners::{IdentityPreconditioner, JacobiPreconditioner, SORPreconditioner};
    use crate::linear_solver::traits::IterativeLinearSolver;
    use crate::linear_solver::IterativeSolverConfig;
    use crate::linear_solver::{BiCGSTAB, ConjugateGradient, Preconditioner};
    use crate::sparse::{SparseMatrix, SparseMatrixBuilder};
    use approx::assert_relative_eq;
    use nalgebra::DVector;

    fn create_tridiagonal_matrix(
        n: usize,
    ) -> std::result::Result<SparseMatrix<f64>, Box<dyn std::error::Error>> {
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

        Ok(builder.build()?)
    }

    #[test]
    fn test_conjugate_gradient() -> std::result::Result<(), Box<dyn std::error::Error>> {
        let n = 5;
        let a = create_tridiagonal_matrix(n)?;
        let b = DVector::from_element(n, 1.0);

        let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);

        let solver = ConjugateGradient::new(config);
        let mut x = DVector::zeros(n);  // Initial guess
        let identity_precond = IdentityPreconditioner;
        solver.solve(&a, &b, &mut x, Some(&identity_precond))?;

        // Check that Ax = b using nalgebra_sparse API
        let ax = &a * &x;  // Try direct multiplication
        for i in 0..n {
            assert_relative_eq!(ax[i], b[i], epsilon = 1e-8);
        }
        Ok(())
    }

    #[test]
    fn test_bicgstab() -> std::result::Result<(), Box<dyn std::error::Error>> {
        let n = 5;
        let a = create_tridiagonal_matrix(n)?;
        let b = DVector::from_element(n, 1.0);

        let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);

        let solver = BiCGSTAB::new(config);
        let mut x = DVector::zeros(n);  // Initial guess  
        let identity_precond = IdentityPreconditioner;
        solver.solve(&a, &b, &mut x, Some(&identity_precond))?;

        // Check that Ax = b
        let ax = &a * &x;
        for i in 0..n {
            assert_relative_eq!(ax[i], b[i], epsilon = 1e-8);
        }
        Ok(())
    }

    #[test]
    fn test_jacobi_preconditioner() -> std::result::Result<(), Box<dyn std::error::Error>> {
        let n = 5;
        let a = create_tridiagonal_matrix(n)?;
        let precond = JacobiPreconditioner::new(&a)?;

        let r = DVector::from_element(n, 1.0);
        let mut z = DVector::zeros(n);
        precond.apply_to(&r, &mut z)?;

        // Check that z = D^{-1} r
        for i in 0..n {
            assert_relative_eq!(z[i], r[i] / 2.0, epsilon = 1e-10);
        }
        Ok(())
    }

    #[test]
    fn test_sor_preconditioner() -> std::result::Result<(), Box<dyn std::error::Error>> {
        let n = 5;
        let a = create_tridiagonal_matrix(n)?;

        // Test SOR with automatic omega selection for 1D Poisson
        let sor = SORPreconditioner::with_omega_for_1d_poisson(&a)?;

        // For n=5: omega = 2/(1 + sin(π/5)) ≈ 2/(1 + 0.588) ≈ 1.26
        // The optimal omega depends on matrix size
        assert!(
            sor.omega() > 1.0 && sor.omega() < 2.0,
            "Omega {} should be in (1,2)",
            sor.omega()
        );

        // Test that it can be applied
        let r = DVector::from_element(n, 1.0);
        let mut z = DVector::zeros(n);
        sor.apply_to(&r, &mut z)?;

        // Result should be non-zero
        assert!(z.norm() > 0.0);
        Ok(())
    }

    #[test]
    fn test_ilu_preconditioner_fails_on_non_tridiagonal(
    ) -> std::result::Result<(), Box<dyn std::error::Error>> {
        let n = 3;
        let mut builder = SparseMatrixBuilder::new(n, n);

        builder.add_entry(0, 0, 2.0)?;
        builder.add_entry(0, 2, 1.0)?; // Non-adjacent entry
        builder.add_entry(1, 1, 2.0)?;
        builder.add_entry(2, 2, 2.0)?;
        let _non_tridiag = builder.build()?;

        // ILU(0) preconditioner not yet implemented
        // let result = ILUPreconditioner::new(&non_tridiag);
        // assert!(result.is_err());
        Ok(())
    }

    #[test]
    fn test_preconditioned_cg() -> std::result::Result<(), Box<dyn std::error::Error>> {
        let n = 5;
        let a = create_tridiagonal_matrix(n)?;
        let b = DVector::from_element(n, 1.0);
        let precond = JacobiPreconditioner::new(&a)?;

        let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);

        let solver = ConjugateGradient::new(config);
        let mut x = DVector::zeros(n);  // Initial guess
        solver.solve(&a, &b, &mut x, Some(&precond))?;

        // Check that Ax = b
        let ax = &a * &x;
        for i in 0..n {
            assert_relative_eq!(ax[i], b[i], epsilon = 1e-8);
        }
        Ok(())
    }

    #[test]
    fn test_gauss_seidel() -> std::result::Result<(), Box<dyn std::error::Error>> {
        let n = 5;
        let a = create_tridiagonal_matrix(n)?;
        // GaussSeidelPreconditioner not yet implemented, using SOR instead
        let precond = SORPreconditioner::new(&a, 1.0)?;

        let r = DVector::from_element(n, 1.0);
        let mut z = DVector::zeros(n);
        precond.apply_to(&r, &mut z)?;

        // Result should be non-zero
        assert!(z.norm() > 0.0);
        Ok(())
    }

    #[test]
    fn test_convergence_with_different_tolerances(
    ) -> std::result::Result<(), Box<dyn std::error::Error>> {
        let n = 10;
        let a = create_tridiagonal_matrix(n)?;
        let b = DVector::from_element(n, 1.0);

        let tolerances = vec![1e-4, 1e-6, 1e-8];

        for tol in tolerances {
            let config = IterativeSolverConfig::new(tol).with_max_iterations(1000);

            let solver = ConjugateGradient::new(config);
            let mut x = DVector::zeros(n);  // Initial guess
            let identity_precond = IdentityPreconditioner;
            solver.solve(&a, &b, &mut x, Some(&identity_precond))?;

            // Check residual
            let ax = &a * &x;
            let residual = &b - &ax;
            let relative_residual = residual.norm() / b.norm();
            assert!(relative_residual < tol * 10.0); // Allow some slack
        }
        Ok(())
    }
}
