//! Tests for linear solvers

#[cfg(test)]
mod tests {
    use super::super::*;
    use approx::assert_relative_eq;
    use nalgebra::{DMatrix, DVector};
    use nalgebra_sparse::{CooMatrix, CsrMatrix};

    fn create_validation_matrix() -> CsrMatrix<f64> {
        // Create a simple 3x3 SPD matrix
        let dense = DMatrix::from_row_slice(3, 3, &[
            4.0, -1.0, 0.0,
            -1.0, 4.0, -1.0,
            0.0, -1.0, 4.0,
        ]);
        CsrMatrix::from(&dense)
    }

    fn create_1d_poisson_matrix(n: usize) -> CsrMatrix<f64> {
        // Create tridiagonal matrix for 1D Poisson: [-1, 2, -1]
        let mut builder = crate::sparse::SparseMatrixBuilder::new(n, n);
        
        for i in 0..n {
            // Diagonal
            builder.add_entry(i, i, 2.0).expect("Failed to add diagonal entry");
            
            // Off-diagonals
            if i > 0 {
                builder.add_entry(i, i-1, -1.0).expect("Failed to add lower diagonal");
            }
            if i < n-1 {
                builder.add_entry(i, i+1, -1.0).expect("Failed to add upper diagonal");
            }
        }
        
        builder.build().expect("Failed to build matrix")
    }

    fn create_test_matrix() -> CsrMatrix<f64> {
        let mut coo = CooMatrix::new(3, 3);
        coo.push(0, 0, 4.0);
        coo.push(0, 1, -1.0);
        coo.push(1, 0, -1.0);
        coo.push(1, 1, 4.0);
        coo.push(1, 2, -1.0);
        coo.push(2, 1, -1.0);
        coo.push(2, 2, 4.0);
        CsrMatrix::from(&coo)
    }

    #[test]
    fn test_cg_solver() {
        let a = create_test_matrix();
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        
        let solver = ConjugateGradient::default();
        let x = solver.solve(&a, &b, None).expect("Failed to solve");
        
        let residual = &b - &a * &x;
        assert!(residual.norm() < 1e-10);
    }

    #[test]
    fn test_bicgstab_solver() {
        let a = create_test_matrix();
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        
        let solver = BiCGSTAB::default();
        let x = solver.solve(&a, &b, None).expect("Failed to solve");
        
        let residual = &b - &a * &x;
        assert!(residual.norm() < 1e-10);
    }

    #[test]
    fn test_jacobi_preconditioner() {
        let a = create_test_matrix();
        let r = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        
        let precond = JacobiPreconditioner::new(&a).expect("Failed to create preconditioner");
        let mut z = DVector::zeros(3);
        precond.apply_to(&r, &mut z).expect("Failed to apply preconditioner");
        
        // Check that z is correctly scaled by inverse diagonal
        assert_relative_eq!(z[0], r[0] / 4.0, epsilon = 1e-12);
        assert_relative_eq!(z[1], r[1] / 4.0, epsilon = 1e-12);
        assert_relative_eq!(z[2], r[2] / 4.0, epsilon = 1e-12);
    }

    #[test]
    fn test_sor_omega_validation() {
        let a = create_test_matrix();
        
        // Test invalid omega values
        assert!(SORPreconditioner::new(&a, 0.0).is_err());
        assert!(SORPreconditioner::new(&a, 2.0).is_err());
        assert!(SORPreconditioner::new(&a, -0.5).is_err());
        
        // Test valid omega
        assert!(SORPreconditioner::new(&a, 1.5).is_ok());
    }

    #[test]
    fn test_1d_poisson_omega_optimization() {
        let a = create_1d_poisson_matrix(10);
        
        // Should successfully create SOR with optimized omega
        let sor = SORPreconditioner::with_omega_for_1d_poisson(&a).expect("Failed to create SOR");
        
        // Omega should be in valid range
        assert!(sor.omega > 0.0 && sor.omega < 2.0);
        
        // For n=10, calculate expected omega
        let n = 10.0_f64;
        let expected_omega = 2.0 / (1.0 + (std::f64::consts::PI / n).sin());
        
        // For n=10, optimal omega ≈ 1.69 (allow larger tolerance)
        assert!((sor.omega - expected_omega).abs() < 0.01);
    }

    #[test]
    fn test_1d_poisson_validation_rejects_invalid_matrix() {
        // Create a non-tridiagonal matrix
        let mut builder = crate::sparse::SparseMatrixBuilder::new(3, 3);
        builder.add_entry(0, 0, 2.0).expect("Failed to add entry");
        builder.add_entry(0, 2, 1.0).expect("Failed to add entry"); // Non-adjacent entry
        builder.add_entry(1, 1, 2.0).expect("Failed to add entry");
        builder.add_entry(2, 2, 2.0).expect("Failed to add entry");
        let non_tridiag = builder.build().expect("Failed to build matrix");
        
        // Should reject non-tridiagonal matrix
        assert!(SORPreconditioner::with_omega_for_1d_poisson(&non_tridiag).is_err());
    }

    #[test]
    fn test_preconditioned_cg() {
        let a = create_test_matrix();
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        
        let solver = ConjugateGradient::default();
        let precond = JacobiPreconditioner::new(&a).expect("Failed to create preconditioner");
        
        let x = solver.solve_preconditioned(&a, &b, &precond, None).expect("Failed to solve");
        
        let residual = &b - &a * &x;
        assert!(residual.norm() < 1e-10);
    }

    #[test]
    fn test_identity_preconditioner() {
        let r = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let precond = IdentityPreconditioner;
        
        let mut z = DVector::zeros(3);
        precond.apply_to(&r, &mut z).expect("Failed to apply preconditioner");
        
        assert_relative_eq!(z, r, epsilon = 1e-12);
    }
}