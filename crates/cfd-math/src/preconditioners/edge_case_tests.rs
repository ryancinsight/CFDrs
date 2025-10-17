//! Comprehensive edge case tests for preconditioners
//!
//! Tests cover positive/negative/zero/boundary cases per SRS requirements
//! Focus on ILU, SSOR, and IncompleteCholesky preconditioners

#[cfg(test)]
mod preconditioner_edge_tests {
    use crate::linear_solver::Preconditioner;
    use crate::preconditioners::cholesky::IncompleteCholesky;
    use crate::preconditioners::ilu::IncompleteLU;
    use crate::preconditioners::ssor::SSOR;
    use crate::sparse::SparseMatrixBuilder;
    use approx::assert_relative_eq;
    use cfd_core::error::Result;
    use nalgebra::DVector;
    use nalgebra_sparse::CsrMatrix;

    /// Helper: Create tridiagonal SPD matrix
    fn create_tridiagonal_spd(n: usize) -> Result<CsrMatrix<f64>> {
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

    /// Test ILU with zero vector input (boundary case)
    #[test]
    fn test_ilu_zero_vector() -> Result<()> {
        let a = create_tridiagonal_spd(10)?;
        let ilu = IncompleteLU::new(&a)?;
        
        let r = DVector::zeros(10);
        let mut z = DVector::zeros(10);
        ilu.apply_to(&r, &mut z)?;
        
        // ILU(M)*0 = 0
        for i in 0..10 {
            assert_relative_eq!(z[i], 0.0, epsilon = 1e-14);
        }
        Ok(())
    }

    /// Test ILU with negative values (numerical robustness)
    #[test]
    fn test_ilu_negative_values() -> Result<()> {
        let a = create_tridiagonal_spd(10)?;
        let ilu = IncompleteLU::new(&a)?;
        
        let r = DVector::from_element(10, -1.0);
        let mut z = DVector::zeros(10);
        ilu.apply_to(&r, &mut z)?;
        
        // Solution should be finite and consistent
        assert!(z.iter().all(|&zi: &f64| zi.is_finite()));
        assert!(z.norm() > 0.0);
        Ok(())
    }

    /// Test SSOR with boundary omega values
    #[test]
    fn test_ssor_boundary_omega() -> Result<()> {
        let a = create_tridiagonal_spd(10)?;
        let r = DVector::from_element(10, 1.0);
        
        // Test omega near boundaries: (0, 2) range
        for omega in [0.1, 0.5, 1.0, 1.5, 1.9] {
            let ssor = SSOR::with_omega(a.clone(), omega)?;
            let mut z = DVector::zeros(10);
            ssor.apply_to(&r, &mut z)?;
            
            // Verify convergence properties
            assert!(z.iter().all(|&zi: &f64| zi.is_finite()), "omega={}: NaN detected", omega);
            assert!(z.norm() > 0.0, "omega={}: Zero output", omega);
        }
        Ok(())
    }

    /// Test SSOR with zero input vector
    #[test]
    fn test_ssor_zero_vector() -> Result<()> {
        let a = create_tridiagonal_spd(10)?;
        let ssor = SSOR::new(a)?;
        
        let r = DVector::zeros(10);
        let mut z = DVector::zeros(10);
        ssor.apply_to(&r, &mut z)?;
        
        // SSOR(M)*0 = 0
        for i in 0..10 {
            assert_relative_eq!(z[i], 0.0, epsilon = 1e-14);
        }
        Ok(())
    }

    /// Test Incomplete Cholesky with positive definite matrix
    #[test]
    fn test_cholesky_spd_matrix() -> Result<()> {
        let a = create_tridiagonal_spd(10)?;
        let chol = IncompleteCholesky::new(&a)?;
        
        let r = DVector::from_element(10, 1.0);
        let mut z = DVector::zeros(10);
        chol.apply_to(&r, &mut z)?;
        
        // Cholesky should produce valid factorization
        assert!(z.iter().all(|&zi: &f64| zi.is_finite()));
        assert!(z.norm() > 0.0);
        Ok(())
    }

    /// Test Incomplete Cholesky with zero vector
    #[test]
    fn test_cholesky_zero_vector() -> Result<()> {
        let a = create_tridiagonal_spd(10)?;
        let chol = IncompleteCholesky::new(&a)?;
        
        let r = DVector::zeros(10);
        let mut z = DVector::zeros(10);
        chol.apply_to(&r, &mut z)?;
        
        // Cholesky(L)*0 = 0
        for i in 0..10 {
            assert_relative_eq!(z[i], 0.0, epsilon = 1e-14);
        }
        Ok(())
    }

    /// Test ILU with pentadiagonal matrix (2D Laplacian pattern)
    #[test]
    fn test_ilu_pentadiagonal() -> Result<()> {
        let n = 25; // 5x5 grid
        let mut builder = SparseMatrixBuilder::new(n, n);
        let nx = 5;
        
        // 2D Laplacian (5-point stencil)
        for i in 0..n {
            let row = i / nx;
            let col = i % nx;
            
            builder.add_entry(i, i, 4.0)?;
            if col > 0 {
                builder.add_entry(i, i - 1, -1.0)?;
            }
            if col < nx - 1 {
                builder.add_entry(i, i + 1, -1.0)?;
            }
            if row > 0 {
                builder.add_entry(i, i - nx, -1.0)?;
            }
            if row < nx - 1 {
                builder.add_entry(i, i + nx, -1.0)?;
            }
        }
        let a = builder.build()?;
        
        let ilu = IncompleteLU::new(&a)?;
        let r = DVector::from_element(n, 1.0);
        let mut z = DVector::zeros(n);
        ilu.apply_to(&r, &mut z)?;
        
        // Verify ILU handles 2D pattern
        assert!(z.iter().all(|&zi: &f64| zi.is_finite()));
        assert!(z.norm() > 0.0);
        Ok(())
    }

    /// Test SSOR with mixed positive/negative residual
    #[test]
    fn test_ssor_mixed_residual() -> Result<()> {
        let a = create_tridiagonal_spd(10)?;
        let ssor = SSOR::new(a)?;
        
        // Create mixed sign residual
        let mut r = DVector::zeros(10);
        for i in 0..10 {
            r[i] = if i % 2 == 0 { 1.0 } else { -1.0 };
        }
        
        let mut z = DVector::zeros(10);
        ssor.apply_to(&r, &mut z)?;
        
        // SSOR should handle mixed signs
        assert!(z.iter().all(|&zi: &f64| zi.is_finite()));
        Ok(())
    }
}
