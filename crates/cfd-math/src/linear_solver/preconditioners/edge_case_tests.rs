//! Comprehensive edge case tests for preconditioners
//!
//! Tests cover positive/negative/zero/boundary cases per SRS requirements
//! Focus on ILU, SSOR, and IncompleteCholesky preconditioners

#[cfg(test)]
mod preconditioner_edge_tests {
    use crate::linear_solver::preconditioners::basic::IdentityPreconditioner;
    use crate::linear_solver::preconditioners::cholesky::IncompleteCholesky;
    use crate::linear_solver::preconditioners::deflation::DeflationPreconditioner;
    use crate::linear_solver::preconditioners::ilu::IncompleteLU;
    use crate::linear_solver::preconditioners::ssor::SSOR;
    use crate::linear_solver::Preconditioner;
    use crate::sparse::SparseMatrixBuilder;
    use approx::assert_relative_eq;
    use cfd_core::error::{Error, Result};
    use leto::Array1;
    use leto_ops::CsrMatrix;

    fn l2_norm(v: &Array1<f64>) -> f64 {
        v.iter().map(|x| x * x).sum::<f64>().sqrt()
    }

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
        builder.build()
    }

    fn assert_invalid_configuration(result: Result<()>, expected_message: &str) {
        match result {
            Err(Error::InvalidConfiguration(message)) => assert_eq!(message, expected_message),
            Err(error) => panic!("expected invalid configuration, got {error:?}"),
            Ok(()) => panic!("expected invalid configuration error"),
        }
    }

    /// Test ILU with zero vector input (boundary case)
    #[test]
    fn test_ilu_zero_vector() -> Result<()> {
        let a = create_tridiagonal_spd(10)?;
        let ilu = IncompleteLU::new(&a)?;

        let r = Array1::zeros([10]);
        let mut z = Array1::zeros([10]);
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

        let r = Array1::from_elem([10], -1.0);
        let mut z = Array1::zeros([10]);
        ilu.apply_to(&r, &mut z)?;

        // Solution should be finite and consistent
        assert!(z.iter().all(|&zi: &f64| zi.is_finite()));
        assert!(l2_norm(&z) > 0.0);
        Ok(())
    }

    /// Test SSOR with boundary omega values
    #[test]
    fn test_ssor_boundary_omega() -> Result<()> {
        let a = create_tridiagonal_spd(10)?;
        let r = Array1::from_elem([10], 1.0);

        // Test omega near boundaries: (0, 2) range
        for omega in [0.1, 0.5, 1.0, 1.5, 1.9] {
            let ssor = SSOR::with_omega(a.clone(), omega)?;
            let mut z = Array1::zeros([10]);
            ssor.apply_to(&r, &mut z)?;

            // Verify convergence properties
            assert!(
                z.iter().all(|&zi: &f64| zi.is_finite()),
                "omega={omega}: NaN detected"
            );
            assert!(l2_norm(&z) > 0.0, "omega={omega}: Zero output");
        }
        Ok(())
    }

    /// Test SSOR with zero input vector
    #[test]
    fn test_ssor_zero_vector() -> Result<()> {
        let a = create_tridiagonal_spd(10)?;
        let ssor = SSOR::new(a.clone())?;

        let r = Array1::zeros([10]);
        let mut z = Array1::zeros([10]);
        ssor.apply_to(&r, &mut z)?;

        // SSOR(M)*0 = 0
        for i in 0..10 {
            assert_relative_eq!(z[i], 0.0, epsilon = 1e-14);
        }
        Ok(())
    }

    /// Test SSOR rejects mismatched Leto vector lengths.
    #[test]
    fn test_ssor_rejects_mismatched_vector_lengths() -> Result<()> {
        let a = create_tridiagonal_spd(10)?;
        let ssor = SSOR::new(a.clone())?;

        let r = Array1::zeros([9]);
        let mut z = Array1::zeros([10]);
        assert_invalid_configuration(
            ssor.apply_to(&r, &mut z),
            "SSOR residual length mismatch: expected 10, got 9",
        );

        let r = Array1::zeros([10]);
        let mut z = Array1::zeros([9]);
        assert_invalid_configuration(
            ssor.apply_to(&r, &mut z),
            "SSOR output length mismatch: expected 10, got 9",
        );
        Ok(())
    }

    /// Test Incomplete Cholesky with positive definite matrix
    #[test]
    fn test_cholesky_spd_matrix() -> Result<()> {
        let a = create_tridiagonal_spd(10)?;
        let chol = IncompleteCholesky::new(&a)?;

        let r = Array1::from_elem([10], 1.0);
        let mut z = Array1::zeros([10]);
        chol.apply_to(&r, &mut z)?;

        // Cholesky should produce valid factorization
        assert!(z.iter().all(|&zi: &f64| zi.is_finite()));
        assert!(l2_norm(&z) > 0.0);
        Ok(())
    }

    /// Test Incomplete Cholesky with zero vector
    #[test]
    fn test_cholesky_zero_vector() -> Result<()> {
        let a = create_tridiagonal_spd(10)?;
        let chol = IncompleteCholesky::new(&a)?;

        let r = Array1::zeros([10]);
        let mut z = Array1::zeros([10]);
        chol.apply_to(&r, &mut z)?;

        // Cholesky(L)*0 = 0
        for i in 0..10 {
            assert_relative_eq!(z[i], 0.0, epsilon = 1e-14);
        }
        Ok(())
    }

    /// Test Incomplete Cholesky rejects mismatched Leto vector lengths.
    #[test]
    fn test_cholesky_rejects_mismatched_vector_lengths() -> Result<()> {
        let a = create_tridiagonal_spd(10)?;
        let chol = IncompleteCholesky::new(&a)?;

        let r = Array1::zeros([9]);
        let mut z = Array1::zeros([10]);
        assert!(chol.apply_to(&r, &mut z).is_err());

        let r = Array1::zeros([10]);
        let mut z = Array1::zeros([9]);
        assert!(chol.apply_to(&r, &mut z).is_err());
        Ok(())
    }

    /// Test deflation projection directly over Leto vector state.
    #[test]
    fn test_deflation_projects_known_eigenvector() -> Result<()> {
        let eigenvector =
            Array1::from_shape_vec([3], vec![1.0, 0.0, 0.0]).expect("valid eigenvector");
        let preconditioner = DeflationPreconditioner::new(
            Box::new(IdentityPreconditioner),
            vec![eigenvector],
            vec![1.0],
        )?;

        let r = Array1::from_shape_vec([3], vec![2.0, 3.0, 4.0]).expect("valid residual");
        let mut z = Array1::zeros([3]);
        preconditioner.apply_to(&r, &mut z)?;

        assert_relative_eq!(z[0], 0.0, epsilon = 1e-14);
        assert_relative_eq!(z[1], 3.0, epsilon = 1e-14);
        assert_relative_eq!(z[2], 4.0, epsilon = 1e-14);
        Ok(())
    }

    /// Test deflation rejects mismatched Leto vector lengths.
    #[test]
    fn test_deflation_rejects_mismatched_vector_lengths() -> Result<()> {
        let eigenvector =
            Array1::from_shape_vec([3], vec![1.0, 0.0, 0.0]).expect("valid eigenvector");
        let preconditioner = DeflationPreconditioner::new(
            Box::new(IdentityPreconditioner),
            vec![eigenvector],
            vec![1.0],
        )?;

        let r = Array1::zeros([3]);
        let mut z = Array1::zeros([2]);
        assert_invalid_configuration(
            preconditioner.apply_to(&r, &mut z),
            "deflation output length mismatch: expected 3, got 2",
        );

        let eigenvector =
            Array1::from_shape_vec([3], vec![1.0, 0.0, 0.0]).expect("valid eigenvector");
        let mut preconditioner = DeflationPreconditioner::new(
            Box::new(IdentityPreconditioner),
            vec![eigenvector],
            vec![1.0],
        )?;
        let bad_eigenvector =
            Array1::from_shape_vec([2], vec![1.0, 0.0]).expect("valid eigenvector");
        assert_invalid_configuration(
            preconditioner.add_eigenpair(bad_eigenvector, 1.0),
            "Eigenvector dimension mismatch: expected 3, got 2",
        );
        Ok(())
    }

    /// Test deflation rejects zero eigenvalues before division.
    #[test]
    fn test_deflation_rejects_zero_eigenvalue() {
        let eigenvector =
            Array1::from_shape_vec([3], vec![1.0, 0.0, 0.0]).expect("valid eigenvector");
        let result = DeflationPreconditioner::new(
            Box::new(IdentityPreconditioner),
            vec![eigenvector],
            vec![0.0],
        );

        match result {
            Err(Error::InvalidConfiguration(message)) => assert_eq!(
                message,
                "Eigenvalue 0 is zero; deflation requires nonzero eigenvalues"
            ),
            Err(error) => panic!("expected invalid configuration, got {error:?}"),
            Ok(_) => panic!("expected invalid configuration error"),
        }
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
        let r = Array1::from_elem([n], 1.0);
        let mut z = Array1::zeros([n]);
        ilu.apply_to(&r, &mut z)?;

        // Verify ILU handles 2D pattern
        assert!(z.iter().all(|&zi: &f64| zi.is_finite()));
        assert!(l2_norm(&z) > 0.0);
        Ok(())
    }

    /// Test SSOR with mixed positive/negative residual
    #[test]
    fn test_ssor_mixed_residual() -> Result<()> {
        let a = create_tridiagonal_spd(10)?;
        let ssor = SSOR::new(a.clone())?;

        // Create mixed sign residual
        let mut r = Array1::zeros([10]);
        for i in 0..10 {
            r[i] = if i % 2 == 0 { 1.0 } else { -1.0 };
        }

        let mut z = Array1::zeros([10]);
        ssor.apply_to(&r, &mut z)?;

        // SSOR should handle mixed signs
        assert!(z.iter().all(|&zi: &f64| zi.is_finite()));
        Ok(())
    }

    /// Test ILU with ill-conditioned matrix (large condition number)
    /// Edge case: κ(A) ≈ 10^6 tests numerical stability
    #[test]
    fn test_ilu_ill_conditioned() -> Result<()> {
        let n = 10;
        let mut builder = SparseMatrixBuilder::new(n, n);

        // Create ill-conditioned tridiagonal: large diagonal, small off-diagonal
        let diag_large = 1e6;
        let offdiag_small = 1.0;

        for i in 0..n {
            builder.add_entry(i, i, diag_large)?;
            if i > 0 {
                builder.add_entry(i, i - 1, -offdiag_small)?;
            }
            if i < n - 1 {
                builder.add_entry(i, i + 1, -offdiag_small)?;
            }
        }
        let a = builder.build()?;

        let ilu = IncompleteLU::new(&a)?;
        let r = Array1::from_elem([n], 1.0);
        let mut z = Array1::zeros([n]);
        ilu.apply_to(&r, &mut z)?;

        // Verify numerical stability: no NaN/Inf
        assert!(z.iter().all(|&zi: &f64| zi.is_finite()));
        assert!(l2_norm(&z) > 0.0);

        // Preconditioner should scale residual
        assert!(l2_norm(&z) < l2_norm(&r) * 1e3); // Scaled by ~1/diag
        Ok(())
    }

    /// Test IncompleteCholesky with nearly diagonal matrix
    /// Edge case: minimal off-diagonal coupling tests preconditioner effectiveness
    #[test]
    fn test_cholesky_nearly_diagonal() -> Result<()> {
        let n = 8;
        let mut builder = SparseMatrixBuilder::new(n, n);

        // Nearly diagonal: large diagonal, tiny off-diagonal
        let diag = 10.0;
        let offdiag = 1e-3;

        for i in 0..n {
            builder.add_entry(i, i, diag)?;
            if i > 0 {
                builder.add_entry(i, i - 1, offdiag)?;
            }
            if i < n - 1 {
                builder.add_entry(i, i + 1, offdiag)?;
            }
        }
        let a = builder.build()?;

        let cholesky = IncompleteCholesky::new(&a)?;
        let r = Array1::from_elem([n], 5.0);
        let mut z = Array1::zeros([n]);
        cholesky.apply_to(&r, &mut z)?;

        // Nearly diagonal → preconditioner ≈ diagonal inverse
        for i in 0..n {
            assert_relative_eq!(z[i], r[i] / diag, epsilon = 1e-2);
        }
        Ok(())
    }

    /// Test ILU with single element matrix (boundary case: n=1)
    /// Edge case: minimal system tests implementation correctness
    #[test]
    fn test_ilu_single_element() -> Result<()> {
        let n = 1;
        let mut builder = SparseMatrixBuilder::new(n, n);
        builder.add_entry(0, 0, 5.0)?;
        let a = builder.build()?;

        let ilu = IncompleteLU::new(&a)?;
        let r = Array1::from_elem([n], 10.0);
        let mut z = Array1::zeros([n]);
        ilu.apply_to(&r, &mut z)?;

        // Single element: ILU(M) * r = A^{-1} * r = r / a[0,0]
        assert_relative_eq!(z[0], 10.0 / 5.0, epsilon = 1e-14);
        Ok(())
    }

    /// Test SSOR with strongly diagonally dominant matrix
    /// Edge case: |a_ii| >> sum|a_ij| tests convergence behavior
    #[test]
    fn test_ssor_strongly_diagonal_dominant() -> Result<()> {
        let n = 10;
        let mut builder = SparseMatrixBuilder::new(n, n);

        // Strong diagonal dominance: diag = 100, off-diag = 1
        for i in 0..n {
            builder.add_entry(i, i, 100.0)?;
            if i > 0 {
                builder.add_entry(i, i - 1, -1.0)?;
            }
            if i < n - 1 {
                builder.add_entry(i, i + 1, -1.0)?;
            }
        }
        let a = builder.build()?;

        let ssor = SSOR::new(a.clone())?;
        let r = Array1::from_elem([n], 1.0);
        let mut z = Array1::zeros([n]);
        ssor.apply_to(&r, &mut z)?;

        // Strong diagonal dominance → preconditioner ≈ diagonal inverse
        for i in 0..n {
            assert_relative_eq!(z[i], r[i] / 100.0, epsilon = 1e-1);
        }
        Ok(())
    }
}
