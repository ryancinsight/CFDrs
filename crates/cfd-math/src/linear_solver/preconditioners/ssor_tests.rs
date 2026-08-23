//! Value-semantic coverage for the successive over-relaxation preconditioner
//! the tiered solver chain applies.
//!
//! Athena owns the preconditioner seam (Atlas ADR 0033), and
//! [`athena_leto::SuccessiveOverRelaxation`] is the implementation
//! [`LinearSolverChain`](crate::linear_solver::LinearSolverChain) reaches for.
//! Its apply is one forward sweep of the lower splitting factor,
//!
//! ```text
//! z = (D/omega + L)^-1 r
//! ```
//!
//! which is the operator these tests pin: the reference vectors below are the
//! closed-form substitution over the tridiagonal operator, not values recorded
//! from a run.

#[cfg(test)]
mod tests {
    use athena_core::Preconditioner;
    use athena_leto::{LetoBackend, LetoBackendError, SuccessiveOverRelaxation};
    use cfd_core::error::Result;
    use leto::{Array1, Array2};
    use leto_ops::CsrMatrix;

    /// `tridiag(-1, 2, -1)` of the requested order: symmetric positive
    /// definite, the one-dimensional Laplacian stencil.
    fn tridiagonal_spd(size: usize) -> CsrMatrix<f64> {
        let mut values = vec![0.0; size * size];
        for row in 0..size {
            values[row * size + row] = 2.0;
            if row > 0 {
                values[row * size + row - 1] = -1.0;
            }
            if row + 1 < size {
                values[row * size + row + 1] = -1.0;
            }
        }

        let dense = Array2::from_shape_vec([size, size], values)
            .expect("tridiagonal construction has a valid shape");
        CsrMatrix::from_dense(&dense.view())
    }

    /// Apply the preconditioner over Athena's borrowed views.
    fn apply(
        preconditioner: &SuccessiveOverRelaxation<f64>,
        residual: &Array1<f64>,
        output: &mut Array1<f64>,
    ) -> std::result::Result<(), LetoBackendError> {
        let backend = LetoBackend::<f64>::default();
        preconditioner.apply(&backend, residual.view(), output.view_mut())
    }

    #[test]
    fn ssor_preserves_zero_and_produces_input_sensitive_output() -> Result<()> {
        let preconditioner = SuccessiveOverRelaxation::from_csr(&tridiagonal_spd(4), 1.0)
            .expect("relaxation 1.0 lies in the convergent range (0, 2)");
        let zero = Array1::zeros([4]);
        let mut zero_output = Array1::zeros([4]);
        apply(&preconditioner, &zero, &mut zero_output).expect("matching lengths");
        assert!(zero_output.iter().all(|value| *value == 0.0));

        // Forward substitution over (D + L) with r = 1:
        //   z0 = 1/2                     = 0.5
        //   z1 = (1 + z0)/2              = 0.75
        //   z2 = (1 + z1)/2              = 0.875
        //   z3 = (1 + z2)/2              = 0.9375
        let input = Array1::from_elem([4], 1.0);
        let mut output = Array1::zeros([4]);
        apply(&preconditioner, &input, &mut output).expect("matching lengths");
        assert_eq!(
            output.as_slice(),
            Some(&[0.5, 0.75, 0.875, 0.9375][..]),
            "Gauss-Seidel sweep must match the closed-form substitution"
        );
        Ok(())
    }

    #[test]
    fn ssor_rejects_mismatched_vector_lengths_at_the_seam() -> Result<()> {
        let preconditioner = SuccessiveOverRelaxation::from_csr(&tridiagonal_spd(4), 1.0)
            .expect("relaxation 1.0 lies in the convergent range (0, 2)");
        let residual = Array1::zeros([3]);
        let mut output = Array1::zeros([4]);
        let result = apply(&preconditioner, &residual, &mut output);
        assert!(
            matches!(
                result,
                Err(LetoBackendError::LengthMismatch { left, right })
                    if left == 3 && right == 4
            ),
            "expected a typed length mismatch, got {result:?}"
        );
        Ok(())
    }

    #[test]
    fn ssor_uses_relaxation_parameter_in_result() -> Result<()> {
        let matrix = tridiagonal_spd(4);
        let input = Array1::from_elem([4], 1.0);
        let unit = SuccessiveOverRelaxation::from_csr(&matrix, 1.0)
            .expect("relaxation 1.0 lies in the convergent range (0, 2)");
        let relaxed = SuccessiveOverRelaxation::from_csr(&matrix, 1.5)
            .expect("relaxation 1.5 lies in the convergent range (0, 2)");
        let mut unit_output = Array1::zeros([4]);
        let mut relaxed_output = Array1::zeros([4]);
        apply(&unit, &input, &mut unit_output).expect("matching lengths");
        apply(&relaxed, &input, &mut relaxed_output).expect("matching lengths");
        assert!(unit_output
            .iter()
            .zip(relaxed_output.iter())
            .any(|(unit, relaxed)| (unit - relaxed).abs() > 1.0e-12));

        // Forward substitution over (D/1.5 + L) with r = 1; the scaled pivot
        // is 2/1.5 = 4/3, so each step multiplies by 3/4:
        //   z0 = 1 * 3/4                 = 0.75
        //   z1 = (1 + z0) * 3/4          = 1.3125
        //   z2 = (1 + z1) * 3/4          = 1.734375
        //   z3 = (1 + z2) * 3/4          = 2.05078125
        assert_eq!(
            relaxed_output.as_slice(),
            Some(&[0.75, 1.3125, 1.734375, 2.050_781_25][..])
        );
        Ok(())
    }
}
