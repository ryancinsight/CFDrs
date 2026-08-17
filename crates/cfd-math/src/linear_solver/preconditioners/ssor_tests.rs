//! Value-semantic coverage for the Leto SSOR preconditioner.

#[cfg(test)]
mod tests {
    use cfd_core::error::Result;
    use leto::{Array1, Array2, LetoError};
    use leto_ops::{CsrMatrix, Preconditioner, SSORPreconditioner};

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

    #[test]
    fn ssor_preserves_zero_and_produces_input_sensitive_output() -> Result<()> {
        let preconditioner = SSORPreconditioner::new(tridiagonal_spd(4))?;
        let zero = Array1::zeros([4]);
        let mut zero_output = Array1::zeros([4]);
        preconditioner.apply_to(&zero, &mut zero_output)?;
        assert!(zero_output.iter().all(|value| *value == 0.0));

        let input = Array1::from_elem([4], 1.0);
        let mut output = Array1::zeros([4]);
        preconditioner.apply_to(&input, &mut output)?;
        assert!(output.iter().all(|value| value.is_finite()));
        assert!(output.iter().any(|value| *value != 0.0));
        Ok(())
    }

    #[test]
    fn ssor_rejects_mismatched_vector_lengths_at_provider_boundary() -> Result<()> {
        let preconditioner = SSORPreconditioner::new(tridiagonal_spd(4))?;
        let residual = Array1::zeros([3]);
        let mut output = Array1::zeros([4]);
        let result = preconditioner.apply_to(&residual, &mut output);
        assert!(matches!(
            result,
            Err(LetoError::InvalidInput(message))
                if message == "SSOR residual length mismatch: expected 4, got 3"
        ));
        Ok(())
    }

    #[test]
    fn ssor_uses_relaxation_parameter_in_result() -> Result<()> {
        let matrix = tridiagonal_spd(4);
        let input = Array1::from_elem([4], 1.0);
        let unit = SSORPreconditioner::with_omega(matrix.clone(), 1.0)?;
        let relaxed = SSORPreconditioner::with_omega(matrix, 1.5)?;
        let mut unit_output = Array1::zeros([4]);
        let mut relaxed_output = Array1::zeros([4]);
        unit.apply_to(&input, &mut unit_output)?;
        relaxed.apply_to(&input, &mut relaxed_output)?;
        assert!(unit_output
            .iter()
            .zip(relaxed_output.iter())
            .any(|(unit, relaxed)| (unit - relaxed).abs() > 1.0e-12));
        assert_eq!(
            unit_output.as_slice(),
            Some(&[1.2109375, 1.421875, 1.34375, 0.9375][..])
        );
        Ok(())
    }
}
