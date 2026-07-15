//! Checkpoint validation utilities

use crate::checkpoint::Checkpoint;
use crate::leto_arrays::all_row_major;
use eunomia::{FloatElement, RealField};
use leto::Array2;

/// Numerical constants for finite difference calculations
mod numerical_constants {
    /// Central difference denominator factor (`h_forward` - `h_backward` = 2h)
    pub const CENTRAL_DIFFERENCE_FACTOR: f64 = 2.0;
}

/// Checkpoint validator for ensuring data integrity
pub struct CheckpointValidator;

impl CheckpointValidator {
    /// Validate checkpoint for physical consistency
    pub fn validate_physics<T: RealField>(checkpoint: &Checkpoint<T>) -> Result<(), String> {
        // Check basic data consistency
        checkpoint.validate()?;

        // Check for NaN or infinite values
        if !Self::is_field_finite(&checkpoint.u_velocity) {
            return Err("U velocity contains non-finite values".to_string());
        }

        if !Self::is_field_finite(&checkpoint.v_velocity) {
            return Err("V velocity contains non-finite values".to_string());
        }

        if !Self::is_field_finite(&checkpoint.pressure) {
            return Err("Pressure contains non-finite values".to_string());
        }

        // Check optional fields
        if let Some(ref temp) = checkpoint.temperature {
            if !Self::is_field_finite(temp) {
                return Err("Temperature contains non-finite values".to_string());
            }

            // Temperature should be positive (in Kelvin)
            if !all_row_major(temp, |value| value > T::ZERO) {
                return Err("Temperature contains non-positive values".to_string());
            }
        }

        if let Some(ref k) = checkpoint.turbulence_k {
            if !Self::is_field_finite(k) {
                return Err("Turbulence k contains non-finite values".to_string());
            }

            // Turbulent kinetic energy should be non-negative
            if !all_row_major(k, |value| value >= T::ZERO) {
                return Err("Turbulence k contains negative values".to_string());
            }
        }

        if let Some(ref eps) = checkpoint.turbulence_epsilon {
            if !Self::is_field_finite(eps) {
                return Err("Turbulence epsilon contains non-finite values".to_string());
            }

            // Dissipation rate should be positive
            if !all_row_major(eps, |value| value > T::ZERO) {
                return Err("Turbulence epsilon contains non-positive values".to_string());
            }
        }

        Ok(())
    }

    /// Check if all values in a matrix are finite
    fn is_field_finite<T: RealField>(field: &Array2<T>) -> bool {
        all_row_major(field, |value| value.is_finite())
    }

    /// Check mass conservation (divergence-free condition)
    pub fn check_mass_conservation<T: RealField>(checkpoint: &Checkpoint<T>, tolerance: T) -> bool {
        let (ny, nx) = checkpoint.dimensions();
        let (domain_x, domain_y) = checkpoint.metadata.domain_size;

        let Some(nx_f64) = usize_to_exact_f64(nx) else {
            return false;
        };
        let Some(ny_f64) = usize_to_exact_f64(ny) else {
            return false;
        };

        let dx = <T as FloatElement>::from_f64(domain_x / nx_f64);
        let dy = <T as FloatElement>::from_f64(domain_y / ny_f64);
        let central_difference =
            <T as FloatElement>::from_f64(numerical_constants::CENTRAL_DIFFERENCE_FACTOR);

        let mut max_divergence = T::ZERO;

        // Check interior points
        for i in 1..(ny.saturating_sub(1)) {
            for j in 1..(nx.saturating_sub(1)) {
                let dudx = (*checkpoint
                    .u_velocity
                    .get([i, j + 1])
                    .expect("invariant: interior u east index is in bounds")
                    - *checkpoint
                        .u_velocity
                        .get([i, j - 1])
                        .expect("invariant: interior u west index is in bounds"))
                    / (central_difference * dx);
                let dvdy = (*checkpoint
                    .v_velocity
                    .get([i + 1, j])
                    .expect("invariant: interior v north index is in bounds")
                    - *checkpoint
                        .v_velocity
                        .get([i - 1, j])
                        .expect("invariant: interior v south index is in bounds"))
                    / (central_difference * dy);

                let divergence = (dudx + dvdy).abs();
                if divergence > max_divergence {
                    max_divergence = divergence;
                }
            }
        }

        max_divergence < tolerance
    }
}

#[allow(clippy::cast_precision_loss)] // invariant: values above f64's exact integer range are rejected before casting.
fn usize_to_exact_f64(value: usize) -> Option<f64> {
    if value == 0 || value as u128 > (1u128 << f64::MANTISSA_DIGITS) {
        None
    } else {
        Some(value as f64)
    }
}
