//! Checkpoint validation utilities

use crate::checkpoint::Checkpoint;
use nalgebra::RealField;

/// Numerical constants for finite difference calculations
mod numerical_constants {
    /// Central difference denominator factor (h_forward - h_backward = 2h)
    pub const CENTRAL_DIFFERENCE_FACTOR: f64 = 2.0;
}

/// Checkpoint validator for ensuring data integrity
pub struct CheckpointValidator;

impl CheckpointValidator {
    /// Validate checkpoint for physical consistency
    pub fn validate_physics<T: RealField + Copy>(checkpoint: &Checkpoint<T>) -> Result<(), String> {
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
            for val in temp.iter() {
                if *val <= T::zero() {
                    return Err("Temperature contains non-positive values".to_string());
                }
            }
        }

        if let Some(ref k) = checkpoint.turbulence_k {
            if !Self::is_field_finite(k) {
                return Err("Turbulence k contains non-finite values".to_string());
            }

            // Turbulent kinetic energy should be non-negative
            for val in k.iter() {
                if *val < T::zero() {
                    return Err("Turbulence k contains negative values".to_string());
                }
            }
        }

        if let Some(ref eps) = checkpoint.turbulence_epsilon {
            if !Self::is_field_finite(eps) {
                return Err("Turbulence epsilon contains non-finite values".to_string());
            }

            // Dissipation rate should be positive
            for val in eps.iter() {
                if *val <= T::zero() {
                    return Err("Turbulence epsilon contains non-positive values".to_string());
                }
            }
        }

        Ok(())
    }

    /// Check if all values in a matrix are finite
    fn is_field_finite<T: RealField>(field: &nalgebra::DMatrix<T>) -> bool {
        field.iter().all(|v| v.is_finite())
    }

    /// Check mass conservation (divergence-free condition)
    pub fn check_mass_conservation<T: RealField + Copy>(
        checkpoint: &Checkpoint<T>,
        tolerance: T,
    ) -> bool {
        let (ny, nx) = checkpoint.dimensions();
        let (domain_x, domain_y) = checkpoint.metadata.domain_size;

        let dx = T::from_f64(domain_x / (nx as f64)).unwrap_or_else(T::one);
        let dy = T::from_f64(domain_y / (ny as f64)).unwrap_or_else(T::one);

        let mut max_divergence = T::zero();

        // Check interior points
        for i in 1..ny - 1 {
            for j in 1..nx - 1 {
                let dudx = (checkpoint.u_velocity[(i, j + 1)] - checkpoint.u_velocity[(i, j - 1)])
                    / (T::from_f64(numerical_constants::CENTRAL_DIFFERENCE_FACTOR)
                        .unwrap_or_else(T::one)
                        * dx);
                let dvdy = (checkpoint.v_velocity[(i + 1, j)] - checkpoint.v_velocity[(i - 1, j)])
                    / (T::from_f64(numerical_constants::CENTRAL_DIFFERENCE_FACTOR)
                        .unwrap_or_else(T::one)
                        * dy);

                let divergence = (dudx + dvdy).abs();
                if divergence > max_divergence {
                    max_divergence = divergence;
                }
            }
        }

        max_divergence < tolerance
    }
}
