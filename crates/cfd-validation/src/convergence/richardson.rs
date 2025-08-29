//! Richardson extrapolation for grid-independent solutions
//!
//! Implements Richardson extrapolation following ASME V&V 20-2009 guidelines.

use cfd_core::error::{Error, Result};
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Richardson extrapolation calculator
///
/// Estimates grid-independent solutions using systematic grid refinement
#[derive(Debug, Clone)]
pub struct RichardsonExtrapolation<T: RealField + Copy> {
    /// Assumed or computed order of accuracy
    pub order: T,
    /// Grid refinement ratio (r = h_coarse / h_fine)
    pub refinement_ratio: T,
}

impl<T: RealField + Copy + FromPrimitive> RichardsonExtrapolation<T> {
    /// Create a new Richardson extrapolation with known order
    pub fn with_order(order: T, refinement_ratio: T) -> Result<Self> {
        if order <= T::zero() {
            return Err(Error::InvalidInput(
                "Order of accuracy must be positive".to_string(),
            ));
        }

        if refinement_ratio <= T::one() {
            return Err(Error::InvalidInput(
                "Refinement ratio must be greater than 1".to_string(),
            ));
        }

        Ok(Self {
            order,
            refinement_ratio,
        })
    }

    /// Create with standard second-order accuracy
    pub fn second_order(refinement_ratio: T) -> Result<Self> {
        let two = T::from_f64(2.0).ok_or_else(|| {
            Error::InvalidInput("Cannot represent 2.0 in numeric type T".to_string())
        })?;
        Self::with_order(two, refinement_ratio)
    }

    /// Extrapolate to zero grid spacing using two solutions
    ///
    /// # Arguments
    /// * `f_fine` - Solution on fine grid
    /// * `f_coarse` - Solution on coarse grid
    ///
    /// # Returns
    /// Extrapolated solution at h→0
    pub fn extrapolate(&self, f_fine: T, f_coarse: T) -> T {
        let r_p = self.refinement_ratio.powf(self.order);
        (r_p * f_fine - f_coarse) / (r_p - T::one())
    }

    /// Compute grid convergence index (GCI) following Roache (1998)
    ///
    /// GCI provides an error band for the grid-converged solution
    pub fn grid_convergence_index(&self, f_fine: T, f_coarse: T, safety_factor: T) -> T {
        let epsilon = (f_fine - f_coarse).abs();
        let r_p = self.refinement_ratio.powf(self.order);

        safety_factor * epsilon / (r_p - T::one())
    }

    /// Estimate order of accuracy from three grid solutions
    ///
    /// Uses the generalized Richardson extrapolation formula
    pub fn estimate_order(f_coarse: T, f_medium: T, f_fine: T, refinement_ratio: T) -> Result<T>
    where
        T: RealField + Copy + FromPrimitive,
    {
        let epsilon_21 = f_medium - f_fine;
        let epsilon_32 = f_coarse - f_medium;

        let epsilon_tolerance = T::from_f64(cfd_core::constants::numerical::solver::EPSILON_TOLERANCE)
            .ok_or_else(|| Error::InvalidInput(
                "Cannot represent EPSILON_TOLERANCE in numeric type T".to_string()
            ))?;
        
        if epsilon_21.abs() < epsilon_tolerance {
            return Err(Error::InvalidInput(
                "Solutions too close to estimate order".to_string(),
            ));
        }

        let ratio = epsilon_32 / epsilon_21;
        let order = ratio.ln() / refinement_ratio.ln();

        Ok(order)
    }

    /// Check if solutions are in asymptotic range
    ///
    /// Returns true if the convergence ratio is consistent with the expected order
    pub fn is_asymptotic(&self, f_coarse: T, f_medium: T, f_fine: T) -> bool {
        let epsilon_21 = (f_medium - f_fine).abs();
        let epsilon_32 = (f_coarse - f_medium).abs();

        let epsilon_tolerance = T::from_f64(cfd_core::constants::numerical::solver::EPSILON_TOLERANCE)
            .unwrap_or_else(|| T::from_f64(1e-10).expect("Failed to represent 1e-10 in numeric type T"));
        
        if epsilon_21 < epsilon_tolerance {
            return false;
        }

        let observed_ratio = epsilon_32 / epsilon_21;
        let expected_ratio = self.refinement_ratio.powf(self.order);

        // Check if within 10% of expected ratio
        let relative_diff = ((observed_ratio - expected_ratio) / expected_ratio).abs();
        let ten_percent = T::from_f64(0.1).expect("Failed to represent 0.1 in numeric type T");
        relative_diff < ten_percent
    }
}

/// Perform Richardson extrapolation with automatic order estimation
///
/// Uses three grid levels to estimate order and extrapolate
pub fn richardson_extrapolate<T>(solutions: &[T], grid_sizes: &[T]) -> Result<(T, T)>
where
    T: RealField + Copy + FromPrimitive,
{
    if solutions.len() < 2 || solutions.len() != grid_sizes.len() {
        return Err(Error::InvalidInput(
            "Need at least 2 solutions with corresponding grid sizes".to_string(),
        ));
    }

    // Sort by grid size (finest first)
    let mut paired: Vec<_> = grid_sizes
        .iter()
        .zip(solutions.iter())
        .map(|(h, f)| (*h, *f))
        .collect();
    paired.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let f_fine = paired[0].1;
    let f_coarse = paired[1].1;
    let h_fine = paired[0].0;
    let h_coarse = paired[1].0;

    let r21 = h_coarse / h_fine;  // Refinement ratio between fine and medium grids

    // Estimate order if we have 3 or more solutions
    let order = if solutions.len() >= 3 {
        let h_medium = paired[1].0;
        let h_coarse_actual = paired[2].0;
        let r32 = h_coarse_actual / h_medium;  // Refinement ratio between medium and coarse
        
        // Check if refinement ratios are uniform (within 1% tolerance)
        let one_percent = T::from_f64(0.01)
            .ok_or_else(|| Error::InvalidInput("Cannot represent 0.01 in numeric type T".to_string()))?;
        
        if ((r21 - r32).abs() / r21) > one_percent {
            return Err(Error::InvalidInput(format!(
                "Non-uniform grid refinement ratios ({:?} and {:?}) detected. \
                 Richardson extrapolation requires constant refinement ratio.",
                r21, r32
            )));
        }
        
        let f_medium = paired[1].1;
        let f_coarse_actual = paired[2].1;
        RichardsonExtrapolation::estimate_order(f_coarse_actual, f_medium, f_fine, r21)?
    } else {
        // Assume second order if not enough data
        T::from_f64(2.0).ok_or_else(|| 
            Error::InvalidInput("Cannot represent 2.0 in numeric type T".to_string())
        )?
    };

    let extrapolator = RichardsonExtrapolation::with_order(order, r21)?;
    let extrapolated = extrapolator.extrapolate(f_fine, f_coarse);

    Ok((extrapolated, order))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_richardson_second_order() {
        // Test with exact second-order convergence
        let extrapolator = RichardsonExtrapolation::<f64>::second_order(2.0).unwrap();

        // Solutions: f(h) = 1 + h²
        let f_fine = 1.0 + 0.01; // h = 0.1
        let f_coarse = 1.0 + 0.04; // h = 0.2

        let extrapolated = extrapolator.extrapolate(f_fine, f_coarse);
        assert_relative_eq!(extrapolated, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_order_estimation() {
        // Test order estimation with known convergence
        let f_coarse = 1.16; // h = 0.4, f = 1 + h²
        let f_medium = 1.04; // h = 0.2
        let f_fine = 1.01; // h = 0.1

        let order = RichardsonExtrapolation::<f64>::estimate_order(f_coarse, f_medium, f_fine, 2.0)
            .unwrap();

        assert_relative_eq!(order, 2.0, epsilon = 0.01);
    }

    #[test]
    fn test_gci_calculation() {
        let extrapolator = RichardsonExtrapolation::<f64>::second_order(2.0).unwrap();
        let f_fine = 1.01;
        let f_coarse = 1.04;
        let safety_factor = 1.25; // Recommended for 3+ grids

        let gci = extrapolator.grid_convergence_index(f_fine, f_coarse, safety_factor);

        // GCI should be small for well-converged solutions
        assert!(gci < 0.02);
    }
}
