//! Wall roughness classification and roughness function computation.
//!
//! ## Roughness Effects on the Log-Law
//!
//! For rough surfaces, the log-law becomes:
//! ```text
//! u⁺ = (1/κ) ln((y + k_s)/k_s) + B − ΔB(k_s⁺)
//! ```
//!
//! where:
//! - k_s = equivalent sand-grain roughness height \[m]
//! - k_s⁺ = k_s u_τ / ν (roughness Reynolds number)
//! - ΔB = roughness function (downward shift of the log-law)
//!
//! ## Reference
//!
//! Jiménez, J. (2004). Turbulent flows over rough walls. *Anna. Rev. Fluid Mech.*, 36, 173–196.

use crate::physics::turbulence::constants::KAPPA;
use eunomia::{FloatElement, RealField};

/// Wall roughness parameters.
#[derive(Debug, Clone)]
pub struct WallRoughness<T: RealField> {
    /// Equivalent sand-grain roughness height (k_s) \[m].
    pub equivalent_sand_grain: T,
    /// Roughness type classification.
    pub roughness_type: RoughnessType,
    /// Roughness function ΔB for log-law shift.
    pub roughness_function: T,
}

impl<T: RealField> WallRoughness<T> {
    /// Create smooth wall (k_s = 0).
    pub fn smooth() -> Self {
        Self {
            equivalent_sand_grain: T::ZERO,
            roughness_type: RoughnessType::Smooth,
            roughness_function: T::ZERO,
        }
    }

    /// Create rough wall with equivalent sand-grain roughness.
    pub fn rough(k_s: T) -> Self {
        let roughness_type = RoughnessType::classify(k_s);
        let roughness_function = Self::calculate_roughness_function(k_s, roughness_type);

        Self {
            equivalent_sand_grain: k_s,
            roughness_type,
            roughness_function,
        }
    }

    /// Calculate roughness function ΔB based on roughness type.
    fn calculate_roughness_function(k_s: T, roughness_type: RoughnessType) -> T {
        match roughness_type {
            RoughnessType::Smooth => T::ZERO,
            RoughnessType::TransitionallyRough => {
                let kappa = T::from_f64(KAPPA);
                (T::ONE / kappa) * <T as FloatElement>::ln(k_s) + T::from_f64(8.5)
                    - (T::ONE / kappa) * <T as FloatElement>::ln(T::from_f64(30.0))
                    - T::from_f64(8.5)
            }
            RoughnessType::FullyRough => {
                let kappa = T::from_f64(KAPPA);
                (T::ONE / kappa) * <T as FloatElement>::ln(k_s) + T::from_f64(8.5)
            }
            RoughnessType::VeryRough => {
                let kappa = T::from_f64(KAPPA);
                (T::ONE / kappa) * <T as FloatElement>::ln(k_s) + T::from_f64(13.0)
            }
        }
    }
}

/// Roughness type classification based on k_s⁺.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RoughnessType {
    /// Smooth surface (k_s⁺ ≈ 0).
    Smooth,
    /// Transitionally rough (0.1 < k_s⁺ < 10).
    TransitionallyRough,
    /// Fully rough (k_s⁺ > 10).
    FullyRough,
    /// Very rough surfaces (k_s⁺ > 1000).
    VeryRough,
}

impl RoughnessType {
    /// Classify roughness based on k_s⁺ value.
    pub fn classify<T: RealField>(k_s_plus: T) -> Self {
        let threshold_trans = T::from_f64(10.0);
        let threshold_very = T::from_f64(1000.0);

        if k_s_plus <= T::from_f64(0.1) {
            RoughnessType::Smooth
        } else if k_s_plus <= threshold_trans {
            RoughnessType::TransitionallyRough
        } else if k_s_plus <= threshold_very {
            RoughnessType::FullyRough
        } else {
            RoughnessType::VeryRough
        }
    }
}
