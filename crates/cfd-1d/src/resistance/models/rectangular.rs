//! Rectangular channel resistance model with exact solution.
//!
//! ## Shah-London Correlations for Rectangular Ducts
//!
//! **Theorem**: For laminar flow in rectangular ducts, the Poiseuille number
//! Po = f * Re varies with aspect ratio α = b/h (where b ≥ h):
//!
//! ### Poiseuille Number Correlations (Shah & London, 1978)
//!
//! **Square Duct (α = 1.0)**: Po = 56.91
//!
//! **Wide Rectangular Duct (α ≥ 1.0)**:
//! Po = 24 * (1 - 1.3553/α + 1.9467/α² - 1.7012/α³ + 0.9564/α⁴ - 0.2537/α⁵)
//!
//! **Tall Rectangular Duct (α ≤ 1.0)**:
//! Po = 24 * (1 - 0.628/α + 0.300/α² - 0.238/α³ + 0.091/α⁴ - 0.013/α⁵)
//!
//! where:
//! - Po is the Poiseuille number (dimensionless)
//! - f is the Fanning friction factor
//! - Re is the Reynolds number based on hydraulic diameter
//! - α is the aspect ratio (width/height, α ≥ 1)
//!
//! ### Hydraulic Resistance
//!
//! R = (Po μ L) / (Re ρ A Dh)
//!
//! where:
//! - μ is dynamic viscosity [Pa·s]
//! - L is channel length [m]
//! - ρ is fluid density [kg/m³]
//! - A is cross-sectional area [m²]
//! - Dh is hydraulic diameter [m]
//!
//! ### Validity Conditions
//!
//! 1. **Laminar Flow**: Re < 2300 (based on hydraulic diameter)
//! 2. **Fully Developed Flow**: L/Dh > 10 (entrance effects negligible)
//! 3. **Aspect Ratio**: 0.1 ≤ α ≤ 10 (recommended range)
//! 4. **Newtonian Fluid**: Constant viscosity
//! 5. **No Surface Roughness**: Smooth walls
//!
//! ### References
//!
//! - Shah, R. K., & London, A. L. (1978). *Laminar Flow Forced Convection in Ducts*.
//!   Academic Press. Chapter 7.
//! - White, F. M. (2006). *Viscous Fluid Flow* (3rd ed.). McGraw-Hill. Section 3.6.

use super::traits::{FlowConditions, ResistanceModel};
use cfd_core::error::Result;
use cfd_core::physics::fluid::Fluid;
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;
use serde::{Deserialize, Serialize};

// Named constants for hydraulic calculations
const HYDRAULIC_DIAMETER_FACTOR: f64 = 4.0;
const PERIMETER_FACTOR: f64 = 2.0;

/// Rectangular channel resistance model with exact solution
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RectangularChannelModel<T: RealField + Copy> {
    /// Channel width [m]
    pub width: T,
    /// Channel height [m]
    pub height: T,
    /// Channel length [m]
    pub length: T,
}

impl<T: RealField + Copy> RectangularChannelModel<T> {
    /// Create a new rectangular channel model
    pub fn new(width: T, height: T, length: T) -> Self {
        Self {
            width,
            height,
            length,
        }
    }
}

impl<T: RealField + Copy + FromPrimitive> ResistanceModel<T>
    for RectangularChannelModel<T>
{
    fn calculate_resistance(&self, fluid: &Fluid<T>, conditions: &FlowConditions<T>) -> Result<T> {
        let (r, k) = self.calculate_coefficients(fluid, conditions)?;
        let q = conditions.flow_rate.unwrap_or_else(T::zero);
        let q_abs = if q >= T::zero() { q } else { -q };
        Ok(r + k * q_abs)
    }

    fn calculate_coefficients(&self, fluid: &Fluid<T>, _conditions: &FlowConditions<T>) -> Result<(T, T)> {
        let viscosity = fluid.viscosity;

        // Calculate hydraulic diameter
        let area = self.width * self.height;
        let perimeter =
            T::from_f64(PERIMETER_FACTOR).unwrap_or_else(|| T::zero()) * (self.width + self.height);
        let dh =
            T::from_f64(HYDRAULIC_DIAMETER_FACTOR).unwrap_or_else(|| T::zero()) * area / perimeter;

        // Calculate aspect ratio (always ≥ 1 for consistency)
        let aspect_ratio =
            RealField::max(self.width, self.height) / RealField::min(self.width, self.height);

        // Calculate Poiseuille number using Shah-London correlations
        // Po = f_darcy * Re (Note: Shah-London often use f_fanning * Re = Po_darcy / 4)
        let poiseuille_number = self.calculate_poiseuille_number(aspect_ratio);

        // For laminar flow, the hydraulic resistance is constant:
        // R = (Po * mu * L) / (2 * A * Dh^2)
        let r = (poiseuille_number * viscosity * self.length)
            / (T::from_f64(2.0).unwrap_or_else(|| T::zero()) * area * dh * dh);

        Ok((r, T::zero()))
    }

    fn model_name(&self) -> &'static str {
        "Rectangular Channel (Exact)"
    }

    fn validate_invariants(&self, fluid: &Fluid<T>, conditions: &FlowConditions<T>) -> Result<()> {
        // Call Mach number validation
        self.validate_mach_number(fluid, conditions)?;

        // Entrance length validation: L/Dh > 10
        let area = self.width * self.height;
        let perimeter =
            T::from_f64(PERIMETER_FACTOR).unwrap_or_else(|| T::zero()) * (self.width + self.height);
        let dh =
            T::from_f64(HYDRAULIC_DIAMETER_FACTOR).unwrap_or_else(|| T::zero()) * area / perimeter;

        let ratio = self.length / dh;
        let limit = T::from_f64(10.0).unwrap_or_else(|| T::zero());

        if ratio < limit {
            return Err(cfd_core::error::Error::PhysicsViolation(format!(
                "Entrance length violation: L/Dh = {:.2} < 10. Flow may not be fully developed for model '{}'",
                ratio,
                self.model_name()
            )));
        }

        Ok(())
    }

    fn reynolds_range(&self) -> (T, T) {
        (
            T::zero(),
            T::from_f64(cfd_core::physics::constants::physics::dimensionless::reynolds::PIPE_LAMINAR_MAX)
                .unwrap_or_else(|| T::zero()),
        )
    }
}

impl<T: RealField + Copy + FromPrimitive> RectangularChannelModel<T> {
    /// Calculate Poiseuille number for rectangular channels using Shah-London correlations
    ///
    /// # Implementation Details
    /// Based on Shah & London (1978) exact series solutions for rectangular ducts:
    /// - Po = f_D * Re where f_D is the Darcy friction factor
    /// - Valid for aspect ratios between 0.1 and 10.0
    /// - Exact solution for laminar flow in rectangular ducts
    ///
    /// # Applicable Range
    /// - Reynolds number: Re < 2300 (laminar flow)
    /// - Aspect ratio: 0.1 ≤ α ≤ 10.0 (recommended)
    /// - Relative roughness: ε/Dh < 0.05
    fn calculate_poiseuille_number(&self, aspect_ratio: T) -> T {
        // Ensure aspect ratio is >= 1 (α = max(w,h)/min(w,h))
        let alpha = RealField::max(aspect_ratio, T::one() / aspect_ratio);

        if alpha == T::one() {
            // Square duct: Po_Darcy = 56.91
            T::from_f64(56.91).unwrap_or_else(|| T::zero())
        } else {
            // Po_Darcy = 4 * Po_Fanning
            // Po_Fanning = 24 * (1 - 1.3553/α + 1.9467/α² - 1.7012/α³ + 0.9564/α⁴ - 0.2537/α⁵)
            let a1 = T::from_f64(1.3553).unwrap_or_else(|| T::zero());
            let a2 = T::from_f64(1.9467).unwrap_or_else(|| T::zero());
            let a3 = T::from_f64(1.7012).unwrap_or_else(|| T::zero());
            let a4 = T::from_f64(0.9564).unwrap_or_else(|| T::zero());
            let a5 = T::from_f64(0.2537).unwrap_or_else(|| T::zero());

            let base = T::from_f64(96.0).unwrap_or_else(|| T::zero()); // 4 * 24.0
            let correction = T::one() - a1 / alpha + a2 / (alpha * alpha)
                - a3 / (alpha * alpha * alpha)
                + a4 / (alpha * alpha * alpha * alpha)
                - a5 / (alpha * alpha * alpha * alpha * alpha);

            base * correction
        }
    }
}
