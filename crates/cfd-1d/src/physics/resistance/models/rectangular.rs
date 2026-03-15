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
use cfd_core::physics::fluid::FluidTrait;
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

impl<T: RealField + Copy + FromPrimitive> ResistanceModel<T> for RectangularChannelModel<T> {
    fn calculate_resistance<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<T> {
        let (r, k) = self.calculate_coefficients(fluid, conditions)?;
        let q = conditions.flow_rate.unwrap_or_else(T::zero);
        let q_abs = if q >= T::zero() { q } else { -q };
        Ok(r + k * q_abs)
    }

    fn calculate_coefficients<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<(T, T)> {
        let area = self.width * self.height;
        let dh = self.hydraulic_diameter();

        let velocity = if let Some(v) = conditions.velocity {
            v
        } else if let Some(q) = conditions.flow_rate {
            q / area
        } else {
            T::zero()
        };
        let v_abs = if velocity >= T::zero() {
            velocity
        } else {
            -velocity
        };

        // Always query the shear-dependent viscosity. Newtonian fluids gracefully
        // ignore the shear_rate argument via the default trait implementation.
        let shear_rate = T::from_f64(8.0).unwrap_or(T::one()) * v_abs / dh;
        let viscosity = fluid
            .viscosity_at_shear(shear_rate, conditions.temperature, conditions.pressure)
            .unwrap_or(
                fluid
                    .properties_at(conditions.temperature, conditions.pressure)?
                    .dynamic_viscosity,
            );

        // Calculate aspect ratio (always ≥ 1 for consistency)
        let aspect_ratio =
            RealField::max(self.width, self.height) / RealField::min(self.width, self.height);

        // Calculate Poiseuille number using Shah-London correlations
        // Po = f_darcy * Re (Note: Shah-London often use f_fanning * Re = Po_darcy / 4)
        let poiseuille_number = self.calculate_poiseuille_number(aspect_ratio);

        // For laminar flow, the hydraulic resistance is constant:
        // R = (Po * mu * L) / (2 * A * Dh^2)
        let r = (poiseuille_number * viscosity * self.length)
            / ((T::one() + T::one()) * area * dh * dh);

        Ok((r, T::zero()))
    }

    fn model_name(&self) -> &'static str {
        "Rectangular Channel (Exact)"
    }

    fn validate_invariants<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<()> {
        // Call Mach number validation
        self.validate_mach_number(fluid, conditions)?;

        // Entrance length validation: L/Dh > 10
        let dh = self.hydraulic_diameter();

        let ratio = self.length / dh;
        let limit = T::from_f64(10.0).expect("Mathematical constant conversion compromised");

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
        (T::zero(), T::from_f64(2300.0).unwrap_or_else(|| T::zero()))
    }
}

impl<T: RealField + Copy + FromPrimitive> RectangularChannelModel<T> {
    /// Hydraulic diameter: Dh = 4A/P = 2wh/(w+h)
    fn hydraulic_diameter(&self) -> T {
        let area = self.width * self.height;
        let perimeter = T::from_f64(PERIMETER_FACTOR)
            .expect("Mathematical constant conversion compromised")
            * (self.width + self.height);
        T::from_f64(HYDRAULIC_DIAMETER_FACTOR)
            .expect("Mathematical constant conversion compromised")
            * area
            / perimeter
    }

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
            T::from_f64(56.91).expect("Mathematical constant conversion compromised")
        } else {
            // Po_Darcy = 4 * Po_Fanning
            // Po_Fanning = 24 * (1 - 1.3553/α + 1.9467/α² - 1.7012/α³ + 0.9564/α⁴ - 0.2537/α⁵)
            let a1 = T::from_f64(1.3553).expect("Mathematical constant conversion compromised");
            let a2 = T::from_f64(1.9467).expect("Mathematical constant conversion compromised");
            let a3 = T::from_f64(1.7012).expect("Mathematical constant conversion compromised");
            let a4 = T::from_f64(0.9564).expect("Mathematical constant conversion compromised");
            let a5 = T::from_f64(0.2537).expect("Mathematical constant conversion compromised");

            let base = T::from_f64(96.0).expect("Mathematical constant conversion compromised"); // 4 * 24.0
            let inv = T::one() / alpha;
            let inv2 = inv * inv;
            let inv3 = inv2 * inv;
            let inv4 = inv3 * inv;
            let inv5 = inv4 * inv;
            let correction =
                T::one() - a1 * inv + a2 * inv2 - a3 * inv3 + a4 * inv4 - a5 * inv5;

            base * correction
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use cfd_core::physics::fluid::ConstantPropertyFluid;

    fn water() -> ConstantPropertyFluid<f64> {
        ConstantPropertyFluid::new(
            "water".to_string(),
            1000.0, // density [kg/m³]
            0.001,  // viscosity [Pa·s]
            4186.0, // specific heat
            0.598,  // thermal conductivity
            1480.0, // speed of sound
        )
    }

    #[test]
    fn square_duct_poiseuille_number() {
        // Square duct (alpha=1.0) should return Po = 56.91
        let model = RectangularChannelModel::new(0.001_f64, 0.001_f64, 0.01_f64);
        let po = model.calculate_poiseuille_number(1.0);
        assert_relative_eq!(po, 56.91, max_relative = 1e-10);
    }

    #[test]
    fn aspect_ratio_symmetry() {
        // width/height = height/width should give same resistance
        let model_wh = RectangularChannelModel::new(0.001_f64, 0.0005_f64, 0.01_f64);
        let model_hw = RectangularChannelModel::new(0.0005_f64, 0.001_f64, 0.01_f64);
        let conditions = FlowConditions::new(0.0);
        let (r_wh, _) = model_wh
            .calculate_coefficients(&water(), &conditions)
            .unwrap();
        let (r_hw, _) = model_hw
            .calculate_coefficients(&water(), &conditions)
            .unwrap();
        assert_relative_eq!(r_wh, r_hw, max_relative = 1e-10);
    }

    #[test]
    fn hydraulic_diameter_rectangle() {
        // Dh = 4*A/P = 4*w*h / (2*(w+h)) = 2*w*h/(w+h)
        // w=0.001, h=0.0005 → Dh = 2*0.001*0.0005/(0.001+0.0005) = 1e-6/1.5e-3 ≈ 6.667e-4
        let model = RectangularChannelModel::new(0.001_f64, 0.0005_f64, 0.01_f64);
        let dh = model.hydraulic_diameter();
        let expected = 2.0 * 0.001 * 0.0005 / (0.001 + 0.0005);
        assert_relative_eq!(dh, expected, max_relative = 1e-10);
        assert_relative_eq!(dh, 6.667e-4, max_relative = 1e-3);
    }
}
