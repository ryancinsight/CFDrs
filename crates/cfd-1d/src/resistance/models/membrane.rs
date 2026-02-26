//! Porous membrane resistance model for organ-on-chip and barrier simulations.
//!
//! ## Theorem: Hagen-Poiseuille Resistance of a Single Cylindrical Pore
//!
//! **Theorem** (exact for laminar, fully-developed Newtonian flow in a cylinder):
//!
//! For a single cylindrical pore of length L and radius r:
//!
//! ```text
//! R_pore = 8 μ L / (π r⁴)    [Pa·s/m³]
//! ```
//!
//! **Derivation**: For Poiseuille flow (parabolic profile), the volumetric flow rate is
//! `Q = π r⁴ ΔP / (8 μ L)`, giving `R = ΔP / Q = 8μL/(πr⁴)`.  ■
//!
//! ## Theorem: Parallel Pore Resistance
//!
//! **Theorem**: For N identical pores in parallel, the total resistance is:
//!
//! ```text
//! R_total = R_pore / N
//! ```
//!
//! The number of pores per unit area is given by the open-area porosity φ
//! (fraction of membrane area occupied by pores):
//!
//! ```text
//! N = φ · A_membrane / (π r²)
//! ```
//!
//! Substituting:
//!
//! ```text
//! R_total = R_pore / N
//!         = [8 μ L / (π r⁴)]  /  [φ A / (π r²)]
//!         = 8 μ L / (φ A r²)    [Pa·s/m³]
//! ```
//!   ■
//!
//! ## Theorem: Continuum Flow Validity (Knudsen Number)
//!
//! The Hagen-Poiseuille model is valid only in the **continuum regime**.
//! The Knudsen number Kn is:
//!
//! ```text
//! Kn = λ / r
//! ```
//!
//! where λ is the mean free path of the fluid molecules.
//! For continuum flow: Kn < 0.01.
//! For slip flow (0.01 < Kn < 0.1), a slip-flow correction applies.
//!
//! In typical organ-on-chip membranes (r ≈ 0.5–1 μm, water), Kn ~ 10⁻⁷ ≪ 0.01,
//! so the continuum assumption is fully justified.
//!
//! ## References
//!
//! - Deen, W. M. (1987). "Hindered transport of large molecules in liquid filled pores."
//!   *AIChE Journal*, 33(9), 1409–1425.
//! - Bungay, P. M., & Brenner, H. (1973). "The motion of a closely-fitting sphere in a
//!   fluid-filled tube." *International Journal of Multiphase Flow*, 1(1), 25–56.
//! - Truskey, G. A., Yuan, F., & Deen, D. F. (2010).
//!   *Transport Phenomena in Biological Systems* (2nd ed.). Prentice Hall. Ch. 3.

use super::traits::{FlowConditions, ResistanceModel};
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid::FluidTrait;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Porous membrane represented by equivalent parallel cylindrical pores.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MembranePoreModel<T: RealField + Copy> {
    /// Membrane thickness [m]
    pub thickness: T,
    /// Membrane width [m]
    pub width: T,
    /// Membrane height [m]
    pub height: T,
    /// Pore radius [m]
    pub pore_radius: T,
    /// Open-area fraction, in [0, 1]
    pub porosity: T,
}

impl<T: RealField + Copy> MembranePoreModel<T> {
    /// Create a new membrane pore model.
    pub fn new(thickness: T, width: T, height: T, pore_radius: T, porosity: T) -> Self {
        Self {
            thickness,
            width,
            height,
            pore_radius,
            porosity,
        }
    }
}

impl<T: RealField + Copy + FromPrimitive> ResistanceModel<T> for MembranePoreModel<T> {
    fn calculate_resistance<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<T> {
        self.validate_invariants(fluid, conditions)?;

        let state = fluid.properties_at(conditions.temperature, conditions.pressure)?;
        let mu = state.dynamic_viscosity;

        let area = self.width * self.height;
        let denom = self.porosity * area * self.pore_radius * self.pore_radius;
        let eight = T::from_f64(8.0).unwrap_or_else(T::zero);
        Ok(eight * mu * self.thickness / denom)
    }

    fn model_name(&self) -> &'static str {
        "MembranePore"
    }

    fn reynolds_range(&self) -> (T, T) {
        (T::zero(), T::from_f64(100.0).unwrap_or_else(T::one))
    }

    fn validate_invariants<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<()> {
        self.validate_mach_number(fluid, conditions)?;

        let zero = T::zero();
        let one = T::one();
        if self.thickness <= zero {
            return Err(Error::InvalidConfiguration(
                "Membrane thickness must be positive".to_string(),
            ));
        }
        if self.width <= zero || self.height <= zero {
            return Err(Error::InvalidConfiguration(
                "Membrane width and height must be positive".to_string(),
            ));
        }
        if self.pore_radius <= zero {
            return Err(Error::InvalidConfiguration(
                "Membrane pore radius must be positive".to_string(),
            ));
        }
        if self.porosity <= zero || self.porosity > one {
            return Err(Error::InvalidConfiguration(
                "Membrane porosity must be in (0, 1]".to_string(),
            ));
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use cfd_core::physics::fluid::database::water_20c;

    fn water() -> impl FluidTrait<f64> {
        water_20c::<f64>().unwrap()
    }

    /// Verify R_total matches the theoretical formula R = 8μL/(φAr²).
    ///
    /// Theorem: R = 8μL/(φAr²) for cylindrical parallel-pore membranes.
    #[test]
    fn test_membrane_resistance_matches_formula() {
        let model = MembranePoreModel::<f64>::new(
            10e-6,   // 10 μm thick
            1e-3,    // 1 mm wide
            1e-3,    // 1 mm high
            0.5e-6,  // 0.5 μm pore radius
            0.2,     // 20% porosity
        );
        let fluid = water();
        let cond = FlowConditions::new(0.0_f64);

        let r = model.calculate_resistance(&fluid, &cond).unwrap();

        // Analytical: mu = 1.002e-3, L=10e-6, phi=0.2, A=1e-6, r=0.5e-6
        // R = 8 * mu * L / (phi * A * r^2)
        let mu = 1.002e-3_f64; // water at 20°C [Pa·s]
        let expected = 8.0 * mu * 10e-6 / (0.2 * 1e-6 * 0.5e-6 * 0.5e-6);
        assert_relative_eq!(r, expected, max_relative = 0.02); // 2% tolerance for tabulated mu
    }

    /// Resistance must be positive (physical invariant).
    #[test]
    fn test_membrane_resistance_positive() {
        let model = MembranePoreModel::<f64>::new(10e-6, 1e-3, 1e-3, 0.5e-6, 0.2);
        let fluid = water();
        let cond = FlowConditions::new(0.0);
        let r = model.calculate_resistance(&fluid, &cond).unwrap();
        assert!(r > 0.0, "Membrane resistance must be positive, got {r}");
        assert!(r.is_finite(), "Membrane resistance must be finite, got {r}");
    }

    /// R ∝ L (linear scaling with thickness).
    ///
    /// Theorem: R = 8μL/(φAr²) is linear in L.
    #[test]
    fn test_membrane_resistance_linear_in_thickness() {
        let base = MembranePoreModel::<f64>::new(10e-6, 1e-3, 1e-3, 0.5e-6, 0.2);
        let double = MembranePoreModel::<f64>::new(20e-6, 1e-3, 1e-3, 0.5e-6, 0.2);
        let fluid = water();
        let cond = FlowConditions::new(0.0);
        let r_base = base.calculate_resistance(&fluid, &cond).unwrap();
        let r_double = double.calculate_resistance(&fluid, &cond).unwrap();
        // R ∝ L: doubling thickness should double resistance
        assert_relative_eq!(r_double / r_base, 2.0, epsilon = 1e-9);
    }

    /// R ∝ 1/φ (inverse porosity scaling).
    #[test]
    fn test_membrane_resistance_inverse_porosity() {
        let phi20 = MembranePoreModel::<f64>::new(10e-6, 1e-3, 1e-3, 0.5e-6, 0.2);
        let phi40 = MembranePoreModel::<f64>::new(10e-6, 1e-3, 1e-3, 0.5e-6, 0.4);
        let fluid = water();
        let cond = FlowConditions::new(0.0);
        let r_phi20 = phi20.calculate_resistance(&fluid, &cond).unwrap();
        let r_phi40 = phi40.calculate_resistance(&fluid, &cond).unwrap();
        // R ∝ 1/φ: doubling porosity should halve resistance
        assert_relative_eq!(r_phi20 / r_phi40, 2.0, epsilon = 1e-9);
    }

    /// R ∝ 1/r² (inverse pore-radius squared). 
    #[test]
    fn test_membrane_resistance_inverse_radius_squared() {
        let r1 = MembranePoreModel::<f64>::new(10e-6, 1e-3, 1e-3, 0.5e-6, 0.2);
        let r2 = MembranePoreModel::<f64>::new(10e-6, 1e-3, 1e-3, 1.0e-6, 0.2); // 2x radius
        let fluid = water();
        let cond = FlowConditions::new(0.0);
        let res1 = r1.calculate_resistance(&fluid, &cond).unwrap();
        let res2 = r2.calculate_resistance(&fluid, &cond).unwrap();
        // R ∝ 1/r² ⇒ R(r1) / R(r2) = (r2/r1)² = 4
        // R ∝ 1/r²: doubling pore radius should quarter resistance → ratio = 4
        assert_relative_eq!(res1 / res2, 4.0, epsilon = 1e-9);
    }

    /// Validate invariant: zero porosity returns error.
    #[test]
    fn test_membrane_zero_porosity_error() {
        let model = MembranePoreModel::<f64>::new(10e-6, 1e-3, 1e-3, 0.5e-6, 0.0);
        let fluid = water();
        let cond = FlowConditions::new(0.0);
        assert!(model.validate_invariants(&fluid, &cond).is_err(),
            "Zero porosity should return error");
    }

    /// Validate invariant: porosity > 1 returns error.
    #[test]
    fn test_membrane_porosity_exceeds_one_error() {
        let model = MembranePoreModel::<f64>::new(10e-6, 1e-3, 1e-3, 0.5e-6, 1.5);
        let fluid = water();
        let cond = FlowConditions::new(0.0);
        assert!(model.validate_invariants(&fluid, &cond).is_err(),
            "Porosity > 1 should return error");
    }

    /// Validate invariant: zero pore radius returns error.
    #[test]
    fn test_membrane_zero_pore_radius_error() {
        let model = MembranePoreModel::<f64>::new(10e-6, 1e-3, 1e-3, 0.0, 0.2);
        let fluid = water();
        let cond = FlowConditions::new(0.0);
        assert!(model.validate_invariants(&fluid, &cond).is_err(),
            "Zero pore radius should return error");
    }
}
