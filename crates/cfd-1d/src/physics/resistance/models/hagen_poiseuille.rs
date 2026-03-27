//! Hagen-Poiseuille resistance model for laminar flow in circular pipes.
//!
//! ## Hagen-Poiseuille Theorem
//!
//! **Theorem**: For steady, fully developed, laminar flow of an incompressible Newtonian fluid
//! in a straight circular pipe of constant cross-section, the pressure drop ΔP over a length L
//! is given by:
//!
//! ΔP = (32 μ L V) / D²
//!
//! where:
//! - μ is the dynamic viscosity [Pa·s]
//! - L is the pipe length [m]
//! - V is the average velocity [m/s]
//! - D is the pipe diameter [m]
//!
//! **Hydraulic Resistance Form**: R = ΔP / Q = (128 μ L) / (π D⁴)
//!
//! where Q is the volumetric flow rate [m³/s].
//!
//! ### Derivation
//!
//! Starting from the Navier-Stokes equations for axisymmetric flow:
//!
//! (1/ρ) ∇P = ν ∇²u + body forces
//!
//! For fully developed laminar flow (∂u/∂z = 0), the axial momentum equation reduces to:
//!
//! dP/dz = μ d²u/dr² + (μ/r) du/dr
//!
//! Integrating twice with boundary conditions u(r=R) = 0 and du/dr(r=0) = 0:
//!
//! u(r) = (1/(4μ)) (-dP/dz) (R² - r²)
//!
//! The Hagen-Poiseuille velocity profile is parabolic:
//!
//! u(r) = u_max (1 - (r/R)²)
//!
//! where u_max = (R²/(4μ)) (-dP/dz)
//!
//! Average velocity: V = (1/2) u_max = (R²/(8μ)) (-dP/dz)
//!
//! Pressure drop: ΔP = - (8μ L V) / R² = (32 μ L V) / D²
//!
//! ### Validity Conditions
//!
//! 1. **Laminar Flow**: Re < 2300 (typically Re < 2000 for safety)
//! 2. **Fully Developed Flow**: L/D > 10 (entrance effects negligible)
//! 3. **Newtonian Fluid**: Constant viscosity, no shear-thinning/thickening
//! 4. **Incompressible**: ρ constant, Mach < 0.3
//! 5. **No Body Forces**: Gravity and centrifugal forces negligible
//! 6. **Constant Properties**: Temperature constant, no property variations
//! 7. **Straight Pipe**: No curvature or bends
//! 8. **Smooth Walls**: No surface roughness effects
//!
//! ### References
//!
//! - Hagen, G. (1839). "Über die Bewegung der Flüssigkeiten in cylindrischen Röhren."
//!   *Poggendorff's Annalen der Physik und Chemie*, 46, 423-442.
//! - Poiseuille, J. L. M. (1840). "Recherches expérimentales sur le mouvement des liquides
//!   dans les tubes de très petits diamètres." *Mémoires présentés par divers savants à
//!   l'Académie Royale des Sciences de l'Institut de France*, 9, 433-544.
//! - White, F. M. (2006). *Viscous Fluid Flow* (3rd ed.). McGraw-Hill. Eq. 3-52.

use super::traits::{FlowConditions, ResistanceModel};
use cfd_core::error::Result;
use cfd_core::physics::fluid::FluidTrait;
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;
use serde::{Deserialize, Serialize};

// Named constants
const HAGEN_POISEUILLE_COEFFICIENT: f64 = 128.0;

/// Hagen-Poiseuille resistance model for circular channels
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HagenPoiseuilleModel<T: RealField + Copy> {
    /// Channel diameter [m]
    pub diameter: T,
    /// Channel length [m]
    pub length: T,
}

impl<T: RealField + Copy> HagenPoiseuilleModel<T> {
    /// Create a new Hagen-Poiseuille model
    pub fn new(diameter: T, length: T) -> Self {
        Self { diameter, length }
    }
}

impl<T: RealField + Copy + FromPrimitive> ResistanceModel<T> for HagenPoiseuilleModel<T> {
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
        // Calculate shear rate if not provided
        let shear_rate = if let Some(sr) = conditions.shear_rate {
            sr
        } else {
            let v = if let Some(vel) = conditions.velocity {
                vel
            } else if let Some(q) = conditions.flow_rate {
                let pi = T::pi();
                let area = pi * self.diameter * self.diameter
                    / (T::one() + T::one() + T::one() + T::one());
                q / area
            } else {
                T::zero()
            };
            T::from_f64(8.0).expect("Mathematical constant conversion compromised") * v
                / self.diameter
        };

        let viscosity =
            fluid.viscosity_at_shear(shear_rate, conditions.temperature, conditions.pressure)?;

        let pi = T::pi();

        let coefficient = T::from_f64(HAGEN_POISEUILLE_COEFFICIENT)
            .expect("Mathematical constant conversion compromised");

        // R = (128 * μ * L) / (π * D^4)
        let d2 = self.diameter * self.diameter;
        let d4 = d2 * d2;
        let r = coefficient * viscosity * self.length / (pi * d4);

        Ok((r, T::zero()))
    }

    fn model_name(&self) -> &'static str {
        "Hagen-Poiseuille"
    }

    fn reynolds_range(&self) -> (T, T) {
        (
            T::zero(),
            T::from_f64(2300.0).expect("Mathematical constant conversion compromised"),
        )
    }

    fn validate_invariants<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<()> {
        // Call Mach number validation
        self.validate_mach_number(fluid, conditions)?;

        if self.diameter <= T::zero() {
            return Err(cfd_core::error::Error::PhysicsViolation(format!(
                "Invalid geometry: diameter = {} <= 0 for model '{}'",
                self.diameter,
                self.model_name()
            )));
        }
        if self.length < T::zero() {
            return Err(cfd_core::error::Error::PhysicsViolation(format!(
                "Invalid geometry: length = {} < 0 for model '{}'",
                self.length,
                self.model_name()
            )));
        }

        // Entrance length validation: L/D > 10
        let ratio = self.length / self.diameter;
        let limit = T::from_f64(10.0).expect("Mathematical constant conversion compromised");

        if ratio < limit {
            return Err(cfd_core::error::Error::PhysicsViolation(format!(
                "Entrance length violation: L/D = {:.2} < 10. Flow may not be fully developed for model '{}'",
                ratio,
                self.model_name()
            )));
        }

        Ok(())
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
    fn resistance_formula_circular_pipe() {
        // R = 128 * mu * L / (pi * D^4)
        // mu=0.001, L=0.01, D=0.001
        // R = 128 * 0.001 * 0.01 / (pi * 1e-12) = 1.28e-3 / (pi * 1e-12)
        let model = HagenPoiseuilleModel::new(0.001_f64, 0.01_f64);
        let conditions = FlowConditions::new(0.0);
        let (r, _k) = model.calculate_coefficients(&water(), &conditions).unwrap();
        let expected = 128.0 * 0.001 * 0.01 / (std::f64::consts::PI * 1e-12);
        assert_relative_eq!(r, expected, max_relative = 1e-10);
        assert_relative_eq!(expected, 4.074e8, max_relative = 1e-3);
    }

    #[test]
    fn coefficients_zero_quadratic_term() {
        // Hagen-Poiseuille is purely linear: k = 0
        let model = HagenPoiseuilleModel::new(0.001_f64, 0.01_f64);
        let conditions = FlowConditions::new(0.0);
        let (_r, k) = model.calculate_coefficients(&water(), &conditions).unwrap();
        assert_relative_eq!(k, 0.0, epsilon = 1e-15);
    }

    #[test]
    fn very_small_diameter_produces_large_resistance() {
        let large = HagenPoiseuilleModel::new(0.001_f64, 0.01_f64);
        let small = HagenPoiseuilleModel::new(0.0001_f64, 0.01_f64);
        let conditions = FlowConditions::new(0.0);
        let (r_large, _) = large.calculate_coefficients(&water(), &conditions).unwrap();
        let (r_small, _) = small.calculate_coefficients(&water(), &conditions).unwrap();
        // D halved by 10x → R increases by 10^4 = 10000x (D^4 in denominator)
        assert!(r_small > r_large * 9000.0, "small diameter should produce much larger resistance");
        assert_relative_eq!(r_small / r_large, 1e4, max_relative = 1e-10);
    }
    
    #[test]
    fn negative_dimensions_rejected() {
        let conditions = FlowConditions::new(0.0);
        
        let neg_diam = HagenPoiseuilleModel::new(-0.001_f64, 0.01_f64);
        assert!(neg_diam.validate_invariants(&water(), &conditions).is_err());
        
        let zero_diam = HagenPoiseuilleModel::new(0.0_f64, 0.01_f64);
        assert!(zero_diam.validate_invariants(&water(), &conditions).is_err());
        
        let neg_length = HagenPoiseuilleModel::new(0.001_f64, -0.01_f64);
        assert!(neg_length.validate_invariants(&water(), &conditions).is_err());
    }
}
