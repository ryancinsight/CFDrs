//! Kreutzer (2005) two-phase gas-liquid slug (Taylor) flow resistance model.
//!
//! ## Theorem — Interfacial Pressure Drop in Taylor Flow (Kreutzer 2005)
//!
//! **Theorem**: In a micro or millifluidic capillary tube containing segmented
//! gas-liquid (slug or Taylor) flow, the apparent friction factor $f_{app}$ of the
//! liquid slug exceeds the theoretical Hagen-Poiseuille single-phase friction factor 
//! ($f = 16/Re$) due to capillary forces and internal recirculation within the slug.
//!
//! Kreutzer et al. (2005) mathematically correlated this enhancement as:
//!
//! $f_{app} = \frac{16}{Re} \left[ 1 + a \frac{D}{L_s} \left(\frac{Re}{Ca}\right)^{1/3} \right]$
//!
//! where:
//! - $a \approx 0.17$ is an empirical geometric constant.
//! - $D$ is the capillary diameter [m].
//! - $L_s$ is the length of the liquid slug [m].
//! - $Re = \rho v D / \mu$ is the Reynolds number.
//! - $Ca = \mu v / \sigma$ is the Capillary number (governing interfacial tension).
//!
//! ### Proof: Linear Invariance of the Resistance Multiplier
//!
//! A foundational property of the Kreutzer correlation is the cancellation of velocity
//! in the exponentiated $(Re/Ca)$ quotient:
//!
//! $\frac{Re}{Ca} = \frac{\rho v D / \mu}{\mu v / \sigma} = \frac{\rho \sigma D}{\mu^2}$
//!
//! This dimensionless quantity characterizes the fluid independent of velocity 
//! (often defined as the Laplace or Suratman number). Consequently, the slug
//! flow resistance multiplier:
//!
//! $M_{slug} = 1 + 0.17 \frac{D}{L_s} \left(\frac{\rho \sigma D}{\mu^2}\right)^{1/3}$
//!
//! evaluates as purely geometric/thermophysical. Thus, the system preserves
//! strict linearity where $R_{total} = R_{linear} \cdot M_{slug}$ and $\Delta P$
//! scales strictly linearly with flow rate $Q$, exactly preserving $k = 0$.
//!
//! ### Validity Conditions
//! 1. Segmented Taylor flow regime.
//! 2. Capillary-dominated flow ($Ca < 0.05$).
//!
//! ### References
//! - Kreutzer, M. T. et al. (2005). "Inertial and interfacial effects on pressure 
//!   drop of Taylor flow in capillaries". *AIChE Journal*, 51(9), 2428-2440.

use super::traits::{FlowConditions, ResistanceModel};
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid::FluidTrait;
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;
use serde::{Deserialize, Serialize};

const KREUTZER_CONSTANT: f64 = 0.17;
const HAGEN_POISEUILLE_COEFFICIENT: f64 = 128.0;

/// Kreutzer (2005) gas-liquid slug flow resistance model
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SlugFlowModel<T: RealField + Copy> {
    /// Capillary tube diameter [m]
    pub diameter: T,
    /// Total channel length [m]
    pub length: T,
    /// Liquid slug length (distance between gas bubbles) [m]
    pub slug_length: T,
    /// Surface tension coefficient of the gas-liquid interface [N/m]
    pub surface_tension: T,
}

impl<T: RealField + Copy> SlugFlowModel<T> {
    /// Create a new Slug Flow resistance model
    pub fn new(diameter: T, length: T, slug_length: T, surface_tension: T) -> Self {
        Self {
            diameter,
            length,
            slug_length,
            surface_tension,
        }
    }
}

impl<T: RealField + Copy + FromPrimitive> ResistanceModel<T> for SlugFlowModel<T> {
    fn calculate_resistance<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<T> {
        let (r, _) = self.calculate_coefficients(fluid, conditions)?;
        Ok(r)
    }

    fn calculate_coefficients<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<(T, T)> {
        let state = fluid.properties_at(conditions.temperature, conditions.pressure)?;
        let density = state.density;
        
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

        // Proof: Re/Ca = (ρ D σ) / (μ²)
        let re_over_ca = (density * self.diameter * self.surface_tension) / (viscosity * viscosity);
        
        // Multiplier = 1.0 + 0.17 * (D / L_s) * (Re/Ca)^(1/3)
        let one = T::one();
        let kreutzer = T::from_f64(KREUTZER_CONSTANT).expect("Mathematical conversion compromised");
        let geometric_ratio = self.diameter / self.slug_length;
        let third = one / T::from_f64(3.0).expect("Mathematical conversion compromised");
        
        // Since Re/Ca is strictly non-negative (density, dia, tension, viscosity are positive)
        let multiplier = one + kreutzer * geometric_ratio * re_over_ca.powf(third);

        // Base Hagen-Poiseuille resistance
        let pi = T::pi();
        let coefficient = T::from_f64(HAGEN_POISEUILLE_COEFFICIENT)
            .expect("Mathematical constant conversion compromised");
        let d2 = self.diameter * self.diameter;
        let d4 = d2 * d2;
        let r_linear = coefficient * viscosity * self.length / (pi * d4);

        // Effective total resistance
        let r_slug = r_linear * multiplier;

        Ok((r_slug, T::zero()))
    }

    fn model_name(&self) -> &'static str {
        "Kreutzer Two-Phase Slug Flow"
    }

    fn reynolds_range(&self) -> (T, T) {
        (
            T::zero(),
            T::from_f64(2000.0).expect("Mathematical conversion compromised"),
        )
    }

    fn validate_invariants<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<()> {
        self.validate_mach_number(fluid, conditions)?;

        if self.diameter <= T::zero() {
            return Err(Error::PhysicsViolation("Diameter must be positive".into()));
        }
        if self.slug_length <= T::zero() {
            return Err(Error::PhysicsViolation("Slug length must be positive".into()));
        }
        if self.length < T::zero() {
            return Err(Error::PhysicsViolation("Length must be non-negative".into()));
        }
        if self.surface_tension <= T::zero() {
            return Err(Error::PhysicsViolation("Surface tension must be positive".into()));
        }

        // Invariant: Slug length must be physically meaningful relative to diameter
        // Typically L_s > D for distinct liquid slugs. If L_s < D, the flow is 
        // likely annular or dispersed bubbly flow, not Taylor flow.
        if self.slug_length < self.diameter {
             return Err(Error::PhysicsViolation(format!(
                "Slug Flow invariant violation: L_s ({}) < D ({}). Taylor flow assumes liquid slug length exceeds diameter.",
                self.slug_length, self.diameter
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

    fn test_fluid() -> ConstantPropertyFluid<f64> {
        ConstantPropertyFluid::new(
            "water_air".to_string(),
            1000.0, // density
            0.001,  // viscosity
            4186.0,
            0.6,
            1500.0,
        )
    }

    #[test]
    fn slug_flow_multiplier_preserves_linear_invariance() {
        // D = 1mm, L = 10mm, L_s = 2mm, surface tension = 0.072 N/m
        let model = SlugFlowModel::new(0.001_f64, 0.01_f64, 0.002_f64, 0.072_f64);
        
        let fluid = test_fluid();
        
        let cond1 = FlowConditions::new(0.01); // v = 0.01 m/s
        let cond2 = FlowConditions::new(10.0); // v = 10.0 m/s
        
        let (r1, k1) = model.calculate_coefficients(&fluid, &cond1).unwrap();
        let (r2, k2) = model.calculate_coefficients(&fluid, &cond2).unwrap();
        
        // Exact mathematical proof validation: The multiplier is velocity independent
        assert_relative_eq!(r1, r2, epsilon = 1e-12);
        
        // Strictly linear behavior
        assert_relative_eq!(k1, 0.0, epsilon = 1e-15);
        assert_relative_eq!(k2, 0.0, epsilon = 1e-15);
        
        // Base resistance comparison
        let pi = std::f64::consts::PI;
        let r_base = (128.0 * 0.001 * 0.01) / (pi * 1e-12);
        
        // Re/Ca ratio = (1000 * 0.001 * 0.072) / (0.001)^2 = 72,000
        // Multiplier = 1 + 0.17 * (1mm / 2mm) * (72000)^(1/3)
        // 72000^(1/3) = ~41.6016
        // Mult = 1 + 0.17 * 0.5 * 41.6016 = 4.536
        let expected_mult = 1.0 + 0.17 * 0.5 * (72000.0_f64).powf(1.0/3.0);
        let expected_r = r_base * expected_mult;
        
        assert_relative_eq!(r1, expected_r, max_relative = 1e-10);
    }
    
    #[test]
    fn invariant_validation_rejects_unphysical_slugs() {
        // L_s = 0.0005 < D = 0.001 -> Dispersed bubbly, not Taylor flow
        let model = SlugFlowModel::new(0.001_f64, 0.01_f64, 0.0005_f64, 0.072_f64);
        let result = model.validate_invariants(&test_fluid(), &FlowConditions::new(0.1));
        assert!(matches!(result, Err(Error::PhysicsViolation(_))));
    }
    
    #[test]
    fn negative_dimensions_rejected() {
        let fluid = test_fluid();
        let cond = FlowConditions::new(0.1);
        
        // Zero diameter
        let m1 = SlugFlowModel::new(0.0, 0.01, 0.002, 0.072);
        assert!(m1.validate_invariants(&fluid, &cond).is_err());
        
        // Negative length
        let m2 = SlugFlowModel::new(0.001, -0.01, 0.002, 0.072);
        assert!(m2.validate_invariants(&fluid, &cond).is_err());
        
        // Zero surface tension
        let m3 = SlugFlowModel::new(0.001, 0.01, 0.002, 0.0);
        assert!(m3.validate_invariants(&fluid, &cond).is_err());
    }
}
