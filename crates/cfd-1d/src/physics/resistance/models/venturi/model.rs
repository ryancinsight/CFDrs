//! Venturi tube resistance model implementation.
//!
//! Contains the [`VenturiModel`] struct with the [`ResistanceModel`] trait
//! implementation, including Bernoulli contraction, throat friction, and
//! Borda-Carnot expansion loss calculations.

use super::traits::{FlowConditions, ResistanceModel};
use super::{ExpansionType, VenturiGeometry, BLASIUS_COEFF, BLASIUS_EXP, LAMINAR_FRICTION_COEFF, LAMINAR_LIMIT_RE};
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid::FluidTrait;
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Venturi tube resistance model for converging-diverging channels
///
/// Computes the net pressure drop across a Venturi tube considering:
/// 1. Bernoulli acceleration (contraction)
/// 2. Discharge coefficient (viscous correction at contraction)
/// 3. Throat friction (Darcy-Weisbach in the throat)
/// 4. Expansion recovery/loss (diffuser section)
///
/// The model supports both Newtonian and non-Newtonian fluids through
/// the `FluidTrait<T>` generic parameter.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VenturiModel<T: RealField + Copy> {
    /// Upstream (inlet) diameter [m]
    pub inlet_diameter: T,
    /// Throat diameter [m]
    pub throat_diameter: T,
    /// Downstream (outlet) diameter [m]
    pub outlet_diameter: T,
    /// Throat length [m]
    pub throat_length: T,
    /// Total device length [m] (for entrance length validation)
    pub total_length: T,
    /// Venturi geometry type (determines discharge coefficient)
    pub geometry_type: VenturiGeometry,
    /// Expansion type (determines recovery efficiency)
    pub expansion_type: ExpansionType,
    /// Surface roughness in throat [m]
    pub throat_roughness: T,
}

impl<T: RealField + Copy + FromPrimitive> VenturiModel<T> {
    /// Create a new Venturi model with specified geometry
    pub fn new(
        inlet_diameter: T,
        throat_diameter: T,
        outlet_diameter: T,
        throat_length: T,
        total_length: T,
    ) -> Self {
        Self {
            inlet_diameter,
            throat_diameter,
            outlet_diameter,
            throat_length,
            total_length,
            geometry_type: VenturiGeometry::MachinedConvergent,
            expansion_type: ExpansionType::Gradual {
                half_angle_deg: 7.0,
            },
            throat_roughness: T::zero(),
        }
    }

    /// Create a symmetric Venturi (outlet = inlet diameter)
    pub fn symmetric(
        inlet_diameter: T,
        throat_diameter: T,
        throat_length: T,
        total_length: T,
    ) -> Self {
        Self::new(
            inlet_diameter,
            throat_diameter,
            inlet_diameter,
            throat_length,
            total_length,
        )
    }

    /// Create a millifluidic Venturi with typical parameters
    pub fn millifluidic(inlet_diameter: T, throat_diameter: T, throat_length: T) -> Self {
        let total_length =
            throat_length * T::from_f64(5.0).expect("Mathematical constant conversion compromised");
        Self {
            inlet_diameter,
            throat_diameter,
            outlet_diameter: inlet_diameter,
            throat_length,
            total_length,
            geometry_type: VenturiGeometry::Custom {
                discharge_coefficient: 0.97,
            },
            expansion_type: ExpansionType::Gradual {
                half_angle_deg: 5.0,
            },
            throat_roughness: T::zero(),
        }
    }

    /// Set the geometry type
    pub fn with_geometry(mut self, geometry_type: VenturiGeometry) -> Self {
        self.geometry_type = geometry_type;
        self
    }

    /// Set the expansion type
    pub fn with_expansion(mut self, expansion_type: ExpansionType) -> Self {
        self.expansion_type = expansion_type;
        self
    }

    /// Set throat roughness
    pub fn with_roughness(mut self, roughness: T) -> Self {
        self.throat_roughness = roughness;
        self
    }

    /// Area ratio β² = (A_throat / A_inlet) = (D_throat / D_inlet)²
    ///
    /// Note: the Bernoulli pressure term requires `(1 − β⁴)` = `(1 − β²·β²)` because:
    ///
    /// ```text
    /// ΔP = ½ρV_t²[1 − (A_t/A_i)²]   where (A_t/A_i)² = (D_t/D_i)⁴ = β²·β²
    /// ```
    ///
    /// Callers must use `one - beta_sq * beta_sq`, **not** `one - beta_sq`.
    pub(crate) fn beta_squared(&self) -> T {
        let ratio = self.throat_diameter / self.inlet_diameter;
        ratio * ratio
    }

    /// Inlet cross-sectional area [m²]
    pub(crate) fn inlet_area(&self) -> T {
        let pi = T::pi();
        pi * self.inlet_diameter * self.inlet_diameter / (T::one() + T::one() + T::one() + T::one())
    }

    /// Throat cross-sectional area [m²]
    pub(crate) fn throat_area(&self) -> T {
        let pi = T::pi();
        pi * self.throat_diameter * self.throat_diameter
            / (T::one() + T::one() + T::one() + T::one())
    }

    /// Outlet cross-sectional area [m²]
    pub(crate) fn outlet_area(&self) -> T {
        let pi = T::pi();
        pi * self.outlet_diameter * self.outlet_diameter
            / (T::one() + T::one() + T::one() + T::one())
    }

    /// Calculate the Darcy friction factor in the throat
    ///
    /// Uses Hagen-Poiseuille (laminar) or Blasius (turbulent) correlations
    pub(crate) fn throat_friction_factor(&self, reynolds_throat: T) -> T {
        let re_lam =
            T::from_f64(LAMINAR_LIMIT_RE).expect("Mathematical constant conversion compromised");

        if reynolds_throat < re_lam {
            // Laminar: f = 64/Re
            T::from_f64(LAMINAR_FRICTION_COEFF)
                .expect("Mathematical constant conversion compromised")
                / reynolds_throat
        } else {
            // Blasius: f = 0.3164 / Re^0.25
            let coeff =
                T::from_f64(BLASIUS_COEFF).expect("Mathematical constant conversion compromised");
            let exp =
                T::from_f64(BLASIUS_EXP).expect("Mathematical constant conversion compromised");
            coeff / reynolds_throat.powf(exp)
        }
    }

    /// Effective discharge coefficient adjusted for Reynolds number
    ///
    /// At low Re (millifluidic regime), C_d decreases due to boundary layer growth.
    /// Empirical correction: C_d_eff = C_d_nominal × min(1, 0.5 + 0.5×(Re/1000)^0.3)
    pub(crate) fn effective_discharge_coefficient(&self, reynolds: T) -> T {
        let c_d_nom =
            T::from_f64(self.geometry_type.discharge_coefficient()).unwrap_or_else(T::one);
        let re_ref = T::from_f64(1000.0).expect("Mathematical constant conversion compromised");
        let ratio = reynolds / re_ref;

        // For Re >> 1000, correction → 1.0
        // For Re << 1000, correction reduces C_d
        let half = T::one() / (T::one() + T::one());
        let exp = T::from_f64(0.3).expect("Mathematical constant conversion compromised");
        let correction = half + half * ratio.powf(exp);
        let one = T::one();
        let correction_clamped = if correction > one { one } else { correction };

        c_d_nom * correction_clamped
    }
}

impl<T: RealField + Copy + FromPrimitive> ResistanceModel<T> for VenturiModel<T> {
    fn calculate_resistance<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<T> {
        let (r, k) = self.calculate_coefficients(fluid, conditions)?;

        // Effective resistance: R_eff = R + k|Q|
        let q_mag = if let Some(q) = conditions.flow_rate {
            if q >= T::zero() {
                q
            } else {
                -q
            }
        } else if let Some(v) = conditions.velocity {
            let v_abs = if v >= T::zero() { v } else { -v };
            v_abs * self.inlet_area()
        } else {
            T::zero()
        };

        Ok(r + k * q_mag)
    }

    fn calculate_coefficients<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<(T, T)> {
        let state = fluid.properties_at(conditions.temperature, conditions.pressure)?;
        let density = state.density;

        // --- Compute throat velocity and Reynolds number ---
        let a_inlet = self.inlet_area();
        let a_throat = self.throat_area();
        let a_outlet = self.outlet_area();

        let v_inlet = if let Some(v) = conditions.velocity {
            v
        } else if let Some(q) = conditions.flow_rate {
            q / a_inlet
        } else {
            T::zero()
        };

        let v_throat = v_inlet * a_inlet / a_throat;
        let v_outlet = v_inlet * a_inlet / a_outlet;

        // Compute shear rate at throat wall: γ̇ = 8V/D for circular pipe
        let eight = T::from_f64(8.0).expect("Mathematical constant conversion compromised");
        let shear_rate_throat = eight * v_throat / self.throat_diameter;

        // Get viscosity at throat shear rate (supports non-Newtonian fluids)
        let viscosity = fluid.viscosity_at_shear(
            shear_rate_throat,
            conditions.temperature,
            conditions.pressure,
        )?;

        // Reynolds number at throat
        let re_throat = density * v_throat * self.throat_diameter / viscosity;

        // Reynolds number at inlet (for C_d correction)
        let viscosity_inlet = fluid.viscosity_at_shear(
            eight * v_inlet / self.inlet_diameter,
            conditions.temperature,
            conditions.pressure,
        )?;
        let re_inlet = density * v_inlet * self.inlet_diameter / viscosity_inlet;

        // --- 1. Contraction loss (Bernoulli + discharge coefficient) ---
        let beta_sq = self.beta_squared(); // = (D_t/D_i)² = A_t/A_i
        let c_d = self.effective_discharge_coefficient(re_inlet);

        // ΔP_contraction = ½ρV_t²(1 − (A_t/A_i)²) / C_d²
        //                = ½ρV_t²(1 − β⁴) / C_d²   where β⁴ = beta_sq·beta_sq
        // Note: (A_t/A_i)² = (D_t/D_i)⁴ = beta_sq², NOT beta_sq.
        let half = T::one() / (T::one() + T::one());
        let one = T::one();
        let dp_contraction =
            half * density * v_throat * v_throat * (one - beta_sq * beta_sq) / (c_d * c_d);

        // --- 2. Throat friction loss ---
        let f = self.throat_friction_factor(re_throat);
        let dp_friction =
            f * (self.throat_length / self.throat_diameter) * half * density * v_throat * v_throat;

        // --- 3. Expansion loss (Borda-Carnot) ---
        let k_exp = T::from_f64(self.expansion_type.loss_coefficient()).unwrap_or_else(T::one);
        let dv = v_throat - v_outlet;
        let dp_expansion_loss = k_exp * half * density * dv * dv;

        // --- 4. Expansion recovery (ideal Bernoulli) ---
        let eta_r = T::from_f64(self.expansion_type.recovery_efficiency()).unwrap_or_else(T::one);
        let dp_recovery = eta_r * half * density * (v_throat * v_throat - v_outlet * v_outlet);

        // --- Net pressure drop ---
        let q = v_inlet * a_inlet;
        let q_sq = q * q;

        let vel_threshold =
            T::from_f64(1e-15).expect("Mathematical constant conversion compromised");
        if v_inlet > vel_threshold {
            // Decompose: laminar part (friction ∝ V → R·Q) and inertial part (∝ V² → k·Q²)
            let r = dp_friction / q; // Friction is approximately linear at low Re
            let k_coeff = if q_sq > T::zero() {
                (dp_contraction + dp_expansion_loss - dp_recovery) / q_sq
            } else {
                T::zero()
            };
            Ok((r, k_coeff))
        } else {
            // Zero flow: return linear estimate from viscous part only
            // R = 128 μ L_throat / (π D_throat⁴) (Hagen-Poiseuille in throat)
            let coeff = T::from_f64(128.0).expect("Mathematical constant conversion compromised");
            let pi = T::pi();
            let d2 = self.throat_diameter * self.throat_diameter;
            let d4 = d2 * d2;
            let r = coeff * viscosity * self.throat_length / (pi * d4);
            Ok((r, T::zero()))
        }
    }

    fn model_name(&self) -> &'static str {
        "Venturi"
    }

    fn reynolds_range(&self) -> (T, T) {
        // Valid across laminar and turbulent regimes
        (
            T::from_f64(0.1).expect("Mathematical constant conversion compromised"),
            T::from_f64(1e7).expect("Mathematical constant conversion compromised"),
        )
    }

    fn validate_invariants<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<()> {
        // Mach number validation
        self.validate_mach_number(fluid, conditions)?;

        // Positive dimensions
        if self.throat_length <= T::zero() || self.total_length <= T::zero() {
            return Err(Error::PhysicsViolation(
                "Venturi lengths must be positive".to_string(),
            ));
        }
        if self.inlet_diameter <= T::zero()
            || self.throat_diameter <= T::zero()
            || self.outlet_diameter <= T::zero()
        {
            return Err(Error::PhysicsViolation(
                "Venturi diameters must be positive".to_string(),
            ));
        }

        // Throat must be narrower than inlet and outlet
        if self.throat_diameter >= self.inlet_diameter {
            return Err(Error::PhysicsViolation(
                "Venturi throat diameter must be less than inlet diameter".to_string(),
            ));
        }

        // Area ratio constraint: β = D_throat/D_inlet typically 0.2-0.75
        let beta = self.throat_diameter / self.inlet_diameter;
        let beta_min = T::from_f64(0.1).expect("Mathematical constant conversion compromised");
        let beta_max = T::from_f64(0.9).expect("Mathematical constant conversion compromised");
        if beta < beta_min || beta > beta_max {
            return Err(Error::PhysicsViolation(format!(
                "Venturi β = D_throat/D_inlet = {beta:.3} outside recommended range [0.1, 0.9]"
            )));
        }

        Ok(())
    }
}
