//! Venturi tube resistance model for converging-diverging channels.
//!
//! ## Venturi Effect
//!
//! **Theorem**: In a converging-diverging (Venturi) tube, the pressure drop
//! between the upstream section (1) and the throat section (2) is given by:
//!
//! ΔP₁₂ = (ρ V₂²/2)(1 − (A₂/A₁)²) / C_d²
//!
//! where (A₂/A₁) is the **area ratio** (throat / inlet).
//! In terms of the diameter ratio β = D₂/D₁:
//!   (A₂/A₁) = β²   so   (A₂/A₁)² = β⁴
//! Therefore: ΔP₁₂ = (ρ V₂²/2)(1 − β⁴) / C_d²
//!
//! where:
//! - ρ is the fluid density [kg/m³]
//! - V₂ is the throat velocity [m/s]
//! - A₁ is the upstream cross-sectional area [m²]
//! - A₂ is the throat cross-sectional area [m²]
//! - C_d is the discharge coefficient (accounts for viscous losses)
//!
//! ### Derivation (Bernoulli + Continuity)
//!
//! Starting from Bernoulli's equation along a streamline:
//!
//! P₁ + ½ρV₁² = P₂ + ½ρV₂²  (ideal, inviscid)
//!
//! Continuity: A₁V₁ = A₂V₂  →  V₁ = V₂(A₂/A₁)
//!
//! Substituting:
//! P₁ - P₂ = ½ρV₂²[1 - (A₂/A₁)²]
//!
//! With discharge coefficient: ΔP = ΔP_ideal / C_d²
//!
//! ### Pressure Recovery (Expansion)
//!
//! The recovery in the diverging section is given by:
//!
//! ΔP₂₃ = η_r · ½ρ(V₂² - V₃²)
//!
//! where η_r is the recovery efficiency (typically 0.6-0.9 depending on
//! diffuser half-angle and area ratio).
//!
//! ### Borda-Carnot Loss
//!
//! For sudden or gradual expansions, the irreversible pressure loss is:
//!
//! ΔP_loss = K_exp · ½ρ(V₂ - V₃)²
//!
//! where K_exp depends on the expansion geometry:
//! - Sudden expansion: K_exp = 1.0 (Borda-Carnot)
//! - Gradual (7° half-angle): K_exp ≈ 0.2
//! - Optimal diffuser (3-4°): K_exp ≈ 0.1
//!
//! ### ISO 5167 Discharge Coefficients
//!
//! For classical Venturi tubes (ISO 5167-4):
//! - Machined convergent: C_d = 0.995 for 2×10⁵ < Re_D < 2×10⁶
//! - Rough-cast convergent: C_d = 0.984 for 2×10⁵ < Re_D < 2×10⁶
//! - Rough-welded convergent: C_d = 0.985
//!
//! For millifluidic Venturi tubes (Re < 2300):
//! C_d depends strongly on Reynolds number and is typically 0.90-0.99.
//!
//! ### Millifluidic/Microfluidic Considerations
//!
//! At low Reynolds numbers (Re < 100), viscous effects dominate:
//! - C_d is significantly less than 1.0
//! - The entrance length in the throat affects the total pressure drop
//! - Non-Newtonian effects (blood) modify the velocity profile
//! - Wall shear rate γ̇ = 8V/D for Newtonian, higher for shear-thinning
//!
//! ### Throat Friction Contribution
//!
//! Additional pressure loss from viscous friction in the throat:
//!
//! ΔP_friction = f · (L_t/D_t) · (ρ V_t²/2)
//!
//! where:
//! - f is the Darcy friction factor in the throat
//! - L_t is the throat length [m]
//! - D_t is the throat diameter [m]
//!
//! ### References
//!
//! - Venturi, G. B. (1797). *Recherches expérimentales sur le principe de la
//!   communication latérale du mouvement dans les fluides.*
//! - ISO 5167-4:2003. "Measurement of fluid flow by means of pressure
//!   differential devices — Part 4: Venturi tubes."
//! - Reader-Harris, M. (2015). *Orifice Plates and Venturi Tubes.* Springer.
//! - Nabavi, M. (2009). "Steady and unsteady flow analysis in microdiffusers
//!   and micropumps: A critical review." *Microfluidics and Nanofluidics*, 7, 599-619.
//! - Idelchik, I. E. (2007). *Handbook of Hydraulic Resistance* (4th ed.).
//!   Begell House. Table 5-1 through 5-8.

use super::traits::{FlowConditions, ResistanceModel};
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid::FluidTrait;
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;
use serde::{Deserialize, Serialize};

// Named constants
const LAMINAR_LIMIT_RE: f64 = 2300.0;
const LAMINAR_FRICTION_COEFF: f64 = 64.0;
const BLASIUS_COEFF: f64 = 0.3164;
const BLASIUS_EXP: f64 = 0.25;

/// Venturi tube geometry specification
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum VenturiGeometry {
    /// Machined convergent section (ISO 5167-4, C_d ≈ 0.995)
    MachinedConvergent,
    /// Rough-cast convergent section (ISO 5167-4, C_d ≈ 0.984)
    RoughCastConvergent,
    /// Rough-welded convergent section (ISO 5167-4, C_d ≈ 0.985)
    RoughWeldedConvergent,
    /// Custom discharge coefficient
    Custom {
        /// Discharge coefficient C_d
        discharge_coefficient: f64,
    },
}

impl VenturiGeometry {
    /// Get the nominal discharge coefficient for this geometry type
    #[must_use]
    pub fn discharge_coefficient(&self) -> f64 {
        match self {
            Self::MachinedConvergent => 0.995,
            Self::RoughCastConvergent => 0.984,
            Self::RoughWeldedConvergent => 0.985,
            Self::Custom { discharge_coefficient } => *discharge_coefficient,
        }
    }
}

/// Expansion geometry type for the diverging section
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum ExpansionType {
    /// Sudden expansion (Borda-Carnot, K_exp = 1.0)
    Sudden,
    /// Gradual expansion with specified half-angle [degrees]
    Gradual {
        /// Diffuser half-angle in degrees
        half_angle_deg: f64,
    },
}

impl ExpansionType {
    /// Get the expansion loss coefficient K_exp
    ///
    /// Based on Idelchik (2007) Table 5-1 for conical diffusers
    #[must_use]
    pub fn loss_coefficient(&self) -> f64 {
        match self {
            Self::Sudden => 1.0,
            Self::Gradual { half_angle_deg } => {
                let angle = *half_angle_deg;
                if angle <= 3.0 {
                    0.10 // Optimal diffuser
                } else if angle <= 5.0 {
                    0.14
                } else if angle <= 7.0 {
                    0.20
                } else if angle <= 10.0 {
                    0.30
                } else if angle <= 15.0 {
                    0.45
                } else if angle <= 20.0 {
                    0.55
                } else if angle <= 30.0 {
                    0.70
                } else if angle <= 45.0 {
                    0.85
                } else {
                    1.0 // Approaches sudden expansion
                }
            }
        }
    }

    /// Get pressure recovery efficiency η_r = 1 - K_exp
    #[must_use]
    pub fn recovery_efficiency(&self) -> f64 {
        1.0 - self.loss_coefficient()
    }
}

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
            expansion_type: ExpansionType::Gradual { half_angle_deg: 7.0 },
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
        Self::new(inlet_diameter, throat_diameter, inlet_diameter, throat_length, total_length)
    }

    /// Create a millifluidic Venturi with typical parameters
    pub fn millifluidic(
        inlet_diameter: T,
        throat_diameter: T,
        throat_length: T,
    ) -> Self {
        let total_length = throat_length * T::from_f64(5.0).unwrap_or_else(T::one);
        Self {
            inlet_diameter,
            throat_diameter,
            outlet_diameter: inlet_diameter,
            throat_length,
            total_length,
            geometry_type: VenturiGeometry::Custom {
                discharge_coefficient: 0.97,
            },
            expansion_type: ExpansionType::Gradual { half_angle_deg: 5.0 },
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
    fn beta_squared(&self) -> T {
        let ratio = self.throat_diameter / self.inlet_diameter;
        ratio * ratio
    }

    /// Inlet cross-sectional area [m²]
    fn inlet_area(&self) -> T {
        let pi = T::from_f64(std::f64::consts::PI).unwrap_or_else(T::zero);
        pi * self.inlet_diameter * self.inlet_diameter / T::from_f64(4.0).unwrap_or_else(T::one)
    }

    /// Throat cross-sectional area [m²]
    fn throat_area(&self) -> T {
        let pi = T::from_f64(std::f64::consts::PI).unwrap_or_else(T::zero);
        pi * self.throat_diameter * self.throat_diameter / T::from_f64(4.0).unwrap_or_else(T::one)
    }

    /// Outlet cross-sectional area [m²]
    fn outlet_area(&self) -> T {
        let pi = T::from_f64(std::f64::consts::PI).unwrap_or_else(T::zero);
        pi * self.outlet_diameter * self.outlet_diameter / T::from_f64(4.0).unwrap_or_else(T::one)
    }

    /// Calculate the Darcy friction factor in the throat
    ///
    /// Uses Hagen-Poiseuille (laminar) or Blasius (turbulent) correlations
    fn throat_friction_factor(&self, reynolds_throat: T) -> T {
        let re_lam = T::from_f64(LAMINAR_LIMIT_RE).unwrap_or_else(T::one);

        if reynolds_throat < re_lam {
            // Laminar: f = 64/Re
            T::from_f64(LAMINAR_FRICTION_COEFF).unwrap_or_else(T::one) / reynolds_throat
        } else {
            // Blasius: f = 0.3164 / Re^0.25
            let coeff = T::from_f64(BLASIUS_COEFF).unwrap_or_else(T::one);
            let exp = T::from_f64(BLASIUS_EXP).unwrap_or_else(T::one);
            coeff / reynolds_throat.powf(exp)
        }
    }

    /// Effective discharge coefficient adjusted for Reynolds number
    ///
    /// At low Re (millifluidic regime), C_d decreases due to boundary layer growth.
    /// Empirical correction: C_d_eff = C_d_nominal × min(1, 0.5 + 0.5×(Re/1000)^0.3)
    fn effective_discharge_coefficient(&self, reynolds: T) -> T {
        let c_d_nom = T::from_f64(self.geometry_type.discharge_coefficient())
            .unwrap_or_else(T::one);
        let re_ref = T::from_f64(1000.0).unwrap_or_else(T::one);
        let ratio = reynolds / re_ref;

        // For Re >> 1000, correction → 1.0
        // For Re << 1000, correction reduces C_d
        let half = T::from_f64(0.5).unwrap_or_else(T::one);
        let exp = T::from_f64(0.3).unwrap_or_else(T::one);
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
            if q >= T::zero() { q } else { -q }
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
        let eight = T::from_f64(8.0).unwrap_or_else(T::one);
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
        let beta_sq = self.beta_squared();  // = (D_t/D_i)² = A_t/A_i
        let c_d = self.effective_discharge_coefficient(re_inlet);

        // ΔP_contraction = ½ρV_t²(1 − (A_t/A_i)²) / C_d²
        //                = ½ρV_t²(1 − β⁴) / C_d²   where β⁴ = beta_sq·beta_sq
        // Note: (A_t/A_i)² = (D_t/D_i)⁴ = beta_sq², NOT beta_sq.
        let half = T::from_f64(0.5).unwrap_or_else(T::one);
        let one = T::one();
        let dp_contraction =
            half * density * v_throat * v_throat * (one - beta_sq * beta_sq) / (c_d * c_d);

        // --- 2. Throat friction loss ---
        let f = self.throat_friction_factor(re_throat);
        let dp_friction = f * (self.throat_length / self.throat_diameter)
            * half * density * v_throat * v_throat;

        // --- 3. Expansion loss (Borda-Carnot) ---
        let k_exp = T::from_f64(self.expansion_type.loss_coefficient())
            .unwrap_or_else(T::one);
        let dv = v_throat - v_outlet;
        let dp_expansion_loss = k_exp * half * density * dv * dv;

        // --- 4. Expansion recovery (ideal Bernoulli) ---
        let eta_r = T::from_f64(self.expansion_type.recovery_efficiency())
            .unwrap_or_else(T::one);
        let dp_recovery = eta_r * half * density * (v_throat * v_throat - v_outlet * v_outlet);

        // --- Net pressure drop ---
        // Total ΔP = contraction + friction + expansion_loss - recovery
        // Note: recovery is pressure regain, so it reduces the total drop
        // Convert to resistance form: ΔP = R·Q + k·Q|Q|
        // For the quadratic component (dominant at higher Re):
        // ΔP ∝ V² ∝ Q² → k = ΔP / Q²
        let q = v_inlet * a_inlet;
        let q_sq = q * q;

        // Use velocity-based check instead of q_sq > epsilon, because
        // microfluidic flow rates (Q ~ 1e-10 m³/s) yield q² ~ 1e-20 which
        // is below f64::EPSILON ≈ 2.2e-16 even though the flow is physically
        // meaningful. Only fall back to analytical limit for truly zero flow.
        let vel_threshold = T::from_f64(1e-15).unwrap_or_else(T::zero);
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
            let coeff = T::from_f64(128.0).unwrap_or_else(T::one);
            let pi = T::from_f64(std::f64::consts::PI).unwrap_or_else(T::one);
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
            T::from_f64(0.1).unwrap_or_else(T::zero),
            T::from_f64(1e7).unwrap_or_else(T::zero),
        )
    }

    fn validate_invariants<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<()> {
        // Mach number validation
        self.validate_mach_number(fluid, conditions)?;

        // Throat must be narrower than inlet and outlet
        if self.throat_diameter >= self.inlet_diameter {
            return Err(Error::PhysicsViolation(
                "Venturi throat diameter must be less than inlet diameter".to_string(),
            ));
        }

        // Area ratio constraint: β = D_throat/D_inlet typically 0.2-0.75
        let beta = self.throat_diameter / self.inlet_diameter;
        let beta_min = T::from_f64(0.1).unwrap_or_else(T::zero);
        let beta_max = T::from_f64(0.9).unwrap_or_else(T::one);
        if beta < beta_min || beta > beta_max {
            return Err(Error::PhysicsViolation(format!(
                "Venturi β = D_throat/D_inlet = {:.3} outside recommended range [0.1, 0.9]",
                beta
            )));
        }

        // Positive dimensions
        if self.throat_length <= T::zero() || self.total_length <= T::zero() {
            return Err(Error::PhysicsViolation(
                "Venturi lengths must be positive".to_string(),
            ));
        }

        Ok(())
    }
}

/// Detailed Venturi flow analysis result
#[derive(Debug, Clone)]
pub struct VenturiAnalysis<T: RealField + Copy> {
    /// Throat velocity [m/s]
    pub throat_velocity: T,
    /// Throat Reynolds number
    pub throat_reynolds: T,
    /// Throat wall shear rate [1/s]
    pub throat_shear_rate: T,
    /// Apparent viscosity at throat [Pa·s]
    pub throat_viscosity: T,
    /// Contraction pressure drop [Pa]
    pub dp_contraction: T,
    /// Throat friction pressure drop [Pa]
    pub dp_friction: T,
    /// Expansion loss [Pa]
    pub dp_expansion_loss: T,
    /// Expansion recovery [Pa]
    pub dp_recovery: T,
    /// Net pressure drop [Pa]
    pub dp_total: T,
    /// Discharge coefficient used
    pub discharge_coefficient: T,
    /// Expansion loss coefficient
    pub expansion_loss_coefficient: T,
    /// Darcy friction factor in throat
    pub friction_factor: T,
}

impl<T: RealField + Copy + FromPrimitive> VenturiModel<T> {
    /// Perform detailed Venturi flow analysis
    ///
    /// Returns comprehensive breakdown of pressure contributions
    pub fn analyze<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<VenturiAnalysis<T>> {
        let state = fluid.properties_at(conditions.temperature, conditions.pressure)?;
        let density = state.density;

        let a_inlet = self.inlet_area();
        let a_throat = self.throat_area();
        let a_outlet = self.outlet_area();

        let v_inlet = if let Some(v) = conditions.velocity {
            v
        } else if let Some(q) = conditions.flow_rate {
            q / a_inlet
        } else {
            return Err(Error::InvalidConfiguration(
                "Venturi analysis requires velocity or flow_rate".to_string(),
            ));
        };

        let v_throat = v_inlet * a_inlet / a_throat;
        let v_outlet = v_inlet * a_inlet / a_outlet;

        let eight = T::from_f64(8.0).unwrap_or_else(T::one);
        let half = T::from_f64(0.5).unwrap_or_else(T::one);
        let one = T::one();

        let shear_rate_throat = eight * v_throat / self.throat_diameter;
        let viscosity = fluid.viscosity_at_shear(
            shear_rate_throat,
            conditions.temperature,
            conditions.pressure,
        )?;

        let re_throat = density * v_throat * self.throat_diameter / viscosity;

        let viscosity_inlet = fluid.viscosity_at_shear(
            eight * v_inlet / self.inlet_diameter,
            conditions.temperature,
            conditions.pressure,
        )?;
        let re_inlet = density * v_inlet * self.inlet_diameter / viscosity_inlet;

        let beta_sq = self.beta_squared();
        let c_d = self.effective_discharge_coefficient(re_inlet);
        let f = self.throat_friction_factor(re_throat);
        let k_exp = T::from_f64(self.expansion_type.loss_coefficient()).unwrap_or_else(T::one);
        let eta_r = T::from_f64(self.expansion_type.recovery_efficiency()).unwrap_or_else(T::one);

        // ΔP_contraction = ½ρV_t²(1 − β⁴) / C_d²  where β⁴ = (A_t/A_i)² = beta_sq·beta_sq
        let dp_contraction =
            half * density * v_throat * v_throat * (one - beta_sq * beta_sq) / (c_d * c_d);
        let dp_friction = f * (self.throat_length / self.throat_diameter)
            * half * density * v_throat * v_throat;
        let dv = v_throat - v_outlet;
        let dp_expansion_loss = k_exp * half * density * dv * dv;
        let dp_recovery = eta_r * half * density * (v_throat * v_throat - v_outlet * v_outlet);
        let dp_total = dp_contraction + dp_friction + dp_expansion_loss - dp_recovery;

        Ok(VenturiAnalysis {
            throat_velocity: v_throat,
            throat_reynolds: re_throat,
            throat_shear_rate: shear_rate_throat,
            throat_viscosity: viscosity,
            dp_contraction,
            dp_friction,
            dp_expansion_loss,
            dp_recovery,
            dp_total,
            discharge_coefficient: c_d,
            expansion_loss_coefficient: k_exp,
            friction_factor: f,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_venturi_bernoulli_limit() -> Result<()> {
        // At high Re, the pressure drop should approach the ideal Bernoulli value
        // ΔP = ½ρV₂²(1 - β²)
        let model = VenturiModel::symmetric(
            0.01,   // 10mm inlet
            0.005,  // 5mm throat (β = 0.5)
            0.01,   // 10mm throat length
            0.05,   // 50mm total
        ).with_geometry(VenturiGeometry::Custom { discharge_coefficient: 1.0 })
         .with_expansion(ExpansionType::Gradual { half_angle_deg: 3.0 });

        let fluid = cfd_core::physics::fluid::database::water_20c::<f64>()?;

        // High velocity (turbulent Re)
        let velocity = 2.0; // m/s
        let mut conditions = FlowConditions::new(velocity);
        conditions.reynolds_number = Some(998.2 * velocity * 0.01 / 1.002e-3);

        let analysis = model.analyze(&fluid, &conditions)?;

        // Check throat velocity: V_throat = V_inlet * (D_inlet/D_throat)² = 2 * 4 = 8 m/s
        assert_relative_eq!(analysis.throat_velocity, 8.0, epsilon = 0.01);

        // β = D_t/D_i = 0.5  →  β⁴ = 0.0625  →  (1 − β⁴) = 0.9375
        // Ideal Bernoulli: ΔP = ½ * 998.2 * 8² * (1 − β⁴) = ½ * 998.2 * 64 * 0.9375 ≈ 29,946 Pa
        // With C_d = 1.0, the contraction drop must equal the ideal Bernoulli value.
        let dp_ideal = 0.5 * 998.2 * 64.0 * 0.9375;
        assert_relative_eq!(analysis.dp_contraction, dp_ideal, epsilon = dp_ideal * 0.05);

        Ok(())
    }

    #[test]
    fn test_venturi_millifluidic_blood() -> Result<()> {
        // Test with blood properties at millifluidic scale
        let model = VenturiModel::millifluidic(
            0.001,   // 1mm inlet
            0.0005,  // 0.5mm throat
            0.002,   // 2mm throat length
        );

        // Create a simple Newtonian approximation of blood
        // (real blood models are non-Newtonian but this tests the flow path)
        let fluid = cfd_core::physics::fluid::database::water_20c::<f64>()?;

        let mut conditions = FlowConditions::new(0.01); // 10 mm/s
        let density = 998.2;
        let viscosity = 1.002e-3;
        let re = density * 0.01 * 0.001 / viscosity;
        conditions.reynolds_number = Some(re);

        let analysis = model.analyze(&fluid, &conditions)?;

        // At low Re, friction should be significant relative to Bernoulli term
        assert!(analysis.dp_friction > 0.0);
        assert!(analysis.dp_contraction > 0.0);
        assert!(analysis.dp_total > 0.0);

        // Verify throat Reynolds number
        // V_throat = V_inlet * (A_inlet/A_throat) = 0.01 * 4 = 0.04 m/s
        assert_relative_eq!(analysis.throat_velocity, 0.04, epsilon = 0.001);

        Ok(())
    }

    #[test]
    fn test_venturi_area_ratio() {
        let model = VenturiModel::<f64>::symmetric(0.01, 0.005, 0.01, 0.05);

        // β² = (D_t/D_i)² = (0.005/0.01)² = 0.25
        assert_relative_eq!(model.beta_squared(), 0.25, epsilon = 1e-10);

        // Areas
        let a_inlet = std::f64::consts::PI * 0.01_f64.powi(2) / 4.0;
        assert_relative_eq!(model.inlet_area(), a_inlet, epsilon = 1e-10);

        let a_throat = std::f64::consts::PI * 0.005_f64.powi(2) / 4.0;
        assert_relative_eq!(model.throat_area(), a_throat, epsilon = 1e-10);
    }

    #[test]
    fn test_venturi_friction_factor() {
        let model = VenturiModel::<f64>::symmetric(0.01, 0.005, 0.01, 0.05);

        // Laminar: f = 64/Re
        let f_lam = model.throat_friction_factor(100.0);
        assert_relative_eq!(f_lam, 0.64, epsilon = 1e-6);

        // Blasius: f = 0.3164 / Re^0.25
        let f_turb = model.throat_friction_factor(10000.0);
        let expected = 0.3164 / 10000.0_f64.powf(0.25);
        assert_relative_eq!(f_turb, expected, epsilon = 1e-6);
    }

    #[test]
    fn test_venturi_expansion_coefficients() {
        assert_relative_eq!(ExpansionType::Sudden.loss_coefficient(), 1.0, epsilon = 1e-10);
        assert_relative_eq!(
            ExpansionType::Gradual { half_angle_deg: 3.0 }.loss_coefficient(),
            0.10,
            epsilon = 1e-10
        );
        assert_relative_eq!(
            ExpansionType::Gradual { half_angle_deg: 7.0 }.loss_coefficient(),
            0.20,
            epsilon = 1e-10
        );
    }

    #[test]
    fn test_venturi_validate_invariants() {
        let model = VenturiModel::<f64>::symmetric(0.01, 0.005, 0.01, 0.05);
        let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().unwrap();
        let conditions = FlowConditions::new(0.1);

        // Should validate OK
        assert!(model.validate_invariants(&fluid, &conditions).is_ok());

        // Invalid: throat >= inlet
        let bad_model = VenturiModel::<f64>::symmetric(0.005, 0.01, 0.01, 0.05);
        assert!(bad_model.validate_invariants(&fluid, &conditions).is_err());
    }

    #[test]
    fn test_venturi_resistance_model_trait() -> Result<()> {
        let model = VenturiModel::symmetric(0.01_f64, 0.005, 0.01, 0.05);
        let fluid = cfd_core::physics::fluid::database::water_20c::<f64>()?;

        let mut conditions = FlowConditions::new(0.1);
        conditions.reynolds_number = Some(500.0);

        // Test that calculate_resistance returns positive value
        let resistance = model.calculate_resistance(&fluid, &conditions)?;
        assert!(resistance > 0.0);

        // Test coefficients
        let (r, k) = model.calculate_coefficients(&fluid, &conditions)?;
        assert!(r >= 0.0 || k >= 0.0); // At least one should be positive

        Ok(())
    }

    #[test]
    fn test_iso5167_discharge_coefficients() {
        assert_relative_eq!(
            VenturiGeometry::MachinedConvergent.discharge_coefficient(),
            0.995,
            epsilon = 1e-10
        );
        assert_relative_eq!(
            VenturiGeometry::RoughCastConvergent.discharge_coefficient(),
            0.984,
            epsilon = 1e-10
        );
        assert_relative_eq!(
            VenturiGeometry::RoughWeldedConvergent.discharge_coefficient(),
            0.985,
            epsilon = 1e-10
        );
    }
}
