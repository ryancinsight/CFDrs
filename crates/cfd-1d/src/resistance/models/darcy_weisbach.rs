//! Darcy-Weisbach resistance model for turbulent flow.
//!
//! ## Darcy-Weisbach Theorem
//!
//! **Theorem**: For steady, fully developed, incompressible flow in a straight pipe,
//! the pressure drop ΔP over a length L is given by:
//!
//! ΔP = f * (L/D) * (ρ V²/2)
//!
//! where:
//! - f is the Darcy friction factor (dimensionless)
//! - L is the pipe length [m]
//! - D is the pipe diameter [m]
//! - ρ is the fluid density [kg/m³]
//! - V is the average velocity [m/s]
//!
//! **Hydraulic Resistance Form**: R = ΔP / Q = (f ρ L) / (2 A D)
//!
//! where:
//! - Q is the volumetric flow rate [m³/s]
//! - A is the pipe cross-sectional area [m²]
//!
//! ### Derivation
//!
//! Starting from momentum balance on a pipe element:
//!
//! τ_w * π D L = ΔP * π D²/4
//!
//! where τ_w is the wall shear stress.
//!
//! Wall shear stress: τ_w = (f ρ V²)/8
//!
//! Substituting: (f ρ V²/8) * π D L = ΔP * π D²/4
//!
//! Solving for ΔP: ΔP = f * (L/D) * (ρ V²/2)
//!
//! ### Colebrook-White Theorem
//!
//! **Theorem**: The Darcy friction factor f for rough pipes is given by the
//! implicit equation:
//!
//! 1/√f = -2.0 log₁₀(ε/(3.7 D) + 2.51/(Re √f))
//!
//! where:
//! - ε is the surface roughness [m]
//! - Re is the Reynolds number (ρ V D / μ)
//!
//! **Validity Conditions**:
//! - Turbulent flow: Re > 4000 (typically Re > 2300 for safety)
//! - Roughness ratio: 0 < ε/D < 0.05
//! - Fully developed flow: L/D > 10
//!
//! ### Special Cases
//!
//! **Laminar Flow (Blasius)**: f = 64/Re for Re < 2300
//!
//! **Smooth Pipes (Prandtl-von Kármán)**: 1/√f ≈ 2.0 log₁₀(Re √f) - 0.8
//!
//! **Fully Rough Pipes (Nikuradse)**: 1/√f ≈ 2.0 log₁₀(3.7 D/ε)
//!
//! ### Haaland Explicit Approximation
//!
//! For practical calculations, Haaland (1983) provides an explicit approximation:
//!
//! 1/√f ≈ -1.8 log₁₀[(ε/(3.7 D))¹.¹ + (6.9/Re)¹.¹]¹/¹.¹
//!
//! **Accuracy**: Within 2% of Colebrook-White for 4000 < Re < 10⁸ and 10⁻⁶ < ε/D < 10⁻²
//!
//! ### References
//!
//! - Darcy, H. (1857). "Recherches expérimentales relatives au mouvement de l'eau
//!   dans les tuyaux." *Mémoires de l'Académie des Sciences*, 15, 141-168.
//! - Weisbach, J. (1845). *Lehrbuch der Ingenieur- und Maschinen-Mechanik*.
//! - Colebrook, C. F. (1939). "Turbulent flow in pipes with particular reference
//!   to the transition region between smooth and rough pipe laws."
//!   *Journal of the Institution of Civil Engineers*, 11(4), 133-156.
//! - Moody, L. F. (1944). "Friction factors for pipe flow."
//!   *Transactions of the ASME*, 66(8), 671-684.
//! - Haaland, S. E. (1983). "Simple and explicit formulas for the friction factor
//!   in turbulent pipe flow." *Journal of Fluids Engineering*, 105(1), 89-90.

use super::traits::{FlowConditions, ResistanceModel};
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid::FluidTrait;
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;
use serde::{Deserialize, Serialize};

// Named constants for friction factor calculations
const LAMINAR_FRICTION_COEFFICIENT: f64 = 64.0;
const LAMINAR_TRANSITION_RE: f64 = 2300.0;
const HAALAND_ROUGHNESS_DIVISOR: f64 = 3.6;
const HAALAND_REYNOLDS_FACTOR: f64 = 6.9;
const LOG_BASE_10: f64 = 10.0;
const HAALAND_EXPONENT_FACTOR: f64 = 1.8;
const COLEBROOK_COEFFICIENT: f64 = 2.0;
const MAX_REYNOLDS: f64 = 1e8;

/// Darcy-Weisbach resistance model for turbulent flow
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DarcyWeisbachModel<T: RealField + Copy> {
    /// Hydraulic diameter [m]
    pub hydraulic_diameter: T,
    /// Channel cross-sectional area [m²]
    pub area: T,
    /// Channel length [m]
    pub length: T,
    /// Surface roughness [m]
    pub roughness: T,
}

impl<T: RealField + Copy + FromPrimitive> DarcyWeisbachModel<T> {
    /// Create a new Darcy-Weisbach model
    pub fn new(hydraulic_diameter: T, area: T, length: T, roughness: T) -> Self {
        Self {
            hydraulic_diameter,
            area,
            length,
            roughness,
        }
    }

    /// Create a new Darcy-Weisbach model for a circular pipe
    pub fn circular(diameter: T, length: T, roughness: T) -> Self {
        let pi = T::from_f64(std::f64::consts::PI).unwrap_or_else(|| T::zero());
        let area = pi * diameter * diameter / T::from_f64(4.0).unwrap_or_else(|| T::zero());
        Self {
            hydraulic_diameter: diameter,
            area,
            length,
            roughness,
        }
    }
}

impl<T: RealField + Copy + FromPrimitive> ResistanceModel<T> for DarcyWeisbachModel<T> {
    fn calculate_resistance<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<T> {
        let (r, k) = self.calculate_coefficients(fluid, conditions)?;

        let k_mag = if k >= T::zero() { k } else { -k };
        if k_mag <= T::default_epsilon() {
            return Ok(r);
        }

        let q_mag = if let Some(q) = conditions.flow_rate {
            if q >= T::zero() {
                q
            } else {
                -q
            }
        } else if let Some(v) = conditions.velocity {
            let v_abs = if v >= T::zero() { v } else { -v };
            v_abs * self.area
        } else {
            return Err(Error::InvalidConfiguration(
                "Flow rate or velocity required for turbulent resistance evaluation".to_string(),
            ));
        };

        Ok(r + k * q_mag)
    }

    fn calculate_coefficients<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<(T, T)> {
        let reynolds = conditions.reynolds_number.ok_or_else(|| {
            Error::InvalidConfiguration(
                "Reynolds number required for Darcy-Weisbach model".to_string(),
            )
        })?;

        let friction_factor = self.calculate_friction_factor(reynolds);

        let props = fluid.properties_at(conditions.temperature, conditions.pressure)?;
        let density = props.density;

        // Shear rate for viscosity
        let shear_rate = if let Some(sr) = conditions.shear_rate {
            sr
        } else {
            // Estimate 8*v/D
            if let Some(v) = conditions.velocity {
                T::from_f64(8.0).unwrap_or_else(T::one) * v.abs() / self.hydraulic_diameter
            } else {
                T::zero()
            }
        };

        let viscosity = fluid.viscosity_at_shear(shear_rate, conditions.temperature, conditions.pressure)?;

        let area = self.area;

        let re_transition = T::from_f64(LAMINAR_TRANSITION_RE).unwrap_or_else(|| T::one());

        if reynolds < re_transition {
            // Laminar regime: ΔP = R·Q
            // R = (f·ρ·L·V) / (2·A·Dh)
            // With f = 64/Re = 64μ / (ρ·V·Dh)
            // R = (64μ / (ρ·V·Dh)) · ρ·L·V / (2·A·Dh) = 32μL / (A·Dh²)
            let r = (T::from_f64(32.0).unwrap_or_else(|| T::zero()) * viscosity * self.length)
                / (area * self.hydraulic_diameter * self.hydraulic_diameter);
            Ok((r, T::zero()))
        } else {
            // Turbulent regime: ΔP = k·Q|Q|
            // k = (f·ρ·L) / (2·A²·Dh)
            let k = (friction_factor * density * self.length)
                / (T::from_f64(2.0).unwrap_or_else(|| T::zero())
                    * area
                    * area
                    * self.hydraulic_diameter);
            Ok((T::zero(), k))
        }
    }

    fn model_name(&self) -> &'static str {
        "Darcy-Weisbach"
    }

    fn reynolds_range(&self) -> (T, T) {
        (
            T::from_f64(
                cfd_core::physics::constants::physics::dimensionless::reynolds::PIPE_LAMINAR_MAX,
            )
            .unwrap_or_else(|| T::zero()),
            T::from_f64(MAX_REYNOLDS).unwrap_or_else(|| T::zero()),
        )
    }

    fn validate_invariants<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<()> {
        // Call Mach number validation
        self.validate_mach_number(fluid, conditions)?;

        // Entrance length validation: L/Dh > 10
        let ratio = self.length / self.hydraulic_diameter;
        let limit = T::from_f64(10.0).unwrap_or_else(|| T::zero());

        if ratio < limit {
            return Err(cfd_core::error::Error::PhysicsViolation(format!(
                "Entrance length violation: L/Dh = {:.2} < 10. Flow may not be fully developed for model '{}'",
                ratio,
                self.model_name()
            )));
        }

        // Roughness ratio validation: ε/Dh < 0.05
        let roughness_ratio = self.roughness / self.hydraulic_diameter;
        let roughness_limit = T::from_f64(0.05).unwrap_or_else(|| T::zero());
        if roughness_ratio > roughness_limit {
            return Err(cfd_core::error::Error::PhysicsViolation(format!(
                "Roughness ratio violation: ε/Dh = {:.4} > 0.05. Darcy-Weisbach model '{}' may be inaccurate",
                roughness_ratio,
                self.model_name()
            )));
        }

        Ok(())
    }
}

impl<T: RealField + Copy + FromPrimitive> DarcyWeisbachModel<T> {
    /// Calculate friction factor using iterative Colebrook-White equation
    fn calculate_friction_factor(&self, reynolds: T) -> T {
        use crate::components::constants::COLEBROOK_TOLERANCE;
        use cfd_core::physics::constants::physics::hydraulics::{
            COLEBROOK_REYNOLDS_NUMERATOR, COLEBROOK_ROUGHNESS_DIVISOR,
        };

        let relative_roughness = self.roughness / self.hydraulic_diameter;

        // Check for laminar flow
        let re_transition = T::from_f64(LAMINAR_TRANSITION_RE).unwrap_or_else(|| T::one());
        if reynolds < re_transition {
            // Laminar flow: f = 64/Re
            return T::from_f64(LAMINAR_FRICTION_COEFFICIENT).unwrap_or_else(|| T::one())
                / reynolds;
        }

        // Initial guess using Haaland explicit formula for convergence
        let mut f = {
            let term = relative_roughness
                / T::from_f64(HAALAND_ROUGHNESS_DIVISOR).unwrap_or_else(|| T::one())
                + T::from_f64(HAALAND_REYNOLDS_FACTOR).unwrap_or_else(|| T::one()) / reynolds;
            let log_term = term.ln() / T::from_f64(LOG_BASE_10.ln()).unwrap_or_else(|| T::one());
            T::one()
                / (T::from_f64(HAALAND_EXPONENT_FACTOR).unwrap_or_else(|| T::one()) * log_term)
                    .powi(2)
        };

        // Iterative solution of Colebrook-White equation
        // 1/sqrt(f) = -2.0 * log10(ε/(3.7*D) + 2.51/(Re*sqrt(f)))
        let tolerance = T::from_f64(COLEBROOK_TOLERANCE)
            .unwrap_or_else(|| T::from_f64(1e-6).unwrap_or_else(|| T::zero()));
        let max_iter = 50;
        let ln10 = T::from_f64(LOG_BASE_10.ln()).unwrap_or_else(|| T::one());

        for _ in 0..max_iter {
            let sqrt_f = f.sqrt();
            let term1 = relative_roughness
                / T::from_f64(COLEBROOK_ROUGHNESS_DIVISOR).unwrap_or_else(|| T::one());
            let term2 = T::from_f64(COLEBROOK_REYNOLDS_NUMERATOR).unwrap_or_else(|| T::one())
                / (reynolds * sqrt_f);

            let log_arg = term1 + term2;
            if log_arg <= T::zero() {
                // Fallback to previous value if log argument invalid
                break;
            }

            let inv_sqrt_f = -T::from_f64(COLEBROOK_COEFFICIENT).unwrap_or_else(|| T::one())
                * (log_arg.ln() / ln10);
            let f_next = T::one() / (inv_sqrt_f * inv_sqrt_f);

            // Check convergence
            let diff = f_next - f;
            let diff_abs = if diff >= T::zero() { diff } else { -diff };
            if diff_abs < tolerance {
                f = f_next;
                break;
            }

            f = f_next;
        }

        f
    }
}
