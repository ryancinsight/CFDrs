//! Entrance effects resistance model for developing flow regions.
//!
//! ## Entrance Effects Theorem
//!
//! **Theorem**: For flow entering a channel or pipe, additional pressure losses occur
//! in the entrance region where the velocity profile develops from uniform (at inlet)
//! to the fully developed parabolic profile.
//!
//! The entrance loss coefficient K_entry depends on the channel geometry and Reynolds number:
//!
//! ### Sudden Contraction (Sharp-Edged Inlet)
//!
//! K_entry = 0.5·(1 − A₂/A₁) · (1 + C/Re_D₂)
//!
//! where:
//! - A₁ is the upstream cross-sectional area (larger)
//! - A₂ is the downstream cross-sectional area (smaller)
//! - C is a constant (typically 0.1)
//! - Re_D₂ is Reynolds number based on downstream hydraulic diameter
//!
//! Physical limits:
//! - A₂/A₁ → 0 (large reservoir → pipe): K → 0.5  (Idelchik, §5)
//! - A₂/A₁ = 1 (no contraction):          K → 0
//!
//! ### Smooth Contraction (Well-Rounded Inlet)
//!
//! K_entry = 0.05 + 0.19 (A₁/A₂) for Re > 10⁴
//!
//! ### Developing Flow Length
//!
//! The hydrodynamic entrance length L_h is given by:
//!
//! L_h/D_h ≈ 0.05 Re_D_h for laminar flow
//! L_h/D_h ≈ 4.4 (Re_D_h)^(1/6) for turbulent flow
//!
//! ### Hydraulic Resistance
//!
//! The entrance resistance is added to the fully developed flow resistance:
//!
//! R_total = R_fully_developed + R_entrance
//!
//! where R_entrance = K_entry * ρ V² / (2 ΔP)
//!
//! ### Validity Conditions
//!
//! 1. **Sharp Inlet**: K_entry ≈ 0.5 for sudden contraction
//! 2. **Reynolds Number**: Effects decrease with increasing Re
//! 3. **Geometry**: Strong dependence on contraction ratio A₁/A₂
//! 4. **Flow Type**: Different correlations for laminar vs turbulent
//!
//! ### References
//!
//! - White, F. M. (2006). *Viscous Fluid Flow* (3rd ed.). McGraw-Hill. Section 4.11.
//! - Idelchik, I. E. (1994). *Handbook of Hydraulic Resistance* (3rd ed.).
//! - Blevins, R. D. (1984). *Applied Fluid Dynamics Handbook*.

use super::traits::{FlowConditions, ResistanceModel};
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid::FluidTrait;
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;
use serde::{Deserialize, Serialize};

// Named constants for entrance effects
const SUDDEN_CONTRACTION_CONSTANT: f64 = 0.1;
const SMOOTH_CONTRACTION_BASE: f64 = 0.05;
const SMOOTH_CONTRACTION_SLOPE: f64 = 0.19;
const TURBULENT_TRANSITION_RE: f64 = 1e4;
const MAX_REYNOLDS_ENTRANCE: f64 = 1e6;

/// Entrance effects resistance model
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EntranceEffectsModel<T: RealField + Copy> {
    /// Upstream cross-sectional area [m²]
    pub upstream_area: T,
    /// Downstream cross-sectional area [m²]
    pub downstream_area: T,
    /// Inlet type coefficient (0.0 = sharp, 1.0 = smooth)
    pub inlet_smoothness: T,
}

impl<T: RealField + Copy> EntranceEffectsModel<T> {
    /// Create a new entrance effects model
    pub fn new(upstream_area: T, downstream_area: T, inlet_smoothness: T) -> Self {
        Self {
            upstream_area,
            downstream_area,
            inlet_smoothness,
        }
    }

    /// Create a model for sudden contraction (sharp-edged inlet)
    pub fn sudden_contraction(upstream_area: T, downstream_area: T) -> Self {
        Self::new(upstream_area, downstream_area, T::zero())
    }

    /// Create a model for smooth contraction (well-rounded inlet)
    pub fn smooth_contraction(upstream_area: T, downstream_area: T) -> Self {
        Self::new(upstream_area, downstream_area, T::one())
    }
}

impl<T: RealField + Copy + FromPrimitive + num_traits::Float> ResistanceModel<T>
    for EntranceEffectsModel<T>
{
    fn calculate_resistance<F: FluidTrait<T>>(
        &self,
        _fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<T> {
        let reynolds = conditions.reynolds_number.ok_or_else(|| {
            Error::InvalidConfiguration(
                "Reynolds number required for entrance effects model".to_string(),
            )
        })?;

        let area_ratio = self.upstream_area / self.downstream_area;

        // Calculate entrance loss coefficient
        let k_entry = if self.inlet_smoothness == T::zero() {
            // Sudden contraction (sharp inlet)
            self.calculate_sudden_contraction_coefficient(area_ratio, reynolds)
        } else if self.inlet_smoothness == T::one() {
            // Smooth contraction (rounded inlet)
            self.calculate_smooth_contraction_coefficient(area_ratio, reynolds)
        } else {
            // Interpolated between sudden and smooth
            let k_sudden = self.calculate_sudden_contraction_coefficient(area_ratio, reynolds);
            let k_smooth = self.calculate_smooth_contraction_coefficient(area_ratio, reynolds);
            k_sudden + self.inlet_smoothness * (k_smooth - k_sudden)
        };

        // Convert loss coefficient to hydraulic resistance
        // R = K_entry / (2 A_downstream²)
        // This gives R such that ΔP = R * Q²
        let resistance = k_entry
            / (T::from_f64(2.0).unwrap_or_else(|| T::zero())
                * self.downstream_area
                * self.downstream_area);

        Ok(resistance)
    }

    fn model_name(&self) -> &'static str {
        "Entrance Effects"
    }

    fn reynolds_range(&self) -> (T, T) {
        (
            T::zero(),
            T::from_f64(MAX_REYNOLDS_ENTRANCE).unwrap_or_else(|| T::zero()),
        )
    }
}

impl<T: RealField + Copy + FromPrimitive + num_traits::Float> EntranceEffectsModel<T> {
    /// Calculate loss coefficient for sudden (sharp-edged) contraction.
    ///
    /// Uses the Idelchik (1994) linear approximation:
    ///
    /// ```text
    /// K = 0.5 · (1 − A₂/A₁) · (1 + C/Re)
    ///   = 0.5 · (1 − 1/area_ratio) · (1 + C/Re)
    /// ```
    ///
    /// where `area_ratio = A₁/A₂ ≥ 1`.  As `area_ratio → ∞` (large reservoir
    /// entering a pipe), `K → 0.5`, matching the well-known sharp-edged inlet
    /// coefficient (White 2011, §6.7; Idelchik 1994, §5).
    fn calculate_sudden_contraction_coefficient(&self, area_ratio: T, reynolds: T) -> T {
        // 1 − A₂/A₁ = 1 − 1/area_ratio  ∈ [0, 1)
        let contraction_ratio = T::one() - T::one() / area_ratio;

        // Base loss coefficient: 0.5·(1 − A₂/A₁)
        let k_base = T::from_f64(0.5).unwrap_or_else(T::one) * contraction_ratio;

        // Low-Reynolds correction (Idelchik): multiply by (1 + C/Re)
        let re_correction =
            T::from_f64(SUDDEN_CONTRACTION_CONSTANT).unwrap_or_else(|| T::zero()) / reynolds;

        k_base * (T::one() + re_correction)
    }

    /// Calculate loss coefficient for smooth contraction
    fn calculate_smooth_contraction_coefficient(&self, area_ratio: T, reynolds: T) -> T {
        let re_transition = T::from_f64(TURBULENT_TRANSITION_RE).unwrap_or_else(|| T::zero());

        if reynolds > re_transition {
            // Turbulent flow correlation (Blevins)
            T::from_f64(SMOOTH_CONTRACTION_BASE).unwrap_or_else(|| T::zero())
                + T::from_f64(SMOOTH_CONTRACTION_SLOPE).unwrap_or_else(|| T::zero()) * area_ratio
        } else {
            // Laminar flow: derive from developing flow theory
            // For smooth contractions in laminar flow, the loss coefficient
            // is derived from the Hagen-Poiseuille developing flow analysis.
            // The entrance loss arises from the momentum deficit in the
            // developing boundary layer region.
            //
            // Theoretical estimate: K ≈ 1.0 for laminar sudden contraction
            // based on the velocity profile development from uniform to parabolic.
            // This is significantly higher than the previous hardcoded 0.01.
            T::from_f64(1.0).unwrap_or_else(|| T::one())
        }
    }
}

/// Methods for combining multiple resistance models
#[derive(Debug, Clone, PartialEq)]
pub enum CombinationMethod {
    /// Add resistances in series
    Series,
    /// Add resistances in parallel
    Parallel,
    /// Custom weighted combination
    Weighted(Vec<f64>),
}
