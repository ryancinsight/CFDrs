//! Entrance effects resistance model for developing flow regions.
//!
//! ## Theorem: Entrance Loss (Idelchik 1994 §5, Blevins 1984)
//!
//! **Problem**: A fluid flowing from a large reservoir (area A₁) into a pipe (area A₂ < A₁)
//! experiences an additional irreversible pressure loss at the inlet, beyond the
//! Hagen-Poiseuille fully-developed loss. This arises from:
//!   1. Flow contraction (vena contracta formation in the jet core).
//!   2. Turbulent mixing and expansion downstream of the vena contracta.
//!
//! **Theorem (Sudden Contraction, Idelchik 1994 §5)**:
//! The entrance loss coefficient K_entry for a sharp-edged inlet is:
//!
//! ```text
//! K_entry = 0.5 · (1 − A₂/A₁) · (1 + C/Re)   [dimensionless]
//! ```
//!
//! where:
//! - A₁ = upstream cross-sectional area [m²]
//! - A₂ = downstream cross-sectional area [m²]  (A₂ < A₁)
//! - C  = empirical low-Re correction constant ≈ 0.1 (Idelchik Table 5-1)
//! - Re = Reynolds number based on downstream hydraulic diameter
//!
//! **Physical limits**:
//! - A₂/A₁ → 0 (large reservoir → pipe): K → 0.5 × 1 = 0.5 (classical sharp-inlet Borda coefficient)
//! - A₂/A₁ = 1 (no contraction): K → 0 (no loss)
//! - Re → ∞: correction → 1 (low-Re factor vanishes)
//!
//! ## Theorem: Bernoulli Pressure Loss to Hydraulic Resistance
//!
//! The entrance pressure drop is:
//!
//! ```text
//! ΔP_entry = K_entry · (½ ρ V₂²)          [Pa]
//!          = K_entry · ρ · Q² / (2 A₂²)    [Pa]  (using V₂ = Q/A₂)
//! ```
//!
//! Converting to the lumped resistance form  ΔP = R · Q  requires knowing Q. The model
//! therefore returns a **quadratic coefficient** k (Pa·s²/m⁶) such that:
//!
//! ```text
//! ΔP = k · Q²    →    k = K_entry · ρ / (2 A₂²)
//! ```
//!
//! and the effective linearized resistance (used by the network solver) is:
//!
//! ```text
//! R_eff = k · |Q| = K_entry · ρ · |Q| / (2 A₂²)
//! ```
//!
//! ## Theorem: Smooth Contraction (Blevins 1984, ch. 6)
//!
//! For a well-rounded inlet at turbulent Reynolds numbers (Re > 10⁴):
//!
//! ```text
//! K_entry = 0.05 + 0.19 · (A₁/A₂)
//! ```
//!
//! For laminar flow through a smooth inlet, the exact Langhaar (1942) analytical
//! solution gives a total entrance coefficient (including kinetic energy change) of 1.25.
//!
//! ## Theorem: Developing Flow Length (Schlichting 1979)
//!
//! The hydrodynamic entrance length L_h required for the velocity profile to become fully developed:
//!
//! ```text
//! L_h / D_h ≈ 0.05 · Re    (laminar)
//! L_h / D_h ≈ 4.4 · Re^(1/6)  (turbulent)
//! ```
//!
//! ### References
//!
//! - Idelchik, I. E. (1994). *Handbook of Hydraulic Resistance* (3rd ed.).
//!   Mir Publishers. §5, Tables 5-1 through 5-3.
//! - Blevins, R. D. (1984). *Applied Fluid Dynamics Handbook*. Van Nostrand Reinhold.
//! - Langhaar, H. L. (1942). "Steady flow in the transition length of a straight tube."
//!   *Journal of Applied Mechanics*, 9(2), A55-A58.
//! - Schlichting, H. (1979). *Boundary Layer Theory* (7th ed.). McGraw-Hill. §9.2.

use super::traits::{FlowConditions, ResistanceModel};
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid::FluidTrait;
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;
use serde::{Deserialize, Serialize};

// Named constants for entrance effects (Idelchik 1994 §5)
const SUDDEN_CONTRACTION_CONSTANT: f64 = 0.1; // Low-Re correction coefficient (Idelchik Table 5-1)
const SMOOTH_CONTRACTION_BASE: f64 = 0.05; // Blevins (1984) smooth inlet base
const SMOOTH_CONTRACTION_SLOPE: f64 = 0.19; // Blevins (1984) smooth inlet area-ratio coefficient
const TURBULENT_TRANSITION_RE: f64 = 1e4; // Turbulent-onset Re for smooth-inlet correlation
const MAX_REYNOLDS_ENTRANCE: f64 = 1e6; // Upper Re limit for model applicability
/// Exact laminar entrance coefficient from Langhaar (1942)
const LAMINAR_SMOOTH_ENTRANCE_COEFFICIENT: f64 = 1.25;

/// Entrance effects resistance model
///
/// Models the irreversible pressure loss at the inlet of a pipe/channel
/// due to flow contraction and mixing. Supports:
/// - Sharp-edged (sudden) contraction: Idelchik (1994) §5 correlation
/// - Smooth (well-rounded) contraction: Blevins (1984) / Langhaar (1942)
/// - Interpolated intermediate geometries
///
/// # Units
///
/// The `calculate_coefficients` method returns `(R=0, k)` where `k` has units
/// of [Pa·s²/m⁶], consistent with: `ΔP = k · Q²`.
/// The `calculate_resistance` method returns the effective linearized
/// `R_eff = k · |Q|` in units of [Pa·s/m³].
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EntranceEffectsModel<T: RealField + Copy> {
    /// Upstream cross-sectional area [m²]
    pub upstream_area: T,
    /// Downstream cross-sectional area [m²]
    pub downstream_area: T,
    /// Inlet type: 0.0 = sharp sudden, 1.0 = smooth rounded
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

    /// Create a model for sudden (sharp-edged) contraction.
    ///
    /// Uses Idelchik (1994) §5: K = 0.5·(1 − A₂/A₁)·(1 + 0.1/Re)
    pub fn sudden_contraction(upstream_area: T, downstream_area: T) -> Self {
        Self::new(upstream_area, downstream_area, T::zero())
    }

    /// Create a model for smooth (well-rounded) contraction.
    ///
    /// Turbulent: K = 0.05 + 0.19·(A₁/A₂)   [Blevins 1984]
    /// Laminar:   K = 1.25                     [Langhaar 1942 exact]
    pub fn smooth_contraction(upstream_area: T, downstream_area: T) -> Self {
        Self::new(upstream_area, downstream_area, T::one())
    }
}

impl<T: RealField + Copy + FromPrimitive + num_traits::Float> ResistanceModel<T>
    for EntranceEffectsModel<T>
{
    /// Calculate effective linearised hydraulic resistance [Pa·s/m³].
    ///
    /// Returns `R_eff = k · |Q|` where `k = K_entry · ρ / (2 A₂²)`.
    /// This is only meaningful when |Q| > 0; callers should prefer
    /// `calculate_coefficients` and handle the quadratic coefficient directly.
    fn calculate_resistance<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<T> {
        let (r, k) = self.calculate_coefficients(fluid, conditions)?;
        let q_mag = if let Some(q) = conditions.flow_rate {
            if q >= T::zero() {
                q
            } else {
                -q
            }
        } else if let Some(v) = conditions.velocity {
            let v_abs = if v >= T::zero() { v } else { -v };
            v_abs * self.downstream_area
        } else {
            T::zero()
        };
        Ok(r + k * q_mag)
    }

    /// Calculate quadratic entrance loss coefficient.
    ///
    /// Returns `(R=0, k)` where:
    ///
    /// ```text
    /// k = K_entry · ρ / (2 · A₂²)    [Pa·s²/m⁶]
    /// ```
    ///
    /// and the pressure drop is `ΔP = k · Q²`.
    fn calculate_coefficients<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<(T, T)> {
        let reynolds = conditions.reynolds_number.ok_or_else(|| {
            Error::InvalidConfiguration(
                "Reynolds number required for entrance effects model".to_string(),
            )
        })?;

        let state = fluid.properties_at(conditions.temperature, conditions.pressure)?;
        let rho = state.density;

        let area_ratio = self.upstream_area / self.downstream_area;

        // K_entry [-]
        let k_entry = if self.inlet_smoothness == T::zero() {
            self.calculate_sudden_contraction_coefficient(area_ratio, reynolds)
        } else if self.inlet_smoothness == T::one() {
            self.calculate_smooth_contraction_coefficient(area_ratio, reynolds)
        } else {
            let k_sudden = self.calculate_sudden_contraction_coefficient(area_ratio, reynolds);
            let k_smooth = self.calculate_smooth_contraction_coefficient(area_ratio, reynolds);
            k_sudden + self.inlet_smoothness * (k_smooth - k_sudden)
        };

        // k = K_entry · ρ / (2 · A₂²)   [Pa·s²/m⁶]
        // Derivation: ΔP = K · ½ρV₂² = K · ρQ²/(2A₂²)
        let two = T::one() + T::one();
        let k_coeff = k_entry * rho / (two * self.downstream_area * self.downstream_area);

        Ok((T::zero(), k_coeff))
    }

    fn model_name(&self) -> &'static str {
        "Entrance Effects"
    }

    fn reynolds_range(&self) -> (T, T) {
        (
            T::zero(),
            T::from_f64(MAX_REYNOLDS_ENTRANCE)
                .expect("Mathematical constant conversion compromised"),
        )
    }

    fn validate_invariants<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<()> {
        self.validate_mach_number(fluid, conditions)?;

        // Downstream must be smaller than upstream (or equal — no-op entrance)
        if self.downstream_area > self.upstream_area {
            return Err(Error::PhysicsViolation(
                "Entrance model: downstream_area must be ≤ upstream_area".to_string(),
            ));
        }
        if self.downstream_area <= T::zero() {
            return Err(Error::PhysicsViolation(
                "Entrance model: downstream_area must be positive".to_string(),
            ));
        }
        if self.upstream_area <= T::zero() {
            return Err(Error::PhysicsViolation(
                "Entrance model: upstream_area must be positive".to_string(),
            ));
        }

        Ok(())
    }
}

impl<T: RealField + Copy + FromPrimitive + num_traits::Float> EntranceEffectsModel<T> {
    /// Compute K for sudden (sharp-edged) contraction (Idelchik 1994 §5).
    ///
    /// ```text
    /// K = 0.5 · (1 − 1/area_ratio) · (1 + C/Re)
    /// ```
    ///
    /// where `area_ratio = A₁/A₂ ≥ 1`, giving:
    /// - `(1 − 1/area_ratio)` = `(A₁ − A₂)/A₁` ∈ [0, 1)
    /// - limit area_ratio→∞: `K → 0.5·(1 + C/Re)` → **0.5** at high Re
    ///   (classical Borda-Carnot sharp-elbow coefficient)
    fn calculate_sudden_contraction_coefficient(&self, area_ratio: T, reynolds: T) -> T {
        // contraction_ratio = 1 − A₂/A₁ = 1 − 1/area_ratio  ∈ [0, 1)
        let contraction_ratio = T::one() - T::one() / area_ratio;
        // base: 0.5·(1 − A₂/A₁)
        let k_base = (T::one() / (T::one() + T::one())) * contraction_ratio;
        // Low-Re correction: multiply by (1 + C/Re)
        let c = T::from_f64(SUDDEN_CONTRACTION_CONSTANT)
            .expect("Mathematical constant conversion compromised");
        let re_correction = c / reynolds;
        k_base * (T::one() + re_correction)
    }

    /// Compute K for smooth (well-rounded) contraction.
    ///
    /// Turbulent (Re > 10⁴): `K = 0.05 + 0.19 · (A₁/A₂)`  [Blevins 1984 ch. 6]
    /// Laminar: `K = 1.25`  [Langhaar 1942 exact analytical solution]
    fn calculate_smooth_contraction_coefficient(&self, area_ratio: T, reynolds: T) -> T {
        let re_transition = T::from_f64(TURBULENT_TRANSITION_RE)
            .expect("Mathematical constant conversion compromised");
        if reynolds > re_transition {
            // Turbulent smooth inlet correlation (Blevins 1984)
            T::from_f64(SMOOTH_CONTRACTION_BASE)
                .expect("Mathematical constant conversion compromised")
                + T::from_f64(SMOOTH_CONTRACTION_SLOPE)
                    .expect("Mathematical constant conversion compromised")
                    * area_ratio
        } else {
            // Laminar: Langhaar (1942) exact total entrance coefficient = 1.25
            T::from_f64(LAMINAR_SMOOTH_ENTRANCE_COEFFICIENT)
                .expect("Mathematical constant conversion compromised")
        }
    }
}

/// Methods for combining multiple resistance models
#[derive(Debug, Clone, PartialEq)]
pub enum CombinationMethod {
    /// Add resistances in series: R_total = R₁ + R₂ + …
    Series,
    /// Combine resistances in parallel: 1/R_total = 1/R₁ + 1/R₂ + …
    Parallel,
    /// Custom weighted combination with explicit weights
    Weighted(Vec<f64>),
}

/// Durst et al. (2005) developing-flow entrance length for rectangular channels.
///
/// ## Theorem — Developing-Flow Entrance Length (Durst et al. 2005)
///
/// For laminar flow in rectangular channels, the hydrodynamic entrance
/// length (distance for the velocity profile to develop to 99% of the
/// Poiseuille parabolic profile) is:
///
/// ```text
/// L_e / D_h = [(0.619)^1.6 + (0.0567·Re)^1.6]^(1/1.6)
/// ```
///
/// This formula unifies the creeping-flow limit (L_e → 0.619·D_h as Re → 0)
/// with the high-Re asymptote (L_e → 0.0567·Re·D_h as Re → ∞).
///
/// The entrance-region resistance correction adds an additional pressure
/// drop beyond the fully-developed Hagen-Poiseuille prediction:
///
/// ```text
/// ΔP_entrance = K_entrance · (1/2)·ρ·u_mean²
/// ```
///
/// where K_entrance ≈ 2.28 + 64/(Re·D_h/L) for short channels (L < L_e).
///
/// **Reference**: Durst, F., Ray, S., Ünsal, B. & Bayoumi, O.A. (2005).
/// "The Development Lengths of Laminar Pipe and Channel Flows",
/// *J. Fluids Eng.* 127(6):1154-1160.
pub fn durst_entrance_length(re: f64, d_h: f64) -> f64 {
    let term1 = 0.619_f64.powf(1.6);
    let term2 = (0.0567 * re).powf(1.6);
    d_h * (term1 + term2).powf(1.0 / 1.6)
}

/// Entrance-region pressure drop coefficient for developing flow.
///
/// Returns K_entrance such that ΔP_entrance = K · (1/2)ρu².
/// Valid for L < L_e (entrance region not fully developed).
///
/// Uses the Langhaar (1942) + Durst (2005) composite formula.
pub fn durst_entrance_k(re: f64, l_over_dh: f64) -> f64 {
    // K ≈ 2.28 + 64/(Re·L/D_h) for laminar entrance
    // Bounded below by K_min = 0 (fully developed)
    let k = 2.28 + 64.0 / (re * l_over_dh).max(1e-30);
    k.max(0.0)
}

/// Durst entrance correction as a fraction of fully-developed resistance.
///
/// Returns a multiplier >= 1.0 for channels where L/D_h < 50 (entrance
/// effects are significant). For fully-developed channels (L/D_h >= 50),
/// returns 1.0 (no correction).
///
/// The fractional increase in pressure drop due to entrance effects is
/// approximately K/(64 * L/D_h), derived from the ratio of entrance
/// loss to Hagen-Poiseuille loss in laminar flow.
///
/// Usage: `R_corrected = R_base * durst_resistance_multiplier(Re, L/D_h)`
///
/// ## Reference
///
/// Durst, F., Ray, S., Unsal, B. & Bayoumi, O.A. (2005).
/// "The Development Lengths of Laminar Pipe and Channel Flows",
/// *J. Fluids Eng.* 127(6):1154-1160.
#[inline]
pub fn durst_resistance_multiplier(re: f64, l_over_dh: f64) -> f64 {
    if l_over_dh >= 50.0 {
        return 1.0;
    }
    let k = durst_entrance_k(re, l_over_dh);
    1.0 + k / (64.0 * l_over_dh).max(1.0)
}

// ─── Tests ─────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use cfd_core::error::Result;

    fn water() -> impl FluidTrait<f64> {
        cfd_core::physics::fluid::database::water_20c::<f64>().unwrap()
    }

    // ─── Dimensional correctness ─────────────────────────────────────────

    /// Verify the quadratic coefficient has correct sign and units.
    ///
    /// For water (ρ ≈ 998.2 kg/m³), A₂ = 1 mm², A₁ = 4 mm², Re = 1000:
    ///   contraction_ratio = 1 − 0.25 = 0.75
    ///   K_base = 0.5 × 0.75 = 0.375
    ///   k_entry = 0.375 × (1 + 0.1/1000) ≈ 0.3754
    ///   k_coeff = 0.3754 × 998.2 / (2 × (1e-6)²) ≈ 1.873 × 10¹⁴ Pa·s²/m⁶
    #[test]
    fn test_sudden_contraction_coefficient_dimensional() -> Result<()> {
        let model = EntranceEffectsModel::sudden_contraction(4e-6_f64, 1e-6_f64);
        let fluid = water();
        let mut conditions = FlowConditions::new(0.1_f64);
        conditions.reynolds_number = Some(1000.0);

        let (r, k) = model.calculate_coefficients(&fluid, &conditions)?;
        assert_eq!(
            r, 0.0,
            "Linear coefficient should be zero for entrance effects"
        );
        assert!(k > 0.0, "Quadratic coefficient must be positive, got {}", k);
        assert!(k.is_finite(), "Quadratic coefficient must be finite");

        // Verify magnitude: k ≈ 0.375 × 998.2 / (2 × 1e-12) ≈ 1.87e14
        let k_entry = 0.375_f64 * (1.0 + 0.1 / 1000.0);
        let expected = k_entry * 998.2 / (2.0 * 1e-6_f64 * 1e-6_f64);
        assert_relative_eq!(k, expected, max_relative = 0.01);
        Ok(())
    }

    /// Verify K → 0.5 for large reservoir → pipe (area_ratio → ∞).
    ///
    /// Theorem (Idelchik §5): at high Re, K_entry → 0.5 (Borda-Carnot sharp inlet)
    #[test]
    fn test_sudden_contraction_borda_carnot_limit() -> Result<()> {
        // Very large upstream area (≈ reservoir)
        let model = EntranceEffectsModel::sudden_contraction(1.0_f64, 1e-6_f64);
        let fluid = water();
        let mut conditions = FlowConditions::new(0.1_f64);
        conditions.reynolds_number = Some(1e6); // High Re → correction ≈ 1

        let (_, k) = model.calculate_coefficients(&fluid, &conditions)?;
        // K_entry ≈ 0.5 at high Re and large area_ratio
        // k = 0.5 × ρ / (2 A₂²)
        let rho = 998.2_f64;
        let expected_k = 0.5 * rho / (2.0 * 1e-6 * 1e-6);
        assert_relative_eq!(k, expected_k, max_relative = 0.01);
        Ok(())
    }

    /// Verify K → 0 when A₂ = A₁ (no contraction).
    #[test]
    fn test_sudden_contraction_no_contraction_limit() -> Result<()> {
        // Equal areas: no contraction loss
        let a = 1e-6_f64;
        let model = EntranceEffectsModel::sudden_contraction(a, a);
        let fluid = water();
        let mut conditions = FlowConditions::new(0.1_f64);
        conditions.reynolds_number = Some(1000.0);

        let (_, k) = model.calculate_coefficients(&fluid, &conditions)?;
        // contraction_ratio = 1 − 1/1 = 0 → K = 0 → k = 0
        assert_relative_eq!(k, 0.0, epsilon = 1e-20);
        Ok(())
    }

    /// Entrance resistance decreases with increasing Re (1/Re correction).
    #[test]
    fn test_sudden_contraction_reynolds_dependence() -> Result<()> {
        let model = EntranceEffectsModel::sudden_contraction(4e-6_f64, 1e-6_f64);
        let fluid = water();

        let make_cond = |re: f64| {
            let mut c = FlowConditions::new(0.1_f64);
            c.reynolds_number = Some(re);
            c
        };

        let (_, k_low) = model.calculate_coefficients(&fluid, &make_cond(100.0))?;
        let (_, k_high) = model.calculate_coefficients(&fluid, &make_cond(1000.0))?;

        // Higher Re → smaller 1/Re correction → smaller K → smaller k
        assert!(
            k_high < k_low,
            "k at Re=1000 ({}) should be < k at Re=100 ({})",
            k_high,
            k_low
        );
        Ok(())
    }

    /// Validate (A₁/A₂)² area scaling of effective resistance.
    ///
    /// With the same K_entry, scaling A₂ by factor s scales k by 1/s⁴.
    #[test]
    fn test_entrance_effects_area_scaling() -> Result<()> {
        let base_down = 1e-6_f64;
        let s = 2.0_f64;
        let model_base = EntranceEffectsModel::sudden_contraction(4.0 * base_down, base_down);
        // Scale both areas by s², keeping area ratio constant
        let model_scaled =
            EntranceEffectsModel::sudden_contraction(4.0 * base_down * s * s, base_down * s * s);
        let fluid = water();
        let mut conditions = FlowConditions::new(0.1_f64);
        conditions.reynolds_number = Some(1000.0);

        let (_, k_base) = model_base.calculate_coefficients(&fluid, &conditions)?;
        let (_, k_scaled) = model_scaled.calculate_coefficients(&fluid, &conditions)?;

        // k_scaled / k_base = (A₂_base / A₂_scaled)² = 1/s⁴
        let expected_ratio = 1.0 / (s * s * s * s);
        assert_relative_eq!(k_scaled / k_base, expected_ratio, epsilon = 1e-6);
        Ok(())
    }

    /// Smooth contraction turbulent: K = 0.05 + 0.19·(A₁/A₂) [Blevins 1984]
    #[test]
    fn test_smooth_contraction_turbulent_blevins() -> Result<()> {
        let upstream = 2e-6_f64;
        let downstream = 1e-6_f64;
        let model = EntranceEffectsModel::smooth_contraction(upstream, downstream);
        let fluid = water();
        let mut conditions = FlowConditions::new(0.1_f64);
        conditions.reynolds_number = Some(50_000.0); // Turbulent

        let (_, k) = model.calculate_coefficients(&fluid, &conditions)?;
        let area_ratio = 2.0_f64;
        let k_expected = (0.05 + 0.19 * area_ratio) * 998.2 / (2.0 * downstream * downstream);
        assert_relative_eq!(k, k_expected, epsilon = 1e-6);
        Ok(())
    }

    /// Validate invariant: downstream > upstream returns error.
    #[test]
    fn test_invariant_downstream_larger_than_upstream() {
        let model = EntranceEffectsModel::<f64>::sudden_contraction(1e-6, 4e-6); // reversed
        let fluid = water();
        let mut conditions = FlowConditions::new(0.1_f64);
        conditions.reynolds_number = Some(1000.0);
        assert!(model.validate_invariants(&fluid, &conditions).is_err());
    }

    /// Validate invariant: zero downstream area returns error.
    #[test]
    fn test_invariant_zero_downstream_area() {
        let model = EntranceEffectsModel::<f64>::sudden_contraction(1e-6, 0.0);
        let fluid = water();
        let mut conditions = FlowConditions::new(0.1_f64);
        conditions.reynolds_number = Some(1000.0);
        assert!(model.validate_invariants(&fluid, &conditions).is_err());
    }

    // ─── Durst et al. (2005) entrance length tests ──────────────────────

    /// Re → 0 gives L_e ≈ 0.619·D_h (creeping-flow limit).
    ///
    /// At Re = 0 the high-Re term vanishes and L_e/D_h → 0.619.
    #[test]
    fn test_durst_entrance_length_creeping_flow() {
        let d_h = 1.0e-3; // 1 mm
        let l_e = durst_entrance_length(0.0, d_h);
        let expected = 0.619 * d_h;
        assert_relative_eq!(l_e, expected, max_relative = 1e-6);
    }

    /// Re = 1000: L_e ≈ 0.0567·Re·D_h = 56.7·D_h (high-Re asymptote dominates).
    #[test]
    fn test_durst_entrance_length_high_re() {
        let d_h = 1.0e-3;
        let re = 1000.0;
        let l_e = durst_entrance_length(re, d_h);
        // At Re=1000 the 0.0567·Re term dominates: 56.7 >> 0.619
        let high_re_approx = 0.0567 * re * d_h;
        // Should be very close to the high-Re asymptote
        assert_relative_eq!(l_e, high_re_approx, max_relative = 0.01);
    }

    /// Short channel (L/D_h = 1): large K (entrance effects dominate).
    #[test]
    fn test_durst_entrance_k_short_channel() {
        let re = 100.0;
        let l_over_dh = 1.0;
        let k = durst_entrance_k(re, l_over_dh);
        // K = 2.28 + 64/(100·1) = 2.28 + 0.64 = 2.92
        assert_relative_eq!(k, 2.92, max_relative = 1e-6);
        assert!(k > 2.0, "Short channel should have large entrance K");
    }

    /// Long channel (L/D_h = 100): small K (fully developed flow).
    #[test]
    fn test_durst_entrance_k_long_channel() {
        let re = 100.0;
        let l_over_dh = 100.0;
        let k = durst_entrance_k(re, l_over_dh);
        // K = 2.28 + 64/(100·100) = 2.28 + 0.0064 ≈ 2.2864
        assert_relative_eq!(k, 2.2864, max_relative = 1e-3);
        // Should be much smaller than the short-channel case
        let k_short = durst_entrance_k(re, 1.0);
        assert!(
            k < k_short,
            "Long channel K ({k}) should be < short channel K ({k_short})"
        );
    }

    // ─── Durst resistance multiplier tests ──────────────────────────────

    /// Long channel (L/D_h = 100 >= 50): multiplier should be exactly 1.0.
    #[test]
    fn test_durst_multiplier_long_channel_unity() {
        let multiplier = durst_resistance_multiplier(100.0, 100.0);
        assert_relative_eq!(multiplier, 1.0, epsilon = 1e-15);
    }

    /// Short channel (L/D_h = 5 < 50): multiplier should be > 1.0.
    #[test]
    fn test_durst_multiplier_short_channel_above_unity() {
        let multiplier = durst_resistance_multiplier(100.0, 5.0);
        assert!(
            multiplier > 1.0,
            "Multiplier for short channel (L/D_h=5) should be > 1.0, got {multiplier}"
        );
        // K at Re=100, L/D_h=5: 2.28 + 64/(100*5) = 2.28 + 0.128 = 2.408
        // multiplier = 1 + 2.408 / (64*5) = 1 + 2.408/320 = 1.007525
        let expected = 1.0 + (2.28 + 64.0 / (100.0 * 5.0)) / (64.0 * 5.0);
        assert_relative_eq!(multiplier, expected, max_relative = 1e-10);
    }
}
