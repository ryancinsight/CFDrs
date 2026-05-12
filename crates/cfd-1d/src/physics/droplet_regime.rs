//! Capillary and Weber number regime classification for two-phase flows.
//!
//! # Theorem — Dimensionless Regime Map (De Menech et al. 2008)
//!
//! Two-phase microfluidic flows (immiscible droplets in a carrier fluid) are
//! governed by the competition between viscous, inertial, and surface tension
//! forces. The flow regime is determined by two dimensionless numbers:
//!
//! ## Capillary number
//!
//! ```text
//! Ca = µ · U / σ
//! ```
//!
//! Ratio of viscous to surface tension forces. Low Ca (< 0.01) → surface
//! tension dominates; high Ca (> 0.1) → viscous forces dominate.
//!
//! ## Weber number
//!
//! ```text
//! We = ρ · U² · L / σ
//! ```
//!
//! Ratio of inertial to surface tension forces. Relevant for high flow rates
//! where inertia cannot be neglected.
//!
//! ## Regime boundaries
//!
//! | Regime     | Ca range         | Physics                                  |
//! |------------|------------------|------------------------------------------|
//! | Squeezing  | Ca < 0.01        | Plug fills channel; backpressure pinch-off |
//! | Dripping   | 0.01 ≤ Ca < 0.1  | Shear-driven neck thinning; regular drops  |
//! | Jetting    | Ca ≥ 0.1         | Continuous thread → Rayleigh-Plateau breakup |
//! | Parallel   | Ca > 1.0         | Co-flow; no droplet formation              |
//!
//! We > 1 indicates inertia-dominated breakup (tip-streaming at high We).
//!
//! # References
//! - De Menech, M., Garstecki, P., Jousse, F. & Stone, H.A. (2008). Transition
//!   from squeezing to dripping in a microfluidic T-junction. *J. Fluid Mech.*
//!   595:141-161.
//! - Christopher, G.F. & Anna, S.L. (2007). Microfluidic methods for generating
//!   continuous droplet streams. *J. Phys. D: Appl. Phys.* 40:R319-R336.
//! - Garstecki, P. et al. (2006). Formation of droplets and bubbles in a
//!   microfluidic T-junction. *Lab Chip* 6:437-446.

use std::fmt;

use cfd_core::error::{Error, Result};

/// Two-phase flow regime classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FlowRegime {
    /// Ca < 0.01: Surface tension dominates. Droplets fill the channel
    /// cross-section; pinch-off is driven by upstream backpressure.
    /// Produces large, monodisperse plugs.
    Squeezing,

    /// 0.01 ≤ Ca < 0.1: Viscous shear thins the neck of the dispersed
    /// phase at the junction. Produces regular, smaller droplets.
    Dripping,

    /// 0.1 ≤ Ca < 1.0: The dispersed phase forms a continuous thread
    /// that breaks up via Rayleigh-Plateau instability downstream of
    /// the junction. Droplet size distribution is broader.
    Jetting,

    /// Ca ≥ 1.0: Viscous forces completely suppress surface-tension-driven
    /// breakup. Both phases flow in parallel (stratified flow).
    Parallel,
}

impl fmt::Display for FlowRegime {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Squeezing => write!(f, "Squeezing"),
            Self::Dripping => write!(f, "Dripping"),
            Self::Jetting => write!(f, "Jetting"),
            Self::Parallel => write!(f, "Parallel"),
        }
    }
}

/// Compute the capillary number.
///
/// ```text
/// Ca = µ · U / σ
/// ```
///
/// # Arguments
/// * `viscosity_pa_s` — Dynamic viscosity of the continuous phase [Pa·s]
/// * `velocity_m_s` — Mean flow velocity [m/s]
/// * `surface_tension_n_m` — Interfacial tension [N/m]
///
/// # Returns
/// Capillary number (dimensionless).
///
/// # Errors
/// Returns [`Error::PhysicsViolation`] when viscosity is negative, surface
/// tension is nonpositive, or any input is nonfinite.
#[inline]
pub fn capillary_number(
    viscosity_pa_s: f64,
    velocity_m_s: f64,
    surface_tension_n_m: f64,
) -> Result<f64> {
    validate_nonnegative_finite("viscosity_pa_s", viscosity_pa_s)?;
    validate_finite("velocity_m_s", velocity_m_s)?;
    validate_positive_finite("surface_tension_n_m", surface_tension_n_m)?;
    Ok(viscosity_pa_s * velocity_m_s.abs() / surface_tension_n_m)
}

/// Compute the Weber number.
///
/// ```text
/// We = ρ · U² · L / σ
/// ```
///
/// # Arguments
/// * `density_kg_m3` — Fluid density [kg/m³]
/// * `velocity_m_s` — Mean flow velocity [m/s]
/// * `length_m` — Characteristic length (channel hydraulic diameter) [m]
/// * `surface_tension_n_m` — Interfacial tension [N/m]
///
/// # Returns
/// Weber number (dimensionless).
///
/// # Errors
/// Returns [`Error::PhysicsViolation`] when density, length, or surface tension
/// is outside its physical domain, or when any input is nonfinite.
#[inline]
pub fn weber_number(
    density_kg_m3: f64,
    velocity_m_s: f64,
    length_m: f64,
    surface_tension_n_m: f64,
) -> Result<f64> {
    validate_nonnegative_finite("density_kg_m3", density_kg_m3)?;
    validate_finite("velocity_m_s", velocity_m_s)?;
    validate_positive_finite("length_m", length_m)?;
    validate_positive_finite("surface_tension_n_m", surface_tension_n_m)?;
    Ok(density_kg_m3 * velocity_m_s * velocity_m_s * length_m / surface_tension_n_m)
}

/// Ohnesorge number: ratio of viscous to inertia-surface-tension forces.
///
/// ```text
/// Oh = µ / √(ρ · σ · L) = √(We) / Re
/// ```
///
/// Oh > 1: viscosity dominates over inertia (small channels, high viscosity).
/// Oh < 1: inertia dominates (large channels, low viscosity).
#[inline]
pub fn ohnesorge_number(
    viscosity_pa_s: f64,
    density_kg_m3: f64,
    surface_tension_n_m: f64,
    length_m: f64,
) -> Result<f64> {
    validate_nonnegative_finite("viscosity_pa_s", viscosity_pa_s)?;
    validate_positive_finite("density_kg_m3", density_kg_m3)?;
    validate_positive_finite("surface_tension_n_m", surface_tension_n_m)?;
    validate_positive_finite("length_m", length_m)?;
    Ok(viscosity_pa_s / (density_kg_m3 * surface_tension_n_m * length_m).sqrt())
}

/// Classify the two-phase flow regime from the capillary number.
///
/// # Theorem — Regime Classification (De Menech 2008)
///
/// The regime boundaries are at $Ca = 0.01$, $Ca = 0.1$, and $Ca = 1.0$:
///
/// ```text
/// Ca < 0.01  → Squeezing
/// Ca ∈ [0.01, 0.1) → Dripping
/// Ca ∈ [0.1, 1.0)  → Jetting
/// Ca ≥ 1.0         → Parallel
/// ```
///
/// **Proof of completeness**: The four regimes partition $[0, \infty)$ via
/// three threshold values, covering all non-negative Ca. ∎
#[inline]
pub fn classify_regime(ca: f64) -> Result<FlowRegime> {
    validate_nonnegative_finite("capillary_number", ca)?;
    if ca < 0.01 {
        Ok(FlowRegime::Squeezing)
    } else if ca < 0.1 {
        Ok(FlowRegime::Dripping)
    } else if ca < 1.0 {
        Ok(FlowRegime::Jetting)
    } else {
        Ok(FlowRegime::Parallel)
    }
}

/// Full regime analysis: compute Ca, We, Oh, and classify.
///
/// # Arguments
/// * `viscosity_pa_s` — Dynamic viscosity [Pa·s]
/// * `density_kg_m3` — Fluid density [kg/m³]
/// * `velocity_m_s` — Mean velocity [m/s]
/// * `length_m` — Characteristic length [m]
/// * `surface_tension_n_m` — Interfacial tension [N/m]
///
/// # Returns
/// `(Ca, We, Oh, FlowRegime)`.
///
/// # Errors
/// Returns [`Error::PhysicsViolation`] when any dimensional input violates the
/// physical domain of the dimensionless group being evaluated.
pub fn regime_analysis(
    viscosity_pa_s: f64,
    density_kg_m3: f64,
    velocity_m_s: f64,
    length_m: f64,
    surface_tension_n_m: f64,
) -> Result<(f64, f64, f64, FlowRegime)> {
    let ca = capillary_number(viscosity_pa_s, velocity_m_s, surface_tension_n_m)?;
    let we = weber_number(density_kg_m3, velocity_m_s, length_m, surface_tension_n_m)?;
    let oh = ohnesorge_number(viscosity_pa_s, density_kg_m3, surface_tension_n_m, length_m)?;
    let regime = classify_regime(ca)?;
    Ok((ca, we, oh, regime))
}

fn validate_finite(name: &str, value: f64) -> Result<()> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(Error::PhysicsViolation(format!("{name} must be finite")))
    }
}

fn validate_nonnegative_finite(name: &str, value: f64) -> Result<()> {
    validate_finite(name, value)?;
    if value >= 0.0 {
        Ok(())
    } else {
        Err(Error::PhysicsViolation(format!(
            "{name} must be nonnegative"
        )))
    }
}

fn validate_positive_finite(name: &str, value: f64) -> Result<()> {
    validate_finite(name, value)?;
    if value > 0.0 {
        Ok(())
    } else {
        Err(Error::PhysicsViolation(format!("{name} must be positive")))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const MU_WATER: f64 = 0.001; // Pa·s
    const RHO_WATER: f64 = 1000.0; // kg/m³
    const SIGMA_WATER_OIL: f64 = 0.035; // N/m (typical water-oil)

    /// Very slow flow → squeezing regime.
    #[test]
    fn slow_flow_squeezing() -> Result<()> {
        let ca = capillary_number(MU_WATER, 1e-4, SIGMA_WATER_OIL)?;
        assert!(ca < 0.01, "Ca={ca:.6} should be < 0.01 for squeezing");
        assert_eq!(classify_regime(ca)?, FlowRegime::Squeezing);
        Ok(())
    }

    /// Moderate flow → dripping regime.
    #[test]
    fn moderate_flow_dripping() -> Result<()> {
        // Ca ≈ 0.001 * 1.0 / 0.035 ≈ 0.029
        let ca = capillary_number(MU_WATER, 1.0, SIGMA_WATER_OIL)?;
        assert!(
            ca >= 0.01 && ca < 0.1,
            "Ca={ca:.6} should be in dripping range [0.01, 0.1)"
        );
        assert_eq!(classify_regime(ca)?, FlowRegime::Dripping);
        Ok(())
    }

    /// Fast flow → jetting regime.
    #[test]
    fn fast_flow_jetting() -> Result<()> {
        let ca = capillary_number(MU_WATER, 5.0, SIGMA_WATER_OIL)?;
        assert!(
            ca >= 0.1 && ca < 1.0,
            "Ca={ca:.6} should be in jetting range [0.1, 1.0)"
        );
        assert_eq!(classify_regime(ca)?, FlowRegime::Jetting);
        Ok(())
    }

    /// Very fast flow → parallel regime.
    #[test]
    fn very_fast_flow_parallel() -> Result<()> {
        let ca = capillary_number(MU_WATER, 50.0, SIGMA_WATER_OIL)?;
        assert!(ca >= 1.0, "Ca={ca:.6} should be >= 1.0 for parallel");
        assert_eq!(classify_regime(ca)?, FlowRegime::Parallel);
        Ok(())
    }

    /// Ca scales linearly with velocity.
    #[test]
    fn ca_scales_with_velocity() -> Result<()> {
        let ca1 = capillary_number(MU_WATER, 1.0, SIGMA_WATER_OIL)?;
        let ca2 = capillary_number(MU_WATER, 2.0, SIGMA_WATER_OIL)?;
        let ratio = ca2 / ca1;
        assert!(
            (ratio - 2.0).abs() < 1e-10,
            "Ca ratio {ratio:.6} should be 2.0"
        );
        Ok(())
    }

    /// We scales with velocity squared.
    #[test]
    fn we_scales_with_velocity_squared() -> Result<()> {
        let l = 100e-6;
        let we1 = weber_number(RHO_WATER, 1.0, l, SIGMA_WATER_OIL)?;
        let we2 = weber_number(RHO_WATER, 2.0, l, SIGMA_WATER_OIL)?;
        let ratio = we2 / we1;
        assert!(
            (ratio - 4.0).abs() < 1e-10,
            "We ratio {ratio:.6} should be 4.0"
        );
        Ok(())
    }

    /// Ohnesorge is always positive.
    #[test]
    fn ohnesorge_positive() -> Result<()> {
        let oh = ohnesorge_number(MU_WATER, RHO_WATER, SIGMA_WATER_OIL, 100e-6)?;
        assert!(oh > 0.0, "Oh should be positive, got {oh:.6}");
        Ok(())
    }

    /// Regime analysis returns consistent results.
    #[test]
    fn regime_analysis_consistent() -> Result<()> {
        let (ca, we, oh, regime) =
            regime_analysis(MU_WATER, RHO_WATER, 1.0, 100e-6, SIGMA_WATER_OIL)?;
        assert!(ca > 0.0 && we > 0.0 && oh > 0.0);
        assert_eq!(classify_regime(ca)?, regime);
        Ok(())
    }

    #[test]
    fn nonpositive_surface_tension_rejected() {
        let err = capillary_number(MU_WATER, 1.0, 0.0)
            .expect_err("zero surface tension makes capillary number undefined");
        assert!(
            err.to_string().contains("surface_tension_n_m"),
            "error must identify invalid surface tension: {err}"
        );

        let err = weber_number(RHO_WATER, 1.0, 100e-6, -SIGMA_WATER_OIL)
            .expect_err("negative surface tension makes Weber number undefined");
        assert!(
            err.to_string().contains("surface_tension_n_m"),
            "error must identify invalid surface tension: {err}"
        );
    }

    #[test]
    fn nonpositive_length_rejected_for_length_dependent_groups() {
        let err = ohnesorge_number(MU_WATER, RHO_WATER, SIGMA_WATER_OIL, 0.0)
            .expect_err("zero characteristic length makes Ohnesorge number undefined");
        assert!(
            err.to_string().contains("length_m"),
            "error must identify invalid characteristic length: {err}"
        );

        let err = weber_number(RHO_WATER, 1.0, -100e-6, SIGMA_WATER_OIL)
            .expect_err("negative characteristic length makes Weber number undefined");
        assert!(
            err.to_string().contains("length_m"),
            "error must identify invalid characteristic length: {err}"
        );
    }

    #[test]
    fn invalid_capillary_number_rejected_for_classification() {
        let err =
            classify_regime(f64::NAN).expect_err("NaN capillary number cannot define a regime");
        assert!(
            err.to_string().contains("capillary_number"),
            "error must identify invalid capillary number: {err}"
        );
    }
}
