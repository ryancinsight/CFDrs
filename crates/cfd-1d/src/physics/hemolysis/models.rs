use aequitas::systems::si::quantities::{Pressure, Time};

/// Giersiepen (1990) haemolysis index: `HI = C · t^α · τ^β`.
///
/// Thin wrapper delegating to [`cfd_core::physics::hemolysis::HemolysisModel::giersiepen_millifluidic`]
/// so that the single source of truth for all model constants lives in cfd-core.
/// Returns `0.0` for negative or zero inputs, including negative infinity.
/// NaN inputs propagate as NaN; positive infinity follows the IEEE-754
/// arithmetic of the provider model.
///
/// # Arguments
///
/// * `shear`    — wall shear stress as an Aequitas `Pressure`
/// * `duration` — exposure duration as an Aequitas `Time`
///
/// # Example
///
/// ```
/// use aequitas::systems::si::quantities::{Pressure, Time};
/// use cfd_1d::giersiepen_hi;
/// let hi = giersiepen_hi(Pressure::from_base(50.0), Time::from_base(0.01));
/// assert!((hi - 2.578_444_864_181_722_3e-3).abs() < 1e-15);
/// ```
#[inline]
#[must_use]
pub fn giersiepen_hi(shear: Pressure, duration: Time) -> f64 {
    let shear_pa = shear.into_base();
    let duration_s = duration.into_base();
    if shear_pa.is_nan() || duration_s.is_nan() {
        return f64::NAN;
    }
    if shear_pa == 0.0 || duration_s == 0.0 {
        return 0.0;
    }
    if shear_pa < 0.0 || duration_s < 0.0 {
        return 0.0;
    }

    cfd_core::physics::hemolysis::HemolysisModel::giersiepen_millifluidic()
        .damage_index(shear_pa, duration_s)
        .expect("invariant: non-negative shear and duration satisfy the hemolysis model")
}

/// Amplify a baseline Giersiepen HI by the local cavitation potential.
///
/// Thin wrapper delegating to [`cfd_core::physics::hemolysis::HemolysisModel::cavitation_amplified`]
/// so that the single source of truth for the amplification slope lives in cfd-core.
///
/// Bubble collapse generates micro-jets (localised shear amplification) and
/// pressure shockwaves that damage RBC membranes independently of macroscopic
/// steady shear.
///
/// Conservative model: 3× amplification at full cavitation potential
/// (`cav_potential = 1.0`).
///
/// ```text
///   HI_amplified = base_hi × (1 + 3 × cav_potential.clamp(0, 1))
/// ```
///
/// # Arguments
///
/// * `base_hi`          — baseline Giersiepen HI (from [`giersiepen_hi`])
/// * `cav_potential`    — cavitation potential ∈ [0, 1]; 0 = no cavitation
///
/// # Example
///
/// ```
/// use aequitas::systems::si::quantities::{Pressure, Time};
/// use cfd_1d::{cavitation_amplified_hi, giersiepen_hi};
/// let base = giersiepen_hi(Pressure::from_base(100.0), Time::from_base(1.0));
/// let amplified = cavitation_amplified_hi(base, 1.0);
/// assert!((amplified - base * 4.0).abs() < 1e-15);  // 1 + 3×1.0 = 4×
/// ```
#[inline]
#[must_use]
pub fn cavitation_amplified_hi(base_hi: f64, cav_potential: f64) -> f64 {
    cfd_core::physics::hemolysis::HemolysisModel::cavitation_amplified(base_hi, cav_potential)
}

/// Taskin (2012) calibration constant `C_T = 1.228 × 10⁻⁵`.
pub const TASKIN_C: f64 = 1.228e-5;

/// Taskin (2012) shear stress exponent `β_T = 1.9918`.
///
/// Slightly different from the Giersiepen value of 1.991, reflecting the
/// re-calibration against Lagrangian particle tracking data.
pub const TASKIN_BETA: f64 = 1.9918;

/// Taskin (2012) strain-based haemolysis index (single-segment evaluation).
///
/// Computes the haemolysis index using the Taskin integral-form model for a
/// single channel segment with constant shear stress:
///
/// ```text
/// HI_Taskin = C_T · τ^β_T · t
/// ```
///
/// This is the single-segment (constant shear) evaluation of the full
/// integral form:
///
/// ```text
/// HI_Taskin = C_T · ∫ τ(t)^β_T dt
/// ```
///
/// ## Theorem: Strain-Rate Path Dependence (Taskin et al. 2012)
///
/// The Taskin model captures cumulative damage along the flow path by integrating
/// the instantaneous shear-stress contribution over time, rather than using a
/// single power-law evaluation as in Giersiepen (1990). This integral form
/// correctly accounts for varying shear exposure histories and yields more
/// accurate hemolysis predictions for complex flow geometries where shear stress
/// varies along particle trajectories.
///
/// For constant shear stress, the integral reduces to `C_T · τ^β_T · t`, which
/// differs from the Giersiepen model (`C · t^α · τ^β`) in that the time
/// dependence is linear rather than sub-linear (`α = 0.765` in Giersiepen).
///
/// ## Reference
///
/// Taskin, M. E. et al. (2012). Evaluation of Eulerian and Lagrangian Models
/// for Hemolysis Estimation. *ASAIO J.*, 58(4), 363–372.
///
/// # Arguments
///
/// * `shear_stress` — wall shear stress \[Pa]
/// * `exposure_time` — exposure duration \[s]
///
/// # Returns
///
/// Haemolysis index (dimensionless). Returns `0.0` for non-positive inputs.
#[inline]
#[must_use]
pub fn taskin_hi(shear_stress: Pressure, exposure_time: Time) -> f64 {
    let shear_pa = shear_stress.into_base();
    let duration_s = exposure_time.into_base();
    if shear_pa <= 0.0 || duration_s <= 0.0 {
        return 0.0;
    }
    TASKIN_C * shear_pa.abs().max(1e-30).powf(TASKIN_BETA) * duration_s
}

/// Activation rate constant for Hematoporphyrin \[s⁻¹\] (Rosenthal 2004).
pub const SENSITIZER_K_ACT_HEMATOPORPHYRIN: f64 = 0.5;

/// Activation rate constant for Chlorin e6 \[s⁻¹\].
pub const SENSITIZER_K_ACT_CHLORIN_E6: f64 = 0.8;

/// Sonosensitizer activation efficiency as a function of transit time.
///
/// ## Theorem — First-Order Activation Kinetics (Rosenthal 2004)
///
/// The fraction of sonosensitizer activated during transit through a
/// cavitation zone is governed by first-order kinetics:
///
/// ```text
/// η_act = 1 − exp(−k_act · I_cav · t_transit)
/// ```
///
/// where:
/// - k_act ≈ 0.1–1.0 s⁻¹ (activation rate constant, drug-dependent)
/// - I_cav = cavitation intensity (dimensionless, 0–1)
/// - t_transit = residence time in cavitation zone \[s\]
///
/// For short transits (k·I·t << 1): η_act ≈ k·I·t (linear regime)
/// For long transits (k·I·t >> 1): η_act → 1 (saturation)
///
/// **Reference**: Rosenthal, I., Sostaric, J.Z. & Riesz, P. (2004).
/// "Sonodynamic therapy – a review of the synergistic effects of drugs
/// and ultrasound", *Ultrason. Sonochem.* 11:349-363.
///
/// # Arguments
///
/// * `k_act` — activation rate constant \[s⁻¹\]
/// * `cavitation_intensity` — I_cav (0–1)
/// * `transit_time_s` — residence time in cavitation zone \[s\]
///
/// # Returns
///
/// Activation fraction η_act ∈ \[0, 1\].
///
/// # Example
///
/// ```
/// use cfd_1d::sonosensitizer_activation_efficiency;
/// let eta = sonosensitizer_activation_efficiency(0.5, 0.8, 2.0);
/// assert!(eta > 0.0 && eta < 1.0);
/// ```
#[inline]
#[must_use]
pub fn sonosensitizer_activation_efficiency(
    k_act: f64,
    cavitation_intensity: f64,
    transit_time_s: f64,
) -> f64 {
    (1.0 - (-k_act * cavitation_intensity * transit_time_s)
        .max(-500.0)
        .exp())
    .clamp(0.0, 1.0)
}
