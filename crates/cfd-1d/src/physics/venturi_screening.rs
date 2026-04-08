//! Fast reduced-order venturi screening helpers.
//!
//! These utilities intentionally stay 1D and inexpensive. They are used for
//! throat-level cavitation and diffuser-recovery screening, not full bubble
//! dynamics or multiphase CFD.
//!
//! # Theorem - Reduced-Order Venturi Throat Balance
//!
//! For an admissible venturi operating point with positive throat geometry,
//! positive fluid properties, and nonnegative section velocities with
//! $u_t \ge u_u$, the reduced-order throat balance can be written as
//!
//! ```text
//! p_t,raw = p_u - Δp_B - Δp_f
//! Δp_B = 1/2 ρ max(u_eff^2 - u_u^2, 0)
//! Δp_f = f_D (L_t / D_h) (1/2 ρ u_eff^2)
//! ```
//!
//! where $u_eff = u_t / C_c,eff$ is the vena-contracta-corrected throat
//! velocity, $Δp_B$ is the Bernoulli contraction drop, and $Δp_f$ is the
//! Darcy-Weisbach throat friction drop. Physical reporting uses
//! $p_t = max(p_t,raw, 0)$, but cavitation inception and margin are evaluated
//! against the unclamped $p_t,raw$ so that sub-vapor excursions are not hidden.
//!
//! **Proof sketch**: the implementation evaluates the throat velocity after the
//! Reynolds-dependent discharge correction, then applies Bernoulli between the
//! upstream section and the effective throat while adding the Darcy-Weisbach
//! loss over the declared throat length. The reported static pressure is a
//! nonnegative projection of the raw balance, so physical output remains valid
//! without erasing the sign information required by the cavitation number.

use crate::physics::resistance::models::VenturiModel;
use cfd_core::error::{Error, Result};
use serde::{Deserialize, Serialize};
use std::fmt;

/// Input state for a throat-level venturi screening evaluation.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct VenturiScreeningInput {
    /// Absolute upstream static pressure [Pa].
    pub upstream_pressure_pa: f64,
    /// Mean velocity in the upstream section [m/s].
    pub upstream_velocity_m_s: f64,
    /// Mean velocity in the geometric throat [m/s].
    pub throat_velocity_m_s: f64,
    /// Hydraulic diameter of the throat [m].
    pub throat_hydraulic_diameter_m: f64,
    /// Throat length [m].
    pub throat_length_m: f64,
    /// Fluid density [kg/m^3].
    pub density_kg_m3: f64,
    /// Dynamic viscosity [Pa·s].
    pub viscosity_pa_s: f64,
    /// Vapor pressure [Pa].
    pub vapor_pressure_pa: f64,
    /// Vena-contracta coefficient.
    pub vena_contracta_coeff: f64,
    /// Diffuser recovery coefficient.
    pub diffuser_recovery_coeff: f64,
    /// Upstream bubble nuclei volume fraction.
    #[serde(default)]
    pub upstream_nuclei_fraction: f64,
}

/// Output of a throat-level venturi screening evaluation.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct VenturiScreeningResult {
    /// Velocity corrected for vena-contracta acceleration [m/s].
    pub effective_throat_velocity_m_s: f64,
    /// Raw static pressure at the throat before clamping [Pa].
    #[serde(default)]
    pub raw_throat_static_pressure_pa: f64,
    /// Bernoulli contraction drop from upstream section to vena contracta [Pa].
    pub bernoulli_drop_pa: f64,
    /// Darcy-Weisbach friction drop through the throat [Pa].
    pub throat_friction_drop_pa: f64,
    /// Estimated static pressure at the throat [Pa].
    pub throat_static_pressure_pa: f64,
    /// Effective vapor-pressure threshold after nuclei amplification [Pa].
    #[serde(default)]
    pub effective_vapor_pressure_pa: f64,
    /// Cavitation number based on throat static pressure and dynamic head.
    pub cavitation_number: f64,
    /// Cavitation margin in pressure units: throat static pressure minus effective vapor pressure [Pa].
    #[serde(default)]
    pub cavitation_margin_pa: f64,
    /// Reynolds number computed from the raw throat velocity estimate.
    #[serde(default)]
    pub raw_throat_reynolds_number: f64,
    /// Reynolds number computed from the corrected throat velocity.
    #[serde(default)]
    pub effective_throat_reynolds_number: f64,
    /// Reynolds-number correction applied to the vena-contracta coefficient.
    #[serde(default)]
    pub discharge_coefficient_correction: f64,
    /// Effective vena-contracta coefficient after Reynolds correction.
    #[serde(default)]
    pub effective_vena_contracta_coeff: f64,
    /// Diffuser static-pressure recovery estimate [Pa].
    pub diffuser_recovery_pa: f64,
    /// Downstream bubble nuclei volume fraction after generation.
    pub outlet_nuclei_fraction: f64,
}

/// Screening classification for venturi throat operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum VenturiScreeningRisk {
    /// No meaningful flow is present.
    Quiescent,
    /// Cavitation margin is comfortably positive.
    Stable,
    /// Operation is close enough to inception that a design review is warranted.
    Watch,
    /// Cavitation inception is likely at the current operating point.
    CavitationLikely,
    /// Static pressure has already dropped below the effective vapor pressure.
    CavitationCritical,
}

impl fmt::Display for VenturiScreeningRisk {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Quiescent => write!(f, "Quiescent"),
            Self::Stable => write!(f, "Stable"),
            Self::Watch => write!(f, "Watch"),
            Self::CavitationLikely => write!(f, "CavitationLikely"),
            Self::CavitationCritical => write!(f, "CavitationCritical"),
        }
    }
}

/// Derived assessment for a venturi screening result.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct VenturiScreeningAssessment {
    /// Screening classification.
    pub risk: VenturiScreeningRisk,
    /// Pressure head above the effective vapor threshold [Pa].
    pub cavitation_margin_pa: f64,
    /// Diffuser recovery divided by Bernoulli contraction drop.
    pub pressure_recovery_ratio: f64,
}

fn validate_venturi_screening_input(input: &VenturiScreeningInput) -> Result<()> {
    let finite = |value: f64| value.is_finite();
    if !finite(input.upstream_pressure_pa)
        || !finite(input.upstream_velocity_m_s)
        || !finite(input.throat_velocity_m_s)
        || !finite(input.throat_hydraulic_diameter_m)
        || !finite(input.throat_length_m)
        || !finite(input.density_kg_m3)
        || !finite(input.viscosity_pa_s)
        || !finite(input.vapor_pressure_pa)
        || !finite(input.vena_contracta_coeff)
        || !finite(input.diffuser_recovery_coeff)
        || !finite(input.upstream_nuclei_fraction)
    {
        return Err(Error::InvalidConfiguration(
            "Venturi screening inputs must be finite".to_string(),
        ));
    }

    if input.upstream_pressure_pa < 0.0 {
        return Err(Error::InvalidConfiguration(
            "Venturi upstream pressure must be nonnegative".to_string(),
        ));
    }
    if input.upstream_velocity_m_s < 0.0 || input.throat_velocity_m_s < 0.0 {
        return Err(Error::InvalidConfiguration(
            "Venturi velocities must be nonnegative".to_string(),
        ));
    }
    if input.throat_velocity_m_s > 0.0 && input.throat_velocity_m_s < input.upstream_velocity_m_s {
        return Err(Error::InvalidConfiguration(
            "Venturi screening assumes the geometric throat accelerates relative to the upstream section".to_string(),
        ));
    }
    if input.throat_hydraulic_diameter_m <= 0.0 || input.throat_length_m <= 0.0 {
        return Err(Error::InvalidConfiguration(
            "Venturi throat geometry must be positive".to_string(),
        ));
    }
    if input.density_kg_m3 <= 0.0 || input.viscosity_pa_s <= 0.0 {
        return Err(Error::InvalidConfiguration(
            "Venturi fluid properties must be positive".to_string(),
        ));
    }
    if input.vapor_pressure_pa < 0.0 {
        return Err(Error::InvalidConfiguration(
            "Venturi vapor pressure must be nonnegative".to_string(),
        ));
    }
    if input.vena_contracta_coeff <= 0.0
        || input.vena_contracta_coeff > 1.0
        || input.diffuser_recovery_coeff < 0.0
        || input.diffuser_recovery_coeff > 1.0
    {
        return Err(Error::InvalidConfiguration(
            "Venturi recovery coefficients must lie in [0, 1] and the vena-contracta coefficient must lie in (0, 1]".to_string(),
        ));
    }
    if !(0.0..=1.0).contains(&input.upstream_nuclei_fraction) {
        return Err(Error::InvalidConfiguration(
            "Venturi upstream nuclei fraction must lie in [0, 1]".to_string(),
        ));
    }

    Ok(())
}

/// Convert a half-angle taper into an equivalent taper length [m].
///
/// # Theorem
///
/// For fixed outer and throat widths with $0 < w_t < w_o$ and a physical
/// convergent half-angle $\theta \in (0, \pi/2)$,
///
/// ```text
/// L = (w_o - w_t) / (2 tan θ)
/// ```
///
/// is strictly positive and strictly decreases as $\theta$ increases.
///
/// **Proof sketch**: the numerator is positive by construction, and
/// $\tan(\theta)$ is positive and strictly increasing on $(0, \pi/2)$, so its
/// reciprocal is strictly decreasing.
#[must_use]
pub fn venturi_taper_length_m(
    outer_width_m: f64,
    throat_width_m: f64,
    half_angle_deg: f64,
) -> Result<f64> {
    if !outer_width_m.is_finite()
        || !throat_width_m.is_finite()
        || !half_angle_deg.is_finite()
    {
        return Err(Error::InvalidConfiguration(
            "Venturi taper geometry must be finite".to_string(),
        ));
    }
    if outer_width_m <= 0.0 || throat_width_m <= 0.0 {
        return Err(Error::InvalidConfiguration(
            "Venturi widths must be positive".to_string(),
        ));
    }
    if throat_width_m >= outer_width_m {
        return Err(Error::InvalidConfiguration(
            "Venturi throat width must be less than the outer width".to_string(),
        ));
    }
    if half_angle_deg <= 0.0 || half_angle_deg >= 90.0 {
        return Err(Error::InvalidConfiguration(
            "Venturi half-angle must lie in the open interval (0, 90) degrees".to_string(),
        ));
    }

    let delta = (outer_width_m - throat_width_m) * 0.5;
    let angle_rad = half_angle_deg.to_radians();
    Ok(delta / angle_rad.tan())
}

/// Reduced-order discharge coefficient from convergent half-angle.
///
/// # Theorem
///
/// On the calibrated domain $\theta \in [1^\circ, 45^\circ]$, the returned
/// piecewise fit is bounded in $(0, 1)$ and is monotone nonincreasing in the
/// half-angle, representing the loss of contraction efficiency as the
/// convergent becomes more abrupt.
///
/// **Proof sketch**: each branch is affine with a nonpositive slope and the
/// branch endpoints are ordered so the piecewise fit does not increase across
/// interval boundaries; the final clamp preserves the same bound and monotonic
/// direction.
#[must_use]
pub fn discharge_coefficient_from_convergent_half_angle_deg(half_angle_deg: f64) -> Result<f64> {
    if !half_angle_deg.is_finite() || !(1.0..=45.0).contains(&half_angle_deg) {
        return Err(Error::InvalidConfiguration(
            "Venturi convergent half-angle must lie in [1, 45] degrees".to_string(),
        ));
    }

    let angle = half_angle_deg;
    let coefficient = if angle <= 7.0 {
        0.985
    } else if angle <= 15.0 {
        0.975 - 0.002 * (angle - 7.0)
    } else if angle <= 25.0 {
        0.959 - 0.003 * (angle - 15.0)
    } else {
        0.929 - 0.002 * (angle - 25.0)
    };

    Ok(coefficient.clamp(0.88, 0.995))
}

/// Classify a venturi screening result into a qualitative risk band.
#[must_use]
pub fn classify_venturi_screening(result: &VenturiScreeningResult) -> VenturiScreeningRisk {
    if result.effective_throat_velocity_m_s == 0.0 {
        VenturiScreeningRisk::Quiescent
    } else if result.cavitation_margin_pa < 0.0 {
        VenturiScreeningRisk::CavitationCritical
    } else if result.cavitation_number < 1.0 {
        VenturiScreeningRisk::CavitationLikely
    } else if result.cavitation_number < 2.0 || result.outlet_nuclei_fraction > 0.05 {
        VenturiScreeningRisk::Watch
    } else {
        VenturiScreeningRisk::Stable
    }
}

/// Build a compact qualitative assessment from a screening result.
#[must_use]
pub fn assess_venturi_screening(result: &VenturiScreeningResult) -> VenturiScreeningAssessment {
    let pressure_recovery_ratio = if result.bernoulli_drop_pa > 0.0 {
        result.diffuser_recovery_pa / result.bernoulli_drop_pa
    } else {
        0.0
    };

    VenturiScreeningAssessment {
        risk: classify_venturi_screening(result),
        cavitation_margin_pa: result.cavitation_margin_pa,
        pressure_recovery_ratio,
    }
}

/// Evaluate venturi throat screening quantities from local 1D state.
#[must_use]
pub fn evaluate_venturi_screening(input: VenturiScreeningInput) -> Result<VenturiScreeningResult> {
    validate_venturi_screening_input(&input)?;

    let effective_vapor_pressure_pa =
        cfd_core::physics::cavitation::nuclei_transport::nuclei_adjusted_vapor_pressure(
            input.vapor_pressure_pa,
            input.upstream_nuclei_fraction,
        );

    if input.throat_velocity_m_s == 0.0 {
        return Ok(VenturiScreeningResult {
            effective_throat_velocity_m_s: 0.0,
            raw_throat_static_pressure_pa: input.upstream_pressure_pa,
            bernoulli_drop_pa: 0.0,
            throat_friction_drop_pa: 0.0,
            throat_static_pressure_pa: input.upstream_pressure_pa,
            effective_vapor_pressure_pa,
            cavitation_number: f64::INFINITY,
            cavitation_margin_pa: input.upstream_pressure_pa - effective_vapor_pressure_pa,
            raw_throat_reynolds_number: 0.0,
            effective_throat_reynolds_number: 0.0,
            discharge_coefficient_correction: 1.0,
            effective_vena_contracta_coeff: input.vena_contracta_coeff,
            diffuser_recovery_pa: 0.0,
            outlet_nuclei_fraction: input.upstream_nuclei_fraction,
        });
    }

    // Reynolds-dependent discharge coefficient correction.
    // Reuse the venturi resistance-model correction so the reduced-order
    // screening path and the pressure-loss model degrade discharge
    // consistently in the millifluidic low-Re regime.
    let v_raw = input.throat_velocity_m_s / input.vena_contracta_coeff;
    let re_raw = input.density_kg_m3 * v_raw * input.throat_hydraulic_diameter_m
        / input.viscosity_pa_s;
    let cd_correction = VenturiModel::<f64>::discharge_coefficient_reynolds_correction(re_raw);
    let effective_cc = input.vena_contracta_coeff * cd_correction;
    let v_eff = input.throat_velocity_m_s / effective_cc;

    let bernoulli_drop = 0.5
        * input.density_kg_m3
        * (v_eff * v_eff - input.upstream_velocity_m_s * input.upstream_velocity_m_s).max(0.0);

    let re_throat = input.density_kg_m3 * v_eff * input.throat_hydraulic_diameter_m
        / input.viscosity_pa_s;
    let f_darcy = if re_throat > 0.0 {
        VenturiModel::<f64>::friction_factor_for_reynolds(re_throat)
    } else {
        64.0
    };
    let friction_drop = VenturiModel::<f64>::throat_friction_pressure_drop(
        f_darcy,
        input.density_kg_m3,
        input.throat_length_m,
        input.throat_hydraulic_diameter_m,
        v_eff,
    );

    // Unclamped static pressure preserves full Bernoulli + friction
    // differentiation in the sigma formula.  The physical (non-negative)
    // value is stored separately in `throat_static_pressure_pa`.
    let throat_static_raw = input.upstream_pressure_pa - bernoulli_drop - friction_drop;
    let throat_static = throat_static_raw.max(0.0);
    let dyn_p = 0.5 * input.density_kg_m3 * v_eff * v_eff;

    let cavitation_number = if dyn_p > 1e-12 {
        (throat_static_raw - effective_vapor_pressure_pa) / dyn_p
    } else {
        f64::INFINITY
    };
    let cavitation_margin_pa = throat_static_raw - effective_vapor_pressure_pa;

    // Nuclei generation at the throat: use the physical model from
    // NucleiTransport.  Cavitation (sigma < 1) generates nuclei from
    // bubble collapse; the generation rate scales with the cavitation
    // intensity (1 - sigma).  The dissolution during throat transit
    // partially offsets generation.
    let throat_transit_s = input.throat_length_m / v_eff;
    let cavitation_source = (1.0 - cavitation_number).max(0.0);
    let transport = cfd_core::physics::cavitation::nuclei_transport::NucleiTransport::new(
        cfd_core::physics::cavitation::nuclei_transport::NucleiTransportConfig::default(),
    );
    let net_rate = transport.calculate_net_reaction_rate(
        input.upstream_nuclei_fraction,
        cavitation_source,
    );
    let outlet_nuclei_fraction =
        (input.upstream_nuclei_fraction + net_rate * throat_transit_s).clamp(0.0, 1.0);

    let re_correction = VenturiModel::<f64>::diffuser_recovery_reynolds_correction(re_throat);
    let diffuser_recovery = input.diffuser_recovery_coeff
        * re_correction
        * 0.5
        * input.density_kg_m3
        * (v_eff * v_eff - input.upstream_velocity_m_s * input.upstream_velocity_m_s).max(0.0);

    Ok(VenturiScreeningResult {
        effective_throat_velocity_m_s: v_eff,
        raw_throat_static_pressure_pa: throat_static_raw,
        bernoulli_drop_pa: bernoulli_drop,
        throat_friction_drop_pa: friction_drop,
        throat_static_pressure_pa: throat_static,
        effective_vapor_pressure_pa,
        cavitation_number,
        cavitation_margin_pa,
        raw_throat_reynolds_number: re_raw,
        effective_throat_reynolds_number: re_throat,
        discharge_coefficient_correction: cd_correction,
        effective_vena_contracta_coeff: effective_cc,
        diffuser_recovery_pa: diffuser_recovery,
        outlet_nuclei_fraction,
    })
}

#[cfg(test)]
mod tests {
    use super::{
        assess_venturi_screening, classify_venturi_screening,
        discharge_coefficient_from_convergent_half_angle_deg, evaluate_venturi_screening,
        venturi_taper_length_m, VenturiScreeningInput, VenturiScreeningRisk,
    };
    use crate::physics::resistance::models::VenturiModel;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn venturi_taper_length_decreases_with_half_angle(
            outer_width_m in 2.0e-4_f64..2.0e-2,
            throat_width_ratio in 0.1_f64..0.95,
            small_half_angle_deg in 1.0_f64..44.0,
            angle_increment_deg in 1.0e-3_f64..1.0,
        ) {
            let throat_width_m = outer_width_m * throat_width_ratio;
            let large_half_angle_deg = (small_half_angle_deg + angle_increment_deg).min(45.0);

            prop_assume!(large_half_angle_deg > small_half_angle_deg);

            let long_taper = venturi_taper_length_m(
                outer_width_m,
                throat_width_m,
                small_half_angle_deg,
            )
            .expect("physical taper geometry should evaluate");
            let short_taper = venturi_taper_length_m(
                outer_width_m,
                throat_width_m,
                large_half_angle_deg,
            )
            .expect("physical taper geometry should evaluate");

            prop_assert!(long_taper.is_sign_positive());
            prop_assert!(short_taper.is_sign_positive());
            prop_assert!(long_taper > short_taper);
        }
    }

    proptest! {
        #[test]
        fn discharge_coefficient_is_monotone_nonincreasing_in_half_angle(
            small_half_angle_deg in 1.0_f64..45.0,
            large_half_angle_deg in 1.0_f64..45.0,
        ) {
            prop_assume!(small_half_angle_deg <= large_half_angle_deg);

            let small = discharge_coefficient_from_convergent_half_angle_deg(small_half_angle_deg)
                .expect("angle is inside calibration range");
            let large = discharge_coefficient_from_convergent_half_angle_deg(large_half_angle_deg)
                .expect("angle is inside calibration range");

            prop_assert!(small >= large);
            prop_assert!((0.88..=0.995).contains(&small));
            prop_assert!((0.88..=0.995).contains(&large));
        }
    }

    #[test]
    fn venturi_taper_length_rejects_invalid_geometry() {
        assert!(venturi_taper_length_m(1.0, 1.0, 7.0).is_err());
        assert!(venturi_taper_length_m(1.0, 0.5, 0.0).is_err());
    }

    #[test]
    fn discharge_coefficient_rejects_out_of_domain_angles() {
        assert!(discharge_coefficient_from_convergent_half_angle_deg(0.5).is_err());
        assert!(discharge_coefficient_from_convergent_half_angle_deg(46.0).is_err());
    }

    #[test]
    fn evaluate_venturi_screening_rejects_nonphysical_inputs() {
        let input = VenturiScreeningInput {
            upstream_pressure_pa: -1.0,
            upstream_velocity_m_s: 0.1,
            throat_velocity_m_s: 0.2,
            throat_hydraulic_diameter_m: 1.0e-3,
            throat_length_m: 1.0e-3,
            density_kg_m3: 1060.0,
            viscosity_pa_s: 3.5e-3,
            vapor_pressure_pa: 3170.0,
            vena_contracta_coeff: 0.9,
            diffuser_recovery_coeff: 0.7,
            upstream_nuclei_fraction: 0.0,
        };

        assert!(evaluate_venturi_screening(input).is_err());
    }

    #[test]
    fn evaluate_venturi_screening_rejects_out_of_regime_coefficients() {
        let input = VenturiScreeningInput {
            upstream_pressure_pa: 101_325.0,
            upstream_velocity_m_s: 0.1,
            throat_velocity_m_s: 0.2,
            throat_hydraulic_diameter_m: 1.0e-3,
            throat_length_m: 1.0e-3,
            density_kg_m3: 1060.0,
            viscosity_pa_s: 3.5e-3,
            vapor_pressure_pa: 3170.0,
            vena_contracta_coeff: 1.05,
            diffuser_recovery_coeff: 0.7,
            upstream_nuclei_fraction: 0.0,
        };

        assert!(evaluate_venturi_screening(input).is_err());

        let input = VenturiScreeningInput {
            diffuser_recovery_coeff: 1.1,
            vena_contracta_coeff: 0.95,
            ..input
        };
        assert!(evaluate_venturi_screening(input).is_err());
    }

    #[test]
    fn evaluate_venturi_screening_rejects_nonaccelerating_throat_state() {
        let input = VenturiScreeningInput {
            upstream_pressure_pa: 101_325.0,
            upstream_velocity_m_s: 0.25,
            throat_velocity_m_s: 0.15,
            throat_hydraulic_diameter_m: 1.0e-3,
            throat_length_m: 1.0e-3,
            density_kg_m3: 1060.0,
            viscosity_pa_s: 3.5e-3,
            vapor_pressure_pa: 3170.0,
            vena_contracta_coeff: 0.9,
            diffuser_recovery_coeff: 0.7,
            upstream_nuclei_fraction: 0.0,
        };

        assert!(evaluate_venturi_screening(input).is_err());
    }

    #[test]
    fn evaluate_venturi_screening_handles_zero_throat_velocity() {
        let input = VenturiScreeningInput {
            upstream_pressure_pa: 101_325.0,
            upstream_velocity_m_s: 0.1,
            throat_velocity_m_s: 0.0,
            throat_hydraulic_diameter_m: 1.0e-3,
            throat_length_m: 1.0e-3,
            density_kg_m3: 1060.0,
            viscosity_pa_s: 3.5e-3,
            vapor_pressure_pa: 3170.0,
            vena_contracta_coeff: 0.9,
            diffuser_recovery_coeff: 0.7,
            upstream_nuclei_fraction: 0.0,
        };

        let result = evaluate_venturi_screening(input).expect("zero-flow limit should succeed");
        assert_eq!(result.effective_throat_velocity_m_s, 0.0);
        assert_eq!(result.bernoulli_drop_pa, 0.0);
        assert_eq!(result.throat_friction_drop_pa, 0.0);
        assert_eq!(result.diffuser_recovery_pa, 0.0);
        assert!(result.cavitation_number.is_infinite());
        assert_eq!(classify_venturi_screening(&result), VenturiScreeningRisk::Quiescent);

        let assessment = assess_venturi_screening(&result);
        assert_eq!(assessment.risk, VenturiScreeningRisk::Quiescent);
        assert_eq!(assessment.pressure_recovery_ratio, 0.0);
    }

    #[test]
    fn evaluate_venturi_screening_exposes_raw_and_effective_metrics() {
        let input = VenturiScreeningInput {
            upstream_pressure_pa: 101_325.0,
            upstream_velocity_m_s: 0.1,
            throat_velocity_m_s: 0.2,
            throat_hydraulic_diameter_m: 1.0e-3,
            throat_length_m: 1.0e-3,
            density_kg_m3: 1060.0,
            viscosity_pa_s: 3.5e-3,
            vapor_pressure_pa: 3170.0,
            vena_contracta_coeff: 0.9,
            diffuser_recovery_coeff: 0.7,
            upstream_nuclei_fraction: 0.0,
        };

        let result = evaluate_venturi_screening(input).expect("screening should succeed");

        assert!(result.raw_throat_static_pressure_pa >= result.throat_static_pressure_pa);
        assert!(result.raw_throat_reynolds_number > 0.0);
        assert!(result.effective_throat_reynolds_number > 0.0);
        assert!(result.discharge_coefficient_correction <= 1.0);
        assert!(result.effective_vena_contracta_coeff > 0.0);
        assert_eq!(classify_venturi_screening(&result), VenturiScreeningRisk::Stable);

        let assessment = assess_venturi_screening(&result);
        assert_eq!(assessment.risk, VenturiScreeningRisk::Stable);
        assert!(assessment.cavitation_margin_pa > 0.0);
        assert!(assessment.pressure_recovery_ratio >= 0.0);
    }

    #[test]
    fn evaluate_venturi_screening_flags_cavitation_critical_operation() {
        let input = VenturiScreeningInput {
            upstream_pressure_pa: 4_000.0,
            upstream_velocity_m_s: 0.1,
            throat_velocity_m_s: 5.0,
            throat_hydraulic_diameter_m: 1.0e-3,
            throat_length_m: 1.0e-3,
            density_kg_m3: 1060.0,
            viscosity_pa_s: 3.5e-3,
            vapor_pressure_pa: 3170.0,
            vena_contracta_coeff: 0.9,
            diffuser_recovery_coeff: 0.7,
            upstream_nuclei_fraction: 0.0,
        };

        let result = evaluate_venturi_screening(input).expect("screening should succeed");

        assert!(result.raw_throat_static_pressure_pa < 0.0);
        assert_eq!(classify_venturi_screening(&result), VenturiScreeningRisk::CavitationCritical);
    }

    #[test]
    fn evaluate_venturi_screening_uses_turbulent_throat_friction_scaling() {
        let input = VenturiScreeningInput {
            upstream_pressure_pa: 101_325.0,
            upstream_velocity_m_s: 1.0,
            throat_velocity_m_s: 10.0,
            throat_hydraulic_diameter_m: 1.0e-3,
            throat_length_m: 2.0e-3,
            density_kg_m3: 1000.0,
            viscosity_pa_s: 1.0e-3,
            vapor_pressure_pa: 3170.0,
            vena_contracta_coeff: 0.98,
            diffuser_recovery_coeff: 0.7,
            upstream_nuclei_fraction: 0.0,
        };

        let result = evaluate_venturi_screening(input).expect("high-Re screening should succeed");

        assert!(result.effective_throat_reynolds_number > 2300.0);

        let expected_f = 0.3164 / result.effective_throat_reynolds_number.powf(0.25);
        let expected_drop = expected_f
            * (input.throat_length_m / input.throat_hydraulic_diameter_m)
            * 0.5
            * input.density_kg_m3
            * result.effective_throat_velocity_m_s.powi(2);
        let laminar_drop = (64.0 / result.effective_throat_reynolds_number)
            * (input.throat_length_m / input.throat_hydraulic_diameter_m)
            * 0.5
            * input.density_kg_m3
            * result.effective_throat_velocity_m_s.powi(2);

        assert!((result.throat_friction_drop_pa - expected_drop).abs() < 1e-9 * expected_drop.max(1.0));
        assert!(result.throat_friction_drop_pa > laminar_drop);
    }

    #[test]
    fn evaluate_venturi_screening_matches_shared_discharge_correction() {
        let input = VenturiScreeningInput {
            upstream_pressure_pa: 101_325.0,
            upstream_velocity_m_s: 0.02,
            throat_velocity_m_s: 0.04,
            throat_hydraulic_diameter_m: 1.0e-3,
            throat_length_m: 1.0e-3,
            density_kg_m3: 1000.0,
            viscosity_pa_s: 1.0e-3,
            vapor_pressure_pa: 3170.0,
            vena_contracta_coeff: 0.95,
            diffuser_recovery_coeff: 0.7,
            upstream_nuclei_fraction: 0.0,
        };

        let result = evaluate_venturi_screening(input).expect("screening should succeed");
        let expected_correction = (0.5
            + 0.5 * (result.raw_throat_reynolds_number / 1000.0).powf(0.3))
            .min(1.0);

        assert!((result.discharge_coefficient_correction - expected_correction).abs() < 1e-12);
        assert!((result.effective_vena_contracta_coeff - input.vena_contracta_coeff * expected_correction).abs() < 1e-12);
    }

    #[test]
    fn evaluate_venturi_screening_matches_shared_recovery_correction() {
        let input = VenturiScreeningInput {
            upstream_pressure_pa: 101_325.0,
            upstream_velocity_m_s: 0.02,
            throat_velocity_m_s: 0.04,
            throat_hydraulic_diameter_m: 1.0e-3,
            throat_length_m: 1.0e-3,
            density_kg_m3: 1000.0,
            viscosity_pa_s: 1.0e-3,
            vapor_pressure_pa: 3170.0,
            vena_contracta_coeff: 0.95,
            diffuser_recovery_coeff: 0.7,
            upstream_nuclei_fraction: 0.0,
        };

        let result = evaluate_venturi_screening(input).expect("screening should succeed");
        let expected_correction =
            VenturiModel::<f64>::diffuser_recovery_reynolds_correction(
                result.effective_throat_reynolds_number,
            );
        let dynamic_recovery = 0.5
            * input.density_kg_m3
            * (result.effective_throat_velocity_m_s.powi(2) - input.upstream_velocity_m_s.powi(2))
                .max(0.0);

        assert!(
            (result.diffuser_recovery_pa
                - input.diffuser_recovery_coeff * expected_correction * dynamic_recovery)
                .abs()
                < 1e-12
        );
    }
}
