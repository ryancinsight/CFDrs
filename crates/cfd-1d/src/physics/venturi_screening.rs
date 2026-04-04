//! Fast reduced-order venturi screening helpers.
//!
//! These utilities intentionally stay 1D and inexpensive. They are used for
//! throat-level cavitation and diffuser-recovery screening, not full bubble
//! dynamics or multiphase CFD.

use cfd_core::error::{Error, Result};
use serde::{Deserialize, Serialize};

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
    /// Bernoulli contraction drop from upstream section to vena contracta [Pa].
    pub bernoulli_drop_pa: f64,
    /// Darcy-Weisbach friction drop through the throat [Pa].
    pub throat_friction_drop_pa: f64,
    /// Estimated static pressure at the throat [Pa].
    pub throat_static_pressure_pa: f64,
    /// Cavitation number based on throat static pressure and dynamic head.
    pub cavitation_number: f64,
    /// Diffuser static-pressure recovery estimate [Pa].
    pub diffuser_recovery_pa: f64,
    /// Downstream bubble nuclei volume fraction after generation.
    pub outlet_nuclei_fraction: f64,
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
    if input.vena_contracta_coeff <= 0.0 || input.diffuser_recovery_coeff < 0.0 {
        return Err(Error::InvalidConfiguration(
            "Venturi recovery coefficients must be nonnegative and the vena-contracta coefficient must be positive".to_string(),
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

/// Evaluate venturi throat screening quantities from local 1D state.
#[must_use]
pub fn evaluate_venturi_screening(input: VenturiScreeningInput) -> Result<VenturiScreeningResult> {
    validate_venturi_screening_input(&input)?;

    if input.throat_velocity_m_s == 0.0 {
        return Ok(VenturiScreeningResult {
            effective_throat_velocity_m_s: 0.0,
            bernoulli_drop_pa: 0.0,
            throat_friction_drop_pa: 0.0,
            throat_static_pressure_pa: input.upstream_pressure_pa,
            cavitation_number: f64::INFINITY,
            diffuser_recovery_pa: 0.0,
            outlet_nuclei_fraction: input.upstream_nuclei_fraction,
        });
    }

    // Reynolds-dependent discharge coefficient correction.
    // At low Re, viscous boundary layers in the convergent section
    // reduce the effective throat velocity beyond what the inviscid
    // vena-contracta coefficient predicts.
    //
    // - Re > 2300 (turbulent): C_d = base vena-contracta coeff (ISO 5167)
    // - Re 500-2300 (laminar): C_d reduced linearly to 92% of base
    // - Re < 500 (creeping): C_d reduced to 85% of base
    let v_raw = input.throat_velocity_m_s / input.vena_contracta_coeff;
    let re_raw = input.density_kg_m3 * v_raw * input.throat_hydraulic_diameter_m
        / input.viscosity_pa_s;
    let cd_correction = if re_raw > 2300.0 {
        1.0
    } else if re_raw > 500.0 {
        0.92 + 0.08 * (re_raw - 500.0) / 1800.0
    } else {
        0.85 + 0.07 * (re_raw / 500.0)
    };
    let effective_cc = input.vena_contracta_coeff * cd_correction;
    let v_eff = input.throat_velocity_m_s / effective_cc;

    let bernoulli_drop = 0.5
        * input.density_kg_m3
        * (v_eff * v_eff - input.upstream_velocity_m_s * input.upstream_velocity_m_s).max(0.0);

    let re_throat = input.density_kg_m3 * v_eff * input.throat_hydraulic_diameter_m
        / input.viscosity_pa_s;
    let f_darcy = if re_throat > 1.0 {
        64.0 / re_throat
    } else {
        64.0
    };
    let friction_drop = f_darcy
        * (input.throat_length_m / input.throat_hydraulic_diameter_m)
        * 0.5
        * input.density_kg_m3
        * v_eff
        * v_eff;

    // Unclamped static pressure preserves full Bernoulli + friction
    // differentiation in the sigma formula.  The physical (non-negative)
    // value is stored separately in `throat_static_pressure_pa`.
    let throat_static_raw = input.upstream_pressure_pa - bernoulli_drop - friction_drop;
    let throat_static = throat_static_raw.max(0.0);
    let dyn_p = 0.5 * input.density_kg_m3 * v_eff * v_eff;

    // Cascade effect: upstream nuclei elevate the effective vapor pressure
    // making inception more likely.
    let effective_vapor_pressure_pa =
        input.vapor_pressure_pa + input.upstream_nuclei_fraction * 10_000.0;

    let cavitation_number = if dyn_p > 1e-12 {
        (throat_static_raw - effective_vapor_pressure_pa) / dyn_p
    } else {
        f64::INFINITY
    };

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

    // Diffuser pressure recovery with Re-dependent efficiency.
    // At low Re (< 500), the boundary layer fills more of the diffuser
    // cross-section, reducing the pressure recovery.  At high Re (> 2000),
    // separation in the diffuser limits recovery to ~60-80% of ideal.
    //
    // Re-dependent correction based on Idelchik (2007), Diagram 5-2:
    // - Re > 2000: eta = base coefficient (turbulent, ~0.7-0.8)
    // - Re 200-2000: eta decreases linearly to 60% of base
    // - Re < 200: eta = 50% of base (viscous losses dominate)
    let re_correction = if re_throat > 2000.0 {
        1.0
    } else if re_throat > 200.0 {
        0.60 + 0.40 * (re_throat - 200.0) / 1800.0
    } else {
        0.50 + 0.10 * (re_throat / 200.0)
    };
    let diffuser_recovery = input.diffuser_recovery_coeff
        * re_correction
        * 0.5
        * input.density_kg_m3
        * (v_eff * v_eff - input.upstream_velocity_m_s * input.upstream_velocity_m_s).max(0.0);

    Ok(VenturiScreeningResult {
        effective_throat_velocity_m_s: v_eff,
        bernoulli_drop_pa: bernoulli_drop,
        throat_friction_drop_pa: friction_drop,
        throat_static_pressure_pa: throat_static,
        cavitation_number,
        diffuser_recovery_pa: diffuser_recovery,
        outlet_nuclei_fraction,
    })
}

#[cfg(test)]
mod tests {
    use super::{
        discharge_coefficient_from_convergent_half_angle_deg, evaluate_venturi_screening,
        venturi_taper_length_m, VenturiScreeningInput,
    };

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
    }
}
