//! Fast reduced-order venturi screening helpers.
//!
//! These utilities intentionally stay 1D and inexpensive. They are used for
//! throat-level cavitation and diffuser-recovery screening, not full bubble
//! dynamics or multiphase CFD.

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

/// Convert a half-angle taper into an equivalent taper length [m].
#[must_use]
pub fn venturi_taper_length_m(outer_width_m: f64, throat_width_m: f64, half_angle_deg: f64) -> f64 {
    let delta = ((outer_width_m - throat_width_m) * 0.5).abs();
    let angle_rad = half_angle_deg
        .to_radians()
        .clamp(1e-6, std::f64::consts::FRAC_PI_2 - 1e-6);
    (delta / angle_rad.tan()).max(delta.max(1e-6))
}

/// Reduced-order discharge coefficient from convergent half-angle.
#[must_use]
pub fn discharge_coefficient_from_convergent_half_angle_deg(half_angle_deg: f64) -> f64 {
    let angle = half_angle_deg.clamp(1.0, 45.0);
    if angle <= 7.0 {
        0.985
    } else if angle <= 15.0 {
        0.975 - 0.002 * (angle - 7.0)
    } else if angle <= 25.0 {
        0.959 - 0.003 * (angle - 15.0)
    } else {
        0.929 - 0.002 * (angle - 25.0)
    }
    .clamp(0.88, 0.995)
}

/// Evaluate venturi throat screening quantities from local 1D state.
#[must_use]
pub fn evaluate_venturi_screening(input: VenturiScreeningInput) -> VenturiScreeningResult {
    // Reynolds-dependent discharge coefficient correction.
    // At low Re, viscous boundary layers in the convergent section
    // reduce the effective throat velocity beyond what the inviscid
    // vena-contracta coefficient predicts.
    //
    // - Re > 2300 (turbulent): C_d = base vena-contracta coeff (ISO 5167)
    // - Re 500-2300 (laminar): C_d reduced linearly to 92% of base
    // - Re < 500 (creeping): C_d reduced to 85% of base
    let v_raw = input.throat_velocity_m_s / input.vena_contracta_coeff.max(1e-9);
    let re_raw = input.density_kg_m3 * v_raw * input.throat_hydraulic_diameter_m.max(1e-18)
        / input.viscosity_pa_s.max(1e-18);
    let cd_correction = if re_raw > 2300.0 {
        1.0
    } else if re_raw > 500.0 {
        0.92 + 0.08 * (re_raw - 500.0) / 1800.0
    } else {
        0.85 + 0.07 * (re_raw / 500.0)
    };
    let effective_cc = input.vena_contracta_coeff * cd_correction;
    let v_eff = input.throat_velocity_m_s / effective_cc.max(1e-9);

    let bernoulli_drop = 0.5
        * input.density_kg_m3
        * (v_eff * v_eff - input.upstream_velocity_m_s * input.upstream_velocity_m_s).max(0.0);

    let re_throat = input.density_kg_m3 * v_eff * input.throat_hydraulic_diameter_m.max(1e-18)
        / input.viscosity_pa_s.max(1e-18);
    let f_darcy = if re_throat > 1.0 {
        64.0 / re_throat
    } else {
        64.0
    };
    let friction_drop = f_darcy
        * (input.throat_length_m / input.throat_hydraulic_diameter_m.max(1e-18))
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
    let effective_vapor_pressure_pa = input.vapor_pressure_pa + input.upstream_nuclei_fraction * 10_000.0;

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
    let throat_transit_s = input.throat_length_m / v_eff.max(1e-9);
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

    VenturiScreeningResult {
        effective_throat_velocity_m_s: v_eff,
        bernoulli_drop_pa: bernoulli_drop,
        throat_friction_drop_pa: friction_drop,
        throat_static_pressure_pa: throat_static,
        cavitation_number,
        diffuser_recovery_pa: diffuser_recovery,
        outlet_nuclei_fraction,
    }
}
