use aequitas::systems::si::quantities::{
    DynamicViscosity, EnergyPerVolume, Length, Power, Pressure, ReciprocalLength, ReciprocalTime,
    TemperatureDifference, Time, Velocity, Volume, VolumetricFlowRate,
};
use aequitas::systems::si::units::{JoulePerCubicMeter, JoulePerMilliliter, Kelvin};
use cfd_1d::{
    acoustic_contrast_factor, acoustic_energy_density, cavitation_amplified_hi,
    cavitation_hemolysis_amplification, giersiepen_hi, sonosensitizer_activation_efficiency,
    KAPPA_CTC, KAPPA_PLASMA, KAPPA_RBC, RHO_CTC, RHO_PLASMA, RHO_RBC,
    SENSITIZER_K_ACT_HEMATOPORPHYRIN,
};
use cfd_schematics::topology::TreatmentActuationMode;
use hyperion::{
    coefficient::{InteractionCoefficient, LinearAttenuation},
    quantity::PathLength,
};

use super::report_math::{
    coefficient_of_variation, cumulative_pass_damage, direct_linear_risk, inverse_linear_risk,
    log_risk, mean, percentile, resonance_match, split_stage_flow_fractions,
};
use crate::constraints::{
    BLOOD_ATTENUATION_405NM_INV_M, BLOOD_DENSITY_KG_M3, BLOOD_VAPOR_PRESSURE_PA,
    BLOOD_VISCOSITY_PA_S, BUBBLE_POLYTROPIC_K, CLOTTING_BFR_CAUTION_ML_MIN,
    CLOTTING_BFR_HIGH_RISK_ML_MIN, CLOTTING_BFR_LOW_RISK_ML_MIN, CLOTTING_BFR_STRICT_10MLS_ML_MIN,
    CLOTTING_RESIDENCE_HIGH_RISK_S, CLOTTING_RESIDENCE_LOW_RISK_S, CLOTTING_SHEAR_HIGH_RISK_INV_S,
    CLOTTING_SHEAR_LOW_RISK_INV_S, C_P_BLOOD_J_KG_K, DEAD_VOLUME_SHEAR_THRESHOLD_INV_S,
    EXPANSION_RATIO_LOW_RISK, FDA_MAX_WALL_SHEAR_PA, FDA_THROAT_TEMP_RISE_LIMIT_K,
    FDA_TRANSIENT_SHEAR_PA, FDA_TRANSIENT_TIME_S, MILESTONE_TREATMENT_DURATION_MIN,
    PATIENT_BLOOD_VOLUME_ML, PEDIATRIC_BLOOD_VOLUME_ML_PER_KG, PEDIATRIC_FLOW_CAUTION_ML_MIN,
    PEDIATRIC_FLOW_EXCESSIVE_ML_MIN, PEDIATRIC_REFERENCE_WEIGHT_KG, PLATE_HEIGHT_MM,
    PLATE_WIDTH_MM, P_ATM_PA, SONO_REF_P_ABS_PA, THERAPEUTIC_WINDOW_REF,
    VENTURI_EXPANSION_RATIO_HIGH_RISK, VENTURI_VEL_RATIO_REF,
};
use crate::domain::BlueprintCandidate;
use crate::error::OptimError;
use crate::metrics::{
    compute_blueprint_separation_metrics, compute_blueprint_venturi_metrics,
    compute_typed_blueprint_safety_metrics, compute_typed_residence_metrics,
    healthy_cell_protection_index as compute_healthy_cell_protection_index,
    solve_blueprint_candidate, ChannelHemolysis, SdtMetrics,
};

const MILLIMETRES_PER_METRE: f64 = 1_000.0;
const MILLILITRES_PER_CUBIC_METRE: f64 = 1_000_000.0;
const SECONDS_PER_MINUTE: f64 = 60.0;

/// Unit-bearing report values before conversion to the legacy serialized DTO.
///
/// `SdtMetrics` is a persisted/reporting contract whose field names encode
/// display units. Keeping the physical values together here makes the
/// computation boundary dimension-safe while [`write_to`] remains the sole
/// explicit conversion into those display units.
#[derive(Debug, Clone, Copy)]
struct TypedReportPhysicalMetrics {
    throat_shear_rate: ReciprocalTime,
    throat_shear: Pressure,
    max_main_channel_shear: Pressure,
    mean_residence_time: Time,
    total_pressure_drop: Pressure,
    total_path_length: Length,
    total_ecv: Volume,
    flow_rate: VolumetricFlowRate,
    main_channel_shear_rate: ReciprocalTime,
    throat_transit_time: Time,
    optical_path_length: Length,
    safety_margin: Pressure,
    wall_shear_p95: Pressure,
    wall_shear_p99: Pressure,
    wall_shear_mean: Pressure,
    diffuser_recovery: Pressure,
    mechanical_power: Power,
    treatment_zone_dwell_time: Time,
    specific_cavitation_energy: EnergyPerVolume,
    acoustic_energy_density: EnergyPerVolume,
    throat_temperature_rise: TemperatureDifference,
}

#[derive(Debug, Clone, Copy)]
struct TypedChannelHemolysis<'a> {
    channel_id: &'a str,
    is_venturi_throat: bool,
    hi_contribution: f64,
    wall_shear: Pressure,
    transit_time: Time,
    flow_fraction: f64,
}

impl TypedChannelHemolysis<'_> {
    fn into_serialized(self) -> ChannelHemolysis {
        ChannelHemolysis {
            channel_id: self.channel_id.to_owned(),
            is_venturi_throat: self.is_venturi_throat,
            hi_contribution: self.hi_contribution,
            wall_shear_pa: self.wall_shear.into_base(),
            transit_time_s: self.transit_time.into_base(),
            flow_fraction: self.flow_fraction,
        }
    }
}

impl TypedReportPhysicalMetrics {
    /// Convert the typed computation result into the serialized report units.
    fn write_to(self, metrics: &mut SdtMetrics) {
        let Self {
            throat_shear_rate,
            throat_shear,
            max_main_channel_shear,
            mean_residence_time,
            total_pressure_drop,
            total_path_length,
            total_ecv,
            flow_rate,
            main_channel_shear_rate,
            throat_transit_time,
            optical_path_length,
            safety_margin,
            wall_shear_p95,
            wall_shear_p99,
            wall_shear_mean,
            diffuser_recovery,
            mechanical_power,
            treatment_zone_dwell_time,
            specific_cavitation_energy,
            acoustic_energy_density,
            throat_temperature_rise,
        } = self;

        metrics.throat_shear_rate_inv_s = throat_shear_rate.into_base();
        metrics.throat_shear_pa = throat_shear.into_base();
        metrics.max_main_channel_shear_pa = max_main_channel_shear.into_base();
        metrics.mean_residence_time_s = mean_residence_time.into_base();
        metrics.total_pressure_drop_pa = total_pressure_drop.into_base();
        metrics.total_path_length_mm = total_path_length.into_base() * MILLIMETRES_PER_METRE;
        metrics.total_ecv_ml = total_ecv.into_base() * MILLILITRES_PER_CUBIC_METRE;
        metrics.flow_rate_ml_min =
            flow_rate.into_base() * MILLILITRES_PER_CUBIC_METRE * SECONDS_PER_MINUTE;
        metrics.main_channel_shear_rate_inv_s = main_channel_shear_rate.into_base();
        metrics.throat_transit_time_s = throat_transit_time.into_base();
        metrics.optical_path_length_405_m = optical_path_length.into_base();
        metrics.safety_margin_pa = safety_margin.into_base();
        metrics.wall_shear_p95_pa = wall_shear_p95.into_base();
        metrics.wall_shear_p99_pa = wall_shear_p99.into_base();
        metrics.wall_shear_mean_pa = wall_shear_mean.into_base();
        metrics.diffuser_recovery_pa = diffuser_recovery.into_base();
        metrics.mechanical_power_w = mechanical_power.into_base();
        metrics.treatment_zone_dwell_time_s = treatment_zone_dwell_time.into_base();
        metrics.specific_cavitation_energy_j_ml =
            specific_cavitation_energy.in_unit::<JoulePerMilliliter>();
        metrics.acoustic_energy_density_j_m3 =
            acoustic_energy_density.in_unit::<JoulePerCubicMeter>();
        metrics.throat_temperature_rise_k = throat_temperature_rise.in_unit::<Kelvin>();
    }
}

fn cavitation_strength_from_sigma(cavitation_number: f64) -> f64 {
    let strength = (1.0 - cavitation_number).max(0.0);
    strength / (1.0 + strength)
}

fn blue_light_delivery_index(
    optical_path_length_m: f64,
    feed_hematocrit: f64,
) -> Result<f64, hyperion::TransportError<f64>> {
    let attenuation = InteractionCoefficient::<f64, LinearAttenuation>::new(
        ReciprocalLength::from_base(BLOOD_ATTENUATION_405NM_INV_M * feed_hematocrit.max(0.0)),
    )?;
    let path = PathLength::new(Length::from_base(optical_path_length_m))?;
    let transmission = attenuation.optical_depth(path)?.transmission();
    Ok(transmission.into_quantity().into_base())
}

fn mechanical_power(pressure_drop_pa: f64, flow_rate_m3_s: f64) -> Power {
    let pressure = Pressure::from_base(pressure_drop_pa);
    let flow_rate = VolumetricFlowRate::from_base(flow_rate_m3_s);
    pressure * flow_rate
}

/// Compute the full set of report-grade SDT metrics for a single blueprint candidate.
///
/// # NOW INTEGRATED acoustic / cavitation physics
///
/// The following validated cfd-1d physics models are called via
/// `compute_sdt_acoustic_metrics()` and their results feed directly into report
/// metrics:
///
/// - **Sonosensitizer activation kinetics (Rosenthal 2004)**: first-order model
///   η = 1 − exp(−k_act · I_cav · t_transit). Modulates `cancer_dose_fraction`
///   — short throat transits (< 1 ms) yield incomplete activation.
///
/// - **Rayleigh-Plesset collapse dynamics (Rayleigh 1917)**: hemolysis amplification
///   A = 1 + α·(R_max/R₀)²·(p/p_ref). Feeds into
///   `hemolysis_index_per_pass_cavitation_amplified` when cavitation_potential > 0.
///
/// - **Acoustic contrast factor / Gor'kov (1962)**: Φ = f₁/3 + f₂/2 for CTCs and
///   RBCs in plasma. Differential radiation force available for narrative reporting
///   via `sdt_acoustic.ctc_contrast_factor` and `sdt_acoustic.rbc_contrast_factor`.
///
/// - **Acoustic energy density (Gor'kov 1962)**: E_ac = p₀²/(4ρc²) at 100 kPa
///   pressure amplitude. Available for narrative reporting via
///   `sdt_acoustic.acoustic_energy_density_j_m3`.
///
/// # Available for 2D/3D refinement
///
/// The following models are implemented and validated but reserved for higher-fidelity
/// cascade pipelines:
///
/// - **Quemada (1978) viscosity**: rouleaux aggregation viscosity for low-shear zones.
///
/// - **Taskin (2012) strain-based hemolysis**: cumulative strain history alternative to
///   the instantaneous-shear Giersiepen (1990) power-law. Relevant for multi-stage
///   venturi designs with repeated high-shear transients.
///
/// - **Fahraeus-Lindqvist (Pries 1992)**: apparent viscosity reduction in D_h < 300 um.
///
/// - **Amini (2014) confinement-dependent lift**: refines focusing for a/D_h > 0.1.
///
/// - **Plasma skimming (Pries 1989)**: hematocrit partitioning at asymmetric
///   bifurcations.
///
/// - **Durst (2005) entrance correction**: developing-flow profile for L/D_h < 20.
///
/// - **Bayat-Rezai (2017) millifluidic Dean**: validated Dean correlation for
///   D_h > 500 um, used in GA serpentine evaluation.
pub fn compute_blueprint_report_metrics(
    candidate: &BlueprintCandidate,
) -> Result<SdtMetrics, OptimError> {
    let topology = candidate.topology_spec()?;
    let solve = solve_blueprint_candidate(candidate)?;
    let residence = compute_typed_residence_metrics(candidate, &solve);
    let separation = compute_blueprint_separation_metrics(candidate)?;
    let venturi = compute_blueprint_venturi_metrics(candidate, &solve, &separation)?;
    let safety = compute_typed_blueprint_safety_metrics(candidate, &solve);
    let max_main_channel_shear_pa = safety.max_main_channel_shear.into_base();
    let max_venturi_shear_pa = safety.max_venturi_shear.into_base();
    let pressure_drop_pa = safety.pressure_drop.into_base();
    let mean_residence_time_s = residence.treatment_residence_time.into_base();

    let mut metrics = SdtMetrics::default();
    let flow_rate_m3_s = candidate.operating_point.flow_rate_m3_s.max(1.0e-18);
    let flow_rate_ml_min = flow_rate_m3_s * 6.0e7;
    let flow_rate = VolumetricFlowRate::from_base(flow_rate_m3_s);
    let total_volume: Volume = Time::from_base(solve.mean_residence_time_s) * flow_rate;
    let total_volume_m3 = total_volume.into_base();
    let total_path_length = Length::from_base(
        solve
            .channel_samples
            .iter()
            .map(|sample| sample.length_m)
            .sum::<f64>(),
    );

    let active_venturi_throat_count = topology
        .venturi_placements
        .iter()
        .map(|placement| usize::from(placement.serial_throat_count.max(1)))
        .sum::<usize>();
    let serial_venturi_stages_per_path =
        topology.venturi_placements.first().map_or(0, |placement| {
            usize::from(placement.serial_throat_count.max(1))
        });
    let strongest_venturi = venturi
        .placements
        .iter()
        .min_by(|a, b| a.cavitation_number.total_cmp(&b.cavitation_number));
    let cavitation_number =
        strongest_venturi.map_or(f64::INFINITY, |placement| placement.cavitation_number);
    let cavitation_potential = if cavitation_number.is_finite() {
        (1.0 - cavitation_number).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let cavitation_strength = if cavitation_number.is_finite() {
        cavitation_strength_from_sigma(cavitation_number)
    } else {
        0.0
    };

    let local_hematocrit_venturi =
        candidate.operating_point.feed_hematocrit * venturi.rbc_exposure_fraction.clamp(0.0, 1.0);
    let hematocrit_factor = if candidate.operating_point.feed_hematocrit > 1.0e-12 {
        (local_hematocrit_venturi / candidate.operating_point.feed_hematocrit).clamp(0.0, 1.0)
    } else {
        0.0
    };

    let n_channels = solve.channel_samples.len();
    let mut main_shears = Vec::with_capacity(n_channels);
    let mut main_transits = Vec::with_capacity(n_channels);
    let mut max_venturi_transit_time_s = 0.0_f64;
    let mut max_venturi_shear_rate_inv_s = 0.0_f64;
    let mut dead_volume_m3 = 0.0_f64;
    let mut bulk_hi = 0.0_f64;
    let mut corrected_hi = 0.0_f64;
    let mut treatment_hi = 0.0_f64;
    let mut bypass_hi_values = Vec::with_capacity(n_channels / 2);
    let mut pai_accumulator = 0.0_f64;
    let mut per_channel_hemolysis = Vec::with_capacity(solve.channel_samples.len());

    for sample in &solve.channel_samples {
        let area_m2 = sample.cross_section.area().max(1.0e-18);
        let velocity_m_s = sample.flow_m3_s.abs() / area_m2;
        let shear_rate_inv_s = sample.cross_section.wall_shear_rate(velocity_m_s);
        let shear = DynamicViscosity::from_base(BLOOD_VISCOSITY_PA_S)
            * ReciprocalTime::from_base(shear_rate_inv_s);
        let shear_pa = shear.into_base();
        let transit_time =
            Length::from_base(sample.length_m) / Velocity::from_base(velocity_m_s.max(1.0e-18));
        let transit_time_s = transit_time.into_base();
        let flow_fraction = (sample.flow_m3_s.abs() / flow_rate_m3_s).clamp(0.0, 1.0);
        let local_cavitation = venturi
            .placements
            .iter()
            .find(|placement| {
                sample.id == placement.target_channel_id
                    || sample.id.starts_with(&placement.target_channel_id)
            })
            .map_or(0.0, |placement| {
                (1.0 - placement.cavitation_number).clamp(0.0, 1.0)
            });
        let base_hi = giersiepen_hi(shear_pa, transit_time_s) * flow_fraction;
        // Hellums 1994: PAI = 1.8e-8 × τ^1.325 × t^0.462, flow-weighted.
        let channel_pai = if shear_pa > 0.0 && transit_time_s > 0.0 {
            1.8e-8 * shear_pa.powf(1.325) * transit_time_s.powf(0.462) * flow_fraction
        } else {
            0.0
        };
        pai_accumulator += channel_pai;
        let corrected_channel_hi = if sample.is_venturi_channel {
            cavitation_amplified_hi(base_hi * hematocrit_factor, local_cavitation)
        } else {
            base_hi
        };

        bulk_hi += base_hi;
        corrected_hi += corrected_channel_hi;
        if sample.is_treatment_channel {
            treatment_hi += corrected_channel_hi;
        } else {
            bypass_hi_values.push(corrected_channel_hi);
        }
        if sample.is_venturi_channel {
            max_venturi_transit_time_s = max_venturi_transit_time_s.max(transit_time_s);
            max_venturi_shear_rate_inv_s = max_venturi_shear_rate_inv_s.max(shear_rate_inv_s);
        } else {
            main_shears.push(shear_pa);
            main_transits.push(transit_time_s);
        }
        if shear_rate_inv_s < DEAD_VOLUME_SHEAR_THRESHOLD_INV_S {
            dead_volume_m3 += sample.length_m * area_m2;
        }

        per_channel_hemolysis.push(TypedChannelHemolysis {
            channel_id: sample.id,
            is_venturi_throat: sample.is_venturi_channel,
            hi_contribution: corrected_channel_hi,
            wall_shear: shear,
            transit_time,
            flow_fraction,
        });
    }
    per_channel_hemolysis.sort_by(|a, b| b.hi_contribution.total_cmp(&a.hi_contribution));

    let wall_shear_p95_pa = percentile(&main_shears, 0.95);
    let wall_shear_p99_pa = percentile(&main_shears, 0.99);
    let wall_shear_mean_pa = mean(&main_shears);
    let wall_shear_cv = coefficient_of_variation(&main_shears, wall_shear_mean_pa);

    let min_main_shear_rate_inv_s = main_shears
        .iter()
        .copied()
        .map(|shear_pa| shear_pa / BLOOD_VISCOSITY_PA_S.max(1.0e-18))
        .filter(|shear_rate| *shear_rate > 0.0)
        .fold(f64::INFINITY, f64::min);
    let max_channel_transit_s = main_transits
        .iter()
        .copied()
        .fold(max_venturi_transit_time_s, f64::max);

    let max_expansion_ratio = topology
        .venturi_placements
        .iter()
        .map(|placement| {
            (placement.throat_geometry.inlet_width_m
                / placement.throat_geometry.throat_width_m.max(1.0e-18))
            .max(1.0)
        })
        .fold(1.0_f64, f64::max);

    let channel_resonance_score = solve
        .channel_samples
        .iter()
        .filter(|sample| sample.is_treatment_channel)
        .map(|sample| resonance_match(sample.cross_section.hydraulic_diameter()))
        .fold(0.0_f64, f64::max);
    let throat_temperature_rise_k = if max_venturi_shear_rate_inv_s > 0.0 {
        max_venturi_shear_pa * max_venturi_shear_rate_inv_s * max_venturi_transit_time_s
            / (BLOOD_DENSITY_KG_M3 * C_P_BLOOD_J_KG_K)
    } else {
        0.0
    };

    let constriction_score = strongest_venturi
        .and_then(|placement| {
            topology
                .venturi_placements
                .iter()
                .find(|spec| spec.placement_id == placement.placement_id)
                .and_then(|spec| {
                    solve
                        .channel_samples
                        .iter()
                        .find(|sample| {
                            sample.id == spec.target_channel_id
                                || sample.id.starts_with(&spec.target_channel_id)
                        })
                        .map(|sample| {
                            let inlet_area = (spec.throat_geometry.inlet_width_m
                                * spec.throat_geometry.throat_height_m)
                                .max(1.0e-18);
                            let upstream_velocity = sample.flow_m3_s.abs() / inlet_area;
                            let ratio = placement.effective_throat_velocity_m_s
                                / upstream_velocity.max(1.0e-18);
                            ratio.max(1.0).ln() / VENTURI_VEL_RATIO_REF.ln()
                        })
                })
        })
        .map_or(0.0, |score| score.clamp(0.0, 1.0));
    let cavitation_intensity =
        cavitation_strength * (0.5 + 0.5 * constriction_score.clamp(0.0, 1.0));

    let treatment_fraction = topology.therapy_channel_fraction().clamp(0.0, 1.0);
    let serial_dose_fraction = if serial_venturi_stages_per_path == 0 {
        treatment_fraction
    } else {
        solve.venturi_flow_fraction
            * (1.0 - (1.0 - cavitation_strength).powi(serial_venturi_stages_per_path as i32))
    }
    .clamp(0.0, 1.0);
    let cancer_targeted_cavitation =
        separation.cancer_center_fraction.clamp(0.0, 1.0) * cavitation_intensity;
    let wbc_targeted_cavitation =
        separation.wbc_center_fraction.clamp(0.0, 1.0) * cavitation_intensity;
    let rbc_venturi_protection = separation.rbc_peripheral_fraction.clamp(0.0, 1.0)
        * (1.0 - cavitation_intensity * venturi.rbc_exposure_fraction.clamp(0.0, 1.0));
    let healthy_cell_protection_index =
        compute_healthy_cell_protection_index(wbc_targeted_cavitation, rbc_venturi_protection);

    let absolute_inlet_pressure = candidate.operating_point.absolute_inlet_pressure_pa();
    let collapse_ratio = (absolute_inlet_pressure / BLOOD_VAPOR_PRESSURE_PA.max(1.0))
        .powf((BUBBLE_POLYTROPIC_K - 1.0) / BUBBLE_POLYTROPIC_K);
    let collapse_ref = (SONO_REF_P_ABS_PA / BLOOD_VAPOR_PRESSURE_PA.max(1.0))
        .powf((BUBBLE_POLYTROPIC_K - 1.0) / BUBBLE_POLYTROPIC_K);
    let sonoluminescence_proxy =
        (cavitation_strength * (collapse_ratio / collapse_ref)).clamp(0.0, 1.0);

    let cancer_dose_fraction = separation.cancer_center_fraction.clamp(0.0, 1.0)
        * if active_venturi_throat_count > 0 {
            serial_dose_fraction
        } else {
            treatment_fraction
        };
    let rbc_load = venturi.rbc_exposure_fraction.clamp(0.0, 1.0)
        * local_hematocrit_venturi
        * cavitation_intensity;
    let cancer_rbc_cavitation_bias_index = if cancer_targeted_cavitation + rbc_load > 1.0e-18 {
        cancer_targeted_cavitation / (cancer_targeted_cavitation + rbc_load)
    } else {
        0.0
    };
    let oncology_selectivity_index = cancer_targeted_cavitation
        * (1.0 - venturi.rbc_exposure_fraction.clamp(0.0, 1.0))
        * (0.5 + 0.5 * separation.cancer_center_fraction.clamp(0.0, 1.0));
    let selective_cavitation_delivery_index =
        oncology_selectivity_index * cancer_rbc_cavitation_bias_index * 1.0;
    let lysis_risk_index = corrected_hi
        * (1.0 + 5.0 * venturi.rbc_exposure_fraction.clamp(0.0, 1.0) * local_hematocrit_venturi);
    let therapeutic_window_score =
        (cancer_targeted_cavitation / (1.0e-6 + lysis_risk_index) / THERAPEUTIC_WINDOW_REF)
            .clamp(0.0, 1.0);

    let pediatric_blood_volume_ml =
        PEDIATRIC_REFERENCE_WEIGHT_KG * PEDIATRIC_BLOOD_VOLUME_ML_PER_KG;
    let projected_passes_15min_pediatric_3kg =
        (flow_rate_ml_min * MILESTONE_TREATMENT_DURATION_MIN) / pediatric_blood_volume_ml;
    let projected_hemolysis_15min_pediatric_3kg =
        cumulative_pass_damage(corrected_hi, projected_passes_15min_pediatric_3kg);
    let projected_passes_15min_adult =
        (flow_rate_ml_min * MILESTONE_TREATMENT_DURATION_MIN) / PATIENT_BLOOD_VOLUME_ML;
    let projected_hemolysis_15min_adult =
        cumulative_pass_damage(corrected_hi, projected_passes_15min_adult);

    let optical_path_length = Length::from_base(
        solve
            .channel_samples
            .iter()
            .filter(|sample| sample.is_treatment_channel)
            .map(|sample| sample.cross_section.dims().1)
            .reduce(f64::max)
            .unwrap_or(0.0),
    );
    let optical_path_length_405_m = optical_path_length.into_base();
    let blue_light_delivery_index_405nm = blue_light_delivery_index(
        optical_path_length_405_m,
        candidate.operating_point.feed_hematocrit,
    )
    .map_err(|source| OptimError::PhysicsError {
        id: candidate.id.clone(),
        reason: format!("405-nm optical transport: {source}"),
    })?;

    metrics.cavitation_number = cavitation_number;
    metrics.cavitation_potential = cavitation_potential;
    metrics.throat_exceeds_fda = max_venturi_shear_pa > FDA_MAX_WALL_SHEAR_PA;
    metrics.fda_main_compliant = max_main_channel_shear_pa <= FDA_MAX_WALL_SHEAR_PA;
    metrics.bulk_hemolysis_index_per_pass = bulk_hi;
    metrics.hemolysis_index_per_pass = corrected_hi;
    metrics.flow_uniformity = solve.flow_uniformity;
    metrics.well_coverage_fraction = treatment_fraction;
    // Use TREATMENT-ZONE residence time, not total device residence.
    // The treatment zone is the subset of channels where therapy (cavitation,
    // ultrasound) is applied.  The total device residence includes all bypass
    // arms and is always longer.
    metrics.pressure_feasible = safety.pressure_feasible;
    metrics.cell_separation_efficiency = separation.separation_efficiency;
    metrics.cancer_center_fraction = separation.cancer_center_fraction;
    metrics.rbc_peripheral_fraction = separation.rbc_peripheral_fraction;
    metrics.three_pop_sep_efficiency = separation.separation_efficiency;
    metrics.wbc_center_fraction = separation.wbc_center_fraction;
    metrics.rbc_peripheral_fraction_three_pop = separation.rbc_peripheral_fraction;
    metrics.wbc_equilibrium_pos = 1.0 - separation.wbc_center_fraction.clamp(0.0, 1.0);
    metrics.cancer_equilibrium_pos = 1.0 - separation.cancer_center_fraction.clamp(0.0, 1.0);
    metrics.rbc_equilibrium_pos = separation.rbc_peripheral_fraction.clamp(0.0, 1.0);
    metrics.wbc_recovery = separation.wbc_center_fraction;
    metrics.rbc_pass_fraction = (1.0 - separation.rbc_peripheral_fraction).clamp(0.0, 1.0);
    metrics.wbc_purity = if metrics.wbc_recovery + metrics.rbc_pass_fraction > 1.0e-18 {
        metrics.wbc_recovery / (metrics.wbc_recovery + metrics.rbc_pass_fraction)
    } else {
        0.0
    };
    metrics.plate_fits =
        topology.box_dims_mm.0 <= PLATE_WIDTH_MM && topology.box_dims_mm.1 <= PLATE_HEIGHT_MM;
    let overlap = candidate.blueprint.channel_overlap_analysis();
    metrics.channel_overlap_fraction = overlap.max_overlap_fraction;
    metrics.overlap_width_ratio = overlap.width_ratio_at_worst;
    metrics.n_outlet_ports = candidate
        .blueprint
        .nodes
        .iter()
        .filter(|node| matches!(node.kind, cfd_schematics::domain::model::NodeKind::Outlet))
        .count();
    metrics.low_flow_stasis_risk = inverse_linear_risk(
        flow_rate_ml_min,
        CLOTTING_BFR_HIGH_RISK_ML_MIN,
        CLOTTING_BFR_LOW_RISK_ML_MIN,
    );
    metrics.min_shear_stasis_risk = direct_linear_risk(
        min_main_shear_rate_inv_s,
        CLOTTING_SHEAR_HIGH_RISK_INV_S,
        CLOTTING_SHEAR_LOW_RISK_INV_S,
    );
    metrics.max_residence_stasis_risk = direct_linear_risk(
        max_channel_transit_s,
        CLOTTING_RESIDENCE_LOW_RISK_S,
        CLOTTING_RESIDENCE_HIGH_RISK_S,
    );
    metrics.expansion_stasis_risk = log_risk(
        max_expansion_ratio,
        EXPANSION_RATIO_LOW_RISK,
        VENTURI_EXPANSION_RATIO_HIGH_RISK,
    );
    metrics.dead_volume_stasis_risk = if total_volume_m3 > 1.0e-18 {
        (dead_volume_m3 / total_volume_m3).clamp(0.0, 1.0)
    } else {
        0.0
    };
    metrics.clotting_risk_index = (metrics.low_flow_stasis_risk
        + metrics.min_shear_stasis_risk
        + metrics.max_residence_stasis_risk
        + metrics.expansion_stasis_risk
        + metrics.dead_volume_stasis_risk)
        / 5.0;
    metrics.clotting_flow_compliant = flow_rate_ml_min >= CLOTTING_BFR_CAUTION_ML_MIN;
    let strict_flow_risk = inverse_linear_risk(
        flow_rate_ml_min,
        CLOTTING_BFR_HIGH_RISK_ML_MIN,
        CLOTTING_BFR_STRICT_10MLS_ML_MIN,
    );
    metrics.clotting_risk_index_10ml_s = (strict_flow_risk
        + metrics.min_shear_stasis_risk
        + metrics.max_residence_stasis_risk
        + metrics.expansion_stasis_risk
        + metrics.dead_volume_stasis_risk)
        / 5.0;
    metrics.clotting_flow_compliant_10ml_s = flow_rate_ml_min >= CLOTTING_BFR_STRICT_10MLS_ML_MIN;
    metrics.pediatric_flow_excess_risk = direct_linear_risk(
        flow_rate_ml_min,
        PEDIATRIC_FLOW_CAUTION_ML_MIN,
        PEDIATRIC_FLOW_EXCESSIVE_ML_MIN,
    );
    metrics.pediatric_flow_compliant = flow_rate_ml_min <= PEDIATRIC_FLOW_CAUTION_ML_MIN;
    metrics.venturi_treatment_enabled = active_venturi_throat_count > 0;
    metrics.treatment_zone_mode = match topology.treatment_mode {
        TreatmentActuationMode::UltrasoundOnly => "UltrasoundOnly",
        TreatmentActuationMode::VenturiCavitation => {
            if topology.has_serpentine() {
                "UltrasoundPlusVenturiThroats"
            } else {
                "VenturiThroats"
            }
        }
    }
    .to_string();
    metrics.active_venturi_throat_count = active_venturi_throat_count;
    metrics.serial_venturi_stages_per_path = serial_venturi_stages_per_path;
    metrics.venturi_flow_fraction = solve.venturi_flow_fraction;
    metrics.rbc_venturi_exposure_fraction = venturi.rbc_exposure_fraction;
    metrics.local_hematocrit_venturi = local_hematocrit_venturi;
    metrics.cancer_dose_fraction = cancer_dose_fraction;
    metrics.cavitation_intensity = cavitation_intensity;
    metrics.cancer_targeted_cavitation = cancer_targeted_cavitation;
    metrics.wbc_targeted_cavitation = wbc_targeted_cavitation;
    metrics.rbc_venturi_protection = rbc_venturi_protection;
    metrics.healthy_cell_protection_index = healthy_cell_protection_index;
    metrics.sonoluminescence_proxy = sonoluminescence_proxy;
    metrics.fda_overall_compliant = metrics.fda_main_compliant
        && (max_venturi_shear_pa <= FDA_MAX_WALL_SHEAR_PA
            || (max_venturi_transit_time_s <= FDA_TRANSIENT_TIME_S
                && max_venturi_shear_pa <= FDA_TRANSIENT_SHEAR_PA));
    metrics.lysis_risk_index = lysis_risk_index;
    metrics.therapeutic_window_score = therapeutic_window_score;
    metrics.oncology_selectivity_index = oncology_selectivity_index;
    metrics.cancer_rbc_cavitation_bias_index = cancer_rbc_cavitation_bias_index;
    metrics.selective_cavitation_delivery_index = selective_cavitation_delivery_index;
    metrics.rbc_lysis_rate_pct_per_h =
        corrected_hi * flow_rate_ml_min * 60.0 / PATIENT_BLOOD_VOLUME_ML * 100.0;
    metrics.projected_passes_15min_pediatric_3kg = projected_passes_15min_pediatric_3kg;
    metrics.projected_hemolysis_15min_pediatric_3kg = projected_hemolysis_15min_pediatric_3kg;
    metrics.projected_hemolysis_15min_adult = projected_hemolysis_15min_adult;
    metrics.cancer_therapy_zone_fraction =
        separation.cancer_center_fraction.clamp(0.0, 1.0) * solve.venturi_flow_fraction;
    metrics.blue_light_delivery_index_405nm = blue_light_delivery_index_405nm;
    metrics.therapy_channel_fraction = treatment_fraction;
    metrics.wall_shear_cv = wall_shear_cv;
    metrics.fda_shear_percentile_compliant =
        wall_shear_p95_pa <= FDA_MAX_WALL_SHEAR_PA && wall_shear_p99_pa <= FDA_TRANSIENT_SHEAR_PA;
    let hydraulic_power = mechanical_power(pressure_drop_pa, flow_rate_m3_s);
    let hydraulic_power_w = hydraulic_power.into_base();
    metrics.acoustic_capture_efficiency =
        (cavitation_potential * (max_venturi_shear_pa / P_ATM_PA)).clamp(0.0, 1.0);
    let specific_cavitation_energy_j_ml = if flow_rate_ml_min > 1.0e-18 {
        metrics.acoustic_capture_efficiency * hydraulic_power_w * mean_residence_time_s
            / (flow_rate_ml_min / 60_000.0)
            * 1.0e3
    } else {
        0.0
    };
    metrics.hemolysis_index_per_pass_cavitation_amplified = corrected_hi;
    metrics.treatment_channel_hi = treatment_hi;
    metrics.bypass_channel_hi_mean = mean(&bypass_hi_values);
    metrics.bypass_channel_hi_max = bypass_hi_values.into_iter().fold(0.0, f64::max);
    metrics.per_channel_hemolysis = per_channel_hemolysis
        .into_iter()
        .map(TypedChannelHemolysis::into_serialized)
        .collect();
    metrics.acoustic_resonance_factor = channel_resonance_score;
    metrics.channel_resonance_score = channel_resonance_score;
    metrics.serial_cavitation_dose_fraction = serial_dose_fraction;
    metrics.fda_thermal_compliant = throat_temperature_rise_k <= FDA_THROAT_TEMP_RISE_LIMIT_K;

    // ── Previously unset metrics ─────────────────────────────────────────────
    // Platelet activation index (Hellums 1994 power-law model).
    metrics.platelet_activation_index = pai_accumulator;

    // Diffuser pressure recovery: maximum across all venturi placements.
    let diffuser_recovery_pa = venturi
        .placements
        .iter()
        .map(|p| p.diffuser_recovery_pa)
        .fold(0.0_f64, f64::max);
    metrics.venturi_total_loss_coefficient = venturi
        .placements
        .iter()
        .map(|p| p.total_loss_coefficient)
        .fold(0.0_f64, f64::max);

    // Outlet tail length and remerge proximity: applicable to all PST topologies.
    // Shorter outlet tails mean treatment-stream remerge happens closer to chip exit,
    // reducing post-therapy dilution and secondary separation of treated cells.
    let outlet_tail_mm = topology.outlet_tail_length_m * 1.0e3;
    metrics.cif_outlet_tail_length_mm = outlet_tail_mm;
    metrics.cif_remerge_proximity_score = if outlet_tail_mm > 0.0 {
        // Exponential decay: 1 mm → 0.95, 5 mm → 0.37, 10 mm → 0.14.
        (-outlet_tail_mm / 5.0_f64).exp().clamp(0.0, 1.0)
    } else {
        0.0
    };

    // Topology-specific split-stage flow fractions via rectangular laminar conductance.
    let (stage_center_fracs, model_frac) = split_stage_flow_fractions(&topology.split_stages);
    // Assign to both CCT and CIF fields since these generalise across topologies.
    metrics.cct_model_venturi_flow_fraction = model_frac.clamp(0.0, 1.0);
    metrics.cif_model_venturi_flow_fraction = model_frac.clamp(0.0, 1.0);
    metrics.cct_solved_venturi_flow_fraction = solve.venturi_flow_fraction;
    metrics.cif_solved_venturi_flow_fraction = solve.venturi_flow_fraction;
    metrics.cct_stage_center_qfrac_mean = if stage_center_fracs.is_empty() {
        0.0
    } else {
        stage_center_fracs.iter().sum::<f64>() / stage_center_fracs.len() as f64
    };
    // CIF pre-trifurcation mean: all stages except the last (which is the terminal).
    metrics.cif_pretri_qfrac_mean = if stage_center_fracs.len() > 1 {
        let pretri = &stage_center_fracs[..stage_center_fracs.len() - 1];
        pretri.iter().sum::<f64>() / pretri.len() as f64
    } else {
        metrics.cct_stage_center_qfrac_mean
    };
    // Terminal trifurcation and bifurcation fractions.
    metrics.cif_terminal_tri_qfrac = stage_center_fracs.last().copied().unwrap_or(0.0);
    metrics.cif_terminal_bi_qfrac = stage_center_fracs.first().copied().unwrap_or(0.0);

    // Update selective_cavitation_delivery_index to include remerge bonus
    // now that cif_remerge_proximity_score is properly computed.
    let remerge_bonus = 0.5 + 0.5 * metrics.cif_remerge_proximity_score;
    metrics.selective_cavitation_delivery_index =
        oncology_selectivity_index * cancer_rbc_cavitation_bias_index * remerge_bonus;

    // ── SDT Acoustic Physics (narrative-only augmentation) ───────────────────
    //
    // Wire the validated cfd-1d acoustic models into the report pipeline.
    let sdt_acoustic = compute_sdt_acoustic_metrics(
        cavitation_intensity,
        max_venturi_transit_time_s,
        pressure_drop_pa,
    );
    metrics.ctc_contrast_factor = sdt_acoustic.ctc_contrast_factor;
    metrics.rbc_contrast_factor = sdt_acoustic.rbc_contrast_factor;

    // Sonosensitizer activation modulates cancer dose (Rosenthal 2004):
    // short throat transits (< 1 ms) yield incomplete activation.
    // Zero cavitation intensity produces zero activation and therefore zero dose.
    let cancer_dose_fraction = cancer_dose_fraction * sdt_acoustic.sensitizer_activation_efficiency;
    metrics.cancer_dose_fraction = cancer_dose_fraction;

    // Rayleigh-Plesset collapse amplification (Rayleigh 1917):
    // bubble collapse micro-jets increase hemolysis beyond the base shear model.
    if cavitation_potential > 0.0 {
        metrics.hemolysis_index_per_pass_cavitation_amplified *=
            sdt_acoustic.rayleigh_plesset_amplification;
    }

    let physical_metrics = TypedReportPhysicalMetrics {
        throat_shear_rate: ReciprocalTime::from_base(max_venturi_shear_rate_inv_s),
        throat_shear: safety.max_venturi_shear,
        max_main_channel_shear: safety.max_main_channel_shear,
        mean_residence_time: residence.treatment_residence_time,
        total_pressure_drop: safety.pressure_drop,
        total_path_length,
        total_ecv: total_volume,
        flow_rate,
        main_channel_shear_rate: ReciprocalTime::from_base(
            max_main_channel_shear_pa / BLOOD_VISCOSITY_PA_S.max(1.0e-18),
        ),
        throat_transit_time: Time::from_base(max_venturi_transit_time_s),
        optical_path_length,
        safety_margin: Pressure::from_base(FDA_MAX_WALL_SHEAR_PA - max_main_channel_shear_pa),
        wall_shear_p95: Pressure::from_base(wall_shear_p95_pa),
        wall_shear_p99: Pressure::from_base(wall_shear_p99_pa),
        wall_shear_mean: Pressure::from_base(wall_shear_mean_pa),
        diffuser_recovery: Pressure::from_base(diffuser_recovery_pa),
        mechanical_power: hydraulic_power,
        treatment_zone_dwell_time: residence.treatment_residence_time,
        specific_cavitation_energy: EnergyPerVolume::from_unit::<JoulePerMilliliter>(
            specific_cavitation_energy_j_ml,
        ),
        acoustic_energy_density: sdt_acoustic.acoustic_energy_density,
        throat_temperature_rise: TemperatureDifference::from_unit::<Kelvin>(
            throat_temperature_rise_k,
        ),
    };
    physical_metrics.write_to(&mut metrics);

    Ok(metrics)
}

// ── SDT Acoustic Metrics (Gor'kov 1962, Rosenthal 2004, Rayleigh 1917) ──────

/// SDT acoustic metrics computed from the validated cfd-1d physics models.
///
/// These AUGMENT the existing cavitation-based metrics.
/// Sensitizer activation and RP amplification are wired into the report pipeline;
/// acoustic energy density and contrast factors are available for narrative reporting.
#[derive(Debug, Clone, Copy)]
pub struct SdtAcousticMetrics {
    /// Sonosensitizer activation fraction (Rosenthal 2004 first-order kinetics).
    /// η_act = 1 − exp(−k_act · I_cav · t_transit), range [0, 1].
    pub sensitizer_activation_efficiency: f64,

    /// Rayleigh-Plesset collapse hemolysis amplification factor (≥ 1.0).
    /// A_collapse = 1 + α · (R_max/R_0)² · (p_inf/p_ref).
    pub rayleigh_plesset_amplification: f64,

    /// Acoustic energy density [J/m³] at 100 kPa pressure amplitude in blood.
    /// E_ac = p₀²/(4ρc²).
    pub acoustic_energy_density: EnergyPerVolume,

    /// Acoustic contrast factor Φ for CTCs in plasma (Gor'kov 1962).
    /// Positive → migrates toward pressure nodes.
    pub ctc_contrast_factor: f64,

    /// Acoustic contrast factor Φ for RBCs in plasma (Gor'kov 1962).
    /// Positive → migrates toward pressure nodes.
    pub rbc_contrast_factor: f64,
}

/// Compute SDT acoustic metrics from the validated cfd-1d physics models.
///
/// This function wires together:
/// - **Sonosensitizer activation** (Rosenthal 2004): first-order kinetics
///   η_act = 1 − exp(−k_act · I_cav · t_transit)
/// - **Rayleigh-Plesset collapse** (Rayleigh 1917): hemolysis amplification
///   from cavitation bubble collapse, A = 1 + α·(R_max/R_0)²·(p/p_ref)
/// - **Acoustic energy density** (Gor'kov 1962): E_ac = p₀²/(4ρc²)
/// - **Acoustic contrast factors**: differential radiation force on CTCs vs RBCs
///
/// # Arguments
///
/// * `cavitation_intensity` — dimensionless cavitation intensity [0, 1]
/// * `throat_transit_time_s` — residence time in the cavitation zone \[s]
/// * `pressure_drop_pa` — total pressure drop across the chip \[Pa]
pub fn compute_sdt_acoustic_metrics(
    cavitation_intensity: f64,
    throat_transit_time_s: f64,
    pressure_drop_pa: f64,
) -> SdtAcousticMetrics {
    // Sonosensitizer activation efficiency (Rosenthal 2004)
    // Guard inputs: transit time must be non-negative, intensity in [0, 1].
    let safe_transit = throat_transit_time_s.max(0.0);
    let safe_intensity = cavitation_intensity.clamp(0.0, 1.0);
    let sensitizer_activation = sonosensitizer_activation_efficiency(
        SENSITIZER_K_ACT_HEMATOPORPHYRIN,
        safe_intensity,
        safe_transit,
    );

    // Rayleigh-Plesset collapse amplification (Rayleigh 1917)
    // Clamp absolute pressure to realistic millifluidic range [50 kPa, 500 kPa]
    // to prevent over-amplification at extreme pressure drops.
    let r_max = 10.0e-6; // 10 µm typical cavitation bubble
    let r_0 = 1.0e-6; // 1 µm equilibrium nucleus
    let p_inf = (pressure_drop_pa + 101_325.0).clamp(50_000.0, 500_000.0);
    let rp_amplification = cavitation_hemolysis_amplification(r_max, r_0, p_inf);

    // Acoustic energy density (Gor'kov 1962)
    // Use actual chip pressure for consistency with RP dynamics.
    let p_acoustic = (pressure_drop_pa + 101_325.0).max(50_000.0);
    let e_acoustic = acoustic_energy_density(
        p_acoustic, // consistent with RP pressure reference
        1060.0,     // blood density [kg/m³]
        1540.0,     // speed of sound in blood [m/s]
    );

    // Acoustic contrast factors for differential radiation force
    let ctc_contrast = acoustic_contrast_factor(RHO_CTC, RHO_PLASMA, KAPPA_CTC, KAPPA_PLASMA);
    let rbc_contrast = acoustic_contrast_factor(RHO_RBC, RHO_PLASMA, KAPPA_RBC, KAPPA_PLASMA);

    SdtAcousticMetrics {
        sensitizer_activation_efficiency: sensitizer_activation,
        rayleigh_plesset_amplification: rp_amplification,
        acoustic_energy_density: EnergyPerVolume::from_unit::<JoulePerCubicMeter>(e_acoustic),
        ctc_contrast_factor: ctc_contrast,
        rbc_contrast_factor: rbc_contrast,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::fixtures::{
        canonical_option1_candidate, operating_point, stage0_venturi_candidate,
    };
    use crate::metrics::{
        compute_blueprint_separation_metrics, compute_blueprint_venturi_metrics,
        solve_blueprint_candidate,
    };
    use cfd_schematics::topology::VenturiPlacementMode;

    #[test]
    fn blue_light_delivery_is_unity_for_zero_path() {
        let transmission = blue_light_delivery_index(0.0, 0.42).expect("valid optical input");

        assert_eq!(transmission, 1.0);
    }

    #[test]
    fn blue_light_delivery_preserves_nonnegative_hematocrit_policy() {
        let transmission = blue_light_delivery_index(0.002, -0.2).expect("valid optical input");

        assert_eq!(transmission, 1.0);
    }

    #[test]
    fn blue_light_delivery_matches_beer_lambert_oracle() {
        let optical_path_length_m = 0.002;
        let feed_hematocrit = 0.42;
        let transmission = blue_light_delivery_index(optical_path_length_m, feed_hematocrit)
            .expect("valid optical input");
        let expected =
            (-BLOOD_ATTENUATION_405NM_INV_M * optical_path_length_m * feed_hematocrit).exp();
        let scaled_epsilon = 32.0 * f64::EPSILON;
        let tolerance = scaled_epsilon / (1.0 - scaled_epsilon) * expected.abs().max(1.0);

        assert!(
            (transmission - expected).abs() <= tolerance,
            "Hyperion transmission {transmission} differs from the f64 Beer-Lambert oracle {expected} by more than {tolerance}"
        );
    }

    #[test]
    fn blue_light_delivery_rejects_negative_path() {
        let error = blue_light_delivery_index(-0.001, 0.42)
            .expect_err("negative optical path must fail validation");

        assert!(matches!(
            error,
            hyperion::TransportError::InvalidValue {
                field: hyperion::ValueKind::PathLength,
                ..
            }
        ));
    }

    #[test]
    fn test_sdt_acoustic_metrics_computable() {
        let m = compute_sdt_acoustic_metrics(0.6, 0.002, 50_000.0);
        assert!(
            m.sensitizer_activation_efficiency.is_finite(),
            "sensitizer_activation_efficiency must be finite"
        );
        assert!(
            m.rayleigh_plesset_amplification.is_finite(),
            "rayleigh_plesset_amplification must be finite"
        );
        assert!(
            m.acoustic_energy_density
                .in_unit::<JoulePerCubicMeter>()
                .is_finite(),
            "acoustic_energy_density must be finite"
        );
        assert!(
            m.ctc_contrast_factor.is_finite(),
            "ctc_contrast_factor must be finite"
        );
        assert!(
            m.rbc_contrast_factor.is_finite(),
            "rbc_contrast_factor must be finite"
        );
        // Activation should be in (0, 1) for non-zero inputs
        assert!(
            m.sensitizer_activation_efficiency > 0.0 && m.sensitizer_activation_efficiency < 1.0,
            "activation = {} should be in (0,1)",
            m.sensitizer_activation_efficiency
        );
        // RP amplification must exceed 1.0
        assert!(
            m.rayleigh_plesset_amplification > 1.0,
            "RP amplification = {} should exceed 1.0",
            m.rayleigh_plesset_amplification
        );
    }

    #[test]
    fn test_sdt_acoustic_metrics_zero_cavitation() {
        let m = compute_sdt_acoustic_metrics(0.0, 0.002, 50_000.0);
        assert_eq!(
            m.sensitizer_activation_efficiency, 0.0,
            "zero cavitation intensity should yield zero activation"
        );
        // RP amplification and energy density are independent of cavitation intensity
        assert!(m.rayleigh_plesset_amplification > 1.0);
        assert!(m.acoustic_energy_density.in_unit::<JoulePerCubicMeter>() > 0.0);
    }

    #[test]
    fn test_sdt_acoustic_metrics_positive_energy() {
        // Energy density must always be positive (p₀² / (4ρc²) > 0)
        for pressure_drop in [0.0, 10_000.0, 100_000.0, 500_000.0] {
            let m = compute_sdt_acoustic_metrics(0.5, 0.001, pressure_drop);
            assert!(
                m.acoustic_energy_density.in_unit::<JoulePerCubicMeter>() > 0.0,
                "energy density must be positive at pressure_drop={pressure_drop}"
            );
        }
    }

    #[test]
    fn report_metrics_propagate_max_venturi_total_loss_coefficient() {
        let candidate = stage0_venturi_candidate(
            "report-loss",
            operating_point(2.0e-6, 30_000.0, 0.18),
            VenturiPlacementMode::StraightSegment,
        );
        let solve = solve_blueprint_candidate(&candidate).expect("solve");
        let separation = compute_blueprint_separation_metrics(&candidate).expect("separation");
        let venturi =
            compute_blueprint_venturi_metrics(&candidate, &solve, &separation).expect("venturi");
        let metrics = compute_blueprint_report_metrics(&candidate).expect("report metrics");
        let expected = venturi
            .placements
            .iter()
            .map(|placement| placement.total_loss_coefficient)
            .fold(0.0_f64, f64::max);

        assert!(
            metrics.venturi_total_loss_coefficient.is_finite()
                && metrics.venturi_total_loss_coefficient >= 0.0,
            "report-level total loss coefficient must be finite and non-negative"
        );
        assert!(
            (metrics.venturi_total_loss_coefficient - expected).abs() < 1.0e-12,
            "expected report total loss coefficient {}, got {}",
            expected,
            metrics.venturi_total_loss_coefficient
        );
    }

    #[test]
    fn report_metrics_zero_cavitation_produces_zero_cancer_dose() {
        let candidate = canonical_option1_candidate(
            "zero-cavitation-report",
            operating_point(2.0e-6, 30_000.0, 0.18),
        );
        let metrics = compute_blueprint_report_metrics(&candidate).expect("report metrics");

        assert_eq!(
            metrics.cancer_dose_fraction, 0.0,
            "zero cavitation intensity must not inject cancer dose"
        );
    }

    #[test]
    fn mechanical_power_composes_pressure_and_volumetric_flow() {
        assert_eq!(
            mechanical_power(12.0, 0.5).into_base().to_bits(),
            6.0_f64.to_bits()
        );
    }

    #[test]
    fn typed_report_physical_metrics_preserve_serialized_units() {
        let mut metrics = SdtMetrics::default();
        TypedReportPhysicalMetrics {
            throat_shear_rate: ReciprocalTime::from_base(2_000.0),
            throat_shear: Pressure::from_base(12.0),
            max_main_channel_shear: Pressure::from_base(8.0),
            mean_residence_time: Time::from_base(0.25),
            total_pressure_drop: Pressure::from_base(40.0),
            total_path_length: Length::from_base(0.012),
            total_ecv: Volume::from_base(3.0e-6),
            flow_rate: VolumetricFlowRate::from_base(2.0e-6),
            main_channel_shear_rate: ReciprocalTime::from_base(800.0),
            throat_transit_time: Time::from_base(0.001),
            optical_path_length: Length::from_base(0.002),
            safety_margin: Pressure::from_base(142.0),
            wall_shear_p95: Pressure::from_base(10.0),
            wall_shear_p99: Pressure::from_base(11.0),
            wall_shear_mean: Pressure::from_base(6.0),
            diffuser_recovery: Pressure::from_base(4.0),
            mechanical_power: Power::from_base(80.0),
            treatment_zone_dwell_time: Time::from_base(0.25),
            specific_cavitation_energy: EnergyPerVolume::from_unit::<JoulePerMilliliter>(0.125),
            acoustic_energy_density: EnergyPerVolume::from_unit::<JoulePerCubicMeter>(2.5),
            throat_temperature_rise: TemperatureDifference::from_unit::<Kelvin>(0.02),
        }
        .write_to(&mut metrics);

        assert_eq!(metrics.total_path_length_mm, 12.0);
        assert_eq!(metrics.total_ecv_ml, 3.0);
        assert_eq!(metrics.flow_rate_ml_min, 120.0);
        assert_eq!(metrics.total_pressure_drop_pa, 40.0);
        assert_eq!(metrics.mechanical_power_w, 80.0);
        assert_eq!(metrics.specific_cavitation_energy_j_ml, 0.125);
        assert_eq!(metrics.acoustic_energy_density_j_m3, 2.5);
        assert_eq!(metrics.throat_temperature_rise_k, 0.02);

        let encoded = serde_json::to_string(&metrics).expect("report metrics serialize");
        let decoded: SdtMetrics =
            serde_json::from_str(&encoded).expect("report metrics deserialize");
        assert_eq!(decoded.total_path_length_mm, 12.0);
        assert_eq!(decoded.total_ecv_ml, 3.0);
        assert_eq!(decoded.flow_rate_ml_min, 120.0);
        assert_eq!(decoded.mechanical_power_w, 80.0);
    }

    #[test]
    fn typed_channel_hemolysis_serializes_provider_values() {
        let typed = TypedChannelHemolysis {
            channel_id: "treatment-0",
            is_venturi_throat: true,
            hi_contribution: 0.25,
            wall_shear: Pressure::from_base(18.0),
            transit_time: Time::from_base(0.004),
            flow_fraction: 0.75,
        };

        let serialized = typed.into_serialized();

        assert_eq!(serialized.channel_id, "treatment-0");
        assert!(serialized.is_venturi_throat);
        assert_eq!(serialized.hi_contribution.to_bits(), 0.25_f64.to_bits());
        assert_eq!(serialized.wall_shear_pa.to_bits(), 18.0_f64.to_bits());
        assert_eq!(serialized.transit_time_s.to_bits(), 0.004_f64.to_bits());
        assert_eq!(serialized.flow_fraction.to_bits(), 0.75_f64.to_bits());
    }

    #[test]
    fn residence_and_safety_adapters_preserve_typed_values() {
        let candidate = stage0_venturi_candidate(
            "typed-residence-safety",
            operating_point(2.0e-6, 30_000.0, 0.18),
            VenturiPlacementMode::StraightSegment,
        );
        let solve = solve_blueprint_candidate(&candidate).expect("solve");
        let typed_residence = compute_typed_residence_metrics(&candidate, &solve);
        let serialized_residence = crate::metrics::compute_residence_metrics(&candidate, &solve);
        let typed_safety = compute_typed_blueprint_safety_metrics(&candidate, &solve);
        let serialized_safety =
            crate::metrics::compute_blueprint_safety_metrics(&candidate, &solve);

        let expected_treatment_volume = solve
            .channel_samples
            .iter()
            .filter(|sample| sample.is_treatment_channel)
            .map(|sample| sample.length_m * sample.cross_section.area())
            .sum::<f64>();

        assert_eq!(
            typed_residence.treatment_volume.into_base().to_bits(),
            expected_treatment_volume.to_bits()
        );
        assert_eq!(
            typed_residence.treatment_volume.into_base().to_bits(),
            serialized_residence.treatment_volume_m3.to_bits()
        );
        assert_eq!(
            typed_residence
                .treatment_residence_time
                .into_base()
                .to_bits(),
            serialized_residence.treatment_residence_time_s.to_bits()
        );
        assert_eq!(
            typed_safety.max_main_channel_shear.into_base().to_bits(),
            serialized_safety.max_main_channel_shear_pa.to_bits()
        );
        assert_eq!(
            typed_safety.pressure_drop.into_base().to_bits(),
            serialized_safety.pressure_drop_pa.to_bits()
        );
        assert_eq!(
            typed_safety
                .mean_device_residence_time
                .into_base()
                .to_bits(),
            serialized_safety.mean_device_residence_time_s.to_bits()
        );
    }
}
