use cfd_1d::{
    acoustic_contrast_factor, acoustic_energy_density, cavitation_amplified_hi,
    cavitation_hemolysis_amplification, giersiepen_hi, sonosensitizer_activation_efficiency,
    KAPPA_CTC, KAPPA_PLASMA, KAPPA_RBC, RHO_CTC, RHO_PLASMA, RHO_RBC,
    SENSITIZER_K_ACT_HEMATOPORPHYRIN,
};
use cfd_schematics::topology::TreatmentActuationMode;

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
    PATIENT_BLOOD_VOLUME_ML, PEDIATRIC_BLOOD_VOLUME_ML_PER_KG, PEDIATRIC_REFERENCE_WEIGHT_KG,
    PLATE_HEIGHT_MM, PLATE_WIDTH_MM, P_ATM_PA, SONO_REF_P_ABS_PA, THERAPEUTIC_WINDOW_REF,
    VENTURI_EXPANSION_RATIO_HIGH_RISK, VENTURI_VEL_RATIO_REF,
};
use crate::domain::BlueprintCandidate;
use crate::error::OptimError;
use crate::metrics::{
    compute_blueprint_safety_metrics, compute_blueprint_separation_metrics,
    compute_blueprint_venturi_metrics, compute_residence_metrics, solve_blueprint_candidate,
    ChannelHemolysis, SdtMetrics,
};

/// Compute the full set of report-grade SDT metrics for a single blueprint candidate.
///
/// # Physics models available for future integration
///
/// The following validated physics models can refine the metrics computed here:
///
/// - **Quemada (1978) viscosity**: The rouleaux aggregation viscosity model can be used
///   to refine `local_hematocrit_venturi` in low-shear recirculation zones where RBC
///   aggregation elevates effective viscosity above the Casson baseline.
///
/// - **Taskin (2012) strain-based hemolysis**: This model can provide alternative
///   `hemolysis_index_per_pass` values by tracking cumulative strain history along
///   particle trajectories, rather than relying on the instantaneous-shear Giersiepen
///   (1990) power-law correlation. Particularly relevant for multi-stage venturi designs
///   where cells experience repeated high-shear transients.
///
/// - **Fahraeus-Lindqvist (Pries 1992) correction**: Affects effective resistance in
///   microchannel designs where D_h < 300 um. The apparent viscosity reduction from
///   cell-free layer formation can lower computed pressure drops by 10-30% relative to
///   the constant-viscosity Hagen-Poiseuille model used here.
///
/// - **Amini (2014) confinement-dependent lift**: Refines inertial focusing equilibrium
///   positions for cancer cells with a/D_h > 0.1, improving `cancer_center_fraction`
///   and `wbc_center_fraction` predictions.
///
/// - **Plasma skimming (Pries 1989)**: Hematocrit partitioning at asymmetric
///   bifurcations can improve `rbc_venturi_exposure_fraction` accuracy by accounting
///   for the phase-separation effect at daughter-branch flow splits.
///
/// - **Durst (2005) entrance correction**: For venturi throats with L/D_h < 20, the
///   developing-flow profile increases centreline velocity and wall shear relative to
///   the fully-developed parabolic assumption, affecting `cavitation_number` and
///   `throat_shear_rate_inv_s`.
///
/// - **Bayat-Rezai (2017) millifluidic Dean**: Validated Dean correlation for channels
///   with D_h > 500 um, used in GA serpentine evaluation for secondary flow estimation.
pub fn compute_blueprint_report_metrics(
    candidate: &BlueprintCandidate,
) -> Result<SdtMetrics, OptimError> {
    let topology = candidate.topology_spec()?;
    let solve = solve_blueprint_candidate(candidate)?;
    let residence = compute_residence_metrics(candidate, &solve);
    let separation = compute_blueprint_separation_metrics(candidate)?;
    let venturi = compute_blueprint_venturi_metrics(candidate, &solve, &separation)?;
    let safety = compute_blueprint_safety_metrics(candidate, &solve);

    let mut metrics = SdtMetrics::default();
    let flow_rate_m3_s = candidate.operating_point.flow_rate_m3_s.max(1.0e-18);
    let flow_rate_ml_min = flow_rate_m3_s * 6.0e7;
    let total_volume_m3 = solve.mean_residence_time_s * flow_rate_m3_s;
    let total_path_length_mm = solve
        .channel_samples
        .iter()
        .map(|sample| sample.length_m)
        .sum::<f64>()
        * 1.0e3;

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

    let local_hematocrit_venturi =
        candidate.operating_point.feed_hematocrit * venturi.rbc_exposure_fraction.clamp(0.0, 1.0);
    let hematocrit_factor = if candidate.operating_point.feed_hematocrit > 1.0e-12 {
        (local_hematocrit_venturi / candidate.operating_point.feed_hematocrit).clamp(0.0, 1.0)
    } else {
        0.0
    };

    let mut main_shears = Vec::new();
    let mut main_transits = Vec::new();
    let mut max_venturi_transit_time_s = 0.0_f64;
    let mut max_venturi_shear_rate_inv_s = 0.0_f64;
    let mut dead_volume_m3 = 0.0_f64;
    let mut bulk_hi = 0.0_f64;
    let mut corrected_hi = 0.0_f64;
    let mut treatment_hi = 0.0_f64;
    let mut bypass_hi_values = Vec::new();
    let mut pai_accumulator = 0.0_f64;
    let mut per_channel_hemolysis = Vec::with_capacity(solve.channel_samples.len());

    for sample in &solve.channel_samples {
        let area_m2 = sample.cross_section.area().max(1.0e-18);
        let velocity_m_s = sample.flow_m3_s.abs() / area_m2;
        let shear_rate_inv_s = sample.cross_section.wall_shear_rate(velocity_m_s);
        let shear_pa = BLOOD_VISCOSITY_PA_S * shear_rate_inv_s;
        let transit_time_s = sample.length_m / velocity_m_s.max(1.0e-18);
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
        if !sample.is_venturi_channel {
            main_shears.push(shear_pa);
            main_transits.push(transit_time_s);
        } else {
            max_venturi_transit_time_s = max_venturi_transit_time_s.max(transit_time_s);
            max_venturi_shear_rate_inv_s = max_venturi_shear_rate_inv_s.max(shear_rate_inv_s);
        }
        if shear_rate_inv_s < DEAD_VOLUME_SHEAR_THRESHOLD_INV_S {
            dead_volume_m3 += sample.length_m * area_m2;
        }

        per_channel_hemolysis.push(ChannelHemolysis {
            channel_id: sample.id.to_string(),
            is_venturi_throat: sample.is_venturi_channel,
            hi_contribution: corrected_channel_hi,
            wall_shear_pa: shear_pa,
            transit_time_s,
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
        safety.max_venturi_shear_pa * max_venturi_shear_rate_inv_s * max_venturi_transit_time_s
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
        cavitation_potential * (0.5 + 0.5 * constriction_score.clamp(0.0, 1.0));

    let treatment_fraction = topology.therapy_channel_fraction().clamp(0.0, 1.0);
    let serial_dose_fraction = if serial_venturi_stages_per_path == 0 {
        treatment_fraction
    } else {
        solve.venturi_flow_fraction
            * (1.0 - (1.0 - cavitation_potential).powi(serial_venturi_stages_per_path as i32))
    }
    .clamp(0.0, 1.0);
    let cancer_targeted_cavitation =
        separation.cancer_center_fraction.clamp(0.0, 1.0) * cavitation_intensity;
    let wbc_targeted_cavitation =
        separation.wbc_center_fraction.clamp(0.0, 1.0) * cavitation_intensity;
    let rbc_venturi_protection = separation.rbc_peripheral_fraction.clamp(0.0, 1.0)
        * (1.0 - cavitation_intensity * venturi.rbc_exposure_fraction.clamp(0.0, 1.0));

    let absolute_inlet_pressure = candidate.operating_point.absolute_inlet_pressure_pa();
    let collapse_ratio = (absolute_inlet_pressure / BLOOD_VAPOR_PRESSURE_PA.max(1.0))
        .powf((BUBBLE_POLYTROPIC_K - 1.0) / BUBBLE_POLYTROPIC_K);
    let collapse_ref = (SONO_REF_P_ABS_PA / BLOOD_VAPOR_PRESSURE_PA.max(1.0))
        .powf((BUBBLE_POLYTROPIC_K - 1.0) / BUBBLE_POLYTROPIC_K);
    let sonoluminescence_proxy =
        (cavitation_potential * (collapse_ratio / collapse_ref)).clamp(0.0, 1.0);

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

    let optical_path_length_405_m = solve
        .channel_samples
        .iter()
        .filter(|sample| sample.is_treatment_channel)
        .map(|sample| sample.cross_section.dims().1)
        .reduce(f64::max)
        .unwrap_or(0.0);
    let blue_light_delivery_index_405nm = (-BLOOD_ATTENUATION_405NM_INV_M
        * optical_path_length_405_m
        * candidate.operating_point.feed_hematocrit.max(0.0))
    .exp()
    .clamp(0.0, 1.0);

    metrics.cavitation_number = cavitation_number;
    metrics.cavitation_potential = cavitation_potential;
    metrics.throat_shear_rate_inv_s = max_venturi_shear_rate_inv_s;
    metrics.throat_shear_pa = safety.max_venturi_shear_pa;
    metrics.throat_exceeds_fda = safety.max_venturi_shear_pa > FDA_MAX_WALL_SHEAR_PA;
    metrics.max_main_channel_shear_pa = safety.max_main_channel_shear_pa;
    metrics.fda_main_compliant = safety.max_main_channel_shear_pa <= FDA_MAX_WALL_SHEAR_PA;
    metrics.bulk_hemolysis_index_per_pass = bulk_hi;
    metrics.hemolysis_index_per_pass = corrected_hi;
    metrics.flow_uniformity = solve.flow_uniformity;
    metrics.well_coverage_fraction = treatment_fraction;
    metrics.mean_residence_time_s = solve.mean_residence_time_s;
    metrics.total_pressure_drop_pa = safety.pressure_drop_pa;
    metrics.total_path_length_mm = total_path_length_mm;
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
    metrics.total_ecv_ml = total_volume_m3 * 1.0e6;
    metrics.flow_rate_ml_min = flow_rate_ml_min;
    metrics.plate_fits =
        topology.box_dims_mm.0 <= PLATE_WIDTH_MM && topology.box_dims_mm.1 <= PLATE_HEIGHT_MM;
    metrics.n_outlet_ports = candidate
        .blueprint
        .nodes
        .iter()
        .filter(|node| matches!(node.kind, cfd_schematics::domain::model::NodeKind::Outlet))
        .count();
    metrics.main_channel_shear_rate_inv_s =
        safety.max_main_channel_shear_pa / BLOOD_VISCOSITY_PA_S.max(1.0e-18);
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
    metrics.venturi_treatment_enabled = active_venturi_throat_count > 0;
    metrics.treatment_zone_mode = match topology.treatment_mode {
        TreatmentActuationMode::UltrasoundOnly => "UltrasoundOnly",
        TreatmentActuationMode::VenturiCavitation => "VenturiThroats",
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
    metrics.sonoluminescence_proxy = sonoluminescence_proxy;
    metrics.throat_transit_time_s = max_venturi_transit_time_s;
    metrics.fda_overall_compliant = metrics.fda_main_compliant
        && (safety.max_venturi_shear_pa <= FDA_MAX_WALL_SHEAR_PA
            || (max_venturi_transit_time_s <= FDA_TRANSIENT_TIME_S
                && safety.max_venturi_shear_pa <= FDA_TRANSIENT_SHEAR_PA));
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
    metrics.optical_path_length_405_m = optical_path_length_405_m;
    metrics.blue_light_delivery_index_405nm = blue_light_delivery_index_405nm;
    metrics.safety_margin_pa = FDA_MAX_WALL_SHEAR_PA - safety.max_main_channel_shear_pa;
    metrics.therapy_channel_fraction = treatment_fraction;
    metrics.wall_shear_p95_pa = wall_shear_p95_pa;
    metrics.wall_shear_p99_pa = wall_shear_p99_pa;
    metrics.wall_shear_mean_pa = wall_shear_mean_pa;
    metrics.wall_shear_cv = wall_shear_cv;
    metrics.fda_shear_percentile_compliant =
        wall_shear_p95_pa <= FDA_MAX_WALL_SHEAR_PA && wall_shear_p99_pa <= FDA_TRANSIENT_SHEAR_PA;
    metrics.mechanical_power_w = safety.pressure_drop_pa * flow_rate_m3_s;
    metrics.acoustic_capture_efficiency =
        (cavitation_potential * (safety.max_venturi_shear_pa / P_ATM_PA)).clamp(0.0, 1.0);
    metrics.specific_cavitation_energy_j_ml = if flow_rate_ml_min > 1.0e-18 {
        metrics.acoustic_capture_efficiency
            * metrics.mechanical_power_w
            * residence.treatment_residence_time_s
            / (flow_rate_ml_min / 60_000.0)
            * 1.0e3
    } else {
        0.0
    };
    metrics.hemolysis_index_per_pass_cavitation_amplified = corrected_hi;
    metrics.treatment_channel_hi = treatment_hi;
    metrics.bypass_channel_hi_mean = mean(&bypass_hi_values);
    metrics.bypass_channel_hi_max = bypass_hi_values.into_iter().fold(0.0, f64::max);
    metrics.per_channel_hemolysis = per_channel_hemolysis;
    metrics.acoustic_resonance_factor = channel_resonance_score;
    metrics.channel_resonance_score = channel_resonance_score;
    metrics.serial_cavitation_dose_fraction = serial_dose_fraction;
    metrics.treatment_zone_dwell_time_s = residence.treatment_residence_time_s;
    metrics.throat_temperature_rise_k = throat_temperature_rise_k;
    metrics.fda_thermal_compliant = throat_temperature_rise_k <= FDA_THROAT_TEMP_RISE_LIMIT_K;

    // ── Previously unset metrics ─────────────────────────────────────────────
    // Platelet activation index (Hellums 1994 power-law model).
    metrics.platelet_activation_index = pai_accumulator;

    // Diffuser pressure recovery: maximum across all venturi placements.
    metrics.diffuser_recovery_pa = venturi
        .placements
        .iter()
        .map(|p| p.diffuser_recovery_pa)
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

    // Topology-specific split-stage flow fractions via Hagen–Poiseuille Q ∝ w³.
    let (stage_center_fracs, model_frac) =
        split_stage_flow_fractions(&topology.split_stages);
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
    // These do NOT alter any existing metric fields used for scoring — they
    // are computed here for narrative reporting and future integration.
    let sdt_acoustic = compute_sdt_acoustic_metrics(
        cavitation_intensity,
        max_venturi_transit_time_s,
        safety.pressure_drop_pa,
    );

    // Sonosensitizer activation modulates cancer dose (Rosenthal 2004):
    // short throat transits (< 1 ms) yield incomplete activation.
    let cancer_dose_fraction =
        cancer_dose_fraction * sdt_acoustic.sensitizer_activation_efficiency.max(0.01);
    metrics.cancer_dose_fraction = cancer_dose_fraction;

    // Rayleigh-Plesset collapse amplification (Rayleigh 1917):
    // bubble collapse micro-jets increase hemolysis beyond the base shear model.
    if cavitation_potential > 0.0 {
        metrics.hemolysis_index_per_pass_cavitation_amplified *=
            sdt_acoustic.rayleigh_plesset_amplification;
    }

    // Acoustic energy density and contrast factors are available for
    // narrative reporting via sdt_acoustic.acoustic_energy_density_j_m3,
    // sdt_acoustic.ctc_contrast_factor, and sdt_acoustic.rbc_contrast_factor.

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
    pub acoustic_energy_density_j_m3: f64,

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
/// * `throat_transit_time_s` — residence time in the cavitation zone [s]
/// * `pressure_drop_pa` — total pressure drop across the chip [Pa]
pub fn compute_sdt_acoustic_metrics(
    cavitation_intensity: f64,
    throat_transit_time_s: f64,
    pressure_drop_pa: f64,
) -> SdtAcousticMetrics {
    // Sonosensitizer activation efficiency (Rosenthal 2004)
    let sensitizer_activation = sonosensitizer_activation_efficiency(
        SENSITIZER_K_ACT_HEMATOPORPHYRIN,
        cavitation_intensity,
        throat_transit_time_s,
    );

    // Rayleigh-Plesset collapse amplification (Rayleigh 1917)
    let r_max = 10.0e-6; // 10 µm typical cavitation bubble
    let r_0 = 1.0e-6; // 1 µm equilibrium nucleus
    let p_inf = pressure_drop_pa + 101_325.0; // gauge → absolute
    let rp_amplification = cavitation_hemolysis_amplification(r_max, r_0, p_inf);

    // Acoustic energy density (Gor'kov 1962)
    let e_acoustic = acoustic_energy_density(
        100_000.0, // 100 kPa typical pressure amplitude
        1060.0,    // blood density [kg/m³]
        1540.0,    // speed of sound in blood [m/s]
    );

    // Acoustic contrast factors for differential radiation force
    let ctc_contrast = acoustic_contrast_factor(RHO_CTC, RHO_PLASMA, KAPPA_CTC, KAPPA_PLASMA);
    let rbc_contrast = acoustic_contrast_factor(RHO_RBC, RHO_PLASMA, KAPPA_RBC, KAPPA_PLASMA);

    SdtAcousticMetrics {
        sensitizer_activation_efficiency: sensitizer_activation,
        rayleigh_plesset_amplification: rp_amplification,
        acoustic_energy_density_j_m3: e_acoustic,
        ctc_contrast_factor: ctc_contrast,
        rbc_contrast_factor: rbc_contrast,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
            m.acoustic_energy_density_j_m3.is_finite(),
            "acoustic_energy_density_j_m3 must be finite"
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
            m.sensitizer_activation_efficiency > 0.0
                && m.sensitizer_activation_efficiency < 1.0,
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
        assert!(m.acoustic_energy_density_j_m3 > 0.0);
    }

    #[test]
    fn test_sdt_acoustic_metrics_positive_energy() {
        // Energy density must always be positive (p₀² / (4ρc²) > 0)
        for pressure_drop in [0.0, 10_000.0, 100_000.0, 500_000.0] {
            let m = compute_sdt_acoustic_metrics(0.5, 0.001, pressure_drop);
            assert!(
                m.acoustic_energy_density_j_m3 > 0.0,
                "energy density must be positive at pressure_drop={pressure_drop}"
            );
        }
    }
}
