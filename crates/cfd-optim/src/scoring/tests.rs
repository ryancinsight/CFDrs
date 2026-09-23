use crate::constraints::HI_PASS_LIMIT;
use super::sigmoid_penalty;
use crate::metrics::SdtMetrics;
use super::score::score_candidate_impl;

use super::*;
use proptest::prelude::*;

fn base_metrics() -> SdtMetrics {
    let mut metrics = SdtMetrics {
        pressure_feasible: true,
        fda_main_compliant: true,
        plate_fits: true,
        hemolysis_index_per_pass: 2.0e-4,
        projected_hemolysis_15min_pediatric_3kg: 0.002,
        wbc_recovery: 0.70,
        rbc_pass_fraction: 0.30,
        wbc_purity: 0.70,
        total_ecv_ml: 0.2,
        flow_rate_ml_min: 120.0,
        three_pop_sep_efficiency: 0.30,
        cancer_targeted_cavitation: 0.20,
        oncology_selectivity_index: 0.10,
        cancer_rbc_cavitation_bias_index: 0.70,
        selective_cavitation_delivery_index: 0.30,
        rbc_venturi_protection: 0.40,
        sonoluminescence_proxy: 0.50,
        wbc_targeted_cavitation: 0.25,
        healthy_cell_protection_index: 0.0,
        ..SdtMetrics::default()
    };
    metrics.healthy_cell_protection_index = (1.0 - metrics.wbc_targeted_cavitation)
        .mul_add(metrics.rbc_venturi_protection, 0.0)
        .sqrt();
    metrics
}

#[test]
fn combined_mode_prefers_high_oncology_selectivity() {
    let mut low = base_metrics();
    low.cancer_targeted_cavitation = 0.04;
    low.oncology_selectivity_index = 0.0;

    let mut high = base_metrics();
    high.cancer_targeted_cavitation = 0.55;
    high.oncology_selectivity_index = 0.40;

    let mode = OptimMode::CombinedSdtLeukapheresis {
        leuka_weight: 0.5,
        sdt_weight: 0.5,
        patient_weight_kg: 3.0,
    };
    let w = SdtWeights::default();
    let s_low = score_candidate(&low, mode, &w);
    let s_high = score_candidate(&high, mode, &w);
    assert!(
        s_high > s_low,
        "combined mode should reward oncology-selective SDT candidates"
    );
}

#[test]
fn hydrosdt_uses_oncology_selectivity_signal() {
    let mut low = base_metrics();
    low.oncology_selectivity_index = 0.05;
    low.cancer_targeted_cavitation = 0.35;

    let mut high = base_metrics();
    high.oncology_selectivity_index = 0.40;
    high.cancer_targeted_cavitation = 0.35;

    let w = SdtWeights::default();
    let mode = OptimMode::HydrodynamicCavitationSDT;
    let s_low = score_candidate(&low, mode, &w);
    let s_high = score_candidate(&high, mode, &w);
    assert!(
        s_high > s_low,
        "HydroSDT score should increase with oncology selectivity at fixed cavitation"
    );
}

#[test]
fn hydrosdt_prefers_higher_cancer_rbc_cavitation_bias() {
    let mut low_bias = base_metrics();
    low_bias.cancer_rbc_cavitation_bias_index = 0.20;
    low_bias.selective_cavitation_delivery_index = 0.12;

    let mut high_bias = base_metrics();
    high_bias.cancer_rbc_cavitation_bias_index = 0.85;
    high_bias.selective_cavitation_delivery_index = 0.52;

    let w = SdtWeights::default();
    let mode = OptimMode::HydrodynamicCavitationSDT;
    let s_low = score_candidate(&low_bias, mode, &w);
    let s_high = score_candidate(&high_bias, mode, &w);
    assert!(
        s_high > s_low,
        "HydroSDT score should increase with stronger cancer-vs-RBC cavitation bias"
    );
}

#[test]
fn combined_mode_prefers_cif_remerge_near_outlet() {
    let mode = OptimMode::CombinedSdtLeukapheresis {
        leuka_weight: 0.5,
        sdt_weight: 0.5,
        patient_weight_kg: 3.0,
    };
    let w = SdtWeights::default();

    let mut far_remerge = base_metrics();
    far_remerge.cif_outlet_tail_length_mm = 5.5;
    far_remerge.cif_remerge_proximity_score = 0.10;

    let mut near_remerge = base_metrics();
    near_remerge.cif_outlet_tail_length_mm = 0.8;
    near_remerge.cif_remerge_proximity_score = 0.92;

    let s_far = score_candidate(&far_remerge, mode, &w);
    let s_near = score_candidate(&near_remerge, mode, &w);
    assert!(
        s_near > s_far,
        "combined mode should reward selective remerge-near-outlet layouts"
    );
}

#[test]
fn clotting_risk_penalizes_score() {
    let mode = OptimMode::CombinedSdtLeukapheresis {
        leuka_weight: 0.5,
        sdt_weight: 0.5,
        patient_weight_kg: 3.0,
    };
    let w = SdtWeights::default();

    let mut low_risk = base_metrics();
    low_risk.clotting_risk_index = 0.0;

    let mut high_risk = base_metrics();
    high_risk.clotting_risk_index = 1.0;

    let s_low = score_candidate(&low_risk, mode, &w);
    let s_high = score_candidate(&high_risk, mode, &w);
    assert!(
        s_low > s_high,
        "higher low-flow clotting risk should reduce score via coagulation penalty"
    );
}

#[test]
fn selective_acoustic_prefers_selective_center_routing() {
    let mut broad = base_metrics();
    broad.cancer_center_fraction = 0.22;
    broad.wbc_center_fraction = 0.18;
    broad.rbc_peripheral_fraction_three_pop = 0.24;
    broad.three_pop_sep_efficiency = 0.14;
    broad.therapy_channel_fraction = 0.16;
    broad.mean_residence_time_s = 0.20;
    broad.treatment_zone_dwell_time_s = 0.12;

    let mut selective = base_metrics();
    selective.cancer_center_fraction = 0.76;
    selective.wbc_center_fraction = 0.66;
    selective.rbc_peripheral_fraction_three_pop = 0.82;
    selective.three_pop_sep_efficiency = 0.61;
    selective.therapy_channel_fraction = 0.34;
    selective.mean_residence_time_s = 1.30;
    selective.treatment_zone_dwell_time_s = 1.05;

    let mode = OptimMode::SdtTherapy;
    let w = SdtWeights::default();
    let broad_score = score_candidate(&broad, mode, &w);
    let selective_score = score_candidate(&selective, mode, &w);
    assert!(
        selective_score > broad_score,
        "sdt therapy mode should reward center-lane enrichment and peripheral RBC routing"
    );
}

#[test]
fn sdt_cavitation_score_uses_physical_components_without_synthetic_floor() {
    let mut metrics = base_metrics();
    metrics.cavitation_number = 1.2;
    metrics.cavitation_potential = 0.0;
    metrics.well_coverage_fraction = 0.42;

    let weights = SdtWeights::default();
    let score = score_candidate(&metrics, OptimMode::SdtCavitation, &weights);

    let hi_ratio = metrics.hemolysis_index_per_pass / HI_PASS_LIMIT.max(1e-12);
    let hi_factor = 1.0_f64 / (1.0 + (hi_ratio / 0.5).powi(2));
    let expected = weights.cav_hemolysis * hi_factor + weights.cav_coverage * 0.42;

    let scale = expected.abs().max(1.0);
    assert!(
        (score - expected).abs() <= 1e-15 * scale,
        "expected {expected}, got {score}"
    );
}

#[test]
fn pediatric_leukapheresis_returns_zero_without_extracorporeal_volume() {
    let mut metrics = base_metrics();
    metrics.total_ecv_ml = 0.0;

    let score = score_candidate(
        &metrics,
        OptimMode::PediatricLeukapheresis {
            patient_weight_kg: 3.0,
        },
        &SdtWeights::default(),
    );

    assert_eq!(score, 0.0, "zero extracorporeal volume must score zero");
}

#[test]
fn hydrodynamic_cavitation_scores_zero_when_physical_terms_vanish() {
    let mut metrics = base_metrics();
    metrics.wbc_recovery = 1.0;
    metrics.three_pop_sep_efficiency = 0.0;
    metrics.cancer_targeted_cavitation = 0.0;
    metrics.oncology_selectivity_index = 0.0;
    metrics.cancer_rbc_cavitation_bias_index = 0.0;
    metrics.selective_cavitation_delivery_index = 0.0;
    metrics.rbc_venturi_protection = 0.0;
    metrics.sonoluminescence_proxy = 0.0;
    metrics.cif_outlet_tail_length_mm = 0.0;
    metrics.cif_remerge_proximity_score = 0.0;
    metrics.blue_light_delivery_index_405nm = 0.0;
    metrics.clotting_risk_index = 0.0;

    let score = score_candidate(
        &metrics,
        OptimMode::HydrodynamicCavitationSDT,
        &SdtWeights::default(),
    );

    assert_eq!(
        score, 0.0,
        "zero physical cavitation signal must score zero"
    );
}

#[test]
fn combined_leukapheresis_scores_zero_when_all_physical_signals_vanish() {
    let mut metrics = base_metrics();
    metrics.total_ecv_ml = 0.0;
    metrics.wbc_recovery = 0.0;
    metrics.rbc_pass_fraction = 1.0;
    metrics.wbc_purity = 0.0;
    metrics.three_pop_sep_efficiency = 0.0;
    metrics.cancer_targeted_cavitation = 0.0;
    metrics.oncology_selectivity_index = 0.0;
    metrics.cancer_rbc_cavitation_bias_index = 0.0;
    metrics.selective_cavitation_delivery_index = 0.0;
    metrics.rbc_venturi_protection = 0.0;
    metrics.sonoluminescence_proxy = 0.0;
    metrics.wbc_targeted_cavitation = 0.0;
    metrics.projected_hemolysis_15min_pediatric_3kg = 0.0;
    metrics.projected_hemolysis_15min_adult = 0.0;
    metrics.cif_outlet_tail_length_mm = 0.0;
    metrics.cif_remerge_proximity_score = 0.0;
    metrics.blue_light_delivery_index_405nm = 0.0;
    metrics.clotting_risk_index = 0.0;
    metrics.therapy_channel_fraction = 0.0;
    metrics.well_coverage_fraction = 0.0;
    metrics.mean_residence_time_s = 0.0;
    metrics.treatment_zone_dwell_time_s = 0.0;

    let score = score_candidate(
        &metrics,
        OptimMode::CombinedSdtLeukapheresis {
            leuka_weight: 0.5,
            sdt_weight: 0.5,
            patient_weight_kg: 3.0,
        },
        &SdtWeights::default(),
    );

    assert_eq!(score, 0.0, "zero physical signal must score zero");
}

#[test]
fn smooth_penalty_plate_overflow_scores_zero() {
    let mut metrics = base_metrics();
    metrics.plate_fits = false;

    let score = score_candidate_impl(
        &metrics,
        OptimMode::SdtCavitation,
        &SdtWeights::default(),
        ScoreMode::SmoothPenalty,
        100.0,
    );

    assert_eq!(
        score, 0.0,
        "plate overflow must score zero under smooth penalty"
    );
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
fn prop_hard_constraints_produce_zero_score(
        cav in 0.0_f64..1.0,
        sep3 in 0.0_f64..1.0,
        wbc_recovery in 0.0_f64..1.0,
        optical_405 in 0.0_f64..1.0
    ) {
        let modes = [
            OptimMode::SdtCavitation,
            OptimMode::UniformExposure,
            OptimMode::Combined { cavitation_weight: 0.5, exposure_weight: 0.5 },
            OptimMode::CellSeparation,
            OptimMode::ThreePopSeparation,
            OptimMode::SdtTherapy,
            OptimMode::HydrodynamicCavitationSDT,
            OptimMode::CombinedSdtLeukapheresis { leuka_weight: 0.5, sdt_weight: 0.5, patient_weight_kg: 3.0 },
            OptimMode::RbcProtectedSdt,
        ];
        let weights = SdtWeights::default();

        for mode in modes {
            let mut m = base_metrics();
            m.cavitation_potential = cav;
            m.three_pop_sep_efficiency = sep3;
            m.wbc_recovery = wbc_recovery;
            m.blue_light_delivery_index_405nm = optical_405;

            m.pressure_feasible = false;
            let s = score_candidate(&m, mode, &weights);
            prop_assert_eq!(s, INFEASIBILITY_SCORE,
                "pressure-infeasible candidate must receive zero score, got {}", s);

            let mut m = base_metrics();
            m.cavitation_potential = cav;
            m.three_pop_sep_efficiency = sep3;
            m.wbc_recovery = wbc_recovery;
            m.blue_light_delivery_index_405nm = optical_405;
            m.fda_main_compliant = false;
            let s = score_candidate(&m, mode, &weights);
            prop_assert_eq!(s, INFEASIBILITY_SCORE,
                "FDA-noncompliant candidate must receive zero score, got {}", s);

            let mut m = base_metrics();
            m.cavitation_potential = cav;
            m.three_pop_sep_efficiency = sep3;
            m.wbc_recovery = wbc_recovery;
            m.blue_light_delivery_index_405nm = optical_405;
            m.plate_fits = false;
            let s = score_candidate(&m, mode, &weights);
            prop_assert_eq!(s, INFEASIBILITY_SCORE,
                "plate-overflow candidate must receive zero score, got {}", s);
        }
    }

    #[test]
    fn prop_score_candidate_bounded_on_feasible_inputs(
        cav in 0.0_f64..1.0,
        sep3 in 0.0_f64..1.0,
        wbc_recovery in 0.0_f64..1.0,
        rbc_pass in 0.0_f64..1.0,
        optical_405 in 0.0_f64..1.0,
        clotting_risk in 0.0_f64..1.0
    ) {
        let modes = [
            OptimMode::SdtCavitation,
            OptimMode::UniformExposure,
            OptimMode::Combined { cavitation_weight: 0.5, exposure_weight: 0.5 },
            OptimMode::CellSeparation,
            OptimMode::ThreePopSeparation,
            OptimMode::SdtTherapy,
            OptimMode::HydrodynamicCavitationSDT,
            OptimMode::CombinedSdtLeukapheresis { leuka_weight: 0.5, sdt_weight: 0.5, patient_weight_kg: 3.0 },
            OptimMode::RbcProtectedSdt,
        ];
        let weights = SdtWeights::default();

        let mut m = base_metrics();
        m.pressure_feasible = true;
        m.fda_main_compliant = true;
        m.plate_fits = true;
        m.cavitation_potential = cav;
        m.three_pop_sep_efficiency = sep3;
        m.wbc_recovery = wbc_recovery;
        m.rbc_pass_fraction = rbc_pass;
        m.blue_light_delivery_index_405nm = optical_405;
        m.clotting_risk_index = clotting_risk;

        for mode in modes {
            let score = score_candidate(&m, mode, &weights);
            prop_assert!(score.is_finite());
            prop_assert!(score >= 0.0, "score = {score}");
            prop_assert!(score <= 1.0, "score = {score}");
        }
    }

    #[test]
    fn prop_sigmoid_penalty_bounded_and_monotone(
        m1 in -2.0_f64..2.0,
        m2 in -2.0_f64..2.0
    ) {
        let s1 = sigmoid_penalty(m1);
        let s2 = sigmoid_penalty(m2);
        prop_assert!((0.0..=1.0).contains(&s1));
        prop_assert!((0.0..=1.0).contains(&s2));
        if m1 <= m2 {
            prop_assert!(s1 <= s2 + 1e-12, "m1={m1}, m2={m2}, s1={s1}, s2={s2}");
        } else {
            prop_assert!(s2 <= s1 + 1e-12, "m1={m1}, m2={m2}, s1={s1}, s2={s2}");
        }
    }
}

/// Sigmoid penalty boundary values.
///
/// `sigmoid_penalty(m) = clamp(0.5 + 5m, 0, 1)`:
/// - At boundary (m=0): 0.5 (50% penalty)
/// - At m=+0.1: 1.0 (fully feasible)
/// - At m=−0.1: 0.0 (fully infeasible)
#[test]
fn sigmoid_penalty_boundary_values() {
    assert_eq!(sigmoid_penalty(0.0), 0.5);
    assert_eq!(sigmoid_penalty(0.1), 1.0);
    assert_eq!(sigmoid_penalty(-0.1), 0.0);
    // Deep interior / exterior
    assert_eq!(sigmoid_penalty(1.0), 1.0);
    assert_eq!(sigmoid_penalty(-1.0), 0.0);
}
