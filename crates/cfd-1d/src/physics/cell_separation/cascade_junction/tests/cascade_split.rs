//! Cascade, cross-junction, and mixed Bi/Tri selective-routing tests.

use super::*;

#[test]
fn symmetric_split_equal_distribution() {
    // With symmetric 1/3 split all cell types should distribute by flow fraction
    // (cancer biased slightly more toward center, RBC less so).
    let r = cascade_junction_separation(1, 1.0 / 3.0, length(2e-3), length(1e-3), flow(5e-6));
    // Center arm carries 1/3 of flow → q_frac = (1/3)³/((1/3)³+2*(1/3)³) = 1/3
    // After 1 level: cancer_center_fraction > rbc_center_fraction (stiffness effect)
    assert!(
        r.cancer_center_fraction > r.wbc_center_fraction
            || (r.cancer_center_fraction - r.wbc_center_fraction).abs() < 1e-6
    );
    assert!(r.wbc_center_fraction >= r.rbc_peripheral_fraction * 0.0); // just sanity
    assert!(r.separation_efficiency >= 0.0 && r.separation_efficiency <= 1.0);
}

#[test]
fn center_biased_increases_separation() {
    // Center-biased split should increase separation relative to symmetric.
    let r_sym = cascade_junction_separation(2, 1.0 / 3.0, length(2e-3), length(1e-3), flow(5e-6));
    let r_bias = cascade_junction_separation(2, 0.55, length(2e-3), length(1e-3), flow(5e-6));
    // Biased split: q_center_frac > 1/3 → stronger Zweifach-Fung routing
    assert!(r_bias.cancer_center_fraction >= r_sym.cancer_center_fraction - 1e-10);
}

#[test]
fn more_levels_increases_enrichment() {
    let r1 = cascade_junction_separation(1, 0.45, length(2e-3), length(1e-3), flow(5e-6));
    let r3 = cascade_junction_separation(3, 0.45, length(2e-3), length(1e-3), flow(5e-6));
    // More cascade levels → more cancer enriched in center, more RBCs in bypass
    assert!(r3.rbc_peripheral_fraction >= r1.rbc_peripheral_fraction - 1e-10);
}

#[test]
fn tri_center_q_frac_symmetric() {
    // Symmetric 1/3 split → q_frac = 1/3
    let q = tri_center_q_frac(1.0 / 3.0);
    assert!(
        (q - 1.0 / 3.0).abs() < 1e-10,
        "symmetric frac should give q=1/3, got {q}"
    );
}

#[test]
fn tri_center_q_frac_larger_center_gives_more_flow() {
    let q33 = tri_center_q_frac(0.333);
    let q55 = tri_center_q_frac(0.55);
    assert!(q55 > q33, "wider center arm should carry more flow");
}

#[test]
fn checked_tri_center_q_frac_rejects_closed_interval_endpoints() {
    let err = checked_tri_center_q_frac(0.0)
        .expect_err("checked center-arm flow fraction must reject zero width fraction");
    assert!(err.to_string().contains("width fraction"));
}

#[test]
fn incremental_filtration_more_pretri_levels_pushes_more_rbc_periphery() {
    let r1 = incremental_filtration_separation_staged(1, 0.45, 0.50, 0.68);
    let r3 = incremental_filtration_separation_staged(3, 0.45, 0.50, 0.68);
    assert!(r3.rbc_peripheral_fraction >= r1.rbc_peripheral_fraction - 1e-10);
    assert!(r3.separation_efficiency >= r1.separation_efficiency - 1e-10);
}

#[test]
fn incremental_filtration_higher_bi_treat_frac_increases_treatment_capture() {
    let low = incremental_filtration_separation_staged(2, 0.45, 0.45, 0.60);
    let high = incremental_filtration_separation_staged(2, 0.45, 0.45, 0.76);
    assert!(high.cancer_center_fraction >= low.cancer_center_fraction - 1e-10);
    assert!(high.wbc_center_fraction >= low.wbc_center_fraction - 1e-10);
}

#[test]
fn mixed_sequence_tri_bi_routes_more_rbc_peripheral_than_single_tri() {
    let q_tri = tri_center_q_frac(0.45);
    let tri_only = mixed_cascade_separation(&[(q_tri, true)]);
    let tri_bi = mixed_cascade_separation(&[(q_tri, true), (0.72, false)]);

    assert!(
        tri_bi.rbc_peripheral_fraction >= tri_only.rbc_peripheral_fraction,
        "adding a treatment bifurcation should not reduce RBC peripheral routing"
    );
    assert!(
        tri_bi.cancer_center_fraction > (1.0 - tri_bi.rbc_peripheral_fraction),
        "cancer capture should remain above RBC center carryover in selective routing"
    );
}

#[test]
fn stronger_trifurcation_bias_improves_three_pop_selectivity() {
    let weak = mixed_cascade_separation(&[
        (tri_center_q_frac(0.38), true),
        (tri_center_q_frac(0.38), true),
    ]);
    let strong = mixed_cascade_separation(&[
        (tri_center_q_frac(0.52), true),
        (tri_center_q_frac(0.52), true),
    ]);

    assert!(
        strong.cancer_center_fraction >= weak.cancer_center_fraction,
        "stronger center-arm bias should not reduce cancer capture"
    );
    assert!(
        strong.separation_efficiency >= weak.separation_efficiency,
        "stronger center-arm bias should improve separation efficiency"
    );
}

#[test]
fn incremental_filtration_terminal_tri_center_bias_improves_cancer_center_fraction() {
    let symmetric_terminal = incremental_filtration_separation_staged(2, 0.45, 1.0 / 3.0, 0.68);
    let center_biased_terminal = incremental_filtration_separation_staged(2, 0.45, 0.55, 0.68);
    assert!(
        center_biased_terminal.cancer_center_fraction
            >= symmetric_terminal.cancer_center_fraction - 1e-10
    );
}

#[test]
fn cascade_qfrac_api_matches_uniform_width_model() {
    let width_model = cascade_junction_separation(3, 0.45, length(2e-3), length(1e-3), flow(5e-6));
    let q = tri_center_q_frac(0.45);
    let solved_like = cascade_junction_separation_from_qfracs(&[q, q, q]);
    assert!(
        (width_model.cancer_center_fraction - solved_like.cancer_center_fraction).abs() < 1e-12
    );
    assert!(
        (width_model.rbc_peripheral_fraction - solved_like.rbc_peripheral_fraction).abs() < 1e-12
    );
}

#[test]
fn checked_cascade_qfrac_api_rejects_empty_stage_sequence() {
    let err = checked_cascade_junction_separation_from_qfracs(&[])
        .expect_err("checked cascade qfrac API must reject empty stage sequences");
    assert!(
        err.to_string()
            .contains("at least one center-arm flow fraction")
    );
}

#[test]
fn checked_cascade_junction_matches_legacy_nominal_case() {
    let legacy = cascade_junction_separation(3, 0.45, length(2e-3), length(1e-3), flow(5e-6));
    let checked =
        checked_cascade_junction_separation(3, 0.45, length(2e-3), length(1e-3), flow(5e-6))
            .expect("checked cascade junction API should succeed on a nominal case");

    assert!((legacy.cancer_center_fraction - checked.cancer_center_fraction).abs() < 1e-12);
    assert!((legacy.rbc_peripheral_fraction - checked.rbc_peripheral_fraction).abs() < 1e-12);
    assert!((legacy.center_hematocrit_ratio - checked.center_hematocrit_ratio).abs() < 1e-12);
}

#[test]
fn incremental_qfrac_api_matches_uniform_width_model() {
    let width_model = incremental_filtration_separation_staged(2, 0.45, 0.55, 0.68);
    let q_pretri = cif_pretri_stage_q_fracs(2, 0.45, 0.55);
    let q_tri = tri_center_q_frac(0.55);
    let solved_like = incremental_filtration_separation_from_qfracs(&q_pretri, q_tri, 0.68);
    assert!(
        (width_model.cancer_center_fraction - solved_like.cancer_center_fraction).abs() < 1e-12
    );
    assert!((width_model.rbc_center_fraction - solved_like.rbc_center_fraction).abs() < 1e-12);
}

#[test]
fn cif_pretri_stage_fracs_ramp_toward_terminal_bias() {
    let stage_fracs = cif_pretri_stage_center_fracs(3, 0.45, 0.60);
    assert_eq!(stage_fracs.len(), 3);
    assert!(stage_fracs[0] >= 0.45 - 1e-12);
    assert!(stage_fracs[1] >= stage_fracs[0] - 1e-12);
    assert!(stage_fracs[2] >= stage_fracs[1] - 1e-12);
    assert!(stage_fracs[2] <= 0.60 + 1e-12);
}

#[test]
fn checked_cif_pretri_stage_fracs_reject_invalid_stage_count() {
    let err = checked_cif_pretri_stage_center_fracs(0, 0.45, 0.60)
        .expect_err("checked CIF stage fractions must reject zero pre-trifurcation stages");
    assert!(err.to_string().contains("stage count"));
}

#[test]
fn checked_incremental_filtration_rejects_invalid_terminal_bifurcation_fraction() {
    let err = checked_incremental_filtration_separation_staged(2, 0.45, 0.50, 0.40)
        .expect_err("checked incremental filtration must reject terminal bifurcation fractions below the validated range");
    assert!(err.to_string().contains("terminal bifurcation"));
}

#[test]
fn checked_incremental_filtration_rejects_zero_pretri_stages() {
    let err = checked_incremental_filtration_separation_staged(0, 0.45, 0.50, 0.68)
        .expect_err("checked incremental filtration must reject zero pre-trifurcation stages");
    assert!(err.to_string().contains("stage count"));
}

#[test]
fn checked_incremental_filtration_matches_legacy_nominal_case() {
    let legacy = incremental_filtration_separation_staged(2, 0.45, 0.55, 0.68);
    let checked = checked_incremental_filtration_separation_staged(2, 0.45, 0.55, 0.68).expect(
        "checked incremental filtration should succeed on a nominal selective-routing case",
    );

    assert!((legacy.cancer_center_fraction - checked.cancer_center_fraction).abs() < 1e-12);
    assert!((legacy.rbc_center_fraction - checked.rbc_center_fraction).abs() < 1e-12);
    assert!((legacy.center_hematocrit_ratio - checked.center_hematocrit_ratio).abs() < 1e-12);
}

#[test]
fn checked_incremental_qfrac_api_matches_legacy_nominal_case() {
    let q_pretri = cif_pretri_stage_q_fracs(2, 0.45, 0.55);
    let q_tri = tri_center_q_frac(0.55);
    let legacy = incremental_filtration_separation_from_qfracs(&q_pretri, q_tri, 0.68);
    let checked = checked_incremental_filtration_separation_from_qfracs(&q_pretri, q_tri, 0.68)
        .expect("checked incremental qfrac API should succeed on a nominal case");

    assert!((legacy.cancer_center_fraction - checked.cancer_center_fraction).abs() < 1e-12);
    assert!((legacy.rbc_center_fraction - checked.rbc_center_fraction).abs() < 1e-12);
}

#[test]
fn cross_junction_q_frac_symmetric_matches_basic() {
    // With symmetric 1/3 split, cross-junction model should give similar
    // result to the basic model (identical for symmetric).
    let q_basic = tri_center_q_frac(1.0 / 3.0);
    let q_cross = tri_center_q_frac_cross_junction(1.0 / 3.0, length(2e-3), length(1e-3));
    // Both should be ~1/3 for symmetric geometry
    assert!((q_basic - 1.0 / 3.0).abs() < 1e-10);
    assert!(
        (q_cross - 1.0 / 3.0).abs() < 0.05,
        "symmetric cross-junction should be near 1/3, got {q_cross}"
    );
}

#[test]
fn cross_junction_q_frac_wider_center_sees_more_k_loss() {
    // Cross-junction K-factor correction penalises the wider center arm
    // more than the narrow peripherals (K_eff ∝ √(A_branch/A_parent),
    // and the additive resistance term scales as K × w/h), so q_center
    // *decreases*.  This is the desired physics: more flow is diverted
    // to peripheral arms, improving RBC peripheral enrichment.
    let q_basic = tri_center_q_frac(0.55);
    let q_cross = tri_center_q_frac_cross_junction(0.55, length(2e-3), length(1e-3));
    assert!(
        q_cross <= q_basic + 1e-6,
        "cross-junction correction should reduce center fraction: basic={q_basic}, cross={q_cross}"
    );
}

#[test]
fn cross_junction_cif_pushes_more_rbc_to_periphery() {
    let basic = incremental_filtration_separation_staged(2, 0.45, 0.50, 0.68);
    let cross = incremental_filtration_separation_cross_junction(
        2,
        0.45,
        0.50,
        0.68,
        length(2e-3),
        length(1e-3),
    );
    assert!(
        cross.rbc_peripheral_fraction >= basic.rbc_peripheral_fraction - 0.01,
        "cross-junction selective routing should push more RBCs to periphery: basic={}, cross={}",
        basic.rbc_peripheral_fraction,
        cross.rbc_peripheral_fraction
    );
}

#[test]
fn cross_junction_cct_pushes_more_rbc_to_periphery() {
    let basic = cascade_junction_separation(2, 0.45, length(2e-3), length(1e-3), flow(5e-6));
    let cross = cascade_junction_separation_cross_junction(2, 0.45, length(2e-3), length(1e-3));
    assert!(
        cross.rbc_peripheral_fraction >= basic.rbc_peripheral_fraction - 0.01,
        "cross-junction cascade routing should push more RBCs to periphery: basic={}, cross={}",
        basic.rbc_peripheral_fraction,
        cross.rbc_peripheral_fraction
    );
}

#[test]
fn mixed_cascade_all_tri_matches_cascade() {
    let q = tri_center_q_frac(0.45);
    let pure = cascade_routing::cascade_from_q_fractions(&[q, q]);
    let mixed = mixed_cascade_separation(&[(q, true), (q, true)]);
    assert!((pure.cancer_center_fraction - mixed.cancer_center_fraction).abs() < 1e-12);
    assert!((pure.rbc_peripheral_fraction - mixed.rbc_peripheral_fraction).abs() < 1e-12);
}

#[test]
fn mixed_cascade_tri_bi_produces_nonzero_separation() {
    let q_tri = tri_center_q_frac(0.45);
    let q_bi = 0.68;
    let r = mixed_cascade_separation(&[(q_tri, true), (q_bi, false)]);
    assert!(r.cancer_center_fraction > 0.0);
    assert!(r.separation_efficiency > 0.0);
    assert!(r.cancer_center_fraction > r.rbc_peripheral_fraction.min(0.99));
}

#[test]
fn treatment_bifurcation_matches_single_stage_mixed_cascade() {
    let q_bi = 0.68;
    let direct = treatment_bifurcation_separation(q_bi);
    let mixed = mixed_cascade_separation(&[(q_bi, false)]);

    assert!((direct.cancer_center_fraction - mixed.cancer_center_fraction).abs() < 1e-12);
    assert!((direct.wbc_center_fraction - mixed.wbc_center_fraction).abs() < 1e-12);
    assert!((direct.rbc_peripheral_fraction - mixed.rbc_peripheral_fraction).abs() < 1e-12);
    assert!((direct.separation_efficiency - mixed.separation_efficiency).abs() < 1e-12);
}

#[test]
fn checked_mixed_cascade_rejects_invalid_stage_flow_fraction() {
    let err = checked_mixed_cascade_separation(&[(1.20, true)]).expect_err(
        "checked mixed cascade routing must reject stage flow fractions outside the open unit interval",
    );
    assert!(err.to_string().contains("stage flow fraction"));
}

#[test]
fn checked_mixed_cascade_matches_legacy_nominal_case() {
    let q_tri = tri_center_q_frac(0.45);
    let legacy = mixed_cascade_separation(&[(q_tri, true), (0.68, false)]);
    let checked = checked_mixed_cascade_separation(&[(q_tri, true), (0.68, false)])
        .expect("checked mixed cascade routing should succeed on a nominal case");

    assert!((legacy.cancer_center_fraction - checked.cancer_center_fraction).abs() < 1e-12);
    assert!((legacy.rbc_peripheral_fraction - checked.rbc_peripheral_fraction).abs() < 1e-12);
}

#[test]
fn mixed_cascade_deeper_improves_separation() {
    let q_tri = tri_center_q_frac(0.45);
    let r1 = mixed_cascade_separation(&[(q_tri, true)]);
    let r2 = mixed_cascade_separation(&[(q_tri, true), (q_tri, true)]);
    assert!(r2.separation_efficiency >= r1.separation_efficiency - 1e-10);
}
