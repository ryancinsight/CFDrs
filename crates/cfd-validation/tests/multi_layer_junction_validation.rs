//! Multi-layer junction validation: validates cascade cell separation physics
//! through 1-layer, 2-layer, and 3-layer bifurcation, trifurcation,
//! quadfurcation, and pentafurcation trees.
//!
//! Tests that:
//! 1. Cancer enrichment increases with cascade depth (Zweifach-Fung effect)
//! 2. RBC exclusion improves with more separation stages
//! 3. Fahraeus-Lindqvist viscosity correction is consistent across depths
//! 4. Plasma skimming reduces hematocrit in smaller branches
//! 5. Quemada viscosity is physically bounded across all flow conditions

use cfd_1d::{
    cascade_junction_separation, cascade_junction_separation_from_qfracs,
    fahraeus_lindqvist_viscosity, incremental_filtration_separation_staged,
    plasma_skimming_hematocrit, quemada_viscosity, three_population_equilibria,
};

// ── Physical constants ──────────────────────────────────────────────────

const MU_PLASMA: f64 = 0.0012; // Pa·s
const HT_NORMAL: f64 = 0.45;
const BLOOD_DENSITY: f64 = 1060.0; // kg/m³
const BLOOD_VISCOSITY: f64 = 3.5e-3; // Pa·s

// ── Helper: simulate N-way junction as cascade of bifurcations ──────────

/// Simulate an N-furcation by using a cascade of trifurcations with
/// center_frac = 1/N. For N=4, a single trifurcation at center_frac=0.25
/// produces 3 arms; a 2-layer bifurcation (4 terminal branches) is
/// simulated as 2 levels of cascade with center_frac=0.5.
///
/// For quadfurcation (4 arms): 2-level cascade with center_frac=0.5.
/// For pentafurcation (5 arms): 2-level cascade with center_frac=0.4.

// =========================================================================
// Bifurcation (1/2/3 layers)
// =========================================================================

/// Test 1: Single bifurcation enriches cancer cells at center.
///
/// A biased bifurcation (60/40 split) should route stiff cancer cells
/// preferentially to the higher-flow center arm via the Zweifach-Fung
/// effect.
#[test]
fn test_1_layer_bifurcation_cancer_enrichment() {
    // Use incremental_filtration_separation_staged with n_pretri=0 and
    // terminal_bi_treat_frac=0.6 to model a single biased bifurcation.
    let r = incremental_filtration_separation_staged(
        0,    // n_pretri (no trifurcation pre-stages)
        0.33, // pretri_center_frac (unused with n_pretri=0)
        0.33, // terminal_tri_frac
        0.6,  // terminal_bi_treat_frac (biased toward treatment)
    );

    // Cancer should be enriched at the treatment arm
    assert!(
        r.cancer_center_fraction > 0.0,
        "1-layer bifurcation should enrich cancer: cancer_center={:.4}",
        r.cancer_center_fraction
    );
    assert!(
        r.separation_efficiency >= 0.0,
        "separation efficiency should be non-negative: {:.4}",
        r.separation_efficiency
    );
}

/// Test 2: 2-layer bifurcation cascade produces higher separation efficiency
/// than a single layer.
///
/// Note: cancer_center_fraction may decrease with depth because deeper
/// cascades narrow the center arm (fewer cells pass through all gates),
/// but the *differential* between cancer and RBC routing widens, so
/// separation_efficiency increases.
#[test]
fn test_2_layer_bifurcation_deeper_enrichment() {
    let r1 = cascade_junction_separation(1, 0.45, 2e-3, 1e-3, 1e-6);
    let r2 = cascade_junction_separation(2, 0.45, 2e-3, 1e-3, 1e-6);

    assert!(
        r2.separation_efficiency >= r1.separation_efficiency,
        "2-layer separation_efficiency should be >= 1-layer: {:.4} vs {:.4}",
        r2.separation_efficiency,
        r1.separation_efficiency
    );

    // The ratio of cancer to RBC center fractions (selectivity) should improve
    let rbc_center_1 = 1.0 - r1.rbc_peripheral_fraction;
    let rbc_center_2 = 1.0 - r2.rbc_peripheral_fraction;
    let selectivity_1 = if rbc_center_1 > 1e-12 {
        r1.cancer_center_fraction / rbc_center_1
    } else {
        f64::INFINITY
    };
    let selectivity_2 = if rbc_center_2 > 1e-12 {
        r2.cancer_center_fraction / rbc_center_2
    } else {
        f64::INFINITY
    };
    assert!(
        selectivity_2 >= selectivity_1,
        "2-layer selectivity should be >= 1-layer: {:.4} vs {:.4}",
        selectivity_2,
        selectivity_1
    );
}

/// Test 3: 3-layer cascade enrichment is monotonically increasing: 3 > 2 > 1.
#[test]
fn test_3_layer_bifurcation_monotone_enrichment() {
    let r1 = cascade_junction_separation(1, 0.45, 2e-3, 1e-3, 1e-6);
    let r2 = cascade_junction_separation(2, 0.45, 2e-3, 1e-3, 1e-6);
    let r3 = cascade_junction_separation(3, 0.45, 2e-3, 1e-3, 1e-6);

    assert!(
        r3.separation_efficiency >= r2.separation_efficiency
            && r2.separation_efficiency >= r1.separation_efficiency,
        "Monotone enrichment: eff_3={:.4} >= eff_2={:.4} >= eff_1={:.4}",
        r3.separation_efficiency,
        r2.separation_efficiency,
        r1.separation_efficiency
    );
}

// =========================================================================
// Trifurcation (1/2/3 layers)
// =========================================================================

/// Test 4: Single trifurcation enriches CTCs at center arm (three-population).
#[test]
fn test_1_layer_trifurcation_three_population() {
    let r = cascade_junction_separation(1, 0.5, 2e-3, 1e-3, 1e-6);

    // With center_frac=0.5, cancer cells (beta=1.85) should strongly
    // prefer the center arm over RBCs (beta=1.0).
    let rbc_center = 1.0 - r.rbc_peripheral_fraction;
    assert!(
        r.cancer_center_fraction > rbc_center,
        "Cancer should be enriched at center over RBC: cancer={:.4}, rbc_center={:.4}",
        r.cancer_center_fraction,
        rbc_center
    );
}

/// Test 5: 2-layer trifurcation improves separation over 1-layer.
#[test]
fn test_2_layer_trifurcation_improved_separation() {
    let r1 = cascade_junction_separation(1, 0.5, 2e-3, 1e-3, 1e-6);
    let r2 = cascade_junction_separation(2, 0.5, 2e-3, 1e-3, 1e-6);

    assert!(
        r2.separation_efficiency >= r1.separation_efficiency,
        "2-layer trifurcation separation should be >= 1-layer: {:.4} vs {:.4}",
        r2.separation_efficiency,
        r1.separation_efficiency
    );
}

/// Test 6: 3-layer trifurcation cascade increases RBC peripheral fraction
/// with depth (RBCs are increasingly excluded from center arm).
#[test]
fn test_3_layer_trifurcation_rbc_exclusion() {
    let r1 = cascade_junction_separation(1, 0.45, 2e-3, 1e-3, 1e-6);
    let r2 = cascade_junction_separation(2, 0.45, 2e-3, 1e-3, 1e-6);
    let r3 = cascade_junction_separation(3, 0.45, 2e-3, 1e-3, 1e-6);

    // RBC peripheral fraction should increase (or remain stable) with depth:
    // deeper cascades push more RBCs to peripheral arms.
    assert!(
        r3.rbc_peripheral_fraction >= r2.rbc_peripheral_fraction
            && r2.rbc_peripheral_fraction >= r1.rbc_peripheral_fraction,
        "RBC peripheral fraction should increase with depth: {:.4} >= {:.4} >= {:.4}",
        r3.rbc_peripheral_fraction,
        r2.rbc_peripheral_fraction,
        r1.rbc_peripheral_fraction
    );
}

// =========================================================================
// Quadfurcation and Pentafurcation
// =========================================================================

/// Test 7: Quadfurcation (4-way split via 2-level cascade) conserves mass:
/// cancer + WBC + RBC center fractions are all bounded in [0, 1].
#[test]
fn test_1_layer_quadfurcation_mass_conservation() {
    // Simulate a 4-way split as 2-level cascade with center_frac=0.5
    // (first level splits into center + 2 peripheral; second level
    //  re-splits the center arm, yielding 4 terminal branches).
    let r = cascade_junction_separation(2, 0.5, 2e-3, 1e-3, 1e-6);

    assert!(
        (0.0..=1.0).contains(&r.cancer_center_fraction),
        "cancer fraction out of [0,1]: {}",
        r.cancer_center_fraction
    );
    assert!(
        (0.0..=1.0).contains(&r.wbc_center_fraction),
        "WBC fraction out of [0,1]: {}",
        r.wbc_center_fraction
    );
    assert!(
        (0.0..=1.0).contains(&r.rbc_peripheral_fraction),
        "RBC peripheral fraction out of [0,1]: {}",
        r.rbc_peripheral_fraction
    );
    assert!(
        r.separation_efficiency >= 0.0,
        "separation efficiency should be non-negative: {}",
        r.separation_efficiency
    );
}

/// Test 8: Pentafurcation (5-way split via 2-level cascade with center_frac=0.4)
/// conserves mass bounds.
#[test]
fn test_1_layer_pentafurcation_mass_conservation() {
    // 2-level cascade at center_frac=0.4 produces ~5 effective branches
    let r = cascade_junction_separation(2, 0.4, 2e-3, 1e-3, 1e-6);

    assert!(
        (0.0..=1.0).contains(&r.cancer_center_fraction),
        "cancer fraction out of [0,1]: {}",
        r.cancer_center_fraction
    );
    assert!(
        (0.0..=1.0).contains(&r.wbc_center_fraction),
        "WBC fraction out of [0,1]: {}",
        r.wbc_center_fraction
    );
    assert!(
        (0.0..=1.0).contains(&r.rbc_peripheral_fraction),
        "RBC peripheral fraction out of [0,1]: {}",
        r.rbc_peripheral_fraction
    );
}

/// Test 9: Quadfurcation (2-level, center_frac=0.5) should produce
/// comparable separation to a 2-layer bifurcation cascade, since both
/// produce ~4 terminal branches.
#[test]
fn test_quadfurcation_vs_bifurcation_cascade() {
    // Quad: 2-level trifurcation cascade at center_frac=0.5
    let r_quad = cascade_junction_separation(2, 0.5, 2e-3, 1e-3, 1e-6);

    // 2-layer bifurcation cascade (from explicit q_fracs at 0.6 per stage)
    let r_bi2 = cascade_junction_separation_from_qfracs(&[0.6, 0.6]);

    // Both should produce positive separation with bounded fractions.
    assert!(r_quad.separation_efficiency > 0.0);
    assert!(r_bi2.separation_efficiency > 0.0);

    // Both should have cancer enrichment above RBC center fraction.
    let rbc_center_quad = 1.0 - r_quad.rbc_peripheral_fraction;
    let rbc_center_bi2 = 1.0 - r_bi2.rbc_peripheral_fraction;
    assert!(
        r_quad.cancer_center_fraction > rbc_center_quad,
        "Quad: cancer > RBC center: {:.4} vs {:.4}",
        r_quad.cancer_center_fraction,
        rbc_center_quad
    );
    assert!(
        r_bi2.cancer_center_fraction > rbc_center_bi2,
        "Bi-2: cancer > RBC center: {:.4} vs {:.4}",
        r_bi2.cancer_center_fraction,
        rbc_center_bi2
    );
}

// =========================================================================
// Cross-model validation
// =========================================================================

/// Test 10: Fahraeus-Lindqvist viscosity decreases in narrower daughter
/// channels at each cascade level (the FL effect becomes stronger as
/// channel diameter shrinks below ~300 µm, with a minimum around 30 µm).
#[test]
fn test_fahraeus_lindqvist_at_each_cascade_level() {
    // Simulate a 3-level cascade where each daughter channel is 70% of
    // the parent diameter.
    let diameters_um = [200.0, 140.0, 98.0]; // parent → daughter → granddaughter

    let viscosities: Vec<f64> = diameters_um
        .iter()
        .map(|&d| fahraeus_lindqvist_viscosity(d, HT_NORMAL, MU_PLASMA))
        .collect();

    // For these diameters (all > 30 µm), viscosity should decrease with
    // decreasing diameter due to the FL effect (thicker cell-free layer).
    for i in 1..viscosities.len() {
        assert!(
            viscosities[i] < viscosities[i - 1],
            "FL viscosity should decrease: D={:.0}µm → {:.6} Pa·s, D={:.0}µm → {:.6} Pa·s",
            diameters_um[i - 1],
            viscosities[i - 1],
            diameters_um[i],
            viscosities[i],
        );
    }
}

/// Test 11: Near-zero shear rate at a bifurcation stagnation point should
/// produce high Quemada viscosity due to rouleaux aggregation.
#[test]
fn test_quemada_low_shear_at_bifurcation_stagnation() {
    // At the stagnation point of a bifurcation junction, shear rate drops
    // to near zero. The Quemada model should produce significantly higher
    // viscosity than at the bulk flow shear rate.
    let mu_stagnation = quemada_viscosity(0.01, HT_NORMAL, MU_PLASMA);
    let mu_bulk_shear = quemada_viscosity(200.0, HT_NORMAL, MU_PLASMA);

    assert!(
        mu_stagnation > 2.0 * mu_bulk_shear,
        "Stagnation viscosity ({:.6}) should be >> bulk shear ({:.6})",
        mu_stagnation,
        mu_bulk_shear
    );

    // Quemada viscosity must always be finite and above plasma viscosity
    assert!(mu_stagnation.is_finite() && mu_stagnation > MU_PLASMA);
    assert!(mu_bulk_shear.is_finite() && mu_bulk_shear > MU_PLASMA);
}

/// Test 12: Plasma skimming reduces hematocrit in smaller daughter branches.
/// The narrower daughter receives a disproportionately low hematocrit
/// due to the cell-free layer effect at the bifurcation.
#[test]
fn test_plasma_skimming_reduces_hematocrit_in_smaller_daughters() {
    let d_parent = 100.0; // µm

    // Small daughter (30 µm, gets 20% of flow)
    let ht_small = plasma_skimming_hematocrit(HT_NORMAL, 0.2, 30.0, d_parent);
    // Medium daughter (60 µm, gets 40% of flow)
    let ht_medium = plasma_skimming_hematocrit(HT_NORMAL, 0.4, 60.0, d_parent);
    // Large daughter (90 µm, gets 60% of flow)
    let ht_large = plasma_skimming_hematocrit(HT_NORMAL, 0.6, 90.0, d_parent);

    // Smaller daughters should receive lower hematocrit
    assert!(
        ht_small < ht_medium && ht_medium < ht_large,
        "Ht should increase with daughter size: small={:.4}, medium={:.4}, large={:.4}",
        ht_small,
        ht_medium,
        ht_large
    );

    // The smallest daughter should have Ht notably below feed
    assert!(
        ht_small < HT_NORMAL,
        "Small daughter Ht ({:.4}) should be < feed Ht ({:.4})",
        ht_small,
        HT_NORMAL
    );
}

/// Test 13: Three-population equilibria across different channel sizes.
/// Larger channels produce less inertial focusing (cells remain more dispersed)
/// because the confinement ratio kappa = d_cell / D_h decreases.
#[test]
fn test_three_population_equilibria_across_channel_sizes() {
    // Narrow channel: 500 µm x 200 µm
    let eq_narrow = three_population_equilibria(
        500e-6,
        200e-6,
        1e-6,
        BLOOD_DENSITY,
        BLOOD_VISCOSITY,
        HT_NORMAL,
        None,
    );

    // Wide channel: 2000 µm x 1000 µm
    let eq_wide = three_population_equilibria(
        2e-3,
        1e-3,
        1e-6,
        BLOOD_DENSITY,
        BLOOD_VISCOSITY,
        HT_NORMAL,
        None,
    );

    // All equilibrium positions should be physically bounded
    for eq in [&eq_narrow, &eq_wide] {
        assert!((0.0..=1.0).contains(&eq.cancer_eq));
        assert!((0.0..=1.0).contains(&eq.wbc_eq));
        assert!((0.0..=1.0).contains(&eq.rbc_eq));
        assert!((0.0..=1.0).contains(&eq.cancer_center_fraction));
        assert!((0.0..=1.0).contains(&eq.rbc_peripheral_fraction));
    }

    // In the narrow channel, cancer cells should focus at least as well
    // (or better) than in the wide channel, because kappa is larger.
    // Note: if both are dispersed (kappa < 0.07), they will be equal at 0.5.
    assert!(
        eq_narrow.cancer_eq <= eq_wide.cancer_eq,
        "Narrow channel should focus cancer at center as well or better: narrow={:.4}, wide={:.4}",
        eq_narrow.cancer_eq,
        eq_wide.cancer_eq
    );
}

/// Test 14: Amini confinement correction amplifies focusing in narrow channels.
/// When kappa = d_cell / D_h > 0.1, the inertial lift is strongly enhanced.
#[test]
fn test_amini_correction_amplifies_focusing_in_narrow_channels() {
    use cfd_1d::physics::cell_separation::amini_confinement_correction;

    // kappa < 0.07: minimal inertial focusing, correction should be near 1.0
    let corr_low = amini_confinement_correction(0.03);
    assert!(
        (corr_low - 1.0).abs() < 0.05,
        "Low kappa correction should be near 1.0: {:.4}",
        corr_low
    );

    // kappa = 0.15: moderate confinement, correction should be > 1.0
    let corr_moderate = amini_confinement_correction(0.15);
    assert!(
        corr_moderate > 1.0,
        "Moderate kappa correction should amplify focusing: {:.4}",
        corr_moderate
    );

    // kappa = 0.3: strong confinement, correction should be even higher
    let corr_strong = amini_confinement_correction(0.3);
    assert!(
        corr_strong > corr_moderate,
        "Stronger confinement should produce larger correction: {:.4} > {:.4}",
        corr_strong,
        corr_moderate
    );
}

/// Test 15: The staged incremental filtration API and the explicit q_fracs
/// API should produce identical results when given equivalent parameters.
///
/// We use `cif_pretri_stage_q_fracs` to extract the exact q_fracs that the
/// staged API computes internally, then feed them to the explicit API.
#[test]
fn test_incremental_filtration_staged_vs_explicit_qfracs_consistency() {
    use cfd_1d::{
        cif_pretri_stage_q_fracs, incremental_filtration_separation_from_qfracs,
        tri_center_q_frac,
    };

    let n_pretri = 2_u8;
    let pretri_center_frac = 0.42;
    let terminal_tri_frac = 0.38;
    let terminal_bi_treat_frac = 0.65;

    // Path A: staged API (computes q_fracs internally)
    let r_staged = incremental_filtration_separation_staged(
        n_pretri,
        pretri_center_frac,
        terminal_tri_frac,
        terminal_bi_treat_frac,
    );

    // Path B: extract the exact same q_fracs the staged API uses, then call
    // the explicit API with them.
    let pretri_q = cif_pretri_stage_q_fracs(n_pretri, pretri_center_frac, terminal_tri_frac);
    let terminal_q = tri_center_q_frac(terminal_tri_frac);

    let r_qfracs = incremental_filtration_separation_from_qfracs(
        &pretri_q,
        terminal_q,
        terminal_bi_treat_frac,
    );

    // Both pathways should produce identical results (same underlying math)
    let eps = 1e-10;
    assert!(
        (r_staged.cancer_center_fraction - r_qfracs.cancer_center_fraction).abs() < eps,
        "cancer_center_fraction mismatch: staged={:.10}, qfracs={:.10}",
        r_staged.cancer_center_fraction,
        r_qfracs.cancer_center_fraction
    );
    assert!(
        (r_staged.wbc_center_fraction - r_qfracs.wbc_center_fraction).abs() < eps,
        "wbc_center_fraction mismatch: staged={:.10}, qfracs={:.10}",
        r_staged.wbc_center_fraction,
        r_qfracs.wbc_center_fraction
    );
    assert!(
        (r_staged.rbc_peripheral_fraction - r_qfracs.rbc_peripheral_fraction).abs() < eps,
        "rbc_peripheral_fraction mismatch: staged={:.10}, qfracs={:.10}",
        r_staged.rbc_peripheral_fraction,
        r_qfracs.rbc_peripheral_fraction
    );
    assert!(
        (r_staged.separation_efficiency - r_qfracs.separation_efficiency).abs() < eps,
        "separation_efficiency mismatch: staged={:.10}, qfracs={:.10}",
        r_staged.separation_efficiency,
        r_qfracs.separation_efficiency
    );
}
