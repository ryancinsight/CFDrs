//! Comprehensive validation tests for the cell separation physics module.
//!
//! These tests validate the Zweifach-Fung effect and inertial lift focusing
//! models used to predict cell-type-specific routing at microfluidic
//! bifurcation and trifurcation junctions.
//!
//! # Physical background
//!
//! At a branching junction in a microfluidic channel, particles (cells)
//! preferentially enter the arm carrying the highest volumetric flow.
//! This is the **Zweifach-Fung effect** (Fung 1969): large, stiff particles
//! (such as CTCs) are routed more strongly toward high-flow arms than
//! small, deformable particles (RBCs).  The routing probability follows
//! a power law P_center = r_c^beta / (r_c^beta + 2*r_p^beta), where beta
//! is a stiffness exponent (cancer ~1.85, WBC ~1.40, RBC ~1.00).
//!
//! **Inertial lift focusing** (Di Carlo 2009) provides a complementary
//! mechanism in straight rectangular channels at finite Reynolds number:
//! shear-gradient and wall-interaction lift forces drive particles to
//! stable lateral equilibrium positions.  Larger, stiffer cells focus
//! near the channel center; smaller, more deformable cells focus near
//! the walls.  In curved channels, Dean flow secondary circulation
//! (De = Re * sqrt(d/R)) enhances this separation.
//!
//! Together, these effects enable selective enrichment of cancer cells
//! in the center (treatment) arm of a cascade junction tree while
//! diverting healthy blood cells to peripheral bypass channels.

use approx::assert_relative_eq;
use cfd_1d::{
    cascade_junction_separation, cascade_junction_separation_from_qfracs,
    incremental_filtration_separation_from_qfracs, incremental_filtration_separation_staged,
    mixed_cascade_separation, mixed_cascade_separation_kappa_aware, three_population_equilibria,
    tri_asymmetric_q_fracs, tri_center_q_frac, CascadeJunctionResult, CascadeStage,
    IncrementalFiltrationResult,
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Assert that all cell fractions and efficiency are physically bounded.
fn assert_cascade_bounds(r: &CascadeJunctionResult) {
    assert!(
        (0.0..=1.0).contains(&r.cancer_center_fraction),
        "cancer_center_fraction out of bounds: {}",
        r.cancer_center_fraction
    );
    assert!(
        (0.0..=1.0).contains(&r.wbc_center_fraction),
        "wbc_center_fraction out of bounds: {}",
        r.wbc_center_fraction
    );
    assert!(
        (0.0..=1.0).contains(&r.rbc_peripheral_fraction),
        "rbc_peripheral_fraction out of bounds: {}",
        r.rbc_peripheral_fraction
    );
    assert!(
        r.separation_efficiency >= 0.0,
        "separation_efficiency negative: {}",
        r.separation_efficiency
    );
}

fn assert_filtration_bounds(r: &IncrementalFiltrationResult) {
    assert!(
        (0.0..=1.0).contains(&r.cancer_center_fraction),
        "cancer_center_fraction out of bounds: {}",
        r.cancer_center_fraction
    );
    assert!(
        (0.0..=1.0).contains(&r.wbc_center_fraction),
        "wbc_center_fraction out of bounds: {}",
        r.wbc_center_fraction
    );
    assert!(
        (0.0..=1.0).contains(&r.rbc_center_fraction),
        "rbc_center_fraction out of bounds: {}",
        r.rbc_center_fraction
    );
    assert!(
        (0.0..=1.0).contains(&r.rbc_peripheral_fraction),
        "rbc_peripheral_fraction out of bounds: {}",
        r.rbc_peripheral_fraction
    );
    assert!(
        r.separation_efficiency >= 0.0,
        "separation_efficiency negative: {}",
        r.separation_efficiency
    );
}

// ---------------------------------------------------------------------------
// Test 1: Cascade junction -- equal trifurcation
// ---------------------------------------------------------------------------

#[test]
fn test_cascade_junction_equal_trifurcation() {
    let r = cascade_junction_separation(
        1,     // n_levels
        0.333, // center_frac (equal widths)
        2e-3,  // channel_width
        1e-3,  // channel_height
        1e-6,  // flow_rate
    );

    assert_cascade_bounds(&r);

    // With equal trifurcation (center_frac ~ 1/3), each arm gets roughly
    // equal flow, so all cell types route ~1/3 to the center arm.
    // Cancer cells (beta=1.85) should still show slight center bias even
    // at near-equal flow due to the power-law exponent.
    assert!(
        r.cancer_center_fraction > 0.0 && r.cancer_center_fraction <= 1.0,
        "cancer fraction should be positive"
    );
}

// ---------------------------------------------------------------------------
// Test 2: Cascade junction -- biased center flow enriches cancer
// ---------------------------------------------------------------------------

#[test]
fn test_cascade_junction_biased_center_increases_cancer_enrichment() {
    let r_equal = cascade_junction_separation(1, 0.333, 2e-3, 1e-3, 1e-6);
    let r_biased = cascade_junction_separation(1, 0.5, 2e-3, 1e-3, 1e-6);

    assert_cascade_bounds(&r_biased);

    // Wider center arm (0.5 vs 0.333) carries more flow, which increases
    // the Zweifach-Fung routing bias for stiff cancer cells.
    assert!(
        r_biased.cancer_center_fraction > r_equal.cancer_center_fraction,
        "biased center (0.5) should yield higher cancer enrichment than equal (0.333): {} vs {}",
        r_biased.cancer_center_fraction,
        r_equal.cancer_center_fraction
    );
}

// ---------------------------------------------------------------------------
// Test 3: Multi-level cascade increases separation efficiency
// ---------------------------------------------------------------------------

#[test]
fn test_multilevel_cascade_increases_separation() {
    let r1 = cascade_junction_separation(1, 0.4, 2e-3, 1e-3, 1e-6);
    let r2 = cascade_junction_separation(2, 0.4, 2e-3, 1e-3, 1e-6);

    assert_cascade_bounds(&r1);
    assert_cascade_bounds(&r2);

    // Deeper cascades apply the Zweifach-Fung routing multiplicatively,
    // widening the gap between cancer (high beta) and RBC (low beta)
    // center fractions. This should increase separation efficiency.
    assert!(
        r2.separation_efficiency > r1.separation_efficiency,
        "2-level cascade should have higher separation efficiency than 1-level: {} vs {}",
        r2.separation_efficiency,
        r1.separation_efficiency
    );
}

// ---------------------------------------------------------------------------
// Test 4: Cascade from explicit q_fracs
// ---------------------------------------------------------------------------

#[test]
fn test_cascade_from_explicit_qfracs() {
    // Two-stage equal trifurcation with explicit center-flow fractions.
    let r = cascade_junction_separation_from_qfracs(&[0.33, 0.33]);

    assert_cascade_bounds(&r);

    // With equal-ish flow fractions, cancer cells should still preferentially
    // route to the center due to beta > 1.
    assert!(
        r.cancer_center_fraction > 0.0,
        "cancer fraction should be positive with equal q_fracs"
    );
}

// ---------------------------------------------------------------------------
// Test 5: Incremental filtration -- staged
// ---------------------------------------------------------------------------

#[test]
fn test_incremental_filtration_staged() {
    let r = incremental_filtration_separation_staged(
        1,    // n_pretri
        0.33, // pretri_center_frac
        0.33, // terminal_tri_frac
        0.68, // terminal_bi_treat_frac
    );

    assert_filtration_bounds(&r);

    // The staged model adds a terminal bifurcation after trifurcation levels,
    // providing additional separation.
    assert!(
        r.separation_efficiency >= 0.0,
        "staged filtration should produce non-negative separation"
    );
}

// ---------------------------------------------------------------------------
// Test 6: Deeper pretri increases cancer enrichment
// ---------------------------------------------------------------------------

#[test]
fn test_deeper_pretri_increases_enrichment() {
    let r1 = incremental_filtration_separation_staged(1, 0.40, 0.33, 0.68);
    let r2 = incremental_filtration_separation_staged(2, 0.40, 0.33, 0.68);

    assert_filtration_bounds(&r1);
    assert_filtration_bounds(&r2);

    // More pre-trifurcation stages provide additional Zweifach-Fung routing
    // opportunities, increasing cancer enrichment in the treatment arm.
    assert!(
        r2.separation_efficiency >= r1.separation_efficiency,
        "2 pretri levels should have >= separation efficiency vs 1: {} vs {}",
        r2.separation_efficiency,
        r1.separation_efficiency
    );
}

// ---------------------------------------------------------------------------
// Test 7: Mixed cascade separation -- trifurcation then bifurcation
// ---------------------------------------------------------------------------

#[test]
fn test_mixed_cascade_separation() {
    // Stage 1: trifurcation with q_frac=0.45 (biased center, is_trifurcation=true)
    // Stage 2: bifurcation with q_frac=0.6 (treatment-biased, is_trifurcation=false)
    let r = mixed_cascade_separation(&[(0.45, true), (0.6, false)]);

    assert_cascade_bounds(&r);

    // With center/treatment-biased flow fractions (> 1/3 for tri, > 0.5 for bi),
    // cancer cells (beta=1.85) should route more strongly to the treatment arm
    // than RBCs (beta=1.0) due to the Zweifach-Fung power-law.
    let rbc_center = 1.0 - r.rbc_peripheral_fraction;
    assert!(
        r.cancer_center_fraction > rbc_center,
        "cancer center fraction should exceed RBC center fraction: cancer={}, rbc_center={}",
        r.cancer_center_fraction,
        rbc_center
    );
}

// ---------------------------------------------------------------------------
// Test 8: Three-population equilibria -- straight channel
// ---------------------------------------------------------------------------

#[test]
fn test_three_population_equilibria_straight() {
    let eq = three_population_equilibria(
        2e-3,   // width
        1e-3,   // height
        1e-6,   // flow_rate
        1060.0, // blood density
        3.5e-3, // viscosity
        0.45,   // hematocrit
        None,   // no bend (straight channel)
    );

    // All equilibrium positions should be in [0, 1] (center=0, wall=1).
    assert!(
        (0.0..=1.0).contains(&eq.cancer_eq),
        "cancer_eq out of bounds: {}",
        eq.cancer_eq
    );
    assert!(
        (0.0..=1.0).contains(&eq.wbc_eq),
        "wbc_eq out of bounds: {}",
        eq.wbc_eq
    );
    assert!(
        (0.0..=1.0).contains(&eq.rbc_eq),
        "rbc_eq out of bounds: {}",
        eq.rbc_eq
    );

    // Center fractions and peripheral fractions should be bounded.
    assert!(
        (0.0..=1.0).contains(&eq.cancer_center_fraction),
        "cancer_center_fraction out of bounds: {}",
        eq.cancer_center_fraction
    );
    assert!(
        (0.0..=1.0).contains(&eq.rbc_peripheral_fraction),
        "rbc_peripheral_fraction out of bounds: {}",
        eq.rbc_peripheral_fraction
    );

    // In a straight channel, large stiff cancer cells (17.5 um) should focus
    // closer to the center (lower x_tilde) than small deformable RBCs (7 um),
    // provided the confinement ratio is sufficient for inertial focusing.
    // If cells don't focus (kappa < 0.07), both return DISPERSED (0.5).
    // We allow equality for the dispersed case.
    assert!(
        eq.cancer_eq <= eq.rbc_eq,
        "cancer should focus at center or equal to RBC: cancer_eq={}, rbc_eq={}",
        eq.cancer_eq,
        eq.rbc_eq
    );
}

// ---------------------------------------------------------------------------
// Test 9: Three-population equilibria -- serpentine (Dean flow)
// ---------------------------------------------------------------------------

#[test]
fn test_three_population_equilibria_serpentine() {
    let eq_straight = three_population_equilibria(
        2e-3,   // width
        1e-3,   // height
        1e-6,   // flow_rate
        1060.0, // blood density
        3.5e-3, // viscosity
        0.45,   // hematocrit
        None,   // straight
    );

    let eq_serpentine = three_population_equilibria(
        2e-3,           // width
        1e-3,           // height
        1e-6,           // flow_rate
        1060.0,         // blood density
        3.5e-3,         // viscosity
        0.45,           // hematocrit
        Some(2e-3),     // bend radius = 2mm (tight serpentine)
    );

    // All equilibrium positions should remain bounded.
    assert!(
        (0.0..=1.0).contains(&eq_serpentine.cancer_eq),
        "serpentine cancer_eq out of bounds: {}",
        eq_serpentine.cancer_eq
    );
    assert!(
        (0.0..=1.0).contains(&eq_serpentine.rbc_eq),
        "serpentine rbc_eq out of bounds: {}",
        eq_serpentine.rbc_eq
    );

    // Dean flow from channel curvature should modify equilibrium positions.
    // The secondary circulation creates an additional drag that shifts focusing
    // positions relative to the straight-channel case.
    // We check that at least one population's equilibrium is different.
    let cancer_changed = (eq_serpentine.cancer_eq - eq_straight.cancer_eq).abs() > 1e-12;
    let wbc_changed = (eq_serpentine.wbc_eq - eq_straight.wbc_eq).abs() > 1e-12;
    let rbc_changed = (eq_serpentine.rbc_eq - eq_straight.rbc_eq).abs() > 1e-12;

    assert!(
        cancer_changed || wbc_changed || rbc_changed,
        "Dean flow should modify at least one population's equilibrium position; \
         straight=({:.4}, {:.4}, {:.4}), serpentine=({:.4}, {:.4}, {:.4})",
        eq_straight.cancer_eq,
        eq_straight.wbc_eq,
        eq_straight.rbc_eq,
        eq_serpentine.cancer_eq,
        eq_serpentine.wbc_eq,
        eq_serpentine.rbc_eq,
    );
}

// ---------------------------------------------------------------------------
// Test 10: Tri center q_frac -- equal widths
// ---------------------------------------------------------------------------

#[test]
fn test_tri_center_q_frac_equal_widths() {
    // With center_frac = 1/3 (equal arm widths), the center arm should
    // carry exactly 1/3 of the total flow.
    let q = tri_center_q_frac(1.0 / 3.0);

    assert_relative_eq!(q, 1.0 / 3.0, epsilon = 1e-6);
}

// ---------------------------------------------------------------------------
// Test 11: Tri center q_frac -- wide center
// ---------------------------------------------------------------------------

#[test]
fn test_tri_center_q_frac_wide_center() {
    // With center_frac = 0.5 (wider center arm), the center arm should carry
    // more than half the total flow due to the cubic width dependence of
    // Hagen-Poiseuille resistance (R ~ 1/w^3).
    let q = tri_center_q_frac(0.5);

    assert!(
        q > 0.5,
        "wider center arm (frac=0.5) should carry > 50% of flow due to cubic width law: q={}",
        q
    );

    // Verify the cubic law amplification: center is 0.5, each peripheral
    // is 0.25. Flow ~ w^3, so q_center = 0.5^3 / (0.5^3 + 2*0.25^3)
    //   = 0.125 / (0.125 + 0.03125) = 0.125 / 0.15625 = 0.8
    let expected = 0.125 / (0.125 + 2.0 * 0.015625);
    assert_relative_eq!(q, expected, epsilon = 1e-6);
}

// ---------------------------------------------------------------------------
// Test 12: Tri asymmetric q_fracs -- sum to 1
// ---------------------------------------------------------------------------

#[test]
fn test_tri_asymmetric_q_fracs_sum_to_one() {
    // Asymmetric trifurcation with widths proportional to 1:2:3.
    // Total width = 6 parts, so center_frac = 2/6, left_periph_frac = 1/6,
    // right_periph = 3/6.
    let fracs = tri_asymmetric_q_fracs(
        2.0 / 6.0, // center_frac
        1.0 / 6.0, // left_periph_frac
        6e-3,      // parent_width_m (6mm so branches are 1mm, 2mm, 3mm)
        1e-3,      // channel_height_m
    );

    assert_eq!(fracs.len(), 3, "should return exactly 3 flow fractions");

    let sum: f64 = fracs.iter().sum();
    assert_relative_eq!(sum, 1.0, epsilon = 1e-6);

    // All fractions should be positive.
    for (i, &f) in fracs.iter().enumerate() {
        assert!(f > 0.0, "fraction[{}] should be positive: {}", i, f);
    }

    // The widest arm (right, 3/6) should carry the most flow.
    assert!(
        fracs[2] > fracs[0] && fracs[2] > fracs[1],
        "widest arm should carry most flow: {:?}",
        fracs
    );
}

// ---------------------------------------------------------------------------
// Test 13: Kappa-aware cascade -- confinement-ratio amplification
// ---------------------------------------------------------------------------

#[test]
fn test_kappa_aware_cascade_selective_enrichment() {
    // Build a single stage with a treatment arm narrow enough that the CTC
    // confinement ratio kappa = d_cancer / Dh is significant (> 0.07),
    // amplifying the stiffness exponent for cancer cells.
    //
    // Treatment arm: 200um x 100um -> Dh = 2*w*h/(w+h) = 2*200e-6*100e-6 / 300e-6 = 133um
    // kappa_cancer = 17.5e-6 / 133e-6 = 0.13 (> 0.07, significant focusing)
    // kappa_rbc = 7e-6 / 133e-6 = 0.053 (< 0.07, minimal focusing)
    let treatment_width = 200e-6;
    let treatment_height = 100e-6;
    let treatment_dh = 2.0 * treatment_width * treatment_height
        / (treatment_width + treatment_height);

    let stage = CascadeStage {
        arm_q_fracs: [0.5, 0.25, 0.25, 0.0, 0.0],
        n_arms: 3,
        treatment_dh_m: treatment_dh,
        parent_v_in_m_s: 0.02,
        peripheral_recoveries: [None, None, None, None],
        n_recoveries: 0,
    };

    let r = mixed_cascade_separation_kappa_aware(&[stage]);

    assert_cascade_bounds(&r);

    // The kappa amplification should increase beta for cancer cells more
    // than for RBCs (whose beta=1.0 is not amplified), leading to stronger
    // preferential routing of cancer cells to the treatment arm.
    let rbc_center_fraction = 1.0 - r.rbc_peripheral_fraction;
    assert!(
        r.cancer_center_fraction > rbc_center_fraction,
        "kappa-aware routing should enrich cancer over RBC in treatment arm: \
         cancer={:.4}, rbc_center={:.4}",
        r.cancer_center_fraction,
        rbc_center_fraction
    );

    // Verify that separation efficiency is non-trivial with confinement effects.
    assert!(
        r.separation_efficiency > 0.0,
        "kappa-aware cascade should produce positive separation efficiency: {}",
        r.separation_efficiency
    );
}

// ---------------------------------------------------------------------------
// Additional physics validation tests
// ---------------------------------------------------------------------------

/// Verify that the incremental filtration from explicit q_fracs produces
/// consistent results with the staged version.
#[test]
fn test_incremental_filtration_from_qfracs_consistency() {
    let r_staged = incremental_filtration_separation_staged(1, 0.40, 0.35, 0.65);

    // Compute the equivalent q_fracs manually:
    // pretri: tri_center_q_frac(0.40) for 1 stage
    // terminal tri: tri_center_q_frac(0.35)
    let pretri_q = tri_center_q_frac(0.40);
    let terminal_tri_q = tri_center_q_frac(0.35);

    let r_qfracs =
        incremental_filtration_separation_from_qfracs(&[pretri_q], terminal_tri_q, 0.65);

    assert_filtration_bounds(&r_qfracs);

    // Both pathways should produce identical results.
    assert_relative_eq!(
        r_staged.cancer_center_fraction,
        r_qfracs.cancer_center_fraction,
        epsilon = 1e-10
    );
    assert_relative_eq!(
        r_staged.rbc_peripheral_fraction,
        r_qfracs.rbc_peripheral_fraction,
        epsilon = 1e-10
    );
    assert_relative_eq!(
        r_staged.separation_efficiency,
        r_qfracs.separation_efficiency,
        epsilon = 1e-10
    );
}

/// Verify that the mixed cascade with a single trifurcation is equivalent
/// to the cascade junction with n_levels=1.
#[test]
fn test_mixed_cascade_matches_single_trifurcation() {
    let q = tri_center_q_frac(0.45);

    let r_cascade = cascade_junction_separation(1, 0.45, 2e-3, 1e-3, 1e-6);
    let r_mixed = mixed_cascade_separation(&[(q, true)]);

    assert_relative_eq!(
        r_cascade.cancer_center_fraction,
        r_mixed.cancer_center_fraction,
        epsilon = 1e-10
    );
    assert_relative_eq!(
        r_cascade.rbc_peripheral_fraction,
        r_mixed.rbc_peripheral_fraction,
        epsilon = 1e-10
    );
    assert_relative_eq!(
        r_cascade.separation_efficiency,
        r_mixed.separation_efficiency,
        epsilon = 1e-10
    );
}

/// Verify that the Zweifach-Fung routing preserves the ordering:
/// cancer > WBC > RBC center fraction for asymmetric trifurcations.
#[test]
fn test_zweifach_fung_cell_type_ordering() {
    let r = cascade_junction_separation(1, 0.5, 2e-3, 1e-3, 1e-6);

    // Cancer (beta=1.85) > WBC (beta=1.40) > RBC (beta=1.00) at the center arm
    // when center_frac > 1/3 (biased flow toward center).
    assert!(
        r.cancer_center_fraction >= r.wbc_center_fraction,
        "cancer should route to center >= WBC: cancer={:.4}, wbc={:.4}",
        r.cancer_center_fraction,
        r.wbc_center_fraction
    );

    let rbc_center = 1.0 - r.rbc_peripheral_fraction;
    assert!(
        r.wbc_center_fraction >= rbc_center,
        "WBC should route to center >= RBC: wbc={:.4}, rbc_center={:.4}",
        r.wbc_center_fraction,
        rbc_center
    );
}
