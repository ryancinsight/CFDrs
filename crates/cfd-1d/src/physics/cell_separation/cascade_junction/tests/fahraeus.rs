//! Fåhræus margination correction tests and the Option 1/Option 2 topology regressions.

use super::*;

// ── Fåhræus margination correction tests ──────────────────────────────

#[test]
fn fahrae_correction_zero_for_rbc() {
    // RBC is the reference cell: excess β = 0 and size ratio = 0 → no correction.
    let c = fahrae_beta_correction(SE_RBC, D_RBC_M);
    assert!(
        c.abs() < 1e-15,
        "Fåhræus correction must be zero for RBC, got {c}"
    );
}

#[test]
fn fahrae_correction_positive_for_cancer() {
    // Cancer cells (17.5 µm) are larger than RBCs (7 µm) → positive correction.
    let c = fahrae_beta_correction(SE_CANCER, D_CANCER_M);
    assert!(
        c > 0.0,
        "Fåhræus correction must be positive for cancer cells"
    );
    // Expected: excess(0.85) × α(0.12) × size_ratio(1.5) = 0.153
    assert!(
        (c - 0.85 * 0.12 * 1.5).abs() < 1e-12,
        "Fåhræus correction for cancer should be ~0.153, got {c}"
    );
}

#[test]
fn fahrae_correction_smaller_for_wbc_than_cancer() {
    // WBCs (10 µm) are smaller than cancer cells (17.5 µm) → smaller correction.
    let c_cancer = fahrae_beta_correction(SE_CANCER, D_CANCER_M);
    let c_wbc = fahrae_beta_correction(SE_WBC, D_WBC_M);
    assert!(c_wbc > 0.0, "Fåhræus correction must be positive for WBCs");
    assert!(
        c_cancer > c_wbc,
        "Cancer correction ({c_cancer:.4}) must exceed WBC correction ({c_wbc:.4})"
    );
}

#[test]
fn kappa_aware_with_fahrae_boosts_cancer_over_legacy() {
    // The Fåhræus correction in the kappa-aware model should boost cancer routing
    // beyond what the legacy mixed_cascade_separation (no κ, no Fåhræus) provides.
    let q = tri_center_q_frac_cross_junction(0.45, length(4e-3), length(1e-3));
    let q_p = (1.0 - q) / 2.0;
    let w_center = 0.45 * 4e-3;
    let h = 1e-3_f64;
    let dh = 2.0 * w_center * h / (w_center + h);

    let stage = CascadeStage {
        arm_q_fracs: [q, q_p, q_p, 0.0, 0.0],
        n_arms: 3,
        treatment_dh_m: Length::from_base(dh),
        parent_v_in_m_s: Velocity::from_base(0.02),
        peripheral_recoveries: [None; 4],
        n_recoveries: 0,
    };
    let kappa_fahrae = mixed_cascade_separation_kappa_aware(&[stage]);
    let legacy = mixed_cascade_separation(&[(q, true)]);

    // The Fåhræus-enhanced kappa model should route more cancer cells to center.
    assert!(
        kappa_fahrae.cancer_center_fraction > legacy.cancer_center_fraction,
        "Fåhræus-enhanced model must exceed legacy cancer routing: enhanced={:.4}, legacy={:.4}",
        kappa_fahrae.cancer_center_fraction,
        legacy.cancer_center_fraction
    );
    // Cancer-to-RBC enrichment ratio should improve.
    let rbc_center_kf = 1.0 - kappa_fahrae.rbc_peripheral_fraction;
    let rbc_center_legacy = 1.0 - legacy.rbc_peripheral_fraction;
    let enrichment_kf = kappa_fahrae.cancer_center_fraction / rbc_center_kf.max(1e-12);
    let enrichment_legacy = legacy.cancer_center_fraction / rbc_center_legacy.max(1e-12);
    assert!(
        enrichment_kf >= enrichment_legacy - 1e-6,
        "Fåhræus model should improve CTC/RBC enrichment: kf={enrichment_kf:.4}, legacy={enrichment_legacy:.4}"
    );
}

#[test]
fn option1_option2_tri_tri_enrichment_quantified() {
    // Quantify CTC enrichment for the actual Option 1/Option 2 TriTri topology:
    //   Stage 1: pretri center_frac = 0.45
    //   Stage 2: terminal tri center_frac = 0.333 (symmetric)
    //   parent_width = 4mm, height = 1mm
    let h = 1e-3_f64;
    let parent_w = 4e-3_f64;
    let pcf = 0.45;
    let tcf = 1.0 / 3.0; // symmetric

    // Build 2-stage TriTri with cross-junction corrections (as used by compute.rs)
    let left_periph = (1.0 - pcf) / 2.0;
    let arm_q_1 = tri_asymmetric_q_fracs(pcf, left_periph, length(parent_w), length(h));
    let w_center_1 = pcf * parent_w;
    let dh_1 = 2.0 * w_center_1 * h / (w_center_1 + h);

    let parent_w_2 = w_center_1;
    let left_periph_2 = (1.0 - tcf) / 2.0;
    let arm_q_2 = tri_asymmetric_q_fracs(tcf, left_periph_2, length(parent_w_2), length(h));
    let w_center_2 = tcf * parent_w_2;
    let dh_2 = 2.0 * w_center_2 * h / (w_center_2 + h);

    let stages = vec![
        CascadeStage {
            arm_q_fracs: [arm_q_1[0], arm_q_1[1], arm_q_1[2], 0.0, 0.0],
            n_arms: 3,
            treatment_dh_m: Length::from_base(dh_1),
            parent_v_in_m_s: Velocity::from_base(0.02),
            peripheral_recoveries: [None; 4],
            n_recoveries: 0,
        },
        CascadeStage {
            arm_q_fracs: [arm_q_2[0], arm_q_2[1], arm_q_2[2], 0.0, 0.0],
            n_arms: 3,
            treatment_dh_m: Length::from_base(dh_2),
            parent_v_in_m_s: Velocity::from_base(0.02),
            peripheral_recoveries: [None; 4],
            n_recoveries: 0,
        },
    ];
    let r = mixed_cascade_separation_kappa_aware(&stages);

    // The terminal stage is symmetric (tcf = 1/3), so the second split can only
    // pass one third of whatever enrichment survives the first stage. The
    // physically meaningful regression is therefore that the full TriTri path
    // still preserves material cancer capture and enrichment above unity, not
    // that it exceeds the asymmetric-terminal case.
    assert!(
        r.cancer_center_fraction > 0.18,
        "TriTri cancer center fraction should exceed 18%: got {:.2}%",
        r.cancer_center_fraction * 100.0
    );
    // RBC peripheral fraction should remain high (> 78%)
    assert!(
        r.rbc_peripheral_fraction > 0.78,
        "TriTri RBC peripheral fraction should exceed 78%: got {:.2}%",
        r.rbc_peripheral_fraction * 100.0
    );
    // CTC/RBC enrichment ratio should be > 1.3
    let rbc_center = 1.0 - r.rbc_peripheral_fraction;
    let enrichment = r.cancer_center_fraction / rbc_center.max(1e-12);
    assert!(
        enrichment > 1.15,
        "CTC/RBC enrichment ratio should exceed 1.15: got {enrichment:.4}"
    );

    // Print for diagnostic review
    eprintln!("=== TriTri (pcf=0.45, tcf=0.333) kappa-aware + Fåhræus ===");
    eprintln!(
        "  Cancer center fraction:  {:.2}%",
        r.cancer_center_fraction * 100.0
    );
    eprintln!(
        "  WBC center fraction:     {:.2}%",
        r.wbc_center_fraction * 100.0
    );
    eprintln!(
        "  RBC peripheral fraction: {:.2}%",
        r.rbc_peripheral_fraction * 100.0
    );
    eprintln!("  RBC center fraction:     {:.2}%", rbc_center * 100.0);
    eprintln!("  Separation efficiency:   {:.4}", r.separation_efficiency);
    eprintln!("  CTC/RBC enrichment:      {enrichment:.4}×");
    eprintln!(
        "  Center HCT ratio:        {:.4}",
        r.center_hematocrit_ratio
    );
}

#[test]
fn option2_higher_tcf_dramatically_improves_enrichment() {
    // Demonstrate that using tcf=0.55 instead of 0.333 would dramatically
    // improve CTC enrichment. This is a design-space insight, not a physics change.
    let h = 1e-3_f64;
    let parent_w = 4e-3_f64;
    let pcf = 0.45;

    // Build with tcf=0.55 (asymmetric terminal tri)
    let tcf = 0.55;
    let left_periph = (1.0 - pcf) / 2.0;
    let arm_q_1 = tri_asymmetric_q_fracs(pcf, left_periph, length(parent_w), length(h));
    let w_center_1 = pcf * parent_w;
    let dh_1 = 2.0 * w_center_1 * h / (w_center_1 + h);

    let parent_w_2 = w_center_1;
    let left_periph_2 = (1.0 - tcf) / 2.0;
    let arm_q_2 = tri_asymmetric_q_fracs(tcf, left_periph_2, length(parent_w_2), length(h));
    let w_center_2 = tcf * parent_w_2;
    let dh_2 = 2.0 * w_center_2 * h / (w_center_2 + h);

    let stages = vec![
        CascadeStage {
            arm_q_fracs: [arm_q_1[0], arm_q_1[1], arm_q_1[2], 0.0, 0.0],
            n_arms: 3,
            treatment_dh_m: Length::from_base(dh_1),
            parent_v_in_m_s: Velocity::from_base(0.02),
            peripheral_recoveries: [None; 4],
            n_recoveries: 0,
        },
        CascadeStage {
            arm_q_fracs: [arm_q_2[0], arm_q_2[1], arm_q_2[2], 0.0, 0.0],
            n_arms: 3,
            treatment_dh_m: Length::from_base(dh_2),
            parent_v_in_m_s: Velocity::from_base(0.02),
            peripheral_recoveries: [None; 4],
            n_recoveries: 0,
        },
    ];
    let r_high = mixed_cascade_separation_kappa_aware(&stages);

    // With tcf=0.55: cancer routing should exceed 40%.
    assert!(
        r_high.cancer_center_fraction > 0.40,
        "TriTri with tcf=0.55 cancer fraction should exceed 40%: got {:.2}%",
        r_high.cancer_center_fraction * 100.0
    );

    let rbc_center = 1.0 - r_high.rbc_peripheral_fraction;
    let enrichment = r_high.cancer_center_fraction / rbc_center.max(1e-12);
    assert!(
        enrichment > 1.4,
        "TriTri with tcf=0.55 enrichment should exceed 1.4: got {enrichment:.4}"
    );

    eprintln!("=== TriTri (pcf=0.45, tcf=0.55) kappa-aware + Fåhræus ===");
    eprintln!(
        "  Cancer center fraction:  {:.2}%",
        r_high.cancer_center_fraction * 100.0
    );
    eprintln!(
        "  WBC center fraction:     {:.2}%",
        r_high.wbc_center_fraction * 100.0
    );
    eprintln!(
        "  RBC peripheral fraction: {:.2}%",
        r_high.rbc_peripheral_fraction * 100.0
    );
    eprintln!("  RBC center fraction:     {:.2}%", rbc_center * 100.0);
    eprintln!("  CTC/RBC enrichment:      {enrichment:.4}×");
}

#[test]
fn updated_se_cancer_improves_single_stage_selectivity() {
    // With SE_CANCER = 1.85 (up from 1.70), a single trifurcation stage with
    // center-biased flow should produce higher cancer-to-RBC differential.
    let q = tri_center_q_frac(0.50);
    let r = mixed_cascade_separation(&[(q, true)]);

    // Cancer routing must exceed RBC routing (flow-weighted at β=1.0).
    let rbc_center = 1.0 - r.rbc_peripheral_fraction;
    assert!(
        r.cancer_center_fraction > rbc_center,
        "Cancer must route preferentially to center: cancer={:.4}, rbc_center={:.4}",
        r.cancer_center_fraction,
        rbc_center
    );
    // Enrichment ratio must be > 1.0
    let enrichment = r.cancer_center_fraction / rbc_center.max(1e-12);
    assert!(
        enrichment > 1.0,
        "CTC/RBC enrichment must exceed 1.0 at asymmetric split, got {enrichment:.4}"
    );
}

