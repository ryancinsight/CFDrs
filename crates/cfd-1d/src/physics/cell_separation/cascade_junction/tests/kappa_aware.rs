//! kappa-aware routing tests: beta amplification and asymmetric arm selection.

use super::*;

// ── kappa-aware tests ─────────────────────────────────────────────────────

#[test]
fn beta_kappa_adjusted_no_change_at_zero_kappa() {
    // At κ = 0 (cell infinitely smaller than channel), β_eff = β_base.
    assert!((beta_kappa_adjusted(SE_CANCER, 0.0, 0.02, false) - SE_CANCER).abs() < 1e-12);
    assert!((beta_kappa_adjusted(SE_WBC, 0.0, 0.02, false) - SE_WBC).abs() < 1e-12);
    assert!((beta_kappa_adjusted(SE_RBC, 0.0, 0.02, true) - SE_RBC).abs() < 1e-12);
}

#[test]
fn beta_kappa_adjusted_rbc_never_amplified() {
    // RBC is fully deformable (β_base = 1.0, excess = 0) — no amplification regardless of κ.
    for kappa in [0.0, 0.03, 0.07, 0.20] {
        // Evaluated below inversion velocity
        let b = beta_kappa_adjusted(SE_RBC, kappa, 0.02, true);
        assert!(
            (b - 1.0).abs() < 1e-12,
            "RBC β should stay 1.0 at κ={kappa}, got {b}"
        );
    }
}

#[test]
fn beta_kappa_adjusted_cancer_increases_with_kappa() {
    let b0 = beta_kappa_adjusted(SE_CANCER, 0.0, 0.02, false);
    let b1 = beta_kappa_adjusted(SE_CANCER, KAPPA_REF, 0.02, false);
    let b2 = beta_kappa_adjusted(SE_CANCER, KAPPA_REF * 2.0, 0.02, false);
    assert!(b1 > b0, "β should increase with κ for stiff cancer cells");
    assert!(b2 >= b1, "β should not decrease as κ grows further");
    assert!(b2 <= 3.0, "β must be capped at 3.0");
}

#[test]
fn p_arm_general_symmetric_trifurcation_matches_p_center() {
    // Symmetric trifurcation: arm_q_fracs = [q, q_p, q_p]
    let q_c = 0.55_f64;
    let q_p = (1.0 - q_c) / 2.0;
    let p_gen = p_arm_general(&[q_c, q_p, q_p], 0, SE_CANCER);
    let p_old = p_center(q_c, SE_CANCER);
    assert!(
        (p_gen - p_old).abs() < 1e-10,
        "generalised p_arm must match p_center for symmetric tri: gen={p_gen}, old={p_old}"
    );
}

#[test]
fn p_arm_general_two_arm_matches_p_treat_bifurcation() {
    let q_t = 0.68_f64;
    let q_b = 1.0 - q_t;
    let p_gen = p_arm_general(&[q_t, q_b], 0, SE_CANCER);
    let p_old = p_treat_bifurcation(q_t, SE_CANCER);
    assert!(
        (p_gen - p_old).abs() < 1e-10,
        "generalised p_arm must match p_treat_bifurcation: gen={p_gen}, old={p_old}"
    );
}

#[test]
fn tri_asymmetric_q_fracs_sum_to_one() {
    let [qc, ql, qr] = tri_asymmetric_q_fracs(0.45, 0.30, length(4e-3), length(1e-3));
    assert!(
        (qc + ql + qr - 1.0).abs() < 1e-10,
        "asymmetric arm fracs must sum to 1: {qc}+{ql}+{qr}={:.6}",
        qc + ql + qr
    );
}

#[test]
fn tri_asymmetric_wider_periph_gets_more_flow() {
    // Left peripheral wider than right → left gets more flow.
    let [_qc, ql, qr] = tri_asymmetric_q_fracs(0.40, 0.40, length(4e-3), length(1e-3));
    // left_frac=0.40, right_frac=0.20 → left should carry more
    assert!(
        ql > qr,
        "wider left peripheral should carry more flow: ql={ql}, qr={qr}"
    );
}

#[test]
fn checked_tri_asymmetric_q_fracs_reject_overfull_width_budget() {
    let err = checked_tri_asymmetric_q_fracs(0.70, 0.35, length(4e-3), length(1e-3))
        .expect_err(
            "checked asymmetric trifurcation flow fractions must reject overfull width budgets",
        );
    assert!(err.to_string().contains("positive right-arm width"));
}

#[test]
fn tri_asymmetric_symmetric_matches_cross_junction() {
    // With equal peripherals, tri_asymmetric_q_fracs should match tri_center_q_frac_cross_junction.
    let center_frac = 0.45_f64;
    let left_frac = (1.0 - center_frac) / 2.0; // symmetric
    let [qc_asym, _, _] =
        tri_asymmetric_q_fracs(center_frac, left_frac, length(4e-3), length(1e-3));
    let qc_sym = tri_center_q_frac_cross_junction(center_frac, length(4e-3), length(1e-3));
    assert!(
        (qc_asym - qc_sym).abs() < 1e-8,
        "symmetric asymmetric must match cross-junction: asym={qc_asym}, sym={qc_sym}"
    );
}

#[test]
fn kappa_aware_higher_beta_than_legacy_for_stiff_cells_in_narrow_channel() {
    // In a narrow channel (Dh = 0.5mm), cancer cell κ ≈ 0.035 (above zero).
    // The kappa-aware model should produce HIGHER cancer_center_fraction than
    // the legacy model using constant β.
    let q = tri_center_q_frac_cross_junction(0.45, length(0.8e-3), length(1e-3));
    let q_p = (1.0 - q) / 2.0;
    let w_center = 0.45 * 0.8e-3;
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

    let kappa_result = mixed_cascade_separation_kappa_aware(&[stage]);
    let legacy_result = mixed_cascade_separation(&[(q, true)]);

    // κ ≈ 17.5µm / Dh → β_cancer > SE_CANCER = 1.85 → stronger cancer routing
    assert!(
        kappa_result.cancer_center_fraction >= legacy_result.cancer_center_fraction - 1e-10,
        "kappa-aware model must not reduce cancer routing vs legacy: kappa={:.4}, legacy={:.4}",
        kappa_result.cancer_center_fraction,
        legacy_result.cancer_center_fraction
    );
    // RBC routing must be unchanged (deformable, β stays 1.0)
    assert!(
        (kappa_result.rbc_peripheral_fraction - legacy_result.rbc_peripheral_fraction).abs()
            < 1e-9,
        "RBC routing must be unchanged by κ correction"
    );
}

#[test]
fn kappa_aware_asymmetric_arms_biases_cells_to_wider_peripheral() {
    // Asymmetric trifurcation: left peripheral wider than right.
    // Cancer cells should preferentially go to the wider LEFT arm (higher q_l).
    let [q_c, q_l, q_r] = tri_asymmetric_q_fracs(0.40, 0.40, length(4e-3), length(1e-3));
    // q_l > q_r → more cells should route to left
    let w_center = 0.40 * 4e-3;
    let h = 1e-3_f64;
    let dh = 2.0 * w_center * h / (w_center + h);
    let stage = CascadeStage {
        arm_q_fracs: [q_c, q_l, q_r, 0.0, 0.0],
        n_arms: 3,
        treatment_dh_m: Length::from_base(dh),
        parent_v_in_m_s: Velocity::from_base(0.02),
        peripheral_recoveries: [None; 4],
        n_recoveries: 0,
    };
    let n = stage.n_arms as usize;
    let dh_s = stage.treatment_dh_m.into_base().max(1e-9);
    let beta = beta_kappa_adjusted(SE_CANCER, D_CANCER_M / dh_s, 0.02, false);
    let p_to_left = p_arm_general(&stage.arm_q_fracs[..n], 1, beta);
    let p_to_right = p_arm_general(&stage.arm_q_fracs[..n], 2, beta);
    assert!(
        p_to_left > p_to_right,
        "wider left peripheral should attract more cancer cells: p_left={p_to_left:.4}, p_right={p_to_right:.4}"
    );
}

#[test]
fn kappa_aware_deeper_cascade_still_improves_separation() {
    // Build a 3-stage TriTriTri with narrowing Dh at each stage.
    let center_frac = 0.45_f64;
    let h = 1e-3_f64;
    let mut parent_w = 4e-3_f64;
    let mut stages = Vec::new();
    for _ in 0..3 {
        let w_c = center_frac * parent_w;
        let dh = 2.0 * w_c * h / (w_c + h);
        let q_c = tri_center_q_frac_cross_junction(center_frac, length(parent_w), length(h));
        let q_p = (1.0 - q_c) / 2.0;
        stages.push(CascadeStage {
            arm_q_fracs: [q_c, q_p, q_p, 0.0, 0.0],
            n_arms: 3,
            treatment_dh_m: Length::from_base(dh),
            parent_v_in_m_s: Velocity::from_base(0.02),
            peripheral_recoveries: [None; 4],
            n_recoveries: 0,
        });
        parent_w *= center_frac;
    }

    let r3 = mixed_cascade_separation_kappa_aware(&stages);
    let r1 = mixed_cascade_separation_kappa_aware(&stages[..1]);

    assert!(
        r3.rbc_peripheral_fraction >= r1.rbc_peripheral_fraction - 1e-10,
        "3-stage TriTriTri should push more RBCs to periphery: r1={:.4}, r3={:.4}",
        r1.rbc_peripheral_fraction,
        r3.rbc_peripheral_fraction
    );
    assert!(
        r3.separation_efficiency >= r1.separation_efficiency - 1e-10,
        "3-stage TriTriTri should improve separation efficiency"
    );
}

