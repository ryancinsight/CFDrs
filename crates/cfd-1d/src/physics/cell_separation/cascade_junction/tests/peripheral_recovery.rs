//! Peripheral recovery sub-split routing tests.

use super::*;

// ── Peripheral recovery routing tests ────────────────────────────────

#[test]
fn recovery_zero_when_no_peripheral_recoveries() {
    let q = tri_center_q_frac(0.45);
    let q_p = (1.0 - q) / 2.0;
    let dh = 1e-3;
    let stage_no_recovery = CascadeStage {
        arm_q_fracs: [q, q_p, q_p, 0.0, 0.0],
        n_arms: 3,
        treatment_dh_m: Length::from_base(dh),
        parent_v_in_m_s: Velocity::from_base(0.02),
        peripheral_recoveries: [None; 4],
        n_recoveries: 0,
    };
    let r = mixed_cascade_separation_kappa_aware(&[stage_no_recovery]);
    let r_legacy = mixed_cascade_separation(&[(q, true)]);
    // With no recovery, kappa-aware should match or exceed legacy
    assert!(r.cancer_center_fraction >= r_legacy.cancer_center_fraction - 1e-10);
}

#[test]
fn recovery_increases_cancer_center_fraction() {
    // 20/60/20 trifurcation with recovery on left peripheral
    let q_c = 0.65; // center gets ~65% of flow (wider channel)
    let q_p = (1.0 - q_c) / 2.0;
    let dh = 1e-3;
    let stage_no_recovery = CascadeStage {
        arm_q_fracs: [q_c, q_p, q_p, 0.0, 0.0],
        n_arms: 3,
        treatment_dh_m: Length::from_base(dh),
        parent_v_in_m_s: Velocity::from_base(0.02),
        peripheral_recoveries: [None; 4],
        n_recoveries: 0,
    };
    // Recovery sub-split on arm 1 (left): 70/30 split, wider arm (index 0) feeds back to treatment
    let recovery = PeripheralRecovery {
        source_arm_idx: 1,
        sub_arm_q_fracs: [0.70, 0.30, 0.0, 0.0, 0.0],
        n_sub_arms: 2,
        recovery_arm_idx: 0,
        recovery_dh_m: Length::from_base(0.5e-3),
    };
    let stage_with_recovery = CascadeStage {
        arm_q_fracs: [q_c, q_p, q_p, 0.0, 0.0],
        n_arms: 3,
        treatment_dh_m: Length::from_base(dh),
        parent_v_in_m_s: Velocity::from_base(0.02),
        peripheral_recoveries: [Some(recovery), None, None, None],
        n_recoveries: 1,
    };
    let r_no = mixed_cascade_separation_kappa_aware(&[stage_no_recovery]);
    let r_yes = mixed_cascade_separation_kappa_aware(&[stage_with_recovery]);
    assert!(
        r_yes.cancer_center_fraction > r_no.cancer_center_fraction,
        "recovery should increase cancer center fraction: without={:.4}, with={:.4}",
        r_no.cancer_center_fraction,
        r_yes.cancer_center_fraction
    );
}

#[test]
fn recovery_bounded_by_one() {
    // Even with aggressive recovery, P_eff must not exceed 1.0
    let q_c = 0.50;
    let q_p = 0.25;
    let recovery_both = [
        Some(PeripheralRecovery {
            source_arm_idx: 1,
            sub_arm_q_fracs: [0.90, 0.10, 0.0, 0.0, 0.0],
            n_sub_arms: 2,
            recovery_arm_idx: 0,
            recovery_dh_m: Length::from_base(0.3e-3),
        }),
        Some(PeripheralRecovery {
            source_arm_idx: 2,
            sub_arm_q_fracs: [0.90, 0.10, 0.0, 0.0, 0.0],
            n_sub_arms: 2,
            recovery_arm_idx: 0,
            recovery_dh_m: Length::from_base(0.3e-3),
        }),
        None,
        None,
    ];
    let stage = CascadeStage {
        arm_q_fracs: [q_c, q_p, q_p, 0.0, 0.0],
        n_arms: 3,
        treatment_dh_m: Length::from_base(1e-3),
        parent_v_in_m_s: Velocity::from_base(0.02),
        peripheral_recoveries: recovery_both,
        n_recoveries: 2,
    };
    let r = mixed_cascade_separation_kappa_aware(&[stage]);
    assert!(
        r.cancer_center_fraction <= 1.0,
        "cancer fraction must be <= 1.0"
    );
    assert!(r.wbc_center_fraction <= 1.0, "wbc fraction must be <= 1.0");
    assert!(
        (1.0 - r.rbc_peripheral_fraction) <= 1.0,
        "rbc center must be <= 1.0"
    );
}

