//! Core Zweifach-Fung probability and beta-amplification unit tests.

use super::*;

// ── Core Zweifach-Fung function unit tests ──────────────────────────

#[test]
fn p_center_equal_trifurcation() {
    // Equal trifurcation: q_center = 1/3, each peripheral = 1/3
    // P_center = (1/3)^β / ((1/3)^β + 2*(1/3)^β) = 1/3 for any β
    use eunomia::assert_relative_eq;
    let p_cancer = p_center(1.0 / 3.0, SE_CANCER);
    let p_wbc = p_center(1.0 / 3.0, SE_WBC);
    let p_rbc = p_center(1.0 / 3.0, SE_RBC);
    assert_relative_eq!(p_cancer, 1.0 / 3.0, max_relative = 1e-9);
    assert_relative_eq!(p_wbc, 1.0 / 3.0, max_relative = 1e-9);
    assert_relative_eq!(p_rbc, 1.0 / 3.0, max_relative = 1e-9);
}

#[test]
fn p_center_biased_flow_stiff_cells_route_more() {
    // With q_center_frac > 1/3, stiffer cells (higher beta) route more to center.
    use eunomia::assert_relative_eq;
    let q = 0.5; // center gets 50% of flow
    let p_cancer = p_center(q, SE_CANCER);
    let p_rbc = p_center(q, SE_RBC);
    // For RBC (beta=1): P = 0.5^1 / (0.5^1 + 2*0.25^1) = 0.5/1.0 = 0.5
    assert_relative_eq!(p_rbc, 0.5, max_relative = 1e-9);
    // Cancer (beta=1.85) should be > 0.5
    assert!(
        p_cancer > 0.5,
        "stiff cancer cells should route more to center: got {p_cancer}"
    );
    assert!(
        p_cancer > p_rbc,
        "cancer p_center should exceed RBC p_center"
    );
}

#[test]
fn tri_center_q_frac_known_values() {
    // For center_frac = 0.5: w_c=0.5, w_p=0.25 each
    // q_frac = 0.5^3 / (0.5^3 + 2*0.25^3) = 0.125 / (0.125 + 0.03125) = 0.125/0.15625 = 0.8
    use eunomia::assert_relative_eq;
    let q = tri_center_q_frac(0.5);
    assert_relative_eq!(q, 0.8, max_relative = 1e-9);
}

#[test]
fn beta_kappa_adjusted_at_kappa_ref() {
    // At kappa = KAPPA_REF: amplification = 1 + 0.07/0.07 = 2
    // beta_eff = 1 + (SE_CANCER - 1) * 2 = 1 + 0.85 * 2 = 2.70
    use eunomia::assert_relative_eq;
    let b = beta_kappa_adjusted(SE_CANCER, KAPPA_REF, 0.02, false);
    assert_relative_eq!(b, 1.0 + 0.85 * 2.0, max_relative = 1e-10);
}

#[test]
fn beta_kappa_adjusted_capped_at_three() {
    // Very large kappa should cap at 3.0
    let b = beta_kappa_adjusted(2.5, KAPPA_REF * 10.0, 0.02, false);
    assert!(b <= 3.0, "beta should be capped at 3.0, got {b}");
    // The clamp on kappa is at 2*KAPPA_REF, so amplification = 1 + 2 = 3
    // excess = 1.5, beta_eff = 1 + 1.5*3 = 5.5, capped to 3.0
    use eunomia::assert_relative_eq;
    assert_relative_eq!(b, 3.0, epsilon = 1e-15);
}

#[test]
fn beta_kappa_adjusted_wbc_moderate_amplification() {
    // WBC at kappa = KAPPA_REF: amplification = 2
    // beta_eff = 1 + (1.40 - 1) * 2 = 1 + 0.80 = 1.80
    use eunomia::assert_relative_eq;
    let b = beta_kappa_adjusted(SE_WBC, KAPPA_REF, 0.02, false);
    assert_relative_eq!(b, 1.0 + 0.40 * 2.0, max_relative = 1e-10);
}
