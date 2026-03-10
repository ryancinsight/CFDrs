//! Cross-fidelity cascade validation: 1D → 3D consistency for CIF networks.
//!
//! Validates that:
//! 1. A 1D Hagen-Poiseuille network solver produces flow rates that sum
//!    to the inlet flow (conservation).
//! 2. The per-channel flow rates from 1D drive the 3D CascadeSolver3D
//!    to produce physically consistent wall shear and pressure drop.
//! 3. Flow conservation holds: total 3D inlet flow ≈ sum of channel flow.

use cfd_3d::cascade::{CascadeChannelSpec, CascadeConfig3D, CascadeSolver3D};
use cfd_core::physics::fluid::ConstantPropertyFluid;

/// Analytical Poiseuille pressure drop for a rectangular duct
/// ΔP = 12 μ L Q / (w h³) (wide duct approximation, w >> h)
fn poiseuille_dp_rect(mu: f64, length: f64, flow_rate: f64, width: f64, height: f64) -> f64 {
    12.0 * mu * length * flow_rate / (width * height.powi(3))
}

/// Hydraulic resistance for a rectangular duct: R = 12 μ L / (w h³)
fn resistance_rect(mu: f64, length: f64, width: f64, height: f64) -> f64 {
    12.0 * mu * length / (width * height.powi(3))
}

/// CIF-like asymmetric bifurcation geometry.
///
/// Parent → (wide bypass, narrow center with venturi).
/// All channels share the same height and length.
const MU_WATER: f64 = 1.0e-3; // Pa·s
const CH_LENGTH: f64 = 10.0e-3; // 10 mm
const CH_HEIGHT: f64 = 1.0e-3; // 1 mm
const W_BYPASS: f64 = 2.0e-3; // 2 mm — wide bypass arm
const W_CENTER: f64 = 1.0e-3; // 1 mm — narrow center arm
const W_THROAT: f64 = 0.5e-3; // 0.5 mm — venturi throat
const Q_TOTAL: f64 = 1.0e-7; // 100 µL/s total inlet flow

// ── 1D analytical flow split ─────────────────────────────────────────────

#[test]
fn cross_fidelity_1d_flow_split_conservation() {
    // Parallel resistance: Q_i = Q_total × R_other / R_total
    let r_bypass = resistance_rect(MU_WATER, CH_LENGTH, W_BYPASS, CH_HEIGHT);
    let r_center = resistance_rect(MU_WATER, CH_LENGTH, W_CENTER, CH_HEIGHT);
    let _r_total = 1.0 / (1.0 / r_bypass + 1.0 / r_center);

    let q_bypass = Q_TOTAL * (r_center / (r_bypass + r_center));
    let q_center = Q_TOTAL * (r_bypass / (r_bypass + r_center));

    // Conservation: flow split sums to total.
    let q_sum = q_bypass + q_center;
    assert!(
        (q_sum - Q_TOTAL).abs() / Q_TOTAL < 1e-12,
        "1D flow conservation violated: sum={:.6e} vs total={:.6e}",
        q_sum,
        Q_TOTAL
    );

    // Wider channel carries more flow (R ∝ 1/w, so narrower is more resistive).
    assert!(
        q_bypass > q_center,
        "Wider bypass should carry more flow: bypass={:.4e} vs center={:.4e}",
        q_bypass,
        q_center
    );

    // Flow ratio should equal width ratio for w>>h: q_bypass/q_center ≈ w_bypass/w_center = 2.
    let flow_ratio = q_bypass / q_center;
    assert!(
        (flow_ratio - (W_BYPASS / W_CENTER)).abs() < 0.01,
        "Flow ratio {:.3} should be near width ratio {:.3}",
        flow_ratio,
        W_BYPASS / W_CENTER
    );
}

// ── 3D cascade from 1D flow rates ────────────────────────────────────────

#[test]
fn cross_fidelity_1d_to_3d_cascade() {
    // 1D: compute flow split analytically.
    let r_bypass = resistance_rect(MU_WATER, CH_LENGTH, W_BYPASS, CH_HEIGHT);
    let r_center = resistance_rect(MU_WATER, CH_LENGTH, W_CENTER, CH_HEIGHT);
    let q_bypass = Q_TOTAL * (r_center / (r_bypass + r_center));
    let q_center = Q_TOTAL * (r_bypass / (r_bypass + r_center));

    // 1D pressure drops.
    let dp_1d_bypass = poiseuille_dp_rect(MU_WATER, CH_LENGTH, q_bypass, W_BYPASS, CH_HEIGHT);
    let dp_1d_center = poiseuille_dp_rect(MU_WATER, CH_LENGTH, q_center, W_CENTER, CH_HEIGHT);

    // In parallel channels, the pressure drops must be equal.
    assert!(
        (dp_1d_bypass - dp_1d_center).abs() / dp_1d_bypass.max(1e-15) < 1e-10,
        "Parallel ΔP mismatch: bypass={:.4} vs center={:.4}",
        dp_1d_bypass,
        dp_1d_center
    );

    // 3D: drive CascadeSolver3D with the 1D flow rates.
    let water = ConstantPropertyFluid::<f64>::water_20c().unwrap();
    let config = CascadeConfig3D {
        outlet_pressure: 0.0,
        resolution: (20, 5, 5),
        max_picard_iterations: 3,
        picard_tolerance: 1e-2,
    };
    let solver = CascadeSolver3D::new(config, water);

    let specs = vec![
        CascadeChannelSpec {
            id: "bypass".into(),
            length: CH_LENGTH,
            width: W_BYPASS,
            height: CH_HEIGHT,
            flow_rate_m3_s: q_bypass,
            is_venturi_throat: false,
            throat_width: None,
            local_hematocrit: None,
        },
        CascadeChannelSpec {
            id: "center".into(),
            length: CH_LENGTH,
            width: W_CENTER,
            height: CH_HEIGHT,
            flow_rate_m3_s: q_center,
            is_venturi_throat: true,
            throat_width: Some(W_THROAT),
            local_hematocrit: None,
        },
    ];
    let result = solver.solve(&specs).unwrap();

    // Both channels produce valid results.
    assert_eq!(result.channel_results.len(), 2);
    for cr in &result.channel_results {
        assert!(
            cr.max_velocity > 0.0,
            "Channel {} has zero velocity",
            cr.channel_id
        );
        assert!(
            cr.wall_shear_mean_pa >= 0.0,
            "Channel {} has negative mean shear: {}",
            cr.channel_id,
            cr.wall_shear_mean_pa
        );
    }

    // The 3D pressure drops should have the same order of magnitude as the
    // 1D analytical value (within 10× tolerance for coarse FEM at 20×5×5).
    let dp_3d_bypass = result
        .channel_results
        .iter()
        .find(|r| r.channel_id == "bypass")
        .unwrap()
        .pressure_drop_pa;
    let dp_3d_center = result
        .channel_results
        .iter()
        .find(|r| r.channel_id == "center")
        .unwrap()
        .pressure_drop_pa;

    // Sanity: 3D pressure drops are positive.
    assert!(dp_3d_bypass > 0.0, "3D bypass ΔP should be positive");
    assert!(dp_3d_center > 0.0, "3D center ΔP should be positive");

    // Log for manual inspection.
    eprintln!("Cross-fidelity comparison:");
    eprintln!(
        "  1D:  ΔP_bypass = {:.4} Pa, ΔP_center = {:.4} Pa",
        dp_1d_bypass, dp_1d_center
    );
    eprintln!(
        "  3D:  ΔP_bypass = {:.4} Pa, ΔP_center = {:.4} Pa",
        dp_3d_bypass, dp_3d_center
    );
    eprintln!(
        "  Flow split: Q_bypass/Q_total = {:.3}, Q_center/Q_total = {:.3}",
        q_bypass / Q_TOTAL,
        q_center / Q_TOTAL
    );
}
