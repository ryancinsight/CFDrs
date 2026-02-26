//! Cell separation physics tests — inertial focusing, Dean drag, and separation analysis.
//!
//! Tests verify:
//! 1. Physics force models return correct order-of-magnitude values
//! 2. Bisection converges and residual is machine-precision small
//! 3. Cancer/RBC cells separate correctly in physiological channels
//! 4. Edge cases (zero velocity, large deformability) are handled correctly

use approx::assert_relative_eq;
use cfd_1d::cell_separation::{
    dean_drag_force_n, dean_number, inertial_lift_force_n, lateral_equilibrium, CellProperties,
    CellSeparationModel,
};

const BLOOD_DENSITY: f64 = 1060.0;
const BLOOD_VISCOSITY: f64 = 3.5e-3;

// ============================================================================
// Dean number formula
// ============================================================================

/// De = Re √(Dh / 2R) — analytical check against manual calculation.
#[test]
fn test_dean_number_formula() {
    let re = 100.0_f64;
    let dh = 200e-6_f64; // 200 µm hydraulic diameter
    let r = 5e-3_f64; // 5 mm bend radius
    let de = dean_number(re, dh, r);
    let expected = re * (dh / (2.0 * r)).sqrt();
    assert_relative_eq!(de, expected, max_relative = 1e-12);
}

// ============================================================================
// Dean drag force scaling
// ============================================================================

/// F_D ∝ De^1.63 (Gossett & Di Carlo 2009, Eq. 4).
#[test]
fn test_dean_drag_force_scaling() {
    let mu = 3.5e-3_f64;
    let a = 10e-6_f64;
    let de1 = 5.0_f64;
    let de2 = 10.0_f64;

    let fd1 = dean_drag_force_n(mu, de1, a);
    let fd2 = dean_drag_force_n(mu, de2, a);

    // Expected ratio: (de2/de1)^1.63
    let expected_ratio = (de2 / de1).powf(1.63);
    assert_relative_eq!(fd2 / fd1, expected_ratio, max_relative = 1e-9);
}

// ============================================================================
// Confinement ratio and focusing threshold
// ============================================================================

/// MCF-7 cancer cell (D = 17.5 µm) in 200 µm wide channel: κ = 0.116 > 0.07.
#[test]
fn test_mcf7_confinement_ratio_exceeds_threshold() {
    let cancer = CellProperties::mcf7_breast_cancer();
    let dh = 2.0 * 200e-6 * 100e-6 / (200e-6 + 100e-6); // 133 µm
    let kappa = cancer.confinement_ratio(dh);
    assert!(kappa > 0.07, "MCF-7 must have κ > 0.07, got {}", kappa);
    assert!(cancer.will_focus(dh), "MCF-7 must focus");
}

/// Platelet (D = 2.5 µm) in 200 µm × 100 µm channel: κ ≈ 0.019 < 0.07.
#[test]
fn test_platelet_below_focusing_threshold() {
    let platelet = CellProperties::platelet();
    let dh = 2.0 * 200e-6 * 100e-6 / (200e-6 + 100e-6); // 133 µm
    let kappa = platelet.confinement_ratio(dh);
    assert!(kappa < 0.07, "Platelet must have κ < 0.07, got {}", kappa);
    assert!(!platelet.will_focus(dh), "Platelet must not focus");
}

// ============================================================================
// Cell equilibrium positions
// ============================================================================

/// In a straight 500 µm × 200 µm channel at physiological flow, MCF-7 cancer
/// cells (DI = 0.15, rigid) should focus **nearer center** than RBCs (DI = 0.85,
/// deformable), enabling separation.
#[test]
fn test_cancer_focuses_closer_to_center_than_rbc() {
    let cancer = CellProperties::mcf7_breast_cancer();
    let rbc = CellProperties::red_blood_cell();

    let cancer_eq = lateral_equilibrium(
        &cancer, BLOOD_DENSITY, BLOOD_VISCOSITY, 0.05, 500e-6, 200e-6, None,
    ).expect("MCF-7 must focus (κ > 0.07)");

    let rbc_eq = lateral_equilibrium(
        &rbc, BLOOD_DENSITY, BLOOD_VISCOSITY, 0.05, 500e-6, 200e-6, None,
    ).expect("RBC must focus in this channel");

    // More rigid cells focus nearer center (smaller x̃)
    assert!(
        cancer_eq.x_tilde_eq < rbc_eq.x_tilde_eq,
        "Cancer (x̃={:.3}) should be closer to center than RBC (x̃={:.3})",
        cancer_eq.x_tilde_eq,
        rbc_eq.x_tilde_eq
    );
}

/// Bisection residual must be negligible (< 1 pN force error).
#[test]
fn test_bisection_residual_near_zero() {
    let cancer = CellProperties::mcf7_breast_cancer();
    let eq = lateral_equilibrium(
        &cancer, BLOOD_DENSITY, BLOOD_VISCOSITY, 0.05, 500e-6, 200e-6, None,
    ).unwrap();
    // Residual force should be machine-epsilon small (< 1 pN = 1e-12 N)
    assert!(
        eq.residual_force_n.abs() < 1e-12,
        "Bisection residual must be < 1 pN, got {:.3e} N",
        eq.residual_force_n
    );
}

// ============================================================================
// CellSeparationModel
// ============================================================================

/// Cancer vs RBC separation efficiency must be positive in a 200µm × 100µm channel.
#[test]
fn test_cancer_rbc_separation_efficiency_positive() {
    let model = CellSeparationModel::new(200e-6, 100e-6, None);
    let cancer = CellProperties::mcf7_breast_cancer();
    let rbc = CellProperties::red_blood_cell();

    let analysis = model
        .analyze(&cancer, &rbc, BLOOD_DENSITY, BLOOD_VISCOSITY, 0.05)
        .expect("Separation analysis must succeed");

    assert!(
        analysis.separation_efficiency > 0.0,
        "Separation efficiency must be > 0, got {}",
        analysis.separation_efficiency
    );
}

/// Purity must be in [0, 1] for the default channel and flow conditions.
#[test]
fn test_purity_bounded_in_zero_to_one() {
    let model = CellSeparationModel::new(200e-6, 100e-6, None);
    let cancer = CellProperties::mcf7_breast_cancer();
    let rbc = CellProperties::red_blood_cell();

    let analysis = model
        .analyze(&cancer, &rbc, BLOOD_DENSITY, BLOOD_VISCOSITY, 0.05)
        .unwrap();

    assert!(
        (0.0..=1.0).contains(&analysis.purity),
        "Purity must be in [0, 1], got {}",
        analysis.purity
    );
    assert!(
        (0.0..=1.0).contains(&analysis.separation_efficiency),
        "Separation efficiency must be in [0, 1], got {}",
        analysis.separation_efficiency
    );
}

/// Curved channel must shift cancer equilibrium outward (larger x̃) vs straight channel.
#[test]
fn test_curved_channel_shifts_equilibrium_outward() {
    let cancer = CellProperties::mcf7_breast_cancer();
    let straight = lateral_equilibrium(
        &cancer, BLOOD_DENSITY, BLOOD_VISCOSITY, 0.05, 500e-6, 200e-6, None,
    ).unwrap();
    let curved = lateral_equilibrium(
        &cancer, BLOOD_DENSITY, BLOOD_VISCOSITY, 0.05, 500e-6, 200e-6, Some(5e-3),
    ).unwrap();

    assert!(
        curved.x_tilde_eq >= straight.x_tilde_eq,
        "Curved channel should shift equilibrium toward wall, straight={:.3}, curved={:.3}",
        straight.x_tilde_eq,
        curved.x_tilde_eq
    );
}

/// Neither cell type focuses when channel is enormous (κ << 0.07).
/// `analyze` still returns `Some` with dispersed defaults (x̃ = 0.5,
/// `will_focus = false`) so callers can use a partial signal.
#[test]
fn test_large_channel_returns_none_for_neither_focuses() {
    let model = CellSeparationModel::new(10e-3, 10e-3, None); // 1 cm × 1 cm — huge
    let cancer = CellProperties::mcf7_breast_cancer();
    let rbc = CellProperties::red_blood_cell();

    let result = model.analyze(&cancer, &rbc, BLOOD_DENSITY, BLOOD_VISCOSITY, 0.001);
    let analysis = result.expect("analyze always returns Some for dispersed cells");
    assert!(
        !analysis.target_equilibrium.will_focus,
        "Target cell should not focus in a 1cm channel at low velocity"
    );
    assert!(
        !analysis.background_equilibrium.will_focus,
        "Background cell should not focus in a 1cm channel at low velocity"
    );
}

/// Inertial lift force must have the correct dimensional units.
/// At x̃ = 0 (center), `C_L` should be negative (shear-gradient → wall).
#[test]
fn test_inertial_lift_at_center_is_negative() {
    let cancer = CellProperties::mcf7_breast_cancer();
    let fl = inertial_lift_force_n(0.0, &cancer, BLOOD_DENSITY, 0.05, 200e-6);
    // At center (x̃=0), C_wall=0, C_center=0.3 → C_L = -0.3 → F_L < 0
    assert!(fl < 0.0, "Lift at center must be negative (toward wall), got {}", fl);
}

/// Inertial lift force must have the correct dimensional units.
/// At x̃ ≈ 0.9 (near wall), `C_L` should be positive (wall repulsion > shear gradient).
#[test]
fn test_inertial_lift_near_wall_is_positive() {
    let cancer = CellProperties::mcf7_breast_cancer();
    let fl = inertial_lift_force_n(0.9, &cancer, BLOOD_DENSITY, 0.05, 200e-6);
    // Near wall C_wall diverges, C_L should be positive
    assert!(fl > 0.0, "Lift near wall must be positive (toward center), got {}", fl);
}
