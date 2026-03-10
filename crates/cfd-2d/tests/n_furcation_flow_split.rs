//! Cross-fidelity validation: 2D N-furcation flow split vs fundamental continuity.
//!
//! Evaluates the generalized NFurcationSolver2D to ensure it accurately resolves
//! 2D velocity fields for junctions with > 2 daughter branches (Quadfurcation/Pentafurcation).

use cfd_2d::solvers::{NFurcationGeometry, NFurcationSolver2D};
use cfd_core::physics::fluid::BloodModel;
use cfd_math::pressure_velocity::SIMPLEConfig;

/// Validates a symmetric quadfurcation (4-way branch: 1 parent, 3 daughters = 4 total arms? 
/// Wait: JunctionFamily defines Quadfurcation = 4 branches total (1 parent, 3 daughters).
/// Pentafurcation = 5 branches total (1 parent, 4 daughters).
/// We'll test 1 parent -> 3 daughters (Quadfurcation) here.)
#[test]
fn symmetric_quadfurcation_mass_conservation() {
    let parent_w = 2.0e-3_f64;
    let parent_l = 4.0e-3;
    let daughter_w = 0.5e-3;
    let daughter_l = 4.0e-3;
    let spread_angle = 1.0_f64; // radians spread between outermost branches
    
    // quadfurcation: 1 parent, 3 daughters
    let geom = NFurcationGeometry::new_symmetric(
        parent_w, parent_l, daughter_w, daughter_l, 3, spread_angle,
    );

    let blood = BloodModel::Newtonian(3.5e-3);
    let density = 1060.0;
    let u_inlet = 0.005;

    let config = SIMPLEConfig {
        max_iterations: 300,
        ..SIMPLEConfig::default()
    };

    let mut solver = NFurcationSolver2D::new(geom, blood, density, 60, 40, config);
    let result = solver
        .solve(u_inlet)
        .expect("quadfurcation solve should succeed");

    // Verify mass conservation.
    assert!(
        result.mass_balance_error < 0.10,
        "mass balance error {:.4} exceeds 10% for Quadfurcation",
        result.mass_balance_error
    );

    // Ensure flow actually left the daughters
    assert!(
        result.q_total_out > 0.0,
        "total out flux should be strictly positive"
    );
}

/// Validates a symmetric pentafurcation (5-way branch: 1 parent, 4 daughters)
#[test]
fn symmetric_pentafurcation_mass_conservation() {
    let parent_w = 2.0e-3_f64;
    let parent_l = 4.0e-3;
    let daughter_w = 0.4e-3;
    let daughter_l = 4.0e-3;
    let spread_angle = 1.2_f64; 
    
    // pentafurcation: 1 parent, 4 daughters
    let geom = NFurcationGeometry::new_symmetric(
        parent_w, parent_l, daughter_w, daughter_l, 4, spread_angle,
    );

    let blood = BloodModel::Newtonian(3.5e-3);
    let density = 1060.0;
    let u_inlet = 0.005;

    let config = SIMPLEConfig {
        max_iterations: 300,
        ..SIMPLEConfig::default()
    };

    let mut solver = NFurcationSolver2D::new(geom, blood, density, 60, 40, config);
    let result = solver
        .solve(u_inlet)
        .expect("pentafurcation solve should succeed");

    // Verify mass conservation.
    assert!(
        result.mass_balance_error < 0.10,
        "mass balance error {:.4} exceeds 10% for Pentafurcation",
        result.mass_balance_error
    );
    
    assert!(
        result.q_total_out > 0.0,
        "total out flux should be strictly positive"
    );
}
