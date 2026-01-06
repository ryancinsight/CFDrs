//! Comprehensive Spalart-Allmaras turbulence model validation tests
//!
//! Tests cover one-equation model, wall distance effects, trip terms, and
//! verification against analytical solutions and literature benchmarks.
//!
//! References:
//! - Spalart, P. R., & Allmaras, S. R. (1994). "A one-equation turbulence model...". AIAA Paper 92-0439

use cfd_2d::physics::turbulence::constants::{SA_CB1, SA_CW1};
use cfd_2d::physics::turbulence::spalart_allmaras::SpalartAllmaras;

/// Test SA model core functionality with public interface
#[test]
fn test_sa_model_basic() {
    let model = SpalartAllmaras::<f64>::new(10, 10);

    let nu_tilde = 1e-4;
    let molecular_viscosity = 1e-5;
    let s_tilde = 100.0;
    let wall_distance = 0.01;

    // Test core functionality with public interface
    let production = model.production(nu_tilde, s_tilde);
    assert!(production > 0.0, "Production should be positive");

    let destruction = model.destruction(nu_tilde, wall_distance, 0.5);
    assert!(destruction > 0.0, "Destruction should be positive");

    let eddy_visc = model.eddy_viscosity(nu_tilde, molecular_viscosity);
    assert!(eddy_visc >= 0.0, "Eddy viscosity should be non-negative");

    let wall_distances = model.wall_distance_field(0.1, 0.1);
    assert_eq!(
        wall_distances.len(),
        10 * 10,
        "Wall distance field size correct"
    );
}

/// Test SA coefficient validation
#[test]
fn test_sa_coefficient_validation() {
    let model = SpalartAllmaras::<f64>::new(5, 5);

    // Test via computational results (literature validation)
    let nu_tilde = 1e-4;
    let s_tilde = 100.0;

    let production = model.production(nu_tilde, s_tilde);
    let expected_production = SA_CB1 * s_tilde * nu_tilde;
    assert!((production - expected_production).abs() < 1e-10);

    let destruction = model.destruction(nu_tilde, 0.01, 0.5);
    let expected_destruction = SA_CW1 * 0.5 * (nu_tilde / 0.01).powi(2);
    assert!((destruction - expected_destruction).abs() < 1e-6);
}

/// Test wall distance computation
#[test]
fn test_sa_wall_distance() {
    let model = SpalartAllmaras::<f64>::new(21, 21);
    let distances = model.wall_distance_field(0.1, 0.1);

    assert_eq!(distances.len(), 21 * 21);
    // Center of domain should have maximum wall distance
    let center_idx = 10 * 21 + 10;
    assert!(distances[center_idx] > distances[0]); // Center > corner
}

/// Test trip term (should be zero)
#[test]
fn test_sa_trip_term() {
    let model = SpalartAllmaras::<f64>::new(5, 5);
    let trip = model.trip_term(1e-4, 0.01);
    assert_eq!(trip, 0.0);
}

/// Test asymptotic behavior of eddy viscosity
#[test]
fn test_eddy_viscosity_asymptotics() {
    let model = SpalartAllmaras::<f64>::new(10, 10);
    let nu = 1e-5;

    // Small ratio (laminar limit)
    let nu_tilde_small = 1e-8;
    let nu_t_small = model.eddy_viscosity(nu_tilde_small, nu);
    assert!(
        nu_t_small < nu_tilde_small,
        "Eddy viscosity should be much smaller in laminar limit"
    );

    // Large ratio (turbulent limit)
    let nu_tilde_large = 1e-3;
    let nu_t_large = model.eddy_viscosity(nu_tilde_large, nu);
    let ratio = nu_t_large / nu_tilde_large;
    assert!(
        ratio > 0.95,
        "Eddy viscosity should approach molecular limit in turbulent limit"
    );
}
