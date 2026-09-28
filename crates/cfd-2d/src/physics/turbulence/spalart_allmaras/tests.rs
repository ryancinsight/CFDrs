use super::*;
use crate::physics::turbulence::constants::{SA_CB1, SA_CW1};
use crate::physics::turbulence::traits::TurbulenceModel;
use eunomia::assert_relative_eq;

#[test]
fn test_spalart_allmaras_creation() {
    let model = SpalartAllmaras::<f64>::new(10, 10);
    assert_eq!(model.nx, 10);
    assert_eq!(model.ny, 10);
}

#[test]
fn test_eddy_viscosity_zero_nu_tilde() {
    let model = SpalartAllmaras::<f64>::new(10, 10);
    let nu_t = model.eddy_viscosity(0.0, 1e-5);
    assert_eq!(nu_t, 0.0);
}

#[test]
fn test_eddy_viscosity_large_nu_tilde() {
    let model = SpalartAllmaras::<f64>::new(10, 10);
    let nu = 1e-5;
    let nu_tilde = 1e-3;
    let nu_t = model.eddy_viscosity(nu_tilde, nu);

    // For χ >> 1, fv1 → 1, so νt → ν̃
    assert!(nu_t > 0.9 * nu_tilde);
    assert!(nu_t <= nu_tilde);
}

#[test]
fn test_vorticity_magnitude() {
    let model = SpalartAllmaras::<f64>::new(10, 10);

    // Test case: uniform rotation
    let velocity_gradient = [[0.0, 1.0], [-1.0, 0.0]];
    let vorticity = model.vorticity_magnitude(&velocity_gradient);

    // Ω = |∂v/∂x - ∂u/∂y| = |-1 - 1| = 2
    assert_relative_eq!(vorticity, 2.0, epsilon = 1e-10);
}

#[test]
fn test_production_term() {
    let model = SpalartAllmaras::<f64>::new(10, 10);
    let nu_tilde = 1e-4;
    let s_tilde = 100.0;

    let production = model.production(nu_tilde, s_tilde);

    // P = Cb1 * S̃ * ν̃
    let expected = SA_CB1 * s_tilde * nu_tilde;
    assert_relative_eq!(production, expected, epsilon = 1e-10);
}

#[test]
fn test_destruction_term() {
    let model = SpalartAllmaras::<f64>::new(10, 10);
    let nu_tilde = 1e-4;
    let wall_distance = 0.01;
    let fw = 1.0;

    let destruction = model.destruction(nu_tilde, wall_distance, fw);

    // D = Cw1 * fw * (ν̃/d)²
    let ratio = nu_tilde / wall_distance;
    let expected = SA_CW1 * fw * ratio * ratio;
    assert_relative_eq!(destruction, expected, epsilon = 1e-10);
}

#[test]
fn test_dissipation_term_mapping() {
    let model = SpalartAllmaras::<f64>::new(10, 10);
    let nu_tilde = 1e-4;
    let wall_distance = 0.01;
    let dissipation = model.dissipation_term(nu_tilde, wall_distance);
    let ratio = nu_tilde / wall_distance;
    let expected = SA_CW1 * ratio * ratio;
    assert_relative_eq!(dissipation, expected, epsilon = 1e-10);
}

#[test]
fn test_cbrt_function() {
    // Test cube root computation via helper
    use super::wall_distance::cbrt;
    assert_relative_eq!(cbrt(8.0), 2.0, epsilon = 1e-8);
    assert_relative_eq!(cbrt(27.0), 3.0, epsilon = 1e-8);
    assert_relative_eq!(cbrt(1.0), 1.0, epsilon = 1e-8);
    assert_eq!(cbrt(0.0), 0.0);
}

#[test]
fn test_wall_distance_field() {
    let model = SpalartAllmaras::<f64>::new(5, 5);
    let distances = model.wall_distance_field(0.1, 0.1);

    assert_eq!(distances.len(), 25);

    // Corner should have minimum distance
    // Center should have maximum distance
    let center_idx = 2 * 5 + 2;
    let corner_idx = 0;

    assert!(distances[center_idx] > distances[corner_idx]);
    assert_relative_eq!(distances[corner_idx], 0.05, epsilon = 1e-10);
}

#[test]
fn test_boundary_conditions() {
    let model = SpalartAllmaras::<f64>::new(5, 5);
    let mut nu_tilde = vec![1.0; 25];

    model.apply_boundary_conditions(&mut nu_tilde);

    // Check walls are zero
    for i in 0..5 {
        assert_eq!(nu_tilde[i], 0.0); // Bottom
        assert_eq!(nu_tilde[4 * 5 + i], 0.0); // Top
    }
    for j in 0..5 {
        assert_eq!(nu_tilde[j * 5], 0.0); // Left
        assert_eq!(nu_tilde[j * 5 + 4], 0.0); // Right
    }
}

#[test]
fn test_trip_term_zero() {
    let model = SpalartAllmaras::<f64>::new(10, 10);
    let trip = model.trip_term(1e-4, 0.01);

    // For fully turbulent flows, trip term should be zero
    assert_eq!(trip, 0.0);
}
