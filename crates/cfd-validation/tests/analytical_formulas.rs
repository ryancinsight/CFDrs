//! Tests for analytical formula calculations
//!
//! These tests verify that the analytical solutions compute the expected values
//! for known cases. They do not validate simulation results against analytical
//! solutions - that would be a true validation test.

#[test]
fn test_calculate_poiseuille_avg_velocity() {
    // Test that Poiseuille flow formula calculates correctly
    // v(r) = (ΔP/4μL) * (R² - r²)
    let radius = 0.01_f64; // 1 cm
    let length = 1.0_f64; // 1 m
    let pressure_drop = 100.0_f64; // Pa
    let viscosity = 0.001_f64; // Pa·s

    // Maximum velocity at centerline
    let v_max = (pressure_drop * radius * radius) / (4.0 * viscosity * length);

    // Average velocity for parabolic profile
    let v_avg = v_max / 2.0;

    // Verify the calculation produces expected result
    assert!(
        (v_avg - 1.25).abs() < 1e-6,
        "Poiseuille average velocity calculation: expected 1.25, got {v_avg}"
    );
}

#[test]
fn test_calculate_reynolds_number() {
    // Test Reynolds number calculation
    let density = 1000.0_f64; // kg/m³
    let velocity = 1.0_f64; // m/s
    let length = 0.1_f64; // m
    let viscosity = 0.001_f64; // Pa·s

    let re = density * velocity * length / viscosity;

    assert!(
        (re - 100000.0).abs() < 1.0,
        "Reynolds number calculation: expected 100000, got {re}"
    );
}
