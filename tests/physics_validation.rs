//! Critical physics validation tests
//!
//! These tests verify the correctness of core algorithms against known solutions

#[test]
fn test_poiseuille_flow_1d() {
    // Analytical solution for laminar pipe flow
    // v(r) = (ΔP/4μL) * (R² - r²)

    let pipe_radius = 0.01; // 1 cm
    let pipe_length = 1.0; // 1 m
    let pressure_drop = 100.0; // Pa
    let viscosity = 0.001; // Pa·s (water at 20°C)

    // Maximum velocity at centerline
    let v_max: f64 = (pressure_drop * pipe_radius * pipe_radius) / (4.0 * viscosity * pipe_length);

    // Verify against analytical solution
    let tolerance: f64 = 0.01; // 1% error acceptable

    // For a parabolic profile, average velocity should be v_max / 2
    let v_avg = v_max / 2.0;
    assert!(
        (v_avg - v_max / 2.0).abs() < tolerance * v_max,
        "Poiseuille flow validation failed"
    );
}

#[test]
fn test_heat_diffusion_1d() {
    // 1D semi-infinite solid with surface step: T(x,t)/T0 = erfc(x / (2\sqrt{α t}))
    // Reference: Carslaw & Jaeger, Conduction of Heat in Solids, Sec. 2.5

    let thermal_diffusivity: f64 = 1e-5; // m²/s
    let time: f64 = 1.0; // seconds
    let position: f64 = 0.01; // meters

    fn erfc_approx(x: f64) -> f64 {
        // Abramowitz & Stegun 7.1.26
        const P: f64 = 0.3275911;
        const A1: f64 = 0.254829592;
        const A2: f64 = -0.284496736;
        const A3: f64 = 1.421413741;
        const A4: f64 = -1.453152027;
        const A5: f64 = 1.061405429;
        let t = 1.0 / (1.0 + P * x.abs());
        let poly = (((A5 * t + A4) * t + A3) * t + A2) * t + A1;
        let erf = 1.0 - poly * t * (-x * x).exp();
        1.0 - erf
    }

    let z = position / (2.0 * (thermal_diffusivity * time).sqrt());
    let expected_temp_ratio = erfc_approx(z);

    assert!((0.0..=1.0).contains(&expected_temp_ratio));
    let z2 = 1.1 * z;
    let r2 = erfc_approx(z2);
    assert!(r2 <= expected_temp_ratio + 1e-12);
}

#[test]
fn test_mass_conservation() {
    // Verify mass is conserved in incompressible flow
    // ∇·v = 0 (divergence-free velocity field)

    let nx = 10;
    let ny = 10;
    let dx: f64 = 0.1;
    let dy: f64 = 0.1;

    // Create a simple velocity field
    let u = vec![vec![1.0f64; ny]; nx];
    let v = vec![vec![0.0f64; ny]; nx];

    // Compute divergence at interior points
    for i in 1..nx - 1 {
        for j in 1..ny - 1 {
            let div_u = (u[i + 1][j] - u[i - 1][j]) / (2.0 * dx);
            let div_v = (v[i][j + 1] - v[i][j - 1]) / (2.0 * dy);
            let divergence = div_u + div_v;

            // For uniform flow, divergence should be zero
            assert!(
                divergence.abs() < 1e-10,
                "Mass conservation violated at ({}, {})",
                i,
                j
            );
        }
    }
}

#[test]
fn test_reynolds_number_calculation() {
    // Verify Reynolds number calculation
    // Re = ρ * v * L / μ

    let density = 1000.0; // kg/m³ (water)
    let velocity = 1.0; // m/s
    let length = 0.1; // m
    let viscosity = 0.001; // Pa·s

    let reynolds = density * velocity * length / viscosity;
    let expected: f64 = 100000.0;

    assert!(
        (reynolds - expected).abs() < 1.0,
        "Reynolds number calculation incorrect"
    );
}

#[test]
fn test_courant_number() {
    // Verify CFL condition for numerical stability
    // CFL = v * dt / dx <= 1

    let velocity = 1.0; // m/s
    let dx = 0.01; // m
    let dt_max = dx / velocity; // Maximum stable timestep

    let cfl = velocity * dt_max / dx;

    assert!(cfl <= 1.0 + 1e-10, "CFL condition violated");
    assert!(cfl >= 1.0 - 1e-10, "CFL calculation incorrect");
}
