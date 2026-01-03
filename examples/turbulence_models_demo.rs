//! Demonstration of turbulence models for 2D CFD simulations
//!
//! This example shows how to use different turbulence models:
//! - k-ε model (standard two-equation model)
//! - k-ω SST model (Menter's Shear Stress Transport)
//! - Spalart-Allmaras model (one-equation model for aerospace)

use cfd_2d::fields::SimulationFields;
use cfd_2d::grid::StructuredGrid2D;
use cfd_2d::physics::turbulence::{
    KEpsilonModel, KOmegaSSTModel, SpalartAllmaras, TurbulenceModel,
};
use nalgebra::Vector2;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Turbulence Models Demonstration");
    println!("================================");

    // Create a simple 2D grid (20x20 cells)
    let grid = StructuredGrid2D::<f64>::new(20, 20, 0.0, 2.0, 0.0, 2.0)?;
    println!(
        "Grid: {}x{} cells, domain: [{:.1}, {:.1}] x [{:.1}, {:.1}]",
        grid.nx, grid.ny, 0.0, 2.0, 0.0, 2.0
    );

    // Initialize simulation fields
    let mut fields = SimulationFields::new(20, 20);

    // Set up test flow conditions (simple shear flow with some turbulence)
    setup_test_conditions(&mut fields);
    println!("Flow conditions: Mean velocity U=1.0, turbulence intensity 5%");

    // Demonstrate different turbulence models
    demonstrate_k_epsilon_model(&grid, &fields)?;
    demonstrate_k_omega_sst_model(&grid, &fields)?;
    demonstrate_spalart_allmaras_model(&grid, &fields)?;

    // Compare model predictions
    compare_model_predictions(&grid, &fields)?;

    println!("Demonstration completed successfully!");
    println!("Turbulence models provide:");
    println!("- Realistic eddy viscosity for momentum transport");
    println!("- Proper near-wall damping for boundary layer flows");
    println!("- Different formulations for different flow regimes");

    Ok(())
}

/// Set up test flow conditions with some turbulence
fn setup_test_conditions(fields: &mut SimulationFields<f64>) {
    let nx = 20;
    let ny = 20;

    // Mean flow: simple shear U = y (parabolic profile would be more realistic)
    for j in 0..ny {
        let y = j as f64 * 0.1; // y from 0 to 2
        let u = y; // Linear shear profile

        for i in 0..nx {
            fields.set_velocity_at(i, j, &Vector2::new(u, 0.0));
        }
    }

    // Turbulent quantities (typical values for moderate turbulence)
    let k_intensity = 0.05; // 5% turbulence intensity
    let u_ref = 1.0; // Reference velocity

    for idx in 0..(nx * ny) {
        let i = idx % nx;
        let j = idx / nx;

        // Turbulent kinetic energy: k = (3/2) * (u'²) ≈ 1.5 * (I*u_ref)²
        let k: f64 = 1.5 * ((k_intensity * u_ref) as f64).powi(2);

        // For k-ε: ε = C_μ^{3/4} * k^{3/2} / l, where l ≈ 0.07 * hydraulic diameter
        // For k-ω: ω = ε / (C_μ * k), with ε ≈ k^{3/2} / l
        let length_scale = 0.1; // Characteristic length
        #[allow(unused_variables)]
        let dissipation = k.powf(1.5) / length_scale; // ε ≈ k^{3/2} / l

        // Store in fields (using pressure for k)
        fields.p[(i, j)] = k;

        // For demo purposes, we'll simulate additional fields
        // In real code, you'd have dedicated fields for turbulence quantities
    }
}

fn demonstrate_k_epsilon_model(
    _grid: &StructuredGrid2D<f64>,
    _fields: &SimulationFields<f64>,
) -> Result<(), Box<dyn std::error::Error>> {
    println!("\n--- k-ε Turbulence Model ---");

    let model = KEpsilonModel::new(20, 20);

    // Typical values for moderate turbulence
    let k = 0.01; // Turbulent kinetic energy [m²/s²]
    let epsilon = 0.1; // Dissipation rate [m²/s³]
    let _density = 1.0; // Density [kg/m³]

    let nu_t = model.turbulent_viscosity(k, epsilon, _density);
    let p_k = model.production_term(&[[0.0, 1.0], [0.0, 0.0]], nu_t); // Simple shear
    let d_k = model.dissipation_term(k, epsilon);

    println!("Input conditions:");
    println!("  k = {:.4} m²/s² (turbulent kinetic energy)", k);
    println!("  ε = {:.4} m²/s³ (dissipation rate)", epsilon);
    println!("Model predictions:");
    println!("  ν_t = {:.6} m²/s (turbulent viscosity)", nu_t);
    println!("  P_k = {:.6} m²/s³ (production)", p_k);
    println!("  ε = {:.6} m²/s³ (dissipation)", d_k);
    println!("Characteristics:");
    println!("  - Two-equation model (k, ε)");
    println!("  - Good for free shear flows");
    println!("  - Requires ε boundary conditions");

    Ok(())
}

fn demonstrate_k_omega_sst_model(
    _grid: &StructuredGrid2D<f64>,
    _fields: &SimulationFields<f64>,
) -> Result<(), Box<dyn std::error::Error>> {
    println!("\n--- k-ω SST Turbulence Model ---");

    let model = KOmegaSSTModel::new(20, 20);

    // Typical values for moderate turbulence
    let k = 0.01; // Turbulent kinetic energy [m²/s²]
    let omega = 100.0; // Specific dissipation rate [1/s]
    let _density = 1.0; // Density [kg/m³]

    let nu_t = model.turbulent_viscosity(k, omega, _density);

    // Test strain rate limiter (important for SST)
    let strain_rate_magnitude = 50.0; // High strain rate [1/s]
    let f2 = 0.5; // Mid-range blending function
    let nu_t_limited =
        model.turbulent_viscosity_with_limiter(k, omega, _density, strain_rate_magnitude, f2);

    let p_k = model.production_term(&[[0.0, 1.0], [0.0, 0.0]], nu_t);
    let d_k = model.dissipation_term(k, omega);

    println!("Input conditions:");
    println!("  k = {:.4} m²/s² (turbulent kinetic energy)", k);
    println!("  ω = {:.1} 1/s (specific dissipation rate)", omega);
    println!("Model predictions:");
    println!("  ν_t = {:.6} m²/s (turbulent viscosity)", nu_t);
    println!(
        "  ν_t_limited = {:.6} m²/s (with strain limiter)",
        nu_t_limited
    );
    println!("  P_k = {:.6} m²/s³ (production)", p_k);
    println!("  D_k = {:.6} m²/s³ (destruction)", d_k);
    println!("Characteristics:");
    println!("  - Two-equation model (k, ω)");
    println!("  - SST blending for boundary layers");
    println!("  - Bradshaw limiter prevents excessive ν_t");
    println!("  - Better near-wall behavior than k-ε");

    Ok(())
}

fn demonstrate_spalart_allmaras_model(
    _grid: &StructuredGrid2D<f64>,
    _fields: &SimulationFields<f64>,
) -> Result<(), Box<dyn std::error::Error>> {
    println!("\n--- Spalart-Allmaras Turbulence Model ---");

    let model = SpalartAllmaras::<f64>::new(20, 20);

    // For SA model: ν̃ is the modified turbulent viscosity
    let nu_tilde = 1e-4; // Modified viscosity [m²/s]
    let molecular_viscosity = 1e-5; // Molecular viscosity [m²/s]
    let _density = 1.0; // Density [kg/m³]

    let nu_t = model.eddy_viscosity(nu_tilde, molecular_viscosity);

    // Test production and destruction terms
    let vorticity_magnitude = 100.0; // Vorticity [1/s]
    let wall_distance = 0.01; // Distance to wall [m]

    let production = model.production(nu_tilde, vorticity_magnitude);
    // For demonstration, use a typical fw value (wall destruction function)
    let fw = 0.5; // Typical value for attached boundary layers
    let destruction = model.destruction(nu_tilde, wall_distance, fw);

    println!("Input conditions:");
    println!("  ν̃ = {:.2e} m²/s (modified turbulent viscosity)", nu_tilde);
    println!(
        "  ν = {:.2e} m²/s (molecular viscosity)",
        molecular_viscosity
    );
    println!("  Ω = {:.1} 1/s (vorticity magnitude)", vorticity_magnitude);
    println!("  d = {:.3} m (wall distance)", wall_distance);
    println!("Model predictions:");
    println!("  ν_t = {:.2e} m²/s (turbulent viscosity)", nu_t);
    println!("  P = {:.2e} m²/s³ (production)", production);
    println!("  D = {:.2e} m²/s³ (destruction)", destruction);
    println!("  fw = {:.4} (wall destruction function)", fw);
    println!("Characteristics:");
    println!("  - One-equation model (ν̃ transport)");
    println!("  - Excellent for aerospace applications");
    println!("  - Robust wall boundary conditions");
    println!("  - Simpler than two-equation models");

    Ok(())
}

fn compare_model_predictions(
    _grid: &StructuredGrid2D<f64>,
    _fields: &SimulationFields<f64>,
) -> Result<(), Box<dyn std::error::Error>> {
    println!("\n--- Model Comparison ---");

    // Common test conditions
    let k = 0.01;
    let density = 1.0;
    let strain_rate = 50.0;

    // k-ε model
    let k_eps_model = KEpsilonModel::new(10, 10);
    let epsilon = k_eps_model.dissipation_term(k, 0.0) * 0.1; // Approximate ε
    let nu_t_k_eps = k_eps_model.turbulent_viscosity(k, epsilon, density);

    // k-ω SST model
    let k_omega_model = KOmegaSSTModel::new(10, 10);
    let omega = epsilon / (0.09 * k); // ω = ε / (C_μ * k)
    let nu_t_k_omega = k_omega_model.turbulent_viscosity(k, omega, density);
    let nu_t_k_omega_limited =
        k_omega_model.turbulent_viscosity_with_limiter(k, omega, density, strain_rate, 0.5);

    // SA model (approximate ν̃ from ν_t)
    let sa_model = SpalartAllmaras::<f64>::new(10, 10);
    let molecular_viscosity = 1e-5;
    // ν̃ ≈ ν_t / fv1, but for comparison we'll use a typical value
    let nu_tilde = 1e-4;
    let nu_t_sa = sa_model.eddy_viscosity(nu_tilde, molecular_viscosity);

    println!("Comparison for moderate turbulence conditions:");
    println!(
        "  k = {:.2e} m²/s², strain rate = {:.1} 1/s",
        k, strain_rate
    );
    println!();
    println!("  k-ε model:         ν_t = {:.2e} m²/s", nu_t_k_eps);
    println!("  k-ω SST model:     ν_t = {:.2e} m²/s", nu_t_k_omega);
    println!(
        "  k-ω SST limited:   ν_t = {:.2e} m²/s",
        nu_t_k_omega_limited
    );
    println!("  SA model:          ν_t = {:.2e} m²/s", nu_t_sa);
    println!();
    println!("Key differences:");
    println!("  - k-ε: Good for free shear flows, may need wall functions");
    println!("  - k-ω SST: Excellent near walls, Bradshaw limiter for stability");
    println!("  - SA: Simple one-equation, robust for complex geometries");

    Ok(())
}
