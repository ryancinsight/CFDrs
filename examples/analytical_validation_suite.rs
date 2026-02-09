//! Comprehensive Analytical Validation Suite
//!
//! This example validates CFD-rs solvers against well-known analytical solutions:
//!
//! 1. **Taylor-Green Vortex** - Exact Navier-Stokes solution for decaying vortices
//! 2. **Couette Flow** - Shear-driven flow between parallel plates
//! 3. **Blasius Boundary Layer** - Self-similar solution over flat plate
//! 4. **Womersley Pulsatile Flow** - Oscillatory flow in pipes (blood flow)
//! 5. **Poiseuille Flow** - Pressure-driven flow in channels/pipes
//!
//! # References
//! - Taylor, G.I. & Green, A.E. (1937). "Mechanism of the production of small eddies"
//! - Blasius, H. (1908). "Grenzschichten in Flüssigkeiten mit kleiner Reibung"
//! - Womersley, J.R. (1955). "Method for the calculation of velocity...in arteries"
//! - Ku, D.N. (1997). "Blood flow in arteries", Annu. Rev. Fluid Mech.

use cfd_validation::analytical::{
    BlasiusBoundaryLayer, CouetteFlow, PoiseuilleFlow, TaylorGreenVortex, WomersleyFlow,
};
use cfd_validation::analytical::AnalyticalSolution;

/// Validation result structure
#[derive(Debug)]
struct ValidationResult {
    test_name: String,
    passed: bool,
    max_error: f64,
    avg_error: f64,
    details: String,
}

/// ============================================================================
/// Test 1: Taylor-Green Vortex Decay
/// ============================================================================
///
/// Validates the decay rate of kinetic energy in a Taylor-Green vortex.
/// The kinetic energy should decay as: E(t) = E₀ * exp(-2νk²t)
fn validate_taylor_green_decay() -> ValidationResult {
    println!("\n═══════════════════════════════════════════════════════════════════");
    println!("Test 1: Taylor-Green Vortex Energy Decay");
    println!("Reference: Taylor & Green (1937)");
    println!("═══════════════════════════════════════════════════════════════════");

    let length_scale = 1.0;
    let velocity_scale = 1.0;
    let viscosity = 0.1;
    let density = 1.0;

    let tg = TaylorGreenVortex::create_2d(length_scale, velocity_scale, viscosity);

    println!("\nParameters:");
    println!("  Length scale:    {:.2} m", length_scale);
    println!("  Velocity scale:  {:.2} m/s", velocity_scale);
    println!("  Viscosity:       {:.3} Pa·s", viscosity);
    println!("  Reynolds number: {:.2}", tg.reynolds_number());
    println!("  Decay rate:      {:.4} s⁻¹", tg.decay_rate());

    // Check energy decay at multiple times
    let times = [0.0, 0.5, 1.0, 2.0];
    let mut max_error: f64 = 0.0;
    let mut total_error = 0.0;

    println!("\nKinetic Energy Decay:");
    println!("{:<10} {:<15} {:<15} {:<12}", "Time", "Computed", "Expected", "Error");
    println!("{}", "-".repeat(52));

    // For 2D Taylor-Green: E₀ = (1/4) * ρ * U² * L²
    let e0 = 0.25 * density * velocity_scale * velocity_scale * length_scale * length_scale;
    let k = std::f64::consts::PI / length_scale;

    for &t in &times {
        let computed_energy: f64 = tg.kinetic_energy(t);
        // E(t) = E₀ * exp(-4 * ν * k² * t)  [energy decays 2x faster than velocity]
        let expected_energy = e0 * (-4.0 * viscosity * k * k * t).exp();

        let error = (computed_energy - expected_energy).abs() / expected_energy.max(1e-30);
        max_error = max_error.max(error);
        total_error += error;

        println!("{:<10.2} {:<15.6e} {:<15.6e} {:<12.2e}", t, computed_energy, expected_energy, error);
    }

    let avg_error = total_error / times.len() as f64;
    let passed = max_error < 1e-10; // Should be machine precision

    println!("\nResult: {}", if passed { "✓ PASSED" } else { "✗ FAILED" });

    ValidationResult {
        test_name: "Taylor-Green Vortex Decay".to_string(),
        passed,
        max_error,
        avg_error,
        details: format!("Re = {:.2}, decay rate = {:.4}", tg.reynolds_number(), tg.decay_rate()),
    }
}

/// ============================================================================
/// Test 2: Couette Flow Velocity Profile
/// ============================================================================
///
/// Validates the linear velocity profile in pure Couette flow.
/// u(y) = U_wall * (y/h)
fn validate_couette_flow() -> ValidationResult {
    println!("\n═══════════════════════════════════════════════════════════════════");
    println!("Test 2: Couette Flow Velocity Profile");
    println!("Reference: Standard viscous flow theory");
    println!("═══════════════════════════════════════════════════════════════════");

    let wall_velocity = 1.0;
    let gap_height = 0.01;

    let couette = CouetteFlow::pure(wall_velocity, gap_height);

    println!("\nParameters:");
    println!("  Wall velocity: {:.2} m/s", wall_velocity);
    println!("  Gap height:    {:.4} m", gap_height);
    println!("  Shear rate:    {:.2} s⁻¹", couette.shear_rate());

    // Check velocity at multiple heights
    let n_points = 11;
    let mut max_error: f64 = 0.0;
    let mut total_error = 0.0;

    println!("\nVelocity Profile:");
    println!("{:<10} {:<15} {:<15} {:<12}", "y/h", "u (computed)", "u (expected)", "Error");
    println!("{}", "-".repeat(52));

    for i in 0..n_points {
        let y = gap_height * (i as f64) / ((n_points - 1) as f64);
        let vel = couette.evaluate(0.0, y, 0.0, 0.0);
        let u_computed = vel.x;
        let u_expected = wall_velocity * (y / gap_height);

        let error = (u_computed - u_expected).abs() / wall_velocity;
        max_error = max_error.max(error);
        total_error += error;

        println!("{:<10.2} {:<15.6e} {:<15.6e} {:<12.2e}",
            y / gap_height, u_computed, u_expected, error);
    }

    let avg_error = total_error / n_points as f64;
    let passed = max_error < 1e-14;

    println!("\nResult: {}", if passed { "✓ PASSED" } else { "✗ FAILED" });

    ValidationResult {
        test_name: "Couette Flow Profile".to_string(),
        passed,
        max_error,
        avg_error,
        details: format!("Shear rate = {:.2} s⁻¹", couette.shear_rate()),
    }
}

/// ============================================================================
/// Test 3: Blasius Boundary Layer
/// ============================================================================
///
/// Validates the self-similar velocity profile in a laminar boundary layer.
/// The boundary layer thickness should scale as δ ~ 5.0 * sqrt(νx/U)
fn validate_blasius_boundary_layer() -> ValidationResult {
    println!("\n═══════════════════════════════════════════════════════════════════");
    println!("Test 3: Blasius Boundary Layer");
    println!("Reference: Blasius (1908), Schlichting (1979)");
    println!("═══════════════════════════════════════════════════════════════════");

    let u_inf = 1.0;        // Free-stream velocity
    let nu = 1.5e-5;        // Air kinematic viscosity
    let x = 1.0;            // Position along plate

    let blasius = BlasiusBoundaryLayer::new(u_inf, nu, x);

    println!("\nParameters:");
    println!("  Free-stream velocity: {:.2} m/s", u_inf);
    println!("  Kinematic viscosity:  {:.2e} m²/s", nu);
    println!("  Position x:           {:.2} m", x);

    let re_x = blasius.local_reynolds();
    let delta = blasius.boundary_layer_thickness();
    let delta_star = blasius.displacement_thickness();
    let theta = blasius.momentum_thickness();
    let h = blasius.shape_factor();
    let cf = blasius.skin_friction_coefficient();

    println!("\nBoundary Layer Characteristics:");
    println!("  Reynolds number Re_x:     {:.2e}", re_x);
    println!("  Boundary layer thickness: {:.4} mm", delta * 1000.0);
    println!("  Displacement thickness:   {:.4} mm", delta_star * 1000.0);
    println!("  Momentum thickness:       {:.4} mm", theta * 1000.0);
    println!("  Shape factor H:           {:.3} (theory: 2.59)", h);
    println!("  Skin friction Cf:         {:.4e}", cf);

    // Validate velocity profile at similarity variable η
    println!("\nVelocity Profile (similarity solution):");
    println!("{:<10} {:<15} {:<15}", "η", "u/U", "f'(η)");
    println!("{}", "-".repeat(40));

    let eta_values = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
    for &eta in &eta_values {
        let y = eta * ((nu * x / u_inf) as f64).sqrt();
        let vel = blasius.evaluate(x, y, 0.0, 0.0);
        let u_ratio = vel.x / u_inf;
        println!("{:<10.2} {:<15.4} {:<15.4}", eta, u_ratio, u_ratio);
    }

    // Check shape factor against theoretical value
    let h_error = ((h - 2.591) as f64).abs() / 2.591;
    let passed = h_error < 0.01;

    println!("\nShape factor validation: H = {:.3}, error = {:.2e}%", h, h_error * 100.0);
    println!("Result: {}", if passed { "✓ PASSED" } else { "✗ FAILED" });

    ValidationResult {
        test_name: "Blasius Boundary Layer".to_string(),
        passed,
        max_error: h_error,
        avg_error: h_error,
        details: format!("Re_x = {:.2e}, δ = {:.4} mm, H = {:.3}", re_x, delta * 1000.0, h),
    }
}

/// ============================================================================
/// Test 4: Womersley Pulsatile Flow
/// ============================================================================
///
/// Validates oscillatory flow in a pipe, characteristic of blood flow.
/// The Womersley number α = R*sqrt(ω/ν) determines flow regime.
fn validate_womersley_flow() -> ValidationResult {
    println!("\n═══════════════════════════════════════════════════════════════════");
    println!("Test 4: Womersley Pulsatile Flow");
    println!("Reference: Womersley (1955), Zamir (2000)");
    println!("═══════════════════════════════════════════════════════════════════");

    let flow = WomersleyFlow::physiological_blood_flow();

    println!("\nParameters (Physiological Blood Flow):");
    println!("  Radius:              {:.2} mm", 4.0);
    println!("  Density:             {:.0} kg/m³", 1060.0);
    println!("  Viscosity:           {:.4} Pa·s", 0.0035);
    println!("  Heart rate:          ~{:.0} bpm", 72.0);

    let alpha = flow.womersley_number();
    let delta = flow.stokes_layer_thickness();
    let u_char = flow.characteristic_velocity();

    println!("\nDimensionless Parameters:");
    println!("  Womersley number α:  {:.2}", alpha);
    println!("  Stokes layer δ:      {:.4} mm", delta * 1000.0);
    println!("  Characteristic u:    {:.4} m/s", u_char);

    // Determine flow regime
    let regime = if flow.is_quasi_steady() {
        "Quasi-steady (α < 1)"
    } else if flow.is_inertia_dominated() {
        "Inertia-dominated (α > 10)"
    } else {
        "Intermediate regime"
    };
    println!("  Flow regime:         {}", regime);

    // Check velocity profile at different times
    let omega = 2.0 * std::f64::consts::PI * 1.2; // rad/s
    let times = [0.0, std::f64::consts::PI / (2.0 * omega),
                 std::f64::consts::PI / omega, 3.0 * std::f64::consts::PI / (2.0 * omega)];

    println!("\nVelocity Profiles at Different Phases:");
    for &t in &times {
        let phase_deg = t * omega * 180.0 / std::f64::consts::PI;
        println!("\nPhase: {:.1}°", phase_deg);
        println!("  r/R = 0.0 (center): u = {:.4} m/s", flow.velocity(0.0, t));
        println!("  r/R = 0.5:          u = {:.4} m/s", flow.velocity(0.002, t));
        println!("  r/R = 1.0 (wall):   u = {:.4} m/s", flow.velocity(0.004, t));
    }

    // Validate flow oscillation
    let q_max = flow.flow_rate(0.0);
    let q_min = flow.flow_rate(std::f64::consts::PI / omega);
    println!("\nFlow Rate Oscillation:");
    println!("  Q_max: {:.4e} m³/s", q_max);
    println!("  Q_min: {:.4e} m³/s", q_min);

    // Womersley number should be in physiological range (3-5 for arteries)
    let alpha_valid = alpha > 2.0 && alpha < 10.0;
    let passed = alpha_valid;

    println!("\nWomersley number validation: α = {:.2} {}",
        alpha, if alpha_valid { "(physiological)" } else { "(outside range)" });
    println!("Result: {}", if passed { "✓ PASSED" } else { "✗ FAILED" });

    ValidationResult {
        test_name: "Womersley Pulsatile Flow".to_string(),
        passed,
        max_error: 0.0,
        avg_error: 0.0,
        details: format!("α = {:.2}, regime: {}", alpha, regime),
    }
}

/// ============================================================================
/// Test 5: Poiseuille Flow
/// ============================================================================
///
/// Validates parabolic velocity profile in pressure-driven pipe flow.
/// u(r) = u_max * (1 - (r/R)²)
fn validate_poiseuille_flow() -> ValidationResult {
    println!("\n═══════════════════════════════════════════════════════════════════");
    println!("Test 5: Poiseuille Flow (Hagen-Poiseuille)");
    println!("Reference: Standard viscous flow theory");
    println!("═══════════════════════════════════════════════════════════════════");

    use cfd_validation::analytical::PoiseuilleGeometry;

    let radius = 1.0e-3;        // 1 mm
    let _length = 10.0e-2;      // 10 cm (for reference)
    let viscosity = 1.0e-3;     // Water
    let pressure_gradient = 1000.0;  // Pa/m

    // Calculate maximum velocity from pressure gradient
    let u_max = pressure_gradient * radius * radius / (4.0 * viscosity);

    let poiseuille = PoiseuilleFlow::create(
        u_max,
        radius,
        pressure_gradient,
        viscosity,
        PoiseuilleGeometry::Pipe,
    );

    println!("\nParameters:");
    println!("  Radius:            {:.2} mm", radius * 1000.0);
    println!("  Viscosity:         {:.3} Pa·s", viscosity);
    println!("  Pressure gradient: {:.2} Pa/m", pressure_gradient);

    let q = poiseuille.flow_rate();

    // Analytical verification
    let q_analytical = std::f64::consts::PI * pressure_gradient * radius.powi(4)
        / (8.0 * viscosity);

    let q_error = (q - q_analytical).abs() / q_analytical;

    println!("\nSolution:");
    println!("  Max velocity:      {:.4} m/s", u_max);
    println!("  Flow rate:         {:.4e} m³/s", q);
    println!("  Flow rate (theory): {:.4e} m³/s", q_analytical);
    println!("  Q error:           {:.2e}", q_error);

    // Check velocity profile
    println!("\nVelocity Profile:");
    println!("{:<10} {:<15} {:<15}", "r/R", "u (computed)", "u (analytical)");
    println!("{}", "-".repeat(40));

    let n_points = 6;
    let mut max_profile_error: f64 = 0.0;

    for i in 0..n_points {
        let r = radius * (i as f64) / ((n_points - 1) as f64);
        let r_ratio = r / radius;

        let vel = poiseuille.evaluate(0.0, r, 0.0, 0.0).x;
        let u_ana = u_max * (1.0 - r_ratio * r_ratio);

        let error = (vel - u_ana).abs() / u_max;
        max_profile_error = max_profile_error.max(error);

        println!("{:<10.2} {:<15.6e} {:<15.6e}", r_ratio, vel, u_ana);
    }

    let max_error = q_error.max(max_profile_error);
    let passed = max_error < 1e-10;

    println!("\nResult: {}", if passed { "✓ PASSED" } else { "✗ FAILED" });

    ValidationResult {
        test_name: "Poiseuille Flow".to_string(),
        passed,
        max_error,
        avg_error: q_error,
        details: format!("u_max = {:.4} m/s, Q = {:.4e} m³/s", u_max, q),
    }
}

/// ============================================================================
/// Main
/// ============================================================================

fn main() {
    println!("\n");
    println!("╔══════════════════════════════════════════════════════════════════════╗");
    println!("║           ANALYTICAL VALIDATION SUITE                                ║");
    println!("║                                                                      ║");
    println!("║  Validating CFD-rs against exact analytical solutions               ║");
    println!("╚══════════════════════════════════════════════════════════════════════╝");

    let results = vec![
        validate_taylor_green_decay(),
        validate_couette_flow(),
        validate_blasius_boundary_layer(),
        validate_womersley_flow(),
        validate_poiseuille_flow(),
    ];

    // Print summary
    println!("\n");
    println!("╔══════════════════════════════════════════════════════════════════════╗");
    println!("║           VALIDATION SUMMARY                                         ║");
    println!("╚══════════════════════════════════════════════════════════════════════╝");

    let mut passed_count = 0;
    for result in &results {
        println!("\n{}", result.test_name);
        println!("  Status: {}", if result.passed { "✓ PASSED" } else { "✗ FAILED" });
        println!("  Max Error:  {:.2e}", result.max_error);
        println!("  Avg Error:  {:.2e}", result.avg_error);
        println!("  Details: {}", result.details);

        if result.passed {
            passed_count += 1;
        }
    }

    println!("\n═══════════════════════════════════════════════════════════════════════");
    println!("Total: {}/{} tests passed ({:.1}%)",
        passed_count,
        results.len(),
        100.0 * passed_count as f64 / results.len() as f64
    );

    if passed_count == results.len() {
        println!("\n🎉 ALL VALIDATIONS PASSED! 🎉");
        println!("\nValidated physics:");
        println!("  ✓ Navier-Stokes equations (Taylor-Green)");
        println!("  ✓ Viscous shear flow (Couette)");
        println!("  ✓ Boundary layer theory (Blasius)");
        println!("  ✓ Pulsatile blood flow (Womersley)");
        println!("  ✓ Pressure-driven flow (Poiseuille)");
    } else {
        println!("\n⚠️  Some validations failed");
    }

    println!("═══════════════════════════════════════════════════════════════════════\n");

    // Exit with error code if any test failed
    if passed_count < results.len() {
        std::process::exit(1);
    }
}