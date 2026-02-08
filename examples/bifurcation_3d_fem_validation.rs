//! 3D Bifurcation FEM Solver Validation
//!
//! This example validates the 3D FEM bifurcation solver against analytical solutions
//! and literature data for blood flow in bifurcating vessels.
//!
//! # Validation Strategy
//!
//! ## 1. Mass Conservation (Fundamental Check)
//! ```text
//! ∫_Ω ∇·u dV = Q_in - Q_out1 - Q_out2 = 0
//! ```
//! For incompressible flow, mass must be conserved exactly.
//!
//! ## 2. Poiseuille Flow Limit (Parent Branch)
//! Far from the bifurcation, the parent branch should exhibit Poiseuille flow:
//! ```text
//! u(r) = u_max (1 - r²/R²)
//! u_max = 2 u_mean (parabolic profile)
//! ```
//!
//! ## 3. Murray's Law Verification
//! For optimal bifurcation (minimum energy dissipation):
//! ```text
//! D₀³ = D₁³ + D₂³
//! ```
//! We verify the flow solver respects this relationship.
//!
//! ## 4. Wall Shear Stress Scaling
//! For Poiseuille flow in cylinder:
//! ```text
//! τ_w = 8μu_mean/D = 4μu_max/D
//! ```
//! Daughters should have higher WSS due to smaller diameter.
//!
//! ## 5. Pressure Drop Validation
//! Using Hagen-Poiseuille for each segment:
//! ```text
//! ΔP = 128μLQ/(πD⁴)
//! ```
//!
//! # References
//!
//! - Ku et al. (1985): Pulsatile flow in human carotid bifurcation
//! - Glagov et al. (1988): Hemodynamics and atherosclerosis
//! - Murray (1926): The physiological principle of minimum work

use cfd_3d::bifurcation::{
    BifurcationConfig3D, BifurcationGeometry3D, BifurcationSolver3D,
    solver::BifurcationSolution3D,
};
use cfd_core::physics::fluid::blood::CassonBlood;
use cfd_core::physics::fluid::{ConstantPropertyFluid, NonNewtonianFluid, FluidTrait};


/// Validation test result
struct ValidationResult {
    name: String,
    passed: bool,
    expected: f64,
    computed: f64,
    relative_error: f64,
    tolerance: f64,
}

impl ValidationResult {
    fn print(&self) {
        let status = if self.passed { "✓ PASS" } else { "✗ FAIL" };
        println!(
            "  {}: {} (expected {:.4e}, got {:.4e}, error {:.2}%, tol {:.2}%)",
            status, self.name, self.expected, self.computed, self.relative_error * 100.0, self.tolerance * 100.0
        );
    }
}

/// Test 1: Mass Conservation
/// 
/// For incompressible flow: Q_in = Q_out1 + Q_out2
fn test_mass_conservation(solution: &BifurcationSolution3D<f64>, tolerance: f64) -> ValidationResult {
    let q_in = solution.q_parent;
    let q_out = solution.q_daughter1 + solution.q_daughter2;
    let error = (q_in - q_out).abs() / q_in;
    
    ValidationResult {
        name: "Mass Conservation".to_string(),
        passed: error < tolerance,
        expected: q_in,
        computed: q_out,
        relative_error: error,
        tolerance,
    }
}

/// Test 2: Flow Split Ratio
///
/// For symmetric bifurcation: Q1/(Q1+Q2) ≈ 0.5
fn test_flow_split(solution: &BifurcationSolution3D<f64>, expected_ratio: f64, tolerance: f64) -> ValidationResult {
    let q_total = solution.q_daughter1 + solution.q_daughter2;
    let split = solution.q_daughter1 / q_total;
    let error = (split - expected_ratio).abs();
    
    ValidationResult {
        name: "Flow Split Ratio".to_string(),
        passed: error < tolerance,
        expected: expected_ratio,
        computed: split,
        relative_error: error,
        tolerance,
    }
}

/// Test 3: Wall Shear Stress Scaling
///
/// For Poiseuille flow: τ_w = 8μu_mean/D
/// WSS ratio: τ_w1/τ_w0 = (D0/D1) * (u1/u0) = (D0/D1)³ (for equal pressure drops)
fn test_wss_scaling(
    solution: &BifurcationSolution3D<f64>,
    d_parent: f64,
    d_daughter: f64,
    tolerance: f64,
) -> ValidationResult {
    // For Murray's law bifurcation with symmetric daughters
    // Expected WSS ratio: τ_daughter/τ_parent = (D_parent/D_daughter)³
    let expected_ratio = (d_parent / d_daughter).powi(3);
    let actual_ratio = solution.wall_shear_stress_daughter1 / solution.wall_shear_stress_parent;
    let error = (actual_ratio - expected_ratio).abs() / expected_ratio;
    
    ValidationResult {
        name: "WSS Scaling (Murray's Law)".to_string(),
        passed: error < tolerance,
        expected: expected_ratio,
        computed: actual_ratio,
        relative_error: error,
        tolerance,
    }
}

/// Test 4: Poiseuille WSS Formula
///
/// Verify: τ_w = 8μu_mean/D for parent branch
fn test_poiseuille_wss(
    solution: &BifurcationSolution3D<f64>,
    viscosity: f64,
    d_parent: f64,
    tolerance: f64,
) -> ValidationResult {
    let expected_wss = 8.0 * viscosity * solution.u_parent_mean / d_parent;
    let actual_wss = solution.wall_shear_stress_parent;
    let error = (actual_wss - expected_wss).abs() / expected_wss;
    
    ValidationResult {
        name: "Poiseuille WSS Formula".to_string(),
        passed: error < tolerance,
        expected: expected_wss,
        computed: actual_wss,
        relative_error: error,
        tolerance,
    }
}

/// Test 5: Pressure Drop Direction
///
/// Pressure should decrease in flow direction: P_in > P_junction > P_out
fn test_pressure_drop(solution: &BifurcationSolution3D<f64>) -> ValidationResult {
    let p_drop_parent = solution.p_inlet - solution.p_junction_mid;
    let p_drop_d1 = solution.p_junction_mid - solution.p_daughter1_outlet;
    
    let parent_positive = p_drop_parent > 0.0;
    let daughter_positive = p_drop_d1 > 0.0;
    let passed = parent_positive && daughter_positive;
    
    ValidationResult {
        name: "Pressure Drop Direction".to_string(),
        passed,
        expected: 1.0, // Positive pressure drop
        computed: if passed { 1.0 } else { 0.0 },
        relative_error: if passed { 0.0 } else { 1.0 },
        tolerance: 0.0,
    }
}

/// Test 6: Reynolds Number Check
///
/// Verify flow is laminar: Re < 2300
fn test_laminar_flow(solver: &BifurcationSolver3D<f64>, fluid: &CassonBlood<f64>) -> ValidationResult {
    let re = solver.reynolds_number(fluid.clone()).unwrap();
    let is_laminar = re < 2300.0;
    
    ValidationResult {
        name: "Laminar Flow (Re < 2300)".to_string(),
        passed: is_laminar,
        expected: 1000.0, // Typical value
        computed: re,
        relative_error: if is_laminar { 0.0 } else { (re - 2300.0) / 2300.0 },
        tolerance: 1.0,
    }
}

/// Run all validation tests for symmetric bifurcation
fn validate_symmetric_bifurcation() -> Vec<ValidationResult> {
    println!("\n{}", "=".repeat(80));
    println!("VALIDATION: Symmetric 3D Bifurcation (100 μm → 79.4 μm each)");
    println!("{}\n", "=".repeat(80));
    
    // Geometry following Murray's law: D₀³ = 2 × D₁³ → D₁ = D₀/2^(1/3)
    let d_parent: f64 = 100e-6; // 100 μm
    let d_daughter: f64 = d_parent / 2.0_f64.powf(1.0/3.0); // Murray's law
    let length: f64 = 1e-3; // 1 mm
    
    println!("Geometry:");
    println!("  Parent diameter: {:.1} μm", d_parent * 1e6);
    println!("  Daughter diameter: {:.1} μm (Murray's law)", d_daughter * 1e6);
    println!("  Vessel length: {:.1} mm", length * 1e3);
    
    // Flow conditions (blood)
    let flow_rate = 1e-9; // 1 μL/s
    let inlet_pressure = 100.0; // Pa
    
    println!("\nFlow Conditions:");
    println!("  Flow rate: {:.1} μL/s", flow_rate * 1e9);
    println!("  Inlet pressure: {:.1} Pa", inlet_pressure);
    
    // Create solver
    let geometry = BifurcationGeometry3D::symmetric(
        d_parent, d_daughter, length, length, 100e-6,
    );
    
    let config = BifurcationConfig3D {
        inlet_flow_rate: flow_rate,
        inlet_pressure,
        outlet_pressure: 0.0,
        ..Default::default()
    };
    
    let solver = BifurcationSolver3D::new(geometry, config);
    
    // Use Casson blood model
    let blood = CassonBlood::<f64>::normal_blood();
    let blood_props = blood.properties_at(310.0, inlet_pressure).unwrap();
    
    println!("\nFluid Properties (Casson Blood):");
    println!("  Density: {:.0} kg/m³", blood_props.density);
    println!("  Viscosity: {:.4} Pa·s", blood_props.dynamic_viscosity);
    println!("  Yield stress: {:.4} Pa", blood.yield_stress().unwrap_or(0.0));
    
    // Solve
    // Debug: check what boundary labels exist in mesh
    println!("\n[Debug] Checking mesh boundary labels...");
    let mesh_builder = cfd_mesh::geometry::branching::BranchingMeshBuilder::bifurcation(
        d_parent, length, d_daughter, length, std::f64::consts::PI / 6.0, 8
    );
    let debug_mesh = mesh_builder.build().expect("Mesh build failed");
    println!("  Total faces: {}", debug_mesh.face_count());
    println!("  Boundary faces: {}", debug_mesh.boundary_faces().len());
    
    // Show inlet/outlet face positions
    for f_idx in debug_mesh.boundary_faces() {
        if let Some(label) = debug_mesh.boundary_label(f_idx) {
            if let Some(face) = debug_mesh.face(f_idx) {
                // Compute face centroid
                let mut centroid = nalgebra::Point3::origin();
                for &v_idx in &face.vertices {
                    if let Some(v) = debug_mesh.vertex(v_idx) {
                        centroid += v.position.coords;
                    }
                }
                centroid /= face.vertices.len() as f64;
                println!("    Face {}: label = '{:10}' at ({:.3e}, {:.3e}, {:.3e})", 
                         f_idx, label, centroid.x, centroid.y, centroid.z);
            }
        }
    }
    
    println!("\nSolving...");
    let solution = solver.solve(blood.clone()).expect("Solver failed");
    
    println!("\nSolution Summary:");
    println!("  Parent flow: {:.3e} m³/s", solution.q_parent);
    println!("  Daughter 1 flow: {:.3e} m³/s", solution.q_daughter1);
    println!("  Daughter 2 flow: {:.3e} m³/s", solution.q_daughter2);
    println!("  Mass conservation error: {:.2e}", solution.mass_conservation_error);
    println!("  Parent WSS: {:.3} Pa", solution.wall_shear_stress_parent);
    println!("  Daughter WSS: {:.3} Pa", solution.wall_shear_stress_daughter1);
    
    // Run validation tests
    let mut results = Vec::new();
    
    results.push(test_mass_conservation(&solution, 1e-6));
    results.push(test_flow_split(&solution, 0.5, 0.01));
    results.push(test_wss_scaling(&solution, d_parent, d_daughter, 0.20));
    results.push(test_poiseuille_wss(&solution, blood_props.dynamic_viscosity, d_parent, 0.10));
    results.push(test_pressure_drop(&solution));
    results.push(test_laminar_flow(&solver, &blood));
    
    results
}

/// Run all validation tests for asymmetric bifurcation
fn validate_asymmetric_bifurcation() -> Vec<ValidationResult> {
    println!("\n{}", "=".repeat(80));
    println!("VALIDATION: Asymmetric 3D Bifurcation (100 μm → 60 μm + 80 μm)");
    println!("{}\n", "=".repeat(80));
    
    // Asymmetric geometry (typical of arterial bifurcations)
    let d_parent: f64 = 100e-6;
    let d_daughter1: f64 = 60e-6; // Smaller branch
    let d_daughter2: f64 = 80e-6; // Larger branch
    let length: f64 = 1e-3;
    
    println!("Geometry:");
    println!("  Parent diameter: {:.1} μm", d_parent * 1e6);
    println!("  Daughter 1 diameter: {:.1} μm (smaller)", d_daughter1 * 1e6);
    println!("  Daughter 2 diameter: {:.1} μm (larger)", d_daughter2 * 1e6);
    
    // Check Murray's law deviation
    let murray_deviation = (d_parent.powi(3) - (d_daughter1.powi(3) + d_daughter2.powi(3))).abs() 
                          / d_parent.powi(3);
    println!("  Murray's law deviation: {:.1}%", murray_deviation * 100.0);
    
    // Flow conditions
    let flow_rate = 1e-9;
    let inlet_pressure = 100.0;
    
    // Create solver
    let geometry = BifurcationGeometry3D::asymmetric(
        d_parent, d_daughter1, d_daughter2,
        length, length, length, 100e-6,
        35.0_f64.to_radians(),
    );
    
    let config = BifurcationConfig3D {
        inlet_flow_rate: flow_rate,
        inlet_pressure,
        outlet_pressure: 0.0,
        ..Default::default()
    };
    
    let solver = BifurcationSolver3D::new(geometry, config);
    
    // Use water for Newtonian comparison
    let water = ConstantPropertyFluid {
        name: "Water".to_string(),
        density: 1000.0,
        viscosity: 0.001,
        specific_heat: 4186.0,
        thermal_conductivity: 0.6,
        speed_of_sound: 1500.0,
    };
    
    println!("\nSolving...");
    let solution = solver.solve(water).expect("Solver failed");
    
    println!("\nSolution Summary:");
    println!("  Parent flow: {:.3e} m³/s", solution.q_parent);
    println!("  Daughter 1 flow: {:.3e} m³/s", solution.q_daughter1);
    println!("  Daughter 2 flow: {:.3e} m³/s", solution.q_daughter2);
    println!("  Flow split: {:.2}% / {:.2}%", 
             100.0 * solution.q_daughter1 / (solution.q_daughter1 + solution.q_daughter2),
             100.0 * solution.q_daughter2 / (solution.q_daughter1 + solution.q_daughter2));
    
    // Expected flow split based on resistance (R ~ 1/D⁴)
    // Q1/Q2 = (D1/D2)⁴ for equal pressure drops
    let expected_split = d_daughter1.powi(4) / (d_daughter1.powi(4) + d_daughter2.powi(4));
    
    let mut results = Vec::new();
    results.push(test_mass_conservation(&solution, 1e-6));
    results.push(test_flow_split(&solution, expected_split, 0.05));
    results.push(test_pressure_drop(&solution));
    
    results
}

fn main() {
    println!("\n{}", "█".repeat(80));
    println!("█{:^78}█", "3D BIFURCATION FEM SOLVER VALIDATION");
    println!("{}", "█".repeat(80));
    
    let start_time = std::time::Instant::now();
    
    // Run validation suites
    let symmetric_results = validate_symmetric_bifurcation();
    let asymmetric_results = validate_asymmetric_bifurcation();
    
    // Combine results
    let all_results: Vec<_> = symmetric_results.into_iter()
        .chain(asymmetric_results.into_iter())
        .collect();
    
    // Print summary
    println!("\n{}", "=".repeat(80));
    println!("VALIDATION SUMMARY");
    println!("{}\n", "=".repeat(80));
    
    let passed = all_results.iter().filter(|r| r.passed).count();
    let total = all_results.len();
    
    for result in &all_results {
        result.print();
    }
    
    let elapsed = start_time.elapsed();
    
    println!("\n{}", "-".repeat(80));
    println!("Results: {}/{} tests passed ({:.1}%)", passed, total, 100.0 * passed as f64 / total as f64);
    println!("Elapsed time: {:.2} seconds", elapsed.as_secs_f64());
    println!("{}", "-".repeat(80));
    
    if passed == total {
        println!("\n✓ ALL VALIDATION TESTS PASSED");
        std::process::exit(0);
    } else {
        println!("\n✗ SOME VALIDATION TESTS FAILED");
        std::process::exit(1);
    }
}
