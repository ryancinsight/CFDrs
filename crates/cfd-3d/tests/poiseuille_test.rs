use cfd_3d::venturi::{VenturiSolver3D, VenturiConfig3D};
use cfd_mesh::geometry::venturi::VenturiMeshBuilder;
use cfd_core::physics::fluid::NewtonianFluid;
use approx::assert_relative_eq;

#[test]
fn validate_poiseuille_flow() {
    // 1. Geometry: Straight Pipe (L=5D)
    // D = 1mm, L = 10mm
    let d = 1.0e-3;
    let l_seg = 2.0e-3; 
    let l_total = 5.0 * l_seg;

    let builder = VenturiMeshBuilder::new(d, d, l_seg, l_seg, l_seg, l_seg, l_seg);

    // 2. Fluid Properties
    // Water-like: rho=1000, mu=0.001
    let rho = 1000.0;
    let mu = 0.001;
    let fluid = NewtonianFluid::new(rho, mu);

    // 3. Flow Conditions
    // Target Re = 10 (Laminar)
    // Re = rho * u * d / mu => u = Re * mu / (rho * d)
    // u = 10 * 0.001 / (1000 * 0.001) = 0.01 m/s
    let u_avg = 0.01;
    let area = std::f64::consts::PI * (d / 2.0).powi(2);
    let q_in = u_avg * area;

    let config = VenturiConfig3D {
        inlet_flow_rate: q_in,
        inlet_pressure: 0.0, // Should be calculated, but initial guess
        outlet_pressure: 0.0,
        resolution: (30, 8), // Coarse mesh for speed, but sufficient for laminar
        circular: true,
        max_nonlinear_iterations: 20,
        nonlinear_tolerance: 1e-5,
    };

    println!("Poiseuille Test Setup:");
    println!("  Diameter: {:.2e} m", d);
    println!("  Length:   {:.2e} m", l_total);
    println!("  Viscosity: {:.2e} Pa.s", mu);
    println!("  Target U: {:.2e} m/s", u_avg);
    println!("  Flow Rate: {:.2e} m^3/s", q_in);

    // 4. Solve
    let solver = VenturiSolver3D::new(builder, config);
    let solution = solver.solve(fluid).unwrap();

    // 5. Validation using Analytical Solution
    // Pressure Drop: Delta P = (128 * mu * L * Q) / (pi * D^4)
    let dp_analytical = (128.0 * mu * l_total * q_in) / (std::f64::consts::PI * d.powi(4));
    
    // Measured Pressure Drop
    // Note: VenturiSolver reports p_inlet and p_outlet
    let dp_measured = solution.p_inlet - solution.p_outlet;
    
    println!("\nResults:");
    println!("  Analytical DP: {:.4} Pa", dp_analytical);
    println!("  Measured DP:   {:.4} Pa", dp_measured);
    println!("  DP Error:      {:.2}%", (dp_measured - dp_analytical).abs() / dp_analytical * 100.0);

    // Centerline Velocity: u_max = 2 * u_avg
    let u_max_analytical = 2.0 * u_avg;
    
    // We need to extract the centerline velocity from the solution.
    // VenturiSolution3D usually exposes u_throat, but that's average. 
    // We can't easily query node values here without 'fem_solution'.
    // However, VenturiSolution3D typically has 'u_throat' (average velocity at throat).
    // For a straight pipe, u_throat should equal u_avg (continuity).
    
    println!("  Target U_avg:  {:.4e} m/s", u_avg);
    println!("  Measured U_th: {:.4e} m/s", solution.u_throat);
    let u_error = (solution.u_throat - u_avg).abs() / u_avg * 100.0;
    println!("  U_avg Error:   {:.2}%", u_error);

    // Check Pressure Drop (Relaxed tolerance for coarse mesh)
    assert_relative_eq!(dp_measured, dp_analytical, epsilon = dp_analytical * 0.20); 
    
    // Check Continuity (U_throat should match U_in)
    assert_relative_eq!(solution.u_throat, u_avg, epsilon = u_avg * 0.05);

    // Check Mass Conservation
    assert!(solution.mass_conservation_error < 1e-6, "Mass conservation failed: {:.2e}", solution.mass_conservation_error);
}
