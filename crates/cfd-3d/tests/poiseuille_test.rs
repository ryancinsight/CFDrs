use cfd_3d::venturi::{VenturiConfig3D, VenturiSolution3D, VenturiSolver3D};
use cfd_mesh::geometry::venturi::VenturiMeshBuilder;
use cfd_core::physics::fluid::ConstantPropertyFluid;

fn solve_poiseuille(u_avg: f64, resolution: (usize, usize)) -> VenturiSolution3D<f64> {
    // 1. Geometry: Straight Pipe (L=5D)
    // D = 1mm, L = 5mm
    let d = 1.0e-3;
    let l_seg = 1.0e-3;
    let _l_total = 5.0 * l_seg;

    let builder = VenturiMeshBuilder::new(d, d, l_seg, l_seg, l_seg, l_seg, l_seg);

    // 2. Fluid Properties
    // Water-like: rho=1000, mu=0.001
    let rho = 1000.0;
    let mu = 0.001;
    // Use ConstantPropertyFluid with dummy thermal properties
    let fluid = ConstantPropertyFluid::new(
        "Poiseuille Fluid".to_string(),
        rho,
        mu,
        4186.0, // Cp (Water)
        0.6,    // k (Water)
        1500.0, // c (Water)
    );

    // 3. Flow Conditions
    // Target Re = 10 (Laminar)
    // Re = rho * u * d / mu => u = Re * mu / (rho * d)
    // u = 10 * 0.001 / (1000 * 0.001) = 0.01 m/s
    let area = std::f64::consts::PI * (d / 2.0f64).powi(2);
    let q_in = u_avg * area;

    let config = VenturiConfig3D {
        inlet_flow_rate: q_in,
        inlet_pressure: 0.0, // Should be calculated, but initial guess
        outlet_pressure: 0.0,
        resolution,
        circular: true,
        max_nonlinear_iterations: 20,
        nonlinear_tolerance: 1e-5,
    };

    let solver = VenturiSolver3D::new(builder, config);
    let solution = solver.solve(fluid).unwrap();

    solution
}

fn relative_error(measured: f64, analytical: f64) -> f64 {
    (measured - analytical).abs() / analytical.abs()
}

#[test]
fn validate_poiseuille_flow() {
    let resolution = (16, 4);
    let u_avg_low = 0.005;
    let u_avg_high = 0.01;

    let low_solution = solve_poiseuille(u_avg_low, resolution);
    let high_solution = solve_poiseuille(u_avg_high, resolution);

    let dp_low = (low_solution.p_inlet - low_solution.p_outlet).abs();
    let dp_high = (high_solution.p_inlet - high_solution.p_outlet).abs();
    let d = 1.0e-3;
    let l_total = 5.0e-3;
    let mu = 0.001;
    let area = std::f64::consts::PI * (d / 2.0f64).powi(2);
    let q_low = u_avg_low * area;
    let q_high = u_avg_high * area;
    let dp_per_q_analytical = (128.0 * mu * l_total) / (std::f64::consts::PI * d.powi(4));
    let dp_per_q_low = dp_low / q_low;
    let dp_per_q_high = dp_high / q_high;

    let u_ratio = high_solution.u_throat / low_solution.u_throat;
    let dp_ratio = dp_high / dp_low;

    println!("Micro-Stokes Scaling Test:");
    println!("  U_avg low:  {:.4e} m/s", u_avg_low);
    println!("  U_avg high: {:.4e} m/s", u_avg_high);
    println!("  DP low:     {:.6} Pa", dp_low);
    println!("  DP high:    {:.6} Pa", dp_high);
    println!("  DP ratio:   {:.3}", dp_ratio);
    println!("  U ratio:    {:.3}", u_ratio);
    println!("  DP/Q low:   {:.3e}", dp_per_q_low);
    println!("  DP/Q high:  {:.3e}", dp_per_q_high);
    println!("  DP/Q ref:   {:.3e}", dp_per_q_analytical);

    assert!(dp_low > 0.0 && dp_high > 0.0, "Pressure drop magnitude should be positive");
    assert!(high_solution.u_throat > low_solution.u_throat, "Throat velocity should increase with flow rate");

    // For Stokes flow, dp and velocity should scale linearly with flow rate.
    assert!(relative_error(dp_ratio, 2.0) < 0.30, "DP scaling deviates from linearity: {:.3}", dp_ratio);
    assert!(relative_error(u_ratio, 2.0) < 0.30, "Velocity scaling deviates from linearity: {:.3}", u_ratio);
    assert!(relative_error(dp_per_q_low, dp_per_q_analytical) < 2.80, "DP/Q magnitude deviates from Hagen–Poiseuille: {:.3e}", dp_per_q_low);
    assert!(relative_error(dp_per_q_high, dp_per_q_analytical) < 2.80, "DP/Q magnitude deviates from Hagen–Poiseuille: {:.3e}", dp_per_q_high);
}
