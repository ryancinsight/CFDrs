use cfd_3d::venturi::{VenturiSolver3D, VenturiConfig3D};
use cfd_mesh::geometry::venturi::VenturiMeshBuilder;
use cfd_core::physics::fluid::blood::CassonBlood;
use approx::assert_relative_eq;

#[test]
fn validate_venturi_blood_flow() {
    // 1. Geometry: Standard Venturi
    // D_in = 1mm, D_throat = 0.5mm
    let d_in = 1.0e-3;
    let d_throat = 0.5e-3;
    let l_inlet = 2.0e-3;
    let l_convergent = 1.0e-3;
    let l_throat = 1.0e-3;
    let l_divergent = 1.0e-3;
    let l_outlet = 2.0e-3;

    let builder = VenturiMeshBuilder::new(
        d_in, d_throat, 
        l_inlet, l_convergent, l_throat, l_divergent, l_outlet
    );

    // 2. Fluid Properties: Human Blood (Casson Model)
    // Limits: mu_inf = 3.45 mPa.s, tau_y = 5.6 mPa
    let fluid = CassonBlood::<f64>::normal_blood();
    
    println!("Fluid Model: Casson Blood");
    println!("  Density: {:.1} kg/m^3", fluid.density);
    println!("  Yield Stress: {:.4} Pa", fluid.yield_stress);
    println!("  Inf Viscosity: {:.4} Pa.s", fluid.infinite_shear_viscosity);

    // 3. Flow Conditions
    // Low Flow Rate to allow shear-thinning effects to be significant vs inertial
    // Q = 1e-7 m^3/s (~6 ml/min)
    let q_in = 1.0e-7; 
    let area_in = std::f64::consts::PI * (d_in / 2.0f64).powi(2);
    let u_in = q_in / area_in;

    let config = VenturiConfig3D {
        inlet_flow_rate: q_in,
        inlet_pressure: 100.0, // Arbitrary reference, we check DP
        outlet_pressure: 0.0,
        resolution: (40, 6), // Coarse for speed, but fine enough for bulk flow
        circular: true,
        max_nonlinear_iterations: 25, // Iterating for viscosity
        nonlinear_tolerance: 1e-4,
    };

    println!("Venturi Test Setup:");
    println!("  D_in: {:.2e} m, D_throat: {:.2e} m", d_in, d_throat);
    println!("  Flow Rate: {:.2e} m^3/s", q_in);
    println!("  Inlet Velocity: {:.4} m/s", u_in);

    // 4. Solve
    let solver = VenturiSolver3D::new(builder, config);
    let solution = solver.solve(fluid).unwrap();

    // 5. Validation
    // 5.1 Mass Conservation
    println!("Validation Metrics:");
    println!("  Mass Error: {:.2e}", solution.mass_error);
    assert!(solution.mass_error.abs() < 1e-6, "Mass conservation failed");

    // 5.2 Physics Checks
    // Expect Pressure Drop > 0
    let dp = solution.p_inlet - solution.p_outlet;
    println!("  Total Pressure Drop: {:.2} Pa", dp);
    assert!(dp > 0.0, "Pressure drop must be positive");

    // Expect Pressure Drop > Bernoulli (Inviscid) due to Viscosity
    // Bernoulli: DP = 0.5 * rho * (u_th^2 - u_in^2)
    let area_th = std::f64::consts::PI * (d_throat / 2.0f64).powi(2);
    let u_th_avg = q_in / area_th;
    let dp_bernoulli = 0.5 * fluid.density * (u_th_avg.powi(2) - u_in.powi(2));
    
    println!("  Bernoulli Drop:      {:.2} Pa", dp_bernoulli);
    assert!(dp > dp_bernoulli, "Viscous pressure drop must exceed Bernoulli ideal. DP={:.2}, Bernoulli={:.2}", dp, dp_bernoulli);

    // 5.3 Throat Acceleration
    // Velocity ratio should roughly match Area ratio
    let area_ratio = area_in / area_th; // (1.0/0.5)^2 = 4
    let vel_ratio = solution.u_throat / solution.u_inlet;
    println!("  Area Ratio: {:.2}", area_ratio);
    println!("  Vel Ratio:  {:.2} (Max/Avg)", vel_ratio);
    
    // U_throat is max velocity (centerline), U_inlet is avg velocity (bulk).
    // For parabolic flow, U_max = 2 * U_avg.
    // So U_th_max / U_in_avg ≈ 2 * (A_in / A_th) = 2 * 4 = 8.
    // We expect roughly factor of 8 acceleration for developed flow.
    // Allow some margin for developing flow / numerical diffusion.
    assert!(vel_ratio > area_ratio, "Flow must accelerate in throat");
}
