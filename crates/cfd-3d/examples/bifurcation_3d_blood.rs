//! 3D Bifurcation with Blood Flow
//!
//! Solves 3D Navier-Stokes in a Y-junction bifurcation using the Casson
//! blood model. Validates Murray's law diameter scaling and mass conservation.
//!
//! # Physics
//!
//! Murray's law: D_parent^3 = D_daughter1^3 + D_daughter2^3
//! For symmetric bifurcation: D_daughter = D_parent / 2^(1/3)
//!
//! # Reference
//!
//! Murray, C.D. (1926). "The Physiological Principle of Minimum Work".
//! *Proc. Natl. Acad. Sci.* 12(3):207-214.

use cfd_3d::bifurcation::{BifurcationConfig3D, BifurcationGeometry3D, BifurcationSolver3D};
use cfd_core::physics::fluid::blood::CassonBlood;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("3D Bifurcation — Blood Flow (Casson Model)");
    println!("============================================\n");

    // Symmetric bifurcation geometry (arteriole scale)
    let d_parent = 100e-6;   // 100 um
    let d_daughter = 80e-6;  // 80 um (slightly larger than Murray optimum)
    let l_parent = 1e-3;     // 1 mm
    let l_daughter = 1e-3;
    let l_transition = 100e-6;

    let geometry = BifurcationGeometry3D::<f64>::symmetric(
        d_parent, d_daughter, l_parent, l_daughter, l_transition,
    );

    // Murray's law analysis
    let d_murray = d_parent / 2.0_f64.powf(1.0 / 3.0);
    let murray_error = ((d_daughter - d_murray) / d_murray).abs() * 100.0;

    println!("Geometry:");
    println!("  Parent diameter:    {:.1} um", d_parent * 1e6);
    println!("  Daughter diameter:  {:.1} um", d_daughter * 1e6);
    println!("  Murray optimal:     {:.1} um", d_murray * 1e6);
    println!("  Murray deviation:   {:.1}%", murray_error);
    println!("  Total volume:       {:.3e} m^3", geometry.total_volume());

    // Solver configuration
    let config = BifurcationConfig3D {
        inlet_flow_rate: 1e-8,   // 10 nL/s
        inlet_pressure: 100.0,   // 100 Pa gauge
        ..BifurcationConfig3D::default()
    };

    let solver = BifurcationSolver3D::new(geometry, config);

    // Solve with Casson blood
    let blood = CassonBlood::<f64>::normal_blood();
    println!("\nFluid: Casson Blood");
    println!("  Yield stress:          {:.4} Pa", blood.yield_stress);
    println!("  Inf-shear viscosity:   {:.4e} Pa.s", blood.infinite_shear_viscosity);
    println!("  Density:               {:.1} kg/m^3", blood.density);

    println!("\nSolving...");
    let solution = solver.solve(blood)?;

    println!("\nFlow Results:");
    println!("  Parent flow:      {:.3e} m^3/s", solution.q_parent);
    println!("  Daughter 1 flow:  {:.3e} m^3/s", solution.q_daughter1);
    println!("  Daughter 2 flow:  {:.3e} m^3/s", solution.q_daughter2);

    println!("\nVelocity:");
    println!("  Parent mean:      {:.4} m/s", solution.u_parent_mean);
    println!("  Daughter 1 mean:  {:.4} m/s", solution.u_daughter1_mean);
    println!("  Daughter 2 mean:  {:.4} m/s", solution.u_daughter2_mean);

    println!("\nPressure:");
    println!("  Inlet:            {:.2} Pa", solution.p_inlet);
    println!("  Junction mid:     {:.2} Pa", solution.p_junction_mid);
    println!("  Daughter 1 out:   {:.2} Pa", solution.p_daughter1_outlet);
    println!("  Daughter 2 out:   {:.2} Pa", solution.p_daughter2_outlet);

    println!("\nWall Shear Stress:");
    println!("  Parent:           {:.4} Pa", solution.wall_shear_stress_parent);
    println!("  Daughter 1:       {:.4} Pa", solution.wall_shear_stress_daughter1);
    println!("  Daughter 2:       {:.4} Pa", solution.wall_shear_stress_daughter2);

    println!("\nValidation:");
    let mass_ok = solution.is_mass_conserved(1e-6);
    println!("  Mass conservation error: {:.2e}", solution.mass_conservation_error);
    println!("  Mass conserved (<1e-6):  {}", if mass_ok { "PASS" } else { "FAIL" });

    // Flow split ratio
    let q_total_out = solution.q_daughter1 + solution.q_daughter2;
    let split_1 = solution.q_daughter1 / q_total_out * 100.0;
    let split_2 = solution.q_daughter2 / q_total_out * 100.0;
    println!("  Flow split: {:.1}% / {:.1}% (expect ~50/50 for symmetric)", split_1, split_2);

    println!("\nBifurcation 3D example completed.");
    Ok(())
}
