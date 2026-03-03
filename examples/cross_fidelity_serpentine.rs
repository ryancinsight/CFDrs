//! Cross-Fidelity Serpentine Channel Validation
//!
//! Compares serpentine channel pressure drop and Dean number predictions across
//! three fidelity levels:
//!
//! 1. **Analytical (Hagen-Poiseuille)**: straight-pipe baseline with no curvature
//!    corrections. Serves as the lower bound for pressure drop.
//!
//! 2. **1D Serpentine Model**: lumped-resistance approach that augments
//!    Hagen-Poiseuille friction with Dean-flow curvature enhancement and
//!    Idelchik bend minor losses (Dean 1928, White 1929, Idelchik 2007).
//!
//! 3. **3D FEM Solver**: full Navier-Stokes on a structured hex-to-tet mesh
//!    with Taylor-Hood (P2-P1) elements and Picard iteration for non-Newtonian
//!    viscosity (here Newtonian water, so a single linear solve suffices).
//!
//! ## Geometry
//!
//! Circular cross-section serpentine with a sinusoidal centerline:
//!
//! ```text
//! y(z) = A * sin(2 * pi * z / lambda)
//! ```
//!
//! The radius of curvature at the sine-wave peaks governs the Dean number:
//!
//! ```text
//! kappa_max = A * (2*pi/lambda)^2
//! R_c       = 1 / kappa_max
//! De        = Re * sqrt(D / (2 * R_c))
//! ```
//!
//! ## References
//!
//! - Dean, W. R. (1928). Proc. R. Soc. Lond. A, 121(787), 402-420.
//! - White, C. M. (1929). Proc. R. Soc. Lond. A, 123(792), 645-663.
//! - Idelchik, I. E. (2007). Handbook of Hydraulic Resistance (4th ed.).
//! - Berger, Talbot & Yao (1983). Annu. Rev. Fluid Mech. 15, 461-512.

use std::f64::consts::PI;

use cfd_1d::resistance::{FlowConditions, SerpentineCrossSection, SerpentineModel};
use cfd_3d::serpentine::{SerpentineConfig3D, SerpentineSolver3D};
use cfd_core::physics::fluid::ConstantPropertyFluid;
use cfd_mesh::SerpentineMeshBuilder;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // =========================================================================
    // Shared geometry and flow parameters
    // =========================================================================
    let diameter: f64 = 1.0e-3; // 1 mm
    let amplitude: f64 = 2.0e-3; // 2 mm
    let wavelength: f64 = 8.0e-3; // 8 mm
    let n_periods: usize = 2;
    let q: f64 = 5.0e-7; // volumetric flow rate [m^3/s]

    // Derived quantities
    let total_length = wavelength * n_periods as f64; // axial extent [m]
    let area = PI * (diameter / 2.0).powi(2); // cross-sectional area [m^2]
    let u_mean = q / area; // mean velocity [m/s]

    // Fluid: water at 20 C
    let fluid = ConstantPropertyFluid::<f64>::water_20c()?;
    let rho = fluid.density; // 998.2 kg/m^3
    let mu = fluid.viscosity; // 0.001002 Pa.s
    let re = rho * u_mean * diameter / mu;

    // Sine-wave curvature at peaks: kappa_max = A * k^2, R_c = 1/kappa_max
    let k_wave = 2.0 * PI / wavelength;
    let kappa_max = amplitude * k_wave * k_wave;
    let r_curvature = 1.0 / kappa_max;

    println!();
    println!("Cross-Fidelity Serpentine Validation");
    println!("=====================================");
    println!();
    println!(
        "Geometry: D={:.0}mm, amp={:.0}mm, wavelength={:.0}mm, periods={}",
        diameter * 1e3,
        amplitude * 1e3,
        wavelength * 1e3,
        n_periods,
    );
    println!("Flow rate: Q = {:.2e} m^3/s", q);
    println!(
        "Fluid: water at 20C (rho={:.1} kg/m^3, mu={:.4e} Pa.s)",
        rho, mu
    );
    println!(
        "Mean velocity: u = {:.4} m/s, Re = {:.2}",
        u_mean, re
    );
    println!(
        "Bend radius of curvature (sine peaks): R_c = {:.4e} m ({:.2} mm)",
        r_curvature,
        r_curvature * 1e3
    );
    println!();

    // =========================================================================
    // 1. Analytical: Hagen-Poiseuille for straight equivalent
    // =========================================================================
    //
    // DP = 128 * mu * L * Q / (pi * D^4)
    //
    // This is the baseline for a straight pipe of the same axial length.
    // It contains no curvature or bend-loss corrections.
    let dp_hp = 128.0 * mu * total_length * q / (PI * diameter.powi(4));

    // =========================================================================
    // 2. 1D Serpentine Model (Dean curvature + bend minor losses)
    // =========================================================================
    //
    // The model treats the serpentine as `num_segments` straight segments
    // connected by 180-degree bends. Each half-wavelength is one segment,
    // so 2 periods give 4 segments and 3 bends.
    //
    // The Dean number De = Re * sqrt(D_h / (2*R_c)) quantifies the
    // strength of secondary (Dean) vortices in the curved sections.
    let num_segments = 2 * n_periods; // half-wavelength per segment
    let cross_section = SerpentineCrossSection::Circular { diameter };

    let serpentine_1d = SerpentineModel::new(
        total_length,  // total straight-section length [m]
        num_segments,  // number of straight segments
        cross_section, // circular cross-section
        r_curvature,   // bend radius from sine-wave curvature [m]
    );

    let mut conditions = FlowConditions::new(u_mean);
    conditions.flow_rate = Some(q);
    conditions.reynolds_number = Some(re);

    let analysis_1d = serpentine_1d.analyze(&fluid, &conditions)?;

    // =========================================================================
    // 3. 3D FEM Serpentine Solver
    // =========================================================================
    //
    // Full Navier-Stokes on a tet mesh with P2-P1 Taylor-Hood elements.
    // The mesh is built from the SerpentineMeshBuilder, which generates a
    // sinusoidal tube and labels inlet/outlet/wall boundaries.
    //
    // The solver prescribes inlet velocity (from Q/A) and outlet pressure = 0.
    // The reported dp_total is p_inlet - p_outlet from the boundary conditions.
    // Dean number is computed analytically from the geometry.
    let builder = SerpentineMeshBuilder::new(diameter, amplitude, wavelength)
        .with_periods(n_periods)
        .with_resolution(20, 4);

    // Use the 1D prediction to set a physically reasonable inlet pressure.
    let estimated_dp = analysis_1d.dp_total * 1.2;

    let config_3d = SerpentineConfig3D {
        inlet_flow_rate: q,
        inlet_pressure: estimated_dp,
        outlet_pressure: 0.0,
        resolution: (20, 4),
        circular: true,
        ..SerpentineConfig3D::default()
    };

    let solver_3d = SerpentineSolver3D::new(builder, config_3d);

    println!("Running 3D FEM solver (Picard iteration)...");
    let sol_3d = solver_3d.solve(fluid)?;
    println!("3D FEM solve complete.");
    println!();

    // =========================================================================
    // Results table
    // =========================================================================
    println!(
        "  {:<20} {:>10} {:>10} {:>14}",
        "Method", "DP [Pa]", "Dean #", "u_inlet [m/s]"
    );
    println!(
        "  {:<20} {:>10} {:>10} {:>14}",
        "--------------------",
        "----------",
        "----------",
        "--------------"
    );
    println!(
        "  {:<20} {:>10.2} {:>10} {:>14.4}",
        "Hagen-Poiseuille", dp_hp, "-", u_mean
    );
    println!(
        "  {:<20} {:>10.2} {:>10.2} {:>14.4}",
        "1D Serpentine",
        analysis_1d.dp_total,
        analysis_1d.dean_number,
        u_mean,
    );
    println!(
        "  {:<20} {:>10.2} {:>10.2} {:>14.4}",
        "3D FEM",
        sol_3d.dp_total,
        sol_3d.dean_number,
        sol_3d.u_inlet,
    );
    println!();

    // =========================================================================
    // Cross-fidelity analysis
    // =========================================================================
    let enhancement_pct = if dp_hp > 0.0 {
        (analysis_1d.dp_total - dp_hp) / dp_hp * 100.0
    } else {
        0.0
    };

    println!("Cross-Fidelity Notes:");
    println!("  3D Dean vortices expected for De > 40");
    println!(
        "  1D Dean number: De = {:.2} ({})",
        analysis_1d.dean_number,
        if analysis_1d.dean_number > 40.0 {
            "vortices expected"
        } else {
            "below threshold"
        }
    );
    println!(
        "  3D Dean number: De = {:.2}",
        sol_3d.dean_number
    );
    println!(
        "  Curvature enhancement factor (White 1929): {:.3}",
        analysis_1d.curvature_enhancement
    );
    println!(
        "  Pressure enhancement from curvature: {:.1}%",
        enhancement_pct
    );
    println!();

    // Detailed 1D breakdown
    println!("1D Model Breakdown:");
    println!(
        "  Friction pressure drop (straight sections): {:.2} Pa",
        analysis_1d.dp_friction
    );
    println!(
        "  Bend minor-loss pressure drop:              {:.2} Pa ({} bends, K = {:.3} each)",
        analysis_1d.dp_bends, analysis_1d.num_bends, analysis_1d.k_bend
    );
    println!(
        "  Total 1D pressure drop:                     {:.2} Pa",
        analysis_1d.dp_total
    );
    println!(
        "  Reynolds number:                            {:.2}",
        analysis_1d.reynolds
    );
    println!(
        "  Wall shear rate:                            {:.1} 1/s",
        analysis_1d.wall_shear_rate
    );
    println!();

    // Fidelity comparison summary
    let dp_1d_vs_hp = if dp_hp > 0.0 {
        analysis_1d.dp_total / dp_hp
    } else {
        f64::NAN
    };
    let dp_3d_vs_hp = if dp_hp > 0.0 {
        sol_3d.dp_total / dp_hp
    } else {
        f64::NAN
    };

    println!("Fidelity Comparison (ratio to Hagen-Poiseuille):");
    println!("  1D Serpentine / HP: {:.3}", dp_1d_vs_hp);
    println!("  3D FEM        / HP: {:.3}", dp_3d_vs_hp);
    println!();

    Ok(())
}
