//! Cross-Fidelity Bifurcation Validation
//!
//! Compares a symmetric Y-junction bifurcation across three fidelity levels:
//!
//! 1. **Analytical** -- Hagen-Poiseuille pressure drop with Murray's law diameter check.
//! 2. **1D lumped-parameter** -- `TwoWayBranchJunction` from `cfd-1d`, which evaluates
//!    the same Poiseuille resistance model per channel segment.
//! 3. **3D FEM** -- `BifurcationSolver3D` from `cfd-3d`, solving the full incompressible
//!    Navier-Stokes equations on a meshed bifurcation domain.
//!
//! 2D is intentionally skipped: the 2D bifurcation solver requires domain decomposition
//! and boundary-fitted mesh configuration that would obscure the cross-fidelity comparison.
//!
//! # Geometry
//!
//! A symmetric bifurcation with parent diameter D_p = 100 um and daughter diameter
//! D_d = 80 um. The Murray-optimal daughter diameter for a symmetric two-way split is
//!
//! ```text
//! D_opt = D_p / 2^(1/3) = 79.37 um
//! ```
//!
//! so the chosen D_d = 80 um deviates by less than 1 % from the optimum.
//!
//! # Murray's Cube Law
//!
//! Murray (1926) showed that the total power dissipation (viscous losses plus metabolic
//! cost of blood volume) in a branching network is minimised when
//!
//! ```text
//! D_0^3 = D_1^3 + D_2^3
//! ```
//!
//! For a symmetric split with D_1 = D_2 = D_d this reduces to D_0^3 = 2 D_d^3, i.e.
//! D_d = D_0 / 2^(1/3).
//!
//! # References
//!
//! - Murray, C.D. (1926). "The Physiological Principle of Minimum Work: I. The Vascular
//!   System and the Cost of Blood Volume." *Proceedings of the National Academy of
//!   Sciences* 12(3):207-214.
//! - Hagen, G.H.L. (1839). *Annalen der Physik und Chemie* 46:423-442.
//! - Poiseuille, J.L.M. (1846). *Memoires des Savants Etrangers* 9:433-544.
//!
//! # Running
//!
//! ```bash
//! cargo run --example cross_fidelity_bifurcation --release --no-default-features
//! ```

use std::f64::consts::PI;

use cfd_1d::channel::{Channel, ChannelGeometry};
use cfd_1d::junctions::branching::TwoWayBranchJunction;
use cfd_3d::bifurcation::{BifurcationConfig3D, BifurcationGeometry3D, BifurcationSolver3D};
use cfd_core::physics::fluid::ConstantPropertyFluid;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Cross-Fidelity Bifurcation Validation");
    println!("=======================================\n");

    // =====================================================================
    // Shared geometry
    // =====================================================================
    let d_parent: f64 = 100e-6; // 100 um
    let d_daughter: f64 = 80e-6; // 80 um (near Murray optimal 79.37 um)
    let l_parent: f64 = 1e-3; // 1 mm
    let l_daughter: f64 = 1e-3; // 1 mm
    let l_transition: f64 = 100e-6; // 100 um conical transition
    let q_parent: f64 = 1e-8; // 10 nL/s
    let p_inlet: f64 = 100.0; // Pa

    // Fluid: water at 20 C (Newtonian, constant properties)
    let mu: f64 = 8.9e-4; // dynamic viscosity [Pa.s]
    let water = ConstantPropertyFluid::<f64>::new(
        "Water (~25 C)".to_string(),
        998.0,  // density [kg/m^3]
        mu,     // dynamic viscosity [Pa.s]
        4186.0, // specific heat [J/(kg.K)]
        0.599,  // thermal conductivity [W/(m.K)]
        1482.0, // speed of sound [m/s]
    );

    println!(
        "Geometry: D_parent = {:.0} um, D_daughter = {:.0} um",
        d_parent * 1e6,
        d_daughter * 1e6
    );
    println!(
        "          L_parent = {:.0} um, L_daughter = {:.0} um, L_transition = {:.0} um",
        l_parent * 1e6,
        l_daughter * 1e6,
        l_transition * 1e6
    );
    println!(
        "Flow:     Q = {:.2e} m^3/s  ({:.0} nL/s)",
        q_parent,
        q_parent * 1e9
    );
    println!(
        "Fluid:    mu = {:.2e} Pa.s  (water)\n",
        mu
    );

    // =====================================================================
    // Murray's law check
    // =====================================================================
    let d_murray = d_parent / 2.0_f64.powf(1.0 / 3.0);
    let murray_error = ((d_daughter - d_murray) / d_murray).abs() * 100.0;

    println!(
        "Murray optimal daughter diameter: {:.1} um",
        d_murray * 1e6
    );
    println!(
        "Actual daughter diameter:         {:.1} um  (deviation: {:.1}%)\n",
        d_daughter * 1e6,
        murray_error
    );

    // =====================================================================
    // 1. Analytical -- Hagen-Poiseuille
    // =====================================================================
    //
    //   dP = 128 mu L Q / (pi D^4)
    //
    // Parent carries Q; each daughter carries Q/2.
    let dp_parent_hp = 128.0 * mu * l_parent * q_parent / (PI * d_parent.powi(4));
    let dp_daughter_hp = 128.0 * mu * l_daughter * (q_parent / 2.0) / (PI * d_daughter.powi(4));
    let dp_total_hp = dp_parent_hp + dp_daughter_hp;

    // =====================================================================
    // 2. 1D lumped-parameter model
    // =====================================================================
    //
    // Build Channel objects and a TwoWayBranchJunction.  For a Newtonian fluid
    // the junction's pressure_drop() is algebraically identical to Hagen-Poiseuille
    // (128 mu Q L / pi D^4), so the 1D result should match the analytical to
    // machine precision.
    //
    // Note: TwoWayBranchJunction::solve() requires F: Fluid<T> + Copy.
    // ConstantPropertyFluid contains a String field and therefore is not Copy.
    // We use the public associated function pressure_drop() directly, which
    // accepts &F and reproduces the same physics.

    let parent_ch = Channel::new(ChannelGeometry::<f64>::circular(l_parent, d_parent, 0.0));
    let daughter1_ch = Channel::new(ChannelGeometry::<f64>::circular(l_daughter, d_daughter, 0.0));
    let daughter2_ch = Channel::new(ChannelGeometry::<f64>::circular(l_daughter, d_daughter, 0.0));

    // Create junction (used for Murray's law deviation report)
    let junction = TwoWayBranchJunction::new(
        parent_ch.clone(),
        daughter1_ch.clone(),
        daughter2_ch,
        0.5, // symmetric 50/50 split
    );

    let murray_dev_1d = junction.murray_law_deviation();

    // Temperature and thermodynamic pressure for viscosity evaluation
    let temperature = 293.15_f64; // 20 C
    let p_thermo = 101_325.0_f64; // 1 atm

    let dp_parent_1d =
        TwoWayBranchJunction::<f64>::pressure_drop(&water, q_parent, &parent_ch, temperature, p_thermo);

    let q_daughter = q_parent / 2.0;
    let dp_daughter_1d =
        TwoWayBranchJunction::<f64>::pressure_drop(&water, q_daughter, &daughter1_ch, temperature, p_thermo);

    let dp_total_1d = dp_parent_1d + dp_daughter_1d;

    // Mass conservation: flow split is prescribed, so error is zero by construction
    let q1_1d = q_daughter;
    let q2_1d = q_parent - q_daughter;
    let mass_err_1d = ((q1_1d + q2_1d) - q_parent).abs();
    let split_1_1d = q1_1d / q_parent * 100.0;
    let split_2_1d = q2_1d / q_parent * 100.0;

    // =====================================================================
    // 3. 3D FEM -- BifurcationSolver3D
    // =====================================================================
    let geometry_3d = BifurcationGeometry3D::<f64>::symmetric(
        d_parent,
        d_daughter,
        l_parent,
        l_daughter,
        l_transition,
    );

    let config_3d = BifurcationConfig3D {
        inlet_flow_rate: q_parent,
        inlet_pressure: p_inlet,
        mesh_resolution: 4,
        ..BifurcationConfig3D::default()
    };

    let solver_3d = BifurcationSolver3D::new(geometry_3d, config_3d);

    println!("Solving 3D FEM (mesh_resolution = 4) ...");

    let sol_3d = solver_3d.solve(water)?;

    let dp_total_3d = sol_3d.dp_parent + sol_3d.dp_daughter1;
    let q_sum_3d = sol_3d.q_daughter1 + sol_3d.q_daughter2;
    let split_1_3d = if q_sum_3d > 0.0 {
        sol_3d.q_daughter1 / q_sum_3d * 100.0
    } else {
        50.0
    };
    let split_2_3d = 100.0 - split_1_3d;
    let mass_err_3d = sol_3d.mass_conservation_error;

    // =====================================================================
    // Results table
    // =====================================================================
    println!();
    println!("  {:<24} {:>14}  {:>12}  {:>10}", "Method", "DP total [Pa]", "Flow split", "Mass err");
    println!("  {:<24} {:>14}  {:>12}  {:>10}", "------------------------", "--------------", "------------", "----------");

    println!(
        "  {:<24} {:>14.4}  {:>12}  {:>10}",
        "Hagen-Poiseuille",
        dp_total_hp,
        "50/50",
        "-"
    );

    println!(
        "  {:<24} {:>14.4}  {:>5.1}/{:.1}%  {:>10.2e}",
        "1D Junction",
        dp_total_1d,
        split_1_1d,
        split_2_1d,
        mass_err_1d
    );

    println!(
        "  {:<24} {:>14.4}  {:>5.1}/{:.1}%  {:>10.2e}",
        "3D FEM (mesh_res=4)",
        dp_total_3d,
        split_1_3d,
        split_2_3d,
        mass_err_3d
    );

    // =====================================================================
    // Component breakdown
    // =====================================================================
    println!("\nPressure-drop breakdown [Pa]:");
    println!(
        "  {:<24} {:>14}  {:>14}",
        "", "DP parent", "DP daughter"
    );
    println!(
        "  {:<24} {:>14.4}  {:>14.4}",
        "Hagen-Poiseuille", dp_parent_hp, dp_daughter_hp
    );
    println!(
        "  {:<24} {:>14.4}  {:>14.4}",
        "1D Junction", dp_parent_1d, dp_daughter_1d
    );
    println!(
        "  {:<24} {:>14.4}  {:>14.4}",
        "3D FEM", sol_3d.dp_parent, sol_3d.dp_daughter1
    );

    // =====================================================================
    // 3D wall shear stress
    // =====================================================================
    println!("\n3D wall shear stress [Pa]:");
    println!(
        "  Parent:    {:.4}",
        sol_3d.wall_shear_stress_parent
    );
    println!(
        "  Daughter1: {:.4}",
        sol_3d.wall_shear_stress_daughter1
    );
    println!(
        "  Daughter2: {:.4}",
        sol_3d.wall_shear_stress_daughter2
    );

    // =====================================================================
    // Murray's law summary
    // =====================================================================
    println!(
        "\nMurray's law deviation (1D junction): {:.4} ({:.2}%)",
        murray_dev_1d,
        murray_dev_1d * 100.0
    );

    // =====================================================================
    // Cross-fidelity consistency check
    // =====================================================================
    //
    // The 1D model uses the same Poiseuille resistance formula as the analytical
    // solution, so they should agree to machine precision.  The 3D result includes
    // entrance/exit effects, the conical transition zone, and numerical discretisation
    // error, so a looser tolerance is appropriate.
    let hp_vs_1d = ((dp_total_1d - dp_total_hp) / dp_total_hp).abs();
    let hp_vs_3d = ((dp_total_3d - dp_total_hp) / dp_total_hp).abs();

    let tolerance_1d = 1e-10; // machine precision
    let tolerance_3d = 0.50; // 50 %: coarse mesh + transition losses

    let pass_1d = hp_vs_1d < tolerance_1d;
    let pass_3d = hp_vs_3d < tolerance_3d;
    let pass_all = pass_1d && pass_3d;

    println!("\nCross-Fidelity Consistency:");
    println!(
        "  HP vs 1D:  {:.2e}  (tol {:.0e})  {}",
        hp_vs_1d,
        tolerance_1d,
        if pass_1d { "PASS" } else { "FAIL" }
    );
    println!(
        "  HP vs 3D:  {:.2e}  (tol {:.0e})  {}",
        hp_vs_3d,
        tolerance_3d,
        if pass_3d { "PASS" } else { "FAIL" }
    );
    println!(
        "\nCross-Fidelity Consistency: {}",
        if pass_all { "PASS" } else { "FAIL" }
    );

    Ok(())
}
