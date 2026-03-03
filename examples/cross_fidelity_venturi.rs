//! Cross-Fidelity Venturi Validation
//!
//! Compares the same Venturi geometry solved at four fidelity levels:
//!
//! 1. **Analytical** -- Bernoulli's equation (inviscid, incompressible)
//! 2. **1D Resistance** -- `VenturiModel` with discharge coefficient, friction,
//!    and expansion losses
//! 3. **2D SIMPLE** -- `VenturiSolver2D` (Navier-Stokes FVM on staggered grid)
//! 4. **3D FEM** -- `VenturiSolver3D` (Navier-Stokes FEM on structured hex mesh)
//!
//! ## Shared Geometry
//!
//! | Parameter       | Value   |
//! |-----------------|---------|
//! | Inlet diameter  | 2.0 mm  |
//! | Throat diameter | 1.0 mm  |
//! | l_inlet         | 3.0 mm  |
//! | l_convergent    | 2.0 mm  |
//! | l_throat        | 1.0 mm  |
//! | l_divergent     | 4.0 mm  |
//! | l_outlet        | 3.0 mm  |
//!
//! ## Fluid
//!
//! Water at 20 degrees C: rho = 998 kg/m^3, mu = 8.9e-4 Pa.s
//!
//! ## References
//!
//! - Venturi, G. B. (1797). *Recherches experimentales sur le principe de la
//!   communication laterale du mouvement dans les fluides.*
//! - ISO 5167-4:2003. "Measurement of fluid flow by means of pressure
//!   differential devices -- Part 4: Venturi tubes."
//! - White, F. M. (2011). *Fluid Mechanics* (7th ed.). McGraw-Hill.
//! - Reader-Harris, M. (2015). *Orifice Plates and Venturi Tubes.* Springer.

use std::f64::consts::PI;

// 1D resistance model
use cfd_1d::resistance::{FlowConditions, VenturiModel};

// 2D FVM solver -- VenturiGeometry from the venturi_flow module is a struct
// describing 2D channel dimensions; it differs from the cfd_1d VenturiGeometry
// enum, so we alias it here.
use cfd_2d::solvers::venturi_flow::{
    VenturiGeometry as VenturiGeometry2D, VenturiSolver2D,
};

// BloodModel with Newtonian(mu) variant, re-exported from cfd_core
use cfd_core::physics::fluid::BloodModel;

// 3D FEM solver
use cfd_3d::venturi::{VenturiConfig3D, VenturiSolver3D};
use cfd_mesh::VenturiMeshBuilder;

// Shared fluid model
use cfd_core::physics::fluid::ConstantPropertyFluid;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // ===================================================================
    // Shared geometry and fluid parameters
    // ===================================================================
    let d_inlet: f64 = 2.0e-3; // 2 mm
    let d_throat: f64 = 1.0e-3; // 1 mm
    let l_inlet: f64 = 3.0e-3;
    let l_convergent: f64 = 2.0e-3;
    let l_throat: f64 = 1.0e-3;
    let l_divergent: f64 = 4.0e-3;
    let l_outlet: f64 = 3.0e-3;
    let total_length = l_inlet + l_convergent + l_throat + l_divergent + l_outlet;

    // Water at 20 degrees C
    let rho: f64 = 998.0; // kg/m^3
    let mu: f64 = 8.9e-4; // Pa.s

    // Volumetric flow rate (same for every fidelity level)
    let q: f64 = 1.0e-7; // m^3/s

    // Circular cross-sectional areas
    let area_inlet = PI * (d_inlet / 2.0).powi(2);
    let area_throat = PI * (d_throat / 2.0).powi(2);

    // Velocities from continuity (A1 * u1 = A2 * u2 = Q)
    let u_inlet_bernoulli = q / area_inlet;
    let u_throat_bernoulli = q / area_throat;

    println!("Cross-Fidelity Venturi Validation");
    println!("==================================\n");
    println!(
        "Geometry: D_inlet={:.0}mm, D_throat={:.0}mm, Q={:.1e} m^3/s\n",
        d_inlet * 1e3,
        d_throat * 1e3,
        q,
    );

    // ===================================================================
    // 0. Analytical -- Bernoulli (inviscid, incompressible)
    //
    //   P1 + 0.5 rho u1^2 = P2 + 0.5 rho u2^2
    //   => dP = 0.5 rho (u_throat^2 - u_inlet^2)
    // ===================================================================
    let dp_bernoulli = 0.5 * rho * (u_throat_bernoulli.powi(2) - u_inlet_bernoulli.powi(2));

    // ===================================================================
    // 1. 1D Resistance -- VenturiModel::symmetric
    //
    //   Includes discharge coefficient, throat friction, and
    //   Borda-Carnot expansion loss with partial pressure recovery.
    // ===================================================================
    let model_1d = VenturiModel::<f64>::symmetric(
        d_inlet,
        d_throat,
        l_throat,
        total_length,
    );

    let conditions = FlowConditions {
        velocity: None,
        flow_rate: Some(q),
        reynolds_number: None,
        shear_rate: None,
        temperature: 293.15,
        pressure: 101325.0,
    };

    let water = ConstantPropertyFluid::<f64>::water_20c()?;
    let analysis_1d = model_1d.analyze(&water, &conditions)?;

    // ===================================================================
    // 2. 2D SIMPLE -- VenturiSolver2D (Navier-Stokes FVM)
    //
    //   The 2D solver uses rectangular channel widths rather than
    //   diameters.  For an approximate cross-fidelity comparison we set
    //   w_inlet = d_inlet, w_throat = d_throat, and height = d_inlet,
    //   giving a square inlet cross-section (d x d) instead of
    //   a circular one (pi d^2 / 4).  The velocity is recomputed
    //   from Q to match this rectangular area.
    //
    //   Note: the 2D geometry does not model a separate outlet section;
    //   the diverging section implicitly recovers to the outlet.
    // ===================================================================
    let w_inlet_2d = d_inlet;
    let w_throat_2d = d_throat;
    let height_2d = d_inlet;

    let geom_2d = VenturiGeometry2D::new(
        w_inlet_2d,
        w_throat_2d,
        l_inlet,
        l_convergent,
        l_throat,
        l_divergent,
        height_2d,
    );

    // Inlet velocity consistent with Q through the rectangular cross-section
    let area_inlet_2d = w_inlet_2d * height_2d;
    let u_inlet_2d = q / area_inlet_2d;

    let blood_model = BloodModel::Newtonian(mu);
    let mut solver_2d = VenturiSolver2D::new(geom_2d, blood_model, rho, 40, 20);
    let sol_2d = solver_2d.solve(u_inlet_2d)?;

    // dp_throat in VenturiFlowSolution is (p_throat - p_inlet), which is
    // negative in a Venturi; negate it to get a positive pressure drop
    // comparable to the Bernoulli and 1D results.
    let dp_2d = sol_2d.dp_throat.abs();
    let u_throat_2d = sol_2d.u_throat;

    // ===================================================================
    // 3. 3D FEM -- VenturiSolver3D (Navier-Stokes on hex mesh)
    //
    //   Uses VenturiMeshBuilder for the geometry and solves the full
    //   3D incompressible Navier-Stokes equations via FEM.
    // ===================================================================
    let builder = VenturiMeshBuilder::new(
        d_inlet,
        d_throat,
        l_inlet,
        l_convergent,
        l_throat,
        l_divergent,
        l_outlet,
    );

    let config_3d = VenturiConfig3D {
        inlet_flow_rate: q,
        inlet_pressure: 200.0,
        outlet_pressure: 0.0,
        resolution: (30, 5),
        ..VenturiConfig3D::default()
    };

    let solver_3d = VenturiSolver3D::new(builder, config_3d);
    let sol_3d = solver_3d.solve(water.clone())?;

    // dp_throat in VenturiSolution3D is (p_inlet - p_throat), already positive
    let dp_3d = sol_3d.dp_throat;
    let u_throat_3d = sol_3d.u_throat;

    // ===================================================================
    // Comparison table
    // ===================================================================
    println!(
        "  {:<16}{:>16}{:>16}  {}",
        "Method", "DP_throat [Pa]", "u_throat [m/s]", "Notes"
    );
    println!(
        "  {:<16}{:>16}{:>16}  {}",
        "----------", "--------------", "--------------", "-----"
    );
    println!(
        "  {:<16}{:>16.2}{:>16.4}  {}",
        "Bernoulli", dp_bernoulli, u_throat_bernoulli, "Inviscid"
    );
    println!(
        "  {:<16}{:>16.2}{:>16.4}  {}",
        "1D Resistance", analysis_1d.dp_total, analysis_1d.throat_velocity, "VenturiModel"
    );
    println!(
        "  {:<16}{:>16.2}{:>16.4}  {}",
        "2D SIMPLE", dp_2d, u_throat_2d, "40x20 grid"
    );
    println!(
        "  {:<16}{:>16.2}{:>16.4}  {}",
        "3D FEM", dp_3d, u_throat_3d, "30x5 grid"
    );

    // ===================================================================
    // Cross-fidelity consistency
    // ===================================================================
    let pct_diff = |a: f64, b: f64| -> f64 { ((a - b) / b).abs() * 100.0 };

    let err_1d = pct_diff(analysis_1d.dp_total, dp_bernoulli);
    let err_2d = pct_diff(dp_2d, dp_bernoulli);
    let err_3d = pct_diff(dp_3d, dp_bernoulli);

    println!("\nCross-Fidelity Consistency:");
    println!("  1D vs Bernoulli: {:.1}%", err_1d);
    println!("  2D vs Bernoulli: {:.1}%", err_2d);
    println!("  3D vs Bernoulli: {:.1}%", err_3d);

    let threshold = 50.0;
    if err_1d < threshold && err_2d < threshold && err_3d < threshold {
        println!(
            "\nAll solvers within {:.0}% of analytical -- PASS",
            threshold
        );
    } else {
        println!(
            "\nWARNING: One or more solvers exceed {:.0}% deviation -- REVIEW",
            threshold
        );
    }

    // ===================================================================
    // Additional diagnostics
    // ===================================================================
    println!("\n--- Additional diagnostics ---");
    println!("  Re_inlet  = {:.1}", rho * u_inlet_bernoulli * d_inlet / mu);
    println!("  Re_throat = {:.1}", analysis_1d.throat_reynolds);
    println!("  1D C_d    = {:.4}", analysis_1d.discharge_coefficient);
    println!("  1D f      = {:.4}", analysis_1d.friction_factor);
    println!("  3D mass err = {:.2e}", sol_3d.mass_error);

    Ok(())
}
