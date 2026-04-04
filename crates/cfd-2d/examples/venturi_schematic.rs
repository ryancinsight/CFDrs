//! Venturi Schematic + 2D CFD Example
//!
//! Demonstrates the two-phase schematic-driven pattern for cfd-2d:
//!
//! **Phase 1 — Design (cfd-schematics)**
//! - Define Venturi topology as `NetworkBlueprint` (5 nodes, 4 channels)
//! - Render schematic PNG → `outputs/venturi_schematic.png`
//!
//! **Phase 2 — Simulation (cfd-2d)**
//! - Run `VenturiSolver2D` (ISO 5167, Newtonian blood, 60×30 grid)
//! - Map pressure results to `AnalysisOverlay` (BlueRed colormap)
//! - Render overlay PNG → `outputs/venturi_pressure_overlay.png`
//!
//! Run with:
//! `cargo run -p cfd-2d --example venturi_schematic`

use cfd_2d::solvers::ns_fvm::BloodModel;
use cfd_2d::solvers::venturi_flow::{BernoulliVenturi, VenturiGeometry, VenturiSolver2D};
use cfd_schematics::domain::model::{ChannelSpec, NetworkBlueprint, NodeKind, NodeSpec};
use cfd_schematics::visualizations::{
    create_plotters_renderer, AnalysisField, AnalysisOverlay, ColormapKind, RenderConfig,
    SchematicRenderer,
};
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🔬 Venturi Schematic + 2D CFD Example");
    println!("======================================");

    // ── 0. Output directory ───────────────────────────────────────────────────
    let out = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("outputs")
        .join("venturi");
    fs::create_dir_all(&out)?;

    // ── Phase 1: Design (cfd-schematics) ─────────────────────────────────────
    println!("\n📐 Phase 1: Schematic Design");

    // ISO 5167 Venturi dimensions (mm for schematic layout)
    // inlet(10mm) | converge(1mm) | throat(2mm) | diverge(3mm)
    let x0 = 0.0_f64;
    let x1 = 10.0; // end of inlet section
    let x2 = 11.0; // end of converging section
    let x3 = 13.0; // end of throat section
    let x4 = 16.0; // end of diverging section (total = 16 mm)
    let y = 5.0; // vertical centre of schematic

    let w_inlet = 10.0; // mm — inlet channel width for schematic
    let w_throat = 7.07; // mm — throat width (√2 ratio)
    let w_mid = (w_inlet + w_throat) / 2.0;

    let mut system = NetworkBlueprint::new_with_explicit_positions("venturi-iso5167");
    system.box_dims = (x4 + 2.0, y * 2.0);

    for (id, px, py, kind) in [
        ("n0", x0, y, NodeKind::Inlet),
        ("n1", x1, y, NodeKind::Junction),
        ("n2", x2, y, NodeKind::Junction),
        ("n3", x3, y, NodeKind::Junction),
        ("n4", x4, y, NodeKind::Outlet),
    ] {
        system.add_node(NodeSpec::new(id, kind).with_point((px, py)));
    }

    let height_m = 1.0 / 1000.0;
    for (id, from, to, width_mm) in [
        ("ch0", "n0", "n1", w_inlet),
        ("ch1", "n1", "n2", w_mid),
        ("ch2", "n2", "n3", w_throat),
        ("ch3", "n3", "n4", w_mid),
    ] {
        let length_m = (width_mm * std::f64::consts::SQRT_2) / 1000.0;
        system.add_channel(ChannelSpec::new_pipe_rect(
            id,
            from,
            to,
            length_m,
            width_mm / 1000.0,
            height_m,
            1.0,
            0.0,
        ));
    }

    let renderer = create_plotters_renderer();
    let schematic_cfg = RenderConfig {
        width: 1200,
        height: 400,
        title: "Venturi — ISO 5167 Topology".to_string(),
        ..Default::default()
    };

    let schematic_path = out.join("venturi_schematic.png");
    renderer.render_system(&system, schematic_path.to_str().unwrap(), &schematic_cfg)?;
    println!("  ✅ Schematic → {}", schematic_path.display());

    // ── Phase 2: Simulation (cfd-2d) ─────────────────────────────────────────
    println!("\n⚙️  Phase 2: 2D CFD Simulation");

    let geom = VenturiGeometry::<f64>::iso_5167_standard();
    let u_inlet = 0.1_f64; // 100 mm/s
    let density = 1060.0_f64; // kg/m³ (blood)

    println!("  Inlet width:  {:.1} mm", geom.w_inlet * 1000.0);
    println!("  Throat width: {:.2} mm", geom.w_throat * 1000.0);
    println!("  Area ratio:   {:.3}", geom.area_ratio());

    let mut solver =
        VenturiSolver2D::new(geom.clone(), BloodModel::Newtonian(1.0e-3), density, 60, 30);
    let sol = solver.solve(u_inlet)?;

    println!("\n  📊 Results:");
    println!("  u_inlet  = {:.4} m/s", sol.u_inlet);
    println!(
        "  u_throat = {:.4} m/s  ({:.1}× acceleration)",
        sol.u_throat,
        sol.u_throat / sol.u_inlet.max(1e-12)
    );
    println!("  p_inlet  = {:.2} Pa", sol.p_inlet);
    println!(
        "  p_throat = {:.2} Pa  (ΔP = {:.2} Pa)",
        sol.p_throat, sol.dp_throat
    );
    println!("  p_outlet = {:.2} Pa", sol.p_outlet);
    println!("  Cp_throat = {:.4}", sol.cp_throat);

    // Bernoulli validation
    let bernoulli = BernoulliVenturi::new(geom, u_inlet, sol.p_inlet, density);
    let u_err = (sol.u_throat - bernoulli.velocity_throat()).abs()
        / bernoulli.velocity_throat().abs()
        * 100.0;
    let p_err = (sol.p_throat - bernoulli.pressure_throat()).abs()
        / bernoulli.pressure_throat().abs().max(1.0)
        * 100.0;
    println!("\n  🔍 Bernoulli Validation:");
    println!(
        "  u_throat error = {:.1}%  (ideal = {:.4} m/s)",
        u_err,
        bernoulli.velocity_throat()
    );
    println!(
        "  p_throat error = {:.1}%  (ideal = {:.2} Pa)",
        p_err,
        bernoulli.pressure_throat()
    );

    // ── Phase 3: Overlay Visualization ───────────────────────────────────────
    println!("\n🎨 Phase 3: Pressure Overlay");

    // Assign pressure values to each node (0=inlet … 4=outlet)
    let node_p: HashMap<usize, f64> = [
        (0, sol.p_inlet),
        (1, sol.p_inlet * 0.99),
        (2, sol.p_throat),
        (3, sol.p_throat + (sol.p_outlet - sol.p_throat) * 0.5),
        (4, sol.p_outlet),
    ]
    .into_iter()
    .collect();

    // Assign mean pressure to each channel edge
    let edge_p: HashMap<usize, f64> = [
        (0, sol.p_inlet),
        (1, (sol.p_inlet + sol.p_throat) / 2.0),
        (2, sol.p_throat),
        (3, (sol.p_throat + sol.p_outlet) / 2.0),
    ]
    .into_iter()
    .collect();

    let overlay = AnalysisOverlay::new(AnalysisField::Pressure, ColormapKind::BlueRed)
        .with_node_data(node_p)
        .with_edge_data(edge_p);

    let overlay_cfg = RenderConfig {
        width: 1200,
        height: 400,
        title: "Venturi — Pressure Distribution [Pa]".to_string(),
        ..Default::default()
    };

    let overlay_path = out.join("venturi_pressure_overlay.png");
    renderer.render_analysis(
        &system,
        overlay_path.to_str().unwrap(),
        &overlay_cfg,
        &overlay,
    )?;
    println!("  ✅ Pressure overlay → {}", overlay_path.display());

    println!("\n✅ Venturi example complete.");
    Ok(())
}
