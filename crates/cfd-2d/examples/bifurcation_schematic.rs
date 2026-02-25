//! Bifurcation Schematic + 2D CFD Example
//!
//! Demonstrates the two-phase schematic-driven pattern for cfd-2d:
//!
//! **Phase 1 — Design (cfd-schematics)**
//! - Define symmetric bifurcation topology as `ChannelSystem`
//! - Murray's law sizing: r_d = r_p / 2^(1/3)
//! - Render schematic PNG → `outputs/bifurcation_schematic.png`
//!
//! **Phase 2 — Simulation (cfd-2d)**
//! - Run `BifurcationSolver2D` (Casson blood, 60×40 grid)
//! - Map flow rate results to `AnalysisOverlay` (Viridis colormap)
//! - Render overlay PNG → `outputs/bifurcation_flow_overlay.png`
//!
//! Run with:
//! `cargo run -p cfd-2d --example bifurcation_schematic`

use cfd_2d::solvers::bifurcation_flow::{BifurcationGeometry, BifurcationSolver2D};
use cfd_2d::solvers::ns_fvm_2d::{BloodModel, SIMPLEConfig};
use cfd_core::physics::fluid::blood::CassonBlood;
use cfd_schematics::geometry::types::{Channel, ChannelSystem, ChannelType, Node};
use cfd_schematics::visualizations::{
    create_plotters_renderer, AnalysisField, AnalysisOverlay, ColormapKind, RenderConfig,
    SchematicRenderer,
};
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🔬 Bifurcation Schematic + 2D CFD Example");
    println!("==========================================");

    // ── 0. Output directory ───────────────────────────────────────────────────
    let out = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("outputs")
        .join("bifurcation");
    fs::create_dir_all(&out)?;

    // ── Phase 1: Design (cfd-schematics) ─────────────────────────────────────
    println!("\n📐 Phase 1: Schematic Design");

    // Murray's law: r_parent³ = r_d1³ + r_d2³ → r_d = r_p / 2^(1/3)
    // Using 1 mm parent radius → 0.794 mm daughter radius
    let r_parent = 1.0_f64;           // mm
    let r_daughter = r_parent / 2.0_f64.powf(1.0 / 3.0); // ≈ 0.794 mm
    let w_parent   = r_parent * 2.0;
    let w_daughter = r_daughter * 2.0;

    // Layout: parent runs left→right, daughters branch up and down at 30°
    // Schematic coordinates in mm
    let angle_deg = 30.0_f64;
    let angle_rad = angle_deg.to_radians();
    let l_parent   = 4.0_f64;  // mm
    let l_daughter = 4.0_f64;  // mm

    let jx = l_parent;
    let jy = 5.0_f64;  // vertical centre

    let d1_end_x = jx + l_daughter * angle_rad.cos();
    let d1_end_y = jy + l_daughter * angle_rad.sin();
    let d2_end_x = jx + l_daughter * angle_rad.cos();
    let d2_end_y = jy - l_daughter * angle_rad.sin();

    let nodes: Vec<Node> = vec![
        Node { id: 0, point: (0.0, jy),       metadata: None }, // inlet
        Node { id: 1, point: (jx,  jy),       metadata: None }, // junction
        Node { id: 2, point: (d1_end_x, d1_end_y), metadata: None }, // daughter 1
        Node { id: 3, point: (d2_end_x, d2_end_y), metadata: None }, // daughter 2
    ];

    let channels: Vec<Channel> = vec![
        Channel { id: 0, from_node: 0, to_node: 1, width: w_parent,   height: 1.0, channel_type: ChannelType::Straight, metadata: None },
        Channel { id: 1, from_node: 1, to_node: 2, width: w_daughter, height: 1.0, channel_type: ChannelType::Straight, metadata: None },
        Channel { id: 2, from_node: 1, to_node: 3, width: w_daughter, height: 1.0, channel_type: ChannelType::Straight, metadata: None },
    ];

    let system = ChannelSystem {
        box_dims: (d1_end_x + 2.0, jy * 2.0),
        nodes,
        channels,
        box_outline: vec![],
    };

    let renderer = create_plotters_renderer();
    let schematic_cfg = RenderConfig {
        width: 900,
        height: 600,
        title: format!(
            "Bifurcation — Murray's Law (r_p={:.1}mm, r_d={:.2}mm, θ={:.0}°)",
            r_parent, r_daughter, angle_deg
        ),
        ..Default::default()
    };

    let schematic_path = out.join("bifurcation_schematic.png");
    renderer.render_system(&system, schematic_path.to_str().unwrap(), &schematic_cfg)?;
    println!("  ✅ Schematic → {}", schematic_path.display());
    println!("  Murray's law: r_parent={:.3} mm, r_daughter={:.3} mm", r_parent, r_daughter);
    println!("  Verification: r_p³ = {:.4}, r_d1³+r_d2³ = {:.4}",
        r_parent.powi(3),
        r_daughter.powi(3) * 2.0
    );

    // ── Phase 2: Simulation (cfd-2d) ─────────────────────────────────────────
    println!("\n⚙️  Phase 2: 2D CFD Simulation");

    // Convert mm to metres for the solver
    let geom = BifurcationGeometry::new_symmetric(
        w_parent   * 1e-3,   // parent width [m]
        l_parent   * 1e-3,   // parent length [m]
        w_daughter * 1e-3,   // daughter width [m]
        l_daughter * 1e-3,   // daughter length [m]
        angle_rad,           // branch angle [rad]
    );

    // Casson blood model
    let casson = CassonBlood::<f64>::normal_blood();
    let density = 1060.0_f64;
    let blood = BloodModel::Casson(casson);
    let u_inlet = 0.1_f64; // 100 mm/s

    let mut config = SIMPLEConfig::default();
    config.max_iterations = 3000;
    config.tolerance = 1e-4;
    config.alpha_u = 0.5;
    config.alpha_p = 0.2;

    println!("  Parent width:   {:.2} mm", w_parent);
    println!("  Daughter width: {:.3} mm  (Murray's law)", w_daughter);
    println!("  Inlet velocity: {:.0} mm/s", u_inlet * 1000.0);
    println!("  Fluid: Casson blood, ρ={:.0} kg/m³", density);

    let mut solver = BifurcationSolver2D::new(geom, blood, density, 60, 40, config);
    let sol = solver.solve(u_inlet)?;

    println!("\n  📊 Results:");
    println!("  Q_parent    = {:.4e} m²/s", sol.q_parent);
    println!("  Q_daughter1 = {:.4e} m²/s", sol.q_daughter1);
    println!("  Q_daughter2 = {:.4e} m²/s", sol.q_daughter2);
    println!("  Mass balance error = {:.2e} ({:.2}%)",
        sol.mass_balance_error,
        sol.mass_balance_error * 100.0
    );

    let symmetry_err = (sol.q_daughter1 - sol.q_daughter2).abs() / sol.q_parent.abs().max(1e-20);
    println!("  Symmetry error = {:.2e} ({:.2}%)", symmetry_err, symmetry_err * 100.0);

    // ── Phase 3: Overlay Visualization ───────────────────────────────────────
    println!("\n🎨 Phase 3: Flow Rate Overlay");

    // Map flow rates to edges
    let edge_q: HashMap<usize, f64> = [
        (0, sol.q_parent),
        (1, sol.q_daughter1),
        (2, sol.q_daughter2),
    ].into_iter().collect();

    // Map derived pressure-like quantity to nodes (proportional to flow)
    let node_q: HashMap<usize, f64> = [
        (0, sol.q_parent),
        (1, sol.q_parent),
        (2, sol.q_daughter1),
        (3, sol.q_daughter2),
    ].into_iter().collect();

    let overlay = AnalysisOverlay::new(AnalysisField::FlowRate, ColormapKind::Viridis)
        .with_edge_data(edge_q)
        .with_node_data(node_q);

    let overlay_cfg = RenderConfig {
        width: 900,
        height: 600,
        title: "Bifurcation — Flow Rate Distribution [m²/s]".to_string(),
        ..Default::default()
    };

    let overlay_path = out.join("bifurcation_flow_overlay.png");
    renderer.render_analysis(&system, overlay_path.to_str().unwrap(), &overlay_cfg, &overlay)?;
    println!("  ✅ Flow rate overlay → {}", overlay_path.display());

    println!("\n✅ Bifurcation example complete.");
    Ok(())
}
