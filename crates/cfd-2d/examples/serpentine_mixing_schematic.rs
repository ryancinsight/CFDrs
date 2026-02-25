//! Serpentine Mixing Schematic + 2D CFD Example
//!
//! Demonstrates the two-phase schematic-driven pattern for cfd-2d:
//!
//! **Phase 1 — Design (cfd-schematics)**
//! - Define serpentine topology as a linear chain of `ChannelSpec` segments
//! - Render schematic PNG → `outputs/serpentine_schematic.png`
//!
//! **Phase 2 — Simulation (cfd-2d)**
//! - Run `SerpentineSolver2D` (microfluidic standard geometry)
//! - Map pressure drop per segment to `AnalysisOverlay` (BlueRed colormap)
//! - Render overlay PNG → `outputs/serpentine_mixing_overlay.png`
//!
//! Run with:
//! `cargo run -p cfd-2d --example serpentine_mixing_schematic`

use cfd_2d::solvers::serpentine_flow::{
    AdvectionDiffusionMixing, SerpentineGeometry, SerpentineSolver2D,
};
use cfd_2d::solvers::ns_fvm_2d::BloodModel;
use cfd_schematics::geometry::types::{Channel, ChannelSystem, ChannelType, Node};
use cfd_schematics::visualizations::{
    create_plotters_renderer, AnalysisField, AnalysisOverlay, ColormapKind, RenderConfig,
    SchematicRenderer,
};
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🔬 Serpentine Mixing Schematic + 2D CFD Example");
    println!("================================================");

    // ── 0. Output directory ───────────────────────────────────────────────────
    let out = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("outputs")
        .join("serpentine_mixing");
    fs::create_dir_all(&out)?;

    // ── Phase 1: Design (cfd-schematics) ─────────────────────────────────────
    println!("\n📐 Phase 1: Schematic Design");

    // Microfluidic standard serpentine: 5 cycles, 200 μm wide, 500 μm straight sections
    let n_cycles = 5_usize;
    let w_mm = 0.2_f64;      // 200 μm channel width in mm
    let l_straight_mm = 0.5; // 500 μm straight section in mm
    let r_mm = 0.2_f64;      // 200 μm turn radius in mm

    // Schematic: lay out as a horizontal chain of 2*n_cycles straight segments
    // Each segment represents one straight section of the serpentine
    // Nodes alternate top/bottom to suggest the serpentine path
    let n_segments = 2 * n_cycles;
    let seg_len = l_straight_mm;
    let y_top = 2.0_f64;
    let y_bot = 0.0_f64;

    let mut nodes: Vec<Node> = Vec::with_capacity(n_segments + 1);
    let mut channels: Vec<Channel> = Vec::with_capacity(n_segments);

    for i in 0..=n_segments {
        let x = i as f64 * seg_len;
        let y = if i % 2 == 0 { y_bot } else { y_top };
        nodes.push(Node { id: i, point: (x, y), metadata: None });
    }

    for i in 0..n_segments {
        channels.push(Channel {
            id: i,
            from_node: i,
            to_node: i + 1,
            width: w_mm,
            height: 0.05, // 50 μm depth
            channel_type: ChannelType::Straight,
            metadata: None,
        });
    }

    let total_x = n_segments as f64 * seg_len + 1.0;
    let system = ChannelSystem {
        box_dims: (total_x, y_top + 1.0),
        nodes,
        channels,
        box_outline: vec![],
    };

    let renderer = create_plotters_renderer();
    let schematic_cfg = RenderConfig {
        width: 1400,
        height: 400,
        title: format!(
            "Serpentine Mixer — {n_cycles} cycles, w={:.0}μm, L_straight={:.0}μm",
            w_mm * 1000.0,
            l_straight_mm * 1000.0
        ),
        ..Default::default()
    };

    let schematic_path = out.join("serpentine_schematic.png");
    renderer.render_system(&system, schematic_path.to_str().unwrap(), &schematic_cfg)?;
    println!("  ✅ Schematic → {}", schematic_path.display());

    // ── Phase 2: Simulation (cfd-2d) ─────────────────────────────────────────
    println!("\n⚙️  Phase 2: 2D CFD Simulation");

    let geom = SerpentineGeometry::<f64>::microfluidic_standard();
    let total_length = geom.total_length();
    println!("  Geometry: microfluidic_standard()");
    println!("  Width:          {:.0} μm", geom.width * 1e6);
    println!("  Straight length:{:.0} μm", geom.l_straight * 1e6);
    println!("  Turn radius:    {:.0} μm", geom.turn_radius * 1e6);
    println!("  Cycles:         {}", geom.n_cycles);
    println!("  Total length:   {:.2} mm", total_length * 1000.0);

    let u_inlet = 0.01_f64;           // 10 mm/s
    let diffusion_coeff = 1e-9_f64;   // 1 nm²/s — typical small molecule in water
    let density = 1000.0_f64;
    let viscosity = 1e-3_f64;         // 1 mPa·s (water)

    // Analytical mixing model for design insight
    let mixing_model = AdvectionDiffusionMixing::new(geom.width, u_inlet, diffusion_coeff);
    let pe = mixing_model.peclet_number();
    let l_mix = mixing_model.mixing_length_90_percent();
    let mixing_frac_analytical = mixing_model.mixing_fraction(total_length);

    println!("\n  📐 Analytical Mixing Estimates:");
    println!("  Peclet number:    {:.1e}", pe);
    println!("  Mixing length:    {:.2} mm  (90% mixing)", l_mix * 1000.0);
    println!("  Mixing fraction:  {:.1}%  at outlet (analytical)", mixing_frac_analytical * 100.0);
    println!("  Well-mixed:       {}", if mixing_frac_analytical > 0.9 { "YES ✅" } else { "NO ❌" });

    println!("\n  Running SerpentineSolver2D (this may take a moment)...");
    let mut solver = SerpentineSolver2D::new(
        geom.clone(),
        BloodModel::Newtonian(viscosity),
        density,
        60, 30,
    );

    let sol = solver.solve(u_inlet, diffusion_coeff, 0.0, 1.0)?;

    println!("\n  📊 Simulation Results:");
    println!("  Peclet number:    {:.1e}", sol.peclet);
    println!("  Mixing fraction:  {:.1}%  at outlet (numerical)", sol.mixing_fraction_outlet * 100.0);
    println!("  Pressure drop:    {:.2} Pa", sol.pressure_drop);
    println!("  Well-mixed:       {}", if sol.is_well_mixed() { "YES ✅" } else { "NO ❌" });

    // Pressure drop per segment (estimated from total)
    let dp_per_segment = sol.pressure_drop / n_segments as f64;
    println!("  ΔP per segment:   {:.3} Pa", dp_per_segment);

    // ── Phase 3: Overlay Visualization ───────────────────────────────────────
    println!("\n🎨 Phase 3: Pressure Drop Overlay");

    // Assign cumulative pressure drop to each node (linear distribution)
    let p_inlet = sol.pressure_drop; // inlet is high pressure
    let node_p: HashMap<usize, f64> = (0..=n_segments)
        .map(|i| {
            let frac = i as f64 / n_segments as f64;
            (i, p_inlet * (1.0 - frac))
        })
        .collect();

    // Assign pressure drop per channel edge (uniform)
    let edge_p: HashMap<usize, f64> = (0..n_segments)
        .map(|i| {
            let frac = (i as f64 + 0.5) / n_segments as f64;
            (i, p_inlet * (1.0 - frac))
        })
        .collect();

    let overlay = AnalysisOverlay::new(AnalysisField::Pressure, ColormapKind::BlueRed)
        .with_node_data(node_p)
        .with_edge_data(edge_p);

    let overlay_cfg = RenderConfig {
        width: 1400,
        height: 400,
        title: format!(
            "Serpentine — Pressure Distribution [Pa]  (ΔP_total = {:.2} Pa, mixing = {:.0}%)",
            sol.pressure_drop,
            sol.mixing_fraction_outlet * 100.0
        ),
        ..Default::default()
    };

    let overlay_path = out.join("serpentine_mixing_overlay.png");
    renderer.render_analysis(&system, overlay_path.to_str().unwrap(), &overlay_cfg, &overlay)?;
    println!("  ✅ Pressure overlay → {}", overlay_path.display());

    println!("\n✅ Serpentine mixing example complete.");
    Ok(())
}
