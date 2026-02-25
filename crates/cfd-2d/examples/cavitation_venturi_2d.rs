//! Cavitation Risk Analysis — 2D Venturi Solver
//!
//! 2D counterpart of `cfd-1d::cavitation_venturi_analysis`. Replaces the
//! lumped-resistance model with a full 2D Navier-Stokes simulation of the
//! Venturi constriction, providing spatially resolved velocity and pressure
//! fields for accurate cavitation assessment.
//!
//! Pipeline:
//! 1. Build Venturi schematic (cfd-schematics `ChannelSystem`)
//! 2. Run `VenturiSolver2D` with water (Newtonian)
//! 3. Validate against Bernoulli analytical solution
//! 4. Compute cavitation numbers using `VenturiCavitation` from cfd-core
//! 5. Render pressure, velocity, cavitation overlays + export JSON
//!
//! Run with:
//! `cargo run -p cfd-2d --example cavitation_venturi_2d`

use cfd_2d::solvers::venturi_flow::{BernoulliVenturi, VenturiGeometry, VenturiSolver2D};
use cfd_2d::solvers::ns_fvm_2d::BloodModel;
use cfd_core::physics::cavitation::{CavitationNumber, VenturiCavitation};
use cfd_schematics::geometry::types::{Channel, ChannelSystem, ChannelType, Node};
use cfd_schematics::visualizations::{
    create_plotters_renderer, AnalysisField, AnalysisOverlay, ColormapKind, RenderConfig,
    SchematicRenderer,
};
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("💧 Cavitation Risk Analysis — 2D Venturi Solver");
    println!("================================================\n");

    // ── 0. Output directory ──────────────────────────────────────────────────
    let out = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("outputs")
        .join("cavitation_venturi");
    fs::create_dir_all(&out)?;

    // ── Phase 1: Design (cfd-schematics) ─────────────────────────────────────
    println!("📐 Phase 1: Schematic Design\n");

    // Millifluidic Venturi sized for 96-well plate SDT device
    // inlet(10mm) → converge(2mm) → throat(3mm) → diverge(5mm)
    let w_inlet_mm = 2.0_f64;   // 2 mm inlet width
    let w_throat_mm = 0.8_f64;  // 800 µm throat — induces high velocity
    let w_mid = (w_inlet_mm + w_throat_mm) / 2.0;

    let x0 = 0.0_f64;
    let x1 = 10.0;  // end of inlet section
    let x2 = 12.0;  // end of converging section
    let x3 = 15.0;  // end of throat section
    let x4 = 20.0;  // end of diverging section
    let y = 3.0;     // vertical centre

    let nodes: Vec<Node> = vec![
        Node { id: 0, point: (x0, y), metadata: None },
        Node { id: 1, point: (x1, y), metadata: None },
        Node { id: 2, point: (x2, y), metadata: None },
        Node { id: 3, point: (x3, y), metadata: None },
        Node { id: 4, point: (x4, y), metadata: None },
    ];

    let channels: Vec<Channel> = vec![
        Channel { id: 0, from_node: 0, to_node: 1, width: w_inlet_mm,  height: 1.0, channel_type: ChannelType::Straight, metadata: None },
        Channel { id: 1, from_node: 1, to_node: 2, width: w_mid,       height: 1.0, channel_type: ChannelType::Straight, metadata: None },
        Channel { id: 2, from_node: 2, to_node: 3, width: w_throat_mm, height: 1.0, channel_type: ChannelType::Straight, metadata: None },
        Channel { id: 3, from_node: 3, to_node: 4, width: w_mid,       height: 1.0, channel_type: ChannelType::Straight, metadata: None },
    ];

    let system = ChannelSystem {
        box_dims: (x4 + 2.0, y * 2.0),
        nodes,
        channels,
        box_outline: vec![],
    };

    let renderer = create_plotters_renderer();
    let schematic_cfg = RenderConfig {
        width: 1200,
        height: 400,
        title: format!(
            "Venturi — Millifluidic (inlet {:.1}mm, throat {:.1}mm)",
            w_inlet_mm, w_throat_mm
        ),
        ..Default::default()
    };

    let schematic_path = out.join("schematic.png");
    renderer.render_system(&system, schematic_path.to_str().unwrap(), &schematic_cfg)?;
    println!("  ✅ Schematic → {}", schematic_path.display());

    // ── Phase 2: 2D CFD Simulation ───────────────────────────────────────────
    println!("\n⚙️  Phase 2: 2D Navier-Stokes Simulation\n");

    // Solver geometry in metres
    let geom = VenturiGeometry::<f64>::new(
        w_inlet_mm * 1e-3,   // inlet width [m]
        w_throat_mm * 1e-3,  // throat width [m]
        10.0e-3,             // inlet section length [m]
        2.0e-3,              // converging section [m]
        3.0e-3,              // throat section [m]
        5.0e-3,              // diverging section [m]
        1.0e-3,              // height [m]
    );

    let u_inlet = 0.2_f64;     // 200 mm/s — aggressive for cavitation study
    let density = 997.0_f64;   // water at 25°C
    let viscosity = 8.9e-4;    // Pa·s
    let vapor_pressure = 3169.0; // Pa at 25°C

    println!("  Geometry:");
    println!("    Inlet width:   {:.1} mm", w_inlet_mm);
    println!("    Throat width:  {:.1} mm", w_throat_mm);
    println!("    Area ratio:    {:.3}", geom.area_ratio());
    println!("    Total length:  {:.1} mm", geom.total_length() * 1e3);
    println!("  Fluid: Water at 25°C, ρ = {:.0} kg/m³, µ = {:.1e} Pa·s", density, viscosity);
    println!("  Inlet velocity:  {:.0} mm/s", u_inlet * 1e3);

    let mut solver = VenturiSolver2D::new(
        geom.clone(),
        BloodModel::Newtonian(viscosity),
        density,
        80, 40,
    );
    let sol = solver.solve(u_inlet)?;

    println!("\n  📊 2D Simulation Results:");
    println!("    u_inlet   = {:.4} m/s", sol.u_inlet);
    println!("    u_throat  = {:.4} m/s  ({:.1}× acceleration)", sol.u_throat, sol.u_throat / sol.u_inlet.max(1e-12));
    println!("    p_inlet   = {:.2} Pa", sol.p_inlet);
    println!("    p_throat  = {:.2} Pa  (ΔP = {:.2} Pa)", sol.p_throat, sol.dp_throat);
    println!("    p_outlet  = {:.2} Pa", sol.p_outlet);
    println!("    Cp_throat = {:.4}", sol.cp_throat);

    // Bernoulli validation
    let bernoulli = BernoulliVenturi::new(geom.clone(), u_inlet, sol.p_inlet, density);
    let u_err = (sol.u_throat - bernoulli.velocity_throat()).abs()
        / bernoulli.velocity_throat().abs().max(1e-12)
        * 100.0;
    let p_err = (sol.p_throat - bernoulli.pressure_throat()).abs()
        / bernoulli.pressure_throat().abs().max(1.0)
        * 100.0;

    println!("\n  🔍 Bernoulli Validation:");
    println!("    u_throat error = {:.1}%  (ideal = {:.4} m/s)", u_err, bernoulli.velocity_throat());
    println!("    p_throat error = {:.1}%  (ideal = {:.2} Pa)", p_err, bernoulli.pressure_throat());

    // ── Phase 3: Cavitation Analysis ─────────────────────────────────────────
    println!("\n🫧 Phase 3: Cavitation Analysis\n");

    // Use VenturiCavitation from cfd-core for detailed cavitation assessment
    // Convert widths to "diameters" for VenturiCavitation (which uses circular cross-section)
    let cav = VenturiCavitation {
        inlet_diameter: w_inlet_mm * 1e-3,
        throat_diameter: w_throat_mm * 1e-3,
        outlet_diameter: w_inlet_mm * 1e-3,
        convergent_angle: 30.0_f64.to_radians(),
        divergent_angle: 15.0_f64.to_radians(),
        inlet_pressure: sol.p_inlet,
        inlet_velocity: sol.u_inlet,
        density,
        vapor_pressure,
    };

    let sigma = cav.cavitation_number();
    let is_cavitating = cav.is_cavitating();
    let throat_v = cav.throat_velocity();
    let throat_p = cav.throat_pressure();

    // Classify regime
    let regime = if sigma > 1.5 {
        "No Cavitation"
    } else if sigma > 0.5 {
        "Inception"
    } else if sigma > 0.1 {
        "Developed"
    } else {
        "Supercavitation"
    };

    println!("  Cavitation number σ = {:.4}", sigma);
    println!("  Is cavitating:       {}", if is_cavitating { "YES ⚠" } else { "NO ✅" });
    println!("  Regime:              {regime}");
    println!("  Throat velocity:     {:.4} m/s", throat_v);
    println!("  Throat pressure:     {:.2} Pa", throat_p);

    if is_cavitating {
        let cavity_len = cav.cavity_length(sigma);
        println!("  Cavity length:       {:.2} mm", cavity_len * 1e3);
    }

    // Per-section cavitation data for overlay
    let section_data: Vec<(usize, f64, f64, f64)> = vec![
        (0, sol.u_inlet,  sol.p_inlet,  u_inlet),
        (1, (sol.u_inlet + sol.u_throat) / 2.0, (sol.p_inlet + sol.p_throat) / 2.0, (sol.u_inlet + sol.u_throat) / 2.0),
        (2, sol.u_throat, sol.p_throat, sol.u_throat),
        (3, (sol.u_throat + sol.u_outlet) / 2.0, (sol.p_throat + sol.p_outlet) / 2.0, (sol.u_throat + sol.u_outlet) / 2.0),
    ];

    let mut edge_sigma = HashMap::<usize, f64>::new();
    let mut edge_velocity = HashMap::<usize, f64>::new();
    let mut edge_pressure = HashMap::<usize, f64>::new();
    let mut node_pressure = HashMap::<usize, f64>::new();

    node_pressure.insert(0, sol.p_inlet);
    node_pressure.insert(1, sol.p_inlet * 0.99);
    node_pressure.insert(2, sol.p_throat);
    node_pressure.insert(3, (sol.p_throat + sol.p_outlet) / 2.0);
    node_pressure.insert(4, sol.p_outlet);

    for &(sec_id, _u, p, v) in &section_data {
        let sec_cav = CavitationNumber {
            reference_pressure: p,
            vapor_pressure,
            density,
            velocity: v,
        };
        let sec_sigma = sec_cav.calculate();
        edge_sigma.insert(sec_id, sec_sigma);
        edge_velocity.insert(sec_id, v);
        edge_pressure.insert(sec_id, p);
    }

    // ── Phase 4: Overlay Visualization ───────────────────────────────────────
    println!("\n🎨 Phase 4: Overlay Visualization\n");

    // Cavitation number overlay (inverted: low σ = red = high risk)
    let max_sigma = edge_sigma.values().copied().fold(0.0_f64, f64::max);
    let inv_sigma: HashMap<usize, f64> = edge_sigma
        .iter()
        .map(|(&k, &v)| (k, max_sigma - v))
        .collect();

    let overlays: Vec<(&str, AnalysisOverlay, String)> = vec![
        (
            "cavitation_number.png",
            AnalysisOverlay::new(
                AnalysisField::Custom("Cavitation Number σ".into()),
                ColormapKind::BlueRed,
            )
            .with_edge_data(inv_sigma)
            .with_node_data(node_pressure.clone()),
            format!("Cavitation Number σ (red = low σ, regime: {regime})"),
        ),
        (
            "pressure.png",
            AnalysisOverlay::new(AnalysisField::Pressure, ColormapKind::Viridis)
                .with_edge_data(edge_pressure.clone())
                .with_node_data(node_pressure.clone()),
            format!("Pressure Distribution (ΔP_throat = {:.2} Pa)", sol.dp_throat),
        ),
        (
            "velocity.png",
            AnalysisOverlay::new(AnalysisField::Velocity, ColormapKind::Viridis)
                .with_edge_data(edge_velocity.clone())
                .with_node_data(node_pressure.clone()),
            format!("Velocity Distribution (u_throat = {:.4} m/s)", sol.u_throat),
        ),
    ];

    for (filename, overlay, title) in &overlays {
        let path = out.join(filename);
        let cfg = RenderConfig {
            width: 1200,
            height: 400,
            title: title.clone(),
            show_axes: true,
            show_grid: false,
            ..Default::default()
        };
        renderer.render_analysis(&system, path.to_str().unwrap(), &cfg, overlay)?;
        println!("  ✓ {filename}");
    }

    // ── Phase 5: JSON Export ─────────────────────────────────────────────────
    println!("\n📄 Phase 5: JSON Export\n");

    let section_json: Vec<serde_json::Value> = section_data
        .iter()
        .map(|&(id, _u, p, v)| {
            serde_json::json!({
                "section_id": id,
                "velocity_m_s": v,
                "pressure_pa": p,
                "cavitation_number": edge_sigma.get(&id),
            })
        })
        .collect();

    let results = serde_json::json!({
        "analysis": "cavitation_venturi_2d",
        "solver": "VenturiSolver2D (Navier-Stokes FVM)",
        "geometry": {
            "inlet_width_mm": w_inlet_mm,
            "throat_width_mm": w_throat_mm,
            "area_ratio": geom.area_ratio(),
            "total_length_mm": geom.total_length() * 1e3
        },
        "fluid": {
            "name": "Water (25°C)",
            "density_kg_m3": density,
            "viscosity_pa_s": viscosity,
            "vapor_pressure_pa": vapor_pressure
        },
        "boundary_conditions": {
            "inlet_velocity_m_s": u_inlet
        },
        "simulation_results": {
            "u_inlet_m_s": sol.u_inlet,
            "u_throat_m_s": sol.u_throat,
            "u_outlet_m_s": sol.u_outlet,
            "p_inlet_pa": sol.p_inlet,
            "p_throat_pa": sol.p_throat,
            "p_outlet_pa": sol.p_outlet,
            "dp_throat_pa": sol.dp_throat,
            "cp_throat": sol.cp_throat,
            "cp_recovery": sol.cp_recovery
        },
        "bernoulli_validation": {
            "u_throat_error_percent": u_err,
            "p_throat_error_percent": p_err
        },
        "cavitation_analysis": {
            "cavitation_number": sigma,
            "is_cavitating": is_cavitating,
            "regime": regime,
            "throat_velocity_m_s": throat_v,
            "throat_pressure_pa": throat_p
        },
        "section_data": section_json
    });

    let json_path = out.join("results.json");
    fs::write(&json_path, serde_json::to_string_pretty(&results)?)?;
    println!("  ✅ JSON → {}", json_path.display());

    println!("\n✅ Cavitation venturi 2D example complete.");
    Ok(())
}
