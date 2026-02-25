//! Hemolysis Analysis in a Serpentine Channel — 2D Navier-Stokes
//!
//! 2D counterpart of `cfd-1d::hemolysis_serpentine_analysis`. Instead of a
//! lumped 1D network, this example resolves the full 2D velocity and pressure
//! fields inside a serpentine micromixer using the `SerpentineSolver2D` Navier-
//! Stokes solver, then back-projects the results onto a `cfd-schematics`
//! overlay for medical-grade hemolysis assessment.
//!
//! Pipeline:
//! 1. Build serpentine schematic (cfd-schematics `ChannelSystem`)
//! 2. Run `SerpentineSolver2D` with Carreau-Yasuda blood
//! 3. Estimate per-segment wall shear stress from 2D pressure drop
//! 4. Compute Giersiepen–Wurzinger hemolysis index per segment
//! 5. Render hemolysis, WSS, pressure overlays + export JSON
//!
//! Run with:
//! `cargo run -p cfd-2d --example hemolysis_serpentine_2d`

use cfd_2d::solvers::serpentine_flow::{AdvectionDiffusionMixing, SerpentineGeometry, SerpentineSolver2D};
use cfd_2d::solvers::ns_fvm_2d::BloodModel;
use cfd_core::physics::cavitation::CavitationNumber;
use cfd_core::physics::fluid::blood::CarreauYasudaBlood;
use cfd_core::physics::hemolysis::HemolysisModel;
use cfd_schematics::geometry::types::{Channel, ChannelSystem, ChannelType, Node};
use cfd_schematics::visualizations::{
    create_plotters_renderer, AnalysisField, AnalysisOverlay, ColormapKind, RenderConfig,
    SchematicRenderer,
};
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🩸 Hemolysis Analysis — 2D Serpentine Solver");
    println!("=============================================\n");

    // ── 0. Output directory ──────────────────────────────────────────────────
    let out = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("outputs")
        .join("hemolysis_serpentine");
    fs::create_dir_all(&out)?;

    // ── Phase 1: Design (cfd-schematics) ─────────────────────────────────────
    println!("📐 Phase 1: Schematic Design\n");

    // Millifluidic serpentine: 4 cycles, 800 µm wide, 2 mm straight sections
    let n_cycles = 4_usize;
    let w_mm = 0.8_f64;          // 800 µm channel width (mm for schematic)
    let l_straight_mm = 2.0;     // 2 mm straight section
    let r_mm = 0.5_f64;          // 500 µm turn radius

    // Schematic layout: zigzag chain of 2*n_cycles segments
    let n_segments = 2 * n_cycles;
    let seg_len = l_straight_mm;
    let y_top = 3.0_f64;
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
            height: 0.5,
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
        height: 500,
        title: format!(
            "Serpentine — {n_cycles} cycles, w={:.0}µm, L_straight={:.1}mm (blood)",
            w_mm * 1000.0,
            l_straight_mm
        ),
        ..Default::default()
    };

    let schematic_path = out.join("schematic.png");
    renderer.render_system(&system, schematic_path.to_str().unwrap(), &schematic_cfg)?;
    println!("  ✅ Schematic → {}", schematic_path.display());

    // ── Phase 2: 2D CFD Simulation ───────────────────────────────────────────
    println!("\n⚙️  Phase 2: 2D Navier-Stokes Simulation\n");

    // Convert mm to metres for the solver
    let w_m = w_mm * 1e-3;
    let l_straight_m = l_straight_mm * 1e-3;
    let r_m = r_mm * 1e-3;
    let h_m = 0.5e-3; // 500 µm height

    let geom = SerpentineGeometry::new(w_m, h_m, l_straight_m, r_m, n_cycles);
    let total_length = geom.total_length();

    let blood = CarreauYasudaBlood::<f64>::normal_blood();
    let density = blood.density;

    // Use Carreau-Yasuda model for non-Newtonian blood
    let blood_model = BloodModel::CarreauYasuda(blood);

    println!("  Geometry:");
    println!("    Width:           {:.0} µm", w_m * 1e6);
    println!("    Straight:        {:.1} mm", l_straight_m * 1e3);
    println!("    Turn radius:     {:.0} µm", r_m * 1e6);
    println!("    Cycles:          {n_cycles}");
    println!("    Total length:    {:.2} mm", total_length * 1e3);
    println!("  Fluid: Carreau-Yasuda blood, ρ = {:.0} kg/m³", density);

    let u_inlet = 0.05_f64;          // 50 mm/s — moderate millifluidic velocity
    let diffusion_coeff = 1e-9_f64;  // typical molecular diffusivity in blood

    // Analytical mixing estimate
    let mixing_model = AdvectionDiffusionMixing::new(w_m, u_inlet, diffusion_coeff);
    let pe = mixing_model.peclet_number();
    println!("\n  Peclet number:     {:.1e}", pe);

    println!("  Running SerpentineSolver2D...");
    let mut solver = SerpentineSolver2D::new(geom.clone(), blood_model, density, 60, 30);
    let sol = solver.solve(u_inlet, diffusion_coeff, 0.0, 1.0)?;

    println!("\n  📊 2D Simulation Results:");
    println!("    Mixing fraction:  {:.1}% at outlet", sol.mixing_fraction_outlet * 100.0);
    println!("    Pressure drop:    {:.2} Pa", sol.pressure_drop);
    println!("    Well-mixed:       {}", if sol.is_well_mixed() { "YES ✅" } else { "NO ❌" });

    // ── Phase 3: Hemolysis & Medical Analysis ────────────────────────────────
    println!("\n🔬 Phase 3: Hemolysis & Medical Analysis\n");

    let hemolysis_model = HemolysisModel::giersiepen_standard();
    let blood_density = density;
    let vapor_pressure = 6300.0; // Pa — blood at 37°C

    // Estimate per-segment wall shear stress from 2D pressure drop
    // For Poiseuille-like flow in a rectangular channel:
    //   τ_w ≈ (ΔP · D_h) / (4 · L)
    // With the 2D solver giving us the actual resolved pressure drop
    let dp_per_seg = sol.pressure_drop / n_segments as f64;
    let d_h = 2.0 * w_m * h_m / (w_m + h_m); // hydraulic diameter
    let _cross_area = w_m * h_m;

    let mut edge_hemolysis = HashMap::<usize, f64>::new();
    let mut edge_shear = HashMap::<usize, f64>::new();
    let mut edge_cavitation = HashMap::<usize, f64>::new();
    let mut node_pressure = HashMap::<usize, f64>::new();

    let mut max_hi = 0.0_f64;
    let mut max_wss = 0.0_f64;

    // Assign pressure distribution across nodes (linear from 2D solution)
    for i in 0..=n_segments {
        let frac = i as f64 / n_segments as f64;
        let p = sol.pressure_drop * (1.0 - frac);
        node_pressure.insert(i, p);
    }

    // Per-segment analysis
    for seg in 0..n_segments {
        // Segment length accounts for turns (every other segment)
        let seg_length = if seg % 2 == 0 {
            l_straight_m
        } else {
            l_straight_m + std::f64::consts::PI * r_m // straight + U-turn
        };

        // Wall shear stress from resolved pressure gradient
        // τ_w = (ΔP / L) · (D_h / 4) — Poiseuille relation from 2D-resolved ΔP
        let dp_seg = dp_per_seg * seg_length / l_straight_m;
        let wall_shear = dp_seg * d_h / (4.0 * seg_length);

        // Mean velocity (from continuity)
        let velocity = u_inlet; // approximately uniform for incompressible flow

        // Residence time in this segment
        let exposure_time = seg_length / velocity.max(1e-12);

        // Hemolysis index (Giersiepen–Wurzinger power law)
        let hi = hemolysis_model.damage_index(wall_shear, exposure_time).unwrap_or(0.0);

        // Cavitation number
        let local_p = sol.pressure_drop * (1.0 - (seg as f64 + 0.5) / n_segments as f64);
        let cav = CavitationNumber {
            reference_pressure: local_p,
            vapor_pressure,
            density: blood_density,
            velocity,
        };
        let sigma = cav.calculate();

        edge_hemolysis.insert(seg, hi);
        edge_shear.insert(seg, wall_shear);
        edge_cavitation.insert(seg, 1.0 / sigma.max(1e-12));

        max_hi = max_hi.max(hi);
        max_wss = max_wss.max(wall_shear);
    }

    println!("  Max hemolysis index:    {:.4e}", max_hi);
    println!("  Max wall shear stress:  {:.2} Pa", max_wss);
    println!("  Hydraulic diameter:     {:.0} µm", d_h * 1e6);

    // ── Phase 4: Overlay Visualization ───────────────────────────────────────
    println!("\n🎨 Phase 4: Overlay Visualization\n");

    let overlays: Vec<(&str, AnalysisOverlay, String)> = vec![
        (
            "hemolysis_index.png",
            AnalysisOverlay::new(
                AnalysisField::Custom("Hemolysis Index (Giersiepen)".into()),
                ColormapKind::BlueRed,
            )
            .with_edge_data(edge_hemolysis.clone())
            .with_node_data(node_pressure.clone()),
            format!("Hemolysis Index — Giersiepen (max HI = {:.2e})", max_hi),
        ),
        (
            "wall_shear_stress.png",
            AnalysisOverlay::new(AnalysisField::WallShearStress, ColormapKind::BlueRed)
                .with_edge_data(edge_shear.clone())
                .with_node_data(node_pressure.clone()),
            format!("Wall Shear Stress (max τ_w = {:.2} Pa)", max_wss),
        ),
        (
            "cavitation_risk.png",
            AnalysisOverlay::new(
                AnalysisField::Custom("Cavitation Risk (1/σ)".into()),
                ColormapKind::BlueRed,
            )
            .with_edge_data(edge_cavitation.clone())
            .with_node_data(node_pressure.clone()),
            "Cavitation Risk (1/σ — red = high risk)".to_string(),
        ),
        (
            "pressure.png",
            AnalysisOverlay::new(AnalysisField::Pressure, ColormapKind::Viridis)
                .with_edge_data(
                    (0..n_segments)
                        .map(|i| {
                            let frac = (i as f64 + 0.5) / n_segments as f64;
                            (i, sol.pressure_drop * (1.0 - frac))
                        })
                        .collect(),
                )
                .with_node_data(node_pressure.clone()),
            format!("Pressure Distribution (ΔP = {:.2} Pa)", sol.pressure_drop),
        ),
    ];

    for (filename, overlay, title) in &overlays {
        let path = out.join(filename);
        let cfg = RenderConfig {
            width: 1400,
            height: 500,
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

    let channel_data: Vec<serde_json::Value> = (0..n_segments)
        .map(|i| {
            serde_json::json!({
                "segment_id": i,
                "hemolysis_index": edge_hemolysis.get(&i),
                "wall_shear_stress_pa": edge_shear.get(&i),
                "cavitation_risk_inv_sigma": edge_cavitation.get(&i),
            })
        })
        .collect();

    let results = serde_json::json!({
        "analysis": "hemolysis_serpentine_2d",
        "solver": "SerpentineSolver2D (Navier-Stokes FVM)",
        "geometry": {
            "type": "serpentine",
            "width_um": w_mm * 1000.0,
            "height_um": 500.0,
            "straight_mm": l_straight_mm,
            "turn_radius_um": r_mm * 1000.0,
            "cycles": n_cycles,
            "total_length_mm": total_length * 1e3,
            "hydraulic_diameter_um": d_h * 1e6
        },
        "fluid": {
            "model": "Carreau-Yasuda blood",
            "density_kg_m3": blood_density,
        },
        "boundary_conditions": {
            "inlet_velocity_m_s": u_inlet,
            "diffusion_coeff_m2_s": diffusion_coeff
        },
        "simulation_results": {
            "pressure_drop_pa": sol.pressure_drop,
            "mixing_fraction_outlet": sol.mixing_fraction_outlet,
            "peclet_number": sol.peclet,
            "well_mixed": sol.is_well_mixed()
        },
        "hemolysis_analysis": {
            "model": "Giersiepen-Wurzinger",
            "max_hemolysis_index": max_hi,
            "max_wall_shear_stress_pa": max_wss
        },
        "channel_data": channel_data
    });

    let json_path = out.join("results.json");
    fs::write(&json_path, serde_json::to_string_pretty(&results)?)?;
    println!("  ✅ JSON → {}", json_path.display());

    println!("\n✅ Hemolysis serpentine 2D example complete.");
    Ok(())
}
