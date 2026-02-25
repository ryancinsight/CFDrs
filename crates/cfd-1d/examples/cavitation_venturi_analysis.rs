//! Cavitation Risk Analysis in a Venturi-Style Millifluidic Device
//!
//! Demonstrates cavitation assessment for Sonodynamic Therapy (SDT) design:
//! 1. Generate a bifurcation geometry with frustum (tapered) channels via `cfd-schematics`
//! 2. Convert to 1D network and solve water flow
//! 3. Compute cavitation number σ = (p − p_v)/(0.5·ρ·v²) per channel
//! 4. Classify cavitation regime (no cavitation / inception / developed / supercavitation)
//! 5. Render colored schematics for cavitation number and pressure distribution
//! 6. Export comprehensive JSON results
//!
//! Cavitation is critical in SDT devices where controlled bubble formation
//! enhances sonosensitizer activation within the 96-well plate footprint.
//!
//! Run with:
//! `cargo run -p cfd-1d --example cavitation_venturi_analysis`

use cfd_1d::network::{EdgeProperties, Network, NetworkBuilder};
use cfd_1d::solver::{NetworkProblem, NetworkSolver, SolverConfig};
use cfd_core::compute::solver::Solver;
use cfd_core::physics::cavitation::CavitationNumber;
use cfd_core::physics::fluid::ConstantPropertyFluid;
use cfd_schematics::config::{ChannelTypeConfig, GeometryConfig};
use cfd_schematics::domain::model::{ChannelSpec, NodeKind, NodeSpec};
use cfd_schematics::geometry::generator::create_geometry;
use cfd_schematics::geometry::{ChannelSystem, SplitType};
use cfd_schematics::plot_geometry;
use cfd_schematics::visualizations::analysis_field::{AnalysisField, AnalysisOverlay, ColormapKind};
use cfd_schematics::visualizations::plotters_backend::create_plotters_renderer;
use cfd_schematics::visualizations::traits::SchematicRenderer;
use cfd_schematics::visualizations::RenderConfig;
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("💧 Cavitation Risk Analysis — Venturi Millifluidic Device");
    println!("=========================================================\n");

    let out = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("outputs");
    fs::create_dir_all(out.join("cavitation"))?;

    // ── 1. Generate Frustum Bifurcation Geometry ─────────────────────────────
    println!("1. Generating bifurcation geometry with frustum channels...");
    let box_dims = (100.0, 50.0); // mm
    let splits = vec![SplitType::Bifurcation, SplitType::Bifurcation];
    let mut geo_config = GeometryConfig::default();
    geo_config.channel_width = 1.5; // 1.5 mm — wider inlet for venturi acceleration
    let channel_type = ChannelTypeConfig::AllStraight;

    let system = create_geometry(box_dims, &splits, &geo_config, &channel_type);
    println!(
        "   {} nodes, {} channels",
        system.nodes.len(),
        system.channels.len()
    );

    // ── 2. Plain schematic ───────────────────────────────────────────────────
    println!("2. Rendering plain schematic...");
    plot_geometry(
        &system,
        out.join("cavitation/channel_schematic.png")
            .to_str()
            .unwrap(),
    )?;

    // ── 3. Convert geometry → 1D specs ───────────────────────────────────────
    println!("3. Converting to simulation specifications...");
    let (node_specs, channel_specs) = convert_geometry_to_specs(&system);

    // ── 4. Build network with water ──────────────────────────────────────────
    println!("4. Building network (water at 25°C)...");
    let water = ConstantPropertyFluid::<f64>::new(
        "Water (25°C)".to_string(),
        997.0,   // kg/m³
        8.9e-4,  // Pa·s
        4186.0,  // J/(kg·K)
        0.606,   // W/(m·K)
        1497.0,  // m/s (speed of sound)
    );

    let mut builder = NetworkBuilder::<f64>::new();
    let mut id_map = HashMap::new();
    let mut edge_id_map = HashMap::new();

    for spec in &node_specs {
        let node = cfd_1d::Node::from(spec);
        let idx = builder.add_node(node);
        id_map.insert(spec.id.as_str().to_string(), idx);
    }
    for spec in &channel_specs {
        let edge = cfd_1d::Edge::from(spec);
        let from_idx = id_map[spec.from.as_str()];
        let to_idx = id_map[spec.to.as_str()];
        let eidx = builder.add_edge(from_idx, to_idx, edge);
        edge_id_map.insert(spec.id.as_str().to_string(), eidx);
    }

    let graph = builder.build()?;
    let mut network = Network::new(graph, water);

    for spec in &channel_specs {
        let eidx = edge_id_map[spec.id.as_str()];
        network.add_edge_properties(eidx, EdgeProperties::from(spec));
    }

    // ── 5. Boundary Conditions (high pressure for cavitation) ────────────────
    println!("5. Applying boundary conditions (high ΔP for cavitation study)...");
    let inlet_pressure = 50_000.0; // 50 kPa — aggressive for millifluidic cavitation
    let outlet_pressure = 0.0;

    for spec in &node_specs {
        let idx = id_map[spec.id.as_str()];
        match spec.kind {
            NodeKind::Inlet => network.set_pressure(idx, inlet_pressure),
            NodeKind::Outlet => network.set_pressure(idx, outlet_pressure),
            _ => {}
        }
    }

    // ── 6. Solve ─────────────────────────────────────────────────────────────
    println!("6. Solving 1D network...");
    let solver = NetworkSolver::with_config(SolverConfig {
        tolerance: 1e-8,
        max_iterations: 200,
    });
    let solution = solver.solve(&NetworkProblem::new(network))?;
    println!("   Solution converged.");

    // ── 7. Compute cavitation metrics ────────────────────────────────────────
    println!("7. Computing cavitation numbers and regime classification...");
    let vapor_pressure_25c = 3169.0; // Pa — water at 25°C

    let mut edge_sigma = HashMap::<usize, f64>::new();
    let mut edge_velocity = HashMap::<usize, f64>::new();
    let mut edge_pressure = HashMap::<usize, f64>::new();
    let mut node_pressure_data = HashMap::<usize, f64>::new();
    let mut regime_counts: HashMap<String, usize> = HashMap::new();

    // Node pressures
    for (nidx, &p) in &solution.pressures {
        if let Some(node) = solution.graph.node_weight(*nidx) {
            if let Ok(id) = node.id.trim_start_matches("node_").parse::<usize>() {
                node_pressure_data.insert(id, p);
            }
        }
    }

    for (eidx, &q) in &solution.flow_rates {
        let Some(edge) = solution.graph.edge_weight(*eidx) else {
            continue;
        };
        let chan_id = match edge.id.trim_start_matches("chan_").parse::<usize>() {
            Ok(id) => id,
            Err(_) => continue,
        };

        if let Some(p) = solution.properties.get(eidx) {
            let velocity = q.abs() / p.area;

            // Local pressure estimate: average of upstream/downstream node pressures
            let (src, dst) = solution.graph.edge_endpoints(*eidx).unwrap();
            let p_src = *solution.pressures.get(&src).unwrap_or(&inlet_pressure);
            let p_dst = *solution.pressures.get(&dst).unwrap_or(&outlet_pressure);
            let local_pressure = (p_src + p_dst) / 2.0;

            let cav_num = CavitationNumber {
                reference_pressure: local_pressure,
                vapor_pressure: vapor_pressure_25c,
                density: 997.0,
                velocity,
            };
            let sigma = cav_num.calculate();

            // Classify regime by sigma thresholds (from literature):
            // σ > 1.5: no cavitation, 0.5 < σ < 1.5: inception, σ < 0.5: developed
            let regime = if sigma > 1.5 {
                "NoCavitation"
            } else if sigma > 0.5 {
                "Inception"
            } else {
                "Developed"
            };

            edge_sigma.insert(chan_id, sigma);
            edge_velocity.insert(chan_id, velocity);
            edge_pressure.insert(chan_id, local_pressure);

            *regime_counts.entry(regime.to_string()).or_insert(0) += 1;

            if sigma < 2.0 {
                println!(
                    "   ⚠ Channel {}: σ = {:.3}, v = {:.2} m/s, regime = {}",
                    chan_id, sigma, velocity, regime
                );
            }
        }
    }

    println!("\n   Regime distribution:");
    for (regime, count) in &regime_counts {
        println!("     {regime}: {count} channels");
    }

    // ── 8. Render cavitation number overlay ──────────────────────────────────
    println!("\n8. Rendering cavitation number overlay...");
    let renderer = create_plotters_renderer();

    // Invert sigma for coloring: low sigma (high risk) → red
    let max_sigma = edge_sigma.values().copied().fold(0.0_f64, f64::max);
    let inverted_sigma: HashMap<usize, f64> = edge_sigma
        .iter()
        .map(|(&k, &v)| (k, max_sigma - v))
        .collect();

    let cav_overlay =
        AnalysisOverlay::new(AnalysisField::Custom("Cavitation Number σ".into()), ColormapKind::BlueRed)
            .with_edge_data(inverted_sigma)
            .with_node_data(node_pressure_data.clone());

    let mut config_cav = RenderConfig::default();
    config_cav.title = format!(
        "Cavitation Number σ (red = low σ = high risk, max σ = {:.1})",
        max_sigma
    );
    config_cav.show_axes = true;
    config_cav.show_grid = false;

    renderer.render_analysis(
        &system,
        out.join("cavitation/cavitation_number.png")
            .to_str()
            .unwrap(),
        &config_cav,
        &cav_overlay,
    )?;

    // ── 9. Render pressure distribution overlay ──────────────────────────────
    println!("9. Rendering pressure distribution overlay...");
    let pressure_overlay = AnalysisOverlay::new(AnalysisField::Pressure, ColormapKind::Viridis)
        .with_edge_data(edge_pressure.clone())
        .with_node_data(node_pressure_data.clone());

    let mut config_p = RenderConfig::default();
    config_p.title = "Pressure Distribution [Pa]".to_string();
    config_p.show_axes = true;
    config_p.show_grid = false;

    renderer.render_analysis(
        &system,
        out.join("cavitation/pressure_distribution.png")
            .to_str()
            .unwrap(),
        &config_p,
        &pressure_overlay,
    )?;

    // ── 10. Render velocity overlay ──────────────────────────────────────────
    println!("10. Rendering velocity overlay...");
    let vel_overlay = AnalysisOverlay::new(AnalysisField::Velocity, ColormapKind::Viridis)
        .with_edge_data(edge_velocity.clone())
        .with_node_data(node_pressure_data.clone());

    let mut config_v = RenderConfig::default();
    config_v.title = "Mean Velocity [m/s]".to_string();
    config_v.show_axes = true;
    config_v.show_grid = false;

    renderer.render_analysis(
        &system,
        out.join("cavitation/velocity_distribution.png")
            .to_str()
            .unwrap(),
        &config_v,
        &vel_overlay,
    )?;

    // ── 11. JSON Export ──────────────────────────────────────────────────────
    println!("11. Exporting results...");
    let results = serde_json::json!({
        "analysis": "cavitation_venturi",
        "geometry": {
            "type": "bifurcation_straight",
            "box_mm": [box_dims.0, box_dims.1],
            "channel_width_mm": 1.5,
            "num_nodes": system.nodes.len(),
            "num_channels": system.channels.len()
        },
        "fluid": {
            "name": "Water (25°C)",
            "density_kg_m3": 997.0,
            "viscosity_pa_s": 8.9e-4,
            "vapor_pressure_pa": vapor_pressure_25c
        },
        "boundary_conditions": {
            "inlet_pressure_pa": inlet_pressure,
            "outlet_pressure_pa": outlet_pressure
        },
        "cavitation_analysis": {
            "regime_distribution": regime_counts,
            "channel_data": edge_sigma.iter().map(|(k, v)| {
                serde_json::json!({
                    "channel_id": k,
                    "cavitation_number": v,
                    "velocity_m_s": edge_velocity.get(k),
                    "local_pressure_pa": edge_pressure.get(k)
                })
            }).collect::<Vec<_>>()
        },
        "outputs": {
            "schematic": "cavitation/channel_schematic.png",
            "cavitation_plot": "cavitation/cavitation_number.png",
            "pressure_plot": "cavitation/pressure_distribution.png",
            "velocity_plot": "cavitation/velocity_distribution.png"
        }
    });

    let json_path = out.join("cavitation/results.json");
    fs::write(&json_path, serde_json::to_string_pretty(&results)?)?;
    println!("   Results exported to {}", json_path.display());

    println!("\n✅ Cavitation analysis complete — all outputs in cavitation/");
    Ok(())
}

/// Convert `ChannelSystem` → `NodeSpec`/`ChannelSpec` with Hagen-Poiseuille rectangular resistance.
fn convert_geometry_to_specs(system: &ChannelSystem) -> (Vec<NodeSpec>, Vec<ChannelSpec>) {
    let min_x = system
        .nodes
        .iter()
        .map(|n| n.point.0)
        .fold(f64::INFINITY, f64::min);
    let max_x = system
        .nodes
        .iter()
        .map(|n| n.point.0)
        .fold(f64::NEG_INFINITY, f64::max);

    let node_specs: Vec<NodeSpec> = system
        .nodes
        .iter()
        .map(|node| {
            let kind = if (node.point.0 - min_x).abs() < 1e-3 {
                NodeKind::Inlet
            } else if (node.point.0 - max_x).abs() < 1e-3 {
                NodeKind::Outlet
            } else {
                NodeKind::Junction
            };
            NodeSpec::new(format!("node_{}", node.id), kind)
        })
        .collect();

    let mu = 8.9e-4_f64; // Water viscosity at 25°C

    let channel_specs: Vec<ChannelSpec> = system
        .channels
        .iter()
        .map(|channel| {
            let from = system.nodes.iter().find(|n| n.id == channel.from_node).unwrap();
            let to = system.nodes.iter().find(|n| n.id == channel.to_node).unwrap();
            let dx = from.point.0 - to.point.0;
            let dy = from.point.1 - to.point.1;
            let length_m = dx.hypot(dy) / 1000.0;

            let width_m = channel.width / 1000.0;
            let height_m = channel.height / 1000.0;

            let (w, h) = if width_m > height_m {
                (width_m, height_m)
            } else {
                (height_m, width_m)
            };
            let resistance = if h > 0.0 {
                (12.0 * mu * length_m) / (w * h.powi(3) * (1.0 - 0.63 * h / w))
            } else {
                f64::INFINITY
            };

            ChannelSpec::new_pipe_rect(
                format!("chan_{}", channel.id),
                format!("node_{}", channel.from_node),
                format!("node_{}", channel.to_node),
                length_m,
                width_m,
                height_m,
                resistance,
                0.0,
            )
        })
        .collect();

    (node_specs, channel_specs)
}
