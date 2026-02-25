//! Hemolysis Analysis in a Serpentine Millifluidic Channel
//!
//! Demonstrates a full-pipeline medical-grade analysis:
//! 1. Generate serpentine geometry via `cfd-schematics`
//! 2. Convert to 1D network specs and solve blood flow (Carreau-Yasuda)
//! 3. Compute wall shear stress and hemolysis index per channel
//! 4. Render colored schematics for hemolysis index and wall shear stress
//! 5. Export comprehensive JSON results
//!
//! The Giersiepen–Wurzinger power-law model computes blood damage:
//!   D = C · τ^α · t^β
//! where τ = wall shear stress, t = exposure time (residence time in channel).
//!
//! Run with:
//! `cargo run -p cfd-1d --example hemolysis_serpentine_analysis`

use cfd_1d::network::{EdgeProperties, Network, NetworkBuilder};
use cfd_1d::solver::{NetworkProblem, NetworkSolver, SolverConfig};
use cfd_core::compute::solver::Solver;
use cfd_core::physics::cavitation::CavitationNumber;
use cfd_core::physics::fluid::non_newtonian::CarreauYasuda;
use cfd_core::physics::fluid::FluidTrait;
use cfd_core::physics::hemolysis::HemolysisModel;
use cfd_schematics::config::presets::smooth_serpentine;
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
    println!("🩸 Hemolysis Analysis — Serpentine Millifluidic Channel");
    println!("=======================================================\n");

    let out = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("outputs");
    fs::create_dir_all(out.join("hemolysis"))?;

    // ── 1. Generate Serpentine Geometry ───────────────────────────────────────
    println!("1. Generating serpentine bifurcation geometry...");
    let box_dims = (80.0, 40.0); // mm — fits within 96-well plate footprint
    let splits = vec![SplitType::Bifurcation, SplitType::Bifurcation];
    let mut geo_config = GeometryConfig::default();
    geo_config.channel_width = 0.8; // 800 µm — millifluidic scale
    let serpentine_config = smooth_serpentine();
    let channel_type = ChannelTypeConfig::AllSerpentine(serpentine_config);

    let system = create_geometry(box_dims, &splits, &geo_config, &channel_type);
    println!(
        "   {} nodes, {} channels generated",
        system.nodes.len(),
        system.channels.len()
    );

    // ── 2. Render plain schematic ────────────────────────────────────────────
    println!("2. Rendering plain channel schematic...");
    plot_geometry(
        &system,
        out.join("hemolysis/channel_schematic.png")
            .to_str()
            .unwrap(),
    )?;

    // ── 3. Convert to 1D simulation specs ────────────────────────────────────
    println!("3. Converting to 1D network specifications...");
    let (node_specs, channel_specs) = convert_geometry_to_specs(&system);

    // ── 4. Build Network with Carreau-Yasuda blood ───────────────────────────
    println!("4. Building 1D network with Carreau-Yasuda blood model...");
    let blood = CarreauYasuda::<f64>::blood();
    println!("   Fluid: {}", blood.name());
    println!(
        "   µ∞ = {:.4} mPa·s, µ₀ = {:.1} mPa·s",
        blood.viscosity_inf * 1e3,
        blood.viscosity_zero * 1e3
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
    let mut network = Network::new(graph, blood.clone());

    for spec in &channel_specs {
        let eidx = edge_id_map[spec.id.as_str()];
        network.add_edge_properties(eidx, EdgeProperties::from(spec));
    }

    // ── 5. Boundary Conditions ───────────────────────────────────────────────
    println!("5. Applying boundary conditions...");
    let inlet_pressure_pa = 2000.0; // 2 kPa — moderate millifluidic pressure
    let outlet_pressure_pa = 0.0;

    for spec in &node_specs {
        let idx = id_map[spec.id.as_str()];
        match spec.kind {
            NodeKind::Inlet => network.set_pressure(idx, inlet_pressure_pa),
            NodeKind::Outlet => network.set_pressure(idx, outlet_pressure_pa),
            _ => {}
        }
    }

    // ── 6. Solve ─────────────────────────────────────────────────────────────
    println!("6. Solving 1D network...");
    let config = SolverConfig {
        tolerance: 1e-8,
        max_iterations: 200,
    };
    let solver = NetworkSolver::with_config(config);
    let solution = solver.solve(&NetworkProblem::new(network))?;
    println!("   Solution converged.");

    // ── 7. Compute hemolysis index & shear stress per channel ────────────────
    println!("7. Computing hemolysis indices (Giersiepen–Wurzinger model)...");
    let hemolysis_model = HemolysisModel::giersiepen_standard();

    let mut edge_hemolysis = HashMap::<usize, f64>::new();
    let mut edge_shear = HashMap::<usize, f64>::new();
    let mut edge_cavitation = HashMap::<usize, f64>::new();
    let mut node_pressure = HashMap::<usize, f64>::new();

    // Extract node pressures
    for (nidx, &p) in &solution.pressures {
        if let Some(node) = solution.graph.node_weight(*nidx) {
            if let Ok(id) = node.id.trim_start_matches("node_").parse::<usize>() {
                node_pressure.insert(id, p);
            }
        }
    }

    let mut max_hemolysis = 0.0_f64;
    let mut max_shear = 0.0_f64;

    for (eidx, &q) in &solution.flow_rates {
        let Some(edge) = solution.graph.edge_weight(*eidx) else {
            continue;
        };
        let props = solution.properties.get(eidx);

        // Parse channel ID
        let chan_id = match edge.id.trim_start_matches("chan_").parse::<usize>() {
            Ok(id) => id,
            Err(_) => continue,
        };

        if let Some(p) = props {
            let velocity = q.abs() / p.area;
            let shear_rate = p
                .hydraulic_diameter
                .map(|d| 8.0 * velocity / d)
                .unwrap_or(0.0);

            // Apparent viscosity at this shear rate
            let viscosity = blood
                .viscosity_at_shear(shear_rate, 310.15, 101_325.0)
                .unwrap_or(3.5e-3);
            let wall_shear = viscosity * shear_rate;

            // Exposure time = channel length / mean velocity
            let exposure_time = if velocity > 1e-12 {
                p.length / velocity
            } else {
                0.0
            };

            // Hemolysis index (Giersiepen power law)
            let hi = hemolysis_model
                .damage_index(wall_shear, exposure_time)
                .unwrap_or(0.0);

            // Cavitation number: σ = (p - p_v) / (0.5 · ρ · v²)
            // For blood at 37°C, vapor pressure ≈ 6.3 kPa
            let cav_num = CavitationNumber {
                reference_pressure: inlet_pressure_pa,
                vapor_pressure: 6300.0,
                density: blood.density,
                velocity,
            };
            let sigma = cav_num.calculate();

            edge_hemolysis.insert(chan_id, hi);
            edge_shear.insert(chan_id, wall_shear);
            // Store inverse cavitation number → higher = more risk
            edge_cavitation.insert(chan_id, 1.0 / sigma);

            max_hemolysis = max_hemolysis.max(hi);
            max_shear = max_shear.max(wall_shear);
        }
    }

    println!(
        "   Max hemolysis index: {:.4e}",
        max_hemolysis
    );
    println!("   Max wall shear stress: {:.2} Pa", max_shear);

    // ── 8. Render hemolysis overlay ──────────────────────────────────────────
    println!("8. Rendering hemolysis index overlay...");
    let renderer = create_plotters_renderer();

    let hemolysis_overlay = AnalysisOverlay::new(
        AnalysisField::Custom("Hemolysis Index (Giersiepen)".into()),
        ColormapKind::BlueRed,
    )
    .with_edge_data(edge_hemolysis.clone())
    .with_node_data(node_pressure.clone());

    let mut config_hi = RenderConfig::default();
    config_hi.title = format!(
        "Hemolysis Index — Giersiepen–Wurzinger (max HI = {:.2e})",
        max_hemolysis
    );
    config_hi.show_axes = true;
    config_hi.show_grid = false;

    renderer.render_analysis(
        &system,
        out.join("hemolysis/hemolysis_index.png")
            .to_str()
            .unwrap(),
        &config_hi,
        &hemolysis_overlay,
    )?;

    // ── 9. Render wall shear stress overlay ──────────────────────────────────
    println!("9. Rendering wall shear stress overlay...");
    let shear_overlay = AnalysisOverlay::new(
        AnalysisField::WallShearStress,
        ColormapKind::BlueRed,
    )
    .with_edge_data(edge_shear.clone())
    .with_node_data(node_pressure.clone());

    let mut config_wss = RenderConfig::default();
    config_wss.title = format!(
        "Wall Shear Stress (max τ_w = {:.2} Pa)",
        max_shear
    );
    config_wss.show_axes = true;
    config_wss.show_grid = false;

    renderer.render_analysis(
        &system,
        out.join("hemolysis/wall_shear_stress.png")
            .to_str()
            .unwrap(),
        &config_wss,
        &shear_overlay,
    )?;

    // ── 10. Render cavitation risk overlay ───────────────────────────────────
    println!("10. Rendering cavitation risk overlay...");
    let cavitation_overlay = AnalysisOverlay::new(
        AnalysisField::Custom("Cavitation Risk (1/σ)".into()),
        ColormapKind::BlueRed,
    )
    .with_edge_data(edge_cavitation.clone())
    .with_node_data(node_pressure.clone());

    let mut config_cav = RenderConfig::default();
    config_cav.title = "Cavitation Risk Index (1/σ)".to_string();
    config_cav.show_axes = true;
    config_cav.show_grid = false;

    renderer.render_analysis(
        &system,
        out.join("hemolysis/cavitation_risk.png")
            .to_str()
            .unwrap(),
        &config_cav,
        &cavitation_overlay,
    )?;

    // ── 11. JSON Export ──────────────────────────────────────────────────────
    println!("11. Exporting results to JSON...");
    let results = serde_json::json!({
        "analysis": "hemolysis_serpentine",
        "geometry": {
            "type": "serpentine_bifurcation",
            "box_mm": [box_dims.0, box_dims.1],
            "channel_width_mm": 0.8,
            "num_nodes": system.nodes.len(),
            "num_channels": system.channels.len()
        },
        "fluid": {
            "model": "Carreau-Yasuda Blood",
            "density_kg_m3": blood.density,
            "viscosity_inf_pa_s": blood.viscosity_inf,
            "viscosity_zero_pa_s": blood.viscosity_zero
        },
        "boundary_conditions": {
            "inlet_pressure_pa": inlet_pressure_pa,
            "outlet_pressure_pa": outlet_pressure_pa
        },
        "hemolysis_model": {
            "type": "Giersiepen-Wurzinger",
            "C": 3.62e-5,
            "alpha": 2.416,
            "beta": 0.785
        },
        "results": {
            "max_hemolysis_index": max_hemolysis,
            "max_wall_shear_stress_pa": max_shear,
            "channel_hemolysis": edge_hemolysis.iter().map(|(k, v)| {
                serde_json::json!({"channel_id": k, "hemolysis_index": v})
            }).collect::<Vec<_>>(),
            "channel_wall_shear": edge_shear.iter().map(|(k, v)| {
                serde_json::json!({"channel_id": k, "wall_shear_stress_pa": v})
            }).collect::<Vec<_>>()
        },
        "outputs": {
            "schematic": "hemolysis/channel_schematic.png",
            "hemolysis_plot": "hemolysis/hemolysis_index.png",
            "shear_plot": "hemolysis/wall_shear_stress.png",
            "cavitation_plot": "hemolysis/cavitation_risk.png"
        }
    });

    let json_path = out.join("hemolysis/results.json");
    fs::write(&json_path, serde_json::to_string_pretty(&results)?)?;
    println!("   Results exported to {}", json_path.display());

    println!("\n✅ Hemolysis analysis complete — all outputs in hemolysis/");
    Ok(())
}

/// Convert a `ChannelSystem` geometry into `NodeSpec`/`ChannelSpec` pairs.
///
/// Nodes at minimum x are inlets; maximum x are outlets; all others are junctions.
/// Channels use the Hele-Shaw rectangular resistance model:
///   R = 12·μ·L / (w·h³·(1 − 0.63·h/w))
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

    // Blood apparent viscosity at mid-shear ≈ 3.5 mPa·s
    let mu = 3.5e-3_f64;

    let channel_specs: Vec<ChannelSpec> = system
        .channels
        .iter()
        .map(|channel| {
            let from = system.nodes.iter().find(|n| n.id == channel.from_node).unwrap();
            let to = system.nodes.iter().find(|n| n.id == channel.to_node).unwrap();
            let dx = from.point.0 - to.point.0;
            let dy = from.point.1 - to.point.1;
            let length_m = dx.hypot(dy) / 1000.0; // mm → m

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
