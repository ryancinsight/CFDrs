//! Venturi Throat Analysis in Parallel Millifluidic Sections
//!
//! Demonstrates frustum (tapered) channels with venturi throats in parallel
//! bifurcation branches. Compares two configurations:
//!
//! 1. **Mild Venturi**: inlet 2.0 mm → throat 1.0 mm → outlet 2.0 mm
//! 2. **Aggressive Venturi**: inlet 2.0 mm → throat 0.5 mm → outlet 2.0 mm
//!
//! At each throat, we evaluate:
//! - Throat velocity (continuity)
//! - Throat pressure (Bernoulli)
//! - Cavitation number σ
//! - Cavity length (Nurick correlation)
//!
//! Pipeline:
//!   cfd-schematics `create_geometry` (AllFrustum) → convert to specs →
//!   cfd-1d NetworkSolver → `VenturiCavitation` analysis →
//!   AnalysisOverlay renders → JSON export
//!
//! Run with:
//! `cargo run -p cfd-1d --example venturi_parallel_analysis`

use cfd_1d::network::{EdgeProperties, Network, NetworkBuilder};
use cfd_1d::resistance::{FlowConditions, ResistanceCalculator, VenturiModel};
use cfd_1d::solver::{NetworkProblem, NetworkSolver, SolverConfig};
use cfd_core::compute::solver::Solver;
use cfd_core::physics::cavitation::VenturiCavitation;
use cfd_core::physics::fluid::ConstantPropertyFluid;
use cfd_schematics::config::{ChannelTypeConfig, FrustumConfig, GeometryConfig, TaperProfile};
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
    println!("🔬 Venturi Throat Analysis — Parallel Millifluidic Sections");
    println!("============================================================\n");

    let out = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("outputs");
    fs::create_dir_all(out.join("venturi_parallel"))?;

    // ── Water at 25 °C ───────────────────────────────────────────────────────
    let water = ConstantPropertyFluid::new(
        "Water".to_string(),
        998.0,   // kg/m³
        8.9e-4,  // Pa·s
        4182.0,  // J/(kg·K)
        0.598,   // W/(m·K)
        1497.0,  // m/s (speed of sound)
    );
    let mu = 8.9e-4_f64;
    let rho = 998.0_f64;
    let p_vapor = 3169.0_f64; // vapor pressure at 25 °C [Pa]

    // ── Configuration A: Mild Venturi (2:1 contraction) ──────────────────────
    let mild_frustum = FrustumConfig {
        inlet_width: 2.0,    // mm
        throat_width: 1.0,   // mm
        outlet_width: 2.0,   // mm
        taper_profile: TaperProfile::Smooth,
        smoothness: 80,
        throat_position: 0.5,
    };

    // ── Configuration B: Aggressive Venturi (4:1 contraction) ────────────────
    let aggressive_frustum = FrustumConfig {
        inlet_width: 2.0,    // mm
        throat_width: 0.5,   // mm
        outlet_width: 2.0,   // mm
        taper_profile: TaperProfile::Smooth,
        smoothness: 80,
        throat_position: 0.5,
    };

    let configs: Vec<(&str, FrustumConfig)> = vec![
        ("mild", mild_frustum),
        ("aggressive", aggressive_frustum),
    ];

    let mut all_results = Vec::new();

    for (label, frustum_config) in &configs {
        println!("━━━ Configuration: {} venturi ━━━", label);

        // ── 1. Generate Frustum Geometry ─────────────────────────────────────
        let box_dims = (90.0, 42.0); // mm — 96-well plate footprint
        let splits = vec![SplitType::Bifurcation];
        let mut geo_config = GeometryConfig::default();
        geo_config.channel_width = frustum_config.inlet_width;
        geo_config.channel_height = 0.5; // 500 µm depth
        let channel_type = ChannelTypeConfig::AllFrustum(*frustum_config);

        let system = create_geometry(box_dims, &splits, &geo_config, &channel_type);
        println!(
            "  Geometry: {} nodes, {} channels",
            system.nodes.len(),
            system.channels.len()
        );

        // Render schematic
        let schematic_path = out.join(format!("venturi_parallel/{label}_schematic.png"));
        plot_geometry(&system, schematic_path.to_str().unwrap())?;
        println!("  Rendered schematic → {label}_schematic.png");

        // ── 2. Convert to Network Specs ──────────────────────────────────────
        let (node_specs, channel_specs) = convert_geometry_to_specs(&system, mu);

        // ── 3. Build Network ─────────────────────────────────────────────────
        let mut builder = NetworkBuilder::<f64>::new();
        let mut id_map = HashMap::new();

        for spec in &node_specs {
            let node = cfd_1d::Node::from(spec);
            let idx = builder.add_node(node);
            id_map.insert(spec.id.as_str().to_string(), idx);
        }

        let mut edge_id_map = HashMap::new();
        for spec in &channel_specs {
            let edge = cfd_1d::Edge::from(spec);
            let from = id_map[spec.from.as_str()];
            let to = id_map[spec.to.as_str()];
            let eidx = builder.add_edge(from, to, edge);
            edge_id_map.insert(spec.id.as_str().to_string(), eidx);
        }

        let graph = builder.build()?;
        let mut network = Network::new(graph, water.clone());

        for spec in &channel_specs {
            let eidx = edge_id_map[spec.id.as_str()];
            network.add_edge_properties(eidx, EdgeProperties::from(spec));
        }

        // ── 4. Boundary Conditions ───────────────────────────────────────────
        let inlet_id = node_specs
            .iter()
            .find(|n| matches!(n.kind, NodeKind::Inlet))
            .unwrap()
            .id
            .as_str()
            .to_string();
        let inlet_idx = id_map[&inlet_id];
        // 2 kPa inlet, 0 Pa outlets — millifluidic-scale driving pressure
        network.set_pressure(inlet_idx, 2_000.0);

        for spec in &node_specs {
            if matches!(spec.kind, NodeKind::Outlet) {
                let idx = id_map[spec.id.as_str()];
                network.set_pressure(idx, 0.0);
            }
        }

        // ── 5. Solve ────────────────────────────────────────────────────────
        let config = SolverConfig { tolerance: 1e-8, max_iterations: 200 };
        let solver = NetworkSolver::<f64, ConstantPropertyFluid<f64>>::with_config(config);
        let solution = solver.solve(&NetworkProblem::new(network))?;
        println!("  Solver converged.");

        // ── 6. Venturi Cavitation Analysis per Channel ───────────────────────
        println!("\n  Channel Results:");
        let mut edge_cavitation_sigma = HashMap::new();
        let mut edge_throat_velocity = HashMap::new();
        let mut edge_throat_pressure = HashMap::new();
        let mut edge_pressure_drop = HashMap::new();
        let mut channel_results = Vec::new();

        let calc = ResistanceCalculator::<f64>::new();

        for spec in &channel_specs {
            let eidx = edge_id_map[spec.id.as_str()];
            let e = solution.graph.edge_weight(eidx).unwrap();
            let q = solution.flow_rates.get(&eidx).copied().unwrap_or(0.0).abs();

            // Upstream/downstream pressures
            let (src, tgt) = solution.graph.edge_endpoints(eidx).unwrap();
            let p_up = solution.pressures.get(&src).copied().unwrap_or(0.0);
            let p_dn = solution.pressures.get(&tgt).copied().unwrap_or(0.0);
            let dp = (p_up - p_dn).abs();

            let props = solution.properties.get(&eidx);
            let length = props.map(|p| p.length).unwrap_or(0.01);

            // Throat geometry from FrustumConfig (SI units)
            let inlet_d = frustum_config.inlet_width / 1000.0; // mm → m
            let throat_d = frustum_config.throat_width / 1000.0;
            let outlet_d = frustum_config.outlet_width / 1000.0;

            // Bulk velocity at channel inlet
            let inlet_area = std::f64::consts::PI * (inlet_d / 2.0).powi(2);
            let v_inlet = if inlet_area > 0.0 { q / inlet_area } else { 0.0 };

            // VenturiCavitation analysis
            let venturi = VenturiCavitation {
                inlet_diameter: inlet_d,
                throat_diameter: throat_d,
                outlet_diameter: outlet_d,
                convergent_angle: 0.2618, // ~15°
                divergent_angle: 0.1222,  // ~7°
                inlet_pressure: p_up,
                inlet_velocity: v_inlet,
                density: rho,
                vapor_pressure: p_vapor,
            };

            let sigma = venturi.cavitation_number();
            let v_throat = venturi.throat_velocity();
            let p_throat = venturi.throat_pressure();
            let is_cav = venturi.is_cavitating();
            let cavity_len = venturi.cavity_length(sigma);

            // Venturi-specific resistance via 1D model
            let venturi_model = VenturiModel::millifluidic(inlet_d, throat_d, length);
            let conditions = FlowConditions {
                reynolds_number: None,
                velocity: None,
                flow_rate: Some(q),
                shear_rate: None,
                temperature: 298.15,
                pressure: p_up,
            };
            let venturi_r = calc
                .calculate_venturi_coefficients(&venturi_model, &water, &conditions)
                .map(|(r, _k)| r)
                .unwrap_or(0.0);

            println!("    {} | Q={:.3e} m³/s | ΔP={:.0} Pa", e.id, q, dp);
            println!(
                "      Throat: v={:.2} m/s  p={:.0} Pa  σ={:.3}  cavitating={}",
                v_throat, p_throat, sigma, is_cav
            );
            if cavity_len > 0.0 {
                println!("      Cavity length: {:.3} mm", cavity_len * 1000.0);
            }
            println!("      Venturi R = {:.3e} Pa·s/m³", venturi_r);

            edge_cavitation_sigma.insert(eidx.index(), sigma);
            edge_throat_velocity.insert(eidx.index(), v_throat);
            edge_throat_pressure.insert(eidx.index(), p_throat);
            edge_pressure_drop.insert(eidx.index(), dp);

            channel_results.push(serde_json::json!({
                "id": e.id,
                "flow_rate_m3s": q,
                "pressure_drop_pa": dp,
                "resistance_solver": e.resistance,
                "resistance_venturi": venturi_r,
                "throat_velocity_ms": v_throat,
                "throat_pressure_pa": p_throat,
                "cavitation_number": sigma,
                "is_cavitating": is_cav,
                "cavity_length_m": cavity_len,
            }));
        }

        // ── 7. Visualization Overlays ────────────────────────────────────────
        let renderer = create_plotters_renderer();

        let overlays: Vec<(&str, AnalysisField, &HashMap<usize, f64>, ColormapKind)> = vec![
            ("cavitation_sigma", AnalysisField::Custom("σ".to_string()), &edge_cavitation_sigma, ColormapKind::BlueRed),
            ("throat_velocity", AnalysisField::Velocity, &edge_throat_velocity, ColormapKind::Viridis),
            ("throat_pressure", AnalysisField::Pressure, &edge_throat_pressure, ColormapKind::BlueRed),
            ("pressure_drop", AnalysisField::Custom("ΔP".to_string()), &edge_pressure_drop, ColormapKind::Viridis),
        ];

        for (filename, field, data, cmap) in &overlays {
            let overlay = AnalysisOverlay::new(field.clone(), *cmap)
                .with_edge_data((*data).clone());
            let mut rc = RenderConfig::default();
            rc.title = format!("{label} venturi — {filename}");
            rc.show_axes = true;
            rc.show_grid = false;
            let path = out.join(format!("venturi_parallel/{label}_{filename}.png"));
            renderer.render_analysis(&system, path.to_str().unwrap(), &rc, &overlay)?;
        }
        println!("  Rendered 4 overlay maps.");

        // ── 8. Collect JSON ──────────────────────────────────────────────────
        let node_results: Vec<_> = solution.graph.node_indices().map(|idx| {
            let n = solution.graph.node_weight(idx).unwrap();
            serde_json::json!({
                "id": n.id,
                "pressure_pa": solution.pressures.get(&idx).unwrap_or(&0.0),
                "type": format!("{:?}", n.node_type)
            })
        }).collect();

        all_results.push(serde_json::json!({
            "configuration": label,
            "frustum": {
                "inlet_width_mm": frustum_config.inlet_width,
                "throat_width_mm": frustum_config.throat_width,
                "outlet_width_mm": frustum_config.outlet_width,
                "contraction_ratio": frustum_config.inlet_width / frustum_config.throat_width,
            },
            "nodes": node_results,
            "channels": channel_results,
        }));

        println!();
    }

    // ── 9. Export Combined JSON ──────────────────────────────────────────────
    let report = serde_json::json!({
        "analysis": "venturi_parallel_comparison",
        "fluid": "Water 25°C",
        "density_kg_m3": rho,
        "viscosity_pa_s": mu,
        "vapor_pressure_pa": p_vapor,
        "configurations": all_results,
    });
    let json_path = out.join("venturi_parallel/results.json");
    fs::write(&json_path, serde_json::to_string_pretty(&report)?)?;
    println!("✅ Exported combined results to {}", json_path.display());

    Ok(())
}

/// Convert a `ChannelSystem` geometry into `NodeSpec`/`ChannelSpec` pairs.
///
/// Nodes at the minimum x-coordinate are classified as inlets; those at the
/// maximum x-coordinate as outlets; all others as junctions.
/// Channels use `ChannelSpec::new_pipe_rect` with the Hele-Shaw rectangular
/// resistance formula: `R = 12·μ·L / (w·h³·(1 − 0.63·h/w))`.
fn convert_geometry_to_specs(system: &ChannelSystem, mu: f64) -> (Vec<NodeSpec>, Vec<ChannelSpec>) {
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

    let channel_specs: Vec<ChannelSpec> = system
        .channels
        .iter()
        .map(|channel| {
            let from_node = system.nodes.iter().find(|n| n.id == channel.from_node).unwrap();
            let to_node = system.nodes.iter().find(|n| n.id == channel.to_node).unwrap();
            let dx = from_node.point.0 - to_node.point.0;
            let dy = from_node.point.1 - to_node.point.1;
            let length_m = dx.hypot(dy) / 1000.0; // mm → m

            let width_m = channel.width / 1000.0;
            let height_m = channel.height / 1000.0;

            let (w, h) = if width_m > height_m { (width_m, height_m) } else { (height_m, width_m) };
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
