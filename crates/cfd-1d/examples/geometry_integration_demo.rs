//! Geometry Integration Demo
//!
//! Demonstrates the full pipeline:
//! 1. Generate a bifurcation geometry via `cfd-schematics`
//! 2. Convert geometry to `ChannelSpec`/`NodeSpec` (no `cfd_1d::domain::channel` imports)
//! 3. Build and solve the 1D network
//! 4. Visualize flow rate distribution using `AnalysisOverlay` (typed CFD field)
//! 5. Export results to JSON
//!
//! Run with:
//! `cargo run -p cfd-1d --example geometry_integration_demo`

use cfd_1d::domain::network::{EdgeProperties, Network, NetworkBuilder};
use cfd_1d::solver::core::{NetworkProblem, NetworkSolver, SolverConfig};
use cfd_core::compute::solver::Solver;
use cfd_schematics::config::{ChannelTypeConfig, GeometryConfig};
use cfd_schematics::domain::model::NodeKind;
use cfd_schematics::geometry::generator::create_geometry;
use cfd_schematics::geometry::SplitType;
use cfd_schematics::visualizations::analysis_field::{
    AnalysisField, AnalysisOverlay, ColormapKind,
};
use cfd_schematics::visualizations::plotters_backend::create_plotters_renderer;
use cfd_schematics::visualizations::traits::SchematicRenderer;
use cfd_schematics::visualizations::RenderConfig;
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🧪 Geometry Integration Demo");
    println!("===========================");

    // ── 1. Generate Geometry ─────────────────────────────────────────────────
    println!("DATA: Generating bifurcation geometry...");
    let box_dims = (100.0, 50.0); // mm
    let splits = vec![SplitType::Bifurcation, SplitType::Bifurcation];
    let mut geo_config = GeometryConfig::default();
    geo_config.channel_width = 1.0; // mm
    let channel_type_config = ChannelTypeConfig::AllStraight;

    let system = create_geometry(box_dims, &splits, &geo_config, &channel_type_config);
    println!(
        "DATA: Generated {} nodes and {} channels",
        system.nodes.len(),
        system.channels.len()
    );

    // ── 2. Convert to Simulation Specs ───────────────────────────────────────
    let node_specs = system.nodes.clone();
    let channel_specs = system.channels.clone();

    // ── 3. Build Network ─────────────────────────────────────────────────────
    println!("DATA: Building network...");
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
        let from_idx = id_map[spec.from.as_str()];
        let to_idx = id_map[spec.to.as_str()];
        let eidx = builder.add_edge(from_idx, to_idx, edge);
        edge_id_map.insert(spec.id.as_str().to_string(), eidx);
    }

    let graph = builder.build()?;

    let fluid = cfd_core::physics::fluid::ConstantPropertyFluid::new(
        "Water".to_string(),
        1000.0, // kg/m³
        1e-3,   // Pa·s
        4186.0, // J/(kg·K)
        0.6,    // W/(m·K)
        1480.0, // m/s
    );
    let mut network = Network::new(graph, fluid);

    // Attach EdgeProperties from ChannelSpec (no cfd_1d::domain::channel imports needed)
    for spec in &channel_specs {
        let eidx = edge_id_map[spec.id.as_str()];
        network.add_edge_properties(eidx, EdgeProperties::from(spec));
    }

    // ── 4. Boundary Conditions ───────────────────────────────────────────────
    println!("DATA: Applying boundary conditions...");
    let inlet_id = node_specs
        .iter()
        .find(|n| matches!(n.kind, NodeKind::Inlet))
        .unwrap()
        .id
        .as_str()
        .to_string();
    let inlet_idx = id_map[&inlet_id];
    network.set_pressure(inlet_idx, 1000.0); // 1 kPa

    for spec in &node_specs {
        if matches!(spec.kind, NodeKind::Outlet) {
            let idx = id_map[spec.id.as_str()];
            network.set_pressure(idx, 0.0);
        }
    }

    // ── 5. Solve ─────────────────────────────────────────────────────────────
    println!("DATA: Solving...");
    let config = SolverConfig {
        tolerance: 1e-8,
        max_iterations: 200,
    };
    let solver = NetworkSolver::with_config(config);
    let solution = solver.solve(&NetworkProblem::new(network))?;
    println!("DATA: Solution converged.");

    // ── 6. Build AnalysisOverlay from Flow Rates ─────────────────────────────
    println!("VISUALIZATION: Building AnalysisOverlay...");

    // Map edge graph indices → channel schematic IDs
    // Channel IDs in the schematic are integers; edge IDs are "chan_<id>"
    let mut edge_flow_data = HashMap::<usize, f64>::new();
    for (eidx, &q) in solution.flow_rates.iter().enumerate() {
        if let Some(edge) = solution
            .graph
            .edge_weight(petgraph::graph::EdgeIndex::new(eidx))
        {
            if let Ok(id) = edge.id.trim_start_matches("chan_").parse::<usize>() {
                edge_flow_data.insert(id, q.abs());
            }
        }
    }

    // Node pressure data
    let mut node_pressure_data = HashMap::<usize, f64>::new();
    for (nidx, &p) in solution.pressures.iter().enumerate() {
        if let Some(node) = solution
            .graph
            .node_weight(petgraph::graph::NodeIndex::new(nidx))
        {
            if let Ok(id) = node.id.trim_start_matches("node_").parse::<usize>() {
                node_pressure_data.insert(id, p);
            }
        }
    }

    let max_flow = edge_flow_data.values().cloned().fold(0.0_f64, f64::max);

    // FlowRate overlay with Viridis colormap
    let overlay = AnalysisOverlay::new(AnalysisField::FlowRate, ColormapKind::Viridis)
        .with_edge_data(edge_flow_data)
        .with_node_data(node_pressure_data);

    // ── 7. Render ────────────────────────────────────────────────────────────
    let output_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("outputs")
        .join("geometry_integration");
    fs::create_dir_all(&output_dir)?;

    let renderer = create_plotters_renderer();
    let mut render_config = RenderConfig::default();
    render_config.title = format!("Flow Distribution (Max Q = {:.2e} m³/s)", max_flow);
    render_config.show_axes = true;
    render_config.show_grid = false;

    renderer.render_analysis(
        &system,
        output_dir.join("flow_analysis.png").to_str().unwrap(),
        &render_config,
        &overlay,
    )?;
    println!("✅ Rendered flow_analysis.png with AnalysisOverlay (Viridis colormap)");

    // ── 8. JSON Export ───────────────────────────────────────────────────────
    let export_data = serde_json::json!({
        "geometry": system,
        "results": {
            "max_flow_m3s": max_flow,
            "boundary_conditions": {
                "inlet_pressure_pa": 1000.0,
                "outlet_pressure_pa": 0.0
            }
        }
    });

    let json_path = output_dir.join("simulation_results.json");
    fs::write(&json_path, serde_json::to_string_pretty(&export_data)?)?;
    println!("✅ Exported results to {}", json_path.display());

    Ok(())
}
