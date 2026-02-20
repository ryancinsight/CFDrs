//! 1D Blood Flow in a Bifurcation (Carreau-Yasuda Model)
//!
//! This example demonstrates simulating blood flow in a bifurcating artery
//! using the Carreau-Yasuda non-Newtonian viscosity model.
//!
//! Network topology and physical geometry are defined entirely via
//! `cfd-schematics` `NodeSpec`/`ChannelSpec`, then converted to solver types
//! via the canonical `From<&ChannelSpec>` bridge. No `cfd_1d::channel` imports
//! are required in this example.
//!
//! Geometry:
//! - Parent vessel: D=4mm, L=50mm
//! - Daughter vessel 1: D=3mm, L=40mm
//! - Daughter vessel 2: D=2.5mm, L=30mm
//!
//! Conditions:
//! - Inlet Flow: 5 mL/s (pulsatile average)
//! - Outlet Pressure: 100 mmHg = 13332.2 Pa

use cfd_1d::network::{EdgeProperties, Network, NetworkBuilder};
use cfd_1d::solver::{NetworkProblem, NetworkSolver, SolverConfig};
use cfd_core::compute::solver::Solver;
use cfd_core::physics::fluid::non_newtonian::CarreauYasuda;
use cfd_core::physics::fluid::FluidTrait;
use scheme::domain::model::{ChannelSpec, NodeKind, NodeSpec};
use scheme::visualizations::{
    analysis_field::{AnalysisField, AnalysisOverlay, ColormapKind},
    traits::{RenderConfig, SchematicRenderer},
    PlottersRenderer,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // ── 1. Fluid Properties ──────────────────────────────────────────────────
    let blood = CarreauYasuda::<f64>::blood();
    println!("Fluid: {}", blood.name());
    println!("Density: {} kg/m³", blood.density);
    println!("Viscosity (inf): {} Pa·s", blood.viscosity_inf);
    println!("Viscosity (0): {} Pa·s", blood.viscosity_zero);

    // ── 2. Network Topology via ChannelSpec/NodeSpec ─────────────────────────
    // All physical geometry lives in ChannelSpec — no cfd_1d::channel imports needed.
    let nodes = vec![
        NodeSpec::new("Inlet", NodeKind::Inlet),
        NodeSpec::new("Bifurcation", NodeKind::Junction),
        NodeSpec::new("Outlet1", NodeKind::Outlet),
        NodeSpec::new("Outlet2", NodeKind::Outlet),
    ];

    // Hagen-Poiseuille resistance: R = 128·μ·L / (π·D⁴)
    // Using apparent viscosity ≈ 3.5e-3 Pa·s (mid-shear blood)
    let mu = 3.5e-3_f64;
    let r_parent = 128.0 * mu * 0.05 / (std::f64::consts::PI * 0.004_f64.powi(4));
    let r_d1 = 128.0 * mu * 0.04 / (std::f64::consts::PI * 0.003_f64.powi(4));
    let r_d2 = 128.0 * mu * 0.03 / (std::f64::consts::PI * 0.0025_f64.powi(4));

    let channels = vec![
        ChannelSpec::new_pipe("Parent", "Inlet", "Bifurcation", 0.05, 0.004, r_parent, 0.0),
        ChannelSpec::new_pipe("Daughter1", "Bifurcation", "Outlet1", 0.04, 0.003, r_d1, 0.0),
        ChannelSpec::new_pipe("Daughter2", "Bifurcation", "Outlet2", 0.03, 0.0025, r_d2, 0.0),
    ];

    // ── 3. Build Network ─────────────────────────────────────────────────────
    let mut builder = NetworkBuilder::<f64>::new();

    let mut node_indices = std::collections::HashMap::new();
    for node_spec in &nodes {
        let idx = match node_spec.kind {
            NodeKind::Inlet => builder.add_inlet(node_spec.id.as_str().to_string()),
            NodeKind::Outlet => builder.add_outlet(node_spec.id.as_str().to_string()),
            NodeKind::Junction | NodeKind::Reservoir => {
                builder.add_junction(node_spec.id.as_str().to_string())
            }
        };
        node_indices.insert(node_spec.id.as_str().to_string(), idx);
    }

    let mut edge_indices = std::collections::HashMap::new();
    for ch in &channels {
        let from = node_indices[ch.from.as_str()];
        let to = node_indices[ch.to.as_str()];
        let eidx = builder.connect_with_pipe(from, to, ch.id.as_str().to_string());
        edge_indices.insert(ch.id.as_str().to_string(), eidx);
    }

    let graph = builder.build()?;
    let mut network = Network::new(graph, blood.clone());

    // ── 4. Assign EdgeProperties from ChannelSpec ────────────────────────────
    for ch in &channels {
        let props = EdgeProperties::from(ch);
        network.add_edge_properties(edge_indices[ch.id.as_str()], props);
    }

    // ── 5. Boundary Conditions ───────────────────────────────────────────────
    let inlet_idx = node_indices["Inlet"];
    let outlet1_idx = node_indices["Outlet1"];
    let outlet2_idx = node_indices["Outlet2"];
    let parent_eidx = edge_indices["Parent"];
    let d1_eidx = edge_indices["Daughter1"];
    let d2_eidx = edge_indices["Daughter2"];

    network.set_neumann_flow(inlet_idx, 5e-6);
    let p_out = 13332.2_f64;
    network.set_pressure(outlet1_idx, p_out);
    network.set_pressure(outlet2_idx, p_out);

    // Initial flow rate guess for non-Newtonian viscosity start
    network.set_flow_rate(parent_eidx, 5e-6);
    network.set_flow_rate(d1_eidx, 3e-6);
    network.set_flow_rate(d2_eidx, 2e-6);
    network.update_resistances()?;

    // ── 6. Solve ─────────────────────────────────────────────────────────────
    let config = SolverConfig { tolerance: 1e-6, max_iterations: 100 };
    let solver = NetworkSolver::<f64, CarreauYasuda<f64>>::with_config(config);
    println!("Solving network...");
    let solution = solver.solve(&NetworkProblem::new(network))?;

    // ── 7. Print Results ─────────────────────────────────────────────────────
    println!("\nNode Pressures:");
    let mut node_pressures = std::collections::HashMap::<usize, f64>::new();
    for idx in solution.graph.node_indices() {
        let n = solution.graph.node_weight(idx).unwrap();
        let p = *solution.pressures.get(&idx).unwrap_or(&0.0);
        println!("  {}: {:.2} Pa ({:.2} mmHg)", n.id, p, p / 133.322);
        node_pressures.insert(idx.index(), p);
    }

    println!("\nEdge Flow Rates & Shear:");
    let mut edge_flow_rates = std::collections::HashMap::<usize, f64>::new();
    let mut edge_shear = std::collections::HashMap::<usize, f64>::new();

    for idx in solution.graph.edge_indices() {
        let e = solution.graph.edge_weight(idx).unwrap();
        let q = *solution.flow_rates.get(&idx).unwrap_or(&0.0);
        let props = solution.properties.get(&idx);

        let (v, shear_rate) = if let Some(p) = props {
            let vel = q.abs() / p.area;
            let sr = p.hydraulic_diameter
                .map(|d| 8.0 * vel / d)
                .unwrap_or(0.0);
            (vel, sr)
        } else {
            (0.0, 0.0)
        };

        let viscosity = blood.viscosity_at_shear(shear_rate, 310.15, 101325.0)?;
        let wall_shear = viscosity * shear_rate;

        println!("  {}: {:.4} mL/s", e.id, q * 1e6);
        println!("     Velocity: {:.2} cm/s", v * 100.0);
        println!("     Shear Rate: {:.1} 1/s", shear_rate);
        println!("     Apparent Viscosity: {:.2} mPa·s", viscosity * 1000.0);
        println!("     Wall Shear Stress: {:.2} Pa", wall_shear);

        edge_flow_rates.insert(idx.index(), q.abs());
        edge_shear.insert(idx.index(), wall_shear);
    }

    // ── 8. Visualization with AnalysisOverlay ────────────────────────────────
    // Build a simple linear schematic for visualization
    // (blood_bifurcation doesn't have a ChannelSystem; we export JSON only here
    //  and note that geometry_integration_demo shows the full schematic path)

    // ── 9. JSON Export ───────────────────────────────────────────────────────
    use std::fs;
    use std::path::Path;

    let output_dir = Path::new("crates/cfd-1d/outputs");
    if !output_dir.exists() {
        fs::create_dir_all(output_dir)?;
    }

    let results = serde_json::json!({
        "fluid": "CarreauYasuda Blood",
        "nodes": solution.graph.node_indices().map(|idx| {
            let n = solution.graph.node_weight(idx).unwrap();
            serde_json::json!({
                "id": n.id,
                "pressure_pa": solution.pressures.get(&idx).unwrap_or(&0.0),
                "pressure_mmhg": solution.pressures.get(&idx).unwrap_or(&0.0) / 133.322,
                "type": format!("{:?}", n.node_type)
            })
        }).collect::<Vec<_>>(),
        "edges": solution.graph.edge_indices().map(|idx| {
            let e = solution.graph.edge_weight(idx).unwrap();
            let q = solution.flow_rates.get(&idx).unwrap_or(&0.0);
            let shear = edge_shear.get(&idx.index()).unwrap_or(&0.0);
            serde_json::json!({
                "id": e.id,
                "flow_rate_m3s": q,
                "flow_rate_mls": q * 1e6,
                "wall_shear_stress_pa": shear,
                "resistance": solution.properties.get(&idx).map(|p| p.resistance).unwrap_or(0.0)
            })
        }).collect::<Vec<_>>()
    });

    let json_path = output_dir.join("blood_bifurcation_results.json");
    fs::write(&json_path, serde_json::to_string_pretty(&results)?)?;
    println!("\n✅ Exported results to {}", json_path.display());

    Ok(())
}
