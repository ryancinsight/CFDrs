//! Medical-Grade Millifluidic CFD Screening
//!
//! Comprehensive example combining all medical-relevant analyses in a single
//! device geometry sized for a 96-well plate (ANSI/SLAS 1-2004):
//!
//! 1. Generate a realistic serpentine-bifurcation geometry
//! 2. Solve blood flow (Carreau-Yasuda non-Newtonian)
//! 3. Compute and visualize:
//!    - Hemolysis index (Giersiepen–Wurzinger)
//!    - Wall shear stress distribution
//!    - FDA shear limit violations
//!    - Cavitation risk assessment
//!    - Flow rate distribution
//!    - Pressure field
//! 4. Render 6 colored overlays + 1 plain schematic
//! 5. Export full JSON report with per-channel data
//!
//! This example targets the Sonodynamic Therapy (SDT) use case where
//! blood-contacting millifluidic devices must simultaneously achieve
//! controlled cavitation while minimizing hemolysis.
//!
//! Run with:
//! `cargo run -p cfd-1d --example medical_millifluidic_screening`

use cfd_1d::network::{EdgeProperties, Network, NetworkBuilder};
use cfd_1d::solver::{NetworkProblem, NetworkSolver, SolverConfig};
use cfd_1d::BloodShearLimits;
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
    println!("🏥 Medical-Grade Millifluidic CFD Screening");
    println!("============================================\n");

    let out = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("outputs");
    fs::create_dir_all(out.join("medical_screening"))?;

    // ── 1. Geometry: serpentine bifurcation within 96-well plate ──────────────
    println!("1. Generating serpentine bifurcation geometry...");
    let box_dims = (90.0, 42.0); // mm — within ANSI/SLAS 1-2004 constraints
    let splits = vec![SplitType::Bifurcation, SplitType::Bifurcation];
    let mut geo_config = GeometryConfig::default();
    geo_config.channel_width = 1.0; // 1 mm — millifluidic scale
    let serpentine = smooth_serpentine();
    let channel_type = ChannelTypeConfig::AllSerpentine(serpentine);

    let system = create_geometry(box_dims, &splits, &geo_config, &channel_type);
    println!(
        "   Topology: {} nodes, {} channels",
        system.nodes.len(),
        system.channels.len()
    );
    println!(
        "   Footprint: {}×{} mm (96-well plate)",
        box_dims.0, box_dims.1
    );

    // ── 2. Plain schematic ───────────────────────────────────────────────────
    println!("2. Rendering plain schematic...");
    plot_geometry(
        &system,
        out.join("medical_screening/schematic.png")
            .to_str()
            .unwrap(),
    )?;

    // ── 3. Convert & build network ───────────────────────────────────────────
    println!("3. Building 1D network with Carreau-Yasuda blood...");
    let (node_specs, channel_specs) = convert_geometry_to_specs(&system);
    let blood = CarreauYasuda::<f64>::blood();

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
        let from = id_map[spec.from.as_str()];
        let to = id_map[spec.to.as_str()];
        let eidx = builder.add_edge(from, to, edge);
        edge_id_map.insert(spec.id.as_str().to_string(), eidx);
    }

    let graph = builder.build()?;
    let mut network = Network::new(graph, blood.clone());

    for spec in &channel_specs {
        let eidx = edge_id_map[spec.id.as_str()];
        network.add_edge_properties(eidx, EdgeProperties::from(spec));
    }

    // ── 4. Boundary conditions ───────────────────────────────────────────────
    println!("4. Applying boundary conditions...");
    let inlet_pressure = 5000.0; // 5 kPa
    let outlet_pressure = 0.0;

    for spec in &node_specs {
        let idx = id_map[spec.id.as_str()];
        match spec.kind {
            NodeKind::Inlet => network.set_pressure(idx, inlet_pressure),
            NodeKind::Outlet => network.set_pressure(idx, outlet_pressure),
            _ => {}
        }
    }

    // ── 5. Solve ─────────────────────────────────────────────────────────────
    println!("5. Solving...");
    let solver = NetworkSolver::with_config(SolverConfig {
        tolerance: 1e-8,
        max_iterations: 200,
    });
    let solution = solver.solve(&NetworkProblem::new(network))?;
    println!("   Converged.");

    // ── 6. Compute all analysis fields ───────────────────────────────────────
    println!("6. Computing medical analysis fields...");
    let hemolysis_model = HemolysisModel::giersiepen_standard();
    let vapor_pressure = 6300.0; // Pa — blood at 37°C

    let mut edge_hemolysis = HashMap::<usize, f64>::new();
    let mut edge_shear = HashMap::<usize, f64>::new();
    let mut edge_flow = HashMap::<usize, f64>::new();
    let mut edge_velocity = HashMap::<usize, f64>::new();
    let mut edge_cavitation = HashMap::<usize, f64>::new();
    let mut edge_pressure = HashMap::<usize, f64>::new();
    let mut node_pressure = HashMap::<usize, f64>::new();

    // FDA shear limit thresholds
    let fda_limits = BloodShearLimits::<f64>::fda_conservative_whole_blood();
    let mut fda_violations = Vec::new();

    // Node pressures
    for (nidx, &p) in &solution.pressures {
        if let Some(node) = solution.graph.node_weight(*nidx) {
            if let Ok(id) = node.id.trim_start_matches("node_").parse::<usize>() {
                node_pressure.insert(id, p);
            }
        }
    }

    let mut max_hi = 0.0_f64;
    let mut max_wss = 0.0_f64;
    let mut total_hi = 0.0_f64;

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
            let shear_rate = p
                .hydraulic_diameter
                .map(|d| 8.0 * velocity / d)
                .unwrap_or(0.0);
            let viscosity = blood
                .viscosity_at_shear(shear_rate, 310.15, 101_325.0)
                .unwrap_or(3.5e-3);
            let wall_shear = viscosity * shear_rate;

            // Hemolysis
            let exposure_time = if velocity > 1e-12 {
                p.length / velocity
            } else {
                0.0
            };
            let hi = hemolysis_model
                .damage_index(wall_shear, exposure_time)
                .unwrap_or(0.0);

            // Cavitation
            let (src, dst) = solution.graph.edge_endpoints(*eidx).unwrap();
            let p_src = *solution.pressures.get(&src).unwrap_or(&inlet_pressure);
            let p_dst = *solution.pressures.get(&dst).unwrap_or(&outlet_pressure);
            let local_p = (p_src + p_dst) / 2.0;

            let cav = CavitationNumber {
                reference_pressure: local_p,
                vapor_pressure,
                density: blood.density,
                velocity,
            };
            let sigma = cav.calculate();

            // FDA check
            if wall_shear > fda_limits.max_wall_shear_stress_pa {
                fda_violations.push((chan_id, wall_shear));
            }

            edge_hemolysis.insert(chan_id, hi);
            edge_shear.insert(chan_id, wall_shear);
            edge_flow.insert(chan_id, q.abs());
            edge_velocity.insert(chan_id, velocity);
            edge_pressure.insert(chan_id, local_p);
            // Cavitation risk: 1/σ (higher = more dangerous)
            edge_cavitation.insert(chan_id, if sigma > 1e-6 { 1.0 / sigma } else { 1e6 });

            max_hi = max_hi.max(hi);
            max_wss = max_wss.max(wall_shear);
            total_hi += hi;
        }
    }

    println!("   Max hemolysis index:      {:.4e}", max_hi);
    println!("   Cumulative hemolysis:     {:.4e}", total_hi);
    println!("   Max wall shear stress:    {:.2} Pa", max_wss);
    println!(
        "   FDA violations:           {} of {} channels",
        fda_violations.len(),
        edge_shear.len()
    );
    println!(
        "   FDA limit (τ_w):          {:.1} Pa",
        fda_limits.max_wall_shear_stress_pa
    );

    if !fda_violations.is_empty() {
        println!("\n   ⚠ FDA Shear Limit Violations:");
        for (id, wss) in &fda_violations {
            println!(
                "     Channel {}: τ_w = {:.2} Pa (limit: {:.1} Pa, {:.1}× exceedance)",
                id,
                wss,
                fda_limits.max_wall_shear_stress_pa,
                wss / fda_limits.max_wall_shear_stress_pa
            );
        }
    }

    // ── 7. Render all overlays ───────────────────────────────────────────────
    println!("\n7. Rendering analysis overlays...");
    let renderer = create_plotters_renderer();

    let overlays: Vec<(&str, AnalysisOverlay, String)> = vec![
        (
            "hemolysis_index.png",
            AnalysisOverlay::new(
                AnalysisField::Custom("Hemolysis Index (Giersiepen)".into()),
                ColormapKind::BlueRed,
            )
            .with_edge_data(edge_hemolysis.clone())
            .with_node_data(node_pressure.clone()),
            format!("Hemolysis Index (max HI = {:.2e})", max_hi),
        ),
        (
            "wall_shear_stress.png",
            AnalysisOverlay::new(AnalysisField::WallShearStress, ColormapKind::BlueRed)
                .with_edge_data(edge_shear.clone())
                .with_node_data(node_pressure.clone()),
            format!("Wall Shear Stress (max = {:.2} Pa)", max_wss),
        ),
        (
            "fda_shear_screening.png",
            {
                // Binary overlay: 0 = safe, 1 = violation
                let fda_binary: HashMap<usize, f64> = edge_shear
                    .iter()
                    .map(|(&k, &v)| {
                        let ratio = v / fda_limits.max_wall_shear_stress_pa;
                        (k, ratio.min(3.0))
                    })
                    .collect();
                AnalysisOverlay::new(
                    AnalysisField::Custom("FDA Shear Exceedance Ratio".into()),
                    ColormapKind::BlueRed,
                )
                .with_edge_data(fda_binary)
                .with_node_data(node_pressure.clone())
            },
            format!(
                "FDA Shear Screening (limit: {:.0} Pa, {} violations)",
                fda_limits.max_wall_shear_stress_pa,
                fda_violations.len()
            ),
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
            "flow_rate.png",
            AnalysisOverlay::new(AnalysisField::FlowRate, ColormapKind::Viridis)
                .with_edge_data(edge_flow.clone())
                .with_node_data(node_pressure.clone()),
            "Flow Rate Distribution [m³/s]".to_string(),
        ),
        (
            "pressure.png",
            AnalysisOverlay::new(AnalysisField::Pressure, ColormapKind::Viridis)
                .with_edge_data(edge_pressure.clone())
                .with_node_data(node_pressure.clone()),
            "Pressure Distribution [Pa]".to_string(),
        ),
    ];

    for (filename, overlay, title) in &overlays {
        let path = out.join(format!("medical_screening/{filename}"));
        let mut rc = RenderConfig::default();
        rc.title = title.clone();
        rc.show_axes = true;
        rc.show_grid = false;
        renderer.render_analysis(&system, path.to_str().unwrap(), &rc, overlay)?;
        println!("   ✓ {filename}");
    }

    // ── 8. JSON Export ───────────────────────────────────────────────────────
    println!("\n8. Exporting comprehensive JSON report...");

    let channel_report: Vec<serde_json::Value> = edge_hemolysis
        .keys()
        .map(|&chan_id| {
            serde_json::json!({
                "channel_id": chan_id,
                "hemolysis_index": edge_hemolysis.get(&chan_id),
                "wall_shear_stress_pa": edge_shear.get(&chan_id),
                "flow_rate_m3s": edge_flow.get(&chan_id),
                "velocity_m_s": edge_velocity.get(&chan_id),
                "local_pressure_pa": edge_pressure.get(&chan_id),
                "cavitation_risk": edge_cavitation.get(&chan_id),
                "fda_violation": edge_shear.get(&chan_id)
                    .map(|&wss| wss > fda_limits.max_wall_shear_stress_pa)
                    .unwrap_or(false)
            })
        })
        .collect();

    let results = serde_json::json!({
        "analysis": "medical_millifluidic_screening",
        "device": {
            "application": "Sonodynamic Therapy (SDT)",
            "standard": "ANSI/SLAS 1-2004 (96-well)",
            "footprint_mm": [box_dims.0, box_dims.1],
            "channel_width_mm": geo_config.channel_width,
            "topology": {
                "num_nodes": system.nodes.len(),
                "num_channels": system.channels.len(),
                "split_pattern": "bifurcation × 2"
            }
        },
        "fluid": {
            "model": "Carreau-Yasuda Blood",
            "density_kg_m3": blood.density,
            "viscosity_inf_pas": blood.viscosity_inf,
            "viscosity_zero_pas": blood.viscosity_zero,
            "temperature_k": 310.15,
            "vapor_pressure_pa": vapor_pressure
        },
        "boundary_conditions": {
            "inlet_pressure_pa": inlet_pressure,
            "outlet_pressure_pa": outlet_pressure
        },
        "hemolysis": {
            "model": "Giersiepen-Wurzinger Power Law",
            "constants": {"C": 3.62e-5, "alpha": 2.416, "beta": 0.785},
            "max_hemolysis_index": max_hi,
            "cumulative_hemolysis": total_hi
        },
        "fda_screening": {
            "max_wall_shear_limit_pa": fda_limits.max_wall_shear_stress_pa,
            "num_violations": fda_violations.len(),
            "max_wall_shear_observed_pa": max_wss,
            "violations": fda_violations.iter().map(|(id, wss)| {
                serde_json::json!({
                    "channel_id": id,
                    "wall_shear_stress_pa": wss,
                    "exceedance_ratio": wss / fda_limits.max_wall_shear_stress_pa
                })
            }).collect::<Vec<_>>()
        },
        "channels": channel_report,
        "outputs": {
            "schematic": "medical_screening/schematic.png",
            "hemolysis": "medical_screening/hemolysis_index.png",
            "wall_shear": "medical_screening/wall_shear_stress.png",
            "fda_screening": "medical_screening/fda_shear_screening.png",
            "cavitation": "medical_screening/cavitation_risk.png",
            "flow_rate": "medical_screening/flow_rate.png",
            "pressure": "medical_screening/pressure.png"
        }
    });

    let json_path = out.join("medical_screening/report.json");
    fs::write(&json_path, serde_json::to_string_pretty(&results)?)?;
    println!("   Report exported to {}", json_path.display());

    println!("\n✅ Medical millifluidic screening complete — 7 plots + 1 JSON report");
    Ok(())
}

/// Convert `ChannelSystem` → (`NodeSpec`, `ChannelSpec`) pairs.
///
/// Inlet/outlet classification by x-coordinate. Hele-Shaw rectangular resistance.
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

    let mu = 3.5e-3_f64; // Blood apparent viscosity at mid-shear

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
