//! Medical-Grade Millifluidic CFD Screening -- 2D Navier-Stokes
//!
//! 2D counterpart of `cfd-1d::medical_millifluidic_screening`. Uses the same
//! millifluidic device topology (inlet -> bifurcation splits -> serpentine
//! parallel branches -> outlets) sized for a 96-well plate (ANSI/SLAS 1-2004),
//! but resolves the 2D velocity, pressure, and concentration fields using
//! `SerpentineSolver2D` instead of a lumped 1D resistance network.
//!
//! Analyses all six medical-relevant fields:
//!   1. Hemolysis index (Giersiepen-Wurzinger)
//!   2. Wall shear stress distribution
//!   3. FDA shear limit violations
//!   4. Cavitation risk assessment
//!   5. Flow rate distribution
//!   6. Pressure field
//!
//! Pipeline:
//!   cfd-schematics `create_geometry` (AllSerpentine) ->
//!   `SerpentineSolver2D` (Carreau-Yasuda blood, 2D N-S) ->
//!   6x AnalysisOverlay renders -> JSON export
//!
//! Run with:
//! `cargo run -p cfd-2d --example medical_screening_2d`

use cfd_2d::solvers::ns_fvm_2d::BloodModel;
use cfd_2d::solvers::serpentine_flow::{SerpentineGeometry, SerpentineSolver2D};
use cfd_core::physics::cavitation::CavitationNumber;
use cfd_core::physics::fluid::blood::CarreauYasudaBlood;
use cfd_core::physics::hemolysis::HemolysisModel;
use cfd_schematics::config::presets::smooth_serpentine;
use cfd_schematics::config::{ChannelTypeConfig, GeometryConfig};
use cfd_schematics::geometry::generator::create_geometry;
use cfd_schematics::geometry::SplitType;
use cfd_schematics::plot_geometry;
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
    println!("Medical-Grade Millifluidic CFD Screening -- 2D N-S");
    println!("===================================================\n");

    let out = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("outputs");
    fs::create_dir_all(out.join("medical_screening"))?;

    // -- 1. Geometry: serpentine bifurcation within 96-well plate --
    println!("1. Generating serpentine bifurcation geometry...");
    let box_dims = (90.0, 42.0); // mm -- within ANSI/SLAS 1-2004 constraints
    let splits = vec![SplitType::Bifurcation, SplitType::Bifurcation];
    let mut geo_config = GeometryConfig::default();
    geo_config.channel_width = 1.0; // 1 mm -- millifluidic scale
    geo_config.channel_height = 0.5; // 500 um depth
    let serpentine = smooth_serpentine();
    let channel_type = ChannelTypeConfig::AllSerpentine(serpentine);

    let system = create_geometry(box_dims, &splits, &geo_config, &channel_type);
    println!(
        "   Topology: {} nodes, {} channels",
        system.nodes.len(),
        system.channels.len()
    );
    println!(
        "   Footprint: {}x{} mm (96-well plate)",
        box_dims.0, box_dims.1
    );

    // -- 2. Plain schematic --
    println!("2. Rendering plain schematic...");
    plot_geometry(
        &system,
        out.join("medical_screening/schematic.png")
            .to_str()
            .unwrap(),
    )?;

    // -- 3. Flow distribution --
    println!("3. Computing bifurcation flow distribution...");
    let channel_levels = compute_channel_levels(&system);
    let max_level = channel_levels.values().copied().max().unwrap_or(0);

    let ch_w = geo_config.channel_width / 1000.0; // m
    let ch_h = geo_config.channel_height / 1000.0; // m
    let cross_area = ch_w * ch_h;
    let d_h = 2.0 * ch_w * ch_h / (ch_w + ch_h);

    // -- 4. Run 2D Navier-Stokes solver --
    println!("4. Running SerpentineSolver2D (Carreau-Yasuda blood)...");
    let blood = CarreauYasudaBlood::<f64>::normal_blood();
    let density = blood.density;
    let blood_model = BloodModel::CarreauYasuda(blood);

    let serp_geom = SerpentineGeometry::new(ch_w, ch_h, 2.0e-3, 0.5e-3, 5);
    let ref_length = serp_geom.total_length();

    let u_inlet = 0.1; // 100 mm/s -- device inlet velocity
    let u_branch = u_inlet / 2.0_f64.powi(max_level as i32);

    let mut solver = SerpentineSolver2D::new(serp_geom, blood_model, density, 60, 30);
    let sol = solver.solve(u_branch, 1e-9, 0.0, 1.0)?;

    let ref_wss = sol.pressure_drop * d_h / (4.0 * ref_length);
    println!(
        "   dP = {:.2} Pa, ref WSS = {:.2} Pa",
        sol.pressure_drop, ref_wss
    );

    // -- 5. Compute all six analysis fields --
    println!("5. Computing medical analysis fields...");
    let hemolysis_model = HemolysisModel::giersiepen_standard();
    let vapor_pressure = 6300.0; // Pa -- blood at 37C
    let fda_limit_wss = 150.0; // Pa -- FDA conservative guidance for whole blood

    let mut edge_hemolysis = HashMap::<usize, f64>::new();
    let mut edge_shear = HashMap::<usize, f64>::new();
    let mut edge_fda = HashMap::<usize, f64>::new();
    let mut edge_cavitation = HashMap::<usize, f64>::new();
    let mut edge_flow = HashMap::<usize, f64>::new();
    let mut edge_pressure = HashMap::<usize, f64>::new();
    let mut node_pressure = HashMap::<usize, f64>::new();

    let mut max_hi = 0.0_f64;
    let mut max_wss = 0.0_f64;
    let mut total_hi = 0.0_f64;
    let mut fda_violations = Vec::new();

    // Node pressure: linear from inlet to outlet
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
    let inlet_pressure = 5000.0; // 5 kPa -- same BC as cfd-1d medical screening
    for node in &system.nodes {
        let frac = (node.point.0 - min_x) / (max_x - min_x).max(1e-12);
        node_pressure.insert(node.id, inlet_pressure * (1.0 - frac));
    }

    for channel in &system.channels {
        let level = channel_levels.get(&channel.id).copied().unwrap_or(0);
        let velocity = u_inlet / 2.0_f64.powi(level as i32);
        let length_m = channel_length_m(channel, &system);
        let flow_rate = velocity * cross_area;

        // WSS scaled from 2D reference
        let wall_shear = ref_wss * (velocity / u_branch);

        // Exposure time
        let exposure_time = if velocity > 1e-12 {
            length_m / velocity
        } else {
            0.0
        };

        // Hemolysis index
        let hi = hemolysis_model
            .damage_index(wall_shear, exposure_time)
            .unwrap_or(0.0);

        // Channel-average pressure
        let p_from = node_pressure
            .get(&channel.from_node)
            .copied()
            .unwrap_or(0.0);
        let p_to = node_pressure.get(&channel.to_node).copied().unwrap_or(0.0);
        let p_avg = (p_from + p_to) / 2.0;

        // Cavitation number
        let cav = CavitationNumber {
            reference_pressure: p_avg.max(100.0),
            vapor_pressure,
            density,
            velocity,
        };
        let sigma = cav.calculate();

        // FDA shear exceedance ratio
        let fda_ratio = (wall_shear / fda_limit_wss).min(3.0);
        if wall_shear > fda_limit_wss {
            fda_violations.push((channel.id, wall_shear));
        }

        edge_hemolysis.insert(channel.id, hi);
        edge_shear.insert(channel.id, wall_shear);
        edge_fda.insert(channel.id, fda_ratio);
        edge_cavitation.insert(channel.id, 1.0 / sigma.max(1e-12));
        edge_flow.insert(channel.id, flow_rate);
        edge_pressure.insert(channel.id, p_avg);

        max_hi = max_hi.max(hi);
        max_wss = max_wss.max(wall_shear);
        total_hi += hi;
    }

    println!("   Max hemolysis index:      {:.4e}", max_hi);
    println!("   Cumulative hemolysis:     {:.4e}", total_hi);
    println!("   Max wall shear stress:    {:.2} Pa", max_wss);
    println!(
        "   FDA violations:           {} of {} channels",
        fda_violations.len(),
        edge_shear.len()
    );
    println!("   FDA limit (tau_w):        {:.1} Pa", fda_limit_wss);

    if !fda_violations.is_empty() {
        println!("\n   FDA Shear Limit Violations:");
        for (id, wss) in &fda_violations {
            println!(
                "     Channel {}: tau_w = {:.2} Pa ({:.1}x limit)",
                id,
                wss,
                wss / fda_limit_wss
            );
        }
    }

    // -- 6. Render all six overlays --
    println!("\n6. Rendering analysis overlays...");
    let renderer = create_plotters_renderer();

    let overlays: Vec<(&str, AnalysisOverlay)> = vec![
        (
            "hemolysis_index",
            AnalysisOverlay::new(
                AnalysisField::Custom("Hemolysis Index (Giersiepen)".into()),
                ColormapKind::BlueRed,
            )
            .with_edge_data(edge_hemolysis.clone())
            .with_node_data(node_pressure.clone()),
        ),
        (
            "wall_shear_stress",
            AnalysisOverlay::new(AnalysisField::WallShearStress, ColormapKind::BlueRed)
                .with_edge_data(edge_shear.clone())
                .with_node_data(node_pressure.clone()),
        ),
        (
            "fda_shear_screening",
            AnalysisOverlay::new(
                AnalysisField::Custom("FDA Shear Exceedance Ratio".into()),
                ColormapKind::BlueRed,
            )
            .with_edge_data(edge_fda.clone())
            .with_node_data(node_pressure.clone()),
        ),
        (
            "cavitation_risk",
            AnalysisOverlay::new(
                AnalysisField::Custom("Cavitation Risk (1/sigma)".into()),
                ColormapKind::BlueRed,
            )
            .with_edge_data(edge_cavitation.clone())
            .with_node_data(node_pressure.clone()),
        ),
        (
            "flow_rate",
            AnalysisOverlay::new(AnalysisField::FlowRate, ColormapKind::Viridis)
                .with_edge_data(edge_flow.clone())
                .with_node_data(node_pressure.clone()),
        ),
        (
            "pressure",
            AnalysisOverlay::new(AnalysisField::Pressure, ColormapKind::Viridis)
                .with_edge_data(edge_pressure.clone())
                .with_node_data(node_pressure.clone()),
        ),
    ];

    for (name, overlay) in &overlays {
        let path = out.join(format!("medical_screening/{name}.png"));
        let mut rc = RenderConfig::default();
        rc.title = format!("Medical 2D -- {name}");
        rc.show_axes = true;
        rc.show_grid = false;
        renderer.render_analysis(&system, path.to_str().unwrap(), &rc, overlay)?;
        println!("   [ok] {name}.png");
    }

    // -- 7. JSON Export --
    println!("\n7. Exporting comprehensive JSON report...");
    let channel_report: Vec<serde_json::Value> = system
        .channels
        .iter()
        .map(|ch| {
            serde_json::json!({
                "channel_id": ch.id,
                "bifurcation_level": channel_levels.get(&ch.id),
                "hemolysis_index": edge_hemolysis.get(&ch.id),
                "wall_shear_stress_pa": edge_shear.get(&ch.id),
                "fda_exceedance_ratio": edge_fda.get(&ch.id),
                "flow_rate_m3s": edge_flow.get(&ch.id),
                "pressure_pa": edge_pressure.get(&ch.id),
                "cavitation_risk": edge_cavitation.get(&ch.id),
                "fda_violation": edge_shear.get(&ch.id)
                    .map(|&wss| wss > fda_limit_wss)
                    .unwrap_or(false)
            })
        })
        .collect();

    let results = serde_json::json!({
        "analysis": "medical_millifluidic_screening_2d",
        "solver": "SerpentineSolver2D (2D Navier-Stokes)",
        "device": {
            "application": "Sonodynamic Therapy (SDT)",
            "standard": "ANSI/SLAS 1-2004 (96-well)",
            "footprint_mm": [box_dims.0, box_dims.1],
            "channel_width_mm": geo_config.channel_width,
            "channel_height_mm": geo_config.channel_height,
            "topology": {
                "num_nodes": system.nodes.len(),
                "num_channels": system.channels.len(),
                "split_pattern": "bifurcation x 2"
            }
        },
        "fluid": {
            "model": "Carreau-Yasuda Blood",
            "density_kg_m3": density,
            "temperature_k": 310.15,
            "vapor_pressure_pa": vapor_pressure
        },
        "solver_2d": {
            "grid": "60 x 30",
            "u_inlet_ms": u_inlet,
            "u_branch_ms": u_branch,
            "pressure_drop_pa": sol.pressure_drop,
            "mixing_fraction": sol.mixing_fraction_outlet,
            "reference_wss_pa": ref_wss,
        },
        "hemolysis": {
            "model": "Giersiepen-Wurzinger Power Law",
            "constants": {"C": 3.62e-5, "alpha": 2.416, "beta": 0.785},
            "max_hemolysis_index": max_hi,
            "cumulative_hemolysis": total_hi
        },
        "fda_screening": {
            "max_wall_shear_limit_pa": fda_limit_wss,
            "num_violations": fda_violations.len(),
            "max_wall_shear_observed_pa": max_wss,
            "violations": fda_violations.iter().map(|(id, wss)| {
                serde_json::json!({
                    "channel_id": id,
                    "wall_shear_stress_pa": wss,
                    "exceedance_ratio": wss / fda_limit_wss
                })
            }).collect::<Vec<_>>()
        },
        "channels": channel_report,
    });

    let json_path = out.join("medical_screening/results.json");
    fs::write(&json_path, serde_json::to_string_pretty(&results)?)?;
    println!("   Exported to {}", json_path.display());

    println!("\nMedical screening complete -- 7 PNGs + JSON in medical_screening/");
    Ok(())
}

/// Compute the straight-line length of a channel in metres.
fn channel_length_m(
    channel: &cfd_schematics::geometry::types::Channel,
    system: &cfd_schematics::geometry::types::ChannelSystem,
) -> f64 {
    let from = system
        .nodes
        .iter()
        .find(|n| n.id == channel.from_node)
        .unwrap();
    let to = system
        .nodes
        .iter()
        .find(|n| n.id == channel.to_node)
        .unwrap();
    let dx = from.point.0 - to.point.0;
    let dy = from.point.1 - to.point.1;
    dx.hypot(dy) / 1000.0 // mm -> m
}

/// Assign bifurcation levels to channels via BFS from inlet nodes.
fn compute_channel_levels(
    system: &cfd_schematics::geometry::types::ChannelSystem,
) -> HashMap<usize, usize> {
    let min_x = system
        .nodes
        .iter()
        .map(|n| n.point.0)
        .fold(f64::INFINITY, f64::min);
    let inlet_ids: Vec<usize> = system
        .nodes
        .iter()
        .filter(|n| (n.point.0 - min_x).abs() < 1e-3)
        .map(|n| n.id)
        .collect();

    let mut levels = HashMap::new();
    let mut current = inlet_ids;
    let mut level = 0;

    while !current.is_empty() {
        let mut next = Vec::new();
        for &node_id in &current {
            for chan in &system.channels {
                if chan.from_node == node_id && !levels.contains_key(&chan.id) {
                    levels.insert(chan.id, level);
                    next.push(chan.to_node);
                }
            }
        }
        level += 1;
        current = next;
    }

    levels
}
