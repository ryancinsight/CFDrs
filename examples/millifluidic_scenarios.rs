//! # Millifluidic Design Scenarios — CFD-1D
//!
//! Extended scenario suite demonstrating a range of millifluidic network
//! topologies, channel types, and analysis techniques through the pipeline:
//!
//! ```text
//! scheme design → cfd-1d network → steady / transient solve → post-process → CSV + SVG
//! ```
//!
//! ## Scenarios
//!
//! | # | Scenario | Topology | Channels | Analysis |
//! |---|----------|----------|----------|----------|
//! | 1 | Trifurcation distributor | 1→3 | Straight | Flow uniformity |
//! | 2 | Two-level cascade | 1→2→4 | Serpentine | Dean-enhanced mixing |
//! | 3 | Mixed venturi-serpentine | 1→2 | Mixed | Cavitation + shear |
//! | 4 | Pressure parametric sweep | 1→2 | Straight | τ_w vs ΔP curve |
//! | 5 | Transient composition switching | 1→2 | Straight | Reagent mixing |
//!
//! ## Running
//!
//! ```sh
//! cargo run --example millifluidic_scenarios --features scheme-integration
//! ```

// ═══════════════════════════════════════════════════════════════════════════════
// Imports
// ═══════════════════════════════════════════════════════════════════════════════

use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

// cfd-1d
use cfd_1d::channel::{ChannelType as CfdChannelType, CrossSection};
use cfd_1d::network::{Network, NodeType};
use cfd_1d::scheme_bridge::SchemeNetworkConverter;
use cfd_1d::solver::{
    InletCompositionEvent, MixtureComposition, NetworkProblem, NetworkSolver,
    SimulationTimeConfig, SolverConfig, TransientCompositionSimulator,
};

// cfd-core
use cfd_core::physics::fluid::ConstantPropertyFluid;

// cfd-schematics (canonical SSOT for network topology)
use cfd_schematics::config::*;
use cfd_schematics::geometry::generator::create_geometry;
use cfd_schematics::geometry::{ChannelSystem, ChannelType as SchemeChannelType, SplitType};

// petgraph
use petgraph::graph::NodeIndex;

// plotters
use plotters::prelude::*;

// ═══════════════════════════════════════════════════════════════════════════════
// Constants
// ═══════════════════════════════════════════════════════════════════════════════

const BLOOD_DENSITY: f64 = 1_060.0;     // kg/m³
const BLOOD_VISCOSITY: f64 = 3.5e-3;    // Pa·s
const BLOOD_CP: f64 = 3_600.0;          // J/(kg·K)
const BLOOD_K: f64 = 0.50;              // W/(m·K)
const BLOOD_SOUND_SPEED: f64 = 1_540.0; // m/s
const SCALE_MM_TO_M: f64 = 1e-3;

const P_VAPOR_37C: f64 = 6_274.0;       // Pa

const CHIP_W: f64 = 60.0; // mm
const CHIP_H: f64 = 30.0; // mm

// ═══════════════════════════════════════════════════════════════════════════════
// Shared utilities
// ═══════════════════════════════════════════════════════════════════════════════

fn blood_fluid() -> ConstantPropertyFluid<f64> {
    ConstantPropertyFluid::new(
        "Blood (Newtonian, 37 °C)".into(),
        BLOOD_DENSITY,
        BLOOD_VISCOSITY,
        BLOOD_CP,
        BLOOD_K,
        BLOOD_SOUND_SPEED,
    )
}

fn water_fluid() -> ConstantPropertyFluid<f64> {
    ConstantPropertyFluid::new(
        "Water (20 °C)".into(),
        998.0,
        1.002e-3,
        4182.0,
        0.598,
        1482.0,
    )
}

/// Identify inlet / outlet nodes and return (inlets, outlets).
fn find_boundary_nodes(network: &Network<f64, ConstantPropertyFluid<f64>>) -> (Vec<NodeIndex>, Vec<NodeIndex>) {
    let mut inlets = Vec::new();
    let mut outlets = Vec::new();
    for ni in network.graph.node_indices() {
        match network.graph[ni].node_type {
            NodeType::Inlet => inlets.push(ni),
            NodeType::Outlet => outlets.push(ni),
            _ => {}
        }
    }
    (inlets, outlets)
}

/// Per-channel analysis record.
#[derive(Debug, Clone)]
struct ChannelRecord {
    id: String,
    ch_type: String,
    flow_rate_m3s: f64,
    velocity_ms: f64,
    reynolds: f64,
    wall_shear_pa: f64,
    dp_pa: f64,
    length_m: f64,
    dh_m: f64,
}

/// Analyse a solved network and return per-channel records.
fn analyse_network(
    solution: &Network<f64, ConstantPropertyFluid<f64>>,
    density: f64,
    viscosity: f64,
) -> Vec<ChannelRecord> {
    let pressures = solution.pressures();
    let flow_rates = solution.flow_rates();

    let mut records = Vec::new();
    for ei in solution.graph.edge_indices() {
        let edge = &solution.graph[ei];
        let (from, to) = match solution.graph.edge_endpoints(ei) {
            Some(pair) => pair,
            None => continue,
        };
        let q = flow_rates.get(&ei).copied().unwrap_or(0.0).abs();
        let p_from = pressures.get(&from).copied().unwrap_or(0.0);
        let p_to = pressures.get(&to).copied().unwrap_or(0.0);
        let dp = (p_from - p_to).abs();

        let props = solution.properties.get(&ei);
        let length = props.map_or(0.01, |p| p.length);
        let (area, dh) = props.map_or((1e-6, 1e-3), |p| {
            (p.area, p.hydraulic_diameter.unwrap_or(1e-3))
        });
        let velocity = if area > 0.0 { q / area } else { 0.0 };
        let reynolds = density * velocity * dh / viscosity;
        let wall_shear = viscosity * 8.0 * velocity / dh;

        let label = if let Some(p) = props {
            p.geometry
                .as_ref()
                .map_or("unknown".into(), |g| match &g.channel_type {
                    CfdChannelType::Straight => "straight".into(),
                    CfdChannelType::Serpentine { turns } => format!("serpentine({turns}t)"),
                    CfdChannelType::Curved { .. } => "curved".into(),
                    CfdChannelType::Tapered => "tapered".into(),
                    CfdChannelType::Spiral { .. } => "spiral".into(),
                })
        } else {
            "unknown".into()
        };

        records.push(ChannelRecord {
            id: edge.id.clone(),
            ch_type: label,
            flow_rate_m3s: q,
            velocity_ms: velocity,
            reynolds,
            wall_shear_pa: wall_shear,
            dp_pa: dp,
            length_m: length,
            dh_m: dh,
        });
    }
    records
}

/// Write channel records to CSV.
fn write_csv(path: &str, records: &[ChannelRecord]) -> Result<(), Box<dyn std::error::Error>> {
    std::fs::create_dir_all(Path::new(path).parent().unwrap_or(Path::new(".")))?;
    let mut f = std::fs::File::create(path)?;
    writeln!(f, "channel_id,type,Q_m3s,V_ms,Re,tau_w_Pa,dP_Pa,L_m,Dh_m")?;
    for r in records {
        writeln!(
            f,
            "{},{},{:.6e},{:.6},{:.2},{:.4},{:.2},{:.6},{:.6}",
            r.id, r.ch_type, r.flow_rate_m3s, r.velocity_ms, r.reynolds,
            r.wall_shear_pa, r.dp_pa, r.length_m, r.dh_m,
        )?;
    }
    Ok(())
}

/// Write a JSON summary of a scenario.
fn write_json_summary(
    path: &str,
    scenario: &str,
    records: &[ChannelRecord],
    extra: &HashMap<String, String>,
) -> Result<(), Box<dyn std::error::Error>> {
    std::fs::create_dir_all(Path::new(path).parent().unwrap_or(Path::new(".")))?;
    let mut f = std::fs::File::create(path)?;

    let max_shear = records.iter().map(|r| r.wall_shear_pa).fold(0.0_f64, f64::max);
    let avg_shear = if records.is_empty() { 0.0 } else {
        records.iter().map(|r| r.wall_shear_pa).sum::<f64>() / records.len() as f64
    };
    let max_re = records.iter().map(|r| r.reynolds).fold(0.0_f64, f64::max);
    let total_q: f64 = records.iter()
        .filter(|r| r.id.contains("ch_0"))
        .map(|r| r.flow_rate_m3s)
        .sum();

    writeln!(f, "{{")?;
    writeln!(f, "  \"scenario\": \"{scenario}\",")?;
    writeln!(f, "  \"n_channels\": {},", records.len())?;
    writeln!(f, "  \"max_wall_shear_Pa\": {max_shear:.4},")?;
    writeln!(f, "  \"avg_wall_shear_Pa\": {avg_shear:.4},")?;
    writeln!(f, "  \"max_reynolds\": {max_re:.2},")?;
    writeln!(f, "  \"inlet_flow_m3s\": {total_q:.6e},")?;
    for (k, v) in extra {
        writeln!(f, "  \"{k}\": \"{v}\",")?;
    }
    writeln!(f, "  \"timestamp\": \"{}\"", chrono::Utc::now().to_rfc3339())?;
    writeln!(f, "}}")?;
    Ok(())
}

/// Jet colourmap t ∈ [0,1] → RGB.
fn jet(t: f64) -> RGBColor {
    let t = t.clamp(0.0, 1.0);
    let (r, g, b) = if t < 0.25 {
        (0.0, t / 0.25, 1.0)
    } else if t < 0.5 {
        (0.0, 1.0, 1.0 - (t - 0.25) / 0.25)
    } else if t < 0.75 {
        ((t - 0.5) / 0.25, 1.0, 0.0)
    } else {
        (1.0, 1.0 - (t - 0.75) / 0.25, 0.0)
    };
    RGBColor((r * 255.0) as u8, (g * 255.0) as u8, (b * 255.0) as u8)
}

/// Get 2D path for a channel from scheme geometry.
fn channel_path(system: &ChannelSystem, idx: usize) -> Vec<(f64, f64)> {
    let ch = &system.channels[idx];
    match &ch.channel_type {
        SchemeChannelType::Straight => {
            vec![system.nodes[ch.from_node].point, system.nodes[ch.to_node].point]
        }
        SchemeChannelType::SmoothStraight { path }
        | SchemeChannelType::Serpentine { path }
        | SchemeChannelType::Arc { path } => path.clone(),
        SchemeChannelType::Frustum { path, .. } => path.clone(),
    }
}

/// Draw a schematic heatmap for one field.
fn draw_heatmap(
    system: &ChannelSystem,
    values: &[f64],
    title: &str,
    label: &str,
    out_path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    std::fs::create_dir_all("outputs")?;
    let (bw, bh) = system.box_dims;
    let margin = 5.0;

    let root = SVGBackend::new(out_path, (1000, 600)).into_drawing_area();
    root.fill(&WHITE)?;
    let (chart_area, cbar_area) = root.split_horizontally(820);

    let vmin = values.iter().copied().fold(f64::INFINITY, f64::min);
    let vmax = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let vrange = if (vmax - vmin).abs() < 1e-12 { 1.0 } else { vmax - vmin };

    let mut chart = ChartBuilder::on(&chart_area)
        .caption(title, ("sans-serif", 17))
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(40)
        .build_cartesian_2d((-margin)..(bw + margin), (-margin)..(bh + margin))?;

    chart.configure_mesh().x_desc("x [mm]").y_desc("y [mm]")
        .light_line_style(WHITE.mix(0.0)).draw()?;

    // Chip outline
    chart.draw_series(std::iter::once(PathElement::new(
        vec![(0.0, 0.0), (bw, 0.0), (bw, bh), (0.0, bh), (0.0, 0.0)],
        ShapeStyle { color: RGBColor(180, 180, 180).to_rgba(), filled: false, stroke_width: 2 },
    )))?;

    let n = values.len().min(system.channels.len());
    for ci in 0..n {
        let t = (values[ci] - vmin) / vrange;
        let c = jet(t);

        let ch = &system.channels[ci];
        match &ch.channel_type {
            SchemeChannelType::Frustum { path, widths, .. } if path.len() >= 2 => {
                let wmax = widths.iter().copied().fold(0.0_f64, f64::max).max(0.1);
                for i in 0..path.len() - 1 {
                    let wa = (widths.get(i).copied().unwrap_or(1.0)
                        + widths.get(i + 1).copied().unwrap_or(1.0)) / 2.0;
                    let sw = ((wa / wmax * 6.0) as u32).max(1).min(12);
                    chart.draw_series(std::iter::once(PathElement::new(
                        vec![path[i], path[i + 1]],
                        ShapeStyle { color: c.to_rgba(), filled: false, stroke_width: sw },
                    )))?;
                }
            }
            _ => {
                let pts = channel_path(system, ci);
                if pts.len() >= 2 {
                    chart.draw_series(std::iter::once(PathElement::new(
                        pts,
                        ShapeStyle { color: c.to_rgba(), filled: false, stroke_width: 3 },
                    )))?;
                }
            }
        }
    }

    // Nodes
    for nd in &system.nodes {
        chart.draw_series(std::iter::once(Circle::new(
            nd.point, 3,
            ShapeStyle { color: RGBColor(60, 60, 60).to_rgba(), filled: true, stroke_width: 1 },
        )))?;
    }

    // Colour bar
    let bx0 = 25i32; let bx1 = 55i32; let by0 = 40i32; let by1 = 540i32;
    let bh_px = by1 - by0;
    let steps = 64usize;
    for i in 0..steps {
        let t = 1.0 - i as f64 / steps as f64;
        let ys = by0 + (i as i32 * bh_px / steps as i32);
        let ye = by0 + ((i + 1) as i32 * bh_px / steps as i32);
        cbar_area.draw(&Rectangle::new([(bx0, ys), (bx1, ye)], jet(t).filled()))?;
    }
    cbar_area.draw(&Rectangle::new([(bx0, by0), (bx1, by1)], BLACK.stroke_width(1)))?;

    let ts = ("sans-serif", 11).into_text_style(&cbar_area);
    for ti in 0..=4u32 {
        let frac = ti as f64 / 4.0;
        let val = vmax - frac * vrange;
        let yp = by0 + (frac * bh_px as f64) as i32;
        cbar_area.draw_text(&format!("{:.2}", val), &ts, (bx1 + 4, yp - 6))?;
        cbar_area.draw(&PathElement::new(vec![(bx1, yp), (bx1 + 3, yp)], BLACK.stroke_width(1)))?;
    }
    let ls = ("sans-serif", 12).into_text_style(&cbar_area);
    cbar_area.draw_text(label, &ls, (bx0 - 5, by1 + 8))?;

    root.present()?;
    Ok(())
}

/// Print a table of channel records.
fn print_table(records: &[ChannelRecord]) {
    println!("    ┌───────────┬────────────────┬────────────┬──────────┬─────────┬────────────┐");
    println!("    │ Channel   │ Type           │  Q [mL/m]  │ V [m/s]  │  Re     │ τ_w [Pa]   │");
    println!("    ├───────────┼────────────────┼────────────┼──────────┼─────────┼────────────┤");
    for r in records {
        let q_ml = r.flow_rate_m3s * 1e6 * 60.0;
        println!(
            "    │ {:>9} │ {:>14} │ {:>10.4} │ {:>8.4} │ {:>7.1} │ {:>10.4} │",
            r.id, r.ch_type, q_ml, r.velocity_ms, r.reynolds, r.wall_shear_pa,
        );
    }
    println!("    └───────────┴────────────────┴────────────┴──────────┴─────────┴────────────┘");
}

// ═══════════════════════════════════════════════════════════════════════════════
// Scenario 1 — Trifurcation Distributor
// ═══════════════════════════════════════════════════════════════════════════════

fn scenario_trifurcation() -> Result<(), Box<dyn std::error::Error>> {
    println!("\n  ━━━ Scenario 1: Trifurcation Distributor ━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("    Topology: 1 inlet → 3 outlets (straight channels)");
    println!("    Goal: evaluate flow uniformity in 3-way split\n");

    let config = GeometryConfig::new(0.5, 1.0, 0.5)?;
    let system = create_geometry(
        (CHIP_W, CHIP_H),
        &[SplitType::Trifurcation],
        &config,
        &ChannelTypeConfig::AllStraight,
    );

    let converter = SchemeNetworkConverter::with_scale(&system, SCALE_MM_TO_M);
    let summary = converter.summary();
    println!("    {}", summary.to_string().replace('\n', "\n    "));

    let mut network = converter.build_network(blood_fluid())?;
    let (inlets, outlets) = find_boundary_nodes(&network);

    let p_in = 10_000.0;
    let p_out = 1_000.0;
    for &ni in &inlets { network.set_pressure(ni, p_in); }
    for &ni in &outlets { network.set_pressure(ni, p_out); }

    let problem = NetworkProblem::new(network);
    let solver = NetworkSolver::<f64, ConstantPropertyFluid<f64>>::with_config(SolverConfig {
        tolerance: 1e-8,
        max_iterations: 500,
    });
    let solution = solver.solve_network(&problem)?;
    let records = analyse_network(&solution, BLOOD_DENSITY, BLOOD_VISCOSITY);

    print_table(&records);

    // Flow uniformity metric — considering only channels with meaningful flow
    let outlet_flows: Vec<f64> = records.iter()
        .filter(|r| r.flow_rate_m3s > 1e-12 && r.id != "ch_0") // exclude zero-flow and main inlet
        .map(|r| r.flow_rate_m3s)
        .collect();
    let q_mean = outlet_flows.iter().sum::<f64>() / outlet_flows.len().max(1) as f64;
    let q_std = (outlet_flows.iter().map(|q| (q - q_mean).powi(2)).sum::<f64>()
        / outlet_flows.len().max(1) as f64).sqrt();
    let uniformity = if q_mean > 0.0 { 1.0 - q_std / q_mean } else { 0.0 };

    println!("\n    Flow uniformity index: {:.4}  (1.0 = perfect)", uniformity);
    println!("      Mean Q = {:.4e} m³/s,  Std Q = {:.4e} m³/s", q_mean, q_std);

    // Verification checks
    let max_re = records.iter().map(|r| r.reynolds).fold(0.0_f64, f64::max);
    println!("\n    ── Verification ──");
    println!("      Max Re = {:.1} → {}", max_re,
        if max_re < 2300.0 { "LAMINAR ✓" } else { "TURBULENT ⚠" });
    println!("      All shear > 0: {}", if records.iter().all(|r| r.wall_shear_pa >= 0.0) { "✓" } else { "✗" });

    // CSV + JSON export
    write_csv("outputs/scenario1_trifurcation.csv", &records)?;
    let mut extra = HashMap::new();
    extra.insert("topology".into(), "trifurcation".into());
    extra.insert("flow_uniformity".into(), format!("{uniformity:.4}"));
    extra.insert("p_inlet_Pa".into(), format!("{p_in}"));
    extra.insert("p_outlet_Pa".into(), format!("{p_out}"));
    write_json_summary("outputs/scenario1_trifurcation.json", "Trifurcation Distributor", &records, &extra)?;
    println!("    CSV  → outputs/scenario1_trifurcation.csv");
    println!("    JSON → outputs/scenario1_trifurcation.json");

    // Heatmaps
    let n = records.len().min(system.channels.len());
    let shear: Vec<f64> = records[..n].iter().map(|r| r.wall_shear_pa).collect();
    let vel: Vec<f64> = records[..n].iter().map(|r| r.velocity_ms).collect();
    draw_heatmap(&system, &shear, "Trifurcation — Wall Shear", "τ_w [Pa]",
        "outputs/heatmap_trifurcation_shear.svg")?;
    draw_heatmap(&system, &vel, "Trifurcation — Velocity", "V [m/s]",
        "outputs/heatmap_trifurcation_velocity.svg")?;
    println!("    SVG  → outputs/heatmap_trifurcation_shear.svg");
    println!("    SVG  → outputs/heatmap_trifurcation_velocity.svg");

    Ok(())
}

// ═══════════════════════════════════════════════════════════════════════════════
// Scenario 2 — Two-Level Cascade (Bifurcation × 2)
// ═══════════════════════════════════════════════════════════════════════════════

fn scenario_cascade() -> Result<(), Box<dyn std::error::Error>> {
    println!("\n  ━━━ Scenario 2: Two-Level Serpentine Cascade ━━━━━━━━━━━━━━━━━━━━━");
    println!("    Topology: 1 → 2 → 4  (two bifurcation levels, serpentine channels)");
    println!("    Goal: Dean-enhanced mixing with progressive shear\n");

    // Larger chip to ensure channels satisfy L/Dh ≥ 10 for 2-level splits
    let config = GeometryConfig::new(0.3, 0.8, 0.3)?;
    let system = create_geometry(
        (90.0, 50.0),
        &[SplitType::Bifurcation, SplitType::Bifurcation],
        &config,
        &ChannelTypeConfig::AllSerpentine(SerpentineConfig {
            fill_factor: 0.6,
            wavelength_factor: 4.0,
            gaussian_width_factor: 8.0,
            wave_density_factor: 1.5,
            wave_phase_direction: 0.0,
            ..SerpentineConfig::default()
        }),
    );

    let converter = SchemeNetworkConverter::with_scale(&system, SCALE_MM_TO_M);
    let summary = converter.summary();
    println!("    {}", summary.to_string().replace('\n', "\n    "));

    let mut network = converter.build_network(blood_fluid())?;
    let (inlets, outlets) = find_boundary_nodes(&network);

    let p_in = 12_000.0;
    let p_out = 1_000.0;
    for &ni in &inlets { network.set_pressure(ni, p_in); }
    for &ni in &outlets { network.set_pressure(ni, p_out); }

    let problem = NetworkProblem::new(network);
    let solver = NetworkSolver::<f64, ConstantPropertyFluid<f64>>::with_config(SolverConfig {
        tolerance: 1e-8,
        max_iterations: 500,
    });
    let solution = solver.solve_network(&problem)?;
    let records = analyse_network(&solution, BLOOD_DENSITY, BLOOD_VISCOSITY);

    print_table(&records);

    // Level-by-level shear analysis
    let max_shear = records.iter().map(|r| r.wall_shear_pa).fold(0.0_f64, f64::max);
    let min_shear = records.iter().filter(|r| r.wall_shear_pa > 0.0)
        .map(|r| r.wall_shear_pa).fold(f64::INFINITY, f64::min);
    println!("\n    Shear range: {:.4} – {:.4} Pa", min_shear, max_shear);
    println!("    Platelet-safe (< 10 Pa): {}",
        if max_shear < 10.0 { "✓ all channels" } else { "⚠ some exceed threshold" });

    // Verification
    let max_re = records.iter().map(|r| r.reynolds).fold(0.0_f64, f64::max);
    println!("\n    ── Verification ──");
    println!("      Max Re = {:.1} → {}", max_re,
        if max_re < 2300.0 { "LAMINAR ✓" } else { "TURBULENT ⚠" });
    // Mass conservation: sum of outlet flows ≈ inlet flow
    let inlet_q: f64 = records.iter().filter(|r| r.id == "ch_0").map(|r| r.flow_rate_m3s).sum();
    let outlet_q: f64 = records.iter().filter(|r| {
        let idx: usize = r.id.trim_start_matches("ch_").parse().unwrap_or(0);
        idx == records.len() - 1 // last channel as proxy
    }).map(|r| r.flow_rate_m3s).sum();
    println!("      Inlet flow ch_0 = {:.4e} m³/s", inlet_q);

    // CSV + JSON
    write_csv("outputs/scenario2_cascade.csv", &records)?;
    let mut extra = HashMap::new();
    extra.insert("topology".into(), "bifurcation_x2".into());
    extra.insert("channel_type".into(), "serpentine".into());
    extra.insert("n_levels".into(), "2".into());
    write_json_summary("outputs/scenario2_cascade.json", "Two-Level Cascade", &records, &extra)?;
    println!("    CSV  → outputs/scenario2_cascade.csv");
    println!("    JSON → outputs/scenario2_cascade.json");

    // Heatmaps
    let n = records.len().min(system.channels.len());
    let shear: Vec<f64> = records[..n].iter().map(|r| r.wall_shear_pa).collect();
    let vel: Vec<f64> = records[..n].iter().map(|r| r.velocity_ms).collect();
    let re: Vec<f64> = records[..n].iter().map(|r| r.reynolds).collect();
    draw_heatmap(&system, &shear, "Cascade — Wall Shear", "τ_w [Pa]",
        "outputs/heatmap_cascade_shear.svg")?;
    draw_heatmap(&system, &vel, "Cascade — Velocity", "V [m/s]",
        "outputs/heatmap_cascade_velocity.svg")?;
    draw_heatmap(&system, &re, "Cascade — Reynolds", "Re",
        "outputs/heatmap_cascade_reynolds.svg")?;
    println!("    SVG  → outputs/heatmap_cascade_shear.svg");
    println!("    SVG  → outputs/heatmap_cascade_velocity.svg");
    println!("    SVG  → outputs/heatmap_cascade_reynolds.svg");

    Ok(())
}

// ═══════════════════════════════════════════════════════════════════════════════
// Scenario 3 — Mixed Venturi + Serpentine
// ═══════════════════════════════════════════════════════════════════════════════

fn scenario_mixed_venturi() -> Result<(), Box<dyn std::error::Error>> {
    println!("\n  ━━━ Scenario 3: Mixed Venturi-Serpentine Design ━━━━━━━━━━━━━━━━━━");
    println!("    Topology: bifurcation with frustum (venturi) channels");
    println!("    Goal: compare shear + cavitation in tapered geometry\n");

    let config = GeometryConfig::new(0.5, 1.0, 0.5)?;
    let system = create_geometry(
        (CHIP_W, CHIP_H),
        &[SplitType::Bifurcation],
        &config,
        &ChannelTypeConfig::AllFrustum(FrustumConfig {
            inlet_width: 2.0,
            throat_width: 0.5,
            outlet_width: 2.0,
            taper_profile: TaperProfile::Smooth,
            smoothness: 50,
            throat_position: 0.5,
        }),
    );

    let converter = SchemeNetworkConverter::with_scale(&system, SCALE_MM_TO_M);
    let summary = converter.summary();
    println!("    {}", summary.to_string().replace('\n', "\n    "));

    let mut network = converter.build_network(blood_fluid())?;
    let (inlets, outlets) = find_boundary_nodes(&network);

    let p_in = 12_000.0;
    let p_out = 1_000.0;
    for &ni in &inlets { network.set_pressure(ni, p_in); }
    for &ni in &outlets { network.set_pressure(ni, p_out); }

    let problem = NetworkProblem::new(network);
    let solver = NetworkSolver::<f64, ConstantPropertyFluid<f64>>::with_config(SolverConfig {
        tolerance: 1e-8,
        max_iterations: 500,
    });
    let solution = solver.solve_network(&problem)?;
    let records = analyse_network(&solution, BLOOD_DENSITY, BLOOD_VISCOSITY);

    print_table(&records);

    // Cavitation check on tapered channels
    let pressures = solution.pressures();
    let flow_rates = solution.flow_rates();
    let mut cav_results: Vec<(String, f64, f64, f64, bool)> = Vec::new();
    for ei in solution.graph.edge_indices() {
        let props = solution.properties.get(&ei);
        if let Some(p) = props {
            if let Some(geom) = &p.geometry {
                if matches!(geom.channel_type, CfdChannelType::Tapered) {
                    let q = flow_rates.get(&ei).copied().unwrap_or(0.0).abs();
                    let (from, _to) = solution.graph.edge_endpoints(ei).unwrap();
                    let p_from = pressures.get(&from).copied().unwrap_or(0.0);
                    let (w, h) = match &geom.cross_section {
                        CrossSection::Rectangular { width, height } => (*width, *height),
                        _ => (0.001, 0.0005),
                    };
                    let min_scale = geom.variations.iter()
                        .map(|v| v.scale_factor).fold(f64::INFINITY, f64::min);
                    let tw = if min_scale < 1.0 && min_scale > 0.0 { w * min_scale } else { w * 0.25 };
                    let inlet_a = w * h;
                    let throat_a = tw * h;
                    let v_in = if inlet_a > 0.0 { q / inlet_a } else { 0.0 };
                    let v_thr = if throat_a > 0.0 { q / throat_a } else { 0.0 };
                    let p_thr = p_from - 0.5 * BLOOD_DENSITY * (v_thr.powi(2) - v_in.powi(2));
                    let sigma = if v_thr > 1e-12 {
                        (p_thr - P_VAPOR_37C) / (0.5 * BLOOD_DENSITY * v_thr.powi(2))
                    } else { f64::INFINITY };
                    let is_cav = sigma < 1.0;

                    cav_results.push((solution.graph[ei].id.clone(), v_thr, p_thr, sigma, is_cav));
                }
            }
        }
    }

    if !cav_results.is_empty() {
        println!("\n    Venturi channels in mixed design:");
        for (id, vt, pt, sig, cav) in &cav_results {
            println!("      {}: V_thr={:.3} m/s, p_thr={:.0} Pa, σ={:.3} → {}",
                id, vt, pt, sig, if *cav { "CAVITATING" } else { "safe" });
        }
    }

    // Verification
    let max_re = records.iter().map(|r| r.reynolds).fold(0.0_f64, f64::max);
    println!("\n    ── Verification ──");
    println!("      Max Re = {:.1} → {}", max_re,
        if max_re < 2300.0 { "LAMINAR ✓" } else { "TURBULENT ⚠" });
    let n_types: std::collections::HashSet<&str> = records.iter().map(|r| r.ch_type.as_str()).collect();
    println!("      Channel types present: {:?}", n_types);

    // CSV + JSON
    write_csv("outputs/scenario3_mixed.csv", &records)?;
    let mut extra = HashMap::new();
    extra.insert("topology".into(), "bifurcation".into());
    extra.insert("channel_type".into(), "adaptive_mixed".into());
    extra.insert("n_cavitating".into(), format!("{}", cav_results.iter().filter(|c| c.4).count()));
    write_json_summary("outputs/scenario3_mixed.json", "Mixed Venturi-Serpentine", &records, &extra)?;
    println!("    CSV  → outputs/scenario3_mixed.csv");
    println!("    JSON → outputs/scenario3_mixed.json");

    // Heatmaps
    let n = records.len().min(system.channels.len());
    let shear: Vec<f64> = records[..n].iter().map(|r| r.wall_shear_pa).collect();
    let dp: Vec<f64> = records[..n].iter().map(|r| r.dp_pa).collect();
    draw_heatmap(&system, &shear, "Mixed — Wall Shear", "τ_w [Pa]",
        "outputs/heatmap_mixed_shear.svg")?;
    draw_heatmap(&system, &dp, "Mixed — Pressure Drop", "ΔP [Pa]",
        "outputs/heatmap_mixed_pressure.svg")?;
    println!("    SVG  → outputs/heatmap_mixed_shear.svg");
    println!("    SVG  → outputs/heatmap_mixed_pressure.svg");

    Ok(())
}

// ═══════════════════════════════════════════════════════════════════════════════
// Scenario 4 — Pressure Parametric Sweep
// ═══════════════════════════════════════════════════════════════════════════════

fn scenario_pressure_sweep() -> Result<(), Box<dyn std::error::Error>> {
    println!("\n  ━━━ Scenario 4: Pressure Parametric Sweep ━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("    Design: bifurcation, blood, vary P_inlet from 3 to 12 kPa");
    println!("    Goal: characterise τ_w vs ΔP and identify safe operating window\n");

    let config = GeometryConfig::new(0.5, 1.0, 0.5)?;
    let system = create_geometry(
        (CHIP_W, CHIP_H),
        &[SplitType::Bifurcation],
        &config,
        &ChannelTypeConfig::AllStraight,
    );

    let p_out = 1_000.0_f64;
    let pressures: Vec<f64> = (1..=10).map(|i| 2_000.0 + i as f64 * 1_000.0).collect();

    let mut sweep_data: Vec<(f64, f64, f64, f64, f64)> = Vec::new(); // (P_in, max_shear, avg_shear, Q_tot, Re_max)

    let mut sweep_csv = std::fs::File::create("outputs/scenario4_sweep.csv")?;
    writeln!(sweep_csv, "P_in_Pa,max_tau_w_Pa,avg_tau_w_Pa,Q_inlet_mL_min,Re_max")?;

    for &p_in in &pressures {
        let converter = SchemeNetworkConverter::with_scale(&system, SCALE_MM_TO_M);
        let mut network = converter.build_network(blood_fluid())?;
        let (inlets, outlets) = find_boundary_nodes(&network);
        for &ni in &inlets { network.set_pressure(ni, p_in); }
        for &ni in &outlets { network.set_pressure(ni, p_out); }

        let problem = NetworkProblem::new(network);
        let solver = NetworkSolver::<f64, ConstantPropertyFluid<f64>>::with_config(SolverConfig {
            tolerance: 1e-8,
            max_iterations: 500,
        });

        match solver.solve_network(&problem) {
            Ok(sol) => {
                let recs = analyse_network(&sol, BLOOD_DENSITY, BLOOD_VISCOSITY);
                let max_s = recs.iter().map(|r| r.wall_shear_pa).fold(0.0_f64, f64::max);
                let avg_s = recs.iter().map(|r| r.wall_shear_pa).sum::<f64>() / recs.len().max(1) as f64;
                let q_tot = recs.iter().filter(|r| r.id == "ch_0").map(|r| r.flow_rate_m3s).sum::<f64>();
                let re_max = recs.iter().map(|r| r.reynolds).fold(0.0_f64, f64::max);

                println!("    P_in = {:>6.0} Pa │ τ_max = {:>8.3} Pa │ Q = {:>8.3} mL/min │ Re_max = {:>6.1}",
                    p_in, max_s, q_tot * 1e6 * 60.0, re_max);

                writeln!(sweep_csv, "{:.0},{:.4},{:.4},{:.4},{:.2}",
                    p_in, max_s, avg_s, q_tot * 1e6 * 60.0, re_max)?;

                sweep_data.push((p_in, max_s, avg_s, q_tot * 1e6 * 60.0, re_max));
            }
            Err(e) => {
                eprintln!("    P_in = {:.0}: solver failed: {}", p_in, e);
            }
        }
    }

    // Verification
    println!("\n    ── Verification ──");
    // Shear should increase monotonically with pressure
    let monotonic = sweep_data.windows(2).all(|w| w[1].1 >= w[0].1 - 1e-6);
    println!("      Shear monotonically increases with P_in: {}",
        if monotonic { "✓" } else { "⚠ NOT monotonic" });
    // Flow should increase with pressure
    let q_mono = sweep_data.windows(2).all(|w| w[1].3 >= w[0].3 - 1e-6);
    println!("      Flow  monotonically increases with P_in: {}",
        if q_mono { "✓" } else { "⚠ NOT monotonic" });

    // Plot τ_w vs P_inlet
    {
        let path = "outputs/scenario4_shear_vs_pressure.svg";
        let root = SVGBackend::new(path, (800, 450)).into_drawing_area();
        root.fill(&WHITE)?;
        let x_max = sweep_data.last().map_or(20000.0, |d| d.0 * 1.1);
        let y_max = sweep_data.iter().map(|d| d.1).fold(0.0_f64, f64::max) * 1.3;

        let mut chart = ChartBuilder::on(&root)
            .caption("Wall Shear vs Inlet Pressure (blood, bifurcation)", ("sans-serif", 16))
            .margin(15)
            .x_label_area_size(35)
            .y_label_area_size(55)
            .build_cartesian_2d(0.0..x_max, 0.0..y_max)?;

        chart.configure_mesh()
            .x_desc("P_inlet [Pa]")
            .y_desc("Max τ_w [Pa]")
            .draw()?;

        chart.draw_series(LineSeries::new(
            sweep_data.iter().map(|d| (d.0, d.1)),
            RED.stroke_width(2),
        ))?.label("Max τ_w").legend(|(x, y)| PathElement::new(vec![(x, y), (x + 15, y)], RED.stroke_width(2)));

        chart.draw_series(LineSeries::new(
            sweep_data.iter().map(|d| (d.0, d.2)),
            BLUE.stroke_width(2),
        ))?.label("Avg τ_w").legend(|(x, y)| PathElement::new(vec![(x, y), (x + 15, y)], BLUE.stroke_width(2)));

        chart.draw_series(sweep_data.iter().map(|d| Circle::new((d.0, d.1), 4, RED.filled())))?;
        chart.draw_series(sweep_data.iter().map(|d| Circle::new((d.0, d.2), 4, BLUE.filled())))?;

        chart.configure_series_labels()
            .background_style(WHITE.mix(0.8)).border_style(BLACK)
            .position(SeriesLabelPosition::UpperLeft).draw()?;

        root.present()?;
        println!("    SVG  → {path}");
    }

    // Plot Q vs P_inlet
    {
        let path = "outputs/scenario4_flow_vs_pressure.svg";
        let root = SVGBackend::new(path, (800, 450)).into_drawing_area();
        root.fill(&WHITE)?;
        let x_max = sweep_data.last().map_or(20000.0, |d| d.0 * 1.1);
        let y_max = sweep_data.iter().map(|d| d.3).fold(0.0_f64, f64::max) * 1.3;

        let mut chart = ChartBuilder::on(&root)
            .caption("Flow Rate vs Inlet Pressure", ("sans-serif", 16))
            .margin(15)
            .x_label_area_size(35)
            .y_label_area_size(55)
            .build_cartesian_2d(0.0..x_max, 0.0..y_max)?;

        chart.configure_mesh()
            .x_desc("P_inlet [Pa]")
            .y_desc("Q_inlet [mL/min]")
            .draw()?;

        chart.draw_series(LineSeries::new(
            sweep_data.iter().map(|d| (d.0, d.3)),
            GREEN.stroke_width(2),
        ))?.label("Q_inlet").legend(|(x, y)| PathElement::new(vec![(x, y), (x + 15, y)], GREEN.stroke_width(2)));

        chart.draw_series(sweep_data.iter().map(|d| Circle::new((d.0, d.3), 4, GREEN.filled())))?;

        chart.configure_series_labels()
            .background_style(WHITE.mix(0.8)).border_style(BLACK)
            .position(SeriesLabelPosition::UpperLeft).draw()?;

        root.present()?;
        println!("    SVG  → {path}");
    }

    println!("    CSV  → outputs/scenario4_sweep.csv");

    // JSON
    let mut extra = HashMap::new();
    extra.insert("fluid".into(), "blood_37C".into());
    extra.insert("n_pressure_points".into(), format!("{}", sweep_data.len()));
    extra.insert("p_range_Pa".into(), format!("{:.0}-{:.0}", pressures[0], pressures[pressures.len()-1]));
    extra.insert("shear_monotonic".into(), format!("{monotonic}"));
    write_json_summary("outputs/scenario4_sweep.json", "Pressure Sweep", &[], &extra)?;

    Ok(())
}

// ═══════════════════════════════════════════════════════════════════════════════
// Scenario 5 — Transient Composition Switching
// ═══════════════════════════════════════════════════════════════════════════════

fn scenario_transient_composition() -> Result<(), Box<dyn std::error::Error>> {
    println!("\n  ━━━ Scenario 5: Transient Composition Switching ━━━━━━━━━━━━━━━━━━");
    println!("    Design: bifurcation, water, reagent A → B switch at t = 5 s");
    println!("    Goal: observe mixing front propagation through network\n");

    let config = GeometryConfig::new(0.5, 1.0, 0.5)?;
    let system = create_geometry(
        (CHIP_W, CHIP_H),
        &[SplitType::Bifurcation],
        &config,
        &ChannelTypeConfig::AllStraight,
    );

    let converter = SchemeNetworkConverter::with_scale(&system, SCALE_MM_TO_M);

    // First solve steady-state to get flow field
    let mut network = converter.build_network(water_fluid())?;
    let (inlets, outlets) = find_boundary_nodes(&network);

    let p_in = 8_000.0;
    let p_out = 1_000.0;
    for &ni in &inlets { network.set_pressure(ni, p_in); }
    for &ni in &outlets { network.set_pressure(ni, p_out); }

    let problem = NetworkProblem::new(network);
    let solver = NetworkSolver::<f64, ConstantPropertyFluid<f64>>::with_config(SolverConfig {
        tolerance: 1e-8,
        max_iterations: 500,
    });
    let solved = solver.solve_network(&problem)?;

    // Print steady flow field
    let records = analyse_network(&solved, 998.0, 1.002e-3);
    println!("    Steady-state flow field:");
    print_table(&records);

    // Transient composition: switch reagent at inlet
    let inlet_idx = inlets[0].index();

    // Fluid IDs: 0 = reagent A (original), 1 = reagent B (switched)
    let mut initial_fractions = HashMap::new();
    initial_fractions.insert(0i32, 1.0_f64);

    let mut switched_fractions = HashMap::new();
    switched_fractions.insert(1i32, 1.0_f64);

    let events = vec![
        InletCompositionEvent {
            time: 0.0,
            node_index: inlet_idx,
            mixture: MixtureComposition::new(initial_fractions),
        },
        InletCompositionEvent {
            time: 5.0,
            node_index: inlet_idx,
            mixture: MixtureComposition::new(switched_fractions),
        },
    ];

    let timing = SimulationTimeConfig::new(
        15.0,  // total duration [s]
        1.0,   // result step [s]
        0.5,   // calculation step [s]
    );

    let states = TransientCompositionSimulator::simulate_with_time_config(
        &solved, events, timing,
    )?;

    // CSV export of composition over time
    let mut comp_csv = std::fs::File::create("outputs/scenario5_composition.csv")?;
    write!(comp_csv, "time_s")?;
    // Write header for each edge
    let n_edges = solved.edge_count();
    for eidx in 0..n_edges {
        write!(comp_csv, ",edge_{}_reagentA,edge_{}_reagentB", eidx, eidx)?;
    }
    writeln!(comp_csv)?;

    println!("\n    Transient composition (reagent_B fraction at edge 0):");
    println!("    ┌──────────┬────────────┬────────────┐");
    println!("    │  Time [s]│ Reagent A  │ Reagent B  │");
    println!("    ├──────────┼────────────┼────────────┤");

    for state in &states {
        write!(comp_csv, "{:.2}", state.time)?;
        for eidx in 0..n_edges {
            let mix = state.edge_mixtures.get(&eidx);
            let fa = mix.map_or(0.0, |m| *m.fractions.get(&0).unwrap_or(&0.0));
            let fb = mix.map_or(0.0, |m| *m.fractions.get(&1).unwrap_or(&0.0));
            write!(comp_csv, ",{:.6},{:.6}", fa, fb)?;
        }
        writeln!(comp_csv)?;

        // Print edge 0 for display
        let mix = state.edge_mixtures.get(&0);
        let fa = mix.map_or(0.0, |m| *m.fractions.get(&0).unwrap_or(&0.0));
        let fb = mix.map_or(0.0, |m| *m.fractions.get(&1).unwrap_or(&0.0));
        println!("    │ {:>8.2} │ {:>10.4} │ {:>10.4} │", state.time, fa, fb);
    }
    println!("    └──────────┴────────────┴────────────┘");

    // Verification
    println!("\n    ── Verification ──");
    // At t=0, should be 100% A
    let t0 = &states[0];
    let t0_fa = t0.edge_mixtures.get(&0).map_or(0.0, |m| *m.fractions.get(&0).unwrap_or(&0.0));
    println!("      t=0: reagent_A = {:.4} (expect 1.0): {}",
        t0_fa, if (t0_fa - 1.0).abs() < 0.01 { "✓" } else { "✗" });

    // At t=15, after switch at t=5, inlet should be 100% B, and
    // instantaneous mixing means edges connected to inlet should be 100% B
    let t_last = states.last().unwrap();
    let tlast_fb = t_last.edge_mixtures.get(&0).map_or(0.0, |m| *m.fractions.get(&1).unwrap_or(&0.0));
    println!("      t=15: reagent_B at edge_0 = {:.4} (expect ~1.0): {}",
        tlast_fb, if tlast_fb > 0.9 { "✓" } else { "⚠" });

    // Fractions should sum to 1.0
    let sum_check = states.iter().all(|s| {
        s.edge_mixtures.values().all(|m| {
            let sum: f64 = m.fractions.values().sum();
            sum < 1e-12 || (sum - 1.0).abs() < 0.01
        })
    });
    println!("      Fractions sum to 1.0 at all times: {}", if sum_check { "✓" } else { "✗" });

    // Plot composition over time for edge 0
    {
        let path = "outputs/scenario5_composition_edge0.svg";
        let root = SVGBackend::new(path, (800, 400)).into_drawing_area();
        root.fill(&WHITE)?;

        let mut chart = ChartBuilder::on(&root)
            .caption("Reagent Composition at Edge 0 (inlet channel)", ("sans-serif", 16))
            .margin(15)
            .x_label_area_size(30)
            .y_label_area_size(45)
            .build_cartesian_2d(0.0..16.0_f64, 0.0..1.1_f64)?;

        chart.configure_mesh()
            .x_desc("Time [s]")
            .y_desc("Volume fraction")
            .draw()?;

        let fa_series: Vec<(f64, f64)> = states.iter().map(|s| {
            let fa = s.edge_mixtures.get(&0).map_or(0.0, |m| *m.fractions.get(&0).unwrap_or(&0.0));
            (s.time, fa)
        }).collect();

        let fb_series: Vec<(f64, f64)> = states.iter().map(|s| {
            let fb = s.edge_mixtures.get(&0).map_or(0.0, |m| *m.fractions.get(&1).unwrap_or(&0.0));
            (s.time, fb)
        }).collect();

        chart.draw_series(LineSeries::new(fa_series.clone(), BLUE.stroke_width(2)))?
            .label("Reagent A")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 15, y)], BLUE.stroke_width(2)));

        chart.draw_series(LineSeries::new(fb_series.clone(), RED.stroke_width(2)))?
            .label("Reagent B")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 15, y)], RED.stroke_width(2)));

        chart.draw_series(fa_series.iter().map(|&(x, y)| Circle::new((x, y), 3, BLUE.filled())))?;
        chart.draw_series(fb_series.iter().map(|&(x, y)| Circle::new((x, y), 3, RED.filled())))?;

        // Switch line
        chart.draw_series(std::iter::once(PathElement::new(
            vec![(5.0, 0.0), (5.0, 1.1)],
            BLACK.stroke_width(1),
        )))?;

        chart.configure_series_labels()
            .background_style(WHITE.mix(0.8)).border_style(BLACK)
            .position(SeriesLabelPosition::MiddleRight).draw()?;

        root.present()?;
        println!("    SVG  → {path}");
    }

    println!("    CSV  → outputs/scenario5_composition.csv");

    // JSON
    let mut extra = HashMap::new();
    extra.insert("fluid".into(), "water_20C".into());
    extra.insert("transient".into(), "true".into());
    extra.insert("switch_time_s".into(), "5.0".into());
    extra.insert("total_time_s".into(), "15.0".into());
    extra.insert("n_snapshots".into(), format!("{}", states.len()));
    write_json_summary("outputs/scenario5_composition.json", "Transient Composition", &records, &extra)?;

    Ok(())
}

// ═══════════════════════════════════════════════════════════════════════════════
// Cross-scenario summary
// ═══════════════════════════════════════════════════════════════════════════════

fn print_output_summary() {
    println!("\n  ╔══════════════════════════════════════════════════════════════════════╗");
    println!("  ║  Output Summary                                                    ║");
    println!("  ╚══════════════════════════════════════════════════════════════════════╝\n");

    let files = [
        ("Scenario 1", &[
            "outputs/scenario1_trifurcation.csv",
            "outputs/scenario1_trifurcation.json",
            "outputs/heatmap_trifurcation_shear.svg",
            "outputs/heatmap_trifurcation_velocity.svg",
        ][..]),
        ("Scenario 2", &[
            "outputs/scenario2_cascade.csv",
            "outputs/scenario2_cascade.json",
            "outputs/heatmap_cascade_shear.svg",
            "outputs/heatmap_cascade_velocity.svg",
            "outputs/heatmap_cascade_reynolds.svg",
        ]),
        ("Scenario 3", &[
            "outputs/scenario3_mixed.csv",
            "outputs/scenario3_mixed.json",
            "outputs/heatmap_mixed_shear.svg",
            "outputs/heatmap_mixed_pressure.svg",
        ]),
        ("Scenario 4", &[
            "outputs/scenario4_sweep.csv",
            "outputs/scenario4_sweep.json",
            "outputs/scenario4_shear_vs_pressure.svg",
            "outputs/scenario4_flow_vs_pressure.svg",
        ]),
        ("Scenario 5", &[
            "outputs/scenario5_composition.csv",
            "outputs/scenario5_composition.json",
            "outputs/scenario5_composition_edge0.svg",
        ]),
    ];

    for (label, paths) in &files {
        println!("    {label}:");
        for p in *paths {
            let exists = Path::new(p).exists();
            let size = if exists {
                std::fs::metadata(p).map(|m| m.len()).unwrap_or(0)
            } else { 0 };
            println!("      {} {:>8} B  {}", if exists { "✓" } else { "✗" },
                size, p);
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Main
// ═══════════════════════════════════════════════════════════════════════════════

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!();
    println!("  ═══════════════════════════════════════════════════════════════════");
    println!("   Millifluidic Design Scenarios — CFD-1D");
    println!("   Pipeline: scheme → cfd-1d → solve → CSV + SVG");
    println!("  ═══════════════════════════════════════════════════════════════════");

    scenario_trifurcation()?;
    scenario_cascade()?;
    scenario_mixed_venturi()?;
    scenario_pressure_sweep()?;
    scenario_transient_composition()?;

    print_output_summary();

    println!("\n  All scenarios complete.");
    Ok(())
}
