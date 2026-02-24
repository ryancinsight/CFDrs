//! # Blood Flow Through Millifluidic Devices
//!
//! Comprehensive simulation demonstrating the full pipeline:
//!
//! ```text
//! scheme design  →  cfd-1d conversion  →  simulation  →  shear & cavitation analysis  →  SVG plots
//! ```
//!
//! ## Scenario
//!
//! Blood is drawn from a brachial artery (arm without AV fistula) and fed into
//! three millifluidic chip designs sharing the same bifurcation topology but
//! different channel types:
//!
//! 1. **Straight** channels — baseline reference
//! 2. **Serpentine** channels — enhanced mixing for cell sorting / lysis
//! 3. **Frustum / Venturi** channels — converging–diverging throats for
//!    hydrodynamic cavitation analysis
//!
//! ## Blood Parameters (arm without AV fistula)
//!
//! | Parameter | Value | Source |
//! |-----------|-------|--------|
//! | Density ρ | 1 060 kg/m³ | Cho & Kensey (1991) |
//! | Dynamic viscosity μ | 3.5 × 10⁻³ Pa·s | High-shear Newtonian approx. |
//! | Specific heat cₚ | 3 600 J/(kg·K) | Duck (1990) |
//! | Thermal conductivity k | 0.50 W/(m·K) | Duck (1990) |
//! | Speed of sound c | 1 540 m/s | Wells (1977) |
//! | Mean arterial pressure | ≈ 93 mmHg ≈ 12 400 Pa | Normal physiology |
//! | Vapour pressure (37 °C) | 6 274 Pa | NIST |
//!
//! ## Running
//!
//! ```sh
//! cargo run --example blood_millifluidic_simulation --features scheme-integration
//! ```

// ═══════════════════════════════════════════════════════════════════════════════
// Imports
// ═══════════════════════════════════════════════════════════════════════════════

use std::collections::HashMap;

// cfd-1d
use cfd_1d::channel::{ChannelType as CfdChannelType, CrossSection};
use cfd_1d::network::{network_from_blueprint, Network, NodeType};
use cfd_1d::solver::{NetworkProblem, NetworkSolver, SolverConfig};

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
// Physical constants — blood from arm without AV fistula
// ═══════════════════════════════════════════════════════════════════════════════

/// Blood density [kg/m³] — Cho & Kensey (1991)
const BLOOD_DENSITY: f64 = 1_060.0;
/// Blood dynamic viscosity [Pa·s] — Newtonian approximation at arterial shear
/// rates (~100–1 000 s⁻¹); infinite-shear plateau of Carreau-Yasuda model.
const BLOOD_VISCOSITY: f64 = 3.5e-3;
/// Blood specific heat [J/(kg·K)] — Duck (1990)
const BLOOD_CP: f64 = 3_600.0;
/// Blood thermal conductivity [W/(m·K)] — Duck (1990)
const BLOOD_K: f64 = 0.50;
/// Blood speed of sound [m/s] — Wells (1977)
const BLOOD_SOUND_SPEED: f64 = 1_540.0;

/// Inlet pressure [Pa] — mean brachial-artery pressure after catheter losses.
/// Normal MAP ≈ 93 mmHg ≈ 12 400 Pa; we use a round 12 000 Pa.
const P_INLET: f64 = 12_000.0;
/// Outlet pressure [Pa] — collection reservoir at slight positive gauge.
const P_OUTLET: f64 = 1_000.0;
/// Water/blood vapour pressure at 37 °C [Pa] — NIST.
const P_VAPOR_37C: f64 = 6_274.0;

/// Chip bounding-box length [mm] — millifluidic scale.
const CHIP_LENGTH_MM: f64 = 60.0;
/// Chip bounding-box width [mm].
const CHIP_WIDTH_MM: f64 = 30.0;
/// Coordinate scale: scheme mm → SI metres.
const SCALE_MM_TO_M: f64 = 1e-3;

// ═══════════════════════════════════════════════════════════════════════════════
// Design specification
// ═══════════════════════════════════════════════════════════════════════════════

struct DesignSpec {
    name: &'static str,
    tag: &'static str,
    description: &'static str,
    channel_config: ChannelTypeConfig,
}

fn design_specs() -> Vec<DesignSpec> {
    vec![
        DesignSpec {
            name: "Straight Bifurcation",
            tag: "straight",
            description: "Baseline — all straight channels; minimal shear enhancement.",
            channel_config: ChannelTypeConfig::AllStraight,
        },
        DesignSpec {
            name: "Serpentine Bifurcation",
            tag: "serpentine",
            description:
                "Serpentine mixing channels; elevated wall shear from Dean-type secondary flows.",
            channel_config: ChannelTypeConfig::AllSerpentine(SerpentineConfig {
                fill_factor: 0.7,
                wavelength_factor: 4.0,
                gaussian_width_factor: 8.0,
                wave_density_factor: 2.0,
                wave_phase_direction: 0.0,
                ..SerpentineConfig::default()
            }),
        },
        DesignSpec {
            name: "Venturi Bifurcation",
            tag: "venturi",
            description:
                "Frustum channels (converging–diverging) for hydrodynamic cavitation analysis.",
            channel_config: ChannelTypeConfig::AllFrustum(FrustumConfig {
                inlet_width: 2.0,    // mm
                throat_width: 0.5,   // mm  — 4 : 1 constriction ratio
                outlet_width: 2.0,   // mm
                taper_profile: TaperProfile::Smooth,
                smoothness: 50,
                throat_position: 0.5,
            }),
        },
    ]
}

// ═══════════════════════════════════════════════════════════════════════════════
// Blood fluid model (constant-property Newtonian approximation)
// ═══════════════════════════════════════════════════════════════════════════════

fn newtonian_blood() -> ConstantPropertyFluid<f64> {
    ConstantPropertyFluid::new(
        "Blood (Newtonian, 37 °C)".into(),
        BLOOD_DENSITY,
        BLOOD_VISCOSITY,
        BLOOD_CP,
        BLOOD_K,
        BLOOD_SOUND_SPEED,
    )
}

// ═══════════════════════════════════════════════════════════════════════════════
// Analysis result types
// ═══════════════════════════════════════════════════════════════════════════════

#[derive(Debug, Clone)]
struct ChannelAnalysis {
    channel_id: String,
    channel_type_label: String,
    flow_rate_m3s: f64,
    velocity_ms: f64,
    reynolds: f64,
    wall_shear_rate: f64,
    wall_shear_stress: f64,
    pressure_drop_pa: f64,
    length_m: f64,
    hydraulic_diameter_m: f64,
}

#[derive(Debug, Clone)]
struct VenturiResult {
    channel_id: String,
    inlet_velocity_ms: f64,
    throat_velocity_ms: f64,
    inlet_pressure_pa: f64,
    throat_pressure_pa: f64,
    cavitation_number: f64,
    is_cavitating: bool,
    throat_wall_shear_pa: f64,
}

#[derive(Debug, Clone)]
struct DesignResult {
    name: String,
    tag: String,
    description: String,
    node_count: usize,
    edge_count: usize,
    total_inlet_flow_m3s: f64,
    channels: Vec<ChannelAnalysis>,
    venturi_results: Vec<VenturiResult>,
    converged: bool,
    /// The original scheme geometry — kept for heatmap visualisation.
    system: Option<ChannelSystem>,
}

// ═══════════════════════════════════════════════════════════════════════════════
// Simulation pipeline
// ═══════════════════════════════════════════════════════════════════════════════

fn run_simulation(spec: &DesignSpec) -> Result<DesignResult, Box<dyn std::error::Error>> {
    // ── 1. Create scheme geometry ────────────────────────────────────────
    let geometry_config = GeometryConfig::new(0.5, 1.0, 0.5, GeometryGenerationConfig::default())?;
    let system = create_geometry(
        (CHIP_LENGTH_MM, CHIP_WIDTH_MM),
        &[SplitType::Bifurcation],
        &geometry_config,
        &spec.channel_config,
    );

    // ── 2. Convert scheme → cfd-1d network ──────────────────────────────
    let system_clone = system.clone();
    let blueprint = system.to_blueprint(SCALE_MM_TO_M).expect("blueprint");
    println!("    Network: {} nodes, {} channels", blueprint.nodes.len(), blueprint.channels.len());

    let blood = newtonian_blood();
    let mut network = network_from_blueprint(&blueprint, blood)?;

    // ── 3. Identify inlets / outlets and apply BCs ──────────────────────
    let mut inlets: Vec<NodeIndex> = Vec::new();
    let mut outlets: Vec<NodeIndex> = Vec::new();

    for ni in network.graph.node_indices() {
        match network.graph[ni].node_type {
            NodeType::Inlet => inlets.push(ni),
            NodeType::Outlet => outlets.push(ni),
            _ => {}
        }
    }

    for &inlet in &inlets {
        network.set_pressure(inlet, P_INLET);
    }
    for &outlet in &outlets {
        network.set_pressure(outlet, P_OUTLET);
    }

    // ── 4. Solve ────────────────────────────────────────────────────────
    let problem = NetworkProblem::new(network);
    let solver = NetworkSolver::<f64, ConstantPropertyFluid<f64>>::with_config(SolverConfig {
        tolerance: 1e-8,
        max_iterations: 500,
    });

    let solution = match solver.solve_network(&problem) {
        Ok(net) => net,
        Err(e) => {
            eprintln!("    Solver failed for {}: {}", spec.name, e);
            return Ok(DesignResult {
                name: spec.name.into(),
                tag: spec.tag.into(),
                description: spec.description.into(),
                node_count: 0,
                edge_count: 0,
                total_inlet_flow_m3s: 0.0,
                channels: vec![],
                venturi_results: vec![],
                converged: false,
                system: None,
            });
        }
    };

    // ── 5. Post-process: shear & cavitation analysis ────────────────────
    let pressures = solution.pressures();
    let flow_rates = solution.flow_rates();

    // Total inlet flow (sum of flow leaving inlet nodes)
    let total_inlet_flow: f64 = solution
        .graph
        .edge_indices()
        .filter_map(|ei| {
            let (from, _) = solution.graph.edge_endpoints(ei)?;
            if inlets.contains(&from) {
                Some(flow_rates.get(&ei).copied().unwrap_or(0.0).abs())
            } else {
                None
            }
        })
        .sum();

    let mut channels = Vec::new();
    let mut venturi_results = Vec::new();

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

        // Retrieve geometry from edge properties
        let props = solution.properties.get(&ei);
        let length = props.map_or(0.01, |p| p.length);
        let (area, dh) = props.map_or((1e-6, 1e-3), |p| {
            (p.area, p.hydraulic_diameter.unwrap_or(1e-3))
        });

        let velocity = if area > 0.0 { q / area } else { 0.0 };
        let reynolds = BLOOD_DENSITY * velocity * dh / BLOOD_VISCOSITY;

        // Wall shear: τ_w = μ · (8 V / D_h)  —  Hagen-Poiseuille analog
        let wall_shear_rate = 8.0 * velocity / dh;
        let wall_shear_stress = BLOOD_VISCOSITY * wall_shear_rate;

        // Classify channel type for labelling
        let ch_type_label = if let Some(p) = props {
            p.geometry
                .as_ref()
                .map_or("unknown".into(), |g| match &g.channel_type {
                    CfdChannelType::Straight => "straight".into(),
                    CfdChannelType::Serpentine { turns } => format!("serpentine({turns}t)"),
                    CfdChannelType::Curved { .. } => "curved".into(),
                    CfdChannelType::Tapered => "tapered/venturi".into(),
                    CfdChannelType::Spiral { .. } => "spiral".into(),
                })
        } else {
            "unknown".into()
        };

        channels.push(ChannelAnalysis {
            channel_id: edge.id.clone(),
            channel_type_label: ch_type_label,
            flow_rate_m3s: q,
            velocity_ms: velocity,
            reynolds,
            wall_shear_rate,
            wall_shear_stress,
            pressure_drop_pa: dp,
            length_m: length,
            hydraulic_diameter_m: dh,
        });

        // ── Venturi / cavitation analysis for tapered channels ──────────
        if let Some(p) = props {
            if let Some(geom) = &p.geometry {
                if matches!(geom.channel_type, CfdChannelType::Tapered) {
                    // Determine throat constriction from geometric variations
                    let min_scale = geom
                        .variations
                        .iter()
                        .map(|v| v.scale_factor)
                        .fold(f64::INFINITY, f64::min);

                    // Reconstruct throat dimensions from cross-section + min scale
                    let (inlet_w, h) = match &geom.cross_section {
                        CrossSection::Rectangular { width, height } => (*width, *height),
                        _ => (0.001, 0.0005),
                    };

                    // The cross-section stored is the *average* width.
                    // Original inlet width = avg / avg_scale ≈ avg * 3/(1+s+1)
                    // With default frustum (2, 0.5, 2 mm → avg = 1.5 mm),
                    // inlet_w stored = avg_width ≈ 1.5 mm.
                    // We recover approximate throat from scale factor.
                    let throat_w = if min_scale < 1.0 && min_scale > 0.0 {
                        inlet_w * min_scale
                    } else {
                        inlet_w * 0.25 // fallback 4:1
                    };

                    let inlet_area_est = inlet_w * h;
                    let throat_area = throat_w * h;

                    // Inlet velocity (through the average cross-section)
                    let v_inlet = if inlet_area_est > 0.0 { q / inlet_area_est } else { 0.0 };
                    // Throat velocity from continuity
                    let v_throat = if throat_area > 0.0 { q / throat_area } else { 0.0 };

                    // Bernoulli: p_throat = p_from − ½ρ(V_throat² − V_inlet²)
                    let p_throat =
                        p_from - 0.5 * BLOOD_DENSITY * (v_throat.powi(2) - v_inlet.powi(2));

                    // Cavitation number σ = (p_throat − p_v) / (½ρV_throat²)
                    let sigma = if v_throat > 1e-12 {
                        (p_throat - P_VAPOR_37C) / (0.5 * BLOOD_DENSITY * v_throat.powi(2))
                    } else {
                        f64::INFINITY
                    };

                    // Wall shear at throat: τ = μ · 8V_throat / D_h,throat
                    let dh_throat = 2.0 * throat_w * h / (throat_w + h);
                    let tau_throat = BLOOD_VISCOSITY * 8.0 * v_throat / dh_throat;

                    venturi_results.push(VenturiResult {
                        channel_id: edge.id.clone(),
                        inlet_velocity_ms: v_inlet,
                        throat_velocity_ms: v_throat,
                        inlet_pressure_pa: p_from,
                        throat_pressure_pa: p_throat,
                        cavitation_number: sigma,
                        is_cavitating: sigma < 1.0,
                        throat_wall_shear_pa: tau_throat,
                    });
                }
            }
        }
    }

    Ok(DesignResult {
        name: spec.name.into(),
        tag: spec.tag.into(),
        description: spec.description.into(),
        node_count: solution.node_count(),
        edge_count: solution.edge_count(),
        total_inlet_flow_m3s: total_inlet_flow,
        channels,
        venturi_results,
        converged: true,
        system: Some(system_clone),
    })
}

// ═══════════════════════════════════════════════════════════════════════════════
// Console output
// ═══════════════════════════════════════════════════════════════════════════════

fn print_results(result: &DesignResult) {
    println!();
    println!(
        "  ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓"
    );
    let pad = 56_usize.saturating_sub(result.name.len());
    println!("  ┃  {}{}┃", result.name, " ".repeat(pad));
    println!(
        "  ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛"
    );
    println!("    {}", result.description);
    println!(
        "    Nodes: {}   Edges: {}   Converged: {}",
        result.node_count, result.edge_count, result.converged
    );
    println!(
        "    Total inlet flow: {:.4e} m³/s  ({:.3} mL/min)",
        result.total_inlet_flow_m3s,
        result.total_inlet_flow_m3s * 1e6 * 60.0
    );

    // ── Shear-stress table ───────────────────────────────────────────────
    println!();
    println!(
        "    ┌───────────┬────────────────┬────────────┬──────────┬─────────┬────────────┬────────────┐"
    );
    println!(
        "    │ Channel   │ Type           │  Q [mL/m]  │ V [m/s]  │  Re     │  γ̇ [1/s]  │ τ_w [Pa]   │"
    );
    println!(
        "    ├───────────┼────────────────┼────────────┼──────────┼─────────┼────────────┼────────────┤"
    );
    for ch in &result.channels {
        let q_ml = ch.flow_rate_m3s * 1e6 * 60.0;
        println!(
            "    │ {:>9} │ {:>14} │ {:>10.4} │ {:>8.4} │ {:>7.1} │ {:>10.0} │ {:>10.4} │",
            ch.channel_id,
            ch.channel_type_label,
            q_ml,
            ch.velocity_ms,
            ch.reynolds,
            ch.wall_shear_rate,
            ch.wall_shear_stress,
        );
    }
    println!(
        "    └───────────┴────────────────┴────────────┴──────────┴─────────┴────────────┴────────────┘"
    );

    // ── Venturi / cavitation table ───────────────────────────────────────
    if !result.venturi_results.is_empty() {
        println!();
        println!("    Venturi Throat Analysis — Hydrodynamic Cavitation Assessment");
        println!(
            "    ┌───────────┬────────────┬────────────┬────────────┬────────────┬──────────────┬────────────┐"
        );
        println!(
            "    │ Channel   │ V_in [m/s] │ V_thr[m/s] │ p_thr [Pa] │   σ        │ Status       │ τ_thr [Pa] │"
        );
        println!(
            "    ├───────────┼────────────┼────────────┼────────────┼────────────┼──────────────┼────────────┤"
        );
        for v in &result.venturi_results {
            let status = if v.is_cavitating {
                "CAVITATING"
            } else if v.cavitation_number < 2.0 {
                "Marginal"
            } else {
                "Safe"
            };
            println!(
                "    │ {:>9} │ {:>10.4} │ {:>10.4} │ {:>10.0} │ {:>10.3} │ {:>12} │ {:>10.3} │",
                v.channel_id,
                v.inlet_velocity_ms,
                v.throat_velocity_ms,
                v.throat_pressure_pa,
                v.cavitation_number,
                status,
                v.throat_wall_shear_pa,
            );
        }
        println!(
            "    └───────────┴────────────┴────────────┴────────────┴────────────┴──────────────┴────────────┘"
        );
    }
}

fn print_comparison(results: &[DesignResult]) {
    println!();
    println!("  ╔══════════════════════════════════════════════════════════════════════╗");
    println!("  ║  Cross-Design Comparison Summary                                   ║");
    println!("  ╚══════════════════════════════════════════════════════════════════════╝");
    println!();

    println!(
        "    {:30} {:>12} {:>12} {:>12} {:>12}",
        "Design", "τ_w max [Pa]", "τ_w avg [Pa]", "Q_tot [mL/m]", "Re_max"
    );
    println!("    {}", "─".repeat(80));

    for r in results {
        if !r.converged || r.channels.is_empty() {
            continue;
        }
        let max_shear = r
            .channels
            .iter()
            .map(|c| c.wall_shear_stress)
            .fold(0.0_f64, f64::max);
        let avg_shear =
            r.channels.iter().map(|c| c.wall_shear_stress).sum::<f64>() / r.channels.len() as f64;
        let q_tot = r.total_inlet_flow_m3s * 1e6 * 60.0;
        let re_max = r
            .channels
            .iter()
            .map(|c| c.reynolds)
            .fold(0.0_f64, f64::max);

        println!(
            "    {:30} {:>12.4} {:>12.4} {:>12.4} {:>12.1}",
            r.name, max_shear, avg_shear, q_tot, re_max
        );
    }

    // Physiological shear thresholds
    println!();
    println!("    Physiological reference thresholds:");
    println!("      Normal arterial wall shear stress:   1 – 7 Pa");
    println!("      Platelet activation threshold:       > 10 Pa");
    println!("      Red-cell damage (haemolysis):        > 150 Pa");
    println!("      Endothelial pathological low shear:  < 0.4 Pa");

    // Cavitation summary
    let venturi_result = results.iter().find(|r| !r.venturi_results.is_empty());
    if let Some(vr) = venturi_result {
        println!();
        println!("    Cavitation Summary — {} design:", vr.name);
        let min_sigma = vr
            .venturi_results
            .iter()
            .map(|v| v.cavitation_number)
            .fold(f64::INFINITY, f64::min);
        let max_v_throat = vr
            .venturi_results
            .iter()
            .map(|v| v.throat_velocity_ms)
            .fold(0.0_f64, f64::max);
        let cavitating = vr
            .venturi_results
            .iter()
            .filter(|v| v.is_cavitating)
            .count();

        println!("      Min cavitation number σ   = {:.3}", min_sigma);
        println!("      Max throat velocity        = {:.4} m/s", max_v_throat);
        println!("      Max throat wall shear      = {:.3} Pa",
            vr.venturi_results.iter().map(|v| v.throat_wall_shear_pa).fold(0.0_f64, f64::max));
        println!(
            "      Channels cavitating        = {} / {}",
            cavitating,
            vr.venturi_results.len()
        );
        println!();

        if min_sigma > 2.0 {
            println!("      Result: No cavitation (sigma > 2.0 for all channels)");
        } else if min_sigma > 1.0 {
            println!("      Result: MARGINAL — onset region (1.0 < sigma < 2.0)");
            println!("      -> Consider widening throat or reducing inlet pressure.");
        } else {
            println!("      Result: ACTIVE CAVITATION (sigma < 1.0)");
            println!("      -> Vapour-pressure threshold exceeded at throat.");
            println!("      -> This may be intentional for lysis / mixing applications.");
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// SVG plot generation
// ═══════════════════════════════════════════════════════════════════════════════

/// Generate an SVG bar chart comparing wall shear stress across designs.
fn plot_shear_comparison(results: &[DesignResult]) -> Result<(), Box<dyn std::error::Error>> {
    let valid: Vec<&DesignResult> = results.iter().filter(|r| r.converged && !r.channels.is_empty()).collect();
    if valid.is_empty() {
        return Ok(());
    }

    let path = "outputs/blood_shear_comparison.svg";
    std::fs::create_dir_all("outputs")?;

    let root = SVGBackend::new(path, (900, 500)).into_drawing_area();
    root.fill(&WHITE)?;

    // Collect data: (design_name, max_shear, avg_shear)
    let data: Vec<(&str, f64, f64)> = valid
        .iter()
        .map(|r| {
            let max_s = r.channels.iter().map(|c| c.wall_shear_stress).fold(0.0_f64, f64::max);
            let avg_s =
                r.channels.iter().map(|c| c.wall_shear_stress).sum::<f64>() / r.channels.len() as f64;
            (r.name.as_str(), max_s, avg_s)
        })
        .collect();

    let y_max = data
        .iter()
        .map(|(_, m, _)| *m)
        .fold(0.0_f64, f64::max)
        * 1.3;

    let mut chart = ChartBuilder::on(&root)
        .caption(
            "Wall Shear Stress — Blood in Millifluidic Designs",
            ("sans-serif", 20),
        )
        .margin(15)
        .x_label_area_size(60)
        .y_label_area_size(60)
        .build_cartesian_2d(0.0..(data.len() as f64), 0.0..y_max)?;

    chart
        .configure_mesh()
        .y_desc("Wall shear stress [Pa]")
        .x_desc("Design")
        .x_label_formatter(&|x| {
            let idx = *x as usize;
            if idx < data.len() {
                data[idx].0.to_string()
            } else {
                String::new()
            }
        })
        .draw()?;

    // Max shear bars
    chart.draw_series(data.iter().enumerate().map(|(i, (_, max_s, _))| {
        let x0 = i as f64 + 0.1;
        let x1 = i as f64 + 0.45;
        Rectangle::new([(x0, 0.0), (x1, *max_s)], RED.mix(0.7).filled())
    }))?
    .label("Max τ_w")
    .legend(|(x, y)| Rectangle::new([(x, y - 5), (x + 15, y + 5)], RED.mix(0.7).filled()));

    // Avg shear bars
    chart.draw_series(data.iter().enumerate().map(|(i, (_, _, avg_s))| {
        let x0 = i as f64 + 0.55;
        let x1 = i as f64 + 0.9;
        Rectangle::new([(x0, 0.0), (x1, *avg_s)], BLUE.mix(0.6).filled())
    }))?
    .label("Avg τ_w")
    .legend(|(x, y)| Rectangle::new([(x, y - 5), (x + 15, y + 5)], BLUE.mix(0.6).filled()));

    // Physiological threshold line: platelet activation = 10 Pa
    let threshold = 10.0_f64;
    if threshold < y_max {
        chart.draw_series(std::iter::once(PathElement::new(
            vec![(0.0, threshold), (data.len() as f64, threshold)],
            RED.stroke_width(2),
        )))?
        .label("Platelet activation (10 Pa)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 15, y)], RED.stroke_width(2)));
    }

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .position(SeriesLabelPosition::UpperRight)
        .draw()?;

    root.present()?;
    println!("    Shear comparison plot saved to {path}");
    Ok(())
}

/// Generate an SVG chart of the venturi cavitation analysis.
fn plot_cavitation(results: &[DesignResult]) -> Result<(), Box<dyn std::error::Error>> {
    let venturi = results
        .iter()
        .find(|r| !r.venturi_results.is_empty() && r.converged);
    let venturi = match venturi {
        Some(v) => v,
        None => return Ok(()),
    };

    let path = "outputs/blood_cavitation_analysis.svg";
    std::fs::create_dir_all("outputs")?;

    let root = SVGBackend::new(path, (800, 500)).into_drawing_area();
    root.fill(&WHITE)?;

    let n = venturi.venturi_results.len();
    if n == 0 {
        return Ok(());
    }

    let max_sigma = venturi
        .venturi_results
        .iter()
        .map(|v| v.cavitation_number)
        .fold(0.0_f64, f64::max)
        * 1.3;
    let y_max = max_sigma.max(3.0);

    let mut chart = ChartBuilder::on(&root)
        .caption(
            "Venturi Cavitation Number — Blood Flow",
            ("sans-serif", 20),
        )
        .margin(15)
        .x_label_area_size(50)
        .y_label_area_size(60)
        .build_cartesian_2d(0.0..(n as f64), 0.0..y_max)?;

    chart
        .configure_mesh()
        .y_desc("Cavitation number σ")
        .x_desc("Venturi channel index")
        .draw()?;

    // Cavitation number bars
    chart.draw_series(venturi.venturi_results.iter().enumerate().map(|(i, v)| {
        let colour = if v.is_cavitating {
            RED.mix(0.8)
        } else if v.cavitation_number < 2.0 {
            YELLOW.mix(0.7)
        } else {
            GREEN.mix(0.6)
        };
        let x0 = i as f64 + 0.15;
        let x1 = i as f64 + 0.85;
        Rectangle::new([(x0, 0.0), (x1, v.cavitation_number)], colour.filled())
    }))?;

    // Threshold lines
    chart.draw_series(std::iter::once(PathElement::new(
        vec![(0.0, 1.0), (n as f64, 1.0)],
        RED.stroke_width(2),
    )))?
    .label("σ = 1.0 (cavitation onset)")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 15, y)], RED.stroke_width(2)));

    chart.draw_series(std::iter::once(PathElement::new(
        vec![(0.0, 2.0), (n as f64, 2.0)],
        YELLOW.stroke_width(2),
    )))?
    .label("σ = 2.0 (marginal)")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 15, y)], YELLOW.stroke_width(2)));

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .position(SeriesLabelPosition::UpperRight)
        .draw()?;

    root.present()?;
    println!("    Cavitation analysis plot saved to {path}");
    Ok(())
}

/// Generate per-channel velocity & shear comparison across designs.
fn plot_velocity_profiles(results: &[DesignResult]) -> Result<(), Box<dyn std::error::Error>> {
    let valid: Vec<&DesignResult> = results.iter().filter(|r| r.converged && !r.channels.is_empty()).collect();
    if valid.is_empty() {
        return Ok(());
    }

    let path = "outputs/blood_velocity_comparison.svg";
    std::fs::create_dir_all("outputs")?;

    let root = SVGBackend::new(path, (900, 500)).into_drawing_area();
    root.fill(&WHITE)?;

    // For each design, collect (channel_index, velocity) pairs
    let max_v: f64 = valid
        .iter()
        .flat_map(|r| r.channels.iter().map(|c| c.velocity_ms))
        .fold(0.0_f64, f64::max)
        * 1.3;

    let max_ch: usize = valid.iter().map(|r| r.channels.len()).max().unwrap_or(1);

    let mut chart = ChartBuilder::on(&root)
        .caption(
            "Channel Velocity — Blood in Millifluidic Designs",
            ("sans-serif", 20),
        )
        .margin(15)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(0usize..max_ch, 0.0..max_v)?;

    chart
        .configure_mesh()
        .y_desc("Mean velocity [m/s]")
        .x_desc("Channel index")
        .draw()?;

    let colors = [&RED, &BLUE, &GREEN];
    for (di, r) in valid.iter().enumerate() {
        let color = colors[di % colors.len()];
        let points: Vec<(usize, f64)> = r
            .channels
            .iter()
            .enumerate()
            .map(|(i, c)| (i, c.velocity_ms))
            .collect();

        chart
            .draw_series(LineSeries::new(points.clone(), color.stroke_width(2)))?
            .label(&r.name)
            .legend(move |(x, y)| {
                PathElement::new(vec![(x, y), (x + 15, y)], color.stroke_width(2))
            });

        chart.draw_series(points.iter().map(|&(x, y)| Circle::new((x, y), 4, color.filled())))?;
    }

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .position(SeriesLabelPosition::UpperRight)
        .draw()?;

    root.present()?;
    println!("    Velocity comparison plot saved to {path}");
    Ok(())
}

// ═══════════════════════════════════════════════════════════════════════════════
// Schematic heatmap — 1D field values overlaid on 2D chip layout
// ═══════════════════════════════════════════════════════════════════════════════

/// Jet-like colormap: blue → cyan → green → yellow → red.
/// Input `t` is clamped to \[0, 1\].
fn jet_colormap(t: f64) -> RGBColor {
    let t = t.clamp(0.0, 1.0);
    let (r, g, b) = if t < 0.25 {
        let s = t / 0.25;
        (0.0, s, 1.0)
    } else if t < 0.5 {
        let s = (t - 0.25) / 0.25;
        (0.0, 1.0, 1.0 - s)
    } else if t < 0.75 {
        let s = (t - 0.5) / 0.25;
        (s, 1.0, 0.0)
    } else {
        let s = (t - 0.75) / 0.25;
        (1.0, 1.0 - s, 0.0)
    };
    RGBColor(
        (r * 255.0) as u8,
        (g * 255.0) as u8,
        (b * 255.0) as u8,
    )
}

/// Extract the 2D centre-line path for a single channel.
fn channel_path_2d(system: &ChannelSystem, ch_idx: usize) -> Vec<(f64, f64)> {
    let ch = &system.channels[ch_idx];
    match &ch.channel_type {
        SchemeChannelType::Straight => {
            vec![
                system.nodes[ch.from_node].point,
                system.nodes[ch.to_node].point,
            ]
        }
        SchemeChannelType::SmoothStraight { path }
        | SchemeChannelType::Serpentine { path }
        | SchemeChannelType::Arc { path } => path.clone(),
        SchemeChannelType::Frustum { path, .. } => path.clone(),
    }
}

/// Draw a schematic heatmap: the chip layout with channels coloured by a
/// per-channel scalar field.
///
/// The resulting SVG shows the 2D chip outline from `scheme`, with each
/// channel drawn using a jet-like colour derived from `channel_values`.
/// A colour bar on the right maps value → colour.
///
/// For frustum (venturi) channels, the stroke width varies along the path
/// to visually represent the converging–diverging geometry.
fn plot_schematic_heatmap(
    system: &ChannelSystem,
    channel_values: &[f64],
    title: &str,
    field_label: &str,
    output_path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    std::fs::create_dir_all("outputs")?;

    let (box_w, box_h) = system.box_dims;
    let margin_mm = 5.0;

    let img_w: u32 = 1000;
    let img_h: u32 = 600;

    let root = SVGBackend::new(output_path, (img_w, img_h)).into_drawing_area();
    root.fill(&WHITE)?;

    // Split: ~82 % chart, ~18 % colour bar
    let (chart_area, cbar_area) = root.split_horizontally(img_w * 82 / 100);

    // ── Value range ─────────────────────────────────────────────────────
    let v_min = channel_values
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    let v_max = channel_values
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let v_range = if (v_max - v_min).abs() < 1e-12 {
        1.0
    } else {
        v_max - v_min
    };

    // ── Build the chart area (chip coordinates, mm) ─────────────────────
    let mut chart = ChartBuilder::on(&chart_area)
        .caption(title, ("sans-serif", 17))
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(40)
        .build_cartesian_2d(
            (-margin_mm)..(box_w + margin_mm),
            (-margin_mm)..(box_h + margin_mm),
        )?;

    chart
        .configure_mesh()
        .x_desc("x [mm]")
        .y_desc("y [mm]")
        .light_line_style(WHITE.mix(0.0)) // hide grid
        .draw()?;

    // ── Draw chip boundary ──────────────────────────────────────────────
    let outline: Vec<(f64, f64)> = vec![
        (0.0, 0.0),
        (box_w, 0.0),
        (box_w, box_h),
        (0.0, box_h),
        (0.0, 0.0),
    ];
    chart.draw_series(std::iter::once(PathElement::new(
        outline,
        ShapeStyle {
            color: RGBColor(180, 180, 180).to_rgba(),
            filled: false,
            stroke_width: 2,
        },
    )))?;

    // ── Draw each channel, coloured by field value ──────────────────────
    let n_ch = channel_values.len().min(system.channels.len());
    for ch_idx in 0..n_ch {
        let value = channel_values[ch_idx];
        let t = (value - v_min) / v_range;
        let color = jet_colormap(t);

        let ch = &system.channels[ch_idx];

        match &ch.channel_type {
            // For frustum channels, vary the stroke width along the path
            // to show the converging–diverging geometry.
            SchemeChannelType::Frustum { path, widths, .. } if path.len() >= 2 => {
                let w_max = widths.iter().copied().fold(0.0_f64, f64::max).max(0.1);
                for i in 0..path.len() - 1 {
                    let w_avg = (widths.get(i).copied().unwrap_or(1.0)
                        + widths.get(i + 1).copied().unwrap_or(1.0))
                        / 2.0;
                    // Scale width to stroke pixels: ~6 px at maximum width
                    let sw = ((w_avg / w_max * 6.0) as u32).max(1).min(12);
                    chart.draw_series(std::iter::once(PathElement::new(
                        vec![path[i], path[i + 1]],
                        ShapeStyle {
                            color: color.to_rgba(),
                            filled: false,
                            stroke_width: sw,
                        },
                    )))?;
                }
            }

            // All other channel types: constant stroke width, full path.
            _ => {
                let pts = channel_path_2d(system, ch_idx);
                if pts.len() >= 2 {
                    chart.draw_series(std::iter::once(PathElement::new(
                        pts,
                        ShapeStyle {
                            color: color.to_rgba(),
                            filled: false,
                            stroke_width: 3,
                        },
                    )))?;
                }
            }
        }
    }

    // ── Draw nodes as small dots ────────────────────────────────────────
    for node in &system.nodes {
        chart.draw_series(std::iter::once(Circle::new(
            node.point,
            3,
            ShapeStyle {
                color: RGBColor(60, 60, 60).to_rgba(),
                filled: true,
                stroke_width: 1,
            },
        )))?;
    }

    // ── Colour bar (right-hand side) ────────────────────────────────────
    let bar_x0: i32 = 25;
    let bar_x1: i32 = 55;
    let bar_y0: i32 = 40;
    let bar_y1: i32 = (img_h as i32) - 60;
    let bar_height = bar_y1 - bar_y0;
    let n_steps: usize = 64;

    for i in 0..n_steps {
        let t = 1.0 - (i as f64 / n_steps as f64);
        let c = jet_colormap(t);
        let y_start = bar_y0 + (i as i32 * bar_height / n_steps as i32);
        let y_end = bar_y0 + ((i + 1) as i32 * bar_height / n_steps as i32);
        cbar_area.draw(&Rectangle::new(
            [(bar_x0, y_start), (bar_x1, y_end)],
            c.filled(),
        ))?;
    }

    // Colour bar outline
    cbar_area.draw(&Rectangle::new(
        [(bar_x0, bar_y0), (bar_x1, bar_y1)],
        BLACK.stroke_width(1),
    ))?;

    // Five evenly-spaced tick labels
    let tick_style = ("sans-serif", 11).into_text_style(&cbar_area);
    for tick_i in 0..=4u32 {
        let frac = tick_i as f64 / 4.0;
        let val = v_max - frac * v_range;
        let y_pos = bar_y0 + (frac * bar_height as f64) as i32;
        cbar_area.draw_text(
            &format!("{:.2}", val),
            &tick_style,
            (bar_x1 + 4, y_pos - 6),
        )?;
        // Small tick mark
        cbar_area.draw(&PathElement::new(
            vec![(bar_x1, y_pos), (bar_x1 + 3, y_pos)],
            BLACK.stroke_width(1),
        ))?;
    }

    // Field label below the bar
    let label_style = ("sans-serif", 12).into_text_style(&cbar_area);
    cbar_area.draw_text(field_label, &label_style, (bar_x0 - 5, bar_y1 + 8))?;

    root.present()?;
    println!("    Heatmap saved to {output_path}");
    Ok(())
}

/// Generate all schematic heatmap views for one design.
fn plot_design_heatmaps(result: &DesignResult) -> Result<(), Box<dyn std::error::Error>> {
    let system = match &result.system {
        Some(s) => s,
        None => return Ok(()),
    };
    let n_ch = result.channels.len().min(system.channels.len());
    if n_ch == 0 {
        return Ok(());
    }

    // ── Wall shear stress ───────────────────────────────────────────────
    let shear: Vec<f64> = result.channels[..n_ch]
        .iter()
        .map(|c| c.wall_shear_stress)
        .collect();
    plot_schematic_heatmap(
        system,
        &shear,
        &format!("{} — Wall Shear Stress", result.name),
        "τ_w [Pa]",
        &format!("outputs/heatmap_{}_shear.svg", result.tag),
    )?;

    // ── Velocity ────────────────────────────────────────────────────────
    let velocity: Vec<f64> = result.channels[..n_ch]
        .iter()
        .map(|c| c.velocity_ms)
        .collect();
    plot_schematic_heatmap(
        system,
        &velocity,
        &format!("{} — Mean Channel Velocity", result.name),
        "V [m/s]",
        &format!("outputs/heatmap_{}_velocity.svg", result.tag),
    )?;

    // ── Pressure drop ───────────────────────────────────────────────────
    let dp: Vec<f64> = result.channels[..n_ch]
        .iter()
        .map(|c| c.pressure_drop_pa)
        .collect();
    plot_schematic_heatmap(
        system,
        &dp,
        &format!("{} — Channel Pressure Drop", result.name),
        "ΔP [Pa]",
        &format!("outputs/heatmap_{}_pressure_drop.svg", result.tag),
    )?;

    // ── Reynolds number ─────────────────────────────────────────────────
    let re: Vec<f64> = result.channels[..n_ch]
        .iter()
        .map(|c| c.reynolds)
        .collect();
    plot_schematic_heatmap(
        system,
        &re,
        &format!("{} — Reynolds Number", result.name),
        "Re",
        &format!("outputs/heatmap_{}_reynolds.svg", result.tag),
    )?;

    // ── Venturi-specific fields ─────────────────────────────────────────
    if !result.venturi_results.is_empty() {
        let venturi_map: HashMap<String, &VenturiResult> = result
            .venturi_results
            .iter()
            .map(|v| (v.channel_id.clone(), v))
            .collect();

        // Cavitation number
        let cav: Vec<f64> = result.channels[..n_ch]
            .iter()
            .map(|c| {
                venturi_map
                    .get(&c.channel_id)
                    .map_or(0.0, |v| v.cavitation_number)
            })
            .collect();
        plot_schematic_heatmap(
            system,
            &cav,
            &format!("{} — Cavitation Number σ", result.name),
            "σ",
            &format!("outputs/heatmap_{}_cavitation.svg", result.tag),
        )?;

        // Throat wall shear
        let throat_shear: Vec<f64> = result.channels[..n_ch]
            .iter()
            .map(|c| {
                venturi_map
                    .get(&c.channel_id)
                    .map_or(0.0, |v| v.throat_wall_shear_pa)
            })
            .collect();
        plot_schematic_heatmap(
            system,
            &throat_shear,
            &format!("{} — Venturi Throat Shear", result.name),
            "τ_thr [Pa]",
            &format!("outputs/heatmap_{}_throat_shear.svg", result.tag),
        )?;
    }

    Ok(())
}

// ═══════════════════════════════════════════════════════════════════════════════
// Main
// ═══════════════════════════════════════════════════════════════════════════════

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!();
    println!("  ═══════════════════════════════════════════════════════════════════");
    println!("   Blood Flow Through Millifluidic Devices");
    println!("   Pipeline: scheme -> cfd-1d -> simulation -> analysis -> SVG");
    println!("   Blood source: brachial artery, arm without AV fistula");
    println!("  ═══════════════════════════════════════════════════════════════════");
    println!();
    println!("  Fluid model : Newtonian blood (high-shear approximation)");
    println!(
        "    rho  = {} kg/m^3    mu = {:.4e} Pa.s",
        BLOOD_DENSITY, BLOOD_VISCOSITY
    );
    println!("  Boundary conditions:");
    println!(
        "    Inlet  : {:.0} Pa  ({:.0} mmHg)",
        P_INLET,
        P_INLET / 133.322
    );
    println!(
        "    Outlet : {:.0} Pa  ({:.0} mmHg)",
        P_OUTLET,
        P_OUTLET / 133.322
    );
    println!(
        "  Chip     : {} x {} mm  (bifurcation topology)",
        CHIP_LENGTH_MM, CHIP_WIDTH_MM
    );
    println!(
        "  Vapour pressure (37 C) : {} Pa",
        P_VAPOR_37C
    );
    println!();

    // ── Run each design ─────────────────────────────────────────────────
    let specs = design_specs();
    let mut results = Vec::with_capacity(specs.len());

    for spec in &specs {
        println!("  ── Simulating: {} ──────────────────────", spec.name);
        match run_simulation(spec) {
            Ok(result) => {
                print_results(&result);
                results.push(result);
            }
            Err(e) => {
                eprintln!("    Error: {}", e);
            }
        }
    }

    // ── Cross-design comparison ─────────────────────────────────────────
    print_comparison(&results);

    // ── Generate SVG plots ──────────────────────────────────────────────
    println!();
    println!("  ── Generating plots ──────────────────────────────────────────");
    if let Err(e) = plot_shear_comparison(&results) {
        eprintln!("    Shear plot error: {}", e);
    }
    if let Err(e) = plot_cavitation(&results) {
        eprintln!("    Cavitation plot error: {}", e);
    }
    if let Err(e) = plot_velocity_profiles(&results) {
        eprintln!("    Velocity plot error: {}", e);
    }

    // ── Generate schematic heatmaps (1D field on 2D layout) ──────────
    println!("  ── Generating schematic heatmaps ────────────────────────────");
    for result in &results {
        if let Err(e) = plot_design_heatmaps(result) {
            eprintln!("    Heatmap error ({}): {}", result.name, e);
        }
    }

    println!();
    println!("  Done.");
    Ok(())
}
