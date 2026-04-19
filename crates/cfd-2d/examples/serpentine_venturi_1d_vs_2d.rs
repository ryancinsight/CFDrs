//! # Serpentine-Venturi 1D vs 2D Cross-Fidelity Comparison
//!
//! Builds a serpentine channel with venturi constrictions at each U-turn
//! bend (Dean-flow apex) using the `serpentine_venturi_rect` preset from
//! `cfd-schematics`, then solves the same geometry at two fidelity levels:
//!
//! 1. **cfd-1d** — lumped Hagen-Poiseuille network via `network_from_blueprint`
//! 2. **cfd-2d** — per-channel 2D SIMPLE N-S via `Network2DSolver::solve_all`
//!
//! The example prints a channel-by-channel comparison table of flow rates,
//! pressure drops, wall shear stresses, and outlet-flow agreement percentages,
//! quantifying the fidelity gap between the lumped 1D model and the resolved
//! 2D field solve.
//!
//! ## Physics
//!
//! At each serpentine U-turn the Dean secondary vortices (De = Re·√(D_h/2R))
//! pre-focus cells toward the centreline.  A venturi throat placed at this
//! Dean apex produces the maximum wall shear rate in the device — ideal for
//! cavitation-enhanced CTC lysis in sonodynamic therapy (SDT).
//!
//! The 1D model captures the resistance network (Kirchhoff pressure balance)
//! but cannot resolve the 2D velocity profiles within throats.  The 2D model
//! resolves intra-channel gradients, revealing the true wall shear which the
//! 1D model can only estimate from the Hagen-Poiseuille formula.
//!
//! ## Run
//!
//! ```sh
//! cargo run -p cfd-2d --example serpentine_venturi_1d_vs_2d --no-default-features
//! ```

#![allow(clippy::cast_precision_loss)]

use std::fs;
use std::path::PathBuf;

// ── cfd-schematics ──────────────────────────────────────────────────────────
use cfd_schematics::application::ports::GraphSink;
use cfd_schematics::interface::presets::serpentine_venturi_rect;

// ── cfd-1d ──────────────────────────────────────────────────────────────────
use cfd_1d::domain::network::network_from_blueprint;
use cfd_1d::{NetworkProblem, NetworkSolver, NodeType, SolverConfig};

// ── cfd-2d ──────────────────────────────────────────────────────────────────
use cfd_2d::network::{Network2DSolver, Network2dBuilderSink};
use cfd_2d::solvers::ns_fvm::BloodModel;

// ── cfd-core ────────────────────────────────────────────────────────────────
use cfd_core::physics::fluid::CassonBlood;

// ── petgraph ────────────────────────────────────────────────────────────────
use petgraph::visit::EdgeRef;

// ─────────────────────────────────────────────────────────────────────────────
// Physical constants
// ─────────────────────────────────────────────────────────────────────────────

/// Blood density [kg/m³]
const BLOOD_DENSITY: f64 = 1_060.0;
/// Blood effective viscosity [Pa·s] (Newtonian approximation for 2D SIMPLE)
const BLOOD_VISCOSITY: f64 = 3.5e-3;
/// Inlet volumetric flow rate [m³/s] — 10 mL/min (millifluidic regime)
const INLET_FLOW_M3_S: f64 = 10.0e-6 / 60.0;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let out = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("outputs")
        .join("serpentine_venturi_1d_vs_2d");
    fs::create_dir_all(&out)?;

    println!("Serpentine-Venturi 1D vs 2D Cross-Fidelity Comparison");
    println!("=====================================================\n");

    // ─────────────────────────────────────────────────────────────────────────
    // Phase 1: Build blueprint (cfd-schematics)
    // ─────────────────────────────────────────────────────────────────────────

    let segments = 4_usize;
    let segment_length_m = 0.015; // 15 mm
    let width_m = 0.002; // 2 mm channel width
    let throat_width_m = 0.001; // 1 mm throat
    let height_m = 0.001; // 1 mm depth
    let throat_length_m = 0.003; // 3 mm throat length
    let bend_radius_m = 0.001; // 1 mm bend radius

    let blueprint = serpentine_venturi_rect(
        "sv_crossfidelity",
        segments,
        segment_length_m,
        width_m,
        throat_width_m,
        height_m,
        throat_length_m,
        bend_radius_m,
    );

    let n_venturis = blueprint.venturi_channels().len();
    let d_h = 2.0 * width_m * height_m / (width_m + height_m);
    let u_mean = INLET_FLOW_M3_S / (width_m * height_m);
    let re = BLOOD_DENSITY * u_mean * d_h / BLOOD_VISCOSITY;
    let de = re * (d_h / (2.0 * bend_radius_m)).sqrt();

    println!("Blueprint: {segments} segments, {n_venturis} venturi throats");
    println!("  Channel width:   {:.1} mm", width_m * 1e3);
    println!("  Throat width:    {:.1} mm", throat_width_m * 1e3);
    println!("  Height:          {:.1} mm", height_m * 1e3);
    println!("  Bend radius:     {:.1} mm", bend_radius_m * 1e3);
    println!("  D_h:             {:.3} mm", d_h * 1e3);
    println!(
        "  Q_inlet:         {:.2} mL/min",
        INLET_FLOW_M3_S * 60.0 * 1e6
    );
    println!("  u_mean:          {:.4} m/s", u_mean);
    println!("  Re:              {:.1}", re);
    println!("  Dean number:     {:.1}", de);
    println!("  Contraction:     {:.1}:1", width_m / throat_width_m);
    println!(
        "  Nodes: {}  Channels: {}",
        blueprint.nodes.len(),
        blueprint.channels.len()
    );

    // ─────────────────────────────────────────────────────────────────────────
    // Phase 2: 1D solve (cfd-1d)
    // ─────────────────────────────────────────────────────────────────────────

    println!("\n--- Phase 2: cfd-1d Lumped Network Solve ---\n");

    let blood_1d = CassonBlood::<f64>::normal_blood();
    let mut network = network_from_blueprint(&blueprint, blood_1d)?;

    // Set boundary conditions: Neumann inlet (fixed flow), Dirichlet outlet (P = 0)
    let mut inlet_nodes = Vec::new();
    let mut outlet_nodes = Vec::new();
    for idx in network.graph.node_indices() {
        if let Some(node) = network.graph.node_weight(idx) {
            match node.node_type {
                NodeType::Inlet => inlet_nodes.push(idx),
                NodeType::Outlet => outlet_nodes.push(idx),
                _ => {}
            }
        }
    }
    let q_per_inlet = INLET_FLOW_M3_S / inlet_nodes.len() as f64;
    for &n in &inlet_nodes {
        network.set_neumann_flow(n, q_per_inlet);
    }
    for &n in &outlet_nodes {
        network.set_pressure(n, 0.0);
    }

    let config_1d = SolverConfig {
        tolerance: 1e-8,
        max_iterations: 200,
    };
    let solver_1d = NetworkSolver::<f64, CassonBlood<f64>>::with_config(config_1d);
    let problem = NetworkProblem::new(network);
    let solved = solver_1d.solve_network(&problem)?;

    // Extract 1D results per channel
    let inlet_p: f64 = inlet_nodes
        .iter()
        .filter_map(|idx| solved.pressures.get(idx.index()).copied())
        .sum::<f64>()
        / inlet_nodes.len() as f64;

    println!("  Inlet pressure:  {:.2} Pa", inlet_p);
    println!(
        "  Total dP:        {:.2} Pa (inlet to outlet = 0 Pa)",
        inlet_p
    );

    // Build a map of channel_id → (flow_rate, pressure_from, pressure_to)
    struct OneDChannelResult {
        id: String,
        flow_rate_m3_s: f64,
        pressure_drop_pa: f64,
        wall_shear_pa: f64,
    }
    let mut results_1d: Vec<OneDChannelResult> = Vec::new();

    for edge_ref in solved.graph.edge_references() {
        let eid = edge_ref.id();
        let edge = edge_ref.weight();
        let (from_idx, to_idx) = solved
            .graph
            .edge_endpoints(eid)
            .expect("edge must have endpoints");
        let q = solved
            .flow_rates
            .get(eid.index())
            .copied()
            .unwrap_or(edge.flow_rate);
        let p_from = solved
            .pressures
            .get(from_idx.index())
            .copied()
            .unwrap_or(0.0);
        let p_to = solved.pressures.get(to_idx.index()).copied().unwrap_or(0.0);
        let dp = (p_from - p_to).abs();

        // Estimate wall shear from 1D Poiseuille: τ_w = 6μQ / (w·h²) (rectangular)
        // Find matching blueprint channel for dimensions
        let (w, h) = blueprint
            .channels
            .iter()
            .find(|c| c.id.as_str() == edge.id)
            .map(|c| c.cross_section.dims())
            .unwrap_or((width_m, height_m));
        let tau_w = 6.0 * BLOOD_VISCOSITY * q.abs() / (w * h * h);

        results_1d.push(OneDChannelResult {
            id: edge.id.clone(),
            flow_rate_m3_s: q,
            pressure_drop_pa: dp,
            wall_shear_pa: tau_w,
        });
    }

    println!(
        "\n  {:<20} {:>12} {:>10} {:>10}",
        "Channel", "Q [mL/min]", "dP [Pa]", "tau_w [Pa]"
    );
    println!("  {}", "-".repeat(56));
    for r in &results_1d {
        println!(
            "  {:<20} {:>12.4} {:>10.2} {:>10.4}",
            r.id,
            r.flow_rate_m3_s * 60.0 * 1e6,
            r.pressure_drop_pa,
            r.wall_shear_pa,
        );
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Phase 3: 2D solve (cfd-2d)
    // ─────────────────────────────────────────────────────────────────────────

    println!("\n--- Phase 3: cfd-2d SIMPLE Network Solve (20x10 grid) ---\n");

    let nx = 20_usize;
    let ny = 10_usize;
    let blood_2d = BloodModel::Newtonian(BLOOD_VISCOSITY);
    let sink = Network2dBuilderSink::new(blood_2d, BLOOD_DENSITY, INLET_FLOW_M3_S, nx, ny);
    let mut net2d: Network2DSolver<f64> = sink.build(&blueprint)?;
    let result_2d = net2d.solve_all(1e-6)?;

    let converged = result_2d.converged_count;
    let total = result_2d.channels.len();

    println!("  Grid:            {nx} x {ny} per channel");
    println!("  Converged:       {converged}/{total}");
    println!("  Total HI:        {:.2e}", result_2d.total_hemolysis_index);
    println!(
        "  Max outlet err:  {:.2}%",
        result_2d.max_field_outlet_flow_error_pct
    );
    println!(
        "  Mean outlet err: {:.2}%",
        result_2d.mean_field_outlet_flow_error_pct
    );

    // Reference trace from the internal 1D solve
    let ref_trace = &result_2d.reference_trace;
    println!(
        "  Ref inlet P:     {:.2} Pa",
        ref_trace.scaled_inlet_pressure_pa
    );
    println!(
        "  Ref total Q_in:  {:.4} mL/min",
        ref_trace.total_inlet_flow_m3_s * 60.0 * 1e6
    );

    // ─────────────────────────────────────────────────────────────────────────
    // Phase 4: Cross-fidelity comparison table
    // ─────────────────────────────────────────────────────────────────────────

    println!("\n--- Phase 4: Cross-Fidelity Comparison ---\n");

    // Header
    println!(
        "  {:<20} {:>10} {:>10} {:>10} {:>10} {:>10} {:>8}",
        "Channel", "Q_1D", "Q_2D", "tau_1D", "tau_2D_max", "tau_2D_mu", "Q_err%"
    );
    println!(
        "  {:<20} {:>10} {:>10} {:>10} {:>10} {:>10} {:>8}",
        "", "[mL/min]", "[mL/min]", "[Pa]", "[Pa]", "[Pa]", ""
    );
    println!("  {}", "-".repeat(82));

    let mut total_q_1d = 0.0_f64;
    let mut total_q_2d_ref = 0.0_f64;
    let mut max_tau_1d = 0.0_f64;
    let mut max_tau_2d = 0.0_f64;

    for ch2d in &result_2d.channels {
        // Find matching 1D result
        let r1d = results_1d.iter().find(|r| r.id == ch2d.channel_id);
        let (q_1d, tau_1d) = r1d
            .map(|r| (r.flow_rate_m3_s, r.wall_shear_pa))
            .unwrap_or((0.0, 0.0));

        let q_2d_ref = ch2d.reference_trace.flow_rate_m3_s;
        let venturi_marker = if ch2d.is_venturi_throat { " [V]" } else { "" };

        println!(
            "  {:<20} {:>10.4} {:>10.4} {:>10.4} {:>10.4} {:>10.4} {:>8.2}",
            format!("{}{}", ch2d.channel_id, venturi_marker),
            q_1d * 60.0 * 1e6,
            q_2d_ref * 60.0 * 1e6,
            tau_1d,
            ch2d.field_wall_shear_max_pa,
            ch2d.field_wall_shear_mean_pa,
            ch2d.field_outlet_flow_error_pct,
        );

        total_q_1d += q_1d.abs();
        total_q_2d_ref += q_2d_ref.abs();
        max_tau_1d = max_tau_1d.max(tau_1d);
        max_tau_2d = max_tau_2d.max(ch2d.field_wall_shear_max_pa);
    }

    println!("  {}", "-".repeat(82));
    println!(
        "  {:<20} {:>10.4} {:>10.4} {:>10.4} {:>10.4}",
        "TOTALS / MAX",
        total_q_1d * 60.0 * 1e6,
        total_q_2d_ref * 60.0 * 1e6,
        max_tau_1d,
        max_tau_2d,
    );

    // ─────────────────────────────────────────────────────────────────────────
    // Phase 5: Venturi throat comparison
    // ─────────────────────────────────────────────────────────────────────────

    let venturi_2d: Vec<_> = result_2d
        .channels
        .iter()
        .filter(|c| c.is_venturi_throat)
        .collect();

    if !venturi_2d.is_empty() {
        println!("\n--- Phase 5: Venturi Throat Summary ---\n");
        println!(
            "  {:<16} {:>10} {:>12} {:>12} {:>10} {:>10}",
            "Throat", "Q [mL/min]", "tau_max [Pa]", "tau_mean [Pa]", "HI", "t_transit [s]"
        );
        println!("  {}", "-".repeat(74));

        for v in &venturi_2d {
            println!(
                "  {:<16} {:>10.4} {:>12.4} {:>12.4} {:>10.2e} {:>10.4e}",
                v.channel_id,
                v.reference_trace.flow_rate_m3_s * 60.0 * 1e6,
                v.field_wall_shear_max_pa,
                v.field_wall_shear_mean_pa,
                v.hemolysis_index,
                v.transit_time_s,
            );
        }

        println!(
            "\n  Total venturi HI:  {:.2e}",
            venturi_2d.iter().map(|v| v.hemolysis_index).sum::<f64>()
        );
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Phase 6: Write report
    // ─────────────────────────────────────────────────────────────────────────

    let report_path = out.join("cross_fidelity_report.md");
    let mut report = String::new();
    report.push_str("# Serpentine-Venturi 1D vs 2D Cross-Fidelity Report\n\n");
    report.push_str(&format!(
        "## Geometry\n\n\
         | Parameter | Value |\n\
         |-----------|-------|\n\
         | Segments | {segments} |\n\
         | Venturi throats | {n_venturis} |\n\
         | Channel width | {:.1} mm |\n\
         | Throat width | {:.1} mm |\n\
         | Height | {:.1} mm |\n\
         | Bend radius | {:.1} mm |\n\
         | D_h | {:.3} mm |\n\
         | Re | {:.1} |\n\
         | De | {:.1} |\n\
         | Contraction ratio | {:.1}:1 |\n\
         | Q_inlet | {:.2} mL/min |\n\n",
        width_m * 1e3,
        throat_width_m * 1e3,
        height_m * 1e3,
        bend_radius_m * 1e3,
        d_h * 1e3,
        re,
        de,
        width_m / throat_width_m,
        INLET_FLOW_M3_S * 60.0 * 1e6,
    ));

    report.push_str("## 1D Results (cfd-1d: Casson blood)\n\n");
    report.push_str(&format!("- Inlet pressure: {:.2} Pa\n", inlet_p));
    report.push_str(&format!("- Total dP: {:.2} Pa\n\n", inlet_p));

    report.push_str(&format!(
        "## 2D Results (cfd-2d: SIMPLE, {nx}x{ny} grid, Newtonian)\n\n"
    ));
    report.push_str(&format!("- Converged: {converged}/{total}\n"));
    report.push_str(&format!(
        "- Total hemolysis index: {:.2e}\n",
        result_2d.total_hemolysis_index
    ));
    report.push_str(&format!(
        "- Max outlet flow error: {:.2}%\n",
        result_2d.max_field_outlet_flow_error_pct
    ));
    report.push_str(&format!(
        "- Mean outlet flow error: {:.2}%\n\n",
        result_2d.mean_field_outlet_flow_error_pct
    ));

    report.push_str("## Cross-Fidelity Comparison\n\n");
    report.push_str("| Channel | Q_1D [mL/min] | Q_2D [mL/min] | tau_1D [Pa] | tau_2D_max [Pa] | tau_2D_mean [Pa] | Q_err% |\n");
    report.push_str("|---------|---------------|---------------|-------------|-----------------|------------------|--------|\n");
    for ch2d in &result_2d.channels {
        let r1d = results_1d.iter().find(|r| r.id == ch2d.channel_id);
        let (q_1d, tau_1d) = r1d
            .map(|r| (r.flow_rate_m3_s, r.wall_shear_pa))
            .unwrap_or((0.0, 0.0));
        let venturi_marker = if ch2d.is_venturi_throat { " [V]" } else { "" };
        report.push_str(&format!(
            "| {}{} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} | {:.2} |\n",
            ch2d.channel_id,
            venturi_marker,
            q_1d * 60.0 * 1e6,
            ch2d.reference_trace.flow_rate_m3_s * 60.0 * 1e6,
            tau_1d,
            ch2d.field_wall_shear_max_pa,
            ch2d.field_wall_shear_mean_pa,
            ch2d.field_outlet_flow_error_pct,
        ));
    }

    if !venturi_2d.is_empty() {
        report.push_str("\n## Venturi Throat Analysis\n\n");
        report.push_str(
            "| Throat | Q [mL/min] | tau_max [Pa] | tau_mean [Pa] | HI | Transit [s] |\n",
        );
        report
            .push_str("|--------|------------|--------------|---------------|----|-----------|\n");
        for v in &venturi_2d {
            report.push_str(&format!(
                "| {} | {:.4} | {:.4} | {:.4} | {:.2e} | {:.4e} |\n",
                v.channel_id,
                v.reference_trace.flow_rate_m3_s * 60.0 * 1e6,
                v.field_wall_shear_max_pa,
                v.field_wall_shear_mean_pa,
                v.hemolysis_index,
                v.transit_time_s,
            ));
        }
    }

    fs::write(&report_path, &report)?;
    println!("\nReport saved to: {}", report_path.display());
    println!("\nCross-fidelity comparison complete.");

    Ok(())
}
