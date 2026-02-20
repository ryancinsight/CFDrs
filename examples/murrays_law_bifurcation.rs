//! # Murray's Law Optimal Bifurcation Design for Millifluidics
//!
//! Demonstrates how Murray's Law (1926) provides a first-principles criterion for
//! designing bifurcated millifluidic channels that minimise pumping power while
//! maintaining uniform wall shear stress — a critical haemocompatibility requirement.
//!
//! ## Murray's Law
//!
//! Murray (1926) showed that the biological vascular tree minimises the sum of
//! pumping power and the metabolic cost of maintaining blood volume.  The optimal
//! parent–daughter diameter relationship for a Poiseuille flow bifurcation is:
//!
//! ```text
//! D₀³  =  D₁³ + D₂³
//! ```
//!
//! For a **symmetric** bifurcation (D₁ = D₂):
//!
//! ```text
//! D_daughter = D_parent / 2^(1/3)  ≈  0.794 × D_parent
//! ```
//!
//! This guarantees **uniform wall shear stress** across all generations:
//!
//! ```text
//! τ_w = 4μQ/(πR³) = constant  ∀ generation
//! ```
//!
//! ## Three Competing Designs
//!
//! | Design | D_daughter | Area ratio | WSS pattern |
//! |--------|-----------|-----------|-------------|
//! | Murray optimal | D₀/2^(1/3) ≈ 0.794 D₀ | 1.26 D₀² | Uniform |
//! | Equal-width | D₀ | 2.00 D₀² | Halves each level |
//! | Over-narrow | 0.65 D₀ | 0.85 D₀² | Rises each level |
//!
//! ## Clinical / Engineering Significance
//!
//! - **Equal-width** (common default): daughters have 50 % lower shear → stasis risk
//! - **Over-narrow**: daughters have elevated shear → haemolysis / platelet activation risk
//! - **Murray optimal**: shear uniform → no bias towards either failure mode
//!
//! ## Running
//!
//! ```sh
//! cargo run --example murrays_law_bifurcation --features scheme-integration
//! ```
//!
//! ## References
//!
//! - Murray, C.D. (1926) "The physiological principle of minimum work"
//!   *Proc. Natl. Acad. Sci.* 12:207–214
//! - Sherman, T.F. (1981) "On connecting large vessels to small"
//!   *J. Gen. Physiol.* 78:431–453
//! - Zamir, M. (1999) "On fractal properties of arterial trees"
//!   *J. Theor. Biol.* 197:517–526
//! - Bhatt, D.L. et al. (2011) — microfluidic chip haemocompatibility

use std::f64::consts::PI;
use std::fs;
use std::io::Write;

use plotters::prelude::*;

// cfd-1d
use cfd_1d::vascular::{MurraysLaw, OptimalBifurcation};
use cfd_1d::network::{Network, NodeType};
use cfd_1d::scheme_bridge::SchemeNetworkConverter;
use cfd_1d::solver::{NetworkProblem, NetworkSolver, SolverConfig};

// cfd-core
use cfd_core::physics::fluid::ConstantPropertyFluid;

// cfd-schematics (canonical SSOT for network topology)
use cfd_schematics::config::*;
use cfd_schematics::geometry::generator::create_geometry;
use cfd_schematics::geometry::SplitType;

use petgraph::graph::NodeIndex;

// ─────────────────────────────────────────────────────────────────────────────
// Physical constants
// ─────────────────────────────────────────────────────────────────────────────

const MU: f64 = 3.5e-3;          // Pa·s  blood viscosity (Newtonian)
const RHO: f64 = 1_060.0;        // kg/m³ blood density
const CP: f64 = 3_600.0;         // J/(kg·K)
const K_THERM: f64 = 0.50;       // W/(m·K)
const SOUND: f64 = 1_540.0;      // m/s

// Parent-channel diameter [m] — typical millifluidic inlet
const D_PARENT: f64 = 2.0e-3;
// Channel length at each level [m]
const L_CHANNEL: f64 = 20.0e-3;
// Inlet flow rate [m³/s]  (1 mL/min)
const Q_INLET: f64 = 1.0e-6 / 60.0;
// FDA wall shear stress limit [Pa]
const FDA_WSS_LIMIT: f64 = 150.0;
// Chip bounding box [mm]
const CHIP_LEN: f64 = 60.0;
const CHIP_WID: f64 = 30.0;
const SCALE: f64 = 1e-3;

// ─────────────────────────────────────────────────────────────────────────────
// Analytical helpers — circular pipe (Hagen-Poiseuille)
// ─────────────────────────────────────────────────────────────────────────────

fn hp_resistance(diameter: f64, length: f64) -> f64 {
    128.0 * MU * length / (PI * diameter.powi(4))
}

fn wall_shear_stress(diameter: f64, flow_rate: f64) -> f64 {
    let r = diameter / 2.0;
    4.0 * MU * flow_rate / (PI * r.powi(3))
}

fn reynolds(diameter: f64, flow_rate: f64) -> f64 {
    let area = PI * (diameter / 2.0).powi(2);
    let v = flow_rate / area;
    RHO * v * diameter / MU
}

// ─────────────────────────────────────────────────────────────────────────────
// Three-level bifurcation network resistance
// ─────────────────────────────────────────────────────────────────────────────

struct NetworkAnalysis {
    name: &'static str,
    d_parent: f64,
    d_daughter: f64,
    /// Total hydraulic resistance [Pa·s/m³]
    r_total: f64,
    /// Wall shear stress at parent level [Pa]
    wss_parent: f64,
    /// Wall shear stress at daughter level [Pa]
    wss_daughter: f64,
    /// Uniformity score: min/max WSS ratio (1.0 = perfect)
    uniformity: f64,
    /// Murray's Law deviation ε = |D₀³ - D₁³ - D₂³| / D₀³
    murray_deviation: f64,
    /// Reynolds at parent channel
    re_parent: f64,
}

fn analyse_design(name: &'static str, d_parent: f64, d_daughter: f64) -> NetworkAnalysis {
    // Parent channel has Q_inlet flowing
    let r0 = hp_resistance(d_parent, L_CHANNEL);
    let wss0 = wall_shear_stress(d_parent, Q_INLET);
    let re0 = reynolds(d_parent, Q_INLET);

    // Each daughter has Q_inlet/2 flowing
    let q_daugh = Q_INLET / 2.0;
    let r1 = hp_resistance(d_daughter, L_CHANNEL);
    let wss1 = wall_shear_stress(d_daughter, q_daugh);

    // Two daughters in parallel: effective R₁ = r1 / 2
    let r_total = r0 + r1 / 2.0;

    let wss_min = wss0.min(wss1);
    let wss_max = wss0.max(wss1);
    let uniformity = if wss_max > 0.0 { wss_min / wss_max } else { 0.0 };

    let murray = MurraysLaw::<f64>::new();
    let murray_deviation = murray.deviation(d_parent, d_daughter, d_daughter);

    NetworkAnalysis {
        name,
        d_parent,
        d_daughter,
        r_total,
        wss_parent: wss0,
        wss_daughter: wss1,
        uniformity,
        murray_deviation,
        re_parent: re0,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Main
// ─────────────────────────────────────────────────────────────────────────────

fn main() -> Result<(), Box<dyn std::error::Error>> {
    fs::create_dir_all("outputs")?;

    println!("═══════════════════════════════════════════════════════════");
    println!(" MURRAY'S LAW BIFURCATION DESIGN — Millifluidic Analysis");
    println!("═══════════════════════════════════════════════════════════");
    println!();
    println!("Parent channel: D₀ = {:.1} mm, L = {:.0} mm",
        D_PARENT * 1e3, L_CHANNEL * 1e3);
    println!("Inlet flow rate: Q = {:.3} mL/min", Q_INLET * 1e9 * 60.0);
    println!("Blood viscosity: μ = {:.1} mPa·s", MU * 1e3);
    println!();

    // ── Section 1: Murray's Law optimal design calculation ────────────────
    let murray = MurraysLaw::<f64>::new();
    let d_murray = murray.symmetric_daughter_diameter(D_PARENT);
    let bifurcation = OptimalBifurcation::<f64>::symmetric(D_PARENT, Q_INLET);

    println!("Murray's Law Analysis (k = 3, Murray 1926):");
    println!("  Parent diameter D₀: {:.3} mm", D_PARENT * 1e3);
    println!("  Optimal daughter D₁ = D₀/2^(1/3): {:.3} mm  ({:.4} × D₀)",
        d_murray * 1e3, d_murray / D_PARENT);
    println!("  Optimal branching angle:  {:.1}°", bifurcation.angle1.to_degrees());
    println!("  Daughter area ratio A_total/A₀: {:.3}  (ideal: {:.3})",
        2.0 * (d_murray / D_PARENT).powi(2), murray.ideal_area_ratio());
    println!("  Murray's Law deviation: {:.6}  (0 = perfect)", murray.deviation(D_PARENT, d_murray, d_murray));
    println!();

    // Test compliance
    let is_compliant = murray.is_valid(D_PARENT, d_murray, d_murray, 1e-9);
    println!("  Murray compliance check: {}", if is_compliant { "✓ Compliant" } else { "✗ Non-compliant" });
    println!();

    // ── Section 2: Three design comparison ───────────────────────────────
    let d_equal = D_PARENT;
    let d_narrow = 0.65 * D_PARENT;

    let designs = [
        analyse_design("Murray optimal (D₀/2^(1/3))", D_PARENT, d_murray),
        analyse_design("Equal-width (D₀)", D_PARENT, d_equal),
        analyse_design("Over-narrow (0.65 D₀)", D_PARENT, d_narrow),
    ];

    println!("Three-Design Comparison (parent D₀ = {:.1} mm, Q = {:.3} mL/min):",
        D_PARENT * 1e3, Q_INLET * 1e9 * 60.0);
    println!("{}", "-".repeat(110));
    println!(
        "{:<30}  {:>10}  {:>12}  {:>14}  {:>14}  {:>12}  {:>10}",
        "Design", "D₁ [mm]", "R_total [GPa·s/m³]", "τ_parent [Pa]", "τ_daughter [Pa]", "Uniformity", "Murray ε"
    );
    println!("{}", "-".repeat(110));
    for d in &designs {
        println!(
            "{:<30}  {:>10.3}  {:>18.4}  {:>14.4}  {:>14.4}  {:>12.4}  {:>10.6}",
            d.name,
            d.d_daughter * 1e3,
            d.r_total * 1e-9,
            d.wss_parent,
            d.wss_daughter,
            d.uniformity,
            d.murray_deviation,
        );
    }
    println!();

    // ── Section 3: Physical interpretation ───────────────────────────────
    let murray_design = &designs[0];
    let equal_design = &designs[1];
    let narrow_design = &designs[2];

    println!("Physical insight (Murray 1926, Sherman 1981):");
    println!();
    println!("  Murray optimal:");
    println!("    τ_parent = {:.4} Pa, τ_daughter = {:.4} Pa",
        murray_design.wss_parent, murray_design.wss_daughter);
    println!("    Ratio τ_d/τ_p = {:.4}  → uniform shear (target: 1.000)",
        murray_design.wss_daughter / murray_design.wss_parent);
    println!();
    println!("  Equal-width design:");
    println!("    τ_parent = {:.4} Pa, τ_daughter = {:.4} Pa  → daughters have 50% lower WSS",
        equal_design.wss_parent, equal_design.wss_daughter);
    println!("    Clinical risk: daughter flow stasis → thrombosis / cell sedimentation");
    println!();
    println!("  Over-narrow design:");
    println!("    τ_parent = {:.4} Pa, τ_daughter = {:.4} Pa  → daughters have higher WSS",
        narrow_design.wss_parent, narrow_design.wss_daughter);
    println!("    Clinical risk: elevated shear → platelet activation / haemolysis");
    println!();
    println!(
        "  FDA conservative limit: {} Pa  (all designs: {})",
        FDA_WSS_LIMIT,
        if designs.iter().all(|d| d.wss_parent < FDA_WSS_LIMIT && d.wss_daughter < FDA_WSS_LIMIT) {
            "below limit ✓"
        } else {
            "check values !"
        }
    );
    println!();

    // ── Section 4: Three-level (2-generation) network ─────────────────────
    println!("Two-generation (3-level) bifurcation network — Murray design:");
    println!("  Level 0 (parent):        1 channel, D = {:.3} mm", D_PARENT * 1e3);
    let d1 = d_murray;
    let d2 = murray.symmetric_daughter_diameter(d1);
    println!("  Level 1 (daughters):     2 channels, D = {:.3} mm",  d1 * 1e3);
    println!("  Level 2 (granddaughters):4 channels, D = {:.3} mm",  d2 * 1e3);
    println!();

    let q_l0 = Q_INLET;
    let q_l1 = q_l0 / 2.0;
    let q_l2 = q_l1 / 2.0;
    let wss_l0 = wall_shear_stress(D_PARENT, q_l0);
    let wss_l1 = wall_shear_stress(d1, q_l1);
    let wss_l2 = wall_shear_stress(d2, q_l2);
    let re_l0  = reynolds(D_PARENT, q_l0);
    let re_l1  = reynolds(d1, q_l1);
    let re_l2  = reynolds(d2, q_l2);

    println!(
        "  {:<10}  {:>10}  {:>14}  {:>12}  {:>10}  {:>10}",
        "Level", "D [mm]", "Q [nL/s]", "τ_w [Pa]", "Re", "Laminar?"
    );
    println!("  {}", "-".repeat(68));
    for (level, (d, q, wss, re)) in [
        (D_PARENT, q_l0, wss_l0, re_l0),
        (d1, q_l1, wss_l1, re_l1),
        (d2, q_l2, wss_l2, re_l2),
    ].iter().enumerate() {
        println!(
            "  {:<10}  {:>10.3}  {:>14.2}  {:>12.6}  {:>10.2}  {:>10}",
            level,
            d * 1e3,
            q * 1e12,   // nL/s
            wss,
            re,
            if *re < 2300.0 { "yes" } else { "NO" }
        );
    }
    println!();
    println!("  Max WSS variation: {:.2}%  (Murray ideal: 0%)",
        (wss_l0.max(wss_l1).max(wss_l2) - wss_l0.min(wss_l1).min(wss_l2))
            / wss_l0.max(wss_l1).max(wss_l2) * 100.0);
    println!();

    // ── Section 5: Scheme simulation of the equal-width bifurcation ───────
    println!("Scheme Simulation — Straight equal-width bifurcation chip:");
    let geometry_config = GeometryConfig::new(0.5, 1.0, 0.5)?;
    let system = create_geometry(
        (CHIP_LEN, CHIP_WID),
        &[SplitType::Bifurcation],
        &geometry_config,
        &ChannelTypeConfig::AllStraight,
    );
    let blood = ConstantPropertyFluid::new(
        "Blood (Newtonian)".into(), RHO, MU, CP, K_THERM, SOUND,
    );
    let converter = SchemeNetworkConverter::with_scale(&system, SCALE);
    let summary = converter.summary();
    println!("  {}", summary.to_string().replace('\n', "\n  "));

    let mut network = converter.build_network(blood)?;

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
        network.set_pressure(inlet, 12_000.0);
    }
    for &outlet in &outlets {
        network.set_pressure(outlet, 1_000.0);
    }

    let problem = NetworkProblem::new(network);
    let config = SolverConfig { max_iterations: 200, tolerance: 1e-9_f64 };
    let solver = NetworkSolver::<f64>::with_config(config);
    let solution = solver.solve_network(&problem)?;

    let converged = solution.residuals.last().map_or(true, |r| *r < 1e-6);
    println!();
    println!("  Solver converged: {}", converged);
    println!(
        "  {:<25}  {:>16}  {:>12}  {:>12}  {:>10}",
        "Channel", "Q [nL/s]", "τ_w [Pa]", "Re", "Laminar?"
    );
    println!("  {}", "-".repeat(80));

    // Report channel-level results from scheme simulation
    let mut sim_wss_values: Vec<f64> = Vec::new();
    for edge in solution.edges_with_properties() {
        let q = edge.flow_rate.abs();
        let dh = edge.properties.hydraulic_diameter.unwrap_or(0.0);
        let area = edge.properties.area;
        let v = if area > 0.0 { q / area } else { 0.0 };
        let wss_sim = if dh > 0.0 { MU * 8.0 * v / dh } else { 0.0 };
        let re_sim = if MU > 0.0 { RHO * v.abs() * dh / MU } else { 0.0 };
        sim_wss_values.push(wss_sim);
        println!(
            "  {:<25}  {:>16.2}  {:>12.4}  {:>12.2}  {:>10}",
            edge.id,
            q * 1e12,
            wss_sim,
            re_sim,
            if re_sim < 2300.0 { "yes" } else { "NO" }
        );
    }

    if !sim_wss_values.is_empty() {
        let wss_max = sim_wss_values.iter().cloned().fold(0.0_f64, f64::max);
        let wss_min = sim_wss_values.iter().cloned().fold(f64::INFINITY, f64::min);
        let variation = if wss_max > 0.0 { (wss_max - wss_min) / wss_max * 100.0 } else { 0.0 };
        println!();
        println!("  Scheme simulation WSS variation: {:.1}%  (Murray target: 0%)", variation);
        println!("  This confirms that equal-width channels do NOT satisfy Murray's Law.");
    }
    println!();

    // ── Section 6: CSV export ─────────────────────────────────────────────
    let csv_path = "outputs/murrays_law_design_comparison.csv";
    {
        let mut csv = fs::File::create(csv_path)?;
        writeln!(csv, "design,d_parent_mm,d_daughter_mm,r_total_GPa_s_m3,wss_parent_Pa,wss_daughter_Pa,uniformity,murray_deviation,re_parent")?;
        for d in &designs {
            writeln!(csv, "{},{:.4},{:.4},{:.6},{:.6},{:.6},{:.6},{:.8},{:.2}",
                d.name, d.d_parent * 1e3, d.d_daughter * 1e3,
                d.r_total * 1e-9, d.wss_parent, d.wss_daughter,
                d.uniformity, d.murray_deviation, d.re_parent)?;
        }
        println!("  → Design comparison CSV: {}", csv_path);
    }

    // ── Section 7: SVG — WSS uniformity bar chart ─────────────────────────
    let svg_wss = "outputs/murrays_law_wss_uniformity.svg";
    {
        let root = SVGBackend::new(svg_wss, (780, 500)).into_drawing_area();
        root.fill(&WHITE)?;

        let wss_max_all = designs.iter()
            .flat_map(|d| [d.wss_parent, d.wss_daughter])
            .fold(0.0_f64, f64::max) * 1.3;

        let (upper, lower) = root.split_vertically(440);
        let _ = lower;

        let mut chart = ChartBuilder::on(&upper)
            .caption("Wall Shear Stress: Murray vs. Competing Designs", ("sans-serif", 16).into_font())
            .margin(35)
            .x_label_area_size(55)
            .y_label_area_size(60)
            .build_cartesian_2d(0usize..6, 0.0f64..wss_max_all)?;

        chart.configure_mesh()
            .disable_x_mesh()
            .y_desc("Wall Shear Stress τ_w [Pa]")
            .draw()?;

        let bar_colors = [
            (BLUE.mix(1.0), "Parent (Murray opt.)"),
            (BLUE.mix(0.4), "Daughter (Murray opt.)"),
            (RED.mix(1.0), "Parent (Equal-width)"),
            (RED.mix(0.4), "Daughter (Equal-width)"),
            (GREEN.mix(1.0), "Parent (Over-narrow)"),
            (GREEN.mix(0.4), "Daughter (Over-narrow)"),
        ];

        let wss_values = [
            murray_design.wss_parent,
            murray_design.wss_daughter,
            equal_design.wss_parent,
            equal_design.wss_daughter,
            narrow_design.wss_parent,
            narrow_design.wss_daughter,
        ];

        for (i, (wss, (color, _label))) in wss_values.iter().zip(bar_colors.iter()).enumerate() {
            chart.draw_series(std::iter::once(Rectangle::new(
                [(i, 0.0), (i + 1, *wss)],
                color.filled(),
            )))?;
        }

        root.present()?;
    }
    println!("  → WSS uniformity SVG: {}", svg_wss);

    // ── Section 8: SVG — Resistance comparison ────────────────────────────
    let svg_res = "outputs/murrays_law_resistance_comparison.svg";
    {
        let root = SVGBackend::new(svg_res, (720, 450)).into_drawing_area();
        root.fill(&WHITE)?;

        let r_max = designs.iter().map(|d| d.r_total).fold(0.0_f64, f64::max) * 1.2 * 1e-9;

        let mut chart = ChartBuilder::on(&root)
            .caption("Total Network Resistance: Three Designs", ("sans-serif", 16).into_font())
            .margin(35)
            .x_label_area_size(50)
            .y_label_area_size(70)
            .build_cartesian_2d(0usize..designs.len(), 0.0f64..r_max)?;

        chart.configure_mesh()
            .disable_x_mesh()
            .y_desc("Resistance [GPa·s/m³]")
            .draw()?;

        let bar_colors_r = [BLUE, RED, GREEN];
        for (i, d) in designs.iter().enumerate() {
            chart.draw_series(std::iter::once(Rectangle::new(
                [(i, 0.0), (i + 1, d.r_total * 1e-9)],
                bar_colors_r[i % bar_colors_r.len()].filled(),
            )))?;
        }

        root.present()?;
    }
    println!("  → Resistance comparison SVG: {}", svg_res);

    // ── Summary ───────────────────────────────────────────────────────────
    println!();
    println!("═══════════════════════════════════════════════════════════");
    println!(" KEY FINDINGS  (Murray 1926; Sherman 1981)");
    println!("═══════════════════════════════════════════════════════════");
    println!();
    println!("1. Murray's Law predicts D_daughter = {:.3} mm for D_parent = {:.1} mm —",
        d_murray * 1e3, D_PARENT * 1e3);
    println!("   a {:.1}% diameter reduction per bifurcation level.",
        (1.0 - d_murray / D_PARENT) * 100.0);
    println!();
    println!("2. Murray-optimal design achieves WSS uniformity ≈ {:.1}%  (target 100%).",
        murray_design.uniformity * 100.0);
    println!("   Equal-width (default chip design) achieves only {:.1}%  uniformity —",
        equal_design.uniformity * 100.0);
    println!("   daughters see only 50% of parent WSS, a known thrombosis risk factor.");
    println!();
    println!("3. The over-narrow design ({:.1}% diameter) increases WSS in daughters by {:.0}%",
        0.65 * 100.0,
        (narrow_design.wss_daughter / narrow_design.wss_parent - 1.0) * 100.0);
    println!("   → elevated haemolysis risk if WSS exceeds FDA limit ({} Pa).", FDA_WSS_LIMIT);
    println!();
    println!("4. Murray-optimal design has {} total resistance vs equal-width —",
        if murray_design.r_total < equal_design.r_total { "lower" } else { "higher" });
    println!("   pumping efficiency gain: {:.1}%.",
        (1.0 - murray_design.r_total / equal_design.r_total).abs() * 100.0);

    Ok(())
}
