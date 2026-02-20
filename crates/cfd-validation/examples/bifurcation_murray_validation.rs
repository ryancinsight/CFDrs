//! Bifurcation Validation Against Murray's Law
//!
//! ## Physical Model
//!
//! Murray's Law (1926) describes the optimal branching geometry in vascular networks
//! that minimizes metabolic cost (Womersley-type metabolic power ∝ vessel volume):
//!
//! ```text
//! D_parent³ = D_daughter1³ + D_daughter2³
//! ```
//!
//! For a symmetric bifurcation this simplifies to:
//!
//! ```text
//! D_daughter = D_parent / ∛2 ≈ 0.7937 * D_parent
//! ```
//!
//! ### Validation Criteria
//!
//! 1. **Murray's Law Diameter Ratio**: D_d / D_p ≈ 0.7937 (symmetric) — geometry check.
//! 2. **Flow Conservation**: Q_parent = Q_d1 + Q_d2 (continuity).
//! 3. **Pressure Continuity**: P at bifurcation node must be equal for both daughters.
//! 4. **Network Solver**: Results agree to within 0.1% of Hagen-Poiseuille analytical.
//!
//! ### Literature Reference
//!
//! - Murray, C. D. (1926). *The physiological principle of minimum work applied to
//!   the angle of branching of arteries.* J. General Physiology, 9, 835-841.
//! - Zamir, M. (1978). *Optimality principles in arterial branching.* J. Theoretical
//!   Biology, 62, 227-251. (Confirms cubic law for minimum power/minimum volume.)
//! - Cassab et al. MMFT comparison: bifurcation pressure correct to <1% vs 1D-lumped.

use cfd_1d::network::{EdgeProperties, Network, NetworkBuilder};
use cfd_1d::resistance::{HagenPoiseuilleModel, FlowConditions, ResistanceModel};
use cfd_1d::solver::{NetworkProblem, NetworkSolver, SolverConfig};
use cfd_core::compute::solver::Solver;
use cfd_core::physics::fluid::non_newtonian::CarreauYasuda;
use cfd_core::physics::fluid::FluidTrait;
use scheme::domain::model::{ChannelSpec, NodeKind, NodeSpec};

/// Validate Murray's Law compliance and flow conservation in a symmetric bifurcation.
fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("═══════════════════════════════════════════════════");
    println!("  CFD Validation: Bifurcation — Murray's Law Check");
    println!("═══════════════════════════════════════════════════");

    // ── Physical Parameters ──────────────────────────────────────────────────
    // Parent: D=4mm, L=50mm (microartery scale)
    // Murray's optimal daughter: D_d = D_p / 2^(1/3) = 3.175 mm
    // Using exact 3mm daughters (slight deviation for validation coverage)
    let d_parent = 4e-3_f64; // 4 mm  → m
    let d_d1     = 3e-3_f64; // 3 mm
    let d_d2     = 3e-3_f64; // 3 mm (symmetric)

    let l_parent = 50e-3_f64; // 50 mm
    let l_d1     = 40e-3_f64; // 40 mm
    let l_d2     = 40e-3_f64; // 40 mm

    let q_inlet  = 5e-6_f64; // 5 mL/s = 5e-6 m³/s
    let p_out    = 13_332.2_f64; // 100 mmHg in Pa

    // ── Blood Properties (Carreau-Yasuda) ────────────────────────────────────
    let blood = CarreauYasuda::<f64>::blood();
    // Apparent viscosity at physiological mid-shear (~100 s⁻¹): ~3.5 mPa·s
    let mu_apparent = 3.5e-3_f64; // Pa·s

    println!("Fluid: {} | density={:.0} kg/m³", blood.name(), 1060.0);
    println!("Apparent viscosity at 100 s⁻¹: {:.2} mPa·s", mu_apparent * 1e3);

    // ── Murray's Law Analytical Check ────────────────────────────────────────
    // Optimal: D_parent³ = Σ D_daughter³
    let murray_lhs = d_parent.powi(3);
    let murray_rhs = d_d1.powi(3) + d_d2.powi(3);
    let murray_deviation_pct = ((murray_lhs - murray_rhs) / murray_lhs).abs() * 100.0;

    // Optimal daughter diameter for Murray's Law
    let d_optimal = d_parent / (2.0_f64).cbrt();
    let murray_ratio_actual = d_d1 / d_parent;
    let murray_ratio_optimal = d_optimal / d_parent;

    println!("\n── Murray's Law Geometry Analysis ──────────────────");
    println!("  D_parent   = {:.2} mm", d_parent * 1e3);
    println!("  D_daughter = {:.2} mm (using)", d_d1 * 1e3);
    println!("  D_optimal  = {:.3} mm (Murray: D_p / ∛2)", d_optimal * 1e3);
    println!("  D_d/D_p actual:  {:.4}", murray_ratio_actual);
    println!("  D_d/D_p optimal: {:.4}", murray_ratio_optimal);
    println!("  Murray ΔP³ deviation: {:.2}%", murray_deviation_pct);

    // ── Hagen-Poiseuille Analytical Resistance ───────────────────────────────
    // R = 128 μ L / (π D⁴)
    let pi = std::f64::consts::PI;
    let r_parent_hp = 128.0 * mu_apparent * l_parent / (pi * d_parent.powi(4));
    let r_d1_hp     = 128.0 * mu_apparent * l_d1     / (pi * d_d1.powi(4));
    let r_d2_hp     = 128.0 * mu_apparent * l_d2     / (pi * d_d2.powi(4));

    // For symmetric split: Q_d = Q_parent / 2
    let q_d = q_inlet / 2.0;
    // Analytical pressure at bifurcation junction
    let dp_parent_analytic = r_parent_hp * q_inlet;
    let dp_d1_analytic     = r_d1_hp     * q_d;
    let p_inlet_analytic   = p_out + dp_d1_analytic + dp_parent_analytic;
    let p_bifurc_analytic  = p_out + dp_d1_analytic;

    println!("\n── Analytical (Hagen-Poiseuille) ────────────────────");
    println!("  R_parent = {:.4e} Pa·s/m³", r_parent_hp);
    println!("  R_d1     = {:.4e} Pa·s/m³", r_d1_hp);
    println!("  R_d2     = {:.4e} Pa·s/m³", r_d2_hp);
    println!("  P_inlet  = {:.2} Pa ({:.2} mmHg)", p_inlet_analytic, p_inlet_analytic / 133.322);
    println!("  P_bifurc = {:.2} Pa ({:.2} mmHg)", p_bifurc_analytic, p_bifurc_analytic / 133.322);

    // ── Network Solver ───────────────────────────────────────────────────────
    let nodes = vec![
        NodeSpec::new("Inlet",      NodeKind::Inlet),
        NodeSpec::new("Bifurcation", NodeKind::Junction),
        NodeSpec::new("Outlet1",    NodeKind::Outlet),
        NodeSpec::new("Outlet2",    NodeKind::Outlet),
    ];
    let channels = vec![
        ChannelSpec::new_pipe("Parent",    "Inlet",       "Bifurcation", l_parent, d_parent, r_parent_hp, 0.0),
        ChannelSpec::new_pipe("Daughter1", "Bifurcation", "Outlet1",     l_d1,     d_d1,     r_d1_hp,     0.0),
        ChannelSpec::new_pipe("Daughter2", "Bifurcation", "Outlet2",     l_d2,     d_d2,     r_d2_hp,     0.0),
    ];

    let mut builder = NetworkBuilder::<f64>::new();
    let mut node_idx = std::collections::HashMap::new();
    for ns in &nodes {
        let idx = match ns.kind {
            NodeKind::Inlet    => builder.add_inlet(ns.id.as_str().to_string()),
            NodeKind::Outlet   => builder.add_outlet(ns.id.as_str().to_string()),
            NodeKind::Junction | NodeKind::Reservoir => {
                builder.add_junction(ns.id.as_str().to_string())
            }
        };
        node_idx.insert(ns.id.as_str().to_string(), idx);
    }

    let mut edge_idx = std::collections::HashMap::new();
    for ch in &channels {
        let from = node_idx[ch.from.as_str()];
        let to   = node_idx[ch.to.as_str()];
        let eidx = builder.connect_with_pipe(from, to, ch.id.as_str().to_string());
        edge_idx.insert(ch.id.as_str().to_string(), eidx);
    }

    let graph   = builder.build()?;
    let mut net = Network::new(graph, blood.clone());

    for ch in &channels {
        net.add_edge_properties(edge_idx[ch.id.as_str()], EdgeProperties::from(ch));
    }

    let inlet_idx   = node_idx["Inlet"];
    let outlet1_idx = node_idx["Outlet1"];
    let outlet2_idx = node_idx["Outlet2"];

    net.set_neumann_flow(inlet_idx, q_inlet);
    net.set_pressure(outlet1_idx, p_out);
    net.set_pressure(outlet2_idx, p_out);

    // Initial guesses for non-Newtonian iterative solver
    net.set_flow_rate(edge_idx["Parent"],    q_inlet);
    net.set_flow_rate(edge_idx["Daughter1"], q_d);
    net.set_flow_rate(edge_idx["Daughter2"], q_d);
    net.update_resistances()?;

    let config   = SolverConfig { tolerance: 1e-9, max_iterations: 200 };
    let solver   = NetworkSolver::<f64, CarreauYasuda<f64>>::with_config(config);
    let solution = solver.solve(&NetworkProblem::new(net))?;

    // ── Extract Results ──────────────────────────────────────────────────────
    let p_inlet  = *solution.pressures.get(&node_idx["Inlet"]).unwrap_or(&0.0);
    let p_bifurc = *solution.pressures.get(&node_idx["Bifurcation"]).unwrap_or(&0.0);
    let p_out1   = *solution.pressures.get(&node_idx["Outlet1"]).unwrap_or(&0.0);
    let p_out2   = *solution.pressures.get(&node_idx["Outlet2"]).unwrap_or(&0.0);

    let q_parent = *solution.flow_rates.get(&edge_idx["Parent"]).unwrap_or(&0.0);
    let q_d1     = *solution.flow_rates.get(&edge_idx["Daughter1"]).unwrap_or(&0.0);
    let q_d2     = *solution.flow_rates.get(&edge_idx["Daughter2"]).unwrap_or(&0.0);

    println!("\n── NetworkSolver Results ────────────────────────────");
    println!("  P_inlet  = {:.2} Pa ({:.2} mmHg)", p_inlet, p_inlet / 133.322);
    println!("  P_bifurc = {:.2} Pa ({:.2} mmHg)", p_bifurc, p_bifurc / 133.322);
    println!("  P_out1   = {:.2} Pa ({:.2} mmHg)", p_out1, p_out1 / 133.322);
    println!("  P_out2   = {:.2} Pa ({:.2} mmHg)", p_out2, p_out2 / 133.322);
    println!("  Q_parent    = {:.4} mL/s", q_parent * 1e6);
    println!("  Q_daughter1 = {:.4} mL/s", q_d1.abs() * 1e6);
    println!("  Q_daughter2 = {:.4} mL/s", q_d2.abs() * 1e6);

    // ── Validation Assertions ────────────────────────────────────────────────
    let mut all_pass = true;
    println!("\n── Validation Results ───────────────────────────────");

    // 1) Flow conservation: Q_parent = Q_d1 + Q_d2
    let q_out_total = q_d1.abs() + q_d2.abs();
    let flow_conservation_err = (q_parent.abs() - q_out_total).abs() / q_parent.abs();
    let fc_pass = flow_conservation_err < 1e-6;
    println!(
        "  [{}] Flow Conservation: |Q_p - (Q_d1+Q_d2)| = {:.2e} ({:.4}%)",
        if fc_pass { "PASS" } else { "FAIL" },
        (q_parent.abs() - q_out_total).abs(),
        flow_conservation_err * 100.0
    );
    all_pass &= fc_pass;

    // 2) Pressure at solver vs Hagen-Poiseuille analytical
    let p_inlet_err  = (p_inlet - p_inlet_analytic).abs() / p_inlet_analytic;
    let p_bifurc_err = (p_bifurc - p_bifurc_analytic).abs() / p_bifurc_analytic;
    let p_inlet_pass  = p_inlet_err  < 0.001; // < 0.1%
    let p_bifurc_pass = p_bifurc_err < 0.001;
    println!(
        "  [{}] P_inlet vs Analytical:  {:.4}% deviation",
        if p_inlet_pass  { "PASS" } else { "FAIL" }, p_inlet_err  * 100.0
    );
    println!(
        "  [{}] P_bifurc vs Analytical: {:.4}% deviation",
        if p_bifurc_pass { "PASS" } else { "FAIL" }, p_bifurc_err * 100.0
    );
    all_pass &= p_inlet_pass && p_bifurc_pass;

    // 3) Pressure continuity: both outlets at p_out ≈ 13332.2 Pa
    let p_out1_err = (p_out1 - p_out).abs() / p_out;
    let p_out2_err = (p_out2 - p_out).abs() / p_out;
    let bc_pass = p_out1_err < 1e-9 && p_out2_err < 1e-9;
    println!(
        "  [{}] Outlet Pressure BCs: P_out1={:.2} Pa, P_out2={:.2} Pa",
        if bc_pass { "PASS" } else { "FAIL" }, p_out1, p_out2
    );
    all_pass &= bc_pass;

    // 4) Murray's Law geometry deviation < 5% (actual vs optimal D_d/D_p)
    let murray_pass = murray_deviation_pct < 5.0;
    println!(
        "  [{}] Murray's Law deviation: {:.2}% (actual D³ sum vs parent D³)",
        if murray_pass { "PASS" } else { "FAIL" },
        murray_deviation_pct
    );
    all_pass &= murray_pass;

    // 5) Symmetric flow split: Q_d1 ≈ Q_d2
    let sym_err = (q_d1.abs() - q_d2.abs()).abs() / q_d1.abs();
    let sym_pass = sym_err < 1e-6;
    println!(
        "  [{}] Symmetric split: Q_d1/Q_d2 = {:.6}",
        if sym_pass { "PASS" } else { "FAIL" },
        q_d1.abs() / q_d2.abs()
    );
    all_pass &= sym_pass;

    println!();
    if all_pass {
        println!("✅ ALL VALIDATIONS PASSED — Bifurcation simulation matches Murray's Law and HP analytical.");
    } else {
        eprintln!("❌ VALIDATION FAILED — See above for failing criteria.");
        std::process::exit(1);
    }

    Ok(())
}
