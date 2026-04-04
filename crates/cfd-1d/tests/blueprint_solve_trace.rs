//! Integration tests: cfd-schematics blueprint → cfd-1d solve → per-node/per-channel trace.
//!
//! These tests validate that:
//! - `network_from_blueprint` correctly converts schematic geometry into a solver-ready network.
//! - `NetworkSolver` produces Kirchhoff-consistent pressures and flow rates.
//! - Mass is conserved at every interior node to machine precision.
//! - Pressure drops across channels are consistent with their computed resistances.
//! - Symmetric topologies produce symmetric flow distributions.

use cfd_1d::{
    domain::network::network_from_blueprint, Network, NetworkProblem, NetworkSolver, SolverConfig,
};
use cfd_core::physics::fluid::database::water_20c;
use cfd_schematics::domain::model::{
    ChannelShape, ChannelSpec, CrossSectionSpec, EdgeId, EdgeKind, NetworkBlueprint, NodeId,
    NodeKind, NodeSpec,
};

type Water = cfd_core::physics::fluid::newtonian::ConstantPropertyFluid<f64>;

// ── Blueprint construction helpers ───────────────────────────────────────────

fn node(id: &str, kind: NodeKind, x: f64, y: f64) -> NodeSpec {
    NodeSpec {
        id: NodeId::new(id),
        kind,
        point: (x, y),
        layout: None,
        junction_geometry: None,
        metadata: None,
    }
}

fn pipe(id: &str, from: &str, to: &str, length_m: f64, diameter_m: f64) -> ChannelSpec {
    ChannelSpec {
        id: EdgeId::new(id),
        kind: EdgeKind::Pipe,
        from: NodeId::new(from),
        to: NodeId::new(to),
        length_m,
        resistance: 0.0, // computed from geometry by network_from_blueprint
        cross_section: CrossSectionSpec::Circular { diameter_m },
        channel_shape: ChannelShape::Straight,
        quad_coeff: 0.0,
        valve_cv: None,
        pump_max_flow: None,
        pump_max_pressure: None,
        path: vec![],
        visual_role: None,
        therapy_zone: None,
        venturi_geometry: None,
        metadata: None,
    }
}

/// Solve a blueprint with Dirichlet pressure at every boundary node.
/// Returns the solved network.
fn solve_dirichlet(
    blueprint: &NetworkBlueprint,
    boundary_pressures: &[(&str, f64)],
) -> Network<f64, Water> {
    let fluid = water_20c::<f64>().expect("water_20c must succeed");
    let mut network = network_from_blueprint(blueprint, fluid).expect("blueprint must be valid");

    // Map node id string → NodeIndex
    let node_ids: Vec<(petgraph::graph::NodeIndex, String)> = network
        .graph
        .node_indices()
        .filter_map(|idx| network.graph.node_weight(idx).map(|n| (idx, n.id.clone())))
        .collect();

    for (id_str, pressure) in boundary_pressures {
        let idx = node_ids
            .iter()
            .find(|(_, id)| id == id_str)
            .map(|(idx, _)| *idx)
            .unwrap_or_else(|| panic!("boundary node '{}' not found in network", id_str));
        network.set_pressure(idx, *pressure);
    }

    let solver = NetworkSolver::<f64, Water>::with_config(SolverConfig {
        tolerance: 1e-10,
        max_iterations: 500,
    });
    let (solved, _diag) = solver
        .solve_network_with_diagnostics(&NetworkProblem::new(network))
        .unwrap_or_else(|e| panic!("solve failed: {e}"));
    solved
}

/// Return the solved pressure at a named node.
fn node_pressure(network: &Network<f64, Water>, id: &str) -> f64 {
    let idx = network
        .graph
        .node_indices()
        .find(|&i| network.graph.node_weight(i).map_or(false, |n| n.id == id))
        .unwrap_or_else(|| panic!("node '{}' not found", id));
    *network
        .pressures()
        .get(idx.index())
        .unwrap_or_else(|| panic!("no pressure for node '{}'", id))
}

/// Return the solved flow rate on a named edge.
fn edge_flow(network: &Network<f64, Water>, id: &str) -> f64 {
    let eidx = network
        .graph
        .edge_indices()
        .find(|&i| network.graph.edge_weight(i).map_or(false, |e| e.id == id))
        .unwrap_or_else(|| panic!("edge '{}' not found", id));
    *network
        .flow_rates()
        .get(eidx.index())
        .unwrap_or_else(|| panic!("no flow rate for edge '{}'", id))
}

/// Return the resistance of a named edge as stored in the graph after network_from_blueprint.
fn edge_resistance(network: &Network<f64, Water>, id: &str) -> f64 {
    let eidx = network
        .graph
        .edge_indices()
        .find(|&i| network.graph.edge_weight(i).map_or(false, |e| e.id == id))
        .unwrap_or_else(|| panic!("edge '{}' not found", id));
    network
        .graph
        .edge_weight(eidx)
        .map(|e| e.resistance)
        .unwrap_or_else(|| panic!("no weight for edge '{}'", id))
}

// ── Test 1: Series circuit ────────────────────────────────────────────────────

/// A single circular pipe with Dirichlet pressure at both ends.
/// Kirchhoff: Q = ΔP / R.
#[test]
fn series_circuit_pressure_drop_matches_kirchhoff() {
    const D_M: f64 = 1e-3; // 1 mm
    const L_M: f64 = 0.1; // 100 mm
    const P_IN: f64 = 1000.0; // Pa
    const P_OUT: f64 = 0.0;

    let mut bp = NetworkBlueprint::new_with_explicit_positions("series");
    bp.add_node(node("inlet", NodeKind::Inlet, 0.0, 0.0));
    bp.add_node(node("outlet", NodeKind::Outlet, 100.0, 0.0));
    bp.add_channel(
        pipe("pipe1", "inlet", "outlet", L_M, D_M).with_path(vec![(0.0, 0.0), (100.0, 0.0)]),
    );

    let network = solve_dirichlet(&bp, &[("inlet", P_IN), ("outlet", P_OUT)]);

    let r = edge_resistance(&network, "pipe1");
    assert!(r > 0.0, "resistance must be positive; got {r}");

    let q = edge_flow(&network, "pipe1");
    let q_expected = (P_IN - P_OUT) / r;

    let rel_err = ((q - q_expected) / q_expected).abs();
    assert!(
        rel_err < 1e-6,
        "series Q = {q:.6e} m³/s, Kirchhoff expects {q_expected:.6e}, rel_err = {rel_err:.2e}"
    );

    // Pressures at boundary nodes must be exactly as prescribed.
    let p_in_solved = node_pressure(&network, "inlet");
    let p_out_solved = node_pressure(&network, "outlet");
    assert!(
        (p_in_solved - P_IN).abs() < 1e-8,
        "inlet pressure not enforced: {p_in_solved}"
    );
    assert!(
        (p_out_solved - P_OUT).abs() < 1e-8,
        "outlet pressure not enforced: {p_out_solved}"
    );
}

// ── Test 2: Symmetric parallel bifurcation ───────────────────────────────────

/// Inlet → trunk → junction → {branch_L, branch_R} → {outlet_L, outlet_R}
/// All channels identical circular pipes.
///
/// Analytical solution (identical R for all channels, P_in=1000 Pa, P_outs=0):
///   P_junction = P_in / 3
///   Q_branch   = P_junction / R_branch = P_in / (3·R)
///   Q_trunk    = 2 · Q_branch           = 2·P_in / (3·R)
///
/// Mass conservation: Q_trunk = Q_branch_L + Q_branch_R.
#[test]
fn symmetric_bifurcation_flow_distribution_and_mass_conservation() {
    const D_M: f64 = 1e-3;
    const L_M: f64 = 0.05; // 50 mm each segment
    const P_IN: f64 = 1000.0;
    const P_OUT: f64 = 0.0;

    let mut bp = NetworkBlueprint::new_with_explicit_positions("bifurcation");
    bp.add_node(node("inlet", NodeKind::Inlet, 0.0, 0.0));
    bp.add_node(node("junction", NodeKind::Junction, 50.0, 0.0));
    bp.add_node(node("out_L", NodeKind::Outlet, 90.0, -30.0));
    bp.add_node(node("out_R", NodeKind::Outlet, 90.0, 30.0));

    bp.add_channel(
        pipe("trunk", "inlet", "junction", L_M, D_M).with_path(vec![(0.0, 0.0), (50.0, 0.0)]),
    );
    bp.add_channel(
        pipe("branch_L", "junction", "out_L", L_M, D_M).with_path(vec![(50.0, 0.0), (90.0, -30.0)]),
    );
    bp.add_channel(
        pipe("branch_R", "junction", "out_R", L_M, D_M).with_path(vec![(50.0, 0.0), (90.0, 30.0)]),
    );

    let network = solve_dirichlet(&bp, &[("inlet", P_IN), ("out_L", P_OUT), ("out_R", P_OUT)]);

    let r_trunk = edge_resistance(&network, "trunk");
    let r_branch_l = edge_resistance(&network, "branch_L");
    let r_branch_r = edge_resistance(&network, "branch_R");

    // All three channels have identical geometry → resistances must be equal.
    let rel_diff_lr = ((r_branch_l - r_branch_r) / r_branch_l).abs();
    assert!(
        rel_diff_lr < 1e-10,
        "branches should have equal resistance; L={r_branch_l:.6e}, R={r_branch_r:.6e}"
    );

    let q_trunk = edge_flow(&network, "trunk");
    let q_branch_l = edge_flow(&network, "branch_L");
    let q_branch_r = edge_flow(&network, "branch_R");

    // ── Mass conservation at junction ────────────────────────────────────────
    // Q_trunk = Q_branch_L + Q_branch_R (signs: trunk in, branches out).
    let mass_error = (q_trunk - (q_branch_l + q_branch_r)).abs();
    let mass_rel = mass_error / q_trunk.abs().max(1e-30);
    assert!(
        mass_rel < 1e-8,
        "mass not conserved at junction: Q_trunk={q_trunk:.6e}, Q_L={q_branch_l:.6e}, Q_R={q_branch_r:.6e}"
    );

    // ── Symmetric branches must carry equal flow ──────────────────────────────
    let sym_rel = ((q_branch_l - q_branch_r) / q_branch_l).abs();
    assert!(
        sym_rel < 1e-6,
        "symmetric branches carry unequal flow: Q_L={q_branch_l:.6e}, Q_R={q_branch_r:.6e}"
    );

    // ── Kirchhoff KVL along each path from inlet to outlet ───────────────────
    // P_in → trunk → P_junction → branch_L → P_out_L
    // ΔP_trunk + ΔP_branch_L = P_in − P_out_L = P_in
    let p_junction = node_pressure(&network, "junction");
    let dp_trunk = P_IN - p_junction;
    let dp_branch_l = p_junction - P_OUT;

    let kvl_trunk = (dp_trunk - q_trunk * r_trunk).abs() / (P_IN + 1e-30);
    let kvl_branch_l = (dp_branch_l - q_branch_l * r_branch_l).abs() / (P_IN + 1e-30);
    assert!(
        kvl_trunk < 1e-6,
        "KVL violation on trunk: ΔP={dp_trunk:.4e}, Q·R={:.4e}",
        q_trunk * r_trunk
    );
    assert!(
        kvl_branch_l < 1e-6,
        "KVL violation on branch_L: ΔP={dp_branch_l:.4e}, Q·R={:.4e}",
        q_branch_l * r_branch_l
    );

    // ── Analytical junction pressure for identical resistances ───────────────
    // P_junction = P_in / 3  (trunk R = branch R, two branches in parallel)
    let p_jn_expected = P_IN / 3.0;
    let p_jn_rel_err = ((p_junction - p_jn_expected) / p_jn_expected).abs();
    assert!(
        p_jn_rel_err < 0.02, // 2% tolerance — junction losses may shift slightly
        "junction pressure {p_junction:.2} Pa deviates from P_in/3={p_jn_expected:.2} Pa by {:.2}%",
        p_jn_rel_err * 100.0
    );
}

// ── Test 3: Asymmetric bifurcation — Kirchhoff ratio ─────────────────────────

/// Inlet → junction → {narrow_branch (D/2), wide_branch (D)} → outlets.
/// Both outlets at P=0; inlet is a flow source node (Dirichlet at junction).
///
/// With R ∝ 1/D⁴ (Hagen-Poiseuille), R_narrow = 16 · R_wide.
/// So Q_narrow / Q_wide = R_wide / R_narrow = 1/16.
/// Both branches see the same ΔP = P_junction − 0.
#[test]
fn asymmetric_bifurcation_flow_ratio_matches_resistance_ratio() {
    const D_WIDE: f64 = 1e-3; // 1 mm
    const D_NARROW: f64 = 5e-4; // 0.5 mm → R_narrow = 2⁴ = 16 × R_wide
    const L_M: f64 = 0.05;
    const P_JN: f64 = 500.0; // Dirichlet at junction, both outlets at 0

    let mut bp = NetworkBlueprint::new_with_explicit_positions("asym_bifurcation");
    bp.add_node(node("junction", NodeKind::Inlet, 0.0, 0.0));
    bp.add_node(node("out_wide", NodeKind::Outlet, 40.0, -30.0));
    bp.add_node(node("out_narr", NodeKind::Outlet, 40.0, 30.0));

    bp.add_channel(
        pipe("branch_wide", "junction", "out_wide", L_M, D_WIDE)
            .with_path(vec![(0.0, 0.0), (40.0, -30.0)]),
    );
    bp.add_channel(
        pipe("branch_narr", "junction", "out_narr", L_M, D_NARROW)
            .with_path(vec![(0.0, 0.0), (40.0, 30.0)]),
    );

    let network = solve_dirichlet(
        &bp,
        &[("junction", P_JN), ("out_wide", 0.0), ("out_narr", 0.0)],
    );

    let r_wide = edge_resistance(&network, "branch_wide");
    let r_narr = edge_resistance(&network, "branch_narr");
    let q_wide = edge_flow(&network, "branch_wide");
    let q_narr = edge_flow(&network, "branch_narr");

    // Both channels have the same pressure drop = P_JN.
    let q_wide_expected = P_JN / r_wide;
    let q_narr_expected = P_JN / r_narr;

    let err_wide = ((q_wide - q_wide_expected) / q_wide_expected).abs();
    let err_narr = ((q_narr - q_narr_expected) / q_narr_expected).abs();

    assert!(
        err_wide < 1e-6,
        "wide branch: Q={q_wide:.4e}, expected P/R={q_wide_expected:.4e}, err={err_wide:.2e}"
    );
    assert!(
        err_narr < 1e-6,
        "narrow branch: Q={q_narr:.4e}, expected P/R={q_narr_expected:.4e}, err={err_narr:.2e}"
    );

    // Resistance ratio: R_narr/R_wide should be (D_wide/D_narrow)^4 = 2^4 = 16.
    let r_ratio = r_narr / r_wide;
    let expected_ratio = (D_WIDE / D_NARROW).powi(4);
    let ratio_err = ((r_ratio - expected_ratio) / expected_ratio).abs();
    assert!(
        ratio_err < 0.01,
        "R_narr/R_wide = {r_ratio:.4}, expected (D_wide/D_narrow)^4 = {expected_ratio:.1}, err = {:.2}%", ratio_err * 100.0
    );

    // Flow ratio is inverse of resistance ratio.
    let q_ratio = q_wide / q_narr;
    let q_ratio_err = ((q_ratio - expected_ratio) / expected_ratio).abs();
    assert!(
        q_ratio_err < 0.01,
        "Q_wide/Q_narr = {q_ratio:.4}, expected {expected_ratio:.1}, err = {:.2}%",
        q_ratio_err * 100.0
    );
}

// ── Test 4: Primitive selective tree — per-node/per-channel trace ─────────────

/// Builds a Tri-Bi primitive selective tree blueprint from cfd-schematics,
/// solves it with a fixed inlet flow rate, and validates:
///   1. All node pressures are finite and positive (monotone decrease toward outlets).
///   2. All channel flow rates are finite and non-zero.
///   3. Mass conservation holds at every interior node.
///   4. Total outlet flow equals total inlet flow.
#[test]
fn primitive_selective_tree_trace_all_nodes_channels() {
    use cfd_1d::NodeType;
    use cfd_schematics::geometry::generator::PrimitiveSelectiveSplitKind;
    use cfd_schematics::interface::presets::primitive_selective_split_tree_rect;

    let blueprint = primitive_selective_split_tree_rect(
        "trace-test",
        (127.76, 85.47),
        &[
            PrimitiveSelectiveSplitKind::Tri,
            PrimitiveSelectiveSplitKind::Bi,
        ],
        4.0e-3,  // main width
        0.50,    // first trifurcation center frac
        0.55,    // later center frac
        0.60,    // bifurcation treatment frac
        35.0e-6, // throat width
        80.0e-6, // throat length
        1.0e-3,  // channel height
        false,   // no venturi (keeps test simple)
        0,
        None,
    );

    const Q_INLET: f64 = 1e-6; // 1 mL/min
    const P_REF: f64 = 0.0;

    let fluid = water_20c::<f64>().expect("water_20c");
    let mut network = network_from_blueprint(&blueprint, fluid).expect("blueprint valid");

    let inlet_nodes: Vec<_> = network
        .graph
        .node_indices()
        .filter(|&i| {
            network
                .graph
                .node_weight(i)
                .map_or(false, |n| n.node_type == NodeType::Inlet)
        })
        .collect();
    let outlet_nodes: Vec<_> = network
        .graph
        .node_indices()
        .filter(|&i| {
            network
                .graph
                .node_weight(i)
                .map_or(false, |n| n.node_type == NodeType::Outlet)
        })
        .collect();

    assert!(
        !inlet_nodes.is_empty(),
        "blueprint must have at least one inlet"
    );
    assert!(
        !outlet_nodes.is_empty(),
        "blueprint must have at least one outlet"
    );

    let q_per_inlet = Q_INLET / inlet_nodes.len() as f64;
    for &i in &inlet_nodes {
        network.set_neumann_flow(i, q_per_inlet);
    }
    for &o in &outlet_nodes {
        network.set_pressure(o, P_REF);
    }

    let solver = NetworkSolver::with_config(SolverConfig {
        tolerance: 1e-9,
        max_iterations: 500,
    });
    let (network, diag) = solver
        .solve_network_with_diagnostics(&NetworkProblem::new(network))
        .expect("solve must succeed");

    assert!(
        diag.failure_detail.is_none(),
        "solver reported failure: {:?}",
        diag.failure_detail
    );

    // ── 1. All node pressures are finite ────────────────────────────────────
    let pressures = network.pressures();
    for (idx, &p) in pressures.iter().enumerate() {
        let id = network
            .graph
            .node_weight(petgraph::graph::NodeIndex::new(idx))
            .map(|n| n.id.as_str())
            .unwrap_or("?");
        assert!(p.is_finite(), "node '{id}' has non-finite pressure {p}");
        // Outlet nodes must be at P_REF; all others above it (flow driven network).
        let is_outlet = network
            .graph
            .node_weight(petgraph::graph::NodeIndex::new(idx))
            .map_or(false, |n| n.node_type == NodeType::Outlet);
        if is_outlet {
            assert!(
                (p - P_REF).abs() < 1e-6,
                "outlet '{id}' pressure {p} deviates from P_REF={P_REF}"
            );
        } else {
            assert!(p >= P_REF, "non-outlet '{id}' pressure {p} below P_REF");
        }
    }

    // ── 2. All channel flow rates are finite and non-zero ────────────────────
    let flow_rates = network.flow_rates();
    for (idx, &q) in flow_rates.iter().enumerate() {
        let id = network
            .graph
            .edge_weight(petgraph::graph::EdgeIndex::new(idx))
            .map(|e| e.id.as_str())
            .unwrap_or("?");
        assert!(q.is_finite(), "edge '{id}' has non-finite flow {q}");
        assert!(q.abs() > 0.0, "edge '{id}' has zero flow");
    }

    // ── 3. Mass conservation at every interior node ──────────────────────────
    // NetworkGraph is a directed graph: edges(node) returns OUTGOING only.
    // Use edges_directed to get both incoming and outgoing, then compute
    // net_outflow = Σ_outgoing Q - Σ_incoming Q = 0 at junctions.
    use petgraph::visit::EdgeRef;
    use petgraph::Direction;
    for node_idx in network.graph.node_indices() {
        let node = match network.graph.node_weight(node_idx) {
            Some(n) => n,
            None => continue,
        };
        if node.node_type == NodeType::Inlet || node.node_type == NodeType::Outlet {
            continue; // BCs applied at boundaries; skip
        }

        let outflow: f64 = network
            .graph
            .edges_directed(node_idx, Direction::Outgoing)
            .map(|eref| flow_rates.get(eref.id().index()).copied().unwrap_or(0.0))
            .sum();
        let inflow: f64 = network
            .graph
            .edges_directed(node_idx, Direction::Incoming)
            .map(|eref| flow_rates.get(eref.id().index()).copied().unwrap_or(0.0))
            .sum();
        let net_flow = outflow - inflow;

        assert!(
            net_flow.abs() < Q_INLET * 1e-6,
            "mass not conserved at junction '{}': net flow = {net_flow:.3e} m³/s (out={outflow:.3e}, in={inflow:.3e})",
            node.id
        );
    }

    // ── 4. Total outlet flow equals total inlet flow ─────────────────────────
    // For each outlet node, sum incoming edge flows (the network pushes flow INTO outlets).
    let total_outlet_flow: f64 = outlet_nodes
        .iter()
        .map(|&outlet_idx| {
            network
                .graph
                .edges_directed(outlet_idx, Direction::Incoming)
                .map(|eref| flow_rates.get(eref.id().index()).copied().unwrap_or(0.0))
                .sum::<f64>()
        })
        .sum::<f64>()
        .abs();

    let global_mass_err = ((total_outlet_flow - Q_INLET) / Q_INLET).abs();
    assert!(
        global_mass_err < 1e-6,
        "global mass balance: Q_in={Q_INLET:.3e}, Q_out={total_outlet_flow:.3e}, rel_err={global_mass_err:.2e}"
    );
}
