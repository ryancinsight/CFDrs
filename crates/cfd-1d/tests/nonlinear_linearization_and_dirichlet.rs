use cfd_1d::network::{Network, NetworkBuilder, Edge, EdgeType, Node, NodeType};
use cfd_1d::solver::{NetworkProblem, NetworkSolver};
use cfd_core::error::Result;
use cfd_core::fluid::newtonian::ConstantPropertyFluid;
use nalgebra::DVector;

fn build_simple_network<T: nalgebra::RealField + Copy + num_traits::FromPrimitive>() -> (Network<T>, petgraph::graph::EdgeIndex, petgraph::graph::NodeIndex, petgraph::graph::NodeIndex) {
    // Build a two-node, one-edge network
    let mut builder = NetworkBuilder::<T>::new();
    let inlet = builder.add_inlet("inlet".to_string());
    let outlet = builder.add_outlet("outlet".to_string());
    let edge_idx = builder.add_edge(inlet, outlet, Edge::new("pipe".to_string(), EdgeType::Pipe));
    let graph = builder.build().unwrap();

    let fluid = ConstantPropertyFluid::<T>::water_20c().unwrap();
    let mut network = Network::new(graph, fluid);
    (network, edge_idx, inlet, outlet)
}

#[test]
fn linearization_effective_resistance_matches_r_plus_2k_abs_q() -> Result<()> {
    type F = f64;
    let (mut network, edge_idx, _, _) = build_simple_network::<F>();

    // Set physical coefficients
    if let Some(edge) = network.graph.edge_weight_mut(edge_idx) {
        edge.resistance = 3.0;
        edge.quad_coeff = 1.5;
        edge.flow_rate = -2.0; // current iterate Q_k
    }

    // Compute parallel edges and check conductance
    let pe: Vec<_> = network.edges_parallel().collect();
    assert_eq!(pe.len(), 1);
    let conductance = pe[0].conductance;

    let q_abs = 2.0_f64;
    let r_eff = 3.0_f64 + 2.0_f64 * 1.5_f64 * q_abs; // R + 2k|Q_k|
    let g_expected = 1.0_f64 / r_eff;
    let rel_err = ((conductance - g_expected) / g_expected).abs();
    assert!(rel_err < 1e-12, "conductance {} vs expected {}", conductance, g_expected);
    Ok(())
}

#[test]
fn update_from_solution_solves_quadratic_with_correct_sign() -> Result<()> {
    type F = f64;
    let (mut network, edge_idx, inlet, outlet) = build_simple_network::<F>();
    // Coefficients
    if let Some(edge) = network.graph.edge_weight_mut(edge_idx) {
        edge.resistance = 2.0;
        edge.quad_coeff = 4.0;
    }

    // Prepare a solution vector with a pressure drop
    let mut x = DVector::<F>::zeros(network.node_count());
    // index mapping: inlet at 0, outlet at 1 in builder insertion order
    x[inlet.index()] = 10.0; // Pa
    x[outlet.index()] = 4.0; // Pa
    network.update_from_solution(&x)?;
    let q = *network.flow_rates.get(&edge_idx).unwrap();
    assert!(q > 0.0, "Flow should align with positive ΔP direction");

    // Closed-form magnitude for Q|Q| model: |Q| = (sqrt(R^2 + 4k|ΔP|) − R)/(2k)
    let dp = 6.0_f64;
    let q_mag = ( (2.0_f64*2.0_f64 + 4.0_f64*4.0_f64*dp).sqrt() - 2.0_f64 ) / (2.0_f64*4.0_f64);
    let rel_err = ((q.abs() - q_mag) / q_mag).abs();
    assert!(rel_err < 1e-12, "Q {} vs expected |Q| {}", q, q_mag);
    Ok(())
}

#[test]
fn dirichlet_enforcement_row_identity_and_rhs() -> Result<()> {
    type F = f64;
    let (mut network, edge_idx, inlet, outlet) = build_simple_network::<F>();
    // Set coefficients to linear to simplify assembly inspection
    if let Some(edge) = network.graph.edge_weight_mut(edge_idx) {
        edge.resistance = 5.0;
        edge.quad_coeff = 0.0;
    }

    // Apply Dirichlet at inlet, free at outlet
    network.set_pressure(inlet, 12.0);

    let assembler = cfd_1d::solver::MatrixAssembler::<F>::new();
    let (a, b) = assembler.assemble(&network)?;

    // Inlet row should be identity
    let row_in = a.row(inlet.index());
    let mut diag_in = 0.0;
    let mut sum_off_in = 0.0;
    for (j, val) in row_in.col_indices().iter().zip(row_in.values()) {
        if *j == inlet.index() { diag_in = *val; } else { sum_off_in += val.abs(); }
    }
    assert_eq!(diag_in, 1.0);
    assert_eq!(sum_off_in, 0.0);
    assert_eq!(b[inlet.index()], 12.0);

    // Outlet diagonal should include conductance contribution, RHS should include Dirichlet value
    let g = 1.0 / 5.0;
    let row_out = a.row(outlet.index());
    let mut diag_out = 0.0;
    let mut off_diag_pos = 0.0;
    for (j, val) in row_out.col_indices().iter().zip(row_out.values()) {
        if *j == outlet.index() { diag_out = *val; } else { off_diag_pos += val; }
    }
    assert!((diag_out - g).abs() < 1e-12);
    assert_eq!(off_diag_pos, 0.0);
    assert!((b[outlet.index()] - g*12.0).abs() < 1e-12);
    Ok(())
}

#[test]
fn dirichlet_enforcement_interior_junction() -> Result<()> {
    type F = f64;
    let mut builder = NetworkBuilder::<F>::new();
    let inlet = builder.add_inlet("inlet".to_string());
    let split = builder.add_junction("split".to_string());
    let outlet = builder.add_outlet("outlet".to_string());

    let mut e1 = Edge::new("b1".to_string(), EdgeType::Pipe);
    let mut e2 = Edge::new("b2".to_string(), EdgeType::Pipe);
    let mut merge = Edge::new("merge".to_string(), EdgeType::Pipe);
    e1.resistance = 4.0;
    e2.resistance = 6.0;
    merge.resistance = 1e-9;

    let _ = builder.add_edge(inlet, split, e1);
    let _ = builder.add_edge(inlet, split, e2);
    let _ = builder.add_edge(split, outlet, merge);
    let graph = builder.build().unwrap();

    let fluid = ConstantPropertyFluid::<F>::water_20c().unwrap();
    let mut network = Network::new(graph, fluid);

    network.set_pressure(split, 10.0);

    let assembler = cfd_1d::solver::MatrixAssembler::<F>::new();
    let (a, b) = assembler.assemble(&network)?;

    let row_split = a.row(split.index());
    let mut diag_split = 0.0;
    let mut sum_off_split = 0.0;
    for (j, val) in row_split.col_indices().iter().zip(row_split.values()) {
        if *j == split.index() { diag_split = *val; } else { sum_off_split += val.abs(); }
    }
    assert_eq!(diag_split, 1.0);
    assert_eq!(sum_off_split, 0.0);
    assert_eq!(b[split.index()], 10.0);

    let g1 = 1.0 / 4.0;
    let g2 = 1.0 / 6.0;
    let g_merge = 1.0 / 1e-9;

    let row_in = a.row(inlet.index());
    let mut rhs_in = 0.0;
    let mut diag_in = 0.0;
    for (j, val) in row_in.col_indices().iter().zip(row_in.values()) {
        if *j == inlet.index() { diag_in = *val; }
    }
    rhs_in = b[inlet.index()];
    assert!((diag_in - (g1 + g2)).abs() < 1e-12);
    assert!((rhs_in - (g1 + g2) * 10.0).abs() < 1e-12);

    let row_out = a.row(outlet.index());
    let mut diag_out = 0.0;
    for (j, val) in row_out.col_indices().iter().zip(row_out.values()) {
        if *j == outlet.index() { diag_out = *val; }
    }
    assert!((diag_out - g_merge).abs() < 1e-12);
    assert!((b[outlet.index()] - g_merge * 10.0).abs() < 1e-12);
    Ok(())
}
