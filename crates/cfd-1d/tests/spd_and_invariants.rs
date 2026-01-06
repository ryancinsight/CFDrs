use cfd_1d::network::{Network, NetworkBuilder, Edge, EdgeType};
use cfd_1d::solver::{NetworkProblem, NetworkSolver, LinearSolverMethod};
use cfd_core::compute::solver::Solver;
use cfd_core::error::Result;
use cfd_core::physics::fluid::ConstantPropertyFluid;

fn network_two_node<T: nalgebra::RealField + Copy + num_traits::FromPrimitive>() -> (Network<T>, petgraph::graph::EdgeIndex, petgraph::graph::NodeIndex, petgraph::graph::NodeIndex) {
    let mut builder = NetworkBuilder::<T>::new();
    let inlet = builder.add_inlet("in".to_string());
    let outlet = builder.add_outlet("out".to_string());
    let edge = builder.add_edge(inlet, outlet, Edge::new("e".to_string(), EdgeType::Pipe));
    let graph = builder.build().unwrap();
    let fluid = ConstantPropertyFluid::<T>::water_20c().unwrap();
    let net = Network::new(graph, fluid);
    (net, edge, inlet, outlet)
}

#[test]
fn spd_heuristic_selects_cg() -> Result<()> {
    type F = f64;
    let (mut net, edge, inlet, _outlet) = network_two_node::<F>();
    if let Some(e) = net.graph.edge_weight_mut(edge) { e.resistance = 4.0; e.quad_coeff = 0.0; }
    net.set_pressure(inlet, 10.0);
    let problem = NetworkProblem::new(net.clone());
    let solver = NetworkSolver::<F>::new();
    let solved = solver.solve(&problem)?;
    match solved.last_solver_method.unwrap() {
        LinearSolverMethod::ConjugateGradient => (),
        _ => panic!("expected CG method selection"),
    }
    Ok(())
}

#[test]
fn negative_coefficients_rejected_in_update() {
    type F = f64;
    let (mut net, edge, inlet, outlet) = network_two_node::<F>();
    if let Some(e) = net.graph.edge_weight_mut(edge) { e.resistance = -1.0; e.quad_coeff = 0.0; }
    let mut x = nalgebra::DVector::<F>::zeros(net.node_count());
    x[inlet.index()] = 5.0; x[outlet.index()] = 0.0;
    let err = net.update_from_solution(&x).unwrap_err();
    match err { cfd_core::error::Error::InvalidConfiguration(_) => (), _ => panic!("unexpected error type") }
}
