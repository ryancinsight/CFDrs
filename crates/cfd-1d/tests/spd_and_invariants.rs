use cfd_1d::domain::network::{
    Edge, EdgeProperties, EdgeType, Network, NetworkBuilder, ResistanceUpdatePolicy,
};
use cfd_1d::solver::core::{LinearSolverMethod, NetworkProblem, NetworkSolver};
use cfd_1d::{ChannelGeometry, CrossSection, SurfaceProperties, Wettability};
use cfd_core::compute::solver::Solver;
use cfd_core::error::Result;
use cfd_core::physics::fluid::ConstantPropertyFluid;

fn network_two_node<T: nalgebra::RealField + Copy + num_traits::FromPrimitive>() -> (
    Network<T>,
    petgraph::graph::EdgeIndex,
    petgraph::graph::NodeIndex,
    petgraph::graph::NodeIndex,
) {
    let mut builder = NetworkBuilder::<T>::new();
    let inlet = builder.add_inlet("in".to_string());
    let outlet = builder.add_outlet("out".to_string());
    let edge = builder.add_edge(inlet, outlet, Edge::new("e".to_string(), EdgeType::Pipe));
    let graph = builder.build().expect("test invariant");
    let fluid = ConstantPropertyFluid::<T>::water_20c().expect("test invariant");
    let net = Network::new(graph, fluid);
    (net, edge, inlet, outlet)
}

#[test]
fn spd_heuristic_selects_cg() -> Result<()> {
    type F = f64;
    let (mut net, edge, inlet, _outlet) = network_two_node::<F>();
    if let Some(e) = net.graph.edge_weight_mut(edge) {
        e.resistance = 4.0;
        e.quad_coeff = 0.0;
    }
    net.set_pressure(inlet, 10.0);
    let problem = NetworkProblem::new(net.clone());
    let solver = NetworkSolver::<F>::new();
    let solved = solver.solve(&problem)?;
    match solved.last_solver_method.expect("test invariant") {
        LinearSolverMethod::ConjugateGradient => (),
        _ => panic!("expected CG method selection"),
    }
    Ok(())
}

#[test]
fn negative_coefficients_rejected_in_update() {
    type F = f64;
    let (mut net, edge, inlet, outlet) = network_two_node::<F>();
    if let Some(e) = net.graph.edge_weight_mut(edge) {
        e.resistance = -1.0;
        e.quad_coeff = 0.0;
    }
    let mut x = nalgebra::DVector::<F>::zeros(net.node_count());
    x[inlet.index()] = 5.0;
    x[outlet.index()] = 0.0;
    let err = net.update_from_solution(&x).unwrap_err();
    match err {
        cfd_core::error::Error::InvalidConfiguration(_) => (),
        _ => panic!("unexpected error type"),
    }
}

#[test]
fn flow_invariant_edges_skip_resistance_recomputation() -> Result<()> {
    type F = f64;
    let (mut net, edge, _inlet, _outlet) = network_two_node::<F>();
    if let Some(e) = net.graph.edge_weight_mut(edge) {
        e.resistance = 7.5;
        e.quad_coeff = 0.0;
    }
    net.set_flow_rate(edge, 3.0e-9);
    net.add_edge_properties(
        edge,
        EdgeProperties {
            id: "e".to_string(),
            component_type: cfd_1d::ComponentType::Pipe,
            length: 1.0e-2,
            area: 1.0e-6,
            hydraulic_diameter: Some(1.0e-3),
            resistance: 7.5,
            geometry: Some(ChannelGeometry {
                channel_type: cfd_1d::ChannelType::Straight,
                length: 1.0e-2,
                cross_section: CrossSection::Circular { diameter: 1.0e-3 },
                surface: SurfaceProperties {
                    roughness: 0.0,
                    contact_angle: None,
                    surface_energy: None,
                    wettability: Wettability::Hydrophilic,
                },
                variations: Vec::new(),
            }),
            resistance_update_policy: ResistanceUpdatePolicy::FlowInvariant,
            properties: std::collections::HashMap::new(),
        },
    );

    net.update_resistances()?;
    let edge_ref = net.graph.edge_weight(edge).expect("edge must exist");
    assert!((edge_ref.resistance - 7.5).abs() < 1.0e-12);
    assert!(edge_ref.quad_coeff.abs() < 1.0e-12);
    Ok(())
}
