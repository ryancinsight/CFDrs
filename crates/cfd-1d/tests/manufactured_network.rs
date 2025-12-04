use cfd_core::error::Result;
use cfd_core::solver::Solver;
use cfd_1d::{NetworkBuilder, Edge, EdgeType, Network, EdgeProperties, NetworkProblem, NetworkSolver};
use cfd_1d::network::ComponentType;
use cfd_core::fluid::database;
use approx::assert_relative_eq;

#[test]
fn test_series_additivity() -> Result<()> {
    type T = f64;
    let mut builder = NetworkBuilder::<T>::new();

    let n_in = builder.add_inlet("in".to_string());
    let n_mid = builder.add_junction("mid".to_string());
    let n_out = builder.add_outlet("out".to_string());

    let mut e1 = Edge::<T>::new("e1".to_string(), EdgeType::Pipe);
    e1.resistance = 2.0;
    let mut e2 = Edge::<T>::new("e2".to_string(), EdgeType::Pipe);
    e2.resistance = 3.0;

    let eidx1 = builder.add_edge(n_in, n_mid, e1);
    let eidx2 = builder.add_edge(n_mid, n_out, e2);

    let graph = builder.build()?;
    let fluid = database::water_20c::<T>()?;
    let mut net = Network::new(graph, fluid);

    let props1 = EdgeProperties { id: "e1".into(), component_type: ComponentType::Pipe,
        length: 1.0, area: 1.0, hydraulic_diameter: None, resistance: 2.0, geometry: None, properties: std::collections::HashMap::new() };
    let props2 = EdgeProperties { id: "e2".into(), component_type: ComponentType::Pipe,
        length: 1.0, area: 1.0, hydraulic_diameter: None, resistance: 3.0, geometry: None, properties: std::collections::HashMap::new() };
    net.add_edge_properties(eidx1, props1);
    net.add_edge_properties(eidx2, props2);

    net.set_pressure(n_in, 1.0);
    net.set_pressure(n_out, 0.0);

    let problem = NetworkProblem::new(net.clone());
    let solver = NetworkSolver::<T>::new();
    let solved = solver.solve(&problem)?;

    let q_expected = 1.0 / 5.0;
    let q1 = *solved.flow_rates().get(&eidx1).unwrap();
    let q2 = *solved.flow_rates().get(&eidx2).unwrap();
    assert_relative_eq!(q1, q_expected, max_relative = 1e-12);
    assert_relative_eq!(q2, q_expected, max_relative = 1e-12);

    Ok(())
}

#[test]
fn test_parallel_quadratic_branches() -> Result<()> {
    type T = f64;
    let mut builder = NetworkBuilder::<T>::new();
    let n_in = builder.add_inlet("in".to_string());
    let n_split = builder.add_junction("split".to_string());
    let n_out = builder.add_outlet("out".to_string());

    let mut e1 = Edge::<T>::new("e1".to_string(), EdgeType::Pipe);
    e1.resistance = 1.0;
    e1.quad_coeff = 0.5;
    let mut e2 = Edge::<T>::new("e2".to_string(), EdgeType::Pipe);
    e2.resistance = 2.0;
    e2.quad_coeff = 1.0;
    let mut merge = Edge::<T>::new("merge".to_string(), EdgeType::Pipe);
    merge.resistance = 1e-9;

    let eidx1 = builder.add_edge(n_in, n_split, e1);
    let eidx2 = builder.add_edge(n_in, n_split, e2);
    let eidx3 = builder.add_edge(n_split, n_out, merge);

    let graph = builder.build()?;
    let fluid = database::water_20c::<T>()?;
    let mut net = Network::new(graph, fluid);

    let props1 = EdgeProperties { id: "e1".into(), component_type: ComponentType::Pipe,
        length: 1.0, area: 1.0, hydraulic_diameter: None, resistance: 1.0, geometry: None, properties: std::collections::HashMap::new() };
    let props2 = EdgeProperties { id: "e2".into(), component_type: ComponentType::Pipe,
        length: 1.0, area: 1.0, hydraulic_diameter: None, resistance: 2.0, geometry: None, properties: std::collections::HashMap::new() };
    let props3 = EdgeProperties { id: "merge".into(), component_type: ComponentType::Pipe,
        length: 1.0, area: 1.0, hydraulic_diameter: None, resistance: 1e-9, geometry: None, properties: std::collections::HashMap::new() };
    net.add_edge_properties(eidx1, props1);
    net.add_edge_properties(eidx2, props2);
    net.add_edge_properties(eidx3, props3);

    if let Some(edge) = net.graph.edge_weight_mut(eidx1) { edge.quad_coeff = 0.5; }
    if let Some(edge) = net.graph.edge_weight_mut(eidx2) { edge.quad_coeff = 1.0; }

    net.set_pressure(n_in, 2.0);
    net.set_pressure(n_out, 0.0);

    let problem = NetworkProblem::new(net.clone());
    let solver = NetworkSolver::<T>::new();
    let solved = solver.solve(&problem)?;

    let dp = 2.0_f64;
    let q1_mag = ((1.0_f64*1.0_f64 + 4.0_f64*0.5_f64*dp).sqrt() - 1.0_f64) / (2.0_f64*0.5_f64);
    let q2_mag = ((2.0_f64*2.0_f64 + 4.0_f64*1.0_f64*dp).sqrt() - 2.0_f64) / (2.0_f64*1.0_f64);

    let q1 = *solved.flow_rates().get(&eidx1).unwrap();
    let q2 = *solved.flow_rates().get(&eidx2).unwrap();
    assert_relative_eq!(q1, q1_mag, max_relative = 1e-9);
    assert_relative_eq!(q2, q2_mag, max_relative = 1e-9);

    let q_total_expected = q1_mag + q2_mag;
    let mut q_total = 0.0;
    for (edge_idx, &q) in solved.flow_rates().iter() {
        let (from, _to) = solved.graph.edge_endpoints(*edge_idx).unwrap();
        if from == n_in { q_total += q; }
    }
    assert_relative_eq!(q_total, q_total_expected, max_relative = 1e-9);

    Ok(())
}

#[test]
fn test_parallel_additivity() -> Result<()> {
    type T = f64;
    let mut builder = NetworkBuilder::<T>::new();
    let n_in = builder.add_inlet("in".to_string());
    let n_split = builder.add_junction("split".to_string());
    let n_out = builder.add_outlet("out".to_string());

    let mut e1 = Edge::<T>::new("e1".to_string(), EdgeType::Pipe);
    e1.resistance = 2.0;
    let mut e2 = Edge::<T>::new("e2".to_string(), EdgeType::Pipe);
    e2.resistance = 3.0;
    let mut merge = Edge::<T>::new("merge".to_string(), EdgeType::Pipe);
    merge.resistance = 1e-9;

    let eidx1 = builder.add_edge(n_in, n_split, e1);
    let eidx2 = builder.add_edge(n_in, n_split, e2);
    let eidx3 = builder.add_edge(n_split, n_out, merge);

    let graph = builder.build()?;
    let fluid = database::water_20c::<T>()?;
    let mut net = Network::new(graph, fluid);

    let props1 = EdgeProperties { id: "e1".into(), component_type: ComponentType::Pipe,
        length: 1.0, area: 1.0, hydraulic_diameter: None, resistance: 2.0, geometry: None, properties: std::collections::HashMap::new() };
    let props2 = EdgeProperties { id: "e2".into(), component_type: ComponentType::Pipe,
        length: 1.0, area: 1.0, hydraulic_diameter: None, resistance: 3.0, geometry: None, properties: std::collections::HashMap::new() };
    let props3 = EdgeProperties { id: "merge".into(), component_type: ComponentType::Pipe,
        length: 1.0, area: 1.0, hydraulic_diameter: None, resistance: 1e-9, geometry: None, properties: std::collections::HashMap::new() };
    net.add_edge_properties(eidx1, props1);
    net.add_edge_properties(eidx2, props2);
    net.add_edge_properties(eidx3, props3);

    net.set_pressure(n_in, 1.0);
    net.set_pressure(n_out, 0.0);

    let problem = NetworkProblem::new(net.clone());
    let solver = NetworkSolver::<T>::new();
    let solved = solver.solve(&problem)?;

    let q_total_expected = 1.0 * (1.0/2.0 + 1.0/3.0);
    let mut q_total = 0.0;
    for (edge_idx, &q) in solved.flow_rates().iter() {
        let (from, _to) = solved.graph.edge_endpoints(*edge_idx).unwrap();
        if from == n_in { q_total += q; }
    }
    assert_relative_eq!(q_total, q_total_expected, max_relative = 1e-9);

    Ok(())
}

#[test]
fn test_conservation_at_junction() -> Result<()> {
    type T = f64;
    let mut builder = NetworkBuilder::<T>::new();
    let n_in = builder.add_inlet("in".to_string());
    let n_j = builder.add_junction("j".to_string());
    let n_out1 = builder.add_outlet("o1".to_string());
    let n_out2 = builder.add_outlet("o2".to_string());

    let mut e_in = Edge::<T>::new("ein".to_string(), EdgeType::Pipe);
    e_in.resistance = 1.0;
    let mut e1 = Edge::<T>::new("e1".to_string(), EdgeType::Pipe);
    e1.resistance = 2.0;
    let mut e2 = Edge::<T>::new("e2".to_string(), EdgeType::Pipe);
    e2.resistance = 3.0;

    let _ein = builder.add_edge(n_in, n_j, e_in);
    let _eidx1 = builder.add_edge(n_j, n_out1, e1);
    let _eidx2 = builder.add_edge(n_j, n_out2, e2);

    let graph = builder.build()?;
    let fluid = database::water_20c::<T>()?;
    let mut net = Network::new(graph, fluid);
    net.set_pressure(n_in, 1.0);
    net.set_pressure(n_out1, 0.0);
    net.set_pressure(n_out2, 0.0);

    let problem = NetworkProblem::new(net.clone());
    let solver = NetworkSolver::<T>::new();
    let solved = solver.solve(&problem)?;

    let mut net = 0.0;
    for (edge_idx, &q) in solved.flow_rates().iter() {
        let (from, to) = solved.graph.edge_endpoints(*edge_idx).unwrap();
        if to == n_j { net += q; }
        if from == n_j { net -= q; }
    }
    assert_relative_eq!(net, 0.0, max_relative = 1e-12);
    Ok(())
}

#[test]
fn test_quadratic_resistance() -> Result<()> {
    // Single pipe with quadratic resistance: ΔP = R·Q + k·Q²
    // R = 1.0, k = 1.0. ΔP = 2.0.
    // 2 = 1·Q + 1·Q² => Q² + Q - 2 = 0 => (Q+2)(Q-1) = 0 => Q = 1.0 (since Q > 0)
    
    type T = f64;
    let mut builder = NetworkBuilder::<T>::new();
    let n_in = builder.add_inlet("in".to_string());
    let n_out = builder.add_outlet("out".to_string());

    let mut e1 = Edge::<T>::new("e1".to_string(), EdgeType::Pipe);
    e1.resistance = 1.0;
    e1.quad_coeff = 1.0; // Non-zero quadratic coefficient

    let eidx1 = builder.add_edge(n_in, n_out, e1);

    let graph = builder.build()?;
    let fluid = database::water_20c::<T>()?;
    let mut net = Network::new(graph, fluid);

    // Set properties
    let props1 = EdgeProperties { 
        id: "e1".into(), 
        component_type: ComponentType::Pipe,
        length: 1.0, 
        area: 1.0, 
        hydraulic_diameter: None, 
        resistance: 1.0, 
        geometry: None, 
        properties: std::collections::HashMap::new() 
    };
    net.add_edge_properties(eidx1, props1);
    
    // We must manually set the quad_coeff in the network graph edge weight because 
    // add_edge_properties might not overwrite it or might not expose it.
    // The Edge struct has quad_coeff. Let's make sure it's set.
    if let Some(edge) = net.graph.edge_weight_mut(eidx1) {
        edge.quad_coeff = 1.0;
    }

    net.set_pressure(n_in, 2.0);
    net.set_pressure(n_out, 0.0);

    let problem = NetworkProblem::new(net.clone());
    let solver = NetworkSolver::<T>::new();
    let solved = solver.solve(&problem)?;

    let q = *solved.flow_rates().get(&eidx1).unwrap();
    
    // Expect Q = 1.0
    assert_relative_eq!(q, 1.0, max_relative = 1e-6);

    Ok(())
}
