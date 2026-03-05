
use cfd_1d::domain::network::{NetworkBuilder, Node, Edge, NodeType, EdgeType};

#[test]
fn test_new_components_integration() {
    let mut builder = NetworkBuilder::<f64>::new();
    
    // Add nodes including new Reservoir type
    // NetworkBuilder returns NodeIndex, which we need for edges
    let inlet = builder.add_node(Node::new("inlet".to_string(), NodeType::Inlet));
    let res1 = builder.add_node(Node::new("res1".to_string(), NodeType::Reservoir));
    let j1 = builder.add_node(Node::new("j1".to_string(), NodeType::Junction));
    let outlet = builder.add_node(Node::new("outlet".to_string(), NodeType::Outlet));
    
    // Add edges including Valve and Pump
    builder.add_edge(inlet, res1, Edge::new("pump1".to_string(), EdgeType::Pump));
    builder.add_edge(res1, j1, Edge::new("valve1".to_string(), EdgeType::Valve));
    builder.add_edge(j1, outlet, Edge::new("pipe1".to_string(), EdgeType::Pipe));
    
    let network = builder.build().expect("Failed to build network with new components");
    
    assert_eq!(network.node_count(), 4);
    assert_eq!(network.edge_count(), 3);
    
    // Verify types by iterating graph weights
    let pump = network.edge_weights().find(|e| e.id == "pump1").unwrap();
    assert_eq!(pump.edge_type, EdgeType::Pump);
    
    let valve = network.edge_weights().find(|e| e.id == "valve1").unwrap();
    assert_eq!(valve.edge_type, EdgeType::Valve);
    
    let res = network.node_weights().find(|n| n.id == "res1").unwrap();
    assert_eq!(res.node_type, NodeType::Reservoir);
}

#[test]
fn test_from_spec_conversion() {
    use cfd_schematics::domain::model::{NodeSpec, ChannelSpec, NodeId, NodeKind, EdgeKind};
    use cfd_1d::domain::network::{Node, Edge};
    
    let n_spec = NodeSpec {
        id: NodeId::new("n1"),
        kind: NodeKind::Reservoir,
    };
    
    let node: Node<f64> = Node::from(&n_spec);
    assert_eq!(node.id, "n1");
    // assert!(matches!(node.node_type, NodeKind::Reservoir)); 
    // Since NodeKind is PartialEq, we can just compare
    assert_eq!(node.node_type, NodeKind::Reservoir);
    
    let c_spec = ChannelSpec::new_valve("c1", "n1", "n2", 0.5);
    
    let edge: Edge<f64> = Edge::from(&c_spec);
    assert_eq!(edge.id, "c1");
    assert_eq!(edge.edge_type, EdgeKind::Valve);
}

#[test]
fn test_blueprint_negative_length_rejected() {
    use cfd_schematics::domain::model::{NetworkBlueprint, NodeSpec, ChannelSpec, CrossSectionSpec, ChannelShape, NodeId, NodeKind, EdgeId, EdgeKind};
    use cfd_core::physics::fluid::database::water_20c;
    let fluid = water_20c::<f64>().unwrap();
    let n1 = NodeSpec { id: NodeId::new("in"), kind: NodeKind::Inlet };
    let n2 = NodeSpec { id: NodeId::new("out"), kind: NodeKind::Outlet };
    let c1 = ChannelSpec {
        id: EdgeId::new("c1"),
        kind: EdgeKind::Pipe,
        from: NodeId::new("in"),
        to: NodeId::new("out"),
        length_m: -0.1, // Invalid negative length!
        resistance: 100.0,
        cross_section: CrossSectionSpec::Circular { diameter_m: 0.001 },
        channel_shape: ChannelShape::Straight,
        quad_coeff: 0.0,
        valve_cv: None,
        pump_max_flow: None,
        pump_max_pressure: None,
        metadata: None,
    };
    
    let blueprint = NetworkBlueprint {
        name: "test_neg_len".to_string(),
        nodes: vec![n1, n2],
        channels: vec![c1],
    };
    
    let result = cfd_1d::domain::network::network_from_blueprint(&blueprint, fluid);
    assert!(result.is_err(), "Blueprint with negative length must be rejected");
}

#[test]
fn test_blueprint_zero_diameter_rejected() {
    use cfd_schematics::domain::model::{NetworkBlueprint, NodeSpec, ChannelSpec, CrossSectionSpec, ChannelShape, NodeId, NodeKind, EdgeId, EdgeKind};
    use cfd_core::physics::fluid::database::water_20c;
    let fluid = water_20c::<f64>().unwrap();
    let n1 = NodeSpec { id: NodeId::new("in"), kind: NodeKind::Inlet };
    let n2 = NodeSpec { id: NodeId::new("out"), kind: NodeKind::Outlet };
    let c1 = ChannelSpec {
        id: EdgeId::new("c1"),
        kind: EdgeKind::Pipe,
        from: NodeId::new("in"),
        to: NodeId::new("out"),
        length_m: 0.1, 
        resistance: 100.0,
        cross_section: CrossSectionSpec::Circular { diameter_m: 0.0 }, // Invalid zero diameter
        channel_shape: ChannelShape::Straight,
        quad_coeff: 0.0,
        valve_cv: None,
        pump_max_flow: None,
        pump_max_pressure: None,
        metadata: None,
    };
    
    let blueprint = NetworkBlueprint {
        name: "test_zero_diam".to_string(),
        nodes: vec![n1, n2],
        channels: vec![c1],
    };
    
    let result = cfd_1d::domain::network::network_from_blueprint(&blueprint, fluid);
    assert!(result.is_err(), "Blueprint with zero diameter must be rejected");
}
