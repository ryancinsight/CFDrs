//! Network topology and graph structures for 1D CFD simulations.
//!
//! This module provides the core network representation for microfluidic and millifluidic
//! systems using graph-based data structures optimized for CFD analysis.

use cfd_core::{Error, Result, Fluid, BoundaryCondition as CoreBoundaryCondition};
use nalgebra::{RealField, Vector3};
use num_traits::cast::FromPrimitive;
use petgraph::{Graph, Directed};
use petgraph::graph::{EdgeIndex, NodeIndex};
use petgraph::visit::EdgeRef;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Type alias for boundary conditions in 1D networks
pub type BoundaryCondition<T> = CoreBoundaryCondition<T>;

/// Network graph type using petgraph
pub type NetworkGraph<N, E> = Graph<N, E, Directed>;

/// Channel properties with clear, self-documenting parameter names
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChannelProperties<T: RealField> {
    /// Hydraulic resistance of the channel
    pub resistance: T,
    /// Physical length of the channel
    pub length: T,
    /// Cross-sectional area of the channel
    pub area: T,
}

impl<T: RealField> ChannelProperties<T> {
    /// Create new channel properties with explicit parameters
    pub fn new(resistance: T, length: T, area: T) -> Self {
        Self {
            resistance,
            length,
            area,
        }
    }
}

/// A node in the microfluidic network
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Node<T: RealField> {
    /// Node identifier
    pub id: String,
    /// Node type
    pub node_type: NodeType,
    /// Position in 2D space (for visualization)
    pub position: Vector3<T>,
    /// Pressure at this node (solution variable)
    pub pressure: Option<T>,
    /// Flow rate through this node
    pub flow_rate: Option<T>,
    /// Node properties
    pub properties: NodeProperties<T>,
}

/// Types of nodes in the network
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum NodeType {
    /// Junction where multiple channels meet
    Junction,
    /// Inlet boundary condition
    Inlet,
    /// Outlet boundary condition
    Outlet,
    /// Component connection point
    Component(String),
    /// Internal node (automatically generated)
    Internal,
}

/// Properties specific to different node types
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NodeProperties<T: RealField> {
    /// Volume of the node (for transient analysis)
    pub volume: Option<T>,
    /// Boundary condition (for inlet/outlet nodes)
    pub boundary_condition: Option<BoundaryCondition<T>>,
    /// Component reference (for component nodes)
    pub component_id: Option<String>,
}



/// An edge representing a channel or component in the network
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Edge<T: RealField> {
    /// Edge identifier
    pub id: String,
    /// Edge type
    pub edge_type: EdgeType,
    /// Hydraulic resistance
    pub resistance: T,
    /// Flow rate through this edge (solution variable)
    pub flow_rate: Option<T>,
    /// Pressure drop across this edge
    pub pressure_drop: Option<T>,
    /// Edge properties
    pub properties: EdgeProperties<T>,
}

/// Types of edges in the network
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum EdgeType {
    /// Physical channel
    Channel,
    /// Pump component
    Pump,
    /// Valve component
    Valve,
    /// Sensor (no resistance)
    Sensor,
    /// Virtual connection
    Virtual,
}

impl EdgeType {
    /// Get static string representation for performance optimization
    pub fn as_str(&self) -> &'static str {
        match self {
            EdgeType::Channel => "Channel",
            EdgeType::Pump => "Pump",
            EdgeType::Valve => "Valve",
            EdgeType::Sensor => "Sensor",
            EdgeType::Virtual => "Virtual",
        }
    }
}

/// Properties specific to different edge types
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EdgeProperties<T: RealField> {
    /// Length of the channel/component
    pub length: Option<T>,
    /// Cross-sectional area
    pub area: Option<T>,
    /// Hydraulic diameter
    pub hydraulic_diameter: Option<T>,
    /// Surface roughness
    pub roughness: Option<T>,
    /// Component-specific parameters
    pub component_params: HashMap<String, T>,
}

/// Main network structure for 1D CFD simulations
#[derive(Debug, Clone)]
pub struct Network<T: RealField> {
    /// The underlying graph structure
    pub graph: NetworkGraph<Node<T>, Edge<T>>,
    /// Fluid properties
    pub fluid: Fluid<T>,
    /// Node index mapping
    node_map: HashMap<String, NodeIndex>,
    /// Edge index mapping
    edge_map: HashMap<String, EdgeIndex>,
    /// Network metadata
    pub metadata: NetworkMetadata<T>,
}

/// Network metadata and global properties
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NetworkMetadata<T: RealField> {
    /// Network name
    pub name: String,
    /// Description
    pub description: String,
    /// Operating temperature
    pub temperature: T,
    /// Operating pressure (reference)
    pub reference_pressure: T,
    /// Gravitational acceleration
    pub gravity: Vector3<T>,
    /// Time step for transient analysis
    pub time_step: Option<T>,
}

impl<T: RealField + FromPrimitive> Default for NodeProperties<T> {
    fn default() -> Self {
        Self {
            volume: None,
            boundary_condition: None,
            component_id: None,
        }
    }
}

impl<T: RealField + FromPrimitive> Default for EdgeProperties<T> {
    fn default() -> Self {
        Self {
            length: None,
            area: None,
            hydraulic_diameter: None,
            roughness: None,
            component_params: HashMap::new(),
        }
    }
}

impl<T: RealField + FromPrimitive> Default for NetworkMetadata<T> {
    fn default() -> Self {
        Self {
            name: "Untitled Network".to_string(),
            description: String::new(),
            temperature: T::from_f64(293.15).unwrap(), // 20°C
            reference_pressure: T::from_f64(101325.0).unwrap(), // 1 atm
            gravity: Vector3::new(T::zero(), T::zero(), -T::from_f64(9.81).unwrap()),
            time_step: None,
        }
    }
}

impl<T: RealField + FromPrimitive> Node<T> {
    /// Create a new node
    pub fn new(id: String, node_type: NodeType, position: Vector3<T>) -> Self {
        Self {
            id,
            node_type,
            position,
            pressure: None,
            flow_rate: None,
            properties: NodeProperties::default(),
        }
    }

    /// Create a junction node
    pub fn junction(id: String, position: Vector3<T>) -> Self {
        Self::new(id, NodeType::Junction, position)
    }

    /// Create an inlet node with pressure boundary condition
    pub fn inlet_pressure(id: String, position: Vector3<T>, pressure: T) -> Self {
        let mut node = Self::new(id, NodeType::Inlet, position);
        node.properties.boundary_condition = Some(BoundaryCondition::pressure_inlet(pressure));
        node
    }

    /// Create an outlet node with pressure boundary condition
    pub fn outlet_pressure(id: String, position: Vector3<T>, pressure: T) -> Self {
        let mut node = Self::new(id, NodeType::Outlet, position);
        node.properties.boundary_condition = Some(BoundaryCondition::pressure_outlet(pressure));
        node
    }

    /// Create an inlet node with flow rate boundary condition
    pub fn inlet_flow_rate(id: String, position: Vector3<T>, flow_rate: T) -> Self {
        let mut node = Self::new(id, NodeType::Inlet, position);
        node.properties.boundary_condition = Some(BoundaryCondition::flow_rate_inlet(flow_rate));
        node
    }

    /// Check if this node has a boundary condition
    pub fn has_boundary_condition(&self) -> bool {
        self.properties.boundary_condition.is_some()
    }

    /// Get the boundary condition if present
    pub fn boundary_condition(&self) -> Option<&BoundaryCondition<T>> {
        self.properties.boundary_condition.as_ref()
    }
}

impl<T: RealField + FromPrimitive> Edge<T> {
    /// Create a new edge
    pub fn new(id: String, edge_type: EdgeType, resistance: T) -> Self {
        Self {
            id,
            edge_type,
            resistance,
            flow_rate: None,
            pressure_drop: None,
            properties: EdgeProperties::default(),
        }
    }

    /// Create a channel edge
    pub fn channel(id: String, resistance: T, length: T, area: T) -> Self {
        let mut edge = Self::new(id, EdgeType::Channel, resistance);
        edge.properties.length = Some(length);
        edge.properties.area = Some(area);
        edge
    }

    /// Create a pump edge
    pub fn pump(id: String, pressure_rise: T) -> Self {
        let mut edge = Self::new(id, EdgeType::Pump, T::zero());
        edge.properties.component_params.insert("pressure_rise".to_string(), pressure_rise);
        edge
    }

    /// Create a valve edge
    pub fn valve(id: String, resistance: T, opening_fraction: T) -> Self {
        let mut edge = Self::new(id, EdgeType::Valve, resistance);
        edge.properties.component_params.insert("opening".to_string(), opening_fraction);
        edge
    }

    /// Get the effective resistance considering component state
    pub fn effective_resistance(&self) -> T {
        match self.edge_type {
            EdgeType::Valve => {
                if let Some(opening) = self.properties.component_params.get("opening") {
                    if *opening <= T::zero() {
                        // Closed valve - infinite resistance
                        T::from_f64(1e12).unwrap()
                    } else {
                        // Resistance inversely proportional to opening
                        self.resistance.clone() / opening.clone()
                    }
                } else {
                    self.resistance.clone()
                }
            },
            EdgeType::Sensor => T::zero(), // No resistance for sensors
            _ => self.resistance.clone(),
        }
    }

    /// Check if this edge represents an active component (pump)
    pub fn is_active_component(&self) -> bool {
        matches!(self.edge_type, EdgeType::Pump)
    }

    /// Get pressure rise for pumps
    pub fn pressure_rise(&self) -> Option<T> {
        if self.edge_type == EdgeType::Pump {
            self.properties.component_params.get("pressure_rise").cloned()
        } else {
            None
        }
    }
}

impl<T: RealField + FromPrimitive> Network<T> {
    /// Create a new empty network
    pub fn new(fluid: Fluid<T>) -> Self {
        Self {
            graph: NetworkGraph::new(),
            fluid,
            node_map: HashMap::new(),
            edge_map: HashMap::new(),
            metadata: NetworkMetadata::default(),
        }
    }

    /// Add a node to the network
    pub fn add_node(&mut self, node: Node<T>) -> Result<NodeIndex> {
        if self.node_map.contains_key(&node.id) {
            return Err(Error::InvalidConfiguration(
                format!("Node with id '{}' already exists", node.id)
            ));
        }

        let node_id = node.id.clone();
        let index = self.graph.add_node(node);
        self.node_map.insert(node_id, index);
        Ok(index)
    }

    /// Add an edge to the network
    pub fn add_edge(&mut self, edge: Edge<T>, from_node: &str, to_node: &str) -> Result<EdgeIndex> {
        let from_index = self.node_map.get(from_node)
            .ok_or_else(|| Error::InvalidConfiguration(
                format!("Source node '{}' not found", from_node)
            ))?;

        let to_index = self.node_map.get(to_node)
            .ok_or_else(|| Error::InvalidConfiguration(
                format!("Target node '{}' not found", to_node)
            ))?;

        if self.edge_map.contains_key(&edge.id) {
            return Err(Error::InvalidConfiguration(
                format!("Edge with id '{}' already exists", edge.id)
            ));
        }

        let edge_id = edge.id.clone();
        let index = self.graph.add_edge(*from_index, *to_index, edge);
        self.edge_map.insert(edge_id, index);
        Ok(index)
    }

    /// Get a node by ID
    pub fn get_node(&self, id: &str) -> Option<&Node<T>> {
        self.node_map.get(id)
            .and_then(|&index| self.graph.node_weight(index))
    }

    /// Get a mutable node by ID
    pub fn get_node_mut(&mut self, id: &str) -> Option<&mut Node<T>> {
        self.node_map.get(id).copied()
            .and_then(move |index| self.graph.node_weight_mut(index))
    }

    /// Get an edge by ID
    pub fn get_edge(&self, id: &str) -> Option<&Edge<T>> {
        self.edge_map.get(id)
            .and_then(|&index| self.graph.edge_weight(index))
    }

    /// Get a mutable edge by ID
    pub fn get_edge_mut(&mut self, id: &str) -> Option<&mut Edge<T>> {
        self.edge_map.get(id).copied()
            .and_then(move |index| self.graph.edge_weight_mut(index))
    }

    /// Get all nodes
    pub fn nodes(&self) -> impl Iterator<Item = &Node<T>> {
        self.graph.node_weights()
    }

    /// Get all edges
    pub fn edges(&self) -> impl Iterator<Item = &Edge<T>> {
        self.graph.edge_weights()
    }

    /// Get the number of nodes
    pub fn node_count(&self) -> usize {
        self.graph.node_count()
    }

    /// Get the number of edges
    pub fn edge_count(&self) -> usize {
        self.graph.edge_count()
    }
    
    /// Get node index by ID
    pub fn get_node_index(&self, id: &str) -> Option<NodeIndex> {
        self.node_map.get(id).copied()
    }
    
    /// Get edge index by ID
    pub fn get_edge_index(&self, id: &str) -> Option<EdgeIndex> {
        self.edge_map.get(id).copied()
    }
    
    /// Get edge ID by index
    pub fn get_edge_id_by_index(&self, edge_idx: EdgeIndex) -> Option<String> {
        self.edge_map.iter()
            .find(|(_, &idx)| idx == edge_idx)
            .map(|(id, _)| id.clone())
    }

    /// Check if the network is connected
    pub fn is_connected(&self) -> bool {
        if self.graph.node_count() <= 1 {
            return true;
        }

        // Use DFS to check connectivity
        use petgraph::visit::Dfs;
        let start_node = self.graph.node_indices().next().unwrap();
        let mut dfs = Dfs::new(&self.graph, start_node);
        let mut visited_count = 0;

        while dfs.next(&self.graph).is_some() {
            visited_count += 1;
        }

        visited_count == self.graph.node_count()
    }

    /// Validate the network structure
    pub fn validate(&self) -> Result<()> {
        // Check for isolated nodes
        for node_index in self.graph.node_indices() {
            let degree = self.graph.edges(node_index).count() +
                        self.graph.edges_directed(node_index, petgraph::Direction::Incoming).count();

            if degree == 0 {
                if let Some(node) = self.graph.node_weight(node_index) {
                    return Err(Error::InvalidConfiguration(
                        format!("Isolated node found: '{}'", node.id)
                    ));
                }
            }
        }

        // Check connectivity
        if !self.is_connected() {
            return Err(Error::InvalidConfiguration(
                "Network is not connected".to_string()
            ));
        }

        // Check for at least one boundary condition
        let has_boundary = self.nodes().any(|node| node.has_boundary_condition());
        if !has_boundary {
            return Err(Error::InvalidConfiguration(
                "Network must have at least one boundary condition".to_string()
            ));
        }

        // Check for valid resistance values
        for edge in self.edges() {
            let resistance = edge.effective_resistance();
            if resistance < T::zero() {
                return Err(Error::InvalidConfiguration(
                    format!("Negative resistance found in edge '{}'", edge.id)
                ));
            }
        }

        Ok(())
    }

    /// Get neighbors of a node
    pub fn neighbors(&self, node_id: &str) -> Result<Vec<&Node<T>>> {
        let node_index = self.node_map.get(node_id)
            .ok_or_else(|| Error::InvalidConfiguration(
                format!("Node '{}' not found", node_id)
            ))?;

        let neighbors: Vec<&Node<T>> = self.graph
            .neighbors(*node_index)
            .filter_map(|idx| self.graph.node_weight(idx))
            .collect();

        Ok(neighbors)
    }

    /// Get edges connected to a node
    pub fn node_edges(&self, node_id: &str) -> Result<Vec<&Edge<T>>> {
        let node_index = self.node_map.get(node_id)
            .ok_or_else(|| Error::InvalidConfiguration(
                format!("Node '{}' not found", node_id)
            ))?;

        let edges: Vec<&Edge<T>> = self.graph
            .edges(*node_index)
            .filter_map(|edge_ref| self.graph.edge_weight(edge_ref.id()))
            .collect();

        Ok(edges)
    }

    /// Clear all solution data (pressures and flow rates)
    pub fn clear_solution(&mut self) {
        for node in self.graph.node_weights_mut() {
            node.pressure = None;
            node.flow_rate = None;
        }

        for edge in self.graph.edge_weights_mut() {
            edge.flow_rate = None;
            edge.pressure_drop = None;
        }
    }

    /// Set fluid properties
    pub fn set_fluid(&mut self, fluid: Fluid<T>) {
        self.fluid = fluid;
    }

    /// Get fluid properties
    pub fn fluid(&self) -> &Fluid<T> {
        &self.fluid
    }
}

/// Builder pattern for constructing networks with error-free configuration
pub struct NetworkBuilder<T: RealField> {
    network: Network<T>,
    /// Pending nodes to be added to the network
    pending_nodes: Vec<Node<T>>,
    /// Pending edges to be added to the network
    pending_edges: Vec<(Edge<T>, String, String)>,
}

impl<T: RealField + FromPrimitive + num_traits::Float> NetworkBuilder<T> {
    /// Create a new network builder
    pub fn new() -> Self {
        Self {
            network: Network::new(Fluid::water()),
            pending_nodes: Vec::new(),
            pending_edges: Vec::new(),
        }
    }

    /// Set the fluid properties
    pub fn with_fluid(mut self, fluid: Fluid<T>) -> Self {
        self.network.set_fluid(fluid);
        self
    }

    /// Set network metadata
    pub fn with_metadata(mut self, metadata: NetworkMetadata<T>) -> Self {
        self.network.metadata = metadata;
        self
    }

    /// Add a junction node - error-free configuration
    pub fn add_junction(mut self, id: &str, x: T, y: T) -> Self {
        let position = Vector3::new(x, y, T::zero());
        let node = Node::junction(id.to_string(), position);
        self.pending_nodes.push(node);
        self
    }

    /// Add an inlet with pressure boundary condition - error-free configuration
    pub fn add_inlet_pressure(mut self, id: &str, x: T, y: T, pressure: T) -> Self {
        let position = Vector3::new(x, y, T::zero());
        let node = Node::inlet_pressure(id.to_string(), position, pressure);
        self.pending_nodes.push(node);
        self
    }

    /// Add an outlet with pressure boundary condition - error-free configuration
    pub fn add_outlet_pressure(mut self, id: &str, x: T, y: T, pressure: T) -> Self {
        let position = Vector3::new(x, y, T::zero());
        let node = Node::outlet_pressure(id.to_string(), position, pressure);
        self.pending_nodes.push(node);
        self
    }

    /// Add an inlet with flow rate boundary condition - error-free configuration
    pub fn add_inlet_flow_rate(mut self, id: &str, x: T, y: T, flow_rate: T) -> Self {
        let position = Vector3::new(x, y, T::zero());
        let node = Node::inlet_flow_rate(id.to_string(), position, flow_rate);
        self.pending_nodes.push(node);
        self
    }

    /// Add a channel between two nodes with explicit properties - error-free configuration
    pub fn add_channel(mut self, id: &str, from: &str, to: &str, props: ChannelProperties<T>) -> Self {
        let edge = Edge::channel(id.to_string(), props.resistance, props.length, props.area);
        self.pending_edges.push((edge, from.to_string(), to.to_string()));
        self
    }

    /// Add a channel between two nodes (legacy method with individual parameters)
    /// 
    /// This method is deprecated in favor of `add_channel` with `ChannelProperties`.
    /// Use `ChannelProperties::new(resistance, length, area)` for clearer parameter naming.
    #[deprecated(since = "0.1.0", note = "Use add_channel with ChannelProperties for clearer parameter naming")]
    pub fn add_channel_legacy(mut self, id: &str, from: &str, to: &str, resistance: T, length: T, area: T) -> Self {
        let edge = Edge::channel(id.to_string(), resistance, length, area);
        self.pending_edges.push((edge, from.to_string(), to.to_string()));
        self
    }

    /// Add a pump between two nodes - error-free configuration
    pub fn add_pump(mut self, id: &str, from: &str, to: &str, pressure_rise: T) -> Self {
        let edge = Edge::pump(id.to_string(), pressure_rise);
        self.pending_edges.push((edge, from.to_string(), to.to_string()));
        self
    }

    /// Add a valve between two nodes - error-free configuration
    pub fn add_valve(mut self, id: &str, from: &str, to: &str, resistance: T, opening: T) -> Self {
        let edge = Edge::valve(id.to_string(), resistance, opening);
        self.pending_edges.push((edge, from.to_string(), to.to_string()));
        self
    }

    /// Build the network and validate it (single point of error handling)
    pub fn build(mut self) -> Result<Network<T>> {
        // Add all pending nodes to the network first
        for node in self.pending_nodes {
            self.network.add_node(node)?;
        }

        // Then add all pending edges to the network
        for (edge, from_id, to_id) in self.pending_edges {
            self.network.add_edge(edge, &from_id, &to_id)?;
        }

        // Perform comprehensive validation at the end
        self.network.validate()?;
        Ok(self.network)
    }

    /// Build the network without validation (for testing)
    pub fn build_unchecked(mut self) -> Network<T> {
        // Add all pending nodes to the network
        for node in self.pending_nodes {
            self.network.add_node(node).unwrap(); // Unwrap because we are building unchecked
        }

        // Add all pending edges to the network
        for (edge, from_id, to_id) in self.pending_edges {
            self.network.add_edge(edge, &from_id, &to_id).unwrap(); // Unwrap because we are building unchecked
        }

        self.network
    }
}

impl<T: RealField + FromPrimitive + num_traits::Float> Default for NetworkBuilder<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use nalgebra::Vector3;

    #[test]
    fn test_node_creation() {
        let position = Vector3::new(0.0, 0.0, 0.0);
        let node = Node::junction("j1".to_string(), position);

        assert_eq!(node.id, "j1");
        assert_eq!(node.node_type, NodeType::Junction);
        assert!(!node.has_boundary_condition());
    }

    #[test]
    fn test_inlet_node_with_pressure() {
        let position = Vector3::new(0.0, 0.0, 0.0);
        let node = Node::inlet_pressure("inlet1".to_string(), position, 1000.0);

        assert_eq!(node.node_type, NodeType::Inlet);
        assert!(node.has_boundary_condition());

        if let Some(bc) = node.boundary_condition() {
            if let Some(p) = bc.pressure_value() {
                assert_relative_eq!(p, 1000.0, epsilon = 1e-10);
            } else {
                panic!("Expected pressure boundary condition");
            }
        } else {
            panic!("Expected boundary condition");
        }
    }

    #[test]
    fn test_edge_creation() {
        let edge = Edge::channel("ch1".to_string(), 100.0, 0.001, 1e-6);

        assert_eq!(edge.id, "ch1");
        assert_eq!(edge.edge_type, EdgeType::Channel);
        assert_relative_eq!(edge.resistance, 100.0, epsilon = 1e-10);
        assert!(!edge.is_active_component());
    }

    #[test]
    fn test_pump_edge() {
        let edge = Edge::pump("pump1".to_string(), 500.0);

        assert_eq!(edge.edge_type, EdgeType::Pump);
        assert!(edge.is_active_component());
        assert_eq!(edge.pressure_rise(), Some(500.0));
    }

    #[test]
    fn test_valve_resistance() {
        let mut edge = Edge::valve("valve1".to_string(), 100.0, 0.5);

        // Half open valve should have 2x resistance
        assert_relative_eq!(edge.effective_resistance(), 200.0, epsilon = 1e-10);

        // Closed valve should have very high resistance
        edge.properties.component_params.insert("opening".to_string(), 0.0);
        assert!(edge.effective_resistance() > 1e10);

        // Fully open valve should have base resistance
        edge.properties.component_params.insert("opening".to_string(), 1.0);
        assert_relative_eq!(edge.effective_resistance(), 100.0, epsilon = 1e-10);
    }

    #[test]
    fn test_network_builder_simple() {
        let network = NetworkBuilder::new()
            .add_inlet_pressure("inlet", 0.0, 0.0, 1000.0)
            .add_outlet_pressure("outlet", 1.0, 0.0, 0.0)
            .add_channel("ch1", "inlet", "outlet", ChannelProperties::new(100.0, 0.001, 1e-6))
            .build().unwrap();

        assert_eq!(network.node_count(), 2);
        assert_eq!(network.edge_count(), 1);
        assert!(network.is_connected());

        // Check nodes exist
        assert!(network.get_node("inlet").is_some());
        assert!(network.get_node("outlet").is_some());

        // Check edge exists
        assert!(network.get_edge("ch1").is_some());
    }

    #[test]
    fn test_network_builder_t_junction() {
        let network = NetworkBuilder::new()
            .add_inlet_pressure("inlet", 0.0, 0.0, 1000.0)
            .add_junction("junction", 0.5, 0.0)
            .add_outlet_pressure("outlet1", 1.0, 0.5, 0.0)
            .add_outlet_pressure("outlet2", 1.0, -0.5, 0.0)
            .add_channel("ch1", "inlet", "junction", ChannelProperties::new(50.0, 0.0005, 1e-6))
            .add_channel("ch2", "junction", "outlet1", ChannelProperties::new(100.0, 0.0005, 5e-7))
            .add_channel("ch3", "junction", "outlet2", ChannelProperties::new(100.0, 0.0005, 5e-7))
            .build().unwrap();

        assert_eq!(network.node_count(), 4);
        assert_eq!(network.edge_count(), 3);
        assert!(network.is_connected());

        // Check junction has correct neighbors
        let neighbors = network.neighbors("junction").unwrap();
        assert!(neighbors.len() >= 2); // Connected to at least inlet and one outlet
    }

    #[test]
    fn test_network_with_pump() {
        let network = NetworkBuilder::new()
            .add_inlet_pressure("inlet", 0.0, 0.0, 0.0)
            .add_outlet_pressure("outlet", 1.0, 0.0, 0.0)
            .add_pump("pump1", "inlet", "outlet", 1000.0)
            .build().unwrap();

        let pump = network.get_edge("pump1").unwrap();
        assert!(pump.is_active_component());
        assert_eq!(pump.pressure_rise(), Some(1000.0));
    }

    #[test]
    fn test_network_validation_isolated_node() {
        let result = NetworkBuilder::new()
            .add_inlet_pressure("inlet", 0.0, 0.0, 1000.0)
            .add_outlet_pressure("outlet", 1.0, 0.0, 0.0)
            .add_junction("isolated", 0.5, 0.5) // Not connected to anything
            .build();

        assert!(result.is_err());
    }

    #[test]
    fn test_network_validation_no_boundary_conditions() {
        let result = NetworkBuilder::new()
            .add_junction("j1", 0.0, 0.0)
            .add_junction("j2", 1.0, 0.0)
            .add_channel("ch1", "j1", "j2", ChannelProperties::new(100.0, 0.001, 1e-6))
            .build();

        assert!(result.is_err());
    }

    #[test]
    fn test_network_validation_disconnected() {
        let result = NetworkBuilder::new()
            .add_inlet_pressure("inlet1", 0.0, 0.0, 1000.0)
            .add_outlet_pressure("outlet1", 0.5, 0.0, 0.0)
            .add_inlet_pressure("inlet2", 1.0, 0.0, 1000.0)
            .add_outlet_pressure("outlet2", 1.5, 0.0, 0.0)
            .add_channel("ch1", "inlet1", "outlet1", ChannelProperties::new(100.0, 0.001, 1e-6))
            .add_channel("ch2", "inlet2", "outlet2", ChannelProperties::new(100.0, 0.001, 1e-6))
            .build();

        assert!(result.is_err());
    }

    #[test]
    fn test_node_edges() {
        let network = NetworkBuilder::new()
            .add_junction("junction", 0.0, 0.0)
            .add_inlet_pressure("inlet", -1.0, 0.0, 1000.0)
            .add_outlet_pressure("outlet1", 1.0, 1.0, 0.0)
            .add_outlet_pressure("outlet2", 1.0, -1.0, 0.0)
            .add_channel("ch1", "inlet", "junction", ChannelProperties::new(50.0, 0.001, 1e-6))
            .add_channel("ch2", "junction", "outlet1", ChannelProperties::new(100.0, 0.0005, 5e-7))
            .add_channel("ch3", "junction", "outlet2", ChannelProperties::new(100.0, 0.0005, 5e-7))
            .build_unchecked(); // Use unchecked to avoid validation

        let junction_edges = network.node_edges("junction").unwrap();
        assert!(junction_edges.len() >= 2); // Should have at least 2 edges

        let inlet_edges = network.node_edges("inlet").unwrap();
        assert_eq!(inlet_edges.len(), 1);
    }

    #[test]
    fn test_clear_solution() {
        let mut network = NetworkBuilder::new()
            .add_inlet_pressure("inlet", 0.0, 0.0, 1000.0)
            .add_outlet_pressure("outlet", 1.0, 0.0, 0.0)
            .add_channel("ch1", "inlet", "outlet", ChannelProperties::new(100.0, 0.001, 1e-6))
            .build().unwrap();

        // Set some solution data
        if let Some(node) = network.get_node_mut("inlet") {
            node.pressure = Some(1000.0);
            node.flow_rate = Some(0.001);
        }

        if let Some(edge) = network.get_edge_mut("ch1") {
            edge.flow_rate = Some(0.001);
            edge.pressure_drop = Some(100.0);
        }

        // Clear solution
        network.clear_solution();

        // Check that solution data is cleared
        assert!(network.get_node("inlet").unwrap().pressure.is_none());
        assert!(network.get_node("inlet").unwrap().flow_rate.is_none());
        assert!(network.get_edge("ch1").unwrap().flow_rate.is_none());
        assert!(network.get_edge("ch1").unwrap().pressure_drop.is_none());
    }

    #[test]
    fn test_fluid_properties() {
        let mut network = NetworkBuilder::new()
            .with_fluid(Fluid::air())
            .add_inlet_pressure("inlet", 0.0, 0.0, 1000.0)
            .add_outlet_pressure("outlet", 1.0, 0.0, 0.0)
            .add_channel("ch1", "inlet", "outlet", ChannelProperties::new(100.0, 0.001, 1e-6))
            .build().unwrap();

        // Check initial fluid
        assert!(network.fluid().name.contains("Air"));

        // Change fluid
        network.set_fluid(Fluid::water());
        assert!(network.fluid().name.contains("Water"));
    }
}
