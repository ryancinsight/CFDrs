//! Network wrapper with convenience methods

use super::{Edge, NetworkGraph, Node};
use crate::channel::ChannelGeometry;
use cfd_core::fluid::Fluid;
use nalgebra::RealField;
use petgraph::graph::{EdgeIndex, NodeIndex};
use petgraph::visit::EdgeRef;
use std::collections::HashMap;

/// Extended network with fluid properties and convenience methods
pub struct Network<T: RealField + Copy> {
    /// The underlying graph
    pub graph: NetworkGraph<T>,
    /// Fluid properties for the network
    pub fluid: Fluid<T>,
    /// Node pressures
    pub pressures: HashMap<NodeIndex, T>,
    /// Edge flow rates
    pub flow_rates: HashMap<EdgeIndex, T>,
    /// Edge properties
    pub properties: HashMap<EdgeIndex, EdgeProperties<T>>,
}

/// Properties for edges in the network
#[derive(Debug, Clone)]
pub struct EdgeProperties<T: RealField + Copy> {
    /// Edge identifier
    pub id: String,
    /// Length of the channel
    pub length: T,
    /// Cross-sectional area
    pub area: T,
    /// Hydraulic diameter
    pub hydraulic_diameter: Option<T>,
    /// Resistance coefficient
    pub resistance: T,
    /// Channel geometry if applicable
    pub geometry: Option<ChannelGeometry<T>>,
    /// Additional properties
    pub properties: HashMap<String, T>,
}

/// Edge with properties for iteration
pub struct EdgeWithProperties<'a, T: RealField + Copy> {
    /// Edge identifier
    pub id: String,
    /// Flow rate
    pub flow_rate: T,
    /// Node indices (from, to)
    pub nodes: (NodeIndex, NodeIndex),
    /// Properties
    pub properties: &'a EdgeProperties<T>,
}

impl<T: RealField + Copy> Network<T> {
    /// Create a new network
    pub fn new(graph: NetworkGraph<T>, fluid: Fluid<T>) -> Self {
        Self {
            graph,
            fluid,
            pressures: HashMap::new(),
            flow_rates: HashMap::new(),
            properties: HashMap::new(),
        }
    }

    /// Get the fluid properties
    pub fn fluid(&self) -> &Fluid<T> {
        &self.fluid
    }

    /// Get edges with their properties
    pub fn edges_with_properties(&self) -> Vec<EdgeWithProperties<T>> {
        self.graph
            .edge_references()
            .filter_map(|edge_ref| {
                let edge_idx = edge_ref.id();
                let edge_data = edge_ref.weight();
                let (from, to) = (edge_ref.source(), edge_ref.target());
                self.properties
                    .get(&edge_idx)
                    .map(|props| EdgeWithProperties {
                        id: edge_data.id.clone(),
                        flow_rate: *self.flow_rates.get(&edge_idx).unwrap_or(&T::zero()),
                        nodes: (from, to),
                        properties: props,
                    })
            })
            .collect()
    }

    /// Get all nodes
    pub fn nodes(&self) -> impl Iterator<Item = &Node<T>> {
        self.graph.node_weights()
    }

    /// Get node pressures
    pub fn pressures(&self) -> &HashMap<NodeIndex, T> {
        &self.pressures
    }

    /// Get edge flow rates
    pub fn flow_rates(&self) -> &HashMap<EdgeIndex, T> {
        &self.flow_rates
    }

    /// Set pressure at a node
    pub fn set_pressure(&mut self, node: NodeIndex, pressure: T) {
        self.pressures.insert(node, pressure);
    }

    /// Set flow rate for an edge
    pub fn set_flow_rate(&mut self, edge: EdgeIndex, flow_rate: T) {
        self.flow_rates.insert(edge, flow_rate);
    }

    /// Add properties for an edge
    pub fn add_edge_properties(&mut self, edge: EdgeIndex, properties: EdgeProperties<T>) {
        self.properties.insert(edge, properties);
    }
}
