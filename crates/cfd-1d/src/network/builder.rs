//! Network builder for constructing fluid networks

use super::{Edge, EdgeType, NetworkGraph, Node, NodeType};
use cfd_core::error::Result;
use nalgebra::RealField;

/// Builder for constructing network graphs
pub struct NetworkBuilder<T: RealField + Copy> {
    graph: NetworkGraph<T>,
}

impl<T: RealField + Copy> NetworkBuilder<T> {
    /// Create a new network builder
    #[must_use]
    pub fn new() -> Self {
        Self {
            graph: NetworkGraph::new(),
        }
    }

    /// Add a node to the network
    pub fn add_node(&mut self, node: Node<T>) -> petgraph::graph::NodeIndex {
        self.graph.add_node(node)
    }

    /// Add an edge between two nodes
    pub fn add_edge(
        &mut self,
        from: petgraph::graph::NodeIndex,
        to: petgraph::graph::NodeIndex,
        edge: Edge<T>,
    ) -> petgraph::graph::EdgeIndex {
        self.graph.add_edge(from, to, edge)
    }

    /// Create an inlet node
    pub fn add_inlet(&mut self, id: String) -> petgraph::graph::NodeIndex {
        self.add_node(Node::new(id, NodeType::Inlet))
    }

    /// Create an outlet node
    pub fn add_outlet(&mut self, id: String) -> petgraph::graph::NodeIndex {
        self.add_node(Node::new(id, NodeType::Outlet))
    }

    /// Create a junction node
    pub fn add_junction(&mut self, id: String) -> petgraph::graph::NodeIndex {
        self.add_node(Node::new(id, NodeType::Junction))
    }

    /// Connect two nodes with a pipe
    pub fn connect_with_pipe(
        &mut self,
        from: petgraph::graph::NodeIndex,
        to: petgraph::graph::NodeIndex,
        id: String,
    ) -> petgraph::graph::EdgeIndex {
        self.add_edge(from, to, Edge::new(id, EdgeType::Pipe))
    }

    /// Build the final network with validation
    pub fn build(self) -> Result<NetworkGraph<T>> {
        // Validate network connectivity
        if self.graph.node_count() == 0 {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "Network has no nodes".to_string(),
            ));
        }

        // Ensure at least one inlet and outlet
        let has_inlet = self
            .graph
            .node_weights()
            .any(|n| matches!(n.node_type, NodeType::Inlet));
        let has_outlet = self
            .graph
            .node_weights()
            .any(|n| matches!(n.node_type, NodeType::Outlet));

        if !has_inlet {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "Network requires at least one inlet".to_string(),
            ));
        }

        if !has_outlet {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "Network requires at least one outlet".to_string(),
            ));
        }

        // Check for disconnected components
        use petgraph::algo::connected_components;
        if connected_components(&self.graph) > 1 {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "Network has disconnected components".to_string(),
            ));
        }

        Ok(self.graph)
    }
}

impl<T: RealField + Copy> Default for NetworkBuilder<T> {
    fn default() -> Self {
        Self::new()
    }
}
