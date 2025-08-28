//! Network graph operations

use super::{Edge, Node};
use nalgebra::RealField;
use petgraph::graph::{EdgeIndex, Graph, NodeIndex};
use petgraph::Directed;

/// Type alias for the network graph
pub type NetworkGraph<T> = Graph<Node<T>, Edge<T>, Directed>;

/// Extension trait for network graph operations
pub trait NetworkGraphExt<T: RealField + Copy> {
    /// Find a node by ID
    fn find_node_by_id(&self, id: &str) -> Option<NodeIndex>;

    /// Find an edge by ID
    fn find_edge_by_id(&self, id: &str) -> Option<EdgeIndex>;

    /// Get all inlet nodes
    fn inlet_nodes(&self) -> Vec<NodeIndex>;

    /// Get all outlet nodes
    fn outlet_nodes(&self) -> Vec<NodeIndex>;
}

impl<T: RealField + Copy> NetworkGraphExt<T> for NetworkGraph<T> {
    fn find_node_by_id(&self, id: &str) -> Option<NodeIndex> {
        self.node_indices().find(|&idx| self[idx].id == id)
    }

    fn find_edge_by_id(&self, id: &str) -> Option<EdgeIndex> {
        self.edge_indices().find(|&idx| self[idx].id == id)
    }

    fn inlet_nodes(&self) -> Vec<NodeIndex> {
        use super::NodeType;
        self.node_indices()
            .filter(|&idx| self[idx].node_type == NodeType::Inlet)
            .collect()
    }

    fn outlet_nodes(&self) -> Vec<NodeIndex> {
        use super::NodeType;
        self.node_indices()
            .filter(|&idx| self[idx].node_type == NodeType::Outlet)
            .collect()
    }
}
