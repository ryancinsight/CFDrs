//! Network topology for 1D CFD simulations
//!
//! Provides graph-based representation of fluid networks

mod builder;
mod edge;
mod graph;
mod metadata;
mod node;

pub use builder::NetworkBuilder;
pub use edge::{ChannelProperties, Edge, EdgeProperties, EdgeType};
pub use graph::NetworkGraph;
pub use metadata::NetworkMetadata;
pub use node::{Node, NodeProperties, NodeType};

// Re-export boundary conditions from core
pub use cfd_core::boundary::BoundaryCondition;

// Type alias for the directed graph
pub type Network<T> = petgraph::graph::Graph<Node<T>, Edge<T>, petgraph::Directed>;
