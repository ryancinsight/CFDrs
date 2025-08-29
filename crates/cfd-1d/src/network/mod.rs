//! Network topology for 1D CFD simulations
//!
//! Provides graph-based representation of fluid networks

mod builder;
mod component_type;
mod edge;
mod graph;
mod metadata;
mod node;
mod wrapper;

pub use builder::NetworkBuilder;
pub use component_type::ComponentType;
pub use edge::{ChannelProperties, Edge, EdgeType};
pub use graph::NetworkGraph;
pub use metadata::NetworkMetadata;
pub use node::{Node, NodeProperties, NodeType};
pub use wrapper::{EdgeProperties, EdgeWithProperties, Network, ParallelEdge};

// Re-export boundary conditions from core
pub use cfd_core::boundary::BoundaryCondition;
