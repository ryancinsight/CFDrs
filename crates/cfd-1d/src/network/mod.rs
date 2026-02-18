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
pub use edge::{ChannelProperties, Edge};
pub use graph::{NetworkGraph, NetworkGraphExt};
pub use metadata::NetworkMetadata;
pub use node::{Node, NodeProperties};
pub use wrapper::{EdgeProperties, EdgeWithProperties, Network, ParallelEdge};
pub use scheme::domain::model::NodeKind as NodeType;
pub use scheme::domain::model::EdgeKind as EdgeType;

// Re-export boundary conditions from core
pub use cfd_core::physics::boundary::BoundaryCondition;
