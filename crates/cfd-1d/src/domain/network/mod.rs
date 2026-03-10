//! Network topology for 1D CFD simulations
//!
//! Provides graph-based representation of fluid networks

mod builder;
mod blueprint_validation;
mod component_type;
mod edge;
mod graph;
mod junction_losses;
mod metadata;
mod node;
pub mod sink;
mod wrapper;

pub use blueprint_validation::validate_blueprint_for_1d_solve;
pub use builder::{network_from_blueprint, NetworkBuilder};
pub use cfd_schematics::domain::model::EdgeKind as EdgeType;
pub use cfd_schematics::domain::model::NodeKind as NodeType;
pub use component_type::ComponentType;
pub use edge::{ChannelProperties, Edge};
pub use graph::{NetworkGraph, NetworkGraphExt};
pub use metadata::NetworkMetadata;
pub use node::{Node, NodeProperties};
pub use sink::NetworkBuilderSink;
pub use wrapper::{EdgeProperties, EdgeWithProperties, Network, ParallelEdge};

// Re-export boundary conditions from core
pub use cfd_core::physics::boundary::BoundaryCondition;
