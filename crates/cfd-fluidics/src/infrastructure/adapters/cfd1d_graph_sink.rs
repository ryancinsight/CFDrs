use crate::application::ports::GraphSink;
use crate::domain::model::{EdgeKind, NetworkBlueprint, NodeKind};
use cfd_1d::{Edge, EdgeType, NetworkBuilder, NetworkGraph, Node, NodeType};
use cfd_core::error::{Error, Result};
use petgraph::graph::NodeIndex;
use std::collections::HashMap;

pub struct Cfd1dGraphSink;

impl GraphSink for Cfd1dGraphSink {
    type Output = NetworkGraph<f64>;

    fn build(&self, blueprint: &NetworkBlueprint) -> Result<Self::Output> {
        let mut builder = NetworkBuilder::<f64>::new();
        let mut node_indices: HashMap<String, NodeIndex> = HashMap::new();

        for node in &blueprint.nodes {
            let node_type = match node.kind {
                NodeKind::Inlet => NodeType::Inlet,
                NodeKind::Outlet => NodeType::Outlet,
                NodeKind::Junction => NodeType::Junction,
            };
            let index = builder.add_node(Node::new(node.id.as_str().to_string(), node_type));
            node_indices.insert(node.id.as_str().to_string(), index);
        }

        for channel in &blueprint.channels {
            let from = node_indices
                .get(channel.from.as_str())
                .copied()
                .ok_or_else(|| {
                    Error::InvalidConfiguration(format!(
                        "Unknown source node '{}' in channel '{}'",
                        channel.from.as_str(),
                        channel.id.as_str()
                    ))
                })?;

            let to = node_indices
                .get(channel.to.as_str())
                .copied()
                .ok_or_else(|| {
                    Error::InvalidConfiguration(format!(
                        "Unknown target node '{}' in channel '{}'",
                        channel.to.as_str(),
                        channel.id.as_str()
                    ))
                })?;

            let edge_type = match channel.kind {
                EdgeKind::Pipe => EdgeType::Pipe,
            };

            let mut edge = Edge::new(channel.id.as_str().to_string(), edge_type);
            edge.resistance = channel.resistance;
            edge.quad_coeff = channel.quad_coeff;

            builder.add_edge(from, to, edge);
        }

        builder.build()
    }
}
