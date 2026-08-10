use crate::application::ports::GraphSink;
use crate::application::use_cases::NetworkGenerationService;
use crate::domain::model::{ChannelSpec, NetworkBlueprint, NodeSpec};
use cfd_core::error::{Error, Result};
use petgraph::graph::{DiGraph, NodeIndex};
use std::collections::HashMap;

/// Petgraph adjacency-graph representation of a network blueprint.
pub type DesignGraph = DiGraph<NodeSpec, ChannelSpec>;

/// `GraphSink` adapter that builds a [`DesignGraph`] from a blueprint.
#[derive(Debug, Clone, Copy, Default)]
pub struct PetgraphGraphSink;

impl GraphSink for PetgraphGraphSink {
    type Output = DesignGraph;

    fn build(&self, blueprint: &NetworkBlueprint) -> Result<Self::Output> {
        let mut graph = DesignGraph::new();
        let mut node_map: HashMap<String, NodeIndex> = HashMap::new();

        for node in &blueprint.nodes {
            let index = graph.add_node(node.clone());
            node_map.insert(node.id.as_str().to_string(), index);
        }

        for channel in &blueprint.channels {
            let from = node_map.get(channel.from.as_str()).ok_or_else(|| {
                Error::InvalidConfiguration(format!(
                    "Channel '{}' references missing source node '{}'",
                    channel.id.as_str(),
                    channel.from.as_str()
                ))
            })?;
            let to = node_map.get(channel.to.as_str()).ok_or_else(|| {
                Error::InvalidConfiguration(format!(
                    "Channel '{}' references missing target node '{}'",
                    channel.id.as_str(),
                    channel.to.as_str()
                ))
            })?;

            graph.add_edge(*from, *to, channel.clone());
        }

        Ok(graph)
    }
}

/// Build a design graph from a blueprint through the default graph sink.
pub fn build_design_graph(blueprint: &NetworkBlueprint) -> Result<DesignGraph> {
    NetworkGenerationService::new(PetgraphGraphSink).generate(blueprint)
}
