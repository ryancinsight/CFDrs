use super::{ChannelSpec, EdgeKind, NodeKind, NodeSpec};
use std::fmt;

#[derive(Debug, Clone)]
pub struct NetworkBlueprint {
    pub name: String,
    pub nodes: Vec<NodeSpec>,
    pub channels: Vec<ChannelSpec>,
}

impl NetworkBlueprint {
    #[must_use]
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            nodes: Vec::new(),
            channels: Vec::new(),
        }
    }

    pub fn add_node(&mut self, node: NodeSpec) {
        self.nodes.push(node);
    }

    pub fn add_channel(&mut self, channel: ChannelSpec) {
        self.channels.push(channel);
    }

    /// Count nodes by kind.
    #[must_use]
    pub fn inlet_count(&self) -> usize {
        self.nodes.iter().filter(|n| matches!(n.kind, NodeKind::Inlet)).count()
    }

    /// Count outlet nodes.
    #[must_use]
    pub fn outlet_count(&self) -> usize {
        self.nodes.iter().filter(|n| matches!(n.kind, NodeKind::Outlet)).count()
    }

    /// Count junction nodes.
    #[must_use]
    pub fn junction_count(&self) -> usize {
        self.nodes.iter().filter(|n| matches!(n.kind, NodeKind::Junction)).count()
    }

    /// Count pipe channels.
    #[must_use]
    pub fn pipe_count(&self) -> usize {
        self.channels.iter().filter(|c| matches!(c.kind, EdgeKind::Pipe)).count()
    }

    /// Total path length across all channels [m].
    #[must_use]
    pub fn total_length_m(&self) -> f64 {
        self.channels.iter().map(|c| c.length_m).sum()
    }

    /// Return a human-readable summary of this blueprint (replaces `ConversionSummary`).
    #[must_use]
    pub fn describe(&self) -> String {
        format!(
            "NetworkBlueprint '{}'\n\
             ─────────────────────────────────────\n\
             Nodes   : {} total  ({} inlets, {} outlets, {} junctions)\n\
             Channels: {}  ({} pipes)\n\
             Length  : {:.6} m total\n",
            self.name,
            self.nodes.len(),
            self.inlet_count(),
            self.outlet_count(),
            self.junction_count(),
            self.channels.len(),
            self.pipe_count(),
            self.total_length_m(),
        )
    }
}

impl fmt::Display for NetworkBlueprint {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.describe())
    }
}

