use super::{ChannelSpec, NodeSpec};

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
}
