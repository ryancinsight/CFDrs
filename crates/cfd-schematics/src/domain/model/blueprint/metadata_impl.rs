use super::super::{ChannelSpec, NodeSpec};
use super::NetworkBlueprint;
use crate::geometry::metadata::{
    BlueprintRenderHints, GeometryAuthoringProvenance, MetadataContainer,
};
use crate::topology::{BlueprintTopologySpec, TopologyLineageMetadata};

impl NetworkBlueprint {
    #[deprecated(
        note = "Nodes created via NetworkBlueprint::new() have no layout positions and fall \
                back to the generic auto-layout, producing incorrect geometry for bifurcation \
                trees. Use `create_geometry()` (cfd_schematics::geometry::generator) for \
                split-tree layouts, or ensure every NodeSpec is created with NodeSpec::new_at()."
    )]
    #[must_use]
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            box_dims: (127.76, 85.47),
            box_outline: Vec::new(),
            nodes: Vec::new(),
            channels: Vec::new(),
            render_hints: None,
            topology: None,
            lineage: None,
            metadata: None,
            geometry_authored: false,
        }
    }

    #[must_use]
    pub fn new_with_explicit_positions(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            box_dims: (127.76, 85.47),
            box_outline: Vec::new(),
            nodes: Vec::new(),
            channels: Vec::new(),
            render_hints: None,
            topology: None,
            lineage: None,
            metadata: None,
            geometry_authored: false,
        }
    }

    #[must_use]
    pub fn with_render_hints(mut self, hints: BlueprintRenderHints) -> Self {
        self.render_hints = Some(hints);
        self
    }

    #[must_use]
    pub fn render_hints(&self) -> Option<&BlueprintRenderHints> {
        self.render_hints
            .as_ref()
            .or_else(|| self.metadata.as_ref()?.get::<BlueprintRenderHints>())
    }

    #[must_use]
    pub fn with_metadata<T: crate::geometry::metadata::Metadata + Clone + 'static>(
        mut self,
        metadata: T,
    ) -> Self {
        self.insert_metadata(metadata);
        self
    }

    pub fn insert_metadata<T: crate::geometry::metadata::Metadata + Clone + 'static>(
        &mut self,
        metadata: T,
    ) {
        if self.metadata.is_none() {
            self.metadata = Some(MetadataContainer::new());
        }
        self.metadata
            .as_mut()
            .expect("metadata container must exist")
            .insert(metadata);
    }

    #[must_use]
    pub fn metadata<T: crate::geometry::metadata::Metadata + 'static>(&self) -> Option<&T> {
        self.metadata.as_ref()?.get::<T>()
    }

    #[must_use]
    pub fn geometry_authoring_provenance(&self) -> Option<&GeometryAuthoringProvenance> {
        self.metadata::<GeometryAuthoringProvenance>()
    }

    #[must_use]
    pub fn is_geometry_authored(&self) -> bool {
        self.geometry_authored || self.geometry_authoring_provenance().is_some()
    }

    #[must_use]
    pub fn with_topology_spec(mut self, topology: BlueprintTopologySpec) -> Self {
        self.topology = Some(topology);
        self
    }

    #[must_use]
    pub fn with_lineage(mut self, lineage: TopologyLineageMetadata) -> Self {
        self.lineage = Some(lineage);
        self
    }

    #[must_use]
    pub fn topology_spec(&self) -> Option<&BlueprintTopologySpec> {
        self.topology.as_ref()
    }

    #[must_use]
    pub fn lineage(&self) -> Option<&TopologyLineageMetadata> {
        self.lineage.as_ref()
    }

    #[must_use]
    pub fn treatment_channel_ids(&self) -> Vec<String> {
        self.topology
            .as_ref()
            .map_or_else(Vec::new, BlueprintTopologySpec::treatment_channel_ids)
    }

    pub fn to_json_pretty(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(self)
    }

    pub fn to_json(&self) -> Result<String, serde_json::Error> {
        self.to_json_pretty()
    }

    pub fn from_json(json: &str) -> Result<Self, serde_json::Error> {
        serde_json::from_str(json)
    }

    pub fn add_node(&mut self, node: NodeSpec) {
        self.nodes.push(node);
    }

    pub fn add_channel(&mut self, channel: ChannelSpec) {
        self.channels.push(channel);
    }
}
