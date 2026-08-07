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
    /// Construct an empty `NetworkBlueprint` with auto-laid-out nodes.
    ///
    /// # Deprecated
    /// Nodes created via this constructor have no layout positions
    /// and fall back to the generic auto-layout, producing incorrect
    /// geometry for bifurcation trees. Use `create_geometry()` from
    /// `cfd_schematics::geometry::generator` for split-tree layouts,
    /// or ensure every `NodeSpec` is built with `NodeSpec::new_at()`.
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

    /// Construct an empty `NetworkBlueprint` that will not run
    /// auto-layout, leaving it the caller's responsibility to position
    /// every node explicitly through `NodeSpec::new_at()`.
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

    /// Attach render hints that visualization passes use to control
    /// annotation density, line styling, and schematic polish.
    /// Returns `self` for chaining.
    #[must_use]
    pub fn with_render_hints(mut self, hints: BlueprintRenderHints) -> Self {
        self.render_hints = Some(hints);
        self
    }

    /// Resolve the blueprint's render hints, preferring the inline
    /// field and falling back to the metadata container when absent.
    #[must_use]
    pub fn render_hints(&self) -> Option<&BlueprintRenderHints> {
        self.render_hints
            .as_ref()
            .or_else(|| self.metadata.as_ref()?.get::<BlueprintRenderHints>())
    }

    /// Attach a metadata entry of type `T` and return `self` for chaining.
    #[must_use]
    pub fn with_metadata<T: crate::geometry::metadata::Metadata + Clone + 'static>(
        mut self,
        metadata: T,
    ) -> Self {
        self.insert_metadata(metadata);
        self
    }

    /// Insert a metadata entry of type `T` into the blueprint's metadata
    /// container, lazy-creating the container on first use.
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

    /// Borrow a metadata entry of type `T` previously inserted into the
    /// blueprint, if present.
    #[must_use]
    pub fn metadata<T: crate::geometry::metadata::Metadata + 'static>(&self) -> Option<&T> {
        self.metadata.as_ref()?.get::<T>()
    }

    /// Borrow the geometry-authoring provenance metadata, if present.
    #[must_use]
    pub fn geometry_authoring_provenance(&self) -> Option<&GeometryAuthoringProvenance> {
        self.metadata::<GeometryAuthoringProvenance>()
    }

    /// Whether geometry has been authored for this blueprint, either by
    /// setting the explicit flag or by carrying an authoring-provenance
    /// record in the metadata.
    #[must_use]
    pub fn is_geometry_authored(&self) -> bool {
        self.geometry_authored || self.geometry_authoring_provenance().is_some()
    }

    /// Attach an authored topology specification and return `self` for
    /// chaining.
    #[must_use]
    pub fn with_topology_spec(mut self, topology: BlueprintTopologySpec) -> Self {
        self.topology = Some(topology);
        self
    }

    /// Attach a topology-lineage metadata block and return `self` for
    /// chaining.
    #[must_use]
    pub fn with_lineage(mut self, lineage: TopologyLineageMetadata) -> Self {
        self.lineage = Some(lineage);
        self
    }

    /// Borrow the topology specification, if one is attached.
    #[must_use]
    pub fn topology_spec(&self) -> Option<&BlueprintTopologySpec> {
        self.topology.as_ref()
    }

    /// Borrow the topology-lineage metadata, if one is attached.
    #[must_use]
    pub fn lineage(&self) -> Option<&TopologyLineageMetadata> {
        self.lineage.as_ref()
    }

    /// Resolve the channel identifiers that lie on the treatment path,
    /// delegating to the attached `BlueprintTopologySpec` when present
    /// and returning an empty vector otherwise.
    #[must_use]
    pub fn treatment_channel_ids(&self) -> Vec<String> {
        self.topology
            .as_ref()
            .map_or_else(Vec::new, BlueprintTopologySpec::treatment_channel_ids)
    }

    /// Serialize the blueprint to a pretty-printed JSON string.
    pub fn to_json_pretty(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(self)
    }

    /// Serialize the blueprint to a JSON string, equivalent to `to_json_pretty`.
    pub fn to_json(&self) -> Result<String, serde_json::Error> {
        self.to_json_pretty()
    }

    /// Parse a `NetworkBlueprint` from a JSON string.
    pub fn from_json(json: &str) -> Result<Self, serde_json::Error> {
        serde_json::from_str(json)
    }

    /// Append a node to the blueprint's node list.
    pub fn add_node(&mut self, node: NodeSpec) {
        self.nodes.push(node);
    }

    /// Append a channel to the blueprint's channel list.
    pub fn add_channel(&mut self, channel: ChannelSpec) {
        self.channels.push(channel);
    }
}
