use super::{ChannelSpec, EdgeKind, NodeKind, NodeSpec};
use crate::geometry::metadata::{
    BlueprintRenderHints, GeometryAuthoringProvenance, MetadataContainer,
};
use crate::topology::{BlueprintTopologyFactory, BlueprintTopologySpec, TopologyLineageMetadata};
use serde::{Deserialize, Serialize};
use std::fmt;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NetworkBlueprint {
    pub name: String,
    pub box_dims: (f64, f64),
    pub box_outline: Vec<((f64, f64), (f64, f64))>,
    pub nodes: Vec<NodeSpec>,
    pub channels: Vec<ChannelSpec>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub render_hints: Option<BlueprintRenderHints>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub topology: Option<BlueprintTopologySpec>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lineage: Option<TopologyLineageMetadata>,
    /// Optional blueprint-level metadata (e.g. rendering hints written by
    /// `cfd-optim` so that a generic renderer can produce fully-annotated
    /// schematics without any domain-specific render logic).
    #[serde(skip)]
    pub metadata: Option<MetadataContainer>,
}

impl NetworkBlueprint {
    /// Create an empty `NetworkBlueprint` with default ANSI/SLAS well-plate dimensions.
    ///
    /// # Warning
    ///
    /// Nodes added via this constructor have **no layout positions** (`point = (0.0, 0.0)`)
    /// and the renderer will fall back to a generic topological-sort column layout
    /// (`auto_layout_positions`), which produces incorrect geometry for bifurcation trees
    /// (diagonal lines, diamond shapes).
    ///
    /// For split-tree / bifurcation layouts, use
    /// [`create_geometry`](crate::geometry::generator::create_geometry) instead — it
    /// computes proper geometric coordinates for every node before returning a
    /// `NetworkBlueprint`.
    ///
    /// This constructor is appropriate only for topologies where all
    /// [`NodeSpec`](crate::domain::model::NodeSpec)s are created with explicit positions
    /// via [`NodeSpec::new_at`](crate::domain::model::NodeSpec::new_at).
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
            box_dims: (127.76, 85.47), // Default to ANSI/SLAS well plate, overrideable later
            box_outline: Vec::new(),
            nodes: Vec::new(),
            channels: Vec::new(),
            render_hints: None,
            topology: None,
            lineage: None,
            metadata: None,
        }
    }

    /// Attach [`BlueprintRenderHints`] and return `self` (builder pattern).
    ///
    /// Hints are subsequently readable via [`Self::render_hints`].
    #[must_use]
    pub fn with_render_hints(mut self, hints: BlueprintRenderHints) -> Self {
        self.render_hints = Some(hints);
        self
    }

    /// Return the [`BlueprintRenderHints`] embedded in this blueprint, if any.
    #[must_use]
    pub fn render_hints(&self) -> Option<&BlueprintRenderHints> {
        self.render_hints
            .as_ref()
            .or_else(|| self.metadata.as_ref()?.get::<BlueprintRenderHints>())
    }

    /// Attach blueprint-level metadata and return `self` (builder pattern).
    #[must_use]
    pub fn with_metadata<T: crate::geometry::metadata::Metadata + Clone + 'static>(
        mut self,
        metadata: T,
    ) -> Self {
        self.insert_metadata(metadata);
        self
    }

    /// Insert blueprint-level metadata in-place.
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

    /// Read blueprint-level metadata by type.
    #[must_use]
    pub fn metadata<T: crate::geometry::metadata::Metadata + 'static>(&self) -> Option<&T> {
        self.metadata.as_ref()?.get::<T>()
    }

    /// Geometry-authoring provenance marker, if present.
    #[must_use]
    pub fn geometry_authoring_provenance(&self) -> Option<&GeometryAuthoringProvenance> {
        self.metadata::<GeometryAuthoringProvenance>()
    }

    /// Whether this blueprint was created by the canonical geometry-authoring
    /// pipeline and is therefore eligible for canonical figure export.
    #[must_use]
    pub fn is_geometry_authored(&self) -> bool {
        self.geometry_authoring_provenance().is_some()
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

    /// Serialize the blueprint, including explicit render/layout state.
    ///
    /// # Errors
    /// Returns a serialization error if the blueprint cannot be encoded as JSON.
    pub fn to_json_pretty(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(self)
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
        self.nodes
            .iter()
            .filter(|n| matches!(n.kind, NodeKind::Inlet))
            .count()
    }

    /// Count outlet nodes.
    #[must_use]
    pub fn outlet_count(&self) -> usize {
        self.nodes
            .iter()
            .filter(|n| matches!(n.kind, NodeKind::Outlet))
            .count()
    }

    /// Count junction nodes.
    #[must_use]
    pub fn junction_count(&self) -> usize {
        self.nodes
            .iter()
            .filter(|n| matches!(n.kind, NodeKind::Junction))
            .count()
    }

    /// Count pipe channels.
    #[must_use]
    pub fn pipe_count(&self) -> usize {
        self.channels
            .iter()
            .filter(|c| matches!(c.kind, EdgeKind::Pipe))
            .count()
    }

    /// Total path length across all channels [m].
    #[must_use]
    pub fn total_length_m(&self) -> f64 {
        self.channels.iter().map(|c| c.length_m).sum()
    }

    // ── Therapy-zone helpers ─────────────────────────────────────────────────

    /// Total channel length [m] in channels tagged with the given therapy zone.
    ///
    /// Returns `0.0` if no channels carry [`TherapyZoneMetadata`] for that zone.
    /// Channels without metadata are not counted.
    #[must_use]
    pub fn length_in_zone(&self, zone: crate::domain::therapy_metadata::TherapyZone) -> f64 {
        use crate::domain::therapy_metadata::TherapyZoneMetadata;
        self.channels
            .iter()
            .filter(|c| {
                c.metadata
                    .as_ref()
                    .and_then(|m| m.get::<TherapyZoneMetadata>())
                    .is_some_and(|tz| tz.zone == zone)
            })
            .map(|c| c.length_m)
            .sum()
    }

    /// All channels tagged with [`VenturiGeometryMetadata`].
    ///
    /// Returns all channels in the blueprint that carry explicit venturi
    /// throat geometry parameters, as set by the venturi preset factories.
    #[must_use]
    pub fn venturi_channels(&self) -> Vec<&ChannelSpec> {
        self.channels
            .iter()
            .filter(|c| {
                c.venturi_geometry.is_some()
                    || c.metadata.as_ref().is_some_and(
                        crate::geometry::metadata::MetadataContainer::contains::<
                            crate::geometry::metadata::VenturiGeometryMetadata,
                        >,
                    )
            })
            .collect()
    }

    /// Serialize this blueprint to a pretty-printed JSON string.
    ///
    /// Alias for [`Self::to_json_pretty`] — the canonical single-source JSON
    /// path for roundtrip serialization.
    pub fn to_json(&self) -> Result<String, serde_json::Error> {
        self.to_json_pretty()
    }

    /// Deserialize a `NetworkBlueprint` from a JSON string.
    pub fn from_json(json: &str) -> Result<Self, serde_json::Error> {
        serde_json::from_str(json)
    }

    /// Detect crossing channel centerlines and insert junction nodes at each
    /// crossing, splitting the two channels into sub-segments.
    ///
    /// Returns a summary of how many intersections were found and the indices
    /// of the newly created junction nodes in [`Self::nodes`].
    pub fn resolve_channel_overlaps(&mut self) -> crate::geometry::IntersectionResult {
        crate::geometry::insert_intersection_nodes(self)
    }

    /// Count unresolved interior channel crossings without mutating the blueprint.
    #[must_use]
    pub fn unresolved_channel_overlap_count(&self) -> usize {
        crate::geometry::unresolved_intersection_count(self)
    }

    /// Report whether the blueprint still contains unresolved interior crossings.
    #[must_use]
    pub fn has_unresolved_channel_overlaps(&self) -> bool {
        self.unresolved_channel_overlap_count() > 0
    }

    /// Check that the blueprint is structurally consistent:
    /// - at least one node and one channel
    /// - every channel's `from`/`to` node ID exists in `self.nodes`
    ///
    /// # Errors
    ///
    /// Returns a descriptive `String` error on the first violation found.
    pub fn validate(&self) -> Result<(), String> {
        if self.nodes.is_empty() {
            return Err("NetworkBlueprint has no nodes".to_string());
        }
        if self.channels.is_empty() {
            return Err("NetworkBlueprint has no channels".to_string());
        }
        let mut node_ids = std::collections::HashSet::with_capacity(self.nodes.len());
        node_ids.extend(self.nodes.iter().map(|n| n.id.as_str()));
        for ch in &self.channels {
            if !node_ids.contains(ch.from.as_str()) {
                return Err(format!(
                    "Channel '{}' references unknown from-node '{}'",
                    ch.id.as_str(),
                    ch.from.as_str()
                ));
            }
            if !node_ids.contains(ch.to.as_str()) {
                return Err(format!(
                    "Channel '{}' references unknown to-node '{}'",
                    ch.id.as_str(),
                    ch.to.as_str()
                ));
            }
        }
        let overlap_count = self.unresolved_channel_overlap_count();
        if overlap_count > 0 {
            return Err(format!(
                "NetworkBlueprint '{}' contains {overlap_count} unresolved interior channel crossing(s)",
                self.name
            ));
        }
        if let Some(topology) = &self.topology {
            BlueprintTopologyFactory::validate_spec(topology)?;
        }
        Ok(())
    }

    /// Return a human-readable summary of this blueprint (replaces `ConversionSummary`).
    #[must_use]
    pub fn describe(&self) -> String {
        // Single pass over nodes and channels instead of four separate iterations.
        let mut inlets = 0usize;
        let mut outlets = 0usize;
        let mut junctions = 0usize;
        for n in &self.nodes {
            match n.kind {
                NodeKind::Inlet => inlets += 1,
                NodeKind::Outlet => outlets += 1,
                NodeKind::Junction => junctions += 1,
                NodeKind::Reservoir => {}
            }
        }
        let mut pipes = 0usize;
        let mut total_len = 0.0f64;
        for c in &self.channels {
            if matches!(c.kind, EdgeKind::Pipe) {
                pipes += 1;
            }
            total_len += c.length_m;
        }
        format!(
            "NetworkBlueprint '{}'\n\
             ─────────────────────────────────────\n\
             Nodes   : {} total  ({} inlets, {} outlets, {} junctions)\n\
             Channels: {}  ({} pipes)\n\
             Length  : {:.6} m total\n",
            self.name,
            self.nodes.len(),
            inlets,
            outlets,
            junctions,
            self.channels.len(),
            pipes,
            total_len,
        )
    }
}

impl fmt::Display for NetworkBlueprint {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.describe())
    }
}
