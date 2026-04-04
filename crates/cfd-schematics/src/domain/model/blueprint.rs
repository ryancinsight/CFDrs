use super::{ChannelSpec, EdgeKind, NodeKind, NodeSpec};
use crate::domain::therapy_metadata::TherapyZone;
use crate::geometry::metadata::{
    BlueprintRenderHints, ChannelVenturiSpec, GeometryAuthoringProvenance, MetadataContainer,
    VenturiGeometryMetadata,
};
use crate::topology::{
    BlueprintTopologyFactory, BlueprintTopologySpec, TopologyLineageMetadata,
    TopologyOptimizationStage, TreatmentActuationMode, VenturiConfig, VenturiPlacementSpec,
};
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
    /// Serializable flag recording that this blueprint was created through
    /// the canonical `create_geometry()` pipeline.  The in-memory
    /// [`MetadataContainer`] cannot survive JSON round-trips (`#[serde(skip)]`),
    /// so this persistent flag ensures deserialized blueprints still pass the
    /// `is_geometry_authored()` validation gate.
    #[serde(default)]
    pub geometry_authored: bool,
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
            geometry_authored: false,
        }
    }

    /// Create an empty `NetworkBlueprint` for callers that will author every
    /// node with an explicit geometric position via
    /// [`NodeSpec::new_at`](crate::domain::model::NodeSpec::new_at).
    ///
    /// Unlike [`Self::new`], this constructor is intended for manually authored
    /// schematics with explicit layout, so it does not carry the generic
    /// auto-layout deprecation warning.
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
    ///
    /// Returns `true` if either the in-memory metadata container carries a
    /// [`GeometryAuthoringProvenance`] marker **or** the serializable
    /// `geometry_authored` flag is set (survives JSON round-trips).
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

    /// Attach venturi augmentation metadata to one or more existing channels.
    ///
    /// When `config.target_channel_ids` is empty, the topology treatment lanes
    /// are used as the augmentation targets.
    ///
    /// # Errors
    ///
    /// Returns an error if the blueprint has no eligible channels or a target
    /// channel ID cannot be found.
    pub fn add_venturi(&mut self, config: &VenturiConfig) -> Result<(), String> {
        let target_channel_ids = if config.target_channel_ids.is_empty() {
            self.treatment_channel_ids()
        } else {
            config.target_channel_ids.clone()
        };
        if target_channel_ids.is_empty() {
            return Err("venturi augmentation requires at least one target channel".to_string());
        }

        let mut placements = Vec::with_capacity(target_channel_ids.len());
        for (index, channel_id) in target_channel_ids.into_iter().enumerate() {
            let channel = self
                .channels
                .iter_mut()
                .find(|channel| channel.id.as_str() == channel_id)
                .ok_or_else(|| {
                    format!(
                        "venturi augmentation target '{}' does not exist in blueprint '{}'",
                        channel_id, self.name
                    )
                })?;

            let resolved_inlet_width_m = if config.throat_geometry.inlet_width_m > 0.0 {
                config.throat_geometry.inlet_width_m
            } else {
                channel.effective_width_m()
            };
            let resolved_outlet_width_m = if config.throat_geometry.outlet_width_m > 0.0 {
                config.throat_geometry.outlet_width_m
            } else {
                channel.effective_width_m()
            };
            let venturi_geometry = VenturiGeometryMetadata {
                throat_width_m: config.throat_geometry.throat_width_m,
                throat_height_m: config.throat_geometry.throat_height_m,
                throat_length_m: config.throat_geometry.throat_length_m,
                inlet_width_m: resolved_inlet_width_m,
                outlet_width_m: resolved_outlet_width_m,
                convergent_half_angle_deg: config.throat_geometry.convergent_half_angle_deg,
                divergent_half_angle_deg: config.throat_geometry.divergent_half_angle_deg,
                throat_position: 0.5,
            };
            channel.venturi_geometry = Some(venturi_geometry.clone());
            if channel.metadata.is_none() {
                channel.metadata = Some(MetadataContainer::new());
            }
            let metadata = channel
                .metadata
                .as_mut()
                .expect("channel metadata container must exist");
            metadata.insert(venturi_geometry);
            metadata.insert(ChannelVenturiSpec {
                n_throats: config.serial_throat_count,
                is_ctc_stream: channel
                    .therapy_zone
                    .is_some_and(|zone| zone == TherapyZone::CancerTarget),
                throat_width_m: config.throat_geometry.throat_width_m,
                height_m: config.throat_geometry.throat_height_m,
                inter_throat_spacing_m: config.throat_geometry.throat_length_m,
            });

            placements.push(VenturiPlacementSpec {
                placement_id: format!("venturi_{index}"),
                target_channel_id: channel_id,
                serial_throat_count: config.serial_throat_count,
                throat_geometry: config.throat_geometry.clone(),
                placement_mode: config.placement_mode,
            });
        }

        if let Some(topology) = &mut self.topology {
            topology.venturi_placements = placements;
            topology.treatment_mode = TreatmentActuationMode::VenturiCavitation;
        }
        if let Some(lineage) = &mut self.lineage {
            lineage.current_stage =
                TopologyOptimizationStage::AsymmetricSplitVenturiCavitationSelectivity;
        }
        Ok(())
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

    /// Analyze physical footprint overlap between all channel pairs.
    ///
    /// Returns [`ChannelOverlapAnalysis`] with the worst-case overlap fraction
    /// and the width ratio at the most-overlapping pair.  When channels of
    /// different widths overlap, the velocity mismatch at the merge zone
    /// creates inertial cell-sorting effects (larger/stiffer CTCs entering
    /// from a narrow high-velocity channel experience different drag/lift
    /// than RBCs from a wide low-velocity bypass).  The width ratio quantifies
    /// this mismatch so downstream scoring can account for unmodeled
    /// separation physics.
    #[must_use]
    pub fn channel_overlap_analysis(&self) -> ChannelOverlapAnalysis {
        let channels_with_paths: Vec<(usize, f64)> = self
            .channels
            .iter()
            .enumerate()
            .filter(|(_, ch)| ch.path.len() >= 2)
            .map(|(i, ch)| {
                let w = match ch.cross_section {
                    super::CrossSectionSpec::Rectangular { width_m, .. } => width_m * 1e3,
                    super::CrossSectionSpec::Circular { diameter_m } => diameter_m * 1e3,
                };
                (i, w)
            })
            .collect();

        let mut worst = ChannelOverlapAnalysis::default();
        for (ai, (idx_a, w_a)) in channels_with_paths.iter().enumerate() {
            let path_a = &self.channels[*idx_a].path;
            for (idx_b, w_b) in &channels_with_paths[ai + 1..] {
                let path_b = &self.channels[*idx_b].path;
                let min_dist = min_path_distance(path_a, path_b);
                let half_sum = (w_a + w_b) * 0.5;
                if min_dist < half_sum {
                    let narrower = w_a.min(*w_b).max(1e-9);
                    let overlap_frac = ((half_sum - min_dist) / narrower).clamp(0.0, 1.0);
                    if overlap_frac > worst.max_overlap_fraction {
                        let wider = w_a.max(*w_b).max(1e-9);
                        worst.max_overlap_fraction = overlap_frac;
                        worst.width_ratio_at_worst = wider / narrower;
                        worst.overlap_pair_count += 1;
                    }
                }
            }
        }
        worst
    }

    /// Convenience: maximum overlap fraction (delegates to full analysis).
    #[must_use]
    pub fn max_channel_overlap_fraction(&self) -> f64 {
        self.channel_overlap_analysis().max_overlap_fraction
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
            if topology.is_selective_routing() && !self.is_geometry_authored() {
                return Err(format!(
                    "NetworkBlueprint '{}' carries selective split-tree topology '{}' but was not authored through create_geometry()",
                    self.name,
                    topology.stage_sequence_label()
                ));
            }
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::therapy_metadata::TherapyZone;
    use crate::geometry::metadata::ChannelVenturiSpec;
    use crate::topology::presets::parallel_path_spec;
    use crate::topology::{ChannelRouteSpec, ParallelChannelSpec, VenturiPlacementMode};

    #[test]
    fn add_venturi_attaches_metadata_to_existing_parallel_channel() {
        let topology = parallel_path_spec(
            "venturi-blueprint",
            2.0e-3,
            2.0e-3,
            12.0e-3,
            12.0e-3,
            vec![ParallelChannelSpec {
                channel_id: "treatment_lane".to_string(),
                route: ChannelRouteSpec {
                    length_m: 10.0e-3,
                    width_m: 1.6e-3,
                    height_m: 1.0e-3,
                    serpentine: None,
                    therapy_zone: TherapyZone::CancerTarget,
                },
            }],
            TreatmentActuationMode::UltrasoundOnly,
        );
        let mut blueprint =
            BlueprintTopologyFactory::build(&topology).expect("parallel topology should build");

        blueprint
            .add_venturi(&VenturiConfig {
                target_channel_ids: vec!["treatment_lane".to_string()],
                serial_throat_count: 2,
                throat_geometry: crate::topology::ThroatGeometrySpec {
                    throat_width_m: 80.0e-6,
                    throat_height_m: 1.0e-3,
                    throat_length_m: 300.0e-6,
                    inlet_width_m: 0.0,
                    outlet_width_m: 0.0,
                    convergent_half_angle_deg: 7.0,
                    divergent_half_angle_deg: 7.0,
                },
                placement_mode: VenturiPlacementMode::StraightSegment,
            })
            .expect("venturi should attach");

        let channel = blueprint
            .channels
            .iter()
            .find(|channel| channel.id.as_str() == "treatment_lane")
            .expect("target channel must exist");
        assert!(channel.venturi_geometry.is_some());
        assert_eq!(
            channel
                .metadata
                .as_ref()
                .and_then(|metadata| metadata.get::<ChannelVenturiSpec>())
                .expect("ChannelVenturiSpec must be inserted")
                .n_throats,
            2
        );
        assert_eq!(
            blueprint
                .topology
                .as_ref()
                .expect("topology metadata must be preserved")
                .venturi_placements
                .len(),
            1
        );
    }
}

/// Summary of physical channel footprint overlap in a blueprint.
#[derive(Debug, Clone, Copy, Default)]
pub struct ChannelOverlapAnalysis {
    /// Maximum overlap fraction across all channel pairs [0, 1].
    /// 0.0 = no overlap; 1.0 = one channel fully inside another.
    pub max_overlap_fraction: f64,

    /// Width ratio (wider / narrower) at the most-overlapping channel pair.
    /// A ratio of 1.0 means equal-width channels (symmetric merge).
    /// A ratio > 2.0 means a narrow channel merging into a channel more than
    /// twice its width, creating a velocity mismatch at the merge zone that
    /// produces inertial cell sorting: larger/stiffer CTCs from the narrow
    /// high-velocity channel experience different drag and lift forces than
    /// RBCs entering from the wide low-velocity bypass.
    pub width_ratio_at_worst: f64,

    /// Number of channel pairs with any physical footprint overlap.
    pub overlap_pair_count: usize,
}

/// Minimum Euclidean distance between two polyline paths, sampled at all
/// vertices of both paths (O(n*m) but paths are short in practice).
fn min_path_distance(path_a: &[(f64, f64)], path_b: &[(f64, f64)]) -> f64 {
    let mut min_d = f64::INFINITY;
    for &(ax, ay) in path_a {
        for &(bx, by) in path_b {
            let d = ((ax - bx).powi(2) + (ay - by).powi(2)).sqrt();
            if d < min_d {
                min_d = d;
            }
        }
    }
    min_d
}
