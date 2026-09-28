use super::{ChannelSpec, NodeSpec};
use crate::geometry::metadata::{BlueprintRenderHints, MetadataContainer};
use crate::topology::{BlueprintTopologySpec, TopologyLineageMetadata};
use serde::{Deserialize, Serialize};
use std::fmt;

mod analysis_impl;
mod metadata_impl;
mod venturi_impl;

/// A network blueprint graph of nodes and channels.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NetworkBlueprint {
    /// Blueprint name.
    pub name: String,
    /// Authored plate envelope in millimeters.
    pub box_dims: (f64, f64),
    /// Plate outline segments in schematic coordinates.
    pub box_outline: Vec<((f64, f64), (f64, f64))>,
    /// Graph nodes.
    pub nodes: Vec<NodeSpec>,
    /// Graph channels.
    pub channels: Vec<ChannelSpec>,
    /// Render hints for visualizers.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub render_hints: Option<BlueprintRenderHints>,
    /// Materialized topology specification, when geometry was authored.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub topology: Option<BlueprintTopologySpec>,
    /// Lineage metadata tracking the blueprint source.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lineage: Option<TopologyLineageMetadata>,
    /// Extension metadata container.
    #[serde(skip)]
    pub metadata: Option<MetadataContainer>,
    /// Whether the geometry was produced by the canonical authoring pipeline.
    #[serde(default)]
    pub geometry_authored: bool,
}

/// Summary of physical channel footprint overlap in a blueprint.
#[derive(Debug, Clone, Copy, Default)]
pub struct ChannelOverlapAnalysis {
    /// Maximum overlap fraction across all channel pairs [0, 1].
    pub max_overlap_fraction: f64,
    /// Width ratio (wider / narrower) at the most-overlapping channel pair.
    pub width_ratio_at_worst: f64,
    /// Number of channel pairs with any physical footprint overlap.
    pub overlap_pair_count: usize,
}

impl fmt::Display for NetworkBlueprint {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.describe())
    }
}

#[cfg(test)]
mod tests;
