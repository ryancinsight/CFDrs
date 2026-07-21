//! Borrow-or-own scalar overlays with precomputed ranges.

use std::{borrow::Cow, collections::HashMap};

use iris::color::NamedColorMap;

use super::{color::colorize_normalized, field::AnalysisField, range::ScalarRange};
use crate::error::VisualizationResult;
use crate::visualizations::traits::Color;

type ScalarMap = HashMap<usize, f64>;

/// Complete scalar overlay for rendering CFD results on a schematic.
///
/// Node and edge maps use [`Cow`] so callers can lend existing solver output or
/// transfer ownership without parallel APIs. Range reduction occurs once when
/// each map enters the overlay; color lookup is then constant time.
#[derive(Debug, Clone)]
pub struct AnalysisOverlay<'a> {
    /// Physical field represented by this overlay.
    pub field: AnalysisField,
    node_data: Cow<'a, ScalarMap>,
    edge_data: Cow<'a, ScalarMap>,
    /// Iris color law used after scalar normalization.
    pub color_map: NamedColorMap,
    node_range: ScalarRange,
    edge_range: ScalarRange,
}

impl<'a> AnalysisOverlay<'a> {
    /// Create an empty overlay using an Iris color law.
    #[must_use]
    pub fn new(field: AnalysisField, color_map: NamedColorMap) -> Self {
        Self {
            field,
            node_data: Cow::Owned(HashMap::new()),
            edge_data: Cow::Owned(HashMap::new()),
            color_map,
            node_range: ScalarRange::default(),
            edge_range: ScalarRange::default(),
        }
    }

    /// Attach borrowed or owned node scalar values.
    ///
    /// # Errors
    ///
    /// Returns an invalid-parameters error when any scalar is non-finite.
    pub fn with_node_data(mut self, data: Cow<'a, ScalarMap>) -> VisualizationResult<Self> {
        self.node_range = ScalarRange::from_values("node data", data.values())?;
        self.node_data = data;
        Ok(self)
    }

    /// Attach borrowed or owned edge scalar values.
    ///
    /// # Errors
    ///
    /// Returns an invalid-parameters error when any scalar is non-finite.
    pub fn with_edge_data(mut self, data: Cow<'a, ScalarMap>) -> VisualizationResult<Self> {
        self.edge_range = ScalarRange::from_values("edge data", data.values())?;
        self.edge_data = data;
        Ok(self)
    }

    /// Borrow node scalar values without copying.
    #[must_use]
    pub fn node_data(&self) -> &ScalarMap {
        &self.node_data
    }

    /// Borrow edge scalar values without copying.
    #[must_use]
    pub fn edge_data(&self) -> &ScalarMap {
        &self.edge_data
    }

    /// Return the precomputed node scalar range.
    #[must_use]
    pub const fn node_range(&self) -> (f64, f64) {
        self.node_range.endpoints()
    }

    /// Return the precomputed edge scalar range.
    #[must_use]
    pub const fn edge_range(&self) -> (f64, f64) {
        self.edge_range.endpoints()
    }

    /// Return the mapped color for an edge channel, if present.
    #[must_use]
    pub fn edge_color(&self, channel_id: usize) -> Option<Color> {
        self.edge_data
            .get(&channel_id)
            .map(|value| colorize_normalized(self.edge_range.normalize(*value), self.color_map))
    }

    /// Return the mapped color for a node, if present.
    #[must_use]
    pub fn node_color(&self, node_id: usize) -> Option<Color> {
        self.node_data
            .get(&node_id)
            .map(|value| colorize_normalized(self.node_range.normalize(*value), self.color_map))
    }
}
