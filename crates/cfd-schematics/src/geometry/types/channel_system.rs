//! `ChannelSystem` — the top-level 2D schematic layout container.

use crate::error::{GeometryError, GeometryResult};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use super::{
    is_finite_point, validate_path, Channel, ChannelType,
    ChannelTypeCategory, Node, Point2D,
};

/// Represents a complete microfluidic channel system
///
/// This is the main data structure that contains all the geometric information
/// needed to represent a 2D microfluidic schematic, including nodes, channels,
/// and the containing boundary box.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChannelSystem {
    /// Dimensions of the containing box (width, height)
    pub box_dims: (f64, f64),
    /// All nodes in the system
    pub nodes: Vec<Node>,
    /// All channels in the system
    pub channels: Vec<Channel>,
    /// Line segments defining the boundary box outline
    pub box_outline: Vec<(Point2D, Point2D)>,
}

impl ChannelSystem {
    pub(crate) const INTERCHANGE_SCHEMA_VERSION: &'static str = "1.0.0";
    pub(crate) const INTERCHANGE_LENGTH_UNITS: &'static str = "mm";
    pub(crate) const INTERCHANGE_PRODUCER_PREFIX: &'static str = "scheme";

    /// Export the channel system to JSON format
    ///
    /// This method serializes the entire channel system to a JSON string,
    /// making it easy to save, load, or transfer channel system data.
    ///
    /// # Returns
    ///
    /// A JSON string representation of the channel system, or an error if
    /// serialization fails.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use cfd_schematics::geometry::generator::create_geometry;
    /// use cfd_schematics::geometry::SplitType;
    /// use cfd_schematics::config::{GeometryConfig, ChannelTypeConfig};
    ///
    /// let system = create_geometry(
    ///     (200.0, 100.0),
    ///     &[SplitType::Bifurcation],
    ///     &GeometryConfig::default(),
    ///     &ChannelTypeConfig::AllStraight,
    /// );
    ///
    /// let json = system.to_json().expect("Failed to serialize");
    /// println!("Exported system: {}", json);
    /// ```
    pub fn to_json(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(self)
    }

    /// Import a channel system from JSON format
    ///
    /// This method deserializes a channel system from a JSON string.
    ///
    /// # Arguments
    ///
    /// * `json` - A JSON string representation of a channel system
    ///
    /// # Returns
    ///
    /// A channel system instance, or an error if deserialization fails.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use cfd_schematics::geometry::ChannelSystem;
    ///
    /// let json = r#"{"box_dims": [200.0, 100.0], "nodes": [], "channels": [], "box_outline": []}"#;
    /// let system = ChannelSystem::from_json(json).expect("Failed to deserialize");
    /// ```
    pub fn from_json(json: &str) -> Result<Self, serde_json::Error> {
        serde_json::from_str(json)
    }

    /// Validate structural integrity for downstream interchange and processing.
    ///
    /// This catches malformed imports before they are consumed by rendering or
    /// downstream 3D/CFD toolchains.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use cfd_schematics::config::{ChannelTypeConfig, GeometryConfig};
    /// use cfd_schematics::geometry::generator::create_geometry;
    /// use cfd_schematics::geometry::SplitType;
    ///
    /// let system = create_geometry(
    ///     (100.0, 50.0),
    ///     &[SplitType::Bifurcation],
    ///     &GeometryConfig::default(),
    ///     &ChannelTypeConfig::AllStraight,
    /// );
    ///
    /// system.validate().expect("system should be structurally valid");
    /// ```
    pub fn validate(&self) -> GeometryResult<()> {
        if self.box_dims.0 <= 0.0
            || self.box_dims.1 <= 0.0
            || !self.box_dims.0.is_finite()
            || !self.box_dims.1.is_finite()
        {
            return Err(GeometryError::invalid_box_dimensions(
                self.box_dims.0,
                self.box_dims.1,
            ));
        }

        for node in &self.nodes {
            if !is_finite_point(node.point) {
                return Err(GeometryError::invalid_point(node.point));
            }
        }

        for channel in &self.channels {
            if channel.from_node >= self.nodes.len() || channel.to_node >= self.nodes.len() {
                return Err(GeometryError::ChannelCreationFailed {
                    from_id: channel.from_node,
                    to_id: channel.to_node,
                    reason: "Channel references node index outside node list".to_string(),
                });
            }

            if channel.width <= 0.0
                || channel.height <= 0.0
                || !channel.width.is_finite()
                || !channel.height.is_finite()
            {
                return Err(GeometryError::ChannelCreationFailed {
                    from_id: channel.from_node,
                    to_id: channel.to_node,
                    reason: "Channel width and height must be finite positive values".to_string(),
                });
            }

            match &channel.channel_type {
                ChannelType::Straight => {}
                ChannelType::SmoothStraight { path }
                | ChannelType::Serpentine { path }
                | ChannelType::Arc { path } => {
                    validate_path(path)?;
                }
                ChannelType::Frustum {
                    path,
                    widths,
                    inlet_width,
                    throat_width,
                    outlet_width,
                    throat_position,
                    ..
                } => {
                    validate_path(path)?;

                    if *inlet_width <= 0.0
                        || *throat_width <= 0.0
                        || *outlet_width <= 0.0
                        || !inlet_width.is_finite()
                        || !throat_width.is_finite()
                        || !outlet_width.is_finite()
                    {
                        return Err(GeometryError::InvalidChannelPath {
                            reason: "Frustum inlet, throat, and outlet widths must be finite positive values"
                                .to_string(),
                        });
                    }

                    if widths.len() != path.len() {
                        return Err(GeometryError::InvalidChannelPath {
                            reason: "Frustum widths must match centerline point count".to_string(),
                        });
                    }

                    if widths.iter().any(|w| *w <= 0.0 || !w.is_finite()) {
                        return Err(GeometryError::InvalidChannelPath {
                            reason: "Frustum width profile must contain finite positive values"
                                .to_string(),
                        });
                    }

                    if *throat_width >= *inlet_width || *throat_width >= *outlet_width {
                        return Err(GeometryError::InvalidChannelPath {
                            reason: "Frustum throat width must be smaller than inlet and outlet"
                                .to_string(),
                        });
                    }

                    if *throat_position < 0.0 || *throat_position > 1.0 {
                        return Err(GeometryError::InvalidChannelPath {
                            reason: "Frustum throat position must be within [0.0, 1.0]".to_string(),
                        });
                    }
                }
            }
        }

        Ok(())
    }

    /// Get all line segments that make up this channel system
    ///
    /// This method extracts all the individual line segments from all channels
    /// in the system, which is useful for rendering and analysis.
    ///
    /// # Returns
    ///
    /// A vector of line segments, where each segment is represented as a tuple
    /// of two points (start, end).
    ///
    /// # Examples
    ///
    /// ```rust
    /// use cfd_schematics::geometry::generator::create_geometry;
    /// use cfd_schematics::geometry::SplitType;
    /// use cfd_schematics::config::{GeometryConfig, ChannelTypeConfig};
    ///
    /// let system = create_geometry(
    ///     (200.0, 100.0),
    ///     &[SplitType::Bifurcation],
    ///     &GeometryConfig::default(),
    ///     &ChannelTypeConfig::AllStraight,
    /// );
    ///
    /// let lines = system.get_lines();
    /// println!("System has {} line segments", lines.len());
    /// ```
    #[must_use]
    pub fn get_lines(&self) -> Vec<(Point2D, Point2D)> {
        let mut lines = self.box_outline.clone();
        for channel in &self.channels {
            match &channel.channel_type {
                ChannelType::Straight => {
                    if let (Some(from), Some(to)) = (
                        self.nodes.get(channel.from_node),
                        self.nodes.get(channel.to_node),
                    ) {
                        lines.push((from.point, to.point));
                    }
                }
                ChannelType::SmoothStraight { path }
                | ChannelType::Serpentine { path }
                | ChannelType::Arc { path }
                | ChannelType::Frustum { path, .. } => {
                    for segment in path.windows(2) {
                        lines.push((segment[0], segment[1]));
                    }
                }
            }
        }
        lines
    }

    /// Get all line segments with their associated channel types for colored rendering
    ///
    /// # Examples
    ///
    /// ```rust
    /// use cfd_schematics::geometry::generator::create_geometry;
    /// use cfd_schematics::geometry::SplitType;
    /// use cfd_schematics::config::{GeometryConfig, ChannelTypeConfig};
    ///
    /// let system = create_geometry(
    ///     (200.0, 100.0),
    ///     &[SplitType::Bifurcation],
    ///     &GeometryConfig::default(),
    ///     &ChannelTypeConfig::AllStraight,
    /// );
    ///
    /// let (boundary_lines, channel_lines) = system.get_lines_by_type();
    /// println!("System has {} boundary lines", boundary_lines.len());
    /// for (channel_type, lines) in channel_lines {
    ///     println!("Channel type {:?} has {} line segments", channel_type, lines.len());
    /// }
    /// ```
    #[must_use]
    #[allow(clippy::type_complexity)]
    pub fn get_lines_by_type(
        &self,
    ) -> (
        Vec<(Point2D, Point2D)>,
        HashMap<ChannelTypeCategory, Vec<(Point2D, Point2D)>>,
    ) {
        let boundary_lines = self.box_outline.clone();
        let mut channel_lines: HashMap<ChannelTypeCategory, Vec<(Point2D, Point2D)>> =
            HashMap::new();

        for channel in &self.channels {
            let category = ChannelTypeCategory::from(&channel.channel_type);
            let lines = channel_lines.entry(category).or_default();

            match &channel.channel_type {
                ChannelType::Straight => {
                    if let (Some(from), Some(to)) = (
                        self.nodes.get(channel.from_node),
                        self.nodes.get(channel.to_node),
                    ) {
                        lines.push((from.point, to.point));
                    }
                }
                ChannelType::SmoothStraight { path }
                | ChannelType::Serpentine { path }
                | ChannelType::Arc { path }
                | ChannelType::Frustum { path, .. } => {
                    for segment in path.windows(2) {
                        lines.push((segment[0], segment[1]));
                    }
                }
            }
        }

        (boundary_lines, channel_lines)
    }

    /// Get the path segments for all channels in the system
    ///
    /// # Examples
    ///
    /// ```rust
    /// use cfd_schematics::geometry::generator::create_geometry;
    /// use cfd_schematics::geometry::SplitType;
    /// use cfd_schematics::config::{GeometryConfig, ChannelTypeConfig, SerpentineConfig};
    ///
    /// let system = create_geometry(
    ///     (200.0, 100.0),
    ///     &[SplitType::Bifurcation],
    ///     &GeometryConfig::default(),
    ///     &ChannelTypeConfig::AllSerpentine(SerpentineConfig::default()),
    /// );
    ///
    /// let paths = system.get_path_segments();
    /// for (i, path) in paths.iter().enumerate() {
    ///     println!("Channel {} has {} points in its path", i, path.len());
    /// }
    /// ```
    #[must_use]
    pub fn get_path_segments(&self) -> Vec<Vec<Point2D>> {
        self.channels
            .iter()
            .filter_map(|c| match &c.channel_type {
                ChannelType::SmoothStraight { path }
                | ChannelType::Serpentine { path }
                | ChannelType::Arc { path }
                | ChannelType::Frustum { path, .. } => Some(path.clone()),
                ChannelType::Straight => None,
            })
            .collect()
    }
}
