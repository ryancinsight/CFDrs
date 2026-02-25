//! Interchange export types for downstream mesh/CFD toolchains.

use crate::config::TaperProfile;
use serde::{Deserialize, Serialize};

use super::{
    centerline_for_channel, hydraulic_diameter, polyline_length, ChannelType, Point2D,
};
use super::channel_system::ChannelSystem;

/// Interchange-ready node payload for downstream meshing/CFD tooling.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InterchangeNode {
    /// Stable node identifier.
    pub id: usize,
    /// Node position in millimeters.
    pub point_mm: Point2D,
}

/// Cross-section profile payload for a channel in interchange export.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "profile_type", rename_all = "snake_case")]
pub enum InterchangeChannelProfile {
    /// Circular cross-section.
    Circular {
        /// Diameter in millimeters.
        diameter_mm: f64,
        /// Cross-sectional area in square millimeters.
        cross_section_area_mm2: f64,
    },
    /// Rounded rectangular cross-section.
    RoundedRectangular {
        /// Width in millimeters.
        width_mm: f64,
        /// Height in millimeters.
        height_mm: f64,
        /// Corner radius in millimeters.
        corner_radius_mm: f64,
        /// Cross-sectional area in square millimeters.
        cross_section_area_mm2: f64,
        /// Hydraulic diameter in millimeters.
        hydraulic_diameter_mm: f64,
    },
    /// Constant-width rectangular channel profile.
    Constant {
        /// Width in millimeters.
        width_mm: f64,
        /// Height in millimeters.
        height_mm: f64,
        /// Cross-sectional area in square millimeters.
        cross_section_area_mm2: f64,
        /// Hydraulic diameter in millimeters.
        hydraulic_diameter_mm: f64,
    },
    /// Variable-width frustum profile with explicit venturi throat data.
    Frustum {
        /// Inlet width in millimeters.
        inlet_width_mm: f64,
        /// Minimum throat width in millimeters.
        throat_width_mm: f64,
        /// Outlet width in millimeters.
        outlet_width_mm: f64,
        /// Channel height in millimeters.
        height_mm: f64,
        /// Taper profile used to build the width transition.
        taper_profile: TaperProfile,
        /// Normalized throat position along channel centerline.
        throat_position: f64,
        /// Width profile sampled along the centerline in millimeters.
        width_profile_mm: Vec<f64>,
        /// Cross-sectional area profile sampled along centerline in square millimeters.
        area_profile_mm2: Vec<f64>,
        /// Explicit flag used by downstream 3D builders.
        has_venturi_throat: bool,
    },
}

/// Interchange-ready channel payload for downstream meshing/CFD tooling.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InterchangeChannel {
    /// Stable channel identifier.
    pub id: usize,
    /// Start node identifier.
    pub from_node_id: usize,
    /// End node identifier.
    pub to_node_id: usize,
    /// Centerline polyline in millimeters.
    pub centerline_mm: Vec<Point2D>,
    /// Total centerline length in millimeters.
    pub centerline_length_mm: f64,
    /// Channel cross-section profile.
    pub profile: InterchangeChannelProfile,
}

/// Interchange export payload designed for external mesh/CFD pipelines.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InterchangeChannelSystem {
    /// Schema version for compatibility checks in downstream tooling.
    pub schema_version: String,
    /// Producer identifier including crate version.
    pub producer: String,
    /// Length units used by all coordinates and dimensions.
    pub length_units: String,
    /// Bounding box dimensions in millimeters.
    pub box_dims_mm: (f64, f64),
    /// All nodes.
    pub nodes: Vec<InterchangeNode>,
    /// All channels with explicit centerlines and profiles.
    pub channels: Vec<InterchangeChannel>,
}

impl ChannelSystem {
    /// Build an explicit interchange payload suitable for mesh/CFD workflows.
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
    /// let interchange = system.to_interchange();
    /// assert_eq!(interchange.length_units, "mm");
    /// ```
    #[must_use]
    pub fn to_interchange(&self) -> InterchangeChannelSystem {
        let nodes = self
            .nodes
            .iter()
            .map(|node| InterchangeNode {
                id: node.id,
                point_mm: node.point,
            })
            .collect();

        let channels = self
            .channels
            .iter()
            .map(|channel| {
                let centerline = centerline_for_channel(channel, &self.nodes);
                let centerline_length = polyline_length(&centerline);

                let profile = if let ChannelType::Frustum {
                        widths,
                        inlet_width,
                        throat_width,
                        outlet_width,
                        taper_profile,
                        throat_position,
                        has_venturi_throat,
                        ..
                    } = &channel.channel_type {
                    InterchangeChannelProfile::Frustum {
                        inlet_width_mm: *inlet_width,
                        throat_width_mm: *throat_width,
                        outlet_width_mm: *outlet_width,
                        height_mm: channel.height,
                        taper_profile: *taper_profile,
                        throat_position: *throat_position,
                        width_profile_mm: widths.clone(),
                        area_profile_mm2: widths.iter().map(|w| w * channel.height).collect(),
                        has_venturi_throat: *has_venturi_throat,
                    }
                    } else {
                        let area = channel.width * channel.height;

                        if let Some(meta) = &channel.metadata {
                            if let Some(geom_meta) = meta.get::<crate::geometry::metadata::ChannelGeometryMetadata>() {
                                InterchangeChannelProfile::Circular {
                                    diameter_mm: geom_meta.channel_diameter_mm,
                                    cross_section_area_mm2: area,
                                }
                            } else {
                                InterchangeChannelProfile::Constant {
                                    width_mm: channel.width,
                                    height_mm: channel.height,
                                    cross_section_area_mm2: area,
                                    hydraulic_diameter_mm: hydraulic_diameter(
                                        channel.width,
                                        channel.height,
                                    ),
                                }
                            }
                        } else {
                            InterchangeChannelProfile::Constant {
                                width_mm: channel.width,
                                height_mm: channel.height,
                                cross_section_area_mm2: area,
                                hydraulic_diameter_mm: hydraulic_diameter(
                                    channel.width,
                                    channel.height,
                                ),
                            }
                        }
                    };

                InterchangeChannel {
                    id: channel.id,
                    from_node_id: channel.from_node,
                    to_node_id: channel.to_node,
                    centerline_mm: centerline,
                    centerline_length_mm: centerline_length,
                    profile,
                }
            })
            .collect();

        InterchangeChannelSystem {
            schema_version: Self::INTERCHANGE_SCHEMA_VERSION.to_string(),
            producer: format!(
                "{}/{}",
                Self::INTERCHANGE_PRODUCER_PREFIX,
                env!("CARGO_PKG_VERSION")
            ),
            length_units: Self::INTERCHANGE_LENGTH_UNITS.to_string(),
            box_dims_mm: self.box_dims,
            nodes,
            channels,
        }
    }

    /// Export the interchange payload as JSON.
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
    /// let json = system
    ///     .to_interchange_json()
    ///     .expect("interchange export should serialize");
    /// assert!(json.contains("\"schema_version\""));
    /// ```
    pub fn to_interchange_json(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(&self.to_interchange())
    }
}
