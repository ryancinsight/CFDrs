//! Interchange export types for downstream mesh/CFD toolchains.

use crate::config::TaperProfile;
use crate::domain::model::CrossSectionSpec;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use super::tpms_fill::TpmsFillSpec;
use super::{hydraulic_diameter, polyline_length, Point2D};
use crate::domain::model::NetworkBlueprint;

const INTERCHANGE_SCHEMA_VERSION: &str = "1.0.0";
const INTERCHANGE_LENGTH_UNITS: &str = "mm";
const INTERCHANGE_PRODUCER_PREFIX: &str = "scheme";

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

impl NetworkBlueprint {
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
        let node_indices: HashMap<&str, usize> = self
            .nodes
            .iter()
            .enumerate()
            .map(|(idx, node)| (node.id.as_str(), idx))
            .collect();

        let nodes = self
            .nodes
            .iter()
            .enumerate()
            .map(|(idx, node)| InterchangeNode {
                id: idx,
                point_mm: node.point,
            })
            .collect();

        let channels = self
            .channels
            .iter()
            .enumerate()
            .map(|(idx, channel)| {
                let centerline = if channel.path.is_empty() {
                    let from_point = self
                        .nodes
                        .iter()
                        .find(|node| node.id == channel.from)
                        .map_or((0.0, 0.0), |node| node.point);
                    let to_point = self
                        .nodes
                        .iter()
                        .find(|node| node.id == channel.to)
                        .map_or(from_point, |node| node.point);
                    vec![from_point, to_point]
                } else {
                    channel.path.clone()
                };
                let centerline_length = polyline_length(&centerline);

                let profile = match (channel.cross_section, channel.venturi_geometry.as_ref()) {
                    (CrossSectionSpec::Rectangular { height_m, .. }, Some(venturi)) => {
                        let inlet_width_mm = venturi.inlet_width_m * 1.0e3;
                        let throat_width_mm = venturi.throat_width_m * 1.0e3;
                        let outlet_width_mm = venturi.outlet_width_m * 1.0e3;
                        let height_mm = height_m * 1.0e3;
                        let width_profile_mm =
                            vec![inlet_width_mm, throat_width_mm, outlet_width_mm];
                        InterchangeChannelProfile::Frustum {
                            inlet_width_mm,
                            throat_width_mm,
                            outlet_width_mm,
                            height_mm,
                            taper_profile: TaperProfile::Linear,
                            throat_position: 0.5,
                            area_profile_mm2: width_profile_mm
                                .iter()
                                .map(|width_mm| width_mm * height_mm)
                                .collect(),
                            width_profile_mm,
                            has_venturi_throat: true,
                        }
                    }
                    (CrossSectionSpec::Circular { diameter_m }, _) => {
                        let diameter_mm = diameter_m * 1.0e3;
                        InterchangeChannelProfile::Circular {
                            diameter_mm,
                            cross_section_area_mm2: std::f64::consts::PI
                                * (diameter_mm * 0.5).powi(2),
                        }
                    }
                    (CrossSectionSpec::Rectangular { width_m, height_m }, _) => {
                        let width_mm = width_m * 1.0e3;
                        let height_mm = height_m * 1.0e3;
                        let area = width_mm * height_mm;
                        InterchangeChannelProfile::Constant {
                            width_mm,
                            height_mm,
                            cross_section_area_mm2: area,
                            hydraulic_diameter_mm: hydraulic_diameter(width_mm, height_mm),
                        }
                    }
                };

                InterchangeChannel {
                    id: idx,
                    from_node_id: *node_indices
                        .get(channel.from.as_str())
                        .expect("blueprint channel source node must exist"),
                    to_node_id: *node_indices
                        .get(channel.to.as_str())
                        .expect("blueprint channel target node must exist"),
                    centerline_mm: centerline,
                    centerline_length_mm: centerline_length,
                    profile,
                }
            })
            .collect();

        InterchangeChannelSystem {
            schema_version: INTERCHANGE_SCHEMA_VERSION.to_string(),
            producer: format!(
                "{}/{}",
                INTERCHANGE_PRODUCER_PREFIX,
                env!("CARGO_PKG_VERSION")
            ),
            length_units: INTERCHANGE_LENGTH_UNITS.to_string(),
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

    /// Build interchange payload after resolving geometric channel overlaps.
    ///
    /// This is useful when imported or procedurally modified schematics contain
    /// crossing centerlines that should become explicit junction nodes for
    /// downstream solver graph construction.
    pub fn to_interchange_resolved(&self) -> InterchangeChannelSystem {
        let mut resolved = self.clone();
        crate::geometry::insert_intersection_nodes(&mut resolved);
        resolved.to_interchange()
    }

    /// Export overlap-resolved interchange payload as JSON.
    pub fn to_interchange_json_resolved(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(&self.to_interchange_resolved())
    }
}

// ── Shell Cuboid interchange types ─────────────────────────────────────────────

/// Interchange payload describing a single port of a [`ShellCuboid`].
///
/// Each port is defined by two points: one on the **outer** wall surface and
/// one on the **inner** cavity surface, forming a short stub through the shell.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InterchangeShellPort {
    /// Human-readable port label (`"inlet"` | `"outlet"`).
    pub label: String,
    /// Midpoint on the outer wall surface (mm).
    pub outer_point_mm: Point2D,
    /// Midpoint on the inner cavity surface (mm).
    pub inner_point_mm: Point2D,
}

/// Interchange export payload for a [`ShellCuboid`].
///
/// Designed to give downstream mesh/CFD toolchains everything they need to
/// reconstruct the device geometry without depending on the schematic's visual
/// representation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InterchangeShellCuboid {
    /// Schema version for compatibility checks in downstream tooling.
    pub schema_version: String,
    /// Producer identifier including crate version.
    pub producer: String,
    /// Length units used by all coordinates and dimensions.
    pub length_units: String,
    /// Outer bounding-box dimensions `(width_mm, height_mm)`.
    pub outer_dims_mm: (f64, f64),
    /// Inner cavity dimensions `(width_mm, height_mm)`.
    pub inner_dims_mm: (f64, f64),
    /// Uniform wall thickness (mm).
    pub shell_thickness_mm: f64,
    /// Port stubs `[inlet, outlet]`.
    pub ports: Vec<InterchangeShellPort>,
    /// Optional TPMS lattice fill for the inner cavity.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tpms_fill: Option<TpmsFillSpec>,
}

impl super::shell_cuboid::ShellCuboid {
    /// Build an explicit interchange payload for mesh/CFD workflows.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use cfd_schematics::geometry::ShellCuboid;
    ///
    /// let sc = ShellCuboid::new((80.0, 40.0), 2.0).expect("structural invariant");
    /// let ix = sc.to_interchange();
    /// assert_eq!(ix.length_units, "mm");
    /// assert_eq!(ix.ports.len(), 2);
    /// ```
    #[must_use]
    pub fn to_interchange(&self) -> InterchangeShellCuboid {
        let (w, h) = self.outer_dims;
        let t = self.shell_thickness_mm;

        let inlet_port = InterchangeShellPort {
            label: "inlet".to_string(),
            outer_point_mm: (0.0, h / 2.0),
            inner_point_mm: (t, h / 2.0),
        };
        let outlet_port = InterchangeShellPort {
            label: "outlet".to_string(),
            outer_point_mm: (w, h / 2.0),
            inner_point_mm: (w - t, h / 2.0),
        };

        InterchangeShellCuboid {
            schema_version: Self::INTERCHANGE_SCHEMA_VERSION.to_string(),
            producer: format!(
                "{}/{}",
                Self::INTERCHANGE_PRODUCER_PREFIX,
                env!("CARGO_PKG_VERSION")
            ),
            length_units: Self::INTERCHANGE_LENGTH_UNITS.to_string(),
            outer_dims_mm: self.outer_dims,
            inner_dims_mm: self.inner_dims,
            shell_thickness_mm: self.shell_thickness_mm,
            ports: vec![inlet_port, outlet_port],
            tpms_fill: self.tpms_fill.clone(),
        }
    }

    /// Export the interchange payload as a pretty-printed JSON string.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use cfd_schematics::geometry::ShellCuboid;
    ///
    /// let sc = ShellCuboid::new((80.0, 40.0), 2.0).expect("structural invariant");
    /// let json = sc.to_interchange_json().expect("should serialize");
    /// assert!(json.contains("\"outer_dims_mm\""));
    /// assert!(json.contains("\"inlet\""));
    /// ```
    pub fn to_interchange_json(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(&self.to_interchange())
    }
}
