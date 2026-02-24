//! geometry/types.rs - Core 2D Geometry Types
//!
//! This module defines the fundamental data structures used throughout
//! the 2D microfluidic schematic design system. It provides types for
//! representing points, nodes, channels, and complete channel systems.
//!
//! The system supports extensible metadata through the metadata module,
//! allowing for easy addition of new tracking variables without breaking
//! existing functionality.

use crate::config::TaperProfile;
use crate::error::{GeometryError, GeometryResult};
use crate::geometry::metadata::MetadataContainer;
use serde::{Deserialize, Serialize};

/// A 2D point represented as (x, y) coordinates
pub type Point2D = (f64, f64);

/// A node represents a connection point in the channel system
///
/// Nodes are used to define the endpoints of channels and serve as
/// junction points where multiple channels can connect.
///
/// The node supports extensible metadata for tracking additional properties
/// like pressure, temperature, or manufacturing tolerances.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Node {
    /// Unique identifier for this node
    pub id: usize,
    /// 2D coordinates of the node
    pub point: Point2D,
    /// Optional metadata container for extensible properties
    #[serde(skip)]
    pub metadata: Option<MetadataContainer>,
}

/// Categories of channel types for visualization and analysis
///
/// This enum groups channel types into categories for consistent coloring
/// and styling in visualizations.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ChannelTypeCategory {
    /// Straight line channels (Straight, `SmoothStraight`)
    Straight,
    /// Curved channels (Serpentine, Arc)
    Curved,
    /// Tapered channels (Frustum)
    Tapered,
}

impl From<&ChannelType> for ChannelTypeCategory {
    fn from(channel_type: &ChannelType) -> Self {
        match channel_type {
            ChannelType::Straight | ChannelType::SmoothStraight { .. } => Self::Straight,
            ChannelType::Serpentine { .. } | ChannelType::Arc { .. } => Self::Curved,
            ChannelType::Frustum { .. } => Self::Tapered,
        }
    }
}

/// Represents the different types of channels that can be generated
///
/// Each channel type has different characteristics:
/// - `Straight`: Direct line between two points
/// - `SmoothStraight`: Straight line with optional smooth transition zones at endpoints
/// - `Serpentine`: Sinusoidal path with Gaussian envelope for smooth transitions
/// - `Arc`: Curved path using quadratic Bezier curves
/// - `Frustum`: Tapered channel with variable width for venturi throat functionality
#[derive(Debug, Default, Clone, Serialize, Deserialize)]
pub enum ChannelType {
    /// A straight line channel between two points
    #[default]
    Straight,
    /// A straight line channel with smooth transition zones at endpoints
    SmoothStraight {
        /// The sequence of points defining the smooth straight path with transitions
        path: Vec<Point2D>,
    },
    /// A serpentine (S-shaped) channel with a predefined path
    Serpentine {
        /// The sequence of points defining the serpentine path
        path: Vec<Point2D>,
    },
    /// An arc channel with a curved path
    Arc {
        /// The sequence of points defining the arc path
        path: Vec<Point2D>,
    },
    /// A frustum (tapered) channel with variable width for venturi throat functionality
    Frustum {
        /// The sequence of points defining the centerline path
        path: Vec<Point2D>,
        /// Width values corresponding to each point in the path
        widths: Vec<f64>,
        /// Inlet width (starting width)
        inlet_width: f64,
        /// Throat width (minimum width at center)
        throat_width: f64,
        /// Outlet width (ending width)
        outlet_width: f64,
        /// Taper profile used to generate the frustum width transition
        #[serde(default = "default_frustum_taper_profile")]
        taper_profile: TaperProfile,
        /// Throat position along channel length (0.0 to 1.0)
        #[serde(default = "default_frustum_throat_position")]
        throat_position: f64,
        /// Explicit tag for downstream 3D reconstruction workflows
        #[serde(default = "default_has_venturi_throat")]
        has_venturi_throat: bool,
    },
}

const fn default_frustum_taper_profile() -> TaperProfile {
    TaperProfile::Linear
}

const fn default_frustum_throat_position() -> f64 {
    0.5
}

const fn default_has_venturi_throat() -> bool {
    true
}

/// Represents a single channel in the microfluidic system
///
/// A channel connects two nodes and has physical properties like width and height.
/// The channel type determines how the path between the nodes is generated.
///
/// The channel supports extensible metadata for tracking additional properties
/// like flow rates, pressure drops, optimization history, or manufacturing data.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Channel {
    /// Unique identifier for this channel
    pub id: usize,
    /// ID of the starting node
    pub from_node: usize,
    /// ID of the ending node
    pub to_node: usize,
    /// Physical width of the channel
    pub width: f64,
    /// Physical height of the channel
    pub height: f64,
    /// The type and path of this channel
    pub channel_type: ChannelType,
    /// Optional metadata container for extensible properties
    #[serde(skip)]
    pub metadata: Option<MetadataContainer>,
}

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
    const INTERCHANGE_SCHEMA_VERSION: &'static str = "1.0.0";
    const INTERCHANGE_LENGTH_UNITS: &'static str = "mm";
    const INTERCHANGE_PRODUCER_PREFIX: &'static str = "scheme";

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
    /// use scheme::geometry::generator::create_geometry;
    /// use scheme::geometry::SplitType;
    /// use scheme::config::{GeometryConfig, ChannelTypeConfig};
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
    /// use scheme::geometry::ChannelSystem;
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
    /// use scheme::config::{ChannelTypeConfig, GeometryConfig};
    /// use scheme::geometry::generator::create_geometry;
    /// use scheme::geometry::SplitType;
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

    /// Build an explicit interchange payload suitable for mesh/CFD workflows.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use scheme::config::{ChannelTypeConfig, GeometryConfig};
    /// use scheme::geometry::generator::create_geometry;
    /// use scheme::geometry::SplitType;
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

                let profile = match &channel.channel_type {
                    ChannelType::Frustum {
                        widths,
                        inlet_width,
                        throat_width,
                        outlet_width,
                        taper_profile,
                        throat_position,
                        has_venturi_throat,
                        ..
                    } => InterchangeChannelProfile::Frustum {
                        inlet_width_mm: *inlet_width,
                        throat_width_mm: *throat_width,
                        outlet_width_mm: *outlet_width,
                        height_mm: channel.height,
                        taper_profile: *taper_profile,
                        throat_position: *throat_position,
                        width_profile_mm: widths.clone(),
                        area_profile_mm2: widths.iter().map(|w| w * channel.height).collect(),
                        has_venturi_throat: *has_venturi_throat,
                    },
                    _ => {
                        let area = channel.width * channel.height;
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
    /// This format is intended for downstream meshing and CFD tools that
    /// require explicit centerlines and cross-section profile data.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use scheme::config::{ChannelTypeConfig, GeometryConfig};
    /// use scheme::geometry::generator::create_geometry;
    /// use scheme::geometry::SplitType;
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

    /// Convert this 2D schematic layout into a [`NetworkBlueprint`] for the 1D solver.
    ///
    /// `scale_m_per_unit` converts schematic coordinate units to SI metres
    /// (e.g. `1e-3` if coordinates are in millimetres, `1e-6` for micrometres).
    ///
    /// # Node classification
    /// - **Inlet** : no incoming channels (source node)
    /// - **Outlet**: no outgoing channels (sink node)
    /// - **Junction**: both incoming and outgoing channels
    ///
    /// # Resistance initialisation
    /// Each `ChannelSpec` is populated with an analytical Hagen-Poiseuille
    /// resistance as a linearisation seed (μ = 0.001 Pa·s, rectangular thin-slit
    /// approximation).  The 1D solver will refine these coefficients with the
    /// actual fluid model on its first Newton step.
    ///
    /// # Errors
    /// Returns [`GeometryError`] when the system has no channels, contains an
    /// out-of-bounds node reference, or produces a network without at least one
    /// inlet and one outlet.
    pub fn to_blueprint(
        &self,
        scale_m_per_unit: f64,
    ) -> GeometryResult<crate::domain::model::NetworkBlueprint> {
        use crate::domain::model::{
            CrossSectionSpec, EdgeKind, NetworkBlueprint, NodeKind, NodeSpec, ChannelSpec,
        };
        use std::collections::HashSet;

        if self.channels.is_empty() {
            return Err(GeometryError::InvalidChannelPath {
                reason: "Cannot build NetworkBlueprint from an empty ChannelSystem".to_string(),
            });
        }

        // ── topology classification ───────────────────────────────────────────
        let max_idx = self.nodes.len();
        for ch in &self.channels {
            if ch.from_node >= max_idx || ch.to_node >= max_idx {
                return Err(GeometryError::ChannelCreationFailed {
                    from_id: ch.from_node,
                    to_id: ch.to_node,
                    reason: format!(
                        "Channel {} references a node index outside the node list (len={})",
                        ch.id, max_idx
                    ),
                });
            }
        }

        let mut has_incoming: HashSet<usize> = HashSet::new();
        let mut has_outgoing: HashSet<usize> = HashSet::new();
        for ch in &self.channels {
            has_outgoing.insert(ch.from_node);
            has_incoming.insert(ch.to_node);
        }

        let classify = |id: usize| -> NodeKind {
            match (has_incoming.contains(&id), has_outgoing.contains(&id)) {
                (false, _) => NodeKind::Inlet,
                (true, false) => NodeKind::Outlet,
                (true, true) => NodeKind::Junction,
            }
        };

        // ── build blueprint ───────────────────────────────────────────────────
        let mut blueprint = NetworkBlueprint::new(format!(
            "from_channel_system_{}x{}",
            self.box_dims.0 as u64, self.box_dims.1 as u64
        ));

        // Validate inlet/outlet presence
        let has_inlet = self.nodes.iter().any(|n| matches!(classify(n.id), NodeKind::Inlet));
        let has_outlet = self.nodes.iter().any(|n| matches!(classify(n.id), NodeKind::Outlet));
        if !has_inlet {
            return Err(GeometryError::InvalidChannelPath {
                reason: "ChannelSystem has no inlet nodes (every node has incoming channels)"
                    .to_string(),
            });
        }
        if !has_outlet {
            return Err(GeometryError::InvalidChannelPath {
                reason: "ChannelSystem has no outlet nodes (every node has outgoing channels)"
                    .to_string(),
            });
        }

        for node in &self.nodes {
            blueprint.add_node(NodeSpec::new(
                format!("node_{}", node.id),
                classify(node.id),
            ));
        }

        // ── channels → ChannelSpec ────────────────────────────────────────────
        for channel in &self.channels {
            let centerline = centerline_for_channel(channel, &self.nodes);
            let length_m = polyline_length(&centerline) * scale_m_per_unit;

            if length_m <= 0.0 {
                return Err(GeometryError::ChannelCreationFailed {
                    from_id: channel.from_node,
                    to_id: channel.to_node,
                    reason: format!("Channel {} has zero or negative path length", channel.id),
                });
            }

            // Average cross-section dims for frustum channels; uniform otherwise.
            let (avg_width, height) = match &channel.channel_type {
                ChannelType::Frustum {
                    inlet_width,
                    throat_width,
                    outlet_width,
                    ..
                } => ((inlet_width + throat_width + outlet_width) / 3.0, channel.height),
                _ => (channel.width, channel.height),
            };
            let width_m = avg_width * scale_m_per_unit;
            let height_m = height * scale_m_per_unit;

            // Hagen-Poiseuille thin-slit approximation: R ≈ 12μL / (wh³)
            // Seeded with water viscosity at 20 °C; refined by solver.
            let mu = 0.001_f64; // Pa·s  (water at 20 °C)
            let resistance = (12.0 * mu * length_m / (width_m * height_m.powi(3))).max(f64::EPSILON);

            let cross_section = CrossSectionSpec::Rectangular { width_m, height_m };

            blueprint.add_channel(ChannelSpec {
                id: crate::domain::model::EdgeId::new(format!("ch_{}", channel.id)),
                kind: EdgeKind::Pipe,
                from: crate::domain::model::NodeId::new(format!("node_{}", channel.from_node)),
                to: crate::domain::model::NodeId::new(format!("node_{}", channel.to_node)),
                length_m,
                cross_section,
                resistance,
                quad_coeff: 0.0,
                valve_cv: None,
                pump_max_flow: None,
                pump_max_pressure: None,
                metadata: None,
            });
        }

        Ok(blueprint)
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
    /// use scheme::geometry::generator::create_geometry;
    /// use scheme::geometry::SplitType;
    /// use scheme::config::{GeometryConfig, ChannelTypeConfig};
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
                | ChannelType::Arc { path } => {
                    for segment in path.windows(2) {
                        lines.push((segment[0], segment[1]));
                    }
                }
                ChannelType::Frustum { path, .. } => {
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
    /// This method extracts line segments along with channel type information,
    /// enabling different colors for different channel types in visualization.
    ///
    /// # Returns
    ///
    /// A tuple containing:
    /// - Boundary lines (for the box outline)
    /// - Channel lines grouped by channel type
    ///
    /// # Examples
    ///
    /// ```rust
    /// use scheme::geometry::generator::create_geometry;
    /// use scheme::geometry::SplitType;
    /// use scheme::config::{GeometryConfig, ChannelTypeConfig};
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
        std::collections::HashMap<ChannelTypeCategory, Vec<(Point2D, Point2D)>>,
    ) {
        use std::collections::HashMap;

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
                | ChannelType::Arc { path } => {
                    for segment in path.windows(2) {
                        lines.push((segment[0], segment[1]));
                    }
                }
                ChannelType::Frustum { path, .. } => {
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
    /// This method returns the complete path information for each channel,
    /// which is particularly useful for serpentine and arc channels that
    /// have complex paths with multiple points.
    ///
    /// # Returns
    ///
    /// A vector where each element is a vector of points representing
    /// the complete path of one channel.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use scheme::geometry::generator::create_geometry;
    /// use scheme::geometry::SplitType;
    /// use scheme::config::{GeometryConfig, ChannelTypeConfig, SerpentineConfig};
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
                | ChannelType::Arc { path } => Some(path.clone()),
                ChannelType::Frustum { path, .. } => Some(path.clone()),
                _ => None,
            })
            .collect()
    }
}

fn is_finite_point(point: Point2D) -> bool {
    point.0.is_finite() && point.1.is_finite()
}

fn validate_path(path: &[Point2D]) -> GeometryResult<()> {
    if path.len() < 2 {
        return Err(GeometryError::InvalidChannelPath {
            reason: "Path must contain at least 2 points".to_string(),
        });
    }

    if path.iter().any(|point| !is_finite_point(*point)) {
        return Err(GeometryError::InvalidChannelPath {
            reason: "Path contains non-finite coordinates".to_string(),
        });
    }

    Ok(())
}

pub fn centerline_for_channel(channel: &Channel, nodes: &[Node]) -> Vec<Point2D> {
    match &channel.channel_type {
        ChannelType::Straight => {
            if let (Some(from), Some(to)) =
                (nodes.get(channel.from_node), nodes.get(channel.to_node))
            {
                vec![from.point, to.point]
            } else {
                Vec::new()
            }
        }
        ChannelType::SmoothStraight { path }
        | ChannelType::Serpentine { path }
        | ChannelType::Arc { path } => path.clone(),
        ChannelType::Frustum { path, .. } => path.clone(),
    }
}

fn polyline_length(points: &[Point2D]) -> f64 {
    points
        .windows(2)
        .map(|segment| {
            let dx = segment[1].0 - segment[0].0;
            let dy = segment[1].1 - segment[0].1;
            dx.hypot(dy)
        })
        .sum()
}

fn hydraulic_diameter(width: f64, height: f64) -> f64 {
    if width <= 0.0 || height <= 0.0 {
        0.0
    } else {
        2.0 * width * height / (width + height)
    }
}

/// Defines the type of channel splitting pattern
///
/// Split types determine how many branches are created at each junction:
/// - `Bifurcation`: Splits into 2 branches
/// - `Trifurcation`: Splits into 3 branches
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub enum SplitType {
    /// Split into two branches
    Bifurcation,
    /// Split into two branches with asymmetric width distribution.
    /// The ratio (0.0-1.0) determines the fraction of available width allocated to the first (top) branch.
    AsymmetricBifurcation {
        /// Ratio of total width allocated to the first branch
        ratio: f64,
    },
    /// Split into three branches
    Trifurcation,
    /// Split into three branches with symmetric width distribution (Center vs Sides).
    /// The center_ratio (0.0-1.0) determines the fraction of total width allocated to the center branch.
    /// Side branches share the remaining width equally.
    SymmetricTrifurcation {
        /// Ratio of total width allocated to the center branch
        center_ratio: f64,
    },
}

impl SplitType {
    /// Returns the number of branches created by this split type
    #[must_use]
    pub const fn branch_count(&self) -> usize {
        match self {
            Self::Bifurcation | Self::AsymmetricBifurcation { .. } => 2,
            Self::Trifurcation | Self::SymmetricTrifurcation { .. } => 3,
        }
    }
}

// CFD functionality removed - Scheme focuses exclusively on 2D schematic design
