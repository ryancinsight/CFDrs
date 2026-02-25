//! geometry/types — Core 2D Geometry Types
//!
//! This module defines the fundamental data structures used throughout
//! the 2D microfluidic schematic design system. It provides types for
//! representing points, nodes, channels, and complete channel systems.
//!
//! The system supports extensible metadata through the metadata module,
//! allowing for easy addition of new tracking variables without breaking
//! existing functionality.

mod blueprint;
mod channel_system;
mod interchange;

pub use self::channel_system::ChannelSystem;
pub use self::interchange::{
    InterchangeChannel, InterchangeChannelProfile, InterchangeChannelSystem, InterchangeNode,
};

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
    /// A curved arc channel with a predefined path
    Arc {
        /// The sequence of points defining the arc path
        path: Vec<Point2D>,
    },
    /// A tapered frustum channel with variable width
    Frustum {
        /// The sequence of centerline points
        path: Vec<Point2D>,
        /// Width values at each centerline point
        widths: Vec<f64>,
        /// Width at the channel inlet
        #[serde(default = "default_frustum_inlet_width")]
        inlet_width: f64,
        /// Width at the throat (narrowest section)
        #[serde(default = "default_frustum_throat_width")]
        throat_width: f64,
        /// Width at the channel outlet
        #[serde(default = "default_frustum_outlet_width")]
        outlet_width: f64,
        /// Taper profile shape
        #[serde(default = "default_frustum_taper_profile")]
        taper_profile: TaperProfile,
        /// Normalized throat position along the channel (0.0 to 1.0)
        #[serde(default = "default_frustum_throat_position")]
        throat_position: f64,
        /// Whether this frustum acts as a venturi throat
        #[serde(default = "default_has_venturi_throat")]
        has_venturi_throat: bool,
    },
}

const fn default_frustum_taper_profile() -> TaperProfile {
    TaperProfile::Smooth
}

const fn default_frustum_throat_position() -> f64 {
    0.5
}

const fn default_has_venturi_throat() -> bool {
    true
}

fn default_frustum_inlet_width() -> f64 {
    1.0
}

fn default_frustum_throat_width() -> f64 {
    0.5
}

fn default_frustum_outlet_width() -> f64 {
    1.0
}

/// A microfluidic channel connecting two nodes
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

// ── Free functions ──────────────────────────────────────────────────────────

pub(crate) fn is_finite_point(point: Point2D) -> bool {
    point.0.is_finite() && point.1.is_finite()
}

pub(crate) fn validate_path(path: &[Point2D]) -> GeometryResult<()> {
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

/// Get the centerline points for a channel, accounting for its type.
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
        | ChannelType::Arc { path }
        | ChannelType::Frustum { path, .. } => path.clone(),
    }
}

pub(crate) fn polyline_length(points: &[Point2D]) -> f64 {
    points
        .windows(2)
        .map(|segment| {
            let dx = segment[1].0 - segment[0].0;
            let dy = segment[1].1 - segment[0].1;
            dx.hypot(dy)
        })
        .sum()
}

pub(crate) fn hydraulic_diameter(width: f64, height: f64) -> f64 {
    if width <= 0.0 || height <= 0.0 {
        0.0
    } else {
        2.0 * width * height / (width + height)
    }
}
