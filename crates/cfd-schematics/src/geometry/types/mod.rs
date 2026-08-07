//! geometry/types — Core 2D Geometry Types
//!
//! This module defines the fundamental data structures used throughout
//! the 2D microfluidic schematic design system. It provides types for
//! representing points, nodes, channels, and complete channel systems.
//!
//! The system supports extensible metadata through the metadata module,
//! allowing for easy addition of new tracking variables without breaking
//! existing functionality.

mod interchange;
mod shell_cuboid;
mod tpms_fill;
mod volume;

pub use self::interchange::{
    InterchangeChannel, InterchangeChannelProfile, InterchangeChannelSystem, InterchangeNode,
    InterchangeShellCuboid, InterchangeShellPort,
};
pub use self::shell_cuboid::ShellCuboid;
pub use self::tpms_fill::{AdaptiveGradient, TpmsFillSpec, TpmsSurfaceKind};
pub use self::volume::{ChannelFluidVolumeSummary, FluidVolumeSummary};

use crate::config::TaperProfile;
use serde::{Deserialize, Serialize};

/// A 2D point represented as (x, y) coordinates
pub type Point2D = (f64, f64);

/// Categories of channel types for visualization and analysis.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ChannelTypeCategory {
    /// Channels without curvature or taper.
    Straight,
    /// Channels with directional change.
    Curved,
    /// Channels with varying width along their length.
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

/// Represents the different types of channels that can be generated.
#[derive(Debug, Default, Clone, Serialize, Deserialize)]
pub enum ChannelType {
    /// Straight channel.
    #[default]
    Straight,
    /// Smoothly curved channel following the given path.
    SmoothStraight {
        /// Centerline waypoints.
        path: Vec<Point2D>,
    },
    /// Serpentine channel following the given path.
    Serpentine {
        /// Centerline waypoints.
        path: Vec<Point2D>,
    },
    /// Circular-arc channel following the given path.
    Arc {
        /// Centerline waypoints.
        path: Vec<Point2D>,
    },
    /// Tapered channel with authored inlet, throat, and outlet widths.
    Frustum {
        /// Centerline waypoints.
        path: Vec<Point2D>,
        /// Per-waypoint widths along the centerline.
        widths: Vec<f64>,
        /// Inlet width of the taper.
        #[serde(default = "default_frustum_inlet_width")]
        inlet_width: f64,
        /// Narrowest width of the taper.
        #[serde(default = "default_frustum_throat_width")]
        throat_width: f64,
        /// Outlet width of the taper.
        #[serde(default = "default_frustum_outlet_width")]
        outlet_width: f64,
        /// Profile governing the transition between widths.
        #[serde(default = "default_frustum_taper_profile")]
        taper_profile: TaperProfile,
        /// Normalized position of the throat along the centerline.
        #[serde(default = "default_frustum_throat_position")]
        throat_position: f64,
        /// Whether the taper represents a venturi throat.
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

/// Defines the type of channel splitting pattern.
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub enum SplitType {
    /// Two-way split.
    Bifurcation,
    /// Two-way split with an asymmetric width ratio.
    AsymmetricBifurcation {
        /// Ratio of the wider branch to the inlet width.
        ratio: f64,
    },
    /// Three-way split.
    Trifurcation,
    /// Three-way split with a symmetric center branch.
    SymmetricTrifurcation {
        /// Center-branch width fraction.
        center_ratio: f64,
    },
    /// Four-way split.
    Quadfurcation,
    /// Five-way split.
    Pentafurcation,
}

impl SplitType {
    /// Number of outlet branches for this split pattern.
    #[must_use]
    pub const fn branch_count(&self) -> usize {
        match self {
            Self::Bifurcation | Self::AsymmetricBifurcation { .. } => 2,
            Self::Trifurcation | Self::SymmetricTrifurcation { .. } => 3,
            Self::Quadfurcation => 4,
            Self::Pentafurcation => 5,
        }
    }
}

// ── Free functions ──────────────────────────────────────────────────────────

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
