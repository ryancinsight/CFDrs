//! Channel geometry definitions for 1D CFD

use super::cross_section::CrossSection;
use super::surface::SurfaceProperties;
use nalgebra::RealField;
use serde::{Deserialize, Serialize};
/// Extended channel geometry representation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChannelGeometry<T: RealField + Copy> {
    /// Channel type
    pub channel_type: ChannelType,
    /// Length [m]
    pub length: T,
    /// Cross-sectional parameters
    pub cross_section: CrossSection<T>,
    /// Surface properties
    pub surface: SurfaceProperties<T>,
    /// Geometric variations along length
    pub variations: Vec<GeometricVariation<T>>,
}
/// Types of channel geometries
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum ChannelType {
    /// Straight channel
    Straight,
    /// Curved channel
    Curved {
        /// Radius of curvature in meters
        radius: f64,
    },
    /// Tapered channel
    Tapered,
    /// Serpentine channel
    Serpentine {
        /// Number of turns in the serpentine path
        turns: usize,
    /// Spiral channel
    Spiral {
        /// Number of turns in the spiral (can be fractional)
        turns: f64,
/// Geometric variation along channel length
pub struct GeometricVariation<T: RealField + Copy> {
    /// Position along channel [0-1]
    pub position: T,
    /// Scale factor for cross-section
    pub scale_factor: T,
    /// Local roughness modification
    pub roughness_factor: T,


}
}
}
}
