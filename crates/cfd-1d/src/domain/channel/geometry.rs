//! Channel geometry definitions for 1D CFD

use super::cross_section::CrossSection;
use super::surface::SurfaceProperties;
use aequitas::systems::si::quantities::Length;
use cfd_core::CfdScalar;
use serde::{Deserialize, Serialize};

/// Extended channel geometry representation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChannelGeometry<T: CfdScalar + Copy> {
    /// Channel type
    pub channel_type: ChannelType,
    /// Length \[m]
    pub length: Length<T>,
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
        radius: Length<f64>,
    },
    /// Tapered channel
    Tapered,
    /// Serpentine channel
    Serpentine {
        /// Number of turns in the serpentine path
        turns: usize,
    },
    /// Spiral channel
    Spiral {
        /// Number of turns in the spiral (can be fractional)
        turns: f64,
    },
}

/// Geometric variation along channel length
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeometricVariation<T: CfdScalar + Copy> {
    /// Position along channel [0-1]
    pub position: T,
    /// Scale factor for cross-section
    pub scale_factor: T,
    /// Local roughness modification
    pub roughness_factor: T,
}

#[cfg(test)]
mod tests {
    use super::ChannelType;
    use aequitas::systems::si::quantities::Length;

    #[test]
    fn curved_channel_radius_preserves_si_length() {
        let channel = ChannelType::Curved {
            radius: Length::from_base(5.0e-3),
        };

        let ChannelType::Curved { radius } = channel else {
            panic!("invariant: channel is curved");
        };
        assert_eq!(radius.into_base(), 5.0e-3);
    }
}
