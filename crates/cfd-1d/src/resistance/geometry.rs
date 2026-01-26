//! Channel geometry definitions for resistance calculations.

use nalgebra::RealField;

/// Geometry enum for resistance model selection
#[derive(Debug, Clone)]
pub enum ChannelGeometry<T: RealField + Copy> {
    /// Circular channel
    Circular {
        /// Diameter of the circular channel
        diameter: T,
        /// Length of the channel
        length: T,
    },
    /// Rectangular channel
    Rectangular {
        /// Width of the rectangular channel
        width: T,
        /// Height of the rectangular channel
        height: T,
        /// Length of the channel
        length: T,
    },
    /// Elliptical channel
    Elliptical {
        /// Major axis length of the ellipse
        major_axis: T,
        /// Minor axis length of the ellipse
        minor_axis: T,
        /// Length of the channel
        length: T,
    },
    /// Trapezoidal channel
    Trapezoidal {
        /// Width at the top of the trapezoid
        top_width: T,
        /// Width at the bottom of the trapezoid
        bottom_width: T,
        /// Height of the trapezoid
        height: T,
        /// Length of the channel
        length: T,
    },
    /// Custom cross-section channel
    Custom {
        /// Cross-sectional area
        area: T,
        /// Hydraulic diameter (4 * area / perimeter)
        hydraulic_diameter: T,
        /// Length of the channel
        length: T,
    },
}
