//! Cross-sectional geometry definitions for channels

use aequitas::systems::si::quantities::{Area, Length};
use cfd_core::CfdScalar;
use serde::{Deserialize, Serialize};

/// Cross-sectional geometry
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum CrossSection<T: CfdScalar + Copy> {
    /// Rectangular cross-section
    Rectangular {
        /// Width of the rectangular channel
        width: Length<T>,
        /// Height of the rectangular channel
        height: Length<T>,
    },
    /// Circular cross-section
    Circular {
        /// Diameter of the circular channel
        diameter: Length<T>,
    },
    /// Elliptical cross-section
    Elliptical {
        /// Major axis length of the ellipse
        major_axis: Length<T>,
        /// Minor axis length of the ellipse
        minor_axis: Length<T>,
    },
    /// Trapezoidal cross-section
    Trapezoidal {
        /// Width at the top of the trapezoid
        top_width: Length<T>,
        /// Width at the bottom of the trapezoid
        bottom_width: Length<T>,
        /// Height of the trapezoid
        height: Length<T>,
    },
    /// Custom cross-section with area and hydraulic diameter
    Custom {
        /// Cross-sectional area
        area: Area<T>,
        /// Hydraulic diameter (4 * area / perimeter)
        hydraulic_diameter: Length<T>,
    },
}
