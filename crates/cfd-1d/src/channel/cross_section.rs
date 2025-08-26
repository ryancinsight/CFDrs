//! Cross-sectional geometry definitions for channels

use nalgebra::RealField;
use serde::{Deserialize, Serialize};
/// Cross-sectional geometry
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum CrossSection<T: RealField + Copy> {
    /// Rectangular cross-section
    Rectangular {
        /// Width of the rectangular channel
        width: T,
        /// Height of the rectangular channel
        height: T,
    },
    /// Circular cross-section
    Circular {
        /// Diameter of the circular channel
        diameter: T,
    /// Elliptical cross-section
    Elliptical {
        /// Major axis length of the ellipse
        major_axis: T,
        /// Minor axis length of the ellipse
        minor_axis: T,
    /// Trapezoidal cross-section
    Trapezoidal {
        /// Width at the top of the trapezoid
        top_width: T,
        /// Width at the bottom of the trapezoid
        bottom_width: T,
        /// Height of the trapezoid
    /// Custom cross-section with area and hydraulic diameter
    Custom {
        /// Cross-sectional area
        area: T,
        /// Hydraulic diameter (4 * area / perimeter)
        hydraulic_diameter: T,
}
