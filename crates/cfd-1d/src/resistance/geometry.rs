//! Channel geometry definitions for resistance calculations.

use nalgebra::RealField;

/// Simplified geometry enum for resistance model selection
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
}
