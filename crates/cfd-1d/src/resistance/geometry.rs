//! Channel geometry definitions for resistance calculations.

use cfd_core::error::{Error, Result};
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;

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

impl<T: RealField + Copy + FromPrimitive> ChannelGeometry<T> {
    /// Cross-sectional area [m²]
    pub fn cross_sectional_area(&self) -> Result<T> {
        let pi = T::from_f64(std::f64::consts::PI).unwrap_or_else(T::zero);
        let four = T::from_f64(4.0).unwrap_or_else(T::one);
        let two = T::from_f64(2.0).unwrap_or_else(T::one);
        match self {
            Self::Circular { diameter, .. } => Ok(pi * *diameter * *diameter / four),
            Self::Rectangular { width, height, .. } => Ok(*width * *height),
            Self::Elliptical { major_axis, minor_axis, .. } => {
                Ok(pi / four * *major_axis * *minor_axis)
            }
            Self::Trapezoidal { top_width, bottom_width, height, .. } => {
                Ok((*top_width + *bottom_width) * *height / two)
            }
            Self::Custom { area, .. } => Ok(*area),
        }
    }

    /// Hydraulic diameter Dh = 4A/P [m]
    pub fn hydraulic_diameter(&self) -> Result<T> {
        let two = T::from_f64(2.0).unwrap_or_else(T::one);
        let four = T::from_f64(4.0).unwrap_or_else(T::one);
        match self {
            Self::Circular { diameter, .. } => Ok(*diameter),
            Self::Rectangular { width, height, .. } => {
                let perimeter = two * (*width + *height);
                let area = *width * *height;
                Ok(four * area / perimeter)
            }
            Self::Elliptical { major_axis, minor_axis, .. } => {
                // Approximation: Dh ≈ 2ab/(a+b) for ellipse (Ramanujan-like)
                Ok(two * *major_axis * *minor_axis / (*major_axis + *minor_axis))
            }
            Self::Trapezoidal { top_width, bottom_width, height, .. } => {
                let area = self.cross_sectional_area()?;
                // Side lengths via Pythagorean theorem
                let half_diff = (*top_width - *bottom_width) / two;
                let side = (*height * *height + half_diff * half_diff).sqrt();
                let perimeter = *top_width + *bottom_width + two * side;
                Ok(four * area / perimeter)
            }
            Self::Custom { hydraulic_diameter, .. } => Ok(*hydraulic_diameter),
        }
    }

    /// Channel length [m]
    pub fn length(&self) -> T {
        match self {
            Self::Circular { length, .. }
            | Self::Rectangular { length, .. }
            | Self::Elliptical { length, .. }
            | Self::Trapezoidal { length, .. }
            | Self::Custom { length, .. } => *length,
        }
    }
}
