//! Cross-section profiles for channels.

use crate::domain::core::constants;
use crate::domain::core::scalar::Real;

/// Cross-section shape for a channel.
#[derive(Clone, Debug)]
pub enum ChannelProfile {
    /// Circular cross-section with a given radius.
    Circular {
        /// Channel radius.
        radius: Real,
        /// Number of segments around the circle.
        segments: usize,
    },
    /// Rectangular cross-section.
    Rectangular {
        /// Width of the channel.
        width: Real,
        /// Height of the channel.
        height: Real,
    },
    /// Rounded rectangular (stadium) cross-section.
    RoundedRectangular {
        /// Width of the channel.
        width: Real,
        /// Height of the channel.
        height: Real,
        /// Corner radius.
        corner_radius: Real,
        /// Segments per corner arc.
        corner_segments: usize,
    },
}

impl ChannelProfile {
    /// Generate 2D profile points in the local XY plane, centered at origin.
    ///
    /// Returns points in CCW winding order.
    pub fn generate_points(&self) -> Vec<[Real; 2]> {
        match self {
            ChannelProfile::Circular { radius, segments } => {
                let n = *segments;
                (0..n)
                    .map(|i| {
                        let angle = constants::TAU * (i as Real) / (n as Real);
                        [radius * angle.cos(), radius * angle.sin()]
                    })
                    .collect()
            }
            ChannelProfile::Rectangular { width, height } => {
                let hw = width / 2.0;
                let hh = height / 2.0;
                vec![[-hw, -hh], [hw, -hh], [hw, hh], [-hw, hh]]
            }
            ChannelProfile::RoundedRectangular {
                width,
                height,
                corner_radius,
                corner_segments,
            } => {
                let hw = width / 2.0 - corner_radius;
                let hh = height / 2.0 - corner_radius;
                let n = *corner_segments;
                let mut points = Vec::with_capacity(4 * (n + 1));

                // Four corners: bottom-right, top-right, top-left, bottom-left
                let corners = [
                    (hw, -hh, 3.0 * constants::FRAC_PI_2), // bottom-right → 270°
                    (hw, hh, 0.0),                         // top-right → 0°
                    (-hw, hh, constants::FRAC_PI_2),       // top-left → 90°
                    (-hw, -hh, constants::PI),             // bottom-left → 180°
                ];

                for (cx, cy, start_angle) in &corners {
                    for i in 0..=n {
                        let angle = start_angle + constants::FRAC_PI_2 * (i as Real) / (n as Real);
                        points.push([
                            cx + corner_radius * angle.cos(),
                            cy + corner_radius * angle.sin(),
                        ]);
                    }
                }

                points
            }
        }
    }

    /// Default circular profile for millifluidics.
    pub fn default_millifluidic() -> Self {
        Self::Circular {
            radius: constants::DEFAULT_CHANNEL_RADIUS,
            segments: 16,
        }
    }
}
