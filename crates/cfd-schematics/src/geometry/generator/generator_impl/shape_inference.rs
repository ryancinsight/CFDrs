//! Reading a serpentine channel's shape back out of its own path.
//!
//! The builder records the path it generated; these read curvature out of
//! that path to classify the result. That is analysis over a finished
//! geometry, not a step in constructing one, so it sits beside the builder
//! rather than inside it.

use super::super::super::types::Point2D;
use super::GeometryGenerator;
use crate::domain::model::ChannelShape;
use aequitas::systems::si::quantities::Length;

impl GeometryGenerator {
    pub(super) fn estimate_local_bend_radius(path: &[Point2D], idx: usize) -> Option<f64> {
        if idx == 0 || idx + 1 >= path.len() {
            return None;
        }

        let a = path[idx - 1];
        let b = path[idx];
        let c = path[idx + 1];
        let ab = ((b.0 - a.0).powi(2) + (b.1 - a.1).powi(2)).sqrt();
        let bc = ((c.0 - b.0).powi(2) + (c.1 - b.1).powi(2)).sqrt();
        let ac = ((c.0 - a.0).powi(2) + (c.1 - a.1).powi(2)).sqrt();
        if ab <= 1e-9 || bc <= 1e-9 || ac <= 1e-9 {
            return None;
        }

        let twice_area = ((b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0)).abs();
        if twice_area <= 1e-12 {
            return Some(f64::INFINITY);
        }

        Some((ab * bc * ac) / (2.0 * twice_area))
    }

    pub(super) fn infer_serpentine_shape(
        path: &[Point2D],
        p1: Point2D,
        p2: Point2D,
        channel_width: f64,
    ) -> ChannelShape {
        if path.len() < 3 {
            return ChannelShape::Serpentine {
                segments: 2,
                bend_radius_m: Length::from_base((channel_width * 0.5) * 1.0e-3),
                wave_type: crate::topology::SerpentineWaveType::Sine,
            };
        }

        let dx = p2.0 - p1.0;
        let dy = p2.1 - p1.1;
        let length = dx.hypot(dy);
        if length <= 1e-9 {
            return ChannelShape::Serpentine {
                segments: 2,
                bend_radius_m: Length::from_base((channel_width * 0.5) * 1.0e-3),
                wave_type: crate::topology::SerpentineWaveType::Sine,
            };
        }

        let nx = -dy / length;
        let ny = dx / length;
        let offsets: Vec<f64> = path
            .iter()
            .map(|point| ((point.0 - p1.0) * nx) + ((point.1 - p1.1) * ny))
            .collect();
        let extrema_threshold = channel_width * 0.15;
        let mut turns = 0usize;
        for idx in 1..offsets.len() - 1 {
            let prev_delta = offsets[idx] - offsets[idx - 1];
            let next_delta = offsets[idx + 1] - offsets[idx];
            if prev_delta.abs() <= 1e-9 || next_delta.abs() <= 1e-9 {
                continue;
            }
            if prev_delta.signum() != next_delta.signum() && offsets[idx].abs() > extrema_threshold
            {
                turns += 1;
            }
        }

        let min_radius_mm = path
            .iter()
            .enumerate()
            .filter_map(|(idx, _)| Self::estimate_local_bend_radius(path, idx))
            .filter(|radius| radius.is_finite())
            .fold(f64::INFINITY, f64::min);
        let bend_radius_mm = if min_radius_mm.is_finite() {
            min_radius_mm.max(channel_width * 0.5)
        } else {
            channel_width * 0.5
        };

        ChannelShape::Serpentine {
            segments: turns.saturating_add(1).max(2),
            bend_radius_m: Length::from_base(bend_radius_mm * 1.0e-3),
            wave_type: crate::topology::SerpentineWaveType::Sine,
        }
    }
}
