//! Shared path-geometry inference for geometry generators.

use super::super::types::Point2D;
use crate::domain::model::ChannelShape;
use aequitas::systems::si::quantities::Length;

/// Estimate the local circumradius at a path vertex.
fn estimate_local_bend_radius(path: &[Point2D], idx: usize) -> Option<f64> {
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

fn default_serpentine_shape(channel_width_mm: f64) -> ChannelShape {
    ChannelShape::Serpentine {
        segments: 2,
        bend_radius_m: Length::from_base((channel_width_mm * 0.5) * 1.0e-3),
        wave_type: crate::topology::SerpentineWaveType::Sine,
    }
}

/// Infer a serpentine shape from a path without materializing offset values.
pub(super) fn infer_serpentine_shape(
    path: &[Point2D],
    start: Point2D,
    end: Point2D,
    channel_width_mm: f64,
) -> ChannelShape {
    if path.len() < 3 {
        return default_serpentine_shape(channel_width_mm);
    }

    let dx = end.0 - start.0;
    let dy = end.1 - start.1;
    let length = dx.hypot(dy);
    if length <= 1e-9 {
        return default_serpentine_shape(channel_width_mm);
    }

    let nx = -dy / length;
    let ny = dx / length;
    let mut offsets = path
        .iter()
        .map(|point| ((point.0 - start.0) * nx) + ((point.1 - start.1) * ny));
    let Some(mut previous_offset) = offsets.next() else {
        return default_serpentine_shape(channel_width_mm);
    };
    let Some(mut current_offset) = offsets.next() else {
        return default_serpentine_shape(channel_width_mm);
    };

    let extrema_threshold = channel_width_mm * 0.15;
    let mut turns = 0usize;
    for next_offset in offsets {
        let previous_delta = current_offset - previous_offset;
        let next_delta = next_offset - current_offset;
        if previous_delta.abs() > 1e-9
            && next_delta.abs() > 1e-9
            && previous_delta.signum() != next_delta.signum()
            && current_offset.abs() > extrema_threshold
        {
            turns += 1;
        }
        previous_offset = current_offset;
        current_offset = next_offset;
    }

    let min_radius_mm = path
        .iter()
        .enumerate()
        .filter_map(|(idx, _)| estimate_local_bend_radius(path, idx))
        .filter(|radius| radius.is_finite())
        .fold(f64::INFINITY, f64::min);
    let bend_radius_mm = if min_radius_mm.is_finite() {
        min_radius_mm.max(channel_width_mm * 0.5)
    } else {
        channel_width_mm * 0.5
    };

    ChannelShape::Serpentine {
        segments: turns.saturating_add(1).max(2),
        bend_radius_m: Length::from_base(bend_radius_mm * 1.0e-3),
        wave_type: crate::topology::SerpentineWaveType::Sine,
    }
}

#[cfg(test)]
mod tests {
    use super::{infer_serpentine_shape, ChannelShape};

    #[test]
    fn rolling_offset_scan_counts_opposing_turns() {
        let path = [(0.0, 0.0), (3.0, 2.0), (5.0, 0.0), (7.0, -2.0), (10.0, 0.0)];
        let shape = infer_serpentine_shape(&path, path[0], path[4], 1.0);

        let ChannelShape::Serpentine {
            segments,
            bend_radius_m,
            ..
        } = shape
        else {
            panic!("shape inference must return a serpentine shape");
        };
        assert_eq!(segments, 3);
        assert!(bend_radius_m.into_base() > 0.0);
    }

    #[test]
    fn short_and_degenerate_paths_keep_the_default_shape() {
        let short_path = [(0.0, 0.0), (10.0, 0.0)];
        let degenerate_path = [(0.0, 0.0), (0.0, 1.0), (0.0, 2.0)];

        assert!(matches!(
            infer_serpentine_shape(&short_path, short_path[0], short_path[1], 2.0),
            ChannelShape::Serpentine { segments: 2, .. }
        ));
        assert!(matches!(
            infer_serpentine_shape(
                &degenerate_path,
                degenerate_path[0],
                degenerate_path[0],
                2.0
            ),
            ChannelShape::Serpentine { segments: 2, .. }
        ));
    }
}
