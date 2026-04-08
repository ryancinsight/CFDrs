//! Pure path and polyline geometry helpers for selective-tree channel routing.
//!
//! All functions in this module are stateless geometric algorithms operating on
//! 2-D point sequences. No blueprint or topology types are referenced here.

use super::super::super::types::Point2D;
use crate::domain::model::ChannelShape;
use super::CenterSerpentinePathSpec;

/// Compute the total path length (mm) of a sequence of 2-D waypoints.
pub(super) fn polyline_length_mm(points: &[Point2D]) -> f64 {
    points
        .windows(2)
        .map(|segment| {
            let dx = segment[1].0 - segment[0].0;
            let dy = segment[1].1 - segment[0].1;
            dx.hypot(dy)
        })
        .sum()
}

/// Remove duplicate consecutive points and collapse three co-linear points into
/// two by dropping the middle one.
pub(super) fn simplify_polyline_points(points: Vec<Point2D>) -> Vec<Point2D> {
    let mut simplified = Vec::with_capacity(points.len());
    for point in points {
        if simplified.last().is_some_and(|last: &Point2D| {
            (last.0 - point.0).abs() < 1e-9 && (last.1 - point.1).abs() < 1e-9
        }) {
            continue;
        }
        simplified.push(point);
    }
    if simplified.len() < 3 {
        return simplified;
    }

    let mut cleaned = Vec::with_capacity(simplified.len());
    for point in simplified {
        cleaned.push(point);
        while cleaned.len() >= 3 {
            let len = cleaned.len();
            if points_are_colinear(cleaned[len - 3], cleaned[len - 2], cleaned[len - 1]) {
                cleaned.remove(len - 2);
            } else {
                break;
            }
        }
    }
    cleaned
}

/// Return `true` when `b` lies on the line segment `a–c` (within floating-point
/// tolerance).
pub(super) fn points_are_colinear(a: Point2D, b: Point2D, c: Point2D) -> bool {
    let cross = (b.0 - a.0) * (c.1 - b.1) - (b.1 - a.1) * (c.0 - b.0);
    cross.abs() < 1e-9
}

/// Test whether two open line segments `p1–p2` and `p3–p4` strictly intersect.
/// Endpoints touching are not counted as intersections.
pub(super) fn segment_intersection_strict(
    p1: Point2D,
    p2: Point2D,
    p3: Point2D,
    p4: Point2D,
) -> bool {
    let dx1 = p2.0 - p1.0;
    let dy1 = p2.1 - p1.1;
    let dx2 = p4.0 - p3.0;
    let dy2 = p4.1 - p3.1;
    let denom = dx1 * dy2 - dy1 * dx2;
    if denom.abs() < 1e-12 {
        return false;
    }

    let t = ((p3.0 - p1.0) * dy2 - (p3.1 - p1.1) * dx2) / denom;
    let u = ((p3.0 - p1.0) * dy1 - (p3.1 - p1.1) * dx1) / denom;
    let eps = 1e-6;
    t > eps && t < (1.0 - eps) && u > eps && u < (1.0 - eps)
}

/// Return `true` when `candidate` intersects any segment in any of
/// `existing_paths`.
pub(super) fn path_intersects_any(
    candidate: &[Point2D],
    existing_paths: &[Vec<Point2D>],
) -> bool {
    candidate.windows(2).any(|segment_a| {
        existing_paths.iter().any(|path| {
            path.windows(2).any(|segment_b| {
                segment_intersection_strict(segment_a[0], segment_a[1], segment_b[0], segment_b[1])
            })
        })
    })
}

/// Build a three-point dogleg polyline from `start` to `end` passing through
/// the intermediate horizontal lane at `lane_y`.
pub(super) fn monotone_dogleg_path(start: Point2D, end: Point2D, lane_y: f64) -> Vec<Point2D> {
    let span_x = end.0 - start.0;
    let lead_frac = if span_x.abs() > 1e-9 { 0.25 } else { 0.0 };
    let x1 = start.0 + span_x * lead_frac;
    let x2 = end.0 - span_x * lead_frac;
    vec![
        start,
        (x1, start.1),
        (x1, lane_y),
        (x2, lane_y),
        (x2, end.1),
        end,
    ]
}

/// Generate a multi-lobe sinusoidal serpentine path between `start` and `end`.
///
/// # Theorem
/// The path visits `spec.segments` lobes whose amplitude is bounded by
/// `amplitude_mm` subject to the physical bend radius `spec.bend_radius_m`.
pub(super) fn build_serpentine_lobe_path(
    start: Point2D,
    end: Point2D,
    width_m: f64,
    spec: CenterSerpentinePathSpec,
    amplitude_mm: f64,
) -> Vec<Point2D> {
    let dx = end.0 - start.0;
    let dy = end.1 - start.1;
    let chord_length_mm = dx.hypot(dy);
    if chord_length_mm <= 1e-9 || spec.segments < 2 {
        return vec![start, end];
    }

    let channel_diameter_mm = (width_m * 1.0e3).max(1.0e-3);
    let bend_radius_mm = (spec.bend_radius_m * 1.0e3).max(channel_diameter_mm * 0.5);
    let requested_turns = spec.segments.saturating_sub(1).max(2);
    // Always materialise at least a mirrored 2-lobe rotated-S when a treatment
    // serpentine is requested. The physical bend radius still modulates the
    // shoulder and amplitude, but the authored path never degrades to a single apex.
    let safe_turns = requested_turns;
    let step_mm = chord_length_mm / safe_turns as f64;
    let shoulder_mm = (step_mm * 0.28)
        .max(channel_diameter_mm * 0.75)
        .min(bend_radius_mm.max(channel_diameter_mm * 1.5));
    let lobe_amplitude_mm = amplitude_mm
        .max(channel_diameter_mm * 1.5)
        .min(step_mm * 0.40);
    let tx = dx / chord_length_mm;
    let ty = dy / chord_length_mm;
    let nx = -ty;
    let ny = tx;

    let mut path = Vec::with_capacity(2 + safe_turns * 5);
    path.push(start);
    let mut last_along_mm = 0.0;
    for turn_idx in 0..safe_turns {
        let center_mm = step_mm * (turn_idx as f64 + 0.5);
        let direction = if turn_idx % 2 == 0 { 1.0 } else { -1.0 };
        for (along_mm, offset_scale) in [
            (center_mm - shoulder_mm, 0.35),
            (center_mm - shoulder_mm * 0.45, 0.72),
            (center_mm, 1.0),
            (center_mm + shoulder_mm * 0.45, 0.72),
            (center_mm + shoulder_mm, 0.35),
        ] {
            if along_mm <= last_along_mm + 1.0e-6 || along_mm >= chord_length_mm - 1.0e-6 {
                continue;
            }
            let offset_mm = direction * lobe_amplitude_mm * offset_scale;
            path.push((
                start.0 + tx * along_mm + nx * offset_mm,
                start.1 + ty * along_mm + ny * offset_mm,
            ));
            last_along_mm = along_mm;
        }
    }
    path.push(end);
    path
}

/// Estimate the local circumradius at vertex `idx` in `path` using adjacent
/// triangle circumscribed-circle formula.  Returns `None` at path endpoints.
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

/// Infer a `ChannelShape::Serpentine` from an existing polyline by counting
/// direction reversals and measuring the minimum local bend radius.
pub(super) fn infer_serpentine_shape(
    path: &[Point2D],
    start: Point2D,
    end: Point2D,
    channel_width_mm: f64,
) -> ChannelShape {
    if path.len() < 3 {
        return ChannelShape::Serpentine {
            segments: 2,
            bend_radius_m: (channel_width_mm * 0.5) * 1.0e-3,
            wave_type: crate::topology::SerpentineWaveType::Sine,
        };
    }

    let dx = end.0 - start.0;
    let dy = end.1 - start.1;
    let length = dx.hypot(dy);
    if length <= 1e-9 {
        return ChannelShape::Serpentine {
            segments: 2,
            bend_radius_m: (channel_width_mm * 0.5) * 1.0e-3,
            wave_type: crate::topology::SerpentineWaveType::Sine,
        };
    }

    let nx = -dy / length;
    let ny = dx / length;
    let offsets: Vec<f64> = path
        .iter()
        .map(|point| ((point.0 - start.0) * nx) + ((point.1 - start.1) * ny))
        .collect();
    let extrema_threshold = channel_width_mm * 0.15;
    let mut turns = 0usize;
    for idx in 1..offsets.len() - 1 {
        let prev_delta = offsets[idx] - offsets[idx - 1];
        let next_delta = offsets[idx + 1] - offsets[idx];
        if prev_delta.abs() <= 1e-9 || next_delta.abs() <= 1e-9 {
            continue;
        }
        if prev_delta.signum() != next_delta.signum() && offsets[idx].abs() > extrema_threshold {
            turns += 1;
        }
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
        bend_radius_m: bend_radius_mm * 1.0e-3,
        wave_type: crate::topology::SerpentineWaveType::Sine,
    }
}

/// Build an overlay serpentine path from the endpoints of an existing polyline.
pub(super) fn serpentine_overlay_path(
    points: &[Point2D],
    width_m: f64,
    spec: CenterSerpentinePathSpec,
    amplitude_mm: f64,
) -> Vec<Point2D> {
    let Some(start) = points.first().copied() else {
        return Vec::new();
    };
    let Some(end) = points.last().copied() else {
        return vec![start];
    };
    build_serpentine_lobe_path(start, end, width_m, spec, amplitude_mm)
}
