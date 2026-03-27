//! Measurement overlay — renders dimension annotations as overlay geometry.
//!
//! Converts `MeasurementResult` entries into `OverlayLine` segments for
//! dimension arrows and leader lines, rendered via the existing `OverlayPipeline`.

use crate::domain::measurement::{MeasurementKind, MeasurementResult, MeasurementStore};
use crate::infrastructure::gpu::overlay_buffer::{OverlayBuffer, OverlayLine};

/// Color for measurement dimension lines.
const MEASURE_COLOR: [f32; 4] = [1.0, 0.8, 0.0, 0.9];

/// Color for measurement extension lines.
const EXTENSION_COLOR: [f32; 4] = [0.7, 0.7, 0.7, 0.6];

/// Perpendicular offset for dimension lines (world units).
const DIM_OFFSET: f32 = 0.15;

/// Generate overlay lines for all active measurements.
#[must_use]
pub fn measurement_overlay_lines(store: &MeasurementStore) -> Vec<OverlayLine> {
    let mut lines = Vec::new();
    for result in store.iter() {
        lines.extend(lines_for_measurement(result));
    }
    lines
}

/// Generate overlay lines for a single measurement.
fn lines_for_measurement(result: &MeasurementResult) -> Vec<OverlayLine> {
    match &result.kind {
        MeasurementKind::PointToPoint { start, end } => {
            point_to_point_lines(start, end)
        }
        MeasurementKind::EdgeLength { v0, v1 } => {
            point_to_point_lines(v0, v1)
        }
        MeasurementKind::FaceAngle { edge_midpoint, .. } => {
            // Small cross-hair at the edge midpoint.
            cross_hair_lines(edge_midpoint)
        }
        MeasurementKind::FaceArea { face_center } => {
            cross_hair_lines(face_center)
        }
        MeasurementKind::MeshVolume { mesh_center } => {
            cross_hair_lines(mesh_center)
        }
    }
}

/// Dimension lines for a two-point measurement (distance or edge length).
fn point_to_point_lines(start: &[f64; 3], end: &[f64; 3]) -> Vec<OverlayLine> {
    let s = arr_f32(start);
    let e = arr_f32(end);

    // Main dimension line.
    let main_line = OverlayLine {
        start: s,
        end: e,
        color: MEASURE_COLOR,
    };

    // Extension lines (small ticks perpendicular to the dimension line).
    let dir = [e[0] - s[0], e[1] - s[1], e[2] - s[2]];
    let len = (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]).sqrt();
    if len < 1e-8 {
        return vec![main_line];
    }

    // Find a perpendicular direction for ticks.
    let perp = perpendicular_direction(&dir);
    let tick_s = [
        s[0] - perp[0] * DIM_OFFSET,
        s[1] - perp[1] * DIM_OFFSET,
        s[2] - perp[2] * DIM_OFFSET,
    ];
    let tick_e_start = [
        s[0] + perp[0] * DIM_OFFSET,
        s[1] + perp[1] * DIM_OFFSET,
        s[2] + perp[2] * DIM_OFFSET,
    ];
    let tick_s2 = [
        e[0] - perp[0] * DIM_OFFSET,
        e[1] - perp[1] * DIM_OFFSET,
        e[2] - perp[2] * DIM_OFFSET,
    ];
    let tick_e2 = [
        e[0] + perp[0] * DIM_OFFSET,
        e[1] + perp[1] * DIM_OFFSET,
        e[2] + perp[2] * DIM_OFFSET,
    ];

    vec![
        main_line,
        OverlayLine { start: tick_s, end: tick_e_start, color: EXTENSION_COLOR },
        OverlayLine { start: tick_s2, end: tick_e2, color: EXTENSION_COLOR },
    ]
}

/// Cross-hair lines at a single point (for area/volume/angle indicators).
fn cross_hair_lines(center: &[f64; 3]) -> Vec<OverlayLine> {
    let c = arr_f32(center);
    let d = DIM_OFFSET;
    vec![
        OverlayLine {
            start: [c[0] - d, c[1], c[2]],
            end: [c[0] + d, c[1], c[2]],
            color: MEASURE_COLOR,
        },
        OverlayLine {
            start: [c[0], c[1] - d, c[2]],
            end: [c[0], c[1] + d, c[2]],
            color: MEASURE_COLOR,
        },
    ]
}

/// Build GPU overlay buffer for measurement lines.
pub fn build_measure_overlay(store: &MeasurementStore, device: &wgpu::Device) -> OverlayBuffer {
    let lines = measurement_overlay_lines(store);
    OverlayBuffer::from_lines(&lines, device)
}

/// Find a unit vector perpendicular to the given direction.
fn perpendicular_direction(dir: &[f32; 3]) -> [f32; 3] {
    let (ax, ay, az) = (dir[0].abs(), dir[1].abs(), dir[2].abs());
    let seed = if ax < ay && ax < az {
        [1.0, 0.0, 0.0]
    } else if ay < az {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    // cross product dir x seed
    let cx = dir[1] * seed[2] - dir[2] * seed[1];
    let cy = dir[2] * seed[0] - dir[0] * seed[2];
    let cz = dir[0] * seed[1] - dir[1] * seed[0];
    let len = (cx * cx + cy * cy + cz * cz).sqrt();
    if len < 1e-8 {
        return [0.0, 1.0, 0.0];
    }
    [cx / len, cy / len, cz / len]
}

/// Convert [f64; 3] to [f32; 3].
fn arr_f32(a: &[f64; 3]) -> [f32; 3] {
    [a[0] as f32, a[1] as f32, a[2] as f32]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn point_to_point_generates_three_lines() {
        let lines = point_to_point_lines(&[0.0, 0.0, 0.0], &[1.0, 0.0, 0.0]);
        assert_eq!(lines.len(), 3); // main + 2 ticks
    }

    #[test]
    fn cross_hair_generates_two_lines() {
        let lines = cross_hair_lines(&[0.0, 0.0, 0.0]);
        assert_eq!(lines.len(), 2);
    }

    #[test]
    fn empty_store_produces_empty_overlay() {
        let store = MeasurementStore::new();
        let lines = measurement_overlay_lines(&store);
        assert!(lines.is_empty());
    }
}
