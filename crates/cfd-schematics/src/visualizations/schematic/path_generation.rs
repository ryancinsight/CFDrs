use crate::geometry::Point2D;

pub(super) fn generated_parallel_path(
    start: Point2D,
    end: Point2D,
    slot: usize,
    count: usize,
    box_dims: (f64, f64),
) -> Vec<Point2D> {
    if count <= 1 {
        return vec![start, end];
    }

    let dx = end.0 - start.0;
    let dy = end.1 - start.1;
    let length = dx.hypot(dy);
    if length <= 1e-9 {
        return vec![start, end];
    }

    let nx = -dy / length;
    let ny = dx / length;
    let base_spacing_mm = (length * 0.18).clamp(6.0, box_dims.1 * 0.18);
    let centered_slot = slot as f64 - (count as f64 - 1.0) * 0.5;
    let spread = if count == 2 { 2.0 } else { 1.0 };
    let offset_mm = centered_slot * base_spacing_mm * spread;
    if offset_mm.abs() <= 1e-9 {
        return vec![start, end];
    }

    let shoulder = 0.28;
    let entry = (
        start.0 + dx * shoulder + nx * offset_mm,
        start.1 + dy * shoulder + ny * offset_mm,
    );
    let exit = (
        end.0 - dx * shoulder + nx * offset_mm,
        end.1 - dy * shoulder + ny * offset_mm,
    );
    vec![start, entry, exit, end]
}

pub(super) fn path_has_visible_serpentine_curvature(path: &[Point2D], bend_radius_mm: f64) -> bool {
    if path.len() < 3 {
        return false;
    }

    let start = path[0];
    let end = *path.last().expect("path has at least one point");
    let dx = end.0 - start.0;
    let dy = end.1 - start.1;
    let chord = dx.hypot(dy);
    if chord <= 1.0e-9 {
        return false;
    }

    let nx = -dy / chord;
    let ny = dx / chord;
    let min_expected_offset_mm = (bend_radius_mm * 0.25).max(0.75);
    let offsets = path
        .iter()
        .map(|point| ((point.0 - start.0) * nx) + ((point.1 - start.1) * ny))
        .filter(|offset| offset.abs() >= min_expected_offset_mm)
        .collect::<Vec<_>>();

    offsets.iter().any(|offset| *offset > 0.0) && offsets.iter().any(|offset| *offset < 0.0)
}

/// Generates a square-wave serpentine path between two endpoints.
///
/// Each segment produces one horizontal run offset alternately +/- amplitude
/// from the center-line (start->end), connected by perpendicular 90-degree steps.
/// This gives a clear square-wave visual distinct from a smooth sine curve.
pub(super) fn generated_serpentine_path(
    start: Point2D,
    end: Point2D,
    segments: usize,
    bend_radius_mm: f64,
) -> Vec<Point2D> {
    if segments < 2 {
        return vec![start, end];
    }

    let dx = end.0 - start.0;
    let dy = end.1 - start.1;
    let length = (dx * dx + dy * dy).sqrt();
    if length < 1e-6 {
        return vec![start, end];
    }

    let px = dx / length;
    let py = dy / length;
    let nx = -py;
    let ny = px;
    let requested_turns = segments.saturating_sub(1).max(2);
    let step = length / requested_turns as f64;
    let shoulder = (step * 0.28).max(length * 0.04);
    let amplitude = (bend_radius_mm * 1.4)
        .max(length / ((requested_turns as f64 + 1.0) * 3.0))
        .min(length * 0.24);

    let mut path = Vec::with_capacity(2 + requested_turns * 5);
    path.push(start);
    for turn_idx in 0..requested_turns {
        let center = step * (turn_idx as f64 + 0.5);
        let sign = if turn_idx % 2 == 0 { 1.0 } else { -1.0 };
        for (along, scale) in [
            (center - shoulder, 0.35),
            (center - shoulder * 0.45, 0.72),
            (center, 1.0),
            (center + shoulder * 0.45, 0.72),
            (center + shoulder, 0.35),
        ] {
            if along <= 1.0e-6 || along >= length - 1.0e-6 {
                continue;
            }
            let offset = sign * amplitude * scale;
            path.push((
                start.0 + px * along + nx * offset,
                start.1 + py * along + ny * offset,
            ));
        }
    }
    path.push(end);
    path
}
