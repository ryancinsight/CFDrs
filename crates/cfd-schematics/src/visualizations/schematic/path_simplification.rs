use crate::geometry::{ChannelTypeCategory, Point2D};

pub(super) fn render_path_for_display(
    path: &[Point2D],
    category: ChannelTypeCategory,
    box_dims: (f64, f64),
) -> Vec<Point2D> {
    if path.len() <= 2 {
        return path.to_vec();
    }

    let tolerance_mm = (box_dims.0.max(box_dims.1) / 1200.0).clamp(0.10, 0.40);
    let simplified = simplify_polyline(path, tolerance_mm);

    if !matches!(category, ChannelTypeCategory::Curved) || simplified.len() < 4 {
        return simplified;
    }

    let smoothed = catmull_rom_resample(&simplified, 5);
    simplify_polyline(&smoothed, tolerance_mm * 0.75)
}

fn simplify_polyline(path: &[Point2D], tolerance_mm: f64) -> Vec<Point2D> {
    if path.len() <= 2 {
        return path.to_vec();
    }

    let mut keep = vec![false; path.len()];
    keep[0] = true;
    keep[path.len() - 1] = true;
    simplify_polyline_range(path, 0, path.len() - 1, tolerance_mm, &mut keep);

    let simplified: Vec<Point2D> = path
        .iter()
        .zip(keep)
        .filter_map(|(point, keep_point)| keep_point.then_some(*point))
        .collect();

    if simplified.len() < 2 {
        vec![path[0], *path.last().expect("path has at least one point")]
    } else {
        simplified
    }
}

fn simplify_polyline_range(
    path: &[Point2D],
    start_idx: usize,
    end_idx: usize,
    tolerance_mm: f64,
    keep: &mut [bool],
) {
    if end_idx <= start_idx + 1 {
        return;
    }

    let start = path[start_idx];
    let end = path[end_idx];
    let mut max_distance = 0.0;
    let mut max_idx = None;

    for (offset, point) in path[start_idx + 1..end_idx].iter().enumerate() {
        let idx = start_idx + 1 + offset;
        let distance = perpendicular_distance_to_segment(*point, start, end);
        if distance > max_distance {
            max_distance = distance;
            max_idx = Some(idx);
        }
    }

    if max_distance > tolerance_mm {
        if let Some(idx) = max_idx {
            keep[idx] = true;
            simplify_polyline_range(path, start_idx, idx, tolerance_mm, keep);
            simplify_polyline_range(path, idx, end_idx, tolerance_mm, keep);
        }
    }
}

fn perpendicular_distance_to_segment(point: Point2D, seg_start: Point2D, seg_end: Point2D) -> f64 {
    let (x0, y0) = point;
    let (x1, y1) = seg_start;
    let (x2, y2) = seg_end;
    let dx = x2 - x1;
    let dy = y2 - y1;

    if dx.abs() < f64::EPSILON && dy.abs() < f64::EPSILON {
        return ((x0 - x1).powi(2) + (y0 - y1).powi(2)).sqrt();
    }

    let t = (((x0 - x1) * dx) + ((y0 - y1) * dy)) / (dx * dx + dy * dy);
    let t = t.clamp(0.0, 1.0);
    let proj_x = x1 + t * dx;
    let proj_y = y1 + t * dy;
    ((x0 - proj_x).powi(2) + (y0 - proj_y).powi(2)).sqrt()
}

fn catmull_rom_resample(path: &[Point2D], samples_per_segment: usize) -> Vec<Point2D> {
    if path.len() < 4 || samples_per_segment == 0 {
        return path.to_vec();
    }

    let mut out = Vec::with_capacity((path.len() - 1) * samples_per_segment + 1);
    out.push(path[0]);

    for idx in 0..path.len() - 1 {
        let p0 = if idx == 0 { path[idx] } else { path[idx - 1] };
        let p1 = path[idx];
        let p2 = path[idx + 1];
        let p3 = if idx + 2 < path.len() {
            path[idx + 2]
        } else {
            path[idx + 1]
        };

        for sample_idx in 1..=samples_per_segment {
            let t = sample_idx as f64 / samples_per_segment as f64;
            out.push(catmull_rom_point(p0, p1, p2, p3, t));
        }
    }

    out
}

fn catmull_rom_point(p0: Point2D, p1: Point2D, p2: Point2D, p3: Point2D, t: f64) -> Point2D {
    let t2 = t * t;
    let t3 = t2 * t;
    let x = 0.5
        * ((2.0 * p1.0)
            + (-p0.0 + p2.0) * t
            + (2.0 * p0.0 - 5.0 * p1.0 + 4.0 * p2.0 - p3.0) * t2
            + (-p0.0 + 3.0 * p1.0 - 3.0 * p2.0 + p3.0) * t3);
    let y = 0.5
        * ((2.0 * p1.1)
            + (-p0.1 + p2.1) * t
            + (2.0 * p0.1 - 5.0 * p1.1 + 4.0 * p2.1 - p3.1) * t2
            + (-p0.1 + 3.0 * p1.1 - 3.0 * p2.1 + p3.1) * t3);
    (x, y)
}