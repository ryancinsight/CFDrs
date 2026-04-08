use crate::domain::model::{ChannelSpec, NetworkBlueprint};
use crate::geometry::Point2D;
use std::collections::HashMap;

/// Test whether two line segments (p1,p2) and (p3,p4) properly intersect.
pub(super) fn segment_intersection(
    p1: Point2D,
    p2: Point2D,
    p3: Point2D,
    p4: Point2D,
) -> Option<(f64, f64, Point2D)> {
    let dx1 = p2.0 - p1.0;
    let dy1 = p2.1 - p1.1;
    let dx2 = p4.0 - p3.0;
    let dy2 = p4.1 - p3.1;

    let denom = dx1 * dy2 - dy1 * dx2;
    if denom.abs() < 1e-12 {
        return None;
    }

    let t = ((p3.0 - p1.0) * dy2 - (p3.1 - p1.1) * dx2) / denom;
    let u = ((p3.0 - p1.0) * dy1 - (p3.1 - p1.1) * dx1) / denom;

    let eps = 1e-6;
    if t > eps && t < (1.0 - eps) && u > eps && u < (1.0 - eps) {
        let ix = t.mul_add(dx1, p1.0);
        let iy = t.mul_add(dy1, p1.1);
        Some((t, u, (ix, iy)))
    } else {
        None
    }
}

/// Detected intersection between two channel segments.
#[derive(Debug, Clone)]
pub(super) struct SegmentCrossing {
    pub(super) channel_a: usize,
    pub(super) seg_a: usize,
    pub(super) t_a: f64,
    pub(super) channel_b: usize,
    pub(super) seg_b: usize,
    pub(super) t_b: f64,
    pub(super) point: Point2D,
}

fn push_distinct_point(points: &mut Vec<Point2D>, point: Point2D) {
    let Some(last) = points.last().copied() else {
        points.push(point);
        return;
    };
    if (last.0 - point.0).abs() > 1e-9 || (last.1 - point.1).abs() > 1e-9 {
        points.push(point);
    }
}

fn channel_centerline_points(
    channel: &ChannelSpec,
    node_points: &HashMap<String, Point2D>,
) -> Vec<Point2D> {
    match channel.path.as_slice() {
        [] => {
            let mut centerline = Vec::with_capacity(2);
            if let Some(start) = node_points.get(channel.from.as_str()).copied() {
                push_distinct_point(&mut centerline, start);
            }
            if let Some(end) = node_points.get(channel.to.as_str()).copied() {
                push_distinct_point(&mut centerline, end);
            }
            centerline
        }
        [midpoint] => {
            let mut centerline = Vec::with_capacity(3);
            if let Some(start) = node_points.get(channel.from.as_str()).copied() {
                push_distinct_point(&mut centerline, start);
            }
            push_distinct_point(&mut centerline, *midpoint);
            if let Some(end) = node_points.get(channel.to.as_str()).copied() {
                push_distinct_point(&mut centerline, end);
            }
            centerline
        }
        points => points.to_vec(),
    }
}

/// Detect all pairwise intersections between channel centerlines.
pub(super) fn detect_crossings(system: &NetworkBlueprint) -> Vec<SegmentCrossing> {
    let mut crossings = Vec::new();
    let node_points: HashMap<String, Point2D> = system
        .nodes
        .iter()
        .map(|node| (node.id.as_str().to_string(), node.point))
        .collect();
    let centerlines: Vec<Vec<Point2D>> = system
        .channels
        .iter()
        .map(|channel| channel_centerline_points(channel, &node_points))
        .collect();

    for (i, ch_a) in system.channels.iter().enumerate() {
        for (j, ch_b) in system.channels.iter().enumerate() {
            if j <= i {
                continue;
            }
            if ch_a.from == ch_b.from
                || ch_a.from == ch_b.to
                || ch_a.to == ch_b.from
                || ch_a.to == ch_b.to
            {
                continue;
            }

            let cl_a = &centerlines[i];
            let cl_b = &centerlines[j];
            if cl_a.len() < 2 || cl_b.len() < 2 {
                continue;
            }

            for seg_a in 0..cl_a.len().saturating_sub(1) {
                for seg_b in 0..cl_b.len().saturating_sub(1) {
                    if let Some((t, u, point)) =
                        segment_intersection(cl_a[seg_a], cl_a[seg_a + 1], cl_b[seg_b], cl_b[seg_b + 1])
                    {
                        crossings.push(SegmentCrossing {
                            channel_a: i,
                            seg_a,
                            t_a: t,
                            channel_b: j,
                            seg_b,
                            t_b: u,
                            point,
                        });
                    }
                }
            }
        }
    }

    crossings
}

/// Count unresolved centerline intersections without mutating the blueprint.
#[must_use]
pub fn unresolved_intersection_count(system: &NetworkBlueprint) -> usize {
    detect_crossings(system).len()
}

/// Report whether the blueprint still contains unresolved centerline crossings.
#[must_use]
pub fn has_unresolved_intersections(system: &NetworkBlueprint) -> bool {
    unresolved_intersection_count(system) > 0
}