//! Treatment-path routing helpers for the primitive selective tree annotator.
//!
//! Provides venturi-path conflict-resolution routing, monotone dogleg routing,
//! and DFS channel-path discovery across a node–adjacency graph.

use std::collections::{HashMap, HashSet};

use super::super::super::types::Point2D;
use crate::domain::model::{ChannelSpec, NodeId};
use super::path_geometry::{
    monotone_dogleg_path, path_intersects_any, polyline_length_mm, simplify_polyline_points,
};
use super::PendingVenturiPath;

/// Compute the midpoint Y of a channel's key waypoints, falling back to
/// `y_fallback` when no points are available.
pub(super) fn preferred_treatment_lane_y(
    points: &[Point2D],
    start: Option<Point2D>,
    end: Option<Point2D>,
    y_fallback: f64,
) -> f64 {
    match (start, end) {
        (Some(a), Some(b)) => f64::midpoint(a.1, b.1),
        (Some(a), None) => a.1,
        (None, Some(b)) => b.1,
        (None, None) => match (points.first(), points.last()) {
            (Some(first), Some(last)) => f64::midpoint(first.1, last.1),
            (Some(point), None) | (None, Some(point)) => point.1,
            (None, None) => y_fallback,
        },
    }
}

/// Assign resolved non-intersecting paths to the pending venturi channels.
pub(super) fn route_pending_venturi_paths(
    channels: &mut [ChannelSpec],
    pending_paths: &[PendingVenturiPath],
    mid_y: f64,
) {
    use std::cmp::Ordering;

    let mut ordered_paths = pending_paths.to_vec();
    ordered_paths.sort_by(|left, right| {
        right
            .preferred_y
            .partial_cmp(&left.preferred_y)
            .unwrap_or(Ordering::Equal)
            .then(left.channel_idx.cmp(&right.channel_idx))
    });

    let mut assigned_paths = Vec::with_capacity(ordered_paths.len());
    for pending in ordered_paths {
        let routed_path = route_monotone_treatment_path(
            pending.start,
            pending.end,
            pending.preferred_y,
            mid_y,
            &assigned_paths,
        );
        let channel = &mut channels[pending.channel_idx];
        channel.path = routed_path;
        channel.length_m = channel_length_from_points_or_endpoints(
            &channel.path,
            pending.start,
            pending.end,
            pending.fallback_length_m,
        );
        assigned_paths.push(channel.path.clone());
    }
}

/// Route a single treatment segment from `start` to `end` without crossing any
/// of the `existing_paths`.  A monotone dogleg is attempted up to 8 times with
/// increasing lateral offsets before falling back to a direct segment.
pub(super) fn route_monotone_treatment_path(
    start: Option<Point2D>,
    end: Option<Point2D>,
    preferred_y: f64,
    mid_y: f64,
    existing_paths: &[Vec<Point2D>],
) -> Vec<Point2D> {
    match (start, end) {
        (Some(a), Some(b)) => {
            let direct = simplify_polyline_points(vec![a, b]);
            if !path_intersects_any(&direct, existing_paths) {
                return direct;
            }

            let direction = if preferred_y >= mid_y { 1.0 } else { -1.0 };
            let offset_step = ((a.1 - b.1).abs() * 0.5).max(2.0);
            for attempt in 0..8 {
                let lane_y = preferred_y + direction * f64::from(attempt) * offset_step;
                let dogleg =
                    simplify_polyline_points(monotone_dogleg_path(a, b, lane_y));
                if !path_intersects_any(&dogleg, existing_paths) {
                    return dogleg;
                }
            }

            direct
        }
        (Some(a), None) => vec![a],
        (None, Some(b)) => vec![b],
        (None, None) => Vec::new(),
    }
}

/// Return the length (m) of a channel, preferring the polyline length over
/// endpoint distance, and falling back to `fallback_m`.
pub(super) fn channel_length_from_points_or_endpoints(
    points: &[Point2D],
    start: Option<Point2D>,
    end: Option<Point2D>,
    fallback_m: f64,
) -> f64 {
    let path_length_m = if points.len() >= 2 {
        polyline_length_mm(points) * 1.0e-3
    } else {
        0.0
    };
    if path_length_m.is_finite() && path_length_m > 0.0 {
        return path_length_m;
    }

    if let (Some(start), Some(end)) = (start, end) {
        let endpoint_length_m = ((end.0 - start.0).hypot(end.1 - start.1)) * 1.0e-3;
        if endpoint_length_m.is_finite() && endpoint_length_m > 0.0 {
            return endpoint_length_m;
        }
    }

    if fallback_m.is_finite() && fallback_m > 0.0 {
        fallback_m
    } else {
        f64::EPSILON
    }
}

/// Determine which leaf indices (0-based, sorted by Y ascending) are treatment
/// path terminals for a given first-split kind.
pub(super) fn primitive_treatment_leaf_indices(
    leaf_count: usize,
    first_stage: Option<super::PrimitiveSelectiveSplitKind>,
) -> Vec<usize> {
    use super::PrimitiveSelectiveSplitKind;
    if leaf_count == 0 {
        return Vec::new();
    }
    match first_stage {
        Some(PrimitiveSelectiveSplitKind::Bi) => {
            // Leaves sorted by Y ascending: index 0 = lower (bypass), index 1 = upper (treatment)
            vec![leaf_count - 1]
        }
        Some(PrimitiveSelectiveSplitKind::Tri) => {
            if leaf_count == 3 {
                vec![1]
            } else {
                let block = leaf_count / 3;
                (block..(2 * block)).collect()
            }
        }
        Some(PrimitiveSelectiveSplitKind::Quad) => {
            // 4 leaves sorted by Y: indices 1,2 are the center pair (treatment)
            let quarter = leaf_count / 4;
            (quarter..(3 * quarter)).collect()
        }
        Some(PrimitiveSelectiveSplitKind::Penta) => {
            // 5 leaves sorted by Y: index 2 is the center (treatment)
            if leaf_count == 5 {
                vec![2]
            } else {
                let fifth = leaf_count / 5;
                (2 * fifth..(3 * fifth)).collect()
            }
        }
        None => Vec::new(),
    }
}

/// Discover all channel indices on the path from `start` to `goal` in the
/// given adjacency graph, optionally restricting traversal to one half-plane
/// of the design box.
pub(super) fn channel_path_between(
    start: &NodeId,
    goal: &NodeId,
    node_pts: &HashMap<NodeId, Point2D>,
    adjacency: &HashMap<NodeId, Vec<(NodeId, usize)>>,
    mid_x: Option<f64>,
    left_half: bool,
) -> Option<Vec<usize>> {
    let mut visited = HashSet::new();
    let mut path = Vec::new();
    if dfs_channel_path(
        start,
        goal,
        node_pts,
        adjacency,
        mid_x,
        left_half,
        &mut visited,
        &mut path,
    ) {
        Some(path)
    } else {
        None
    }
}

fn dfs_channel_path(
    current: &NodeId,
    goal: &NodeId,
    node_pts: &HashMap<NodeId, Point2D>,
    adjacency: &HashMap<NodeId, Vec<(NodeId, usize)>>,
    mid_x: Option<f64>,
    left_half: bool,
    visited: &mut HashSet<NodeId>,
    path: &mut Vec<usize>,
) -> bool {
    if current == goal {
        return true;
    }
    visited.insert(current.clone());
    let Some(neighbors) = adjacency.get(current) else {
        return false;
    };
    for (next, channel_idx) in neighbors {
        if visited.contains(next) {
            continue;
        }
        if let Some(mid_x) = mid_x {
            let from = node_pts[current].0;
            let to = node_pts[next].0;
            let keep = if left_half {
                from <= mid_x + 1e-6 && to <= mid_x + 1e-6
            } else {
                from >= mid_x - 1e-6 && to >= mid_x - 1e-6
            };
            if !keep {
                continue;
            }
        }
        path.push(*channel_idx);
        if dfs_channel_path(
            next, goal, node_pts, adjacency, mid_x, left_half, visited, path,
        ) {
            return true;
        }
        path.pop();
    }
    visited.remove(current);
    false
}
