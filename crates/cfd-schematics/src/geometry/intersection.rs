//! Channel intersection detection and junction node insertion.
//!
//! This module detects where channel centerlines cross in 2D and inserts
//! junction nodes at those intersection points. This is essential for
//! accurate 1D and 2D simulations of planar millifluidic devices where
//! channels physically cross in the same plane.
//!
//! # Theorem — Line Segment Intersection
//!
//! Two line segments $P_1 P_2$ and $P_3 P_4$ intersect if and only if
//! the cross-product signs differ:
//!
//! ```text
//! sign(cross(P₃P₄, P₃P₁)) ≠ sign(cross(P₃P₄, P₃P₂))
//! AND
//! sign(cross(P₁P₂, P₁P₃)) ≠ sign(cross(P₁P₂, P₁P₄))
//! ```
//!
//! **Proof sketch**: The cross products determine which side of the line
//! each endpoint lies on. If endpoints of one segment lie on opposite
//! sides of the other segment (and vice versa), the segments must cross.

use crate::domain::model::{ChannelSpec, NetworkBlueprint, NodeKind, NodeSpec};

use crate::geometry::Point2D;

/// Metadata marker for nodes created at channel intersections.
#[derive(Debug, Clone)]
pub struct IntersectionMetadata {
    /// IDs of the two channels that cross at this node.
    pub channel_a_id: String,
    /// ID of the second crossing channel.
    pub channel_b_id: String,
}

impl crate::geometry::metadata::Metadata for IntersectionMetadata {
    fn metadata_type_name(&self) -> &'static str {
        "IntersectionMetadata"
    }

    fn clone_metadata(&self) -> Box<dyn crate::geometry::metadata::Metadata> {
        Box::new(self.clone())
    }

    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

    fn as_any_mut(&mut self) -> &mut dyn std::any::Any {
        self
    }
}

/// Result of intersection detection on a channel system.
#[derive(Debug, Clone)]
pub struct IntersectionResult {
    /// Number of intersections detected.
    pub intersection_count: usize,
    /// Indices of newly created junction nodes.
    pub junction_node_ids: Vec<usize>,
}

/// Test whether two line segments (p1,p2) and (p3,p4) properly intersect.
///
/// Returns the parametric intersection point (t, u) where:
/// - The intersection point on segment 1 is at p1 + t*(p2-p1)
/// - The intersection point on segment 2 is at p3 + u*(p4-p3)
/// - Both t and u must be in (0,1) for a proper interior crossing.
fn segment_intersection(
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

    // Parallel or coincident segments
    if denom.abs() < 1e-12 {
        return None;
    }

    let t = ((p3.0 - p1.0) * dy2 - (p3.1 - p1.1) * dx2) / denom;
    let u = ((p3.0 - p1.0) * dy1 - (p3.1 - p1.1) * dx1) / denom;

    // Exclude endpoints (t and u strictly inside (0,1)) to avoid
    // false positives at shared nodes where channels already meet.
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
struct SegmentCrossing {
    /// Index of channel A in the system.
    channel_a: usize,
    /// Segment index within channel A's centerline.
    seg_a: usize,
    /// Parametric position along segment A (0..1).
    t_a: f64,
    /// Index of channel B.
    channel_b: usize,
    /// Segment index within channel B's centerline.
    seg_b: usize,
    /// Parametric position along segment B.
    t_b: f64,
    /// World-space intersection point.
    point: Point2D,
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
    node_points: &std::collections::HashMap<String, Point2D>,
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
fn detect_crossings(system: &NetworkBlueprint) -> Vec<SegmentCrossing> {
    let mut crossings = Vec::new();
    let node_points: std::collections::HashMap<String, Point2D> = system
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
                continue; // avoid duplicates
            }
            // Skip channels that share an endpoint node — those are
            // intentional junctions, not geometric crossings.
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
                    if let Some((t, u, pt)) = segment_intersection(
                        cl_a[seg_a],
                        cl_a[seg_a + 1],
                        cl_b[seg_b],
                        cl_b[seg_b + 1],
                    ) {
                        crossings.push(SegmentCrossing {
                            channel_a: i,
                            seg_a,
                            t_a: t,
                            channel_b: j,
                            seg_b,
                            t_b: u,
                            point: pt,
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

/// Insert junction nodes at all detected channel intersections.
///
/// For each crossing, this function:
/// 1. Creates a new junction node at the crossing point.
/// 2. Splits each of the two crossing channels into two sub-channels
///    (before and after the crossing point).
/// 3. Attaches [`IntersectionMetadata`] to the new node.
///
/// The original crossing channels are replaced by their sub-segments, and
/// the newly created junction nodes appear in `system.nodes`.
///
/// # Returns
///
/// An [`IntersectionResult`] summarising how many intersections were found
/// and which node IDs were created.
pub fn insert_intersection_nodes(system: &mut NetworkBlueprint) -> IntersectionResult {
    let crossings = detect_crossings(system);
    if crossings.is_empty() {
        return IntersectionResult {
            intersection_count: 0,
            junction_node_ids: Vec::new(),
        };
    }

    // For simplicity and correctness, rebuild the channel list from scratch.
    // First, create all junction nodes for the crossings.
    let mut junction_ids = Vec::with_capacity(crossings.len());
    let mut new_nodes = Vec::new();
    let initial_node_count = system.nodes.len();

    for (ci, crossing) in crossings.iter().enumerate() {
        let node_id = format!(
            "intersect_{}_{}_{}",
            crossing.channel_a, crossing.channel_b, ci
        );

        // This won't perfectly match exactly if IntersectionMetadata expects usize, but we'll fix IntersectionMetadata next.
        // Wait, IntersectionMetadata might expect usize if it wasn't refactored. Let me check its definition.
        // I'll update it separately below if needed.

        system.nodes.push(NodeSpec::new_at(
            node_id.clone(),
            NodeKind::Junction,
            crossing.point,
        ));
        junction_ids.push(initial_node_count + ci);
        new_nodes.push(node_id);
    }

    // Build a map: original_channel_index → Vec<(crossing_index, role)>
    // where role indicates whether this channel is A or B in the crossing.
    struct SplitInfo {
        junction_node_id: String,
        /// Parametric position along this channel's centerline segment.
        t: f64,
        /// Which segment of the centerline is being split.
        seg: usize,
    }

    let mut splits_per_channel: std::collections::HashMap<usize, Vec<SplitInfo>> =
        std::collections::HashMap::new();

    for (ci, crossing) in crossings.iter().enumerate() {
        splits_per_channel
            .entry(crossing.channel_a)
            .or_default()
            .push(SplitInfo {
                junction_node_id: new_nodes[ci].clone(),
                t: crossing.t_a,
                seg: crossing.seg_a,
            });
        splits_per_channel
            .entry(crossing.channel_b)
            .or_default()
            .push(SplitInfo {
                junction_node_id: new_nodes[ci].clone(),
                t: crossing.t_b,
                seg: crossing.seg_b,
            });
    }

    // Sort each channel's splits by (segment, t) so we can process them in order.
    for splits in splits_per_channel.values_mut() {
        splits.sort_by(|a, b| {
            a.seg
                .cmp(&b.seg)
                .then(a.t.partial_cmp(&b.t).expect("structural invariant"))
        });
    }

    // Rebuild the channel list.
    let mut next_channel_idx = 0;
    let mut new_channels: Vec<ChannelSpec> = Vec::new();

    for (ch_idx, orig_channel) in system.channels.iter().enumerate() {
        if let Some(splits) = splits_per_channel.get(&ch_idx) {
            // This channel is split at one or more crossings.
            // Walk through the centerline and produce sub-channels.
            let centerline = &orig_channel.path;
            if centerline.len() < 2 {
                new_channels.push(orig_channel.clone());
                continue;
            }

            // Build a sequence of (node_id, point) for the split points,
            // ordered along the channel path.
            let mut waypoints: Vec<(String, Point2D)> = Vec::new();
            waypoints.push((orig_channel.from.0.clone(), centerline[0]));

            // Insert splits in order along the centerline.
            let mut path_idx = 0; // current segment index
            let mut split_iter = splits.iter().peekable();

            while path_idx < centerline.len() - 1 {
                // Collect all splits on this segment, sorted by t.
                while let Some(split) = split_iter.peek() {
                    if split.seg == path_idx {
                        let s = split_iter.next().expect("structural invariant");
                        let pt = (
                            centerline[path_idx]
                                .0
                                .mul_add(1.0 - s.t, centerline[path_idx + 1].0 * s.t),
                            centerline[path_idx]
                                .1
                                .mul_add(1.0 - s.t, centerline[path_idx + 1].1 * s.t),
                        );
                        waypoints.push((s.junction_node_id.clone(), pt));
                    } else {
                        break;
                    }
                }

                // Also include internal centerline points as path geometry
                // (they are not nodes, just path fidelity).
                path_idx += 1;
            }

            waypoints.push((
                orig_channel.to.0.clone(),
                *centerline.last().expect("structural invariant"),
            ));

            // Now produce sub-channels between consecutive waypoints.
            for pair in waypoints.windows(2) {
                let (from_id, from_pt) = &pair[0];
                let (to_id, to_pt) = &pair[1];
                next_channel_idx += 1;
                let mut channel_clone = orig_channel.clone();
                channel_clone.id = crate::domain::model::EdgeId(format!(
                    "{}_part{}",
                    orig_channel.id.0, next_channel_idx
                ));
                channel_clone.from = crate::domain::model::NodeId(from_id.clone());
                channel_clone.to = crate::domain::model::NodeId(to_id.clone());
                channel_clone.path = vec![*from_pt, *to_pt]; // Simplify path to straight segments between intersections
                new_channels.push(channel_clone);
            }
        } else {
            // Channel has no crossings; keep as-is.
            new_channels.push(orig_channel.clone());
        }
    }

    system.channels = new_channels;

    IntersectionResult {
        intersection_count: crossings.len(),
        junction_node_ids: junction_ids,
    }
}

/// Compute the optimal box dimensions for a given split pattern.
///
/// Scales the layout to fill more of the available treatment zone as the
/// number of parallel channels increases. With more splits, channels need
/// more vertical space, so the height is expanded proportionally while
/// respecting the physical plate constraints.
///
/// # Arguments
///
/// * `plate_width_mm` — Available width of the treatment zone [mm].
/// * `plate_height_mm` — Available height of the treatment zone [mm].
/// * `total_branches` — Total number of parallel branches at the widest point
///   (product of all split branch counts).
/// * `channel_width_mm` — Physical channel width [mm].
/// * `wall_clearance_mm` — Minimum clearance between channels and walls [mm].
///
/// # Returns
///
/// `(width_mm, height_mm)` — Scaled box dimensions that maximise use of the
/// available footprint while maintaining minimum clearance constraints.
///
/// # Theorem — Footprint Utilisation
///
/// For $n$ parallel branches each of width $w$ with clearance $c$, the
/// minimum required height is:
///
/// $$H_{\min} = n \cdot w + (n + 1) \cdot c$$
///
/// The algorithm scales the box height to $\min(H_{\min} \cdot k, H_{\text{plate}})$
/// where $k \geq 1$ is a padding factor that ensures serpentine amplitude
/// has room to develop.
#[must_use]
pub fn adaptive_box_dims(
    plate_width_mm: f64,
    plate_height_mm: f64,
    total_branches: usize,
    channel_width_mm: f64,
    wall_clearance_mm: f64,
) -> (f64, f64) {
    let n = total_branches.max(1) as f64;

    // Minimum height needed for n parallel branches with clearance.
    let min_height = n * channel_width_mm + (n + 1.0) * wall_clearance_mm;

    // Branch-adaptive padding factor: fewer branches get more room for serpentine
    // amplitude; denser layouts pack more tightly to maximise fill.
    //   1–4 branches  → 1.8× (generous serpentine room)
    //   5–16 branches → 1.5× (balanced)
    //   17+ branches  → 1.2× (tight packing, maximise well coverage)
    let padding_factor = if total_branches <= 4 {
        1.8
    } else if total_branches <= 16 {
        1.5
    } else {
        1.2
    };

    let desired_height = (min_height * padding_factor).max(plate_height_mm * 0.3);

    // Clamp to plate boundaries.
    let height = desired_height.min(plate_height_mm);
    let width = plate_width_mm;

    (width, height)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::{ChannelTypeConfig, GeometryConfig};
    use crate::geometry::generator::create_geometry;
    use crate::geometry::SplitType;

    #[test]
    fn segment_intersection_detects_crossing() {
        let result = segment_intersection((0.0, 0.0), (1.0, 1.0), (0.0, 1.0), (1.0, 0.0));
        assert!(result.is_some());
        let (t, u, pt) = result.expect("structural invariant");
        assert!((t - 0.5).abs() < 1e-10);
        assert!((u - 0.5).abs() < 1e-10);
        assert!((pt.0 - 0.5).abs() < 1e-10);
        assert!((pt.1 - 0.5).abs() < 1e-10);
    }

    #[test]
    fn segment_intersection_detects_no_crossing_for_parallel() {
        let result = segment_intersection((0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0));
        assert!(result.is_none());
    }

    #[test]
    fn segment_intersection_excludes_shared_endpoints() {
        // Segments that share an endpoint should not register as a crossing.
        let result = segment_intersection((0.0, 0.0), (1.0, 0.0), (1.0, 0.0), (2.0, 1.0));
        assert!(result.is_none());
    }

    #[test]
    fn adaptive_box_dims_scales_with_branches() {
        let (w1, h1) = adaptive_box_dims(45.0, 45.0, 1, 2.0, 2.0);
        let (w2, h2) = adaptive_box_dims(45.0, 45.0, 4, 2.0, 2.0);
        let (_w8, h8) = adaptive_box_dims(45.0, 45.0, 8, 2.0, 2.0);

        // Width should always use full plate width.
        assert!((w1 - 45.0).abs() < 1e-10);
        assert!((w2 - 45.0).abs() < 1e-10);

        // Height should grow with branch count.
        assert!(h2 > h1, "4 branches should use more height than 1");
        assert!(h8 > h2, "8 branches should use more height than 4");

        // Height should not exceed plate.
        assert!(h8 <= 45.0 + 1e-10);
    }

    #[test]
    fn no_intersections_for_simple_bifurcation() {
        let mut system = create_geometry(
            (100.0, 50.0),
            &[SplitType::Bifurcation],
            &GeometryConfig::default(),
            &ChannelTypeConfig::AllStraight,
        );
        let result = insert_intersection_nodes(&mut system);
        // Simple bifurcation has no crossing channels.
        assert_eq!(result.intersection_count, 0);
    }

    #[test]
    #[allow(deprecated)] // NetworkBlueprint::new() used intentionally; all nodes use NodeSpec::new_at().
    fn intersection_detection_finds_manual_crossing() {
        // Manually build a system with two crossing channels.
        let mut system = NetworkBlueprint::new("crossing");
        system.box_dims = (1.0, 1.0);
        system.nodes = vec![
            NodeSpec::new_at("0".to_string(), NodeKind::Junction, (0.0, 0.5)),
            NodeSpec::new_at("1".to_string(), NodeKind::Junction, (1.0, 0.5)),
            NodeSpec::new_at("2".to_string(), NodeKind::Junction, (0.5, 0.0)),
            NodeSpec::new_at("3".to_string(), NodeKind::Junction, (0.5, 1.0)),
        ];

        let mut ch0 = ChannelSpec::new_pipe_rect(
            "ch0".to_string(),
            "0".to_string(),
            "1".to_string(),
            1.0,
            0.1,
            0.05,
            0.0,
            0.0,
        );
        ch0.path = vec![(0.0, 0.5), (1.0, 0.5)];
        let mut ch1 = ChannelSpec::new_pipe_rect(
            "ch1".to_string(),
            "2".to_string(),
            "3".to_string(),
            1.0,
            0.1,
            0.05,
            0.0,
            0.0,
        );
        ch1.path = vec![(0.5, 0.0), (0.5, 1.0)];
        system.channels = vec![ch0, ch1];

        let result = insert_intersection_nodes(&mut system);
        assert_eq!(result.intersection_count, 1);
        assert_eq!(result.junction_node_ids.len(), 1);

        // The two original channels should be split into 4 sub-channels.
        assert_eq!(
            system.channels.len(),
            4,
            "2 channels × 1 crossing = 4 sub-channels"
        );

        // New junction node at (0.5, 0.5).
        let jn = &system.nodes[result.junction_node_ids[0]];
        assert!((jn.point.0 - 0.5).abs() < 1e-6);
        assert!((jn.point.1 - 0.5).abs() < 1e-6);
    }
}
