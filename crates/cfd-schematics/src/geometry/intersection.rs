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

use super::types::{centerline_for_channel, Channel, ChannelSystem, ChannelType, Node, Point2D};
use crate::geometry::metadata::MetadataContainer;

/// Metadata marker for nodes created at channel intersections.
#[derive(Debug, Clone)]
pub struct IntersectionMetadata {
    /// IDs of the two channels that cross at this node.
    pub channel_a_id: usize,
    /// ID of the second crossing channel.
    pub channel_b_id: usize,
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

/// Detect all pairwise intersections between channel centerlines.
fn detect_crossings(system: &ChannelSystem) -> Vec<SegmentCrossing> {
    let centerlines: Vec<Vec<Point2D>> = system
        .channels
        .iter()
        .map(|ch| centerline_for_channel(ch, &system.nodes))
        .collect();

    let mut crossings = Vec::new();

    for (i, cl_a) in centerlines.iter().enumerate() {
        for (j, cl_b) in centerlines.iter().enumerate() {
            if j <= i {
                continue; // avoid duplicates
            }
            // Skip channels that share an endpoint node — those are
            // intentional junctions, not geometric crossings.
            let ch_a = &system.channels[i];
            let ch_b = &system.channels[j];
            if ch_a.from_node == ch_b.from_node
                || ch_a.from_node == ch_b.to_node
                || ch_a.to_node == ch_b.from_node
                || ch_a.to_node == ch_b.to_node
            {
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
pub fn insert_intersection_nodes(system: &mut ChannelSystem) -> IntersectionResult {
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
    for crossing in &crossings {
        let node_id = system.nodes.len();
        let mut meta_container = MetadataContainer::new();
        meta_container.insert(IntersectionMetadata {
            channel_a_id: system.channels[crossing.channel_a].id,
            channel_b_id: system.channels[crossing.channel_b].id,
        });
        system.nodes.push(Node {
            id: node_id,
            point: crossing.point,
            metadata: Some(meta_container),
        });
        junction_ids.push(node_id);
    }

    // Build a map: original_channel_index → Vec<(crossing_index, role)>
    // where role indicates whether this channel is A or B in the crossing.
    struct SplitInfo {
        junction_node_id: usize,
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
                junction_node_id: junction_ids[ci],
                t: crossing.t_a,
                seg: crossing.seg_a,
            });
        splits_per_channel
            .entry(crossing.channel_b)
            .or_default()
            .push(SplitInfo {
                junction_node_id: junction_ids[ci],
                t: crossing.t_b,
                seg: crossing.seg_b,
            });
    }

    // Sort each channel's splits by (segment, t) so we can process them in order.
    for (_ch, splits) in splits_per_channel.iter_mut() {
        splits.sort_by(|a, b| a.seg.cmp(&b.seg).then(a.t.partial_cmp(&b.t).unwrap()));
    }

    // Rebuild the channel list.
    let mut next_channel_id = system.channels.iter().map(|c| c.id).max().unwrap_or(0) + 1;
    let mut new_channels: Vec<Channel> = Vec::new();

    for (ch_idx, orig_channel) in system.channels.iter().enumerate() {
        if let Some(splits) = splits_per_channel.get(&ch_idx) {
            // This channel is split at one or more crossings.
            // Walk through the centerline and produce sub-channels.
            let centerline = centerline_for_channel(orig_channel, &system.nodes);
            if centerline.len() < 2 {
                new_channels.push(orig_channel.clone());
                continue;
            }

            // Build a sequence of (node_id, point) for the split points,
            // ordered along the channel path.
            let mut waypoints: Vec<(usize, Point2D)> = Vec::new();
            waypoints.push((orig_channel.from_node, centerline[0]));

            // Insert splits in order along the centerline.
            let mut path_idx = 0; // current segment index
            let mut split_iter = splits.iter().peekable();

            while path_idx < centerline.len() - 1 {
                // Collect all splits on this segment, sorted by t.
                while let Some(split) = split_iter.peek() {
                    if split.seg == path_idx {
                        let s = split_iter.next().unwrap();
                        let pt = (
                            centerline[path_idx]
                                .0
                                .mul_add(1.0 - s.t, centerline[path_idx + 1].0 * s.t),
                            centerline[path_idx]
                                .1
                                .mul_add(1.0 - s.t, centerline[path_idx + 1].1 * s.t),
                        );
                        waypoints.push((s.junction_node_id, pt));
                    } else {
                        break;
                    }
                }

                // Also include internal centerline points as path geometry
                // (they are not nodes, just path fidelity).
                path_idx += 1;
            }

            waypoints.push((orig_channel.to_node, *centerline.last().unwrap()));

            // Now produce sub-channels between consecutive waypoints.
            for pair in waypoints.windows(2) {
                let (from_id, _from_pt) = pair[0];
                let (to_id, _to_pt) = pair[1];
                new_channels.push(Channel {
                    id: next_channel_id,
                    from_node: from_id,
                    to_node: to_id,
                    width: orig_channel.width,
                    height: orig_channel.height,
                    channel_type: ChannelType::Straight,
                    metadata: orig_channel.metadata.clone(),
                });
                next_channel_id += 1;
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
    use crate::geometry::types::SplitType;

    #[test]
    fn segment_intersection_detects_crossing() {
        let result = segment_intersection((0.0, 0.0), (1.0, 1.0), (0.0, 1.0), (1.0, 0.0));
        assert!(result.is_some());
        let (t, u, pt) = result.unwrap();
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
        let (w8, h8) = adaptive_box_dims(45.0, 45.0, 8, 2.0, 2.0);

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
    fn intersection_detection_finds_manual_crossing() {
        // Manually build a system with two crossing channels.
        let nodes = vec![
            Node {
                id: 0,
                point: (0.0, 0.5),
                metadata: None,
            },
            Node {
                id: 1,
                point: (1.0, 0.5),
                metadata: None,
            },
            Node {
                id: 2,
                point: (0.5, 0.0),
                metadata: None,
            },
            Node {
                id: 3,
                point: (0.5, 1.0),
                metadata: None,
            },
        ];
        let channels = vec![
            Channel {
                id: 0,
                from_node: 0,
                to_node: 1,
                width: 0.1,
                height: 0.05,
                channel_type: ChannelType::Straight,
                metadata: None,
            },
            Channel {
                id: 1,
                from_node: 2,
                to_node: 3,
                width: 0.1,
                height: 0.05,
                channel_type: ChannelType::Straight,
                metadata: None,
            },
        ];
        let mut system = ChannelSystem {
            box_dims: (1.0, 1.0),
            nodes,
            channels,
            box_outline: vec![],
        };

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
