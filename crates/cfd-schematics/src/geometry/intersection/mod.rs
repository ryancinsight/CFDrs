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

mod adaptive;
mod detection;
mod insertion;

pub use adaptive::adaptive_box_dims;
pub use detection::{has_unresolved_intersections, unresolved_intersection_count};
pub use insertion::insert_intersection_nodes;

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

#[cfg(test)]
mod tests {
    use super::detection::segment_intersection;
    use super::*;
    use crate::config::{ChannelTypeConfig, GeometryConfig};
    use crate::domain::model::{ChannelSpec, NetworkBlueprint, NodeKind, NodeSpec};
    use crate::geometry::generator::create_geometry;
    use crate::geometry::SplitType;

    #[test]
    fn segment_intersection_detects_crossing() {
        let result = segment_intersection((0.0, 0.0), (1.0, 1.0), (0.0, 1.0), (1.0, 0.0));
        assert!(result.is_some());
        let (t, u, point) = result.expect("structural invariant");
        assert!((t - 0.5).abs() < 1e-10);
        assert!((u - 0.5).abs() < 1e-10);
        assert!((point.0 - 0.5).abs() < 1e-10);
        assert!((point.1 - 0.5).abs() < 1e-10);
    }

    #[test]
    fn segment_intersection_detects_no_crossing_for_parallel() {
        let result = segment_intersection((0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0));
        assert!(result.is_none());
    }

    #[test]
    fn segment_intersection_excludes_shared_endpoints() {
        let result = segment_intersection((0.0, 0.0), (1.0, 0.0), (1.0, 0.0), (2.0, 1.0));
        assert!(result.is_none());
    }

    #[test]
    fn adaptive_box_dims_scales_with_branches() {
        let (w1, h1) = adaptive_box_dims(45.0, 45.0, 1, 2.0, 2.0);
        let (w2, h2) = adaptive_box_dims(45.0, 45.0, 4, 2.0, 2.0);
        let (_w8, h8) = adaptive_box_dims(45.0, 45.0, 8, 2.0, 2.0);

        assert!((w1 - 45.0).abs() < 1e-10);
        assert!((w2 - 45.0).abs() < 1e-10);
        assert!(h2 > h1, "4 branches should use more height than 1");
        assert!(h8 > h2, "8 branches should use more height than 4");
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
        assert_eq!(result.intersection_count, 0);
    }

    #[test]
    #[allow(deprecated)]
    fn intersection_detection_finds_manual_crossing() {
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
        assert_eq!(
            system.channels.len(),
            4,
            "2 channels × 1 crossing = 4 sub-channels"
        );

        let junction = &system.nodes[result.junction_node_ids[0]];
        assert!((junction.point.0 - 0.5).abs() < 1e-6);
        assert!((junction.point.1 - 0.5).abs() < 1e-6);
    }
}
