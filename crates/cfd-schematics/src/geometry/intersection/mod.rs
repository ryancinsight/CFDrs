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
mod tests;
