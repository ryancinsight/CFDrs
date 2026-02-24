//! Types for triangle-triangle intersection

use crate::domain::core::scalar::Point3r;

/// A 3-D intersection segment to be snapped onto a face during CDT co-refinement.
#[derive(Debug, Clone, Copy)]
pub struct SnapSegment {
    /// One endpoint of the intersection segment.
    pub start: Point3r,
    /// The other endpoint of the intersection segment.
    pub end: Point3r,
}

/// The result of an exact triangle–triangle intersection test.
#[derive(Debug, Clone)]
pub enum IntersectionType {
    /// Triangles are disjoint or only touch at a single point / edge.
    ///
    /// No Boolean splitting is required for this pair.
    None,

    /// Triangles are coplanar (all vertices of one lie on the plane of the other).
    ///
    /// Coplanar overlap is handled separately via [`crate::application::csg::clip`].
    Coplanar,

    /// Triangles properly straddle each other and intersect along a segment.
    ///
    /// Both triangles must be split at this segment before Boolean
    /// classification.
    Segment {
        /// One endpoint of the intersection segment.
        start: Point3r,
        /// The other endpoint of the intersection segment.
        end: Point3r,
    },
}
