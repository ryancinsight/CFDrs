//! Polygon classification against a splitting plane.

use crate::core::index::VertexId;
use crate::core::scalar::{Real, Point3r, TOLERANCE};
use crate::geometry::plane::{Plane, PointClassification};
use crate::storage::vertex_pool::VertexPool;

/// Classification of a polygon relative to a plane.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PolygonClassification {
    /// All vertices are in front of (or coplanar with) the plane.
    Front,
    /// All vertices are behind (or coplanar with) the plane.
    Back,
    /// All vertices are coplanar with the plane.
    Coplanar,
    /// Polygon straddles the plane.
    Spanning,
}

/// Classify a triangle relative to a plane.
///
/// Returns the classification and the individual vertex classifications.
pub fn classify_triangle(
    verts: [VertexId; 3],
    pool: &VertexPool,
    plane: &Plane,
) -> (PolygonClassification, [PointClassification; 3]) {
    let classifications = [
        plane.classify_point(&pool.position(verts[0])),
        plane.classify_point(&pool.position(verts[1])),
        plane.classify_point(&pool.position(verts[2])),
    ];

    let mut front = false;
    let mut back = false;
    for &c in &classifications {
        match c {
            PointClassification::Front => front = true,
            PointClassification::Back => back = true,
            PointClassification::Coplanar => {}
        }
    }

    let poly_class = match (front, back) {
        (true, true) => PolygonClassification::Spanning,
        (true, false) => PolygonClassification::Front,
        (false, true) => PolygonClassification::Back,
        (false, false) => PolygonClassification::Coplanar,
    };

    (poly_class, classifications)
}
