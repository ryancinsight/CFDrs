//! Triangle splitting against a plane.
//!
//! When a triangle straddles a plane, it is split into up to 3 sub-triangles,
//! some in front and some behind.

use crate::core::index::{VertexId, RegionId};
use crate::core::scalar::{Point3r, Vector3r};
use crate::geometry::plane::{Plane, PointClassification};
use crate::storage::face_store::FaceData;
use crate::storage::vertex_pool::VertexPool;
use crate::csg::bsp::BSP_PLANE_TOLERANCE;

/// Result of splitting a triangle against a plane.
pub struct SplitResult {
    /// Triangles on the front side.
    pub front: Vec<FaceData>,
    /// Triangles on the back side.
    pub back: Vec<FaceData>,
}

/// Split a triangle by a plane.
///
/// New vertices are inserted into `pool` via welding. The result is a set of
/// sub-triangles partitioned into front and back.
pub fn split_triangle(
    face: &FaceData,
    classifications: [PointClassification; 3],
    plane: &Plane,
    pool: &mut VertexPool,
) -> SplitResult {
    let mut front_verts: Vec<VertexId> = Vec::with_capacity(4);
    let mut back_verts: Vec<VertexId> = Vec::with_capacity(4);

    for i in 0..3 {
        let j = (i + 1) % 3;
        let vi = face.vertices[i];
        let vj = face.vertices[j];
        let ci = classifications[i];
        let cj = classifications[j];

        match ci {
            PointClassification::Front => {
                front_verts.push(vi);
            }
            PointClassification::Back => {
                back_verts.push(vi);
            }
            PointClassification::Coplanar => {
                front_verts.push(vi);
                back_verts.push(vi);
            }
        }

        // If the edge crosses the plane, compute the intersection point
        let crosses = matches!(
            (ci, cj),
            (PointClassification::Front, PointClassification::Back)
                | (PointClassification::Back, PointClassification::Front)
        );

        if crosses {
            let pi = pool.position(vi);
            let pj = pool.position(vj);
            if let Some(t) = plane.intersect_segment_with_eps(&pi, &pj, BSP_PLANE_TOLERANCE) {
                let intersection = Point3r::from(pi.coords.lerp(&pj.coords, t));
                // Compute interpolated normal
                let ni = pool.normal(vi);
                let nj = pool.normal(vj);
                let interp_normal = Vector3r::from(ni.lerp(&nj, t));
                let new_vid = pool.insert_or_weld(intersection, interp_normal);
                front_verts.push(new_vid);
                back_verts.push(new_vid);
            }
        }
    }

    let region = face.region;

    SplitResult {
        front: fan_triangulate(&front_verts, region),
        back: fan_triangulate(&back_verts, region),
    }
}

/// Fan-triangulate a convex polygon.
fn fan_triangulate(verts: &[VertexId], region: RegionId) -> Vec<FaceData> {
    if verts.len() < 3 {
        return Vec::new();
    }
    let mut faces = Vec::with_capacity(verts.len() - 2);
    for i in 1..(verts.len() - 1) {
        faces.push(FaceData {
            vertices: [verts[0], verts[i], verts[i + 1]],
            region,
        });
    }
    faces
}
