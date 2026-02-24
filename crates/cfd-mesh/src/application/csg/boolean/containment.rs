//! Containment detection
//!
//! A pure BSP-only approach fails when one mesh is entirely enclosed within the
//! other (no surface intersection exists), because no splitting plane from one
//! mesh ever splits the other's triangles. The functions here detect this case
//! via exact Generalized Winding Number (GWN) containment (`point_in_mesh`) and handle it analytically:

use crate::domain::core::scalar::Point3r;
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;

/// Spatial relationship between two face soups.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Containment {
    /// The meshes' surfaces actually cross (or cannot be determined otherwise).
    Intersecting,
    /// Mesh B is entirely inside mesh A — no surface crossing.
    BInsideA,
    /// Mesh A is entirely inside mesh B — no surface crossing.
    AInsideB,
    /// The meshes are entirely disjoint — no overlap of any kind.
    Disjoint,
}

/// Test whether `query` is inside a closed triangle mesh using the mathematically
/// exact Generalized Winding Number (GWN).
///
/// Returns `true` if the GWN evaluated at `query` is > 0.5. Since the query point
/// selected by the caller (the AABB centre) is strictly interior to its own bounded volume,
/// it will not lie exactly on any face of the opposing mesh, making the GWN result sharp
/// (close to 0.0 or 1.0) without needing arbitrary raycast nudges.
pub(crate) fn point_in_mesh(query: &Point3r, faces: &[FaceData], pool: &VertexPool) -> bool {
    crate::application::csg::arrangement::classify::gwn(query, faces, pool).abs() > 0.5
}

/// Determine the containment relationship between two face soups in `pool`.
///
/// ## Strategy
///
/// 1. **Disjoint**: AABBs do not overlap → return immediately.
/// 2. **Partial overlap**: neither AABB strictly contains the other → the
///    surfaces must cross → return `Intersecting` immediately (no ray cast
///    needed and avoids sample-point boundary ambiguity).
/// 3. **One AABB inside the other**: do a point-in-mesh ray cast using the
///    centre of the inner AABB as the sample point.  The AABB centre of a
///    convex mesh is strictly interior, so it never lies on a face plane of
///    the enclosing mesh.
pub(crate) fn containment(faces_a: &[FaceData], faces_b: &[FaceData], pool: &VertexPool) -> Containment {
    use crate::domain::geometry::aabb::Aabb;

    let aabb_of = |faces: &[FaceData]| -> Option<Aabb> {
        let mut bb = Aabb::empty();
        let mut any = false;
        for f in faces {
            for &vid in &f.vertices {
                bb.expand(pool.position(vid));
                any = true;
            }
        }
        if any {
            Some(bb)
        } else {
            None
        }
    };

    let aabb_a = match aabb_of(faces_a) {
        Some(a) => a,
        None => return Containment::Disjoint,
    };
    let aabb_b = match aabb_of(faces_b) {
        Some(b) => b,
        None => return Containment::Disjoint,
    };

    if !aabb_a.intersects(&aabb_b) {
        return Containment::Disjoint;
    }

    // Check strict AABB containment: if neither AABB is fully inside the other,
    // the two surfaces must geometrically intersect → use the intersecting-mesh pipeline.
    // `aabb_x.contains_point(min)` and `contains_point(max)` of the other covers all 8 corners.
    let b_in_aabb_a = aabb_a.contains_point(&aabb_b.min) && aabb_a.contains_point(&aabb_b.max);
    let a_in_aabb_b = aabb_b.contains_point(&aabb_a.min) && aabb_b.contains_point(&aabb_a.max);

    if !b_in_aabb_a && !a_in_aabb_b {
        // Partial AABB overlap → surfaces intersect → use the intersecting-mesh pipeline.
        return Containment::Intersecting;
    }

    // One AABB is inside the other.  Ray-cast using the inner AABB's centre:
    // it is guaranteed strictly interior to any convex mesh whose AABB contains it,
    // so it will never sit exactly on a face plane of the outer mesh.
    if b_in_aabb_a {
        let b_sample = aabb_b.center();
        if point_in_mesh(&b_sample, faces_a, pool) {
            return Containment::BInsideA;
        }
    }
    if a_in_aabb_b {
        let a_sample = aabb_a.center();
        if point_in_mesh(&a_sample, faces_b, pool) {
            return Containment::AInsideB;
        }
    }

    // AABB containment but ray cast says otherwise (concave mesh, etc.) → intersecting pipeline.
    Containment::Intersecting
}
