//! Low-level face-soup boolean operations

use crate::domain::core::error::{MeshError, MeshResult};
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;
use super::containment::{containment, Containment};

/// CSG boolean operation type.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BooleanOp {
    /// A ∪ B
    Union,
    /// A ∩ B
    Intersection,
    /// A \ B
    Difference,
}

/// Perform a CSG boolean operation on two sets of faces.
///
/// Both sets share the same `VertexPool`; new intersection vertices are
/// inserted into the pool during splitting.
pub fn csg_boolean(
    op: BooleanOp,
    faces_a: &[FaceData],
    faces_b: &[FaceData],
    pool: &mut VertexPool,
) -> MeshResult<Vec<FaceData>> {
    // ── Coplanar flat-mesh pre-check ─────────────────────────────────────────
    // Flat (zero-volume) meshes like disks must be handled before the 3-D
    // containment/ray-cast logic, which is unreliable on zero-volume geometry.
    // If both operands lie in the same plane, dispatch directly to the 2-D
    // Sutherland-Hodgman pipeline regardless of AABB containment.
    if let (Some(basis), Some(_)) = (
        crate::application::csg::coplanar::detect_flat_plane(faces_a, pool),
        crate::application::csg::coplanar::detect_flat_plane(faces_b, pool),
    ) {
        let result = crate::application::csg::coplanar::boolean_coplanar(op, faces_a, faces_b, pool, &basis);
        if result.is_empty() {
            return Err(MeshError::EmptyBooleanResult {
                op: format!("{:?}", op),
            });
        }
        return Ok(result);
    }

    let result = match op {
        BooleanOp::Union => boolean_union(faces_a, faces_b, pool),
        BooleanOp::Intersection => boolean_intersection(faces_a, faces_b, pool),
        BooleanOp::Difference => boolean_difference(faces_a, faces_b, pool),
    };

    if result.is_empty() {
        return Err(MeshError::EmptyBooleanResult {
            op: format!("{:?}", op),
        });
    }

    Ok(result)
}

/// Reverse the winding of every face (flips the surface normal direction).
fn flip_faces(faces: &[FaceData]) -> Vec<FaceData> {
    faces
        .iter()
        .map(|f| {
            let mut ff = *f;
            ff.flip();
            ff
        })
        .collect()
}

/// A ∪ B — union of two face soups.
fn boolean_union(
    faces_a: &[FaceData],
    faces_b: &[FaceData],
    pool: &mut VertexPool,
) -> Vec<FaceData> {
    match containment(faces_a, faces_b, pool) {
        Containment::BInsideA => faces_a.to_vec(),
        Containment::AInsideB => faces_b.to_vec(),
        Containment::Disjoint => {
            let mut r = faces_a.to_vec();
            r.extend_from_slice(faces_b);
            r
        }
        Containment::Intersecting => {
            // ── Mesh Arrangement pipeline (flat or curved — always watertight).
            // Coplanar case is already handled by the outer csg_boolean guard.
            crate::application::csg::arrangement::boolean_intersecting_arrangement(
                BooleanOp::Union,
                faces_a,
                faces_b,
                pool,
            )
        }
    }
}

/// A ∩ B — intersection of two face soups.
fn boolean_intersection(
    faces_a: &[FaceData],
    faces_b: &[FaceData],
    pool: &mut VertexPool,
) -> Vec<FaceData> {
    match containment(faces_a, faces_b, pool) {
        Containment::BInsideA => faces_b.to_vec(),
        Containment::AInsideB => faces_a.to_vec(),
        Containment::Disjoint => Vec::new(),
        Containment::Intersecting => {
            // ── Mesh Arrangement pipeline (flat or curved — always watertight).
            // Coplanar case is already handled by the outer csg_boolean guard.
            crate::application::csg::arrangement::boolean_intersecting_arrangement(
                BooleanOp::Intersection,
                faces_a,
                faces_b,
                pool,
            )
        }
    }
}

/// A \ B — subtract B from A.
fn boolean_difference(
    faces_a: &[FaceData],
    faces_b: &[FaceData],
    pool: &mut VertexPool,
) -> Vec<FaceData> {
    match containment(faces_a, faces_b, pool) {
        Containment::BInsideA => {
            // B is a cavity inside A: keep A's outer surface + B's surface flipped inward.
            let mut r = faces_a.to_vec();
            r.extend(flip_faces(faces_b));
            r
        }
        Containment::AInsideB => Vec::new(), // A is swallowed by B → empty
        Containment::Disjoint => faces_a.to_vec(), // B doesn't touch A → A unchanged
        Containment::Intersecting => {
            // ── Mesh Arrangement pipeline (flat or curved — always watertight).
            // Coplanar case is already handled by the outer csg_boolean guard.
            crate::application::csg::arrangement::boolean_intersecting_arrangement(
                BooleanOp::Difference,
                faces_a,
                faces_b,
                pool,
            )
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::core::scalar::{Point3r, Vector3r};
    use super::super::containment::point_in_mesh;

    /// Build a closed 2×2×2 axis-aligned cube centred at (1,1,1) with outward CCW faces.
    fn make_unit_cube(pool: &mut VertexPool) -> Vec<FaceData> {
        let mut p =
            |x: f64, y: f64, z: f64| pool.insert_or_weld(Point3r::new(x, y, z), Vector3r::zeros());
        let (p000, p100, p110, p010) = (p(0., 0., 0.), p(2., 0., 0.), p(2., 2., 0.), p(0., 2., 0.));
        let (p001, p101, p111, p011) = (p(0., 0., 2.), p(2., 0., 2.), p(2., 2., 2.), p(0., 2., 2.));
        vec![
            // -Z
            FaceData::untagged(p000, p010, p110),
            FaceData::untagged(p000, p110, p100),
            // +Z
            FaceData::untagged(p001, p101, p111),
            FaceData::untagged(p001, p111, p011),
            // -Y
            FaceData::untagged(p000, p100, p101),
            FaceData::untagged(p000, p101, p001),
            // +Y
            FaceData::untagged(p010, p011, p111),
            FaceData::untagged(p010, p111, p110),
            // -X
            FaceData::untagged(p000, p001, p011),
            FaceData::untagged(p000, p011, p010),
            // +X
            FaceData::untagged(p100, p110, p111),
            FaceData::untagged(p100, p111, p101),
        ]
    }

    /// Build a 0.5×0.5×0.5 cube at offset (o,o,o).
    fn make_inner_cube(pool: &mut VertexPool, o: f64) -> Vec<FaceData> {
        let s = 0.5_f64;
        let mut p =
            |x: f64, y: f64, z: f64| pool.insert_or_weld(Point3r::new(x, y, z), Vector3r::zeros());
        let (p000, p100, p110, p010) = (
            p(o, o, o),
            p(o + s, o, o),
            p(o + s, o + s, o),
            p(o, o + s, o),
        );
        let (p001, p101, p111, p011) = (
            p(o, o, o + s),
            p(o + s, o, o + s),
            p(o + s, o + s, o + s),
            p(o, o + s, o + s),
        );
        vec![
            FaceData::untagged(p000, p010, p110),
            FaceData::untagged(p000, p110, p100),
            FaceData::untagged(p001, p101, p111),
            FaceData::untagged(p001, p111, p011),
            FaceData::untagged(p000, p100, p101),
            FaceData::untagged(p000, p101, p001),
            FaceData::untagged(p010, p011, p111),
            FaceData::untagged(p010, p111, p110),
            FaceData::untagged(p000, p001, p011),
            FaceData::untagged(p000, p011, p010),
            FaceData::untagged(p100, p110, p111),
            FaceData::untagged(p100, p111, p101),
        ]
    }

    #[test]
    fn point_in_cube_inside() {
        let mut pool = VertexPool::default_millifluidic();
        let faces = make_unit_cube(&mut pool);
        assert!(point_in_mesh(&Point3r::new(1.0, 1.0, 1.0), &faces, &pool));
    }

    #[test]
    fn point_in_cube_outside() {
        let mut pool = VertexPool::default_millifluidic();
        let faces = make_unit_cube(&mut pool);
        assert!(!point_in_mesh(&Point3r::new(5.0, 5.0, 5.0), &faces, &pool));
    }

    #[test]
    fn containment_b_inside_a_detected() {
        let mut pool = VertexPool::default_millifluidic();
        let faces_a = make_unit_cube(&mut pool);
        let faces_b = make_inner_cube(&mut pool, 0.75); // centred inside big cube
        assert_eq!(
            containment(&faces_a, &faces_b, &pool),
            Containment::BInsideA
        );
    }

    #[test]
    fn containment_disjoint_detected() {
        let mut pool = VertexPool::default_millifluidic();
        let faces_a = make_unit_cube(&mut pool);
        let mut p =
            |x: f64, y: f64, z: f64| pool.insert_or_weld(Point3r::new(x, y, z), Vector3r::zeros());
         // Tiny cube far away
        let (p0, p1, p2, p3) = (
             p(10., 10., 10.),
             p(11., 10., 10.),
             p(11., 11., 10.),
             p(10., 11., 10.),
        );
        let (p4, p5, p6, p7) = (
            p(10., 10., 11.),
            p(11., 10., 11.),
            p(11., 11., 11.),
            p(10., 11., 11.),
        );
        let faces_b = vec![
            FaceData::untagged(p0, p2, p1),
            FaceData::untagged(p0, p3, p2),
            FaceData::untagged(p4, p5, p6),
            FaceData::untagged(p4, p6, p7),
            FaceData::untagged(p0, p1, p5),
            FaceData::untagged(p0, p5, p4),
            FaceData::untagged(p3, p7, p6),
            FaceData::untagged(p3, p6, p2),
            FaceData::untagged(p0, p4, p7),
            FaceData::untagged(p0, p7, p3),
            FaceData::untagged(p1, p2, p6),
            FaceData::untagged(p1, p6, p5),
        ];
        assert_eq!(
            containment(&faces_a, &faces_b, &pool),
            Containment::Disjoint
        );
    }

    #[test]
    fn difference_enclosed_produces_outer_plus_flipped_inner() {
        let mut pool = VertexPool::default_millifluidic();
        let faces_a = make_unit_cube(&mut pool);
        let faces_b = make_inner_cube(&mut pool, 0.75);
        let result = boolean_difference(&faces_a, &faces_b, &mut pool);
        // Outer (12) + inner flipped (12) = 24
        assert_eq!(result.len(), 24);
    }

    #[test]
    fn union_enclosed_returns_outer() {
        let mut pool = VertexPool::default_millifluidic();
        let faces_a = make_unit_cube(&mut pool);
        let faces_b = make_inner_cube(&mut pool, 0.75);
        let result = boolean_union(&faces_a, &faces_b, &mut pool);
        assert_eq!(result.len(), faces_a.len());
    }

    #[test]
    fn intersection_enclosed_returns_inner() {
        let mut pool = VertexPool::default_millifluidic();
        let faces_a = make_unit_cube(&mut pool);
        let faces_b = make_inner_cube(&mut pool, 0.75);
        let result = boolean_intersection(&faces_a, &faces_b, &mut pool);
        assert_eq!(result.len(), faces_b.len());
    }

    #[test]
    fn difference_disjoint_returns_a() {
        let mut pool = VertexPool::default_millifluidic();
        let faces_a = make_unit_cube(&mut pool);
        let mut p =
            |x: f64, y: f64, z: f64| pool.insert_or_weld(Point3r::new(x, y, z), Vector3r::zeros());
        let (p0, p1, p2, p3) = (
            p(10., 10., 10.),
            p(11., 10., 10.),
            p(11., 11., 10.),
            p(10., 11., 10.),
        );
        let (p4, p5, p6, p7) = (
            p(10., 10., 11.),
            p(11., 10., 11.),
            p(11., 11., 11.),
            p(10., 11., 11.),
        );
        let faces_b = vec![
            FaceData::untagged(p0, p2, p1),
            FaceData::untagged(p0, p3, p2),
            FaceData::untagged(p4, p5, p6),
            FaceData::untagged(p4, p6, p7),
            FaceData::untagged(p0, p1, p5),
            FaceData::untagged(p0, p5, p4),
            FaceData::untagged(p3, p7, p6),
            FaceData::untagged(p3, p6, p2),
            FaceData::untagged(p0, p4, p7),
            FaceData::untagged(p0, p7, p3),
            FaceData::untagged(p1, p2, p6),
            FaceData::untagged(p1, p6, p5),
        ];
        let result = boolean_difference(&faces_a, &faces_b, &mut pool);
        assert_eq!(result.len(), faces_a.len());
    }
}
