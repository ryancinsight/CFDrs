//! 2-D plane basis projection and exact flat-plane detection.
//!
//! ## Algorithm — Exact Coplanar Plane Detection
//!
//! 1. Pick a representative non-degenerate triangle `(a,b,c)` from `faces`.
//! 2. Build a projection basis from `(a,b,c)` (`PlaneBasis::from_triangle`).
//! 3. For every vertex `p` in every face, evaluate `orient3d(a,b,c,p)`.
//! 4. Return `Some(basis)` iff all orientations are exactly `Sign::Zero`.
//!
//! ## Theorem — Coplanarity Equivalence
//!
//! Let `(a,b,c)` be a non-collinear triangle. A point set `P` is coplanar with
//! the plane of `(a,b,c)` iff `orient3d(a,b,c,p) == 0` for every `p ∈ P`.
//!
//! **Proof sketch.**
//! The signed tetrahedral volume determinant `orient3d(a,b,c,p)` is zero
//! exactly when `p` lies in the affine span of `(a,b,c)`. Because `(a,b,c)` is
//! non-collinear, that span is a unique plane. Therefore all-zero determinants
//! are equivalent to global coplanarity. ∎

use crate::domain::core::scalar::{Point3r, Real, Vector3r};
use crate::domain::geometry::predicates::{orient_2d_arr, Orientation};
use crate::domain::topology::predicates::{orient3d, Sign};
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;

pub(crate) struct PlaneBasis {
    pub(crate) origin: Point3r,
    pub(crate) u: Vector3r,
    pub(crate) v: Vector3r,
    pub(crate) normal: Vector3r,
}

impl PlaneBasis {
    pub(crate) fn from_triangle(a: &Point3r, b: &Point3r, c: &Point3r) -> Option<Self> {
        let ab = b - a;
        let ac = c - a;
        let n = ab.cross(&ac);
        let nl = n.norm();
        if nl < 1e-20 {
            return None;
        }
        let ul = ab.norm();
        if ul < 1e-20 {
            return None;
        }
        let u = ab / ul;
        let normal = n / nl;
        let v = normal.cross(&u).normalize();
        Some(Self {
            origin: *a,
            u,
            v,
            normal,
        })
    }

    #[inline]
    pub(crate) fn project(&self, p: &Point3r) -> [Real; 2] {
        let d = p - self.origin;
        [d.dot(&self.u), d.dot(&self.v)]
    }

    /// Lift a 2-D point (u,v) back to 3-D.
    #[inline]
    pub(crate) fn lift(&self, u: Real, v: Real) -> Point3r {
        self.origin + self.u * u + self.v * v
    }
}

#[inline]
fn triangle_is_degenerate_exact(a: &Point3r, b: &Point3r, c: &Point3r) -> bool {
    orient_2d_arr([a.x, a.y], [b.x, b.y], [c.x, c.y]) == Orientation::Degenerate
        && orient_2d_arr([a.x, a.z], [b.x, b.z], [c.x, c.z]) == Orientation::Degenerate
        && orient_2d_arr([a.y, a.z], [b.y, b.z], [c.y, c.z]) == Orientation::Degenerate
}

pub(crate) fn detect_flat_plane(faces: &[FaceData], pool: &VertexPool) -> Option<PlaneBasis> {
    let mut basis: Option<PlaneBasis> = None;
    let mut rep_tri: Option<(Point3r, Point3r, Point3r)> = None;
    for face in faces {
        let a = *pool.position(face.vertices[0]);
        let b = *pool.position(face.vertices[1]);
        let c = *pool.position(face.vertices[2]);
        if triangle_is_degenerate_exact(&a, &b, &c) {
            continue;
        }
        if let Some(b0) = PlaneBasis::from_triangle(&a, &b, &c) {
            basis = Some(b0);
            rep_tri = Some((a, b, c));
            break;
        }
    }
    let basis = basis?;
    let (ra, rb, rc) = rep_tri?;
    for face in faces {
        for &vid in &face.vertices {
            if orient3d(&ra, &rb, &rc, pool.position(vid)) != Sign::Zero {
                return None;
            }
        }
    }
    Some(basis)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::infrastructure::storage::face_store::FaceData;

    fn p(x: Real, y: Real, z: Real) -> Point3r {
        Point3r::new(x, y, z)
    }

    #[test]
    fn detect_flat_plane_accepts_exactly_coplanar_faces() {
        let mut pool = VertexPool::for_csg();
        let n = Vector3r::new(0.0, 0.0, 1.0);

        let v0 = pool.insert_unique(p(0.0, 0.0, 0.0), n);
        let v1 = pool.insert_unique(p(1.0, 0.0, 0.0), n);
        let v2 = pool.insert_unique(p(1.0, 1.0, 0.0), n);
        let v3 = pool.insert_unique(p(0.0, 1.0, 0.0), n);

        let faces = vec![
            FaceData::untagged(v0, v1, v2),
            FaceData::untagged(v0, v2, v3),
        ];
        assert!(detect_flat_plane(&faces, &pool).is_some());
    }

    #[test]
    fn detect_flat_plane_rejects_algebraically_non_coplanar_vertex() {
        let mut pool = VertexPool::for_csg();
        let n = Vector3r::new(0.0, 0.0, 1.0);

        let v0 = pool.insert_unique(p(0.0, 0.0, 0.0), n);
        let v1 = pool.insert_unique(p(1.0, 0.0, 0.0), n);
        let v2 = pool.insert_unique(p(0.0, 1.0, 0.0), n);
        let v3 = pool.insert_unique(p(0.25, 0.25, 1.0e-12), n);

        let faces = vec![
            FaceData::untagged(v0, v1, v2),
            FaceData::untagged(v0, v2, v3),
        ];
        assert!(detect_flat_plane(&faces, &pool).is_none());
    }

    #[test]
    fn detect_flat_plane_ignores_degenerate_seed_face() {
        let mut pool = VertexPool::for_csg();
        let n = Vector3r::new(0.0, 0.0, 1.0);

        let v0 = pool.insert_unique(p(0.0, 0.0, 0.0), n);
        let v1 = pool.insert_unique(p(1.0, 0.0, 0.0), n);
        let v2 = pool.insert_unique(p(2.0, 0.0, 0.0), n); // collinear with v0,v1
        let v3 = pool.insert_unique(p(0.0, 1.0, 0.0), n);

        let faces = vec![
            FaceData::untagged(v0, v1, v2), // degenerate
            FaceData::untagged(v0, v1, v3), // valid plane seed
        ];
        assert!(detect_flat_plane(&faces, &pool).is_some());
    }
}
