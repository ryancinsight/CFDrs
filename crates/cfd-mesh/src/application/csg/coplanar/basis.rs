//! 2-D Plane Basis projection and flat-plane detection.

use crate::domain::core::scalar::{Point3r, Real, Vector3r};
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

pub(crate) fn detect_flat_plane(faces: &[FaceData], pool: &VertexPool) -> Option<PlaneBasis> {
    let mut basis: Option<PlaneBasis> = None;
    for face in faces {
        let a = pool.position(face.vertices[0]);
        let b = pool.position(face.vertices[1]);
        let c = pool.position(face.vertices[2]);
        if let Some(b0) = PlaneBasis::from_triangle(a, b, c) {
            basis = Some(b0);
            break;
        }
    }
    let basis = basis?;
    const TOL: Real = 1e-6;
    for face in faces {
        for &vid in &face.vertices {
            if (pool.position(vid) - basis.origin).dot(&basis.normal).abs() > TOL {
                return None;
            }
        }
    }
    Some(basis)
}
