//! UV-parametric sphere primitive.

use std::f64::consts::{PI, TAU};

use crate::core::index::RegionId;
use crate::core::scalar::{Point3r, Vector3r};
use crate::mesh::IndexedMesh;
use super::{PrimitiveMesh, PrimitiveError};

/// Builds a UV-parametric sphere.
///
/// ## Parametrisation
///
/// ```text
/// x(θ, φ) = cx + r · sin(φ) · cos(θ)
/// y(θ, φ) = cy + r · cos(φ)
/// z(θ, φ) = cz + r · sin(φ) · sin(θ)
///
/// θ ∈ [0, 2π)   — longitude
/// φ ∈ [0, π]    — colatitude from +Y
/// ```
///
/// Pole rows collapse to triangles; interior rows are quads → 2 triangles.
///
/// ## Output
///
/// - `segments × (stacks − 2) × 2 + 2 × segments` faces
/// - `signed_volume ≈ 4π r³ / 3`  (error < 1% for segments ≥ 32)
#[derive(Clone, Debug)]
pub struct UvSphere {
    /// Sphere radius [mm].
    pub radius: f64,
    /// Centre position.
    pub center: Point3r,
    /// Number of longitude subdivisions (≥ 3).
    pub segments: usize,
    /// Number of latitude subdivisions (≥ 2).
    pub stacks: usize,
}

impl Default for UvSphere {
    fn default() -> Self {
        Self {
            radius: 1.0,
            center: Point3r::origin(),
            segments: 32,
            stacks: 16,
        }
    }
}

impl PrimitiveMesh for UvSphere {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build(self)
    }
}

fn build(s: &UvSphere) -> Result<IndexedMesh, PrimitiveError> {
    if s.radius <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "radius must be > 0, got {}", s.radius
        )));
    }
    if s.segments < 3 {
        return Err(PrimitiveError::TooFewSegments(s.segments));
    }
    if s.stacks < 2 {
        return Err(PrimitiveError::InvalidParam(format!(
            "stacks must be ≥ 2, got {}", s.stacks
        )));
    }

    let region = RegionId::new(1);
    let mut mesh = IndexedMesh::new();
    let r = s.radius;
    let cx = s.center.x;
    let cy = s.center.y;
    let cz = s.center.z;

    // Sample a (position, outward_normal) pair at (θ, φ).
    let vertex_at = |theta: f64, phi: f64| -> (Point3r, Vector3r) {
        let sp = phi.sin();
        let cp = phi.cos();
        let ct = theta.cos();
        let st = theta.sin();
        let n = Vector3r::new(sp * ct, cp, sp * st);
        let pos = Point3r::new(cx + r * n.x, cy + r * n.y, cz + r * n.z);
        (pos, n)
    };

    for i in 0..s.segments {
        let theta0 = i as f64 / s.segments as f64 * TAU;
        let theta1 = (i + 1) as f64 / s.segments as f64 * TAU;

        for j in 0..s.stacks {
            let phi0 = j as f64 / s.stacks as f64 * PI;
            let phi1 = (j + 1) as f64 / s.stacks as f64 * PI;

            let (pos00, n00) = vertex_at(theta0, phi0);
            let (pos10, n10) = vertex_at(theta1, phi0);
            let (pos11, n11) = vertex_at(theta1, phi1);
            let (pos01, n01) = vertex_at(theta0, phi1);

            let v00 = mesh.add_vertex(pos00, n00);
            let v10 = mesh.add_vertex(pos10, n10);
            let v11 = mesh.add_vertex(pos11, n11);
            let v01 = mesh.add_vertex(pos01, n01);

            if j == 0 {
                // North-pole row: collapse top two vertices → triangle
                // Outward CCW from outside: v10 → v11 → v01
                mesh.add_face_with_region(v10, v11, v01, region);
            } else if j == s.stacks - 1 {
                // South-pole row: collapse bottom two vertices → triangle
                // Outward CCW from outside: v00 → v10 → v01
                mesh.add_face_with_region(v00, v10, v01, region);
            } else {
                // Interior quad → 2 outward CCW triangles
                mesh.add_face_with_region(v00, v10, v11, region);
                mesh.add_face_with_region(v00, v11, v01, region);
            }
        }
    }

    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::storage::edge_store::EdgeStore;
    use crate::watertight::check::check_watertight;
    use std::f64::consts::PI;

    #[test]
    fn sphere_is_watertight() {
        let mesh = UvSphere::default().build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.is_watertight, "sphere must be watertight");
        assert_eq!(report.euler_characteristic, Some(2));
    }

    #[test]
    fn sphere_volume_approximately_correct() {
        let r = 1.0_f64;
        let mesh = UvSphere { radius: r, segments: 64, stacks: 32, ..UvSphere::default() }
            .build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        let expected = 4.0 * PI / 3.0 * r * r * r;
        let error = (report.signed_volume - expected).abs() / expected;
        assert!(
            error < 0.005,
            "volume error {:.4}% should be < 0.5% with 64×32 mesh",
            error * 100.0
        );
    }

    #[test]
    fn sphere_volume_positive() {
        let mesh = UvSphere::default().build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.signed_volume > 0.0, "outward normals → positive volume");
    }

    #[test]
    fn sphere_too_few_segments() {
        let result = UvSphere { segments: 2, ..UvSphere::default() }.build();
        assert!(result.is_err());
    }
}
