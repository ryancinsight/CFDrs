//! Triaxial ellipsoid primitive.

use std::f64::consts::{PI, TAU};

use super::{PrimitiveError, PrimitiveMesh};
use crate::domain::core::index::RegionId;
use crate::domain::core::scalar::{Point3r, Vector3r};
use crate::domain::mesh::IndexedMesh;

/// Builds a triaxial ellipsoid.
///
/// ## Parametrisation
///
/// ```text
/// x(θ, φ) = cx + a · sin(φ) · cos(θ)
/// y(θ, φ) = cy + b · cos(φ)
/// z(θ, φ) = cz + c · sin(φ) · sin(θ)
///
/// θ ∈ [0, 2π)   — longitude
/// φ ∈ [0, π]    — colatitude from +Y
/// ```
///
/// The outward unit normal at `(θ, φ)` is the gradient of the implicit
/// equation `(x/a)² + (y/b)² + (z/c)² = 1`, which gives:
///
/// ```text
/// n = normalize( sin(φ)·cos(θ)/a,  cos(φ)/b,  sin(φ)·sin(θ)/c )
/// ```
///
/// This differs from the position-vector normal when `a ≠ b ≠ c`, making
/// the `Ellipsoid` non-interchangeable with a scaled [`UvSphere`].
///
/// ## Output
///
/// - `segments × (stacks − 2) × 2 + 2 × segments` faces
/// - `signed_volume = (4/3)·π·a·b·c` (exact formula)
/// - All faces tagged `RegionId(1)`
///
/// ## Uses
///
/// Droplet deformation in extensional/shear flow, red blood cell oblate/prolate
/// states, Taylor deformation parameter `D = (L−B)/(L+B)` reference shapes,
/// non-spherical SDT fibre inclusions.
#[derive(Clone, Debug)]
pub struct Ellipsoid {
    /// Semi-axis along X [mm].
    pub semi_x: f64,
    /// Semi-axis along Y [mm].
    pub semi_y: f64,
    /// Semi-axis along Z [mm].
    pub semi_z: f64,
    /// Centre position.
    pub center: Point3r,
    /// Number of longitude subdivisions (≥ 3).
    pub segments: usize,
    /// Number of latitude subdivisions (≥ 2).
    pub stacks: usize,
}

impl Default for Ellipsoid {
    fn default() -> Self {
        Self {
            semi_x: 1.0,
            semi_y: 0.5,
            semi_z: 0.75,
            center: Point3r::origin(),
            segments: 32,
            stacks: 16,
        }
    }
}

impl PrimitiveMesh for Ellipsoid {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build(self)
    }
}

fn build(e: &Ellipsoid) -> Result<IndexedMesh, PrimitiveError> {
    if e.semi_x <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "semi_x must be > 0, got {}",
            e.semi_x
        )));
    }
    if e.semi_y <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "semi_y must be > 0, got {}",
            e.semi_y
        )));
    }
    if e.semi_z <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "semi_z must be > 0, got {}",
            e.semi_z
        )));
    }
    if e.segments < 3 {
        return Err(PrimitiveError::TooFewSegments(e.segments));
    }
    if e.stacks < 2 {
        return Err(PrimitiveError::InvalidParam(format!(
            "stacks must be ≥ 2, got {}",
            e.stacks
        )));
    }

    let region = RegionId::new(1);
    let mut mesh = IndexedMesh::new();
    let a = e.semi_x;
    let b = e.semi_y;
    let c = e.semi_z;
    let cx = e.center.x;
    let cy = e.center.y;
    let cz = e.center.z;

    // Sample (position, outward_normal) at (theta, phi).
    // Normal = gradient of implicit F = (x/a)²+(y/b)²+(z/c)²−1, normalised.
    let vertex_at = |theta: f64, phi: f64| -> (Point3r, Vector3r) {
        let sp = phi.sin();
        let cp = phi.cos();
        let ct = theta.cos();
        let st = theta.sin();
        let pos = Point3r::new(cx + a * sp * ct, cy + b * cp, cz + c * sp * st);
        // Raw gradient direction: (sp*ct/a, cp/b, sp*st/c)
        let nx = sp * ct / a;
        let ny = cp / b;
        let nz_c = sp * st / c;
        let len = (nx * nx + ny * ny + nz_c * nz_c).sqrt();
        let n = if len > 1e-15 {
            Vector3r::new(nx / len, ny / len, nz_c / len)
        } else {
            // Degenerate (pole with b=0 pathological case — prevented by validation)
            Vector3r::y()
        };
        (pos, n)
    };

    for i in 0..e.segments {
        let theta0 = i as f64 / e.segments as f64 * TAU;
        let theta1 = (i + 1) as f64 / e.segments as f64 * TAU;

        for j in 0..e.stacks {
            let phi0 = j as f64 / e.stacks as f64 * PI;
            let phi1 = (j + 1) as f64 / e.stacks as f64 * PI;

            let (pos00, n00) = vertex_at(theta0, phi0);
            let (pos10, n10) = vertex_at(theta1, phi0);
            let (pos11, n11) = vertex_at(theta1, phi1);
            let (pos01, n01) = vertex_at(theta0, phi1);

            let v00 = mesh.add_vertex(pos00, n00);
            let v10 = mesh.add_vertex(pos10, n10);
            let v11 = mesh.add_vertex(pos11, n11);
            let v01 = mesh.add_vertex(pos01, n01);

            if j == 0 {
                // North-pole row: v00 == v10 (both at north pole) → triangle
                mesh.add_face_with_region(v10, v11, v01, region);
            } else if j == e.stacks - 1 {
                // South-pole row: v01 == v11 (both at south pole) → triangle
                mesh.add_face_with_region(v00, v10, v01, region);
            } else {
                // Interior quad → 2 CCW triangles
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
    use crate::application::watertight::check::check_watertight;
    use crate::infrastructure::storage::edge_store::EdgeStore;
    use std::f64::consts::PI;

    #[test]
    fn ellipsoid_is_watertight() {
        let mesh = Ellipsoid::default().build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.is_watertight, "ellipsoid must be watertight");
        assert_eq!(report.euler_characteristic, Some(2));
    }

    #[test]
    fn ellipsoid_volume_positive_and_approximately_correct() {
        let (a, b, c) = (2.0_f64, 1.0_f64, 1.5_f64);
        let mesh = Ellipsoid {
            semi_x: a,
            semi_y: b,
            semi_z: c,
            segments: 64,
            stacks: 32,
            ..Ellipsoid::default()
        }
        .build()
        .unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.signed_volume > 0.0);
        let expected = 4.0 / 3.0 * PI * a * b * c;
        let error = (report.signed_volume - expected).abs() / expected;
        assert!(
            error < 0.005,
            "volume error {:.4}% should be < 0.5% at 64×32",
            error * 100.0
        );
    }

    #[test]
    fn ellipsoid_sphere_case_matches_uvsphere() {
        // When a=b=c=r the ellipsoid should match the UvSphere volume exactly.
        let r = 1.0_f64;
        let mesh = Ellipsoid {
            semi_x: r,
            semi_y: r,
            semi_z: r,
            segments: 64,
            stacks: 32,
            ..Ellipsoid::default()
        }
        .build()
        .unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        let expected = 4.0 * PI / 3.0 * r * r * r;
        let error = (report.signed_volume - expected).abs() / expected;
        assert!(error < 0.005, "sphere-case error {:.4}%", error * 100.0);
    }

    #[test]
    fn ellipsoid_rejects_invalid_params() {
        assert!(Ellipsoid {
            semi_x: 0.0,
            ..Ellipsoid::default()
        }
        .build()
        .is_err());
        assert!(Ellipsoid {
            semi_y: -1.0,
            ..Ellipsoid::default()
        }
        .build()
        .is_err());
        assert!(Ellipsoid {
            semi_z: 0.0,
            ..Ellipsoid::default()
        }
        .build()
        .is_err());
        assert!(Ellipsoid {
            segments: 2,
            ..Ellipsoid::default()
        }
        .build()
        .is_err());
        assert!(Ellipsoid {
            stacks: 1,
            ..Ellipsoid::default()
        }
        .build()
        .is_err());
    }
}
