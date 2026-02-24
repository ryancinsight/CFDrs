//! Ring torus primitive.

use std::f64::consts::TAU;

use super::{PrimitiveError, PrimitiveMesh};
use crate::domain::core::index::RegionId;
use crate::domain::core::scalar::{Point3r, Vector3r};
use crate::domain::mesh::IndexedMesh;

/// Builds a ring torus centred at the origin, lying in the XZ plane.
///
/// ## Parametrisation
///
/// ```text
/// x(φ, θ) = (R + r·cos θ) · cos φ
/// y(φ, θ) =  r · sin θ
/// z(φ, θ) = (R + r·cos θ) · sin φ
///
/// φ ∈ [0, 2π)  — major angle (around the ring)
/// θ ∈ [0, 2π)  — minor angle (around the tube cross-section)
/// ```
///
/// ## Topology
///
/// The torus has genus 1, so `V − E + F = 0` (not 2).
/// [`check_watertight`][crate::application::watertight::check::check_watertight] will
/// report `euler_characteristic = 0` and `is_watertight = true` (the mesh
/// *is* closed and oriented; only the Euler check differs from a sphere).
///
/// ## Output
///
/// - `2 × major_segments × minor_segments` faces
/// - `signed_volume ≈ 2π² R r²` (positive for outward normals)
/// - All faces tagged `RegionId(1)`
///
/// ## Example
///
/// ```rust,ignore
/// let torus = Torus { major_radius: 3.0, minor_radius: 1.0,
///     major_segments: 48, minor_segments: 24 };
/// let mesh = torus.build().unwrap();
/// ```
#[derive(Clone, Debug)]
pub struct Torus {
    /// Distance from the centre of the tube to the centre of the torus [mm].
    pub major_radius: f64,
    /// Radius of the tube cross-section [mm].
    pub minor_radius: f64,
    /// Number of divisions around the ring (≥ 3).
    pub major_segments: usize,
    /// Number of divisions around the tube cross-section (≥ 3).
    pub minor_segments: usize,
}

impl Default for Torus {
    fn default() -> Self {
        Self {
            major_radius: 3.0,
            minor_radius: 1.0,
            major_segments: 48,
            minor_segments: 24,
        }
    }
}

impl PrimitiveMesh for Torus {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build(self)
    }
}

fn build(t: &Torus) -> Result<IndexedMesh, PrimitiveError> {
    if t.major_radius <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "major_radius must be > 0, got {}",
            t.major_radius
        )));
    }
    if t.minor_radius <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "minor_radius must be > 0, got {}",
            t.minor_radius
        )));
    }
    if t.minor_radius >= t.major_radius {
        return Err(PrimitiveError::InvalidParam(format!(
            "minor_radius ({}) must be < major_radius ({}) for a ring torus",
            t.minor_radius, t.major_radius
        )));
    }
    if t.major_segments < 3 {
        return Err(PrimitiveError::TooFewSegments(t.major_segments));
    }
    if t.minor_segments < 3 {
        return Err(PrimitiveError::TooFewSegments(t.minor_segments));
    }

    let region = RegionId::new(1);
    let mut mesh = IndexedMesh::new();
    let r_maj = t.major_radius;
    let r_min = t.minor_radius;

    // Sample (position, outward_normal) at (φ, θ).
    //
    // The tube-normal at angle θ on the cross-section at ring angle φ is:
    //   n = (cos θ · cos φ,  sin θ,  cos θ · sin φ)
    // (The "outward" direction from the tube axis toward the tube surface.)
    let vertex_at = |phi: f64, theta: f64| -> (Point3r, Vector3r) {
        let cp = phi.cos();
        let sp = phi.sin();
        let ct = theta.cos();
        let st = theta.sin();
        // Outward normal from the tube centre
        let n = Vector3r::new(ct * cp, st, ct * sp);
        // Position on torus surface
        let rho = r_maj + r_min * ct; // distance from torus axis (Y)
        let pos = Point3r::new(rho * cp, r_min * st, rho * sp);
        (pos, n)
    };

    // Build a grid of quads.  Index layout:
    //   (i, j) → vertex at (φ_i, θ_j)
    //
    // Outward CCW quad (viewed from outside the tube):
    //   (i, j) → (i+1, j) → (i+1, j+1) → (i, j+1)
    // which in indices wraps via modulo.
    for i in 0..t.major_segments {
        let phi0 = i as f64 / t.major_segments as f64 * TAU;
        let phi1 = (i + 1) as f64 / t.major_segments as f64 * TAU;

        for j in 0..t.minor_segments {
            let theta0 = j as f64 / t.minor_segments as f64 * TAU;
            let theta1 = (j + 1) as f64 / t.minor_segments as f64 * TAU;

            let (pos00, n00) = vertex_at(phi0, theta0);
            let (pos10, n10) = vertex_at(phi1, theta0);
            let (pos11, n11) = vertex_at(phi1, theta1);
            let (pos01, n01) = vertex_at(phi0, theta1);

            let v00 = mesh.add_vertex(pos00, n00);
            let v10 = mesh.add_vertex(pos10, n10);
            let v11 = mesh.add_vertex(pos11, n11);
            let v01 = mesh.add_vertex(pos01, n01);

            // CCW from outside the tube.
            // The outward-facing quad (φ increases, θ increases) winds as:
            //   v00 → v01 → v11 → v10   (θ increases first, then φ)
            // Split into two CCW triangles:
            //   tri 1: v00 → v01 → v11
            //   tri 2: v00 → v11 → v10
            mesh.add_face_with_region(v00, v01, v11, region);
            mesh.add_face_with_region(v00, v11, v10, region);
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
    fn torus_is_closed_and_oriented() {
        let mesh = Torus::default().build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        // Torus is closed (no boundary edges) and consistently oriented
        assert!(report.is_closed, "torus must be closed");
        assert!(
            report.orientation_consistent,
            "torus must be consistently oriented"
        );
        // Euler characteristic for a torus (genus 1) = 0
        assert_eq!(
            report.euler_characteristic,
            Some(0),
            "torus Euler characteristic must be 0 (genus 1)"
        );
        // The is_watertight flag checks closed && oriented — both should be true
        assert!(report.is_watertight, "torus passes watertight check");
    }

    #[test]
    fn torus_volume_positive_and_approximately_correct() {
        let r_maj = 3.0_f64;
        let r_min = 1.0_f64;
        let t = Torus {
            major_radius: r_maj,
            minor_radius: r_min,
            major_segments: 96,
            minor_segments: 48,
        };
        let mesh = t.build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(
            report.signed_volume > 0.0,
            "outward normals → positive volume"
        );
        let expected = 2.0 * PI * PI * r_maj * r_min * r_min;
        let error = (report.signed_volume - expected).abs() / expected;
        assert!(
            error < 0.01,
            "volume error {:.4}% should be < 1.0% at 96×48",
            error * 100.0
        );
    }

    #[test]
    fn torus_minor_radius_exceeds_major() {
        let result = Torus {
            major_radius: 1.0,
            minor_radius: 2.0,
            ..Torus::default()
        }
        .build();
        assert!(result.is_err());
    }
}
