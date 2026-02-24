//! Closed right circular cylinder primitive.

use std::f64::consts::TAU;

use super::{PrimitiveError, PrimitiveMesh};
use crate::domain::core::index::RegionId;
use crate::domain::core::scalar::{Point3r, Vector3r};
use crate::domain::mesh::IndexedMesh;

/// Builds a closed right circular cylinder.
///
/// The cylinder is aligned with the +Y axis, with the base centre at
/// `base_center` and the top at `base_center + (0, height, 0)`.
///
/// ## Mesh structure
///
/// - **Lateral surface**: `segments` quads → `2 × segments` triangles
/// - **Bottom cap** (−Y normal): `segments` triangles
/// - **Top cap** (+Y normal): `segments` triangles
/// - Total: `4 × segments` faces
///
/// ## Output
///
/// - `signed_volume ≈ π r² h` (error < 0.5% for segments ≥ 32)
/// - All three parts carry `RegionId(1)`
///   (use [`CylinderRegions`] variant for named wall/inlet/outlet regions)
#[derive(Clone, Debug)]
pub struct Cylinder {
    /// Base circle centre.
    pub base_center: Point3r,
    /// Cylinder radius [mm].
    pub radius: f64,
    /// Cylinder height [mm] (extends along +Y from `base_center`).
    pub height: f64,
    /// Number of angular subdivisions (≥ 3).
    pub segments: usize,
}

impl Default for Cylinder {
    fn default() -> Self {
        Self {
            base_center: Point3r::origin(),
            radius: 1.0,
            height: 2.0,
            segments: 32,
        }
    }
}

impl PrimitiveMesh for Cylinder {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build(self)
    }
}

fn build(c: &Cylinder) -> Result<IndexedMesh, PrimitiveError> {
    if c.radius <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "radius must be > 0, got {}",
            c.radius
        )));
    }
    if c.height <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "height must be > 0, got {}",
            c.height
        )));
    }
    if c.segments < 3 {
        return Err(PrimitiveError::TooFewSegments(c.segments));
    }

    let region = RegionId::new(1);
    let mut mesh = IndexedMesh::new();
    let r = c.radius;
    let h = c.height;
    let bx = c.base_center.x;
    let by = c.base_center.y;
    let bz = c.base_center.z;

    // ── Lateral surface ──────────────────────────────────────────────────────
    // Viewed from outside (+radial), the quad winds CCW as:
    //   bot0 → top0 → top1 → bot1
    for i in 0..c.segments {
        let a0 = i as f64 / c.segments as f64 * TAU;
        let a1 = (i + 1) as f64 / c.segments as f64 * TAU;

        let (c0, s0) = (a0.cos(), a0.sin());
        let (c1, s1) = (a1.cos(), a1.sin());

        let n0 = Vector3r::new(c0, 0.0, s0);
        let n1 = Vector3r::new(c1, 0.0, s1);

        let p_bot0 = Point3r::new(bx + r * c0, by, bz + r * s0);
        let p_bot1 = Point3r::new(bx + r * c1, by, bz + r * s1);
        let p_top0 = Point3r::new(bx + r * c0, by + h, bz + r * s0);
        let p_top1 = Point3r::new(bx + r * c1, by + h, bz + r * s1);

        let vb0 = mesh.add_vertex(p_bot0, n0);
        let vb1 = mesh.add_vertex(p_bot1, n1);
        let vt0 = mesh.add_vertex(p_top0, n0);
        let vt1 = mesh.add_vertex(p_top1, n1);

        // CCW from outside: bot0 → top0 → top1, then bot0 → top1 → bot1
        mesh.add_face_with_region(vb0, vt0, vt1, region);
        mesh.add_face_with_region(vb0, vt1, vb1, region);
    }

    // ── Bottom cap (y = by, normal −Y) ───────────────────────────────────────
    {
        let n_down = -Vector3r::y();
        let center = Point3r::new(bx, by, bz);
        let vc = mesh.add_vertex(center, n_down);

        for i in 0..c.segments {
            let a0 = i as f64 / c.segments as f64 * TAU;
            let a1 = (i + 1) as f64 / c.segments as f64 * TAU;
            let p0 = Point3r::new(bx + r * a0.cos(), by, bz + r * a0.sin());
            let p1 = Point3r::new(bx + r * a1.cos(), by, bz + r * a1.sin());
            let v0 = mesh.add_vertex(p0, n_down);
            let v1 = mesh.add_vertex(p1, n_down);
            // CCW from below (−Y outward): centre → p0 → p1
            mesh.add_face_with_region(vc, v0, v1, region);
        }
    }

    // ── Top cap (y = by + h, normal +Y) ─────────────────────────────────────
    {
        let n_up = Vector3r::y();
        let center = Point3r::new(bx, by + h, bz);
        let vc = mesh.add_vertex(center, n_up);

        for i in 0..c.segments {
            let a0 = i as f64 / c.segments as f64 * TAU;
            let a1 = (i + 1) as f64 / c.segments as f64 * TAU;
            let p0 = Point3r::new(bx + r * a0.cos(), by + h, bz + r * a0.sin());
            let p1 = Point3r::new(bx + r * a1.cos(), by + h, bz + r * a1.sin());
            let v0 = mesh.add_vertex(p0, n_up);
            let v1 = mesh.add_vertex(p1, n_up);
            // CCW from above (+Y outward): centre → p1 → p0
            mesh.add_face_with_region(vc, v1, v0, region);
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
    fn cylinder_is_watertight() {
        let mesh = Cylinder::default().build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.is_watertight, "cylinder must be watertight");
        assert_eq!(report.euler_characteristic, Some(2));
    }

    #[test]
    fn cylinder_volume_positive_and_approximately_correct() {
        let c = Cylinder {
            radius: 1.0,
            height: 2.0,
            segments: 64,
            ..Cylinder::default()
        };
        let mesh = c.build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.signed_volume > 0.0);
        let expected = PI * 1.0_f64.powi(2) * 2.0;
        let error = (report.signed_volume - expected).abs() / expected;
        assert!(
            error < 0.005,
            "volume error {:.4}% should be < 0.5%",
            error * 100.0
        );
    }

    #[test]
    fn cylinder_too_few_segments() {
        let result = Cylinder {
            segments: 2,
            ..Cylinder::default()
        }
        .build();
        assert!(result.is_err());
    }
}
