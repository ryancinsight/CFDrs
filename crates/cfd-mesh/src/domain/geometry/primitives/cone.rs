//! Right circular cone primitive.

use std::f64::consts::TAU;

use super::{PrimitiveError, PrimitiveMesh};
use crate::domain::core::index::RegionId;
use crate::domain::core::scalar::{Point3r, Vector3r};
use crate::domain::mesh::IndexedMesh;

/// Builds a closed right circular cone.
///
/// The cone is aligned with the +Y axis: the base circle is at `base_center`
/// and the apex is at `base_center + (0, height, 0)`.
///
/// ## Lateral normal formula
///
/// At a base point `P = (r cosθ, 0, r sinθ)`, the outward lateral surface
/// normal is:
///
/// ```text
/// n = normalize(h · r̂ + r · ĵ)
///   where r̂ = (cosθ, 0, sinθ)  and  ĵ = (0, 1, 0)
/// ```
///
/// (The slant length `l = √(r² + h²)` is the denominator.)
///
/// ## Output
///
/// - `2 × segments` faces (segments lateral + segments base)
/// - `signed_volume ≈ π r² h / 3` (positive for outward normals)
#[derive(Clone, Debug)]
pub struct Cone {
    /// Base circle centre.
    pub base_center: Point3r,
    /// Base circle radius [mm].
    pub radius: f64,
    /// Height from base to apex [mm] (extends along +Y).
    pub height: f64,
    /// Number of angular subdivisions (≥ 3).
    pub segments: usize,
}

impl Default for Cone {
    fn default() -> Self {
        Self {
            base_center: Point3r::origin(),
            radius: 1.0,
            height: 2.0,
            segments: 32,
        }
    }
}

impl PrimitiveMesh for Cone {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build(self)
    }
}

fn build(c: &Cone) -> Result<IndexedMesh, PrimitiveError> {
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

    // Slant length used to normalise lateral normals.
    let slant = (r * r + h * h).sqrt();

    let apex_pos = Point3r::new(bx, by + h, bz);

    // ── Lateral surface ──────────────────────────────────────────────────────
    // Outward winding (CCW from outside = apex → base1 → base0):
    // When facing the cone from outside at angle midpoint, the base ring
    // goes CCW as: base0 → base1 (increasing θ), so the outward-facing
    // triangle apex → base1 → base0 follows the right-hand rule correctly.
    for i in 0..c.segments {
        let a0 = i as f64 / c.segments as f64 * TAU;
        let a1 = (i + 1) as f64 / c.segments as f64 * TAU;

        let (c0, s0) = (a0.cos(), a0.sin());
        let (c1, s1) = (a1.cos(), a1.sin());

        let n0 = Vector3r::new(h * c0 / slant, r / slant, h * s0 / slant);
        let n1 = Vector3r::new(h * c1 / slant, r / slant, h * s1 / slant);
        let n_apex = ((n0 + n1) * 0.5).normalize();

        let p_base0 = Point3r::new(bx + r * c0, by, bz + r * s0);
        let p_base1 = Point3r::new(bx + r * c1, by, bz + r * s1);

        let vapex = mesh.add_vertex(apex_pos, n_apex);
        let vb0 = mesh.add_vertex(p_base0, n0);
        let vb1 = mesh.add_vertex(p_base1, n1);

        // CCW from outside: apex → base1 → base0
        mesh.add_face_with_region(vapex, vb1, vb0, region);
    }

    // ── Base cap (y = by, normal −Y) ─────────────────────────────────────────
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

    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::application::watertight::check::check_watertight;
    use crate::infrastructure::storage::edge_store::EdgeStore;
    use std::f64::consts::PI;

    #[test]
    fn cone_is_watertight() {
        let mesh = Cone::default().build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.is_watertight, "cone must be watertight");
        assert_eq!(report.euler_characteristic, Some(2));
    }

    #[test]
    fn cone_volume_positive_and_approximately_correct() {
        let c = Cone {
            radius: 1.0,
            height: 2.0,
            segments: 64,
            ..Cone::default()
        };
        let mesh = c.build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.signed_volume > 0.0);
        let expected = PI * 1.0_f64.powi(2) * 2.0 / 3.0;
        let error = (report.signed_volume - expected).abs() / expected;
        assert!(
            error < 0.005,
            "volume error {:.4}% should be < 0.5%",
            error * 100.0
        );
    }
}
