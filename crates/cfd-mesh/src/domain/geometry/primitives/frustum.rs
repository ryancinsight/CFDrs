//! Truncated right circular cone (frustum) primitive.

use std::f64::consts::TAU;

use super::{PrimitiveError, PrimitiveMesh};
use crate::domain::core::index::RegionId;
use crate::domain::core::scalar::{Point3r, Vector3r};
use crate::domain::mesh::IndexedMesh;

/// Builds a closed truncated right circular cone (frustum).
///
/// The frustum is aligned with the +Y axis: the bottom circle of radius
/// `bottom_radius` is centred at `base_center`, and the top circle of radius
/// `top_radius` is at `base_center + (0, height, 0)`.
///
/// ## Special cases
///
/// | Condition | Equivalent shape |
/// |-----------|-----------------|
/// | `top_radius == bottom_radius` | Cylinder |
/// | `top_radius = 0` | Invalid — use [`Cone`] instead |
///
/// ## Lateral normal formula
///
/// At angle `θ` on the lateral surface, the outward normal is:
///
/// ```text
/// n = normalize( h·cos(θ),  r0 − r1,  h·sin(θ) )
/// ```
///
/// where `h = height`, `r0 = bottom_radius`, `r1 = top_radius`.
///
/// ## Output
///
/// - `4 × segments` faces  (`2×segments` lateral + `segments` bottom + `segments` top)
/// - `signed_volume = (π·h/3)·(r0² + r0·r1 + r1²)` (Prismatoid formula, exact)
/// - All faces tagged `RegionId(1)`
#[derive(Clone, Debug)]
pub struct Frustum {
    /// Base circle centre.
    pub base_center: Point3r,
    /// Bottom circle radius [mm] (at `y = base_center.y`).
    pub bottom_radius: f64,
    /// Top circle radius [mm] (at `y = base_center.y + height`). Must be > 0.
    pub top_radius: f64,
    /// Height [mm] (extends along +Y from `base_center`).
    pub height: f64,
    /// Number of angular subdivisions (≥ 3).
    pub segments: usize,
}

impl Default for Frustum {
    fn default() -> Self {
        Self {
            base_center: Point3r::origin(),
            bottom_radius: 1.0,
            top_radius: 0.5,
            height: 2.0,
            segments: 32,
        }
    }
}

impl PrimitiveMesh for Frustum {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build(self)
    }
}

fn build(f: &Frustum) -> Result<IndexedMesh, PrimitiveError> {
    if f.bottom_radius <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "bottom_radius must be > 0, got {}",
            f.bottom_radius
        )));
    }
    if f.top_radius <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "top_radius must be > 0, got {} — for a cone use `Cone` instead",
            f.top_radius
        )));
    }
    if f.height <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "height must be > 0, got {}",
            f.height
        )));
    }
    if f.segments < 3 {
        return Err(PrimitiveError::TooFewSegments(f.segments));
    }

    let region = RegionId::new(1);
    let mut mesh = IndexedMesh::new();
    let r0 = f.bottom_radius;
    let r1 = f.top_radius;
    let h = f.height;
    let bx = f.base_center.x;
    let by = f.base_center.y;
    let bz = f.base_center.z;

    // Slant normal denominator: ||(h·r̂, r0−r1, h·ŝ)|| = √(h² + (r0−r1)²)
    let slant = (h * h + (r0 - r1) * (r0 - r1)).sqrt();

    // ── Lateral surface ──────────────────────────────────────────────────────
    for i in 0..f.segments {
        let a0 = i as f64 / f.segments as f64 * TAU;
        let a1 = (i + 1) as f64 / f.segments as f64 * TAU;

        let (c0, s0) = (a0.cos(), a0.sin());
        let (c1, s1) = (a1.cos(), a1.sin());

        // Outward slant normal: n = (h·cosθ, r0−r1, h·sinθ) / slant
        let n0 = Vector3r::new(h * c0 / slant, (r0 - r1) / slant, h * s0 / slant);
        let n1 = Vector3r::new(h * c1 / slant, (r0 - r1) / slant, h * s1 / slant);

        let p_bot0 = Point3r::new(bx + r0 * c0, by, bz + r0 * s0);
        let p_bot1 = Point3r::new(bx + r0 * c1, by, bz + r0 * s1);
        let p_top0 = Point3r::new(bx + r1 * c0, by + h, bz + r1 * s0);
        let p_top1 = Point3r::new(bx + r1 * c1, by + h, bz + r1 * s1);

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
        for i in 0..f.segments {
            let a0 = i as f64 / f.segments as f64 * TAU;
            let a1 = (i + 1) as f64 / f.segments as f64 * TAU;
            let p0 = Point3r::new(bx + r0 * a0.cos(), by, bz + r0 * a0.sin());
            let p1 = Point3r::new(bx + r0 * a1.cos(), by, bz + r0 * a1.sin());
            let v0 = mesh.add_vertex(p0, n_down);
            let v1 = mesh.add_vertex(p1, n_down);
            // CCW from below: centre → p0 → p1
            mesh.add_face_with_region(vc, v0, v1, region);
        }
    }

    // ── Top cap (y = by + h, normal +Y) ─────────────────────────────────────
    {
        let n_up = Vector3r::y();
        let center = Point3r::new(bx, by + h, bz);
        let vc = mesh.add_vertex(center, n_up);
        for i in 0..f.segments {
            let a0 = i as f64 / f.segments as f64 * TAU;
            let a1 = (i + 1) as f64 / f.segments as f64 * TAU;
            let p0 = Point3r::new(bx + r1 * a0.cos(), by + h, bz + r1 * a0.sin());
            let p1 = Point3r::new(bx + r1 * a1.cos(), by + h, bz + r1 * a1.sin());
            let v0 = mesh.add_vertex(p0, n_up);
            let v1 = mesh.add_vertex(p1, n_up);
            // CCW from above: centre → p1 → p0
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
    fn frustum_is_watertight() {
        let mesh = Frustum::default().build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.is_watertight, "frustum must be watertight");
        assert_eq!(report.euler_characteristic, Some(2));
    }

    #[test]
    fn frustum_volume_positive_and_approximately_correct() {
        let (r0, r1, h) = (1.0_f64, 0.5_f64, 2.0_f64);
        let mesh = Frustum {
            bottom_radius: r0,
            top_radius: r1,
            height: h,
            segments: 64,
            ..Frustum::default()
        }
        .build()
        .unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.signed_volume > 0.0);
        // Prismatoid formula: V = (π·h/3)·(r0² + r0·r1 + r1²)
        let expected = PI * h / 3.0 * (r0 * r0 + r0 * r1 + r1 * r1);
        let error = (report.signed_volume - expected).abs() / expected;
        assert!(
            error < 0.005,
            "volume error {:.4}% should be < 0.5%",
            error * 100.0
        );
    }

    #[test]
    fn frustum_cylinder_case() {
        // r0 == r1 → degenerate frustum = cylinder
        let r = 1.0_f64;
        let h = 2.0_f64;
        let mesh = Frustum {
            bottom_radius: r,
            top_radius: r,
            height: h,
            segments: 64,
            ..Frustum::default()
        }
        .build()
        .unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        let expected = PI * r * r * h;
        let error = (report.signed_volume - expected).abs() / expected;
        assert!(error < 0.005, "cylinder-case error {:.4}%", error * 100.0);
    }

    #[test]
    fn frustum_rejects_invalid_params() {
        assert!(Frustum {
            bottom_radius: 0.0,
            ..Frustum::default()
        }
        .build()
        .is_err());
        assert!(Frustum {
            top_radius: 0.0,
            ..Frustum::default()
        }
        .build()
        .is_err());
        assert!(Frustum {
            top_radius: -1.0,
            ..Frustum::default()
        }
        .build()
        .is_err());
        assert!(Frustum {
            height: 0.0,
            ..Frustum::default()
        }
        .build()
        .is_err());
        assert!(Frustum {
            segments: 2,
            ..Frustum::default()
        }
        .build()
        .is_err());
    }
}
