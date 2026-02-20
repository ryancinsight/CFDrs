//! Hollow cylinder (pipe) primitive.

use std::f64::consts::TAU;

use crate::core::index::RegionId;
use crate::core::scalar::{Point3r, Vector3r};
use crate::mesh::IndexedMesh;
use super::{PrimitiveMesh, PrimitiveError};

/// Builds a hollow right circular cylinder (pipe / annular tube).
///
/// The pipe is aligned with the +Y axis: the bottom annulus is centred at
/// `base_center` and the top annulus is at `base_center + (0, height, 0)`.
///
/// ## Mesh structure
///
/// | Section | Faces | Normal direction |
/// |---------|-------|-----------------|
/// | Outer lateral | `2 × segments` | Radially outward |
/// | Inner lateral | `2 × segments` | Toward bore axis (inward) |
/// | Bottom annular cap | `2 × segments` | −Y |
/// | Top annular cap | `2 × segments` | +Y |
///
/// ## Topology
///
/// The bore creates a topological genus-1 handle, so `V − E + F = 0` (like
/// [`Torus`]).  [`check_watertight`][crate::watertight::check::check_watertight]
/// reports `euler_characteristic = 0` and `is_watertight = true`.
///
/// ## Uses
///
/// PDMS soft-lithography tubing walls, FSI vascular models, concentric-cylinder
/// Couette flow devices, hollow microsphere shells.
///
/// ## Output
///
/// - `signed_volume = π (r_outer² − r_inner²) h`
/// - All faces tagged `RegionId(1)`
#[derive(Clone, Debug)]
pub struct Pipe {
    /// Base annulus centre (bottom face of pipe).
    pub base_center: Point3r,
    /// Bore (inner) radius [mm].
    pub inner_radius: f64,
    /// Outer wall radius [mm]. Must be > `inner_radius`.
    pub outer_radius: f64,
    /// Height [mm] (extends along +Y from `base_center`).
    pub height: f64,
    /// Number of angular subdivisions (≥ 3).
    pub segments: usize,
}

impl Default for Pipe {
    fn default() -> Self {
        Self {
            base_center: Point3r::origin(),
            inner_radius: 0.5,
            outer_radius: 1.0,
            height: 2.0,
            segments: 32,
        }
    }
}

impl PrimitiveMesh for Pipe {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build(self)
    }
}

fn build(p: &Pipe) -> Result<IndexedMesh, PrimitiveError> {
    if p.inner_radius <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "inner_radius must be > 0, got {}", p.inner_radius
        )));
    }
    if p.outer_radius <= p.inner_radius {
        return Err(PrimitiveError::InvalidParam(format!(
            "outer_radius ({}) must be > inner_radius ({})",
            p.outer_radius, p.inner_radius
        )));
    }
    if p.height <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "height must be > 0, got {}", p.height
        )));
    }
    if p.segments < 3 {
        return Err(PrimitiveError::TooFewSegments(p.segments));
    }

    let region = RegionId::new(1);
    let mut mesh = IndexedMesh::new();
    let ri = p.inner_radius;
    let ro = p.outer_radius;
    let h  = p.height;
    let bx = p.base_center.x;
    let by = p.base_center.y;
    let bz = p.base_center.z;

    // ── Outer lateral (normal = +radial) ─────────────────────────────────────
    // Same winding as Cylinder: bot0 → top0 → top1, bot0 → top1 → bot1
    for i in 0..p.segments {
        let a0 = i as f64 / p.segments as f64 * TAU;
        let a1 = (i + 1) as f64 / p.segments as f64 * TAU;
        let (c0, s0) = (a0.cos(), a0.sin());
        let (c1, s1) = (a1.cos(), a1.sin());
        let n0 = Vector3r::new(c0, 0.0, s0);
        let n1 = Vector3r::new(c1, 0.0, s1);
        let vb0 = mesh.add_vertex(Point3r::new(bx + ro * c0, by,     bz + ro * s0), n0);
        let vb1 = mesh.add_vertex(Point3r::new(bx + ro * c1, by,     bz + ro * s1), n1);
        let vt0 = mesh.add_vertex(Point3r::new(bx + ro * c0, by + h, bz + ro * s0), n0);
        let vt1 = mesh.add_vertex(Point3r::new(bx + ro * c1, by + h, bz + ro * s1), n1);
        mesh.add_face_with_region(vb0, vt0, vt1, region);
        mesh.add_face_with_region(vb0, vt1, vb1, region);
    }

    // ── Inner lateral (normal = −radial, pointing toward bore axis) ──────────
    // Reversed winding so that the outward normal points into the bore.
    for i in 0..p.segments {
        let a0 = i as f64 / p.segments as f64 * TAU;
        let a1 = (i + 1) as f64 / p.segments as f64 * TAU;
        let (c0, s0) = (a0.cos(), a0.sin());
        let (c1, s1) = (a1.cos(), a1.sin());
        // Inward-facing normal (toward axis)
        let n0 = Vector3r::new(-c0, 0.0, -s0);
        let n1 = Vector3r::new(-c1, 0.0, -s1);
        let vb0 = mesh.add_vertex(Point3r::new(bx + ri * c0, by,     bz + ri * s0), n0);
        let vb1 = mesh.add_vertex(Point3r::new(bx + ri * c1, by,     bz + ri * s1), n1);
        let vt0 = mesh.add_vertex(Point3r::new(bx + ri * c0, by + h, bz + ri * s0), n0);
        let vt1 = mesh.add_vertex(Point3r::new(bx + ri * c1, by + h, bz + ri * s1), n1);
        // Reversed winding compared to outer lateral
        mesh.add_face_with_region(vb0, vt1, vt0, region);
        mesh.add_face_with_region(vb0, vb1, vt1, region);
    }

    // ── Bottom annular cap (y = by, normal −Y) ───────────────────────────────
    // CCW from below: inner_{i+1} → inner_i → outer_i, inner_{i+1} → outer_i → outer_{i+1}
    {
        let n_down = -Vector3r::y();
        for i in 0..p.segments {
            let a0 = i as f64 / p.segments as f64 * TAU;
            let a1 = (i + 1) as f64 / p.segments as f64 * TAU;
            let oi  = mesh.add_vertex(Point3r::new(bx + ro * a0.cos(), by, bz + ro * a0.sin()), n_down);
            let oi1 = mesh.add_vertex(Point3r::new(bx + ro * a1.cos(), by, bz + ro * a1.sin()), n_down);
            let ii  = mesh.add_vertex(Point3r::new(bx + ri * a0.cos(), by, bz + ri * a0.sin()), n_down);
            let ii1 = mesh.add_vertex(Point3r::new(bx + ri * a1.cos(), by, bz + ri * a1.sin()), n_down);
            mesh.add_face_with_region(ii1, ii, oi, region);
            mesh.add_face_with_region(ii1, oi, oi1, region);
        }
    }

    // ── Top annular cap (y = by + h, normal +Y) ──────────────────────────────
    // CCW from above: outer_i → inner_i → inner_{i+1}, outer_i → inner_{i+1} → outer_{i+1}
    {
        let n_up = Vector3r::y();
        for i in 0..p.segments {
            let a0 = i as f64 / p.segments as f64 * TAU;
            let a1 = (i + 1) as f64 / p.segments as f64 * TAU;
            let oi  = mesh.add_vertex(Point3r::new(bx + ro * a0.cos(), by + h, bz + ro * a0.sin()), n_up);
            let oi1 = mesh.add_vertex(Point3r::new(bx + ro * a1.cos(), by + h, bz + ro * a1.sin()), n_up);
            let ii  = mesh.add_vertex(Point3r::new(bx + ri * a0.cos(), by + h, bz + ri * a0.sin()), n_up);
            let ii1 = mesh.add_vertex(Point3r::new(bx + ri * a1.cos(), by + h, bz + ri * a1.sin()), n_up);
            mesh.add_face_with_region(oi, ii, ii1, region);
            mesh.add_face_with_region(oi, ii1, oi1, region);
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
    fn pipe_is_closed_and_oriented() {
        let mesh = Pipe::default().build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.is_closed, "pipe must be closed");
        assert!(report.orientation_consistent, "pipe must be consistently oriented");
        // Genus 1 (through-bore) → Euler characteristic = 0
        assert_eq!(report.euler_characteristic, Some(0),
            "pipe Euler characteristic must be 0 (genus 1)");
        assert!(report.is_watertight, "pipe passes watertight check");
    }

    #[test]
    fn pipe_volume_positive_and_approximately_correct() {
        let (ri, ro, h) = (0.5_f64, 1.0_f64, 2.0_f64);
        let mesh = Pipe {
            inner_radius: ri, outer_radius: ro, height: h,
            segments: 64, ..Pipe::default()
        }.build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.signed_volume > 0.0);
        let expected = PI * (ro * ro - ri * ri) * h;
        let error = (report.signed_volume - expected).abs() / expected;
        assert!(error < 0.005, "volume error {:.4}% should be < 0.5%", error * 100.0);
    }

    #[test]
    fn pipe_rejects_invalid_params() {
        assert!(Pipe { inner_radius: 0.0, ..Pipe::default() }.build().is_err());
        assert!(Pipe { outer_radius: 0.4, inner_radius: 0.5, ..Pipe::default() }.build().is_err());
        assert!(Pipe { outer_radius: 0.5, inner_radius: 0.5, ..Pipe::default() }.build().is_err());
        assert!(Pipe { height: 0.0, ..Pipe::default() }.build().is_err());
        assert!(Pipe { segments: 2, ..Pipe::default() }.build().is_err());
    }
}
