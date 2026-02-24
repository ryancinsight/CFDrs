//! Hollow sphere (spherical shell) primitive.

use std::f64::consts::PI;
use std::f64::consts::TAU;

use super::{PrimitiveError, PrimitiveMesh};
use crate::domain::core::index::RegionId;
use crate::domain::core::scalar::{Point3r, Vector3r};
use crate::domain::mesh::IndexedMesh;

/// Builds a hollow sphere — two concentric sphere surfaces connected by
/// annular polar caps.
///
/// The shell is centred at `center`.  The outer surface has outward normals
/// pointing *away* from the centre; the inner surface has outward normals
/// pointing *toward* the centre (i.e., outward from the enclosed solid
/// material).
///
/// ## Construction
///
/// Both sphere surfaces are UV-parametrised over φ ∈ [φ₁, π − φ₁] where
/// φ₁ = π / stacks (one step from each pole).  This leaves a small polar
/// opening at each end that is closed by an annular quad strip — exactly
/// the same topology as [`Pipe`].
///
/// ## Topology
///
/// The through-bore (the polar hole connecting outer and inner surfaces)
/// creates a genus-1 handle, so `V − E + F = 0`  (χ = 0), identical to
/// [`Torus`] and [`Pipe`].
///
/// ## Output
///
/// - `signed_volume = (4/3)·π·(r_outer³ − r_inner³)`
/// - All faces tagged `RegionId(1)`
///
/// [`Pipe`]: super::pipe::Pipe
#[derive(Clone, Debug)]
pub struct SphericalShell {
    /// Centre of the shell.
    pub center: Point3r,
    /// Outer sphere radius [mm].
    pub outer_radius: f64,
    /// Inner sphere (cavity) radius [mm]. Must be < `outer_radius`.
    pub inner_radius: f64,
    /// Angular subdivisions around the equator (≥ 3).
    pub segments: usize,
    /// Latitude subdivisions (≥ 3). Determines shell resolution.
    pub stacks: usize,
}

impl Default for SphericalShell {
    fn default() -> Self {
        Self {
            center: Point3r::origin(),
            outer_radius: 1.0,
            inner_radius: 0.9,
            segments: 32,
            stacks: 16,
        }
    }
}

impl PrimitiveMesh for SphericalShell {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build(self)
    }
}

fn build(s: &SphericalShell) -> Result<IndexedMesh, PrimitiveError> {
    if s.inner_radius <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "inner_radius must be > 0, got {}",
            s.inner_radius
        )));
    }
    if s.outer_radius <= s.inner_radius {
        return Err(PrimitiveError::InvalidParam(format!(
            "outer_radius ({}) must be > inner_radius ({})",
            s.outer_radius, s.inner_radius
        )));
    }
    if s.segments < 3 {
        return Err(PrimitiveError::TooFewSegments(s.segments));
    }
    if s.stacks < 3 {
        return Err(PrimitiveError::InvalidParam("stacks must be ≥ 3".into()));
    }

    let region = RegionId::new(1);
    let mut mesh = IndexedMesh::new();

    let ro = s.outer_radius;
    let ri = s.inner_radius;
    let cx = s.center.x;
    let cy = s.center.y;
    let cz = s.center.z;
    let ns = s.segments;
    let nk = s.stacks;

    // phi ranges over [phi1, pi - phi1] where phi1 = pi/nk (one step from each pole).
    // This creates nk-1 latitudinal rings (indices 1 .. nk-1 inclusive).

    // Build outer ring VertexId arrays (phi rings 1 to nk-1)
    let mut outer_rings: Vec<Vec<crate::domain::core::index::VertexId>> =
        Vec::with_capacity(nk - 1);
    for k in 0..nk - 1 {
        let phi = (k + 1) as f64 / nk as f64 * PI;
        let sp = phi.sin();
        let cp = phi.cos();
        let y = cy + ro * cp;
        let row: Vec<_> = (0..ns)
            .map(|j| {
                let theta = j as f64 / ns as f64 * TAU;
                let ct = theta.cos();
                let st = theta.sin();
                let pos = Point3r::new(cx + ro * sp * ct, y, cz + ro * sp * st);
                let n = Vector3r::new(sp * ct, cp, sp * st);
                mesh.add_vertex(pos, n)
            })
            .collect();
        outer_rings.push(row);
    }

    // Build inner ring VertexId arrays
    let mut inner_rings: Vec<Vec<crate::domain::core::index::VertexId>> =
        Vec::with_capacity(nk - 1);
    for k in 0..nk - 1 {
        let phi = (k + 1) as f64 / nk as f64 * PI;
        let sp = phi.sin();
        let cp = phi.cos();
        let y = cy + ri * cp;
        let row: Vec<_> = (0..ns)
            .map(|j| {
                let theta = j as f64 / ns as f64 * TAU;
                let ct = theta.cos();
                let st = theta.sin();
                let pos = Point3r::new(cx + ri * sp * ct, y, cz + ri * sp * st);
                let n = Vector3r::new(-sp * ct, -cp, -sp * st);
                mesh.add_vertex(pos, n)
            })
            .collect();
        inner_rings.push(row);
    }

    // Outer sphere lateral (outward normals)
    for k in 0..outer_rings.len() - 1 {
        for j in 0..ns {
            let j1 = (j + 1) % ns;
            let vu0 = outer_rings[k][j];
            let vu1 = outer_rings[k][j1];
            let vl0 = outer_rings[k + 1][j];
            let vl1 = outer_rings[k + 1][j1];
            // CCW from outside: upper[j]->lower[j]->lower[j+1], upper[j]->lower[j+1]->upper[j+1]
            mesh.add_face_with_region(vu0, vl0, vl1, region);
            mesh.add_face_with_region(vu0, vl1, vu1, region);
        }
    }

    // Inner sphere lateral (reversed winding for inward surface normals)
    for k in 0..inner_rings.len() - 1 {
        for j in 0..ns {
            let j1 = (j + 1) % ns;
            let vu0 = inner_rings[k][j];
            let vu1 = inner_rings[k][j1];
            let vl0 = inner_rings[k + 1][j];
            let vl1 = inner_rings[k + 1][j1];
            // Reversed winding (normals point inward = toward cavity centre):
            mesh.add_face_with_region(vu0, vu1, vl1, region);
            mesh.add_face_with_region(vu0, vl1, vl0, region);
        }
    }

    // North polar cap (connects outer_rings[0] and inner_rings[0])
    // Outer lateral at north provides: outer[0][j1]->outer[0][j] (decreasing j direction)
    // Inner lateral at north provides: inner[0][j]->inner[0][j1] (increasing j direction)
    // Cap must provide the opposite directions for manifold topology.
    {
        for j in 0..ns {
            let j1 = (j + 1) % ns;
            let oi = outer_rings[0][j];
            let oi1 = outer_rings[0][j1];
            let ii = inner_rings[0][j];
            let ii1 = inner_rings[0][j1];
            // Face 1: (oi, oi1, ii1) provides outer[j]->outer[j1] (opposite of lateral)
            //         and inner[j1]->inner[j] via face 2 below... wait let me recalculate
            // Need: oi->oi1 (opposite of lateral oi1->oi) and ii1->ii (opposite of lateral ii->ii1)
            // Face 1 (oi, ii1, ii): edges oi->ii1, ii1->ii (correct for inner), ii->oi
            // Face 2 (oi, oi1, ii1): edges oi->oi1 (correct for outer), oi1->ii1, ii1->oi
            mesh.add_face_with_region(oi, ii1, ii, region);
            mesh.add_face_with_region(oi, oi1, ii1, region);
        }
    }

    // South polar cap (connects outer_rings[last] and inner_rings[last])
    // Outer lateral at south provides: outer[last][j]->outer[last][j1] (increasing j direction)
    // Inner lateral at south provides: inner[last][j1]->inner[last][j] (decreasing j direction)
    // Cap must provide opposite directions.
    {
        let outer_last_k = outer_rings.len() - 1;
        let inner_last_k = inner_rings.len() - 1;
        for j in 0..ns {
            let j1 = (j + 1) % ns;
            let oi = outer_rings[outer_last_k][j];
            let oi1 = outer_rings[outer_last_k][j1];
            let ii = inner_rings[inner_last_k][j];
            let ii1 = inner_rings[inner_last_k][j1];
            // Need: oi1->oi (opposite of lateral oi->oi1) and ii->ii1 (opposite of lateral ii1->ii)
            // Face 1 (ii, ii1, oi1): edges ii->ii1 (correct), ii1->oi1, oi1->ii
            // Face 2 (ii, oi1, oi): edges ii->oi1, oi1->oi (correct), oi->ii
            mesh.add_face_with_region(ii, ii1, oi1, region);
            mesh.add_face_with_region(ii, oi1, oi, region);
        }
    }

    // All sections (outer lateral, inner lateral, north/south polar caps) are built
    // with consistent inward winding. Flip all faces to obtain outward-pointing normals
    // and positive signed volume.
    mesh.flip_faces();

    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::application::watertight::check::check_watertight;
    use crate::infrastructure::storage::edge_store::EdgeStore;
    use std::f64::consts::PI;

    #[test]
    fn spherical_shell_is_closed_and_oriented() {
        let mesh = SphericalShell::default().build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.is_closed, "shell must be closed");
        assert!(
            report.orientation_consistent,
            "shell must be consistently oriented"
        );
        // Genus 1 (through-bore at poles) → χ = 0
        assert_eq!(
            report.euler_characteristic,
            Some(0),
            "spherical shell Euler characteristic must be 0 (genus 1)"
        );
        assert!(report.is_watertight, "shell passes watertight check");
    }

    #[test]
    fn spherical_shell_volume_positive_and_approximately_correct() {
        let (ri, ro) = (0.9_f64, 1.0_f64);
        let mesh = SphericalShell {
            outer_radius: ro,
            inner_radius: ri,
            segments: 64,
            stacks: 32,
            ..SphericalShell::default()
        }
        .build()
        .unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.signed_volume > 0.0);
        let expected = 4.0 / 3.0 * PI * (ro * ro * ro - ri * ri * ri);
        let error = (report.signed_volume - expected).abs() / expected;
        assert!(
            error < 0.01,
            "volume error {:.4}% should be < 1%",
            error * 100.0
        );
    }

    #[test]
    fn spherical_shell_rejects_invalid_params() {
        assert!(SphericalShell {
            inner_radius: 0.0,
            ..SphericalShell::default()
        }
        .build()
        .is_err());
        assert!(SphericalShell {
            outer_radius: 0.8,
            inner_radius: 0.9,
            ..SphericalShell::default()
        }
        .build()
        .is_err());
        assert!(SphericalShell {
            segments: 2,
            ..SphericalShell::default()
        }
        .build()
        .is_err());
        assert!(SphericalShell {
            stacks: 2,
            ..SphericalShell::default()
        }
        .build()
        .is_err());
    }
}
