//! Antiprism primitive — n-gon top and bottom rotated by π/n, connected by triangles.

use std::f64::consts::TAU;

use super::{PrimitiveError, PrimitiveMesh};
use crate::domain::core::index::RegionId;
use crate::domain::core::scalar::{Point3r, Vector3r};
use crate::domain::geometry::normal::triangle_normal;
use crate::domain::mesh::IndexedMesh;

/// Builds a right antiprism with two parallel regular n-gon faces.
///
/// An antiprism is like a prism but with the top face rotated by `π / sides`
/// relative to the bottom face, and the lateral strip filled with `2·sides`
/// equilateral triangles rather than rectangles.  The triangular `sides=3`
/// case is an octahedron.
///
/// The antiprism axis is +Y; the base is centred at `base_center` and the top
/// is at `base_center + (0, height, 0)`.
///
/// ## Topology
///
/// - V = 2·sides, E = 4·sides, F = 2·sides + 2
/// - `V − E + F = 2`  (χ = 2, genus 0)
///
/// ## Output
///
/// - `RegionId(1)` on all faces
/// - For `sides = 3`: degenerate octahedron (same geometry as the Octahedron
///   primitive when `height = r·√2`)
#[derive(Clone, Debug)]
pub struct Antiprism {
    /// Centre of the bottom polygon.
    pub base_center: Point3r,
    /// Circumradius of both polygons [mm].
    pub base_radius: f64,
    /// Height of the antiprism [mm].
    pub height: f64,
    /// Number of sides on each face (≥ 3).
    pub sides: usize,
}

impl Default for Antiprism {
    fn default() -> Self {
        Self {
            base_center: Point3r::origin(),
            base_radius: 1.0,
            height: 1.0,
            sides: 6,
        }
    }
}

impl PrimitiveMesh for Antiprism {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build(self)
    }
}

fn build(ap: &Antiprism) -> Result<IndexedMesh, PrimitiveError> {
    if ap.base_radius <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "base_radius must be > 0, got {}",
            ap.base_radius
        )));
    }
    if ap.height <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "height must be > 0, got {}",
            ap.height
        )));
    }
    if ap.sides < 3 {
        return Err(PrimitiveError::TooFewSegments(ap.sides));
    }

    let region = RegionId::new(1);
    let mut mesh = IndexedMesh::new();

    let r = ap.base_radius;
    let bx = ap.base_center.x;
    let by = ap.base_center.y;
    let bz = ap.base_center.z;
    let ns = ap.sides;

    // Bottom polygon vertices (y = by)
    let bot: Vec<Point3r> = (0..ns)
        .map(|i| {
            let angle = i as f64 / ns as f64 * TAU;
            Point3r::new(bx + r * angle.cos(), by, bz + r * angle.sin())
        })
        .collect();

    // Top polygon vertices (y = by + height, rotated by π/n)
    let top: Vec<Point3r> = (0..ns)
        .map(|i| {
            let angle = (i as f64 + 0.5) / ns as f64 * TAU;
            Point3r::new(bx + r * angle.cos(), by + ap.height, bz + r * angle.sin())
        })
        .collect();

    // ── Lateral strip (2·sides triangles) ─────────────────────────────────────
    // For each i:
    //   "Up" triangle:   top[i], bot[i], bot[(i+1)%n]
    //   "Down" triangle: top[i], bot[(i+1)%n], top[(i+1)%n]
    for i in 0..ns {
        let j = (i + 1) % ns;
        {
            let n = triangle_normal(&top[i], &bot[i], &bot[j]).unwrap_or(Vector3r::zeros());
            let vt = mesh.add_vertex(top[i], n);
            let vb0 = mesh.add_vertex(bot[i], n);
            let vb1 = mesh.add_vertex(bot[j], n);
            mesh.add_face_with_region(vt, vb0, vb1, region);
        }
        {
            let n = triangle_normal(&top[i], &bot[j], &top[j]).unwrap_or(Vector3r::zeros());
            let vt0 = mesh.add_vertex(top[i], n);
            let vb = mesh.add_vertex(bot[j], n);
            let vt1 = mesh.add_vertex(top[j], n);
            mesh.add_face_with_region(vt0, vb, vt1, region);
        }
    }

    // ── Bottom cap (normal −Y) ─────────────────────────────────────────────────
    {
        let n_down = -Vector3r::y();
        let vc = mesh.add_vertex(Point3r::new(bx, by, bz), n_down);
        for i in 0..ns {
            let j = (i + 1) % ns;
            let v0 = mesh.add_vertex(bot[i], n_down);
            let v1 = mesh.add_vertex(bot[j], n_down);
            mesh.add_face_with_region(vc, v1, v0, region);
        }
    }

    // ── Top cap (normal +Y) ────────────────────────────────────────────────────
    {
        let n_up = Vector3r::y();
        let vc = mesh.add_vertex(Point3r::new(bx, by + ap.height, bz), n_up);
        for i in 0..ns {
            let j = (i + 1) % ns;
            let v0 = mesh.add_vertex(top[i], n_up);
            let v1 = mesh.add_vertex(top[j], n_up);
            mesh.add_face_with_region(vc, v0, v1, region);
        }
    }

    // Lateral strip winding produces inward normals with this vertex ordering.
    // Flip all faces to restore outward-pointing normals.
    mesh.flip_faces();

    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::application::watertight::check::check_watertight;
    use crate::infrastructure::storage::edge_store::EdgeStore;

    #[test]
    fn antiprism_is_watertight() {
        let mesh = Antiprism::default().build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.is_watertight, "antiprism must be watertight");
        assert_eq!(report.euler_characteristic, Some(2));
    }

    #[test]
    fn antiprism_volume_positive() {
        let mesh = Antiprism {
            sides: 8,
            height: 2.0,
            base_radius: 1.5,
            ..Antiprism::default()
        }
        .build()
        .unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.is_watertight);
        assert!(report.signed_volume > 0.0);
    }

    #[test]
    fn antiprism_rejects_invalid_params() {
        assert!(Antiprism {
            base_radius: 0.0,
            ..Antiprism::default()
        }
        .build()
        .is_err());
        assert!(Antiprism {
            height: 0.0,
            ..Antiprism::default()
        }
        .build()
        .is_err());
        assert!(Antiprism {
            sides: 2,
            ..Antiprism::default()
        }
        .build()
        .is_err());
    }
}
