//! Regular octahedron primitive.

use crate::core::index::RegionId;
use crate::core::scalar::{Point3r, Vector3r};
use crate::geometry::normal::triangle_normal;
use crate::mesh::IndexedMesh;
use super::{PrimitiveMesh, PrimitiveError};

/// Builds a regular octahedron inscribed in a sphere of the given `radius`.
///
/// ## Geometry
///
/// Six vertices are placed at ±`radius` along each Cartesian axis:
///
/// ```text
///   +Y apex  (0,  R, 0)
///   −Y apex  (0, −R, 0)
///   +X       ( R, 0, 0)
///   −X       (−R, 0, 0)
///   +Z       (0, 0,  R)
///   −Z       (0, 0, −R)
/// ```
///
/// All edges have length `R√2`. The solid is centred at `center`.
///
/// ## Topology
///
/// - 6 vertices, 12 edges, 8 equilateral triangular faces
/// - `V − E + F = 6 − 12 + 8 = 2`  (genus 0, χ = 2)
/// - `signed_volume = (4/3) R³`
///
/// ## Output
///
/// - `RegionId(1)` on all faces
/// - `signed_volume > 0` (outward CCW winding)
#[derive(Clone, Debug)]
pub struct Octahedron {
    /// Circumsphere radius [mm].
    pub radius: f64,
    /// Centre position.
    pub center: Point3r,
}

impl Default for Octahedron {
    fn default() -> Self {
        Self { radius: 1.0, center: Point3r::origin() }
    }
}

impl PrimitiveMesh for Octahedron {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build(self)
    }
}

fn build(o: &Octahedron) -> Result<IndexedMesh, PrimitiveError> {
    if o.radius <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "radius must be > 0, got {}", o.radius
        )));
    }

    let region = RegionId::new(1);
    let mut mesh = IndexedMesh::new();
    let r = o.radius;
    let cx = o.center.x;
    let cy = o.center.y;
    let cz = o.center.z;

    // Six axis vertices.
    let px = Point3r::new(cx + r, cy,     cz    );
    let nx = Point3r::new(cx - r, cy,     cz    );
    let py = Point3r::new(cx,     cy + r, cz    );
    let ny = Point3r::new(cx,     cy - r, cz    );
    let pz = Point3r::new(cx,     cy,     cz + r);
    let nz = Point3r::new(cx,     cy,     cz - r);

    // Eight faces with outward CCW winding.
    // Each normal is the unit centroid direction (verified analytically).
    //
    // Top cap (py apex) — ring order CCW from above: pz, px, nz, nx
    //   Face (py, pz, px):  normal (+1,+1,+1)/√3
    //   Face (py, px, nz):  normal (+1,+1,−1)/√3
    //   Face (py, nz, nx):  normal (−1,+1,−1)/√3
    //   Face (py, nx, pz):  normal (−1,+1,+1)/√3
    //
    // Bottom cap (ny apex) — reversed ring order:
    //   Face (ny, px, pz):  normal (+1,−1,+1)/√3
    //   Face (ny, nz, px):  normal (+1,−1,−1)/√3
    //   Face (ny, nx, nz):  normal (−1,−1,−1)/√3
    //   Face (ny, pz, nx):  normal (−1,−1,+1)/√3
    let face_verts: [(&Point3r, &Point3r, &Point3r); 8] = [
        (&py, &pz, &px),
        (&py, &px, &nz),
        (&py, &nz, &nx),
        (&py, &nx, &pz),
        (&ny, &px, &pz),
        (&ny, &nz, &px),
        (&ny, &nx, &nz),
        (&ny, &pz, &nx),
    ];

    for (p0, p1, p2) in face_verts {
        let n = triangle_normal(p0, p1, p2).unwrap_or(Vector3r::zeros());
        let v0 = mesh.add_vertex(*p0, n);
        let v1 = mesh.add_vertex(*p1, n);
        let v2 = mesh.add_vertex(*p2, n);
        mesh.add_face_with_region(v0, v1, v2, region);
    }

    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::storage::edge_store::EdgeStore;
    use crate::watertight::check::check_watertight;
    use approx::assert_relative_eq;

    #[test]
    fn octahedron_is_watertight() {
        let mesh = Octahedron::default().build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.is_watertight, "octahedron must be watertight");
        assert_eq!(report.euler_characteristic, Some(2));
    }

    #[test]
    fn octahedron_volume_positive_and_exact() {
        let r = 1.5_f64;
        let mesh = Octahedron { radius: r, ..Octahedron::default() }.build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.signed_volume > 0.0, "outward normals → positive volume");
        // V = (4/3) R³  (exact for octahedron, no discretisation error)
        let expected = 4.0 / 3.0 * r * r * r;
        assert_relative_eq!(report.signed_volume, expected, epsilon = 1e-10);
    }

    #[test]
    fn octahedron_invalid_radius() {
        assert!(Octahedron { radius: 0.0, ..Octahedron::default() }.build().is_err());
        assert!(Octahedron { radius: -1.0, ..Octahedron::default() }.build().is_err());
    }
}
