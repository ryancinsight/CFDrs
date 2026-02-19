//! Regular tetrahedron primitive.

use crate::core::index::RegionId;
use crate::core::scalar::{Point3r, Vector3r};
use crate::geometry::normal::triangle_normal;
use crate::mesh::IndexedMesh;
use super::{PrimitiveMesh, PrimitiveError};

/// Builds a regular tetrahedron inscribed in a sphere of the given `radius`.
///
/// ## Geometry
///
/// The four vertices are placed at alternating corners of a cube with
/// half-side `r`:
///
/// ```text
///   A = ( r,  r,  r)
///   B = (-r, -r,  r)
///   C = (-r,  r, -r)
///   D = ( r, -r, -r)
/// ```
///
/// All six edge lengths equal `2r√2`.  The solid is centred at the origin.
///
/// ## Output
///
/// - 4 vertices, 4 faces
/// - `RegionId(1)` on all faces
/// - `signed_volume = (8/3) r³ > 0`
#[derive(Clone, Debug)]
pub struct Tetrahedron {
    /// Circumsphere radius [mm].
    pub radius: f64,
}

impl Default for Tetrahedron {
    fn default() -> Self { Self { radius: 1.0 } }
}

impl PrimitiveMesh for Tetrahedron {
    fn build(&self) -> Result<crate::mesh::IndexedMesh, PrimitiveError> {
        build(self)
    }
}

fn build(t: &Tetrahedron) -> Result<IndexedMesh, PrimitiveError> {
    let r = t.radius;
    if r <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "radius must be > 0, got {r}"
        )));
    }

    let region = RegionId::new(1);
    let mut mesh = IndexedMesh::new();

    // Canonical vertices: alternating corners of a cube of half-side r.
    // Centroid = origin.
    let a = Point3r::new( r,  r,  r);
    let b = Point3r::new(-r, -r,  r);
    let c = Point3r::new(-r,  r, -r);
    let d = Point3r::new( r, -r, -r);

    // Outward-facing winding (outward CCW, signed_volume > 0):
    //  face ACB, face ABD, face ADC, face BCD
    // Each face normal is computed from the winding order via cross-product,
    // which is then used as the vertex normal.
    let faces: [(&Point3r, &Point3r, &Point3r); 4] = [
        (&a, &c, &b), // face ACB
        (&a, &b, &d), // face ABD
        (&a, &d, &c), // face ADC
        (&b, &c, &d), // face BCD
    ];

    for (p0, p1, p2) in faces {
        let n = triangle_normal(p0, p1, p2)
            .unwrap_or(Vector3r::zeros());
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
    fn tetrahedron_is_watertight() {
        let mesh = Tetrahedron { radius: 1.0 }.build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.is_watertight, "tetrahedron must be watertight");
        assert_eq!(report.euler_characteristic, Some(2));
    }

    #[test]
    fn tetrahedron_volume_positive() {
        let r = 1.0_f64;
        let mesh = Tetrahedron { radius: r }.build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.signed_volume > 0.0, "outward normals → positive volume");
    }

    #[test]
    fn tetrahedron_invalid_radius() {
        let result = Tetrahedron { radius: -1.0 }.build();
        assert!(result.is_err());
    }
}
