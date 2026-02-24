//! Regular icosahedron primitive.

use super::{PrimitiveError, PrimitiveMesh};
use crate::domain::core::index::RegionId;
use crate::domain::core::scalar::{Point3r, Vector3r};
use crate::domain::geometry::normal::triangle_normal;
use crate::domain::mesh::IndexedMesh;

/// Builds a regular icosahedron inscribed in a sphere of the given `radius`.
///
/// ## Geometry
///
/// Twelve vertices are constructed from all even permutations of `(0, ±1, ±φ)`,
/// where `φ = (1 + √5) / 2` is the golden ratio, then scaled to the desired
/// circumradius.  The 20 equilateral triangular faces use outward CCW winding
/// (verified analytically via centroid dot-product).
///
/// ## Topology
///
/// - 12 vertices, 30 edges, 20 faces
/// - `V − E + F = 12 − 30 + 20 = 2`  (genus 0, χ = 2)
/// - No pole singularities; all triangles are equilateral
///
/// ## Uses
///
/// - Base mesh for [`GeodesicSphere`] subdivision
/// - More isotropic sphere approximation than [`UvSphere`] (no polar crowding)
/// - CSG robustness testing with 5-valent vertices
///
/// ## Output
///
/// - `RegionId(1)` on all faces
/// - `signed_volume = (5(3+√5)/12) · a³`  where `a = 4R / √(10+2√5)`
#[derive(Clone, Debug)]
pub struct Icosahedron {
    /// Circumsphere radius [mm].
    pub radius: f64,
    /// Centre position.
    pub center: Point3r,
}

impl Default for Icosahedron {
    fn default() -> Self {
        Self {
            radius: 1.0,
            center: Point3r::origin(),
        }
    }
}

impl PrimitiveMesh for Icosahedron {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build(self)
    }
}

/// Compute the 12 scaled icosahedron vertices for the given circumradius,
/// centered at the origin (caller offsets by center).
pub(crate) fn icosahedron_vertices(radius: f64) -> [[f64; 3]; 12] {
    let phi = (1.0_f64 + 5.0_f64.sqrt()) / 2.0;
    // Circumradius of the unit icosahedron (vertices at (0,±1,±φ)):
    let unit_r = (1.0_f64 + phi * phi).sqrt();
    let scale = radius / unit_r;

    [
        [0.0, scale, phi * scale],   // v0
        [0.0, scale, -phi * scale],  // v1
        [0.0, -scale, phi * scale],  // v2
        [0.0, -scale, -phi * scale], // v3
        [scale, phi * scale, 0.0],   // v4
        [scale, -phi * scale, 0.0],  // v5
        [-scale, phi * scale, 0.0],  // v6
        [-scale, -phi * scale, 0.0], // v7
        [phi * scale, 0.0, scale],   // v8
        [phi * scale, 0.0, -scale],  // v9
        [-phi * scale, 0.0, scale],  // v10
        [-phi * scale, 0.0, -scale], // v11
    ]
}

/// The 20 icosahedron face triplets with outward CCW winding.
///
/// Each triple is an index into `icosahedron_vertices()`.
/// Winding verified analytically: `centroid · normal > 0` for all 20 faces.
pub(crate) const ICOSAHEDRON_FACES: [[usize; 3]; 20] = [
    // Top cap (around v0)
    [0, 2, 8],  //  0
    [0, 8, 4],  //  1
    [0, 4, 6],  //  2
    [0, 6, 10], //  3
    [0, 10, 2], //  4
    // Middle upper (connecting top equator to bottom vertices)
    [2, 5, 8],   //  5
    [8, 9, 4],   //  6
    [4, 1, 6],   //  7
    [6, 11, 10], //  8
    [10, 7, 2],  //  9
    // Middle lower (connecting bottom equator to top vertices)
    [5, 9, 8],   // 10
    [7, 5, 2],   // 11
    [11, 7, 10], // 12
    [1, 11, 6],  // 13
    [9, 1, 4],   // 14
    // Bottom cap (around v3)
    [3, 9, 5],  // 15
    [3, 5, 7],  // 16
    [3, 7, 11], // 17
    [3, 11, 1], // 18
    [3, 1, 9],  // 19
];

fn build(ico: &Icosahedron) -> Result<IndexedMesh, PrimitiveError> {
    if ico.radius <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "radius must be > 0, got {}",
            ico.radius
        )));
    }

    let region = RegionId::new(1);
    let mut mesh = IndexedMesh::new();

    let raw = icosahedron_vertices(ico.radius);
    let cx = ico.center.x;
    let cy = ico.center.y;
    let cz = ico.center.z;

    let pts: [Point3r; 12] =
        std::array::from_fn(|i| Point3r::new(cx + raw[i][0], cy + raw[i][1], cz + raw[i][2]));

    for [a, b, c] in ICOSAHEDRON_FACES {
        let n = triangle_normal(&pts[a], &pts[b], &pts[c]).unwrap_or(Vector3r::zeros());
        let va = mesh.add_vertex(pts[a], n);
        let vb = mesh.add_vertex(pts[b], n);
        let vc = mesh.add_vertex(pts[c], n);
        mesh.add_face_with_region(va, vb, vc, region);
    }

    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::application::watertight::check::check_watertight;
    use crate::infrastructure::storage::edge_store::EdgeStore;

    #[test]
    fn icosahedron_is_watertight() {
        let mesh = Icosahedron::default().build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.is_watertight, "icosahedron must be watertight");
        assert_eq!(report.euler_characteristic, Some(2));
    }

    #[test]
    fn icosahedron_volume_positive_and_approximately_correct() {
        let r = 1.0_f64;
        let mesh = Icosahedron {
            radius: r,
            ..Icosahedron::default()
        }
        .build()
        .unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(
            report.signed_volume > 0.0,
            "outward normals → positive volume"
        );

        // Analytic volume: V = (5/12)(3+√5) · a³,  a = 4R / √(10+2√5)
        let sqrt5 = 5.0_f64.sqrt();
        let a = 4.0 * r / (10.0 + 2.0 * sqrt5).sqrt();
        let expected = (5.0 / 12.0) * (3.0 + sqrt5) * a.powi(3);
        let error = (report.signed_volume - expected).abs() / expected;
        assert!(
            error < 1e-9,
            "volume error {:.2e} should be < 1e-9 (exact geometry)",
            error
        );
    }

    #[test]
    fn icosahedron_invalid_radius() {
        assert!(Icosahedron {
            radius: 0.0,
            ..Icosahedron::default()
        }
        .build()
        .is_err());
        assert!(Icosahedron {
            radius: -1.0,
            ..Icosahedron::default()
        }
        .build()
        .is_err());
    }

    #[test]
    fn icosahedron_centered_at_offset() {
        let mesh = Icosahedron {
            radius: 1.0,
            center: Point3r::new(10.0, 5.0, -3.0),
        }
        .build()
        .unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.is_watertight);
        assert!(report.signed_volume > 0.0);
    }
}
