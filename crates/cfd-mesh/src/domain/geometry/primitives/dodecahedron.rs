//! Regular dodecahedron primitive — 12 pentagonal faces.

use super::{PrimitiveError, PrimitiveMesh};
use crate::domain::core::index::RegionId;
use crate::domain::core::scalar::{Point3r, Vector3r};
use crate::domain::geometry::normal::triangle_normal;
use crate::domain::mesh::IndexedMesh;

/// Builds a regular dodecahedron inscribed in a sphere of the given `radius`.
///
/// ## Geometry
///
/// The dodecahedron has 20 vertices, 30 edges, and 12 regular pentagonal
/// faces.  Vertices are taken from:
/// - All 8 cube vertices `(±1, ±1, ±1)`, and
/// - All 12 even permutations of `(0, ±1/φ, ±φ)` where `φ = (1+√5)/2`.
///
/// All 20 raw vertices lie on a sphere of radius `√3`.  They are scaled to
/// the requested circumsphere radius before building.
///
/// ## Topology
///
/// - V = 20, E = 30, F = 12 pentagons → 60 triangles (fan-triangulated)
/// - `V − E + F_orig = 2`  (χ = 2, genus 0)
/// - `signed_volume > 0` (outward CCW winding)
///
/// ## Output
///
/// - `RegionId(1)` on all faces
/// - `signed_volume = (1/4)(15 + 7√5) a³` where `a` = edge length
#[derive(Clone, Debug)]
pub struct Dodecahedron {
    /// Circumsphere radius [mm].
    pub radius: f64,
    /// Centre position.
    pub center: Point3r,
}

impl Default for Dodecahedron {
    fn default() -> Self {
        Self {
            radius: 1.0,
            center: Point3r::origin(),
        }
    }
}

impl PrimitiveMesh for Dodecahedron {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build(self)
    }
}

/// Golden ratio φ = (1 + √5) / 2.
const PHI: f64 = 1.618_033_988_749_895;

/// 20 raw vertices (circumradius = √3, centred at origin).
///
/// Index layout:
/// - 0–7   cube corners  (±1, ±1, ±1)
/// - 8–11  (0, ±1/φ, ±φ) permutation
/// - 12–15 (±1/φ, ±φ, 0) permutation
/// - 16–19 (±φ, 0, ±1/φ) permutation
fn raw_vertices() -> [[f64; 3]; 20] {
    let p = PHI; // φ
    let ip = 1.0 / PHI; // 1/φ
    [
        // cube corners
        [1.0, 1.0, 1.0],    // 0
        [1.0, 1.0, -1.0],   // 1
        [1.0, -1.0, 1.0],   // 2
        [1.0, -1.0, -1.0],  // 3
        [-1.0, 1.0, 1.0],   // 4
        [-1.0, 1.0, -1.0],  // 5
        [-1.0, -1.0, 1.0],  // 6
        [-1.0, -1.0, -1.0], // 7
        // (0, ±1/φ, ±φ) family
        [0.0, ip, p],   // 8
        [0.0, ip, -p],  // 9
        [0.0, -ip, p],  // 10
        [0.0, -ip, -p], // 11
        // (±1/φ, ±φ, 0) family
        [ip, p, 0.0],   // 12
        [ip, -p, 0.0],  // 13
        [-ip, p, 0.0],  // 14
        [-ip, -p, 0.0], // 15
        // (±φ, 0, ±1/φ) family
        [p, 0.0, ip],   // 16
        [p, 0.0, -ip],  // 17
        [-p, 0.0, ip],  // 18
        [-p, 0.0, -ip], // 19
    ]
}

/// 12 pentagonal faces in CCW outward winding order.
///
/// Derivation: standard (CW-when-viewed-outside) dodecahedron face rings are
/// reversed to obtain CCW outward winding consistent with the positive
/// signed-volume convention used throughout cfd-mesh.
///
/// Each vertex appears in exactly 3 faces (3 pentagons meet at every vertex
/// of a dodecahedron), and every directed edge appears once in each direction
/// across adjacent faces — confirmed by edge-pair enumeration.
const DODECAHEDRON_FACES: [[usize; 5]; 12] = [
    [0, 8, 10, 2, 16],  // face  1 (+z cluster top)
    [0, 12, 14, 4, 8],  // face  2 (+y front)
    [0, 16, 17, 1, 12], // face  3 (+x top)
    [8, 4, 18, 6, 10],  // face  4 (-x front)
    [4, 14, 5, 19, 18], // face  5 (-x back)
    [14, 12, 1, 9, 5],  // face  6 (+y back)
    [16, 2, 13, 3, 17], // face  7 (+x bottom)
    [2, 10, 6, 15, 13], // face  8 (-x bottom)
    [6, 18, 19, 7, 15], // face  9 (-z cluster bottom)
    [1, 17, 3, 11, 9],  // face 10 (+x back)
    [3, 13, 15, 7, 11], // face 11 (-y bottom)
    [5, 9, 11, 7, 19],  // face 12 (-y back)
];

fn build(d: &Dodecahedron) -> Result<IndexedMesh, PrimitiveError> {
    if d.radius <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "radius must be > 0, got {}",
            d.radius
        )));
    }

    let region = RegionId::new(1);
    let mut mesh = IndexedMesh::new();

    // All 20 raw vertices lie on sphere of radius √3; scale to desired radius.
    let scale = d.radius / 3.0_f64.sqrt();
    let raw = raw_vertices();
    let cx = d.center.x;
    let cy = d.center.y;
    let cz = d.center.z;

    let verts: Vec<Point3r> = raw
        .iter()
        .map(|v| Point3r::new(cx + v[0] * scale, cy + v[1] * scale, cz + v[2] * scale))
        .collect();

    for face_indices in &DODECAHEDRON_FACES {
        // Face centroid — used to orient per-triangle normals outward.
        let centroid_coords: Vector3r = face_indices
            .iter()
            .map(|&i| verts[i].coords)
            .fold(Vector3r::zeros(), |acc, v| acc + v)
            / 5.0;
        let centroid = Point3r::from(centroid_coords);

        // Fan-triangulate from the first vertex of the face pentagon.
        let p0 = &verts[face_indices[0]];
        for k in 1..4 {
            let p1 = &verts[face_indices[k]];
            let p2 = &verts[face_indices[k + 1]];

            let n = match triangle_normal(p0, p1, p2) {
                Some(n) => {
                    // Flip if pointing toward the centre rather than outward.
                    if n.dot(&centroid.coords) < 0.0 {
                        -n
                    } else {
                        n
                    }
                }
                None => Vector3r::zeros(),
            };

            let vi0 = mesh.add_vertex(*p0, n);
            let vi1 = mesh.add_vertex(*p1, n);
            let vi2 = mesh.add_vertex(*p2, n);
            mesh.add_face_with_region(vi0, vi1, vi2, region);
        }
    }

    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::application::watertight::check::check_watertight;
    use crate::infrastructure::storage::edge_store::EdgeStore;

    #[test]
    fn dodecahedron_is_watertight() {
        let mesh = Dodecahedron::default().build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(
            report.is_watertight,
            "dodecahedron must be watertight: {:?}",
            report
        );
        assert_eq!(report.euler_characteristic, Some(2));
    }

    #[test]
    fn dodecahedron_volume_positive() {
        let mesh = Dodecahedron {
            radius: 2.0,
            ..Dodecahedron::default()
        }
        .build()
        .unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.signed_volume > 0.0, "volume must be positive");
    }

    #[test]
    fn dodecahedron_invalid_radius() {
        assert!(Dodecahedron {
            radius: 0.0,
            ..Dodecahedron::default()
        }
        .build()
        .is_err());
        assert!(Dodecahedron {
            radius: -1.0,
            ..Dodecahedron::default()
        }
        .build()
        .is_err());
    }
}
