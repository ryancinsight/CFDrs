//! Cuboctahedron primitive — Archimedean solid with 8 triangles + 6 squares.

use super::{PrimitiveError, PrimitiveMesh};
use crate::domain::core::index::RegionId;
use crate::domain::core::scalar::{Point3r, Vector3r};
use crate::domain::geometry::normal::triangle_normal;
use crate::domain::mesh::IndexedMesh;

/// Builds a regular cuboctahedron inscribed in a sphere of the given `radius`.
///
/// ## Geometry
///
/// The cuboctahedron has 12 vertices at all permutations of `(0, ±1, ±1)`,
/// scaled to the desired circumradius.  It has:
/// - 8 equilateral triangular faces
/// - 6 square faces (each fan-triangulated into 2 triangles)
///
/// It is the intersection of a cube and an octahedron of equal edge length,
/// and corresponds to the FCC / BCC Wigner-Seitz unit cell.
///
/// ## Topology
///
/// - V = 12, E = 24, F = 14  →  `V − E + F = 2`  (χ = 2, genus 0)
/// - After square fan-triangulation: 20 triangular faces.
///
/// ## Output
///
/// - `RegionId(1)` on all faces
/// - `signed_volume = (5√2 / 3) · a³` where `a` = edge length
#[derive(Clone, Debug)]
pub struct Cuboctahedron {
    /// Circumsphere radius [mm].
    pub radius: f64,
    /// Centre position.
    pub center: Point3r,
}

impl Default for Cuboctahedron {
    fn default() -> Self {
        Self {
            radius: 1.0,
            center: Point3r::origin(),
        }
    }
}

impl PrimitiveMesh for Cuboctahedron {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build(self)
    }
}

/// 12 vertices: all permutations of (0, ±1, ±1), circumradius = √2.
fn raw_vertices() -> [[f64; 3]; 12] {
    [
        // (0, ±1, ±1) permutations
        [0.0, 1.0, 1.0],   // 0
        [0.0, 1.0, -1.0],  // 1
        [0.0, -1.0, 1.0],  // 2
        [0.0, -1.0, -1.0], // 3
        // (±1, 0, ±1) permutations
        [1.0, 0.0, 1.0],   // 4
        [1.0, 0.0, -1.0],  // 5
        [-1.0, 0.0, 1.0],  // 6
        [-1.0, 0.0, -1.0], // 7
        // (±1, ±1, 0) permutations
        [1.0, 1.0, 0.0],   // 8
        [1.0, -1.0, 0.0],  // 9
        [-1.0, 1.0, 0.0],  // 10
        [-1.0, -1.0, 0.0], // 11
    ]
}

/// 8 triangular faces — one per sign-octant (sx,sy,sz).
/// Each face connects (0,sy,sz), (sx,0,sz), (sx,sy,0) for the given sign triple.
/// Winding checked via centroid dot-product in `build`.
///
/// Vertex key:
///  0=(0,+1,+1), 1=(0,+1,−1), 2=(0,−1,+1), 3=(0,−1,−1)
///  4=(+1,0,+1), 5=(+1,0,−1), 6=(−1,0,+1), 7=(−1,0,−1)
///  8=(+1,+1,0), 9=(+1,−1,0),10=(−1,+1,0),11=(−1,−1,0)
const TRI_FACES: [[usize; 3]; 8] = [
    [0, 4, 8],  // (+,+,+): (0,1,1),(1,0,1),(1,1,0)
    [1, 5, 8],  // (+,+,−): (0,1,−1),(1,0,−1),(1,1,0)
    [2, 4, 9],  // (+,−,+): (0,−1,1),(1,0,1),(1,−1,0)
    [3, 5, 9],  // (+,−,−): (0,−1,−1),(1,0,−1),(1,−1,0)
    [0, 6, 10], // (−,+,+): (0,1,1),(−1,0,1),(−1,1,0)
    [1, 7, 10], // (−,+,−): (0,1,−1),(−1,0,−1),(−1,1,0)
    [2, 6, 11], // (−,−,+): (0,−1,1),(−1,0,1),(−1,−1,0)
    [3, 7, 11], // (−,−,−): (0,−1,−1),(−1,0,−1),(−1,−1,0)
];

/// 6 square faces as vertex quads (fan-triangulate to 2 triangles each).
/// Each quad lists vertices CCW as seen from outside (outward normal direction).
///
/// Vertex positions:
///  0=(0,1,1), 1=(0,1,-1), 2=(0,-1,1), 3=(0,-1,-1)
///  4=(1,0,1), 5=(1,0,-1), 6=(-1,0,1), 7=(-1,0,-1)
///  8=(1,1,0), 9=(1,-1,0),10=(-1,1,0),11=(-1,-1,0)
const QUAD_FACES: [[usize; 4]; 6] = [
    // +Z face (z=1): verts 0,2,4,6. CCW from +Z: 0→4→2→6.
    [0, 4, 2, 6],
    // −Z face (z=-1): verts 1,3,5,7. CCW from −Z: 1→7→3→5.
    [1, 7, 3, 5],
    // +Y face (y=1): verts 0,1,8,10. CCW from +Y: 0→8→1→10.
    [0, 8, 1, 10],
    // −Y face (y=-1): verts 2,3,9,11. CCW from −Y: 2→11→3→9.
    [2, 11, 3, 9],
    // +X face (x=1): verts 4,5,8,9. CCW from +X: 8→4→9→5.
    [8, 4, 9, 5],
    // −X face (x=-1): verts 6,7,10,11. CCW from −X: 10→6→11→7.
    [10, 6, 11, 7],
];

fn build(co: &Cuboctahedron) -> Result<IndexedMesh, PrimitiveError> {
    if co.radius <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "radius must be > 0, got {}",
            co.radius
        )));
    }

    let region = RegionId::new(1);
    let mut mesh = IndexedMesh::new();

    // Raw circumradius = √2; scale to desired radius.
    let scale = co.radius / 2.0_f64.sqrt();
    let raw = raw_vertices();
    let cx = co.center.x;
    let cy = co.center.y;
    let cz = co.center.z;

    let verts: Vec<Point3r> = raw
        .iter()
        .map(|v| Point3r::new(cx + v[0] * scale, cy + v[1] * scale, cz + v[2] * scale))
        .collect();

    // Triangular faces
    for [ia, ib, ic] in TRI_FACES {
        let centroid = (verts[ia].coords + verts[ib].coords + verts[ic].coords) / 3.0;
        let raw_n = triangle_normal(&verts[ia], &verts[ib], &verts[ic]);
        let n = match raw_n {
            Some(n) => {
                if n.dot(&centroid) < 0.0 {
                    -n
                } else {
                    n
                }
            }
            None => Vector3r::zeros(),
        };
        let vi0 = mesh.add_vertex(verts[ia], n);
        let vi1 = mesh.add_vertex(verts[ib], n);
        let vi2 = mesh.add_vertex(verts[ic], n);
        // Ensure face winding matches outward normal direction.
        if raw_n.map(|n| n.dot(&centroid) < 0.0).unwrap_or(false) {
            mesh.add_face_with_region(vi0, vi2, vi1, region); // flipped
        } else {
            mesh.add_face_with_region(vi0, vi1, vi2, region); // normal
        }
    }

    // Square faces -- fan-triangulated from vertex 0 of each quad.
    for [ia, ib, ic, id] in QUAD_FACES {
        let centroid =
            (verts[ia].coords + verts[ib].coords + verts[ic].coords + verts[id].coords) / 4.0;
        // Triangle 1: ia, ib, ic
        {
            let raw_n = triangle_normal(&verts[ia], &verts[ib], &verts[ic]);
            let n = match raw_n {
                Some(n) => {
                    if n.dot(&centroid) < 0.0 {
                        -n
                    } else {
                        n
                    }
                }
                None => Vector3r::zeros(),
            };
            let vi0 = mesh.add_vertex(verts[ia], n);
            let vi1 = mesh.add_vertex(verts[ib], n);
            let vi2 = mesh.add_vertex(verts[ic], n);
            if raw_n.map(|n| n.dot(&centroid) < 0.0).unwrap_or(false) {
                mesh.add_face_with_region(vi0, vi2, vi1, region);
            } else {
                mesh.add_face_with_region(vi0, vi1, vi2, region);
            }
        }
        // Triangle 2: ia, ic, id
        {
            let raw_n = triangle_normal(&verts[ia], &verts[ic], &verts[id]);
            let n = match raw_n {
                Some(n) => {
                    if n.dot(&centroid) < 0.0 {
                        -n
                    } else {
                        n
                    }
                }
                None => Vector3r::zeros(),
            };
            let vi0 = mesh.add_vertex(verts[ia], n);
            let vi1 = mesh.add_vertex(verts[ic], n);
            let vi2 = mesh.add_vertex(verts[id], n);
            if raw_n.map(|n| n.dot(&centroid) < 0.0).unwrap_or(false) {
                mesh.add_face_with_region(vi0, vi2, vi1, region);
            } else {
                mesh.add_face_with_region(vi0, vi1, vi2, region);
            }
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
    fn cuboctahedron_is_watertight() {
        let mesh = Cuboctahedron::default().build().unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(
            report.is_watertight,
            "cuboctahedron must be watertight: {:?}",
            report
        );
        assert_eq!(report.euler_characteristic, Some(2));
    }

    #[test]
    fn cuboctahedron_volume_positive() {
        let mesh = Cuboctahedron {
            radius: 2.0,
            ..Cuboctahedron::default()
        }
        .build()
        .unwrap();
        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
        assert!(report.signed_volume > 0.0);
    }

    #[test]
    fn cuboctahedron_invalid_radius() {
        assert!(Cuboctahedron {
            radius: 0.0,
            ..Cuboctahedron::default()
        }
        .build()
        .is_err());
        assert!(Cuboctahedron {
            radius: -1.0,
            ..Cuboctahedron::default()
        }
        .build()
        .is_err());
    }
}
