//! Convert a Delaunay / CDT triangulation to an [`IndexedMesh`].
//!
//! The 2-D triangulation is embedded in 3-D by setting `z = 0` for all
//! vertices, with normals pointing in the `+z` direction.
//!
//! # Vertex Deduplication
//!
//! Since triangulation vertices are exact (no floating-point mesh welding),
//! we insert them directly into the `VertexPool` with a tight tolerance.
//! The spatial-hash deduplication in `VertexPool` ensures no duplicates.

use nalgebra::{Point3, Vector3};

use crate::application::delaunay::triangulation::bowyer_watson::DelaunayTriangulation;
use crate::domain::core::index::VertexId;
use crate::domain::core::scalar::Real;
use crate::domain::mesh::IndexedMesh;

/// Convert a `DelaunayTriangulation` into an `IndexedMesh<f64>`.
///
/// Only interior triangles (not touching super-triangle vertices) are
/// included.  The triangulation's 2-D vertices are embedded at `z = 0`
/// with normals `(0, 0, 1)`.
///
/// # Arguments
///
/// - `dt`: the Delaunay triangulation (or CDT) to convert
///
/// # Returns
///
/// A new `IndexedMesh<f64>` containing the triangulation.
pub fn to_indexed_mesh(dt: &DelaunayTriangulation) -> IndexedMesh<f64> {
    let mut mesh = IndexedMesh::<f64>::new();

    let normal = Vector3::new(0.0, 0.0, 1.0);

    // Map from PslgVertexId → VertexId in the IndexedMesh.
    let n = dt.vertices().len();
    let mut id_map: Vec<Option<VertexId>> = vec![None; n];

    // Only insert vertices referenced by interior triangles.
    for (_, tri) in dt.interior_triangles() {
        for &vid in &tri.vertices {
            if id_map[vid.idx()].is_none() {
                let v = dt.vertex(vid);
                let pos = Point3::new(v.x, v.y, 0.0);
                let mid = mesh.add_vertex(pos, normal);
                id_map[vid.idx()] = Some(mid);
            }
        }
    }

    // Insert faces.
    for (_, tri) in dt.interior_triangles() {
        let v0 = id_map[tri.vertices[0].idx()].expect("vertex not mapped");
        let v1 = id_map[tri.vertices[1].idx()].expect("vertex not mapped");
        let v2 = id_map[tri.vertices[2].idx()].expect("vertex not mapped");
        mesh.add_face(v0, v1, v2);
    }

    mesh
}

/// Convert a `DelaunayTriangulation` into an `IndexedMesh<f64>` with a
/// prescribed z-coordinate for the embedding plane.
pub fn to_indexed_mesh_at_z(dt: &DelaunayTriangulation, z: Real) -> IndexedMesh<f64> {
    let mut mesh = IndexedMesh::<f64>::new();
    let normal = Vector3::new(0.0, 0.0, 1.0);

    let n = dt.vertices().len();
    let mut id_map: Vec<Option<VertexId>> = vec![None; n];

    for (_, tri) in dt.interior_triangles() {
        for &vid in &tri.vertices {
            if id_map[vid.idx()].is_none() {
                let v = dt.vertex(vid);
                let pos = Point3::new(v.x, v.y, z);
                let mid = mesh.add_vertex(pos, normal);
                id_map[vid.idx()] = Some(mid);
            }
        }
    }

    for (_, tri) in dt.interior_triangles() {
        let v0 = id_map[tri.vertices[0].idx()].expect("vertex not mapped");
        let v1 = id_map[tri.vertices[1].idx()].expect("vertex not mapped");
        let v2 = id_map[tri.vertices[2].idx()].expect("vertex not mapped");
        mesh.add_face(v0, v1, v2);
    }

    mesh
}
