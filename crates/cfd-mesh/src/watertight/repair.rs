//! Mesh repair utilities.

use crate::core::index::FaceId;
use crate::storage::edge_store::EdgeStore;
use crate::storage::face_store::FaceStore;
use crate::storage::vertex_pool::VertexPool;
use crate::topology::orientation;

/// Mesh repair operations.
pub struct MeshRepair;

impl MeshRepair {
    /// Fix inconsistent winding orientations.
    ///
    /// Returns the number of faces flipped.
    pub fn fix_orientations(
        face_store: &mut FaceStore,
        edge_store: &EdgeStore,
    ) -> usize {
        orientation::fix_orientation(face_store, edge_store)
    }

    /// Remove degenerate faces (zero-area triangles or faces with duplicate vertices).
    ///
    /// Returns the IDs of removed faces.
    pub fn remove_degenerate_faces(
        face_store: &FaceStore,
        vertex_pool: &VertexPool,
    ) -> Vec<FaceId> {
        let mut degenerate = Vec::new();

        for (fid, face) in face_store.iter_enumerated() {
            // Check for duplicate vertex references
            if face.vertices[0] == face.vertices[1]
                || face.vertices[1] == face.vertices[2]
                || face.vertices[2] == face.vertices[0]
            {
                degenerate.push(fid);
                continue;
            }

            // Check for zero area
            let a = vertex_pool.position(face.vertices[0]);
            let b = vertex_pool.position(face.vertices[1]);
            let c = vertex_pool.position(face.vertices[2]);
            let area = 0.5 * (b - a).cross(&(c - a)).norm();

            if area < crate::core::scalar::TOLERANCE {
                degenerate.push(fid);
            }
        }

        degenerate
    }
}
