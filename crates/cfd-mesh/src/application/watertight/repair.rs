//! Mesh repair utilities.

use crate::application::csg::arrangement::stitch;
use crate::domain::core::index::FaceId;
use crate::domain::topology::orientation;
use crate::infrastructure::storage::edge_store::EdgeStore;
use crate::infrastructure::storage::face_store::{FaceData, FaceStore};
use crate::infrastructure::storage::vertex_pool::VertexPool;

/// Mesh repair operations.
pub struct MeshRepair;

impl MeshRepair {
    /// Fix inconsistent winding orientations.
    ///
    /// Returns the number of faces flipped.
    pub fn fix_orientations(face_store: &mut FaceStore, edge_store: &EdgeStore) -> usize {
        orientation::fix_orientation(face_store, edge_store)
    }

    /// Remove degenerate faces (zero-area triangles or faces with duplicate vertices).
    ///
    /// Returns the IDs of removed faces.
    #[must_use]
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

            if area < crate::domain::core::scalar::TOLERANCE {
                degenerate.push(fid);
            }
        }

        degenerate
    }

    /// Apply bounded iterative snap-rounding and loop filling to reduce
    /// residual boundary edges without unbounded repair churn.
    ///
    /// Returns the number of repair passes that strictly reduced the boundary
    /// edge count.
    pub fn iterative_boundary_stitch(
        face_store: &mut FaceStore,
        vertex_pool: &VertexPool,
        max_passes: usize,
    ) -> usize {
        let mut boundary_edges = EdgeStore::from_face_store(face_store).boundary_edges().len();
        let mut improved_passes = 0usize;

        for _ in 0..max_passes {
            if boundary_edges == 0 {
                break;
            }

            let mut repaired_faces: Vec<FaceData> = face_store.iter().copied().collect();
            stitch::snap_round_tjunctions(&mut repaired_faces, vertex_pool);
            stitch::fill_boundary_loops(&mut repaired_faces, vertex_pool);

            let mut repaired_store = FaceStore::new();
            for face in repaired_faces {
                repaired_store.push(face);
            }

            let next_boundary_edges =
                EdgeStore::from_face_store(&repaired_store).boundary_edges().len();
            if next_boundary_edges >= boundary_edges {
                break;
            }

            *face_store = repaired_store;
            boundary_edges = next_boundary_edges;
            improved_passes += 1;
        }

        improved_passes
    }
}
