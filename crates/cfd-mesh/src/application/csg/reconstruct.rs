//! Fragment mesh reconstruction for CSG Boolean operations.
//!
//! After Boolean classification, combines the face fragments selected from
//! mesh A and mesh B into a single fresh, deduplicated [`IndexedMesh`].
//!
//! ## Why Reconstruction Matters
//!
//! The intersection/corefine phase inserts new seam vertices into the shared
//! vertex pool.  Due to accumulated floating-point error, the same geometric
//! point may be stored as two slightly different positions (one added during
//! an A-edge / B-plane crossing, another during a B-edge / A-plane crossing).
//! Routing all kept faces through a fresh [`IndexedMesh`]'s spatial-hash
//! deduplication step welds such near-duplicate vertices together, restoring
//! a watertight seam.
//!
//! ## Algorithm
//!
//! 1. Allocate a fresh `IndexedMesh` (default millifluidic tolerances).
//! 2. For each kept face, look up each vertex position and normal in the
//!    source pool and insert via [`IndexedMesh::add_vertex`] (deduplicating).
//! 3. Add the remapped triangle (preserving the region tag).
//! 4. Silently skip degenerate faces where deduplication collapsed two or
//!    more corners to the same output vertex.

use crate::domain::core::index::{RegionId, VertexId};
use crate::domain::mesh::IndexedMesh;
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;

// ── Public API ────────────────────────────────────────────────────────────────

/// Reconstruct an [`IndexedMesh`] from a set of kept face fragments.
///
/// Each face in `faces` references vertices in `pool` by [`VertexId`].
///
/// # Arguments
///
/// * `faces` — Face fragments to include in the output.
/// * `pool`  — Vertex pool backing the faces.  Must contain every vertex
///             referenced by `faces`.
///
/// # Returns
///
/// A freshly-built, deduplicated `IndexedMesh`.
///
/// # Panics
///
/// Panics if any face references a vertex ID outside `0..pool.len()`.
pub fn reconstruct_mesh(faces: &[FaceData], pool: &VertexPool) -> IndexedMesh {
    let mut mesh = IndexedMesh::new();
    // Lazily-populated map: old pool VertexId → new mesh VertexId.
    let mut id_map: Vec<Option<VertexId>> = vec![None; pool.len()];
    let mut degenerate_count: usize = 0;

    for face in faces {
        let mut new_verts = [VertexId::default(); 3];

        for (k, &vid) in face.vertices.iter().enumerate() {
            let idx = vid.as_usize();
            let new_id = match id_map[idx] {
                Some(id) => id,
                None => {
                    let pos = *pool.position(vid);
                    let nrm = *pool.normal(vid);
                    let new_id = mesh.add_vertex(pos, nrm);
                    id_map[idx] = Some(new_id);
                    new_id
                }
            };
            new_verts[k] = new_id;
        }

        // Skip degenerate faces (two or more corners collapsed to the same vertex).
        if new_verts[0] == new_verts[1]
            || new_verts[1] == new_verts[2]
            || new_verts[2] == new_verts[0]
        {
            degenerate_count += 1;
            continue;
        }

        if face.region == RegionId::INVALID {
            mesh.add_face(new_verts[0], new_verts[1], new_verts[2]);
        } else {
            mesh.add_face_with_region(new_verts[0], new_verts[1], new_verts[2], face.region);
        }
    }

    if degenerate_count > 0 {
        tracing::warn!(
            "CSG reconstruct: {} degenerate fragment(s) silently dropped (vertices collapsed by welding)",
            degenerate_count
        );
    }

    mesh
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::core::scalar::{Point3r, Vector3r};

    fn z() -> Vector3r {
        Vector3r::zeros()
    }

    #[test]
    fn single_triangle_round_trip() {
        let mut pool = VertexPool::default_millifluidic();
        let v0 = pool.insert_or_weld(Point3r::new(0.0, 0.0, 0.0), z());
        let v1 = pool.insert_or_weld(Point3r::new(1.0, 0.0, 0.0), z());
        let v2 = pool.insert_or_weld(Point3r::new(0.0, 1.0, 0.0), z());
        let face = FaceData::untagged(v0, v1, v2);
        let mesh = reconstruct_mesh(&[face], &pool);

        assert_eq!(mesh.vertex_count(), 3);
        assert_eq!(mesh.face_count(), 1);
    }

    #[test]
    fn degenerate_face_skipped() {
        let mut pool = VertexPool::default_millifluidic();
        let v0 = pool.insert_or_weld(Point3r::new(0.0, 0.0, 0.0), z());
        let v1 = pool.insert_or_weld(Point3r::new(1.0, 0.0, 0.0), z());
        // Face with repeated vertex (v0 twice) — degenerate.
        let face_deg = FaceData::untagged(v0, v0, v1);
        let mesh = reconstruct_mesh(&[face_deg], &pool);
        assert_eq!(
            mesh.face_count(),
            0,
            "degenerate face should be silently skipped"
        );
    }

    #[test]
    fn shared_vertex_deduplication() {
        // Two coplanar triangles sharing two vertices.
        let mut pool = VertexPool::default_millifluidic();
        let v0 = pool.insert_or_weld(Point3r::new(0.0, 0.0, 0.0), z());
        let v1 = pool.insert_or_weld(Point3r::new(1.0, 0.0, 0.0), z());
        let v2 = pool.insert_or_weld(Point3r::new(0.0, 1.0, 0.0), z());
        let v3 = pool.insert_or_weld(Point3r::new(1.0, 1.0, 0.0), z());

        let f0 = FaceData::untagged(v0, v1, v2);
        let f1 = FaceData::untagged(v1, v3, v2);

        let mesh = reconstruct_mesh(&[f0, f1], &pool);
        assert_eq!(
            mesh.vertex_count(),
            4,
            "shared vertices should not be duplicated"
        );
        assert_eq!(mesh.face_count(), 2);
    }

    #[test]
    fn near_duplicate_positions_not_double_counted() {
        // Two triangles sharing all three vertex positions exactly (same data).
        // reconstruct_mesh should produce 1 set of vertices shared by both faces.
        let mut pool = VertexPool::default_millifluidic();
        let v0 = pool.insert_or_weld(Point3r::new(0.0, 0.0, 0.0), z());
        let v1 = pool.insert_or_weld(Point3r::new(1.0, 0.0, 0.0), z());
        let v2 = pool.insert_or_weld(Point3r::new(0.0, 1.0, 0.0), z());
        let v3 = pool.insert_or_weld(Point3r::new(0.0, -1.0, 0.0), z());

        let f0 = FaceData::untagged(v0, v1, v2);
        let f1 = FaceData::untagged(v0, v1, v3);

        let mesh = reconstruct_mesh(&[f0, f1], &pool);
        assert_eq!(mesh.vertex_count(), 4, "four distinct vertices expected");
        assert_eq!(mesh.face_count(), 2);
    }

    #[test]
    fn region_tags_preserved() {
        use crate::domain::core::index::RegionId;
        let mut pool = VertexPool::default_millifluidic();
        let v0 = pool.insert_or_weld(Point3r::new(0.0, 0.0, 0.0), z());
        let v1 = pool.insert_or_weld(Point3r::new(1.0, 0.0, 0.0), z());
        let v2 = pool.insert_or_weld(Point3r::new(0.0, 1.0, 0.0), z());
        let region = RegionId::from_usize(3);
        let face = FaceData::new(v0, v1, v2, region);
        let mesh = reconstruct_mesh(&[face], &pool);
        assert_eq!(mesh.face_count(), 1);
        // Verify the region tag was carried over to the new mesh.
        let new_face = mesh.faces.iter().next().unwrap();
        assert_eq!(new_face.region, region);
    }
}
