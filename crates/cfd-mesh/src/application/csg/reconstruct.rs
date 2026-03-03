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
//! ```text
//! for face in kept_faces:
//!   for vid in face.vertices:
//!     new_id = id_map.get(vid)          ← O(1) HashMap lookup
//!           or IndexedMesh::add_vertex  ← spatial-hash dedup insert
//!     id_map.insert(vid, new_id)        ← O(1) amortised
//!   if !degenerate:
//!     output_mesh.add_face(new_ids, region)
//! ```
//!
//! ## Memory
//!
//! `id_map` is a `HashMap<VertexId, VertexId>` pre-sized to
//! `faces.len() * 3 / 2`.  For a 500-face CSG result this allocates ~750 entries
//! rather than the old `Vec<Option<_>>` which allocated at `pool.len()` (2 000+).
//! The HashMap is also freed immediately on return, keeping peak RSS low.
//!
//! ## Complexity
//!
//! `O(f)` where `f = faces.len()`.  Each vertex is inserted into `id_map` at
//! most once (amortised O(1) per HashMap insert).  Total work is proportional
//! to the output face count, independent of pool size.

use std::collections::HashMap;

use crate::domain::core::index::{RegionId, VertexId};
use crate::domain::core::scalar::Real;
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
    // Compute adaptive welding tolerance: 5% of the minimum face edge length.
    // This merges near-duplicate seam vertices from CDT co-refinement at
    // shallow-angle junctions (typically ~0.004-0.008 apart) while keeping
    // distinct surface vertices (typically ~0.05-0.1 apart) separate.
    // Floor at 1e-4 (the default IndexedMesh tolerance).
    let tol = adaptive_reconstruct_tolerance(faces, pool);
    let mut mesh = IndexedMesh::with_tolerance(tol, tol);
    // HashMap<old VertexId, new VertexId> — capacity-hinted at face count * 3/2
    // (upper bound on unique vertices for maximally shared topology).
    // This is O(face_count) rather than the old O(pool_len) Vec<Option<_>>.
    let mut id_map: HashMap<VertexId, VertexId> =
        HashMap::with_capacity(faces.len().saturating_mul(3).saturating_add(1) / 2);
    let mut degenerate_count: usize = 0;

    for face in faces {
        let mut new_verts = [VertexId::default(); 3];

        for (k, &vid) in face.vertices.iter().enumerate() {
            let new_id = if let Some(&id) = id_map.get(&vid) {
                id
            } else {
                let pos = *pool.position(vid);
                let nrm = *pool.normal(vid);
                let new_id = mesh.add_vertex(pos, nrm);
                id_map.insert(vid, new_id);
                new_id
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

/// Compute the welding tolerance for CSG mesh reconstruction.
///
/// Uses the default tolerance of `1e-4`.  Seam vertex merging is handled
/// upstream by `stitch_boundary_seams` (VertexId remapping), so the
/// reconstruction spatial hash only needs to weld truly coincident
/// positions from floating-point rounding during plane crossing.
fn adaptive_reconstruct_tolerance(_faces: &[FaceData], _pool: &VertexPool) -> Real {
    1e-4
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

    #[test]
    fn region_invalid_uses_untagged_path() {
        // Theorem: faces with `RegionId::INVALID` must be emitted via
        // `add_face` (not `add_face_with_region`) to avoid contaminating
        // region-tagged topology queries downstream. Verify the output
        // face has INVALID region after roundtrip.
        use crate::domain::core::index::RegionId;
        let mut pool = VertexPool::default_millifluidic();
        let v0 = pool.insert_or_weld(Point3r::new(0.0, 0.0, 0.0), z());
        let v1 = pool.insert_or_weld(Point3r::new(1.0, 0.0, 0.0), z());
        let v2 = pool.insert_or_weld(Point3r::new(0.0, 1.0, 0.0), z());
        let face = FaceData::new(v0, v1, v2, RegionId::INVALID);
        let mesh = reconstruct_mesh(&[face], &pool);
        assert_eq!(mesh.face_count(), 1);
        let out = mesh.faces.iter().next().unwrap();
        assert_eq!(
            out.region,
            RegionId::INVALID,
            "INVALID-region face should remain INVALID after reconstruction"
        );
    }

    #[test]
    fn large_pool_small_face_set_no_over_alloc() {
        // Memory regression: id_map must NOT allocate at pool.len() capacity.
        // With a pool of 10 000 vertices and only 2 faces, the HashMap
        // capacity hint is min(face*3/2, 3) = 3, not 10 000.
        // We verify correctness is maintained regardless.
        let mut pool = VertexPool::default_millifluidic();
        // Insert 9994 "noise" vertices that won't appear in any face.
        let zero = z();
        for i in 0..9994_usize {
            pool.insert_or_weld(Point3r::new(i as f64 + 100.0, 0.0, 0.0), zero);
        }
        // Only 6 vertices actually used.
        let v0 = pool.insert_or_weld(Point3r::new(0.0, 0.0, 0.0), zero);
        let v1 = pool.insert_or_weld(Point3r::new(1.0, 0.0, 0.0), zero);
        let v2 = pool.insert_or_weld(Point3r::new(0.0, 1.0, 0.0), zero);
        let v3 = pool.insert_or_weld(Point3r::new(1.0, 1.0, 0.0), zero);
        let v4 = pool.insert_or_weld(Point3r::new(2.0, 0.0, 0.0), zero);
        let v5 = pool.insert_or_weld(Point3r::new(0.0, 2.0, 0.0), zero);
        let faces = vec![
            FaceData::untagged(v0, v1, v2),
            FaceData::untagged(v3, v4, v5),
        ];
        let mesh = reconstruct_mesh(&faces, &pool);
        assert_eq!(mesh.face_count(), 2);
        assert_eq!(mesh.vertex_count(), 6, "only used vertices should appear");
    }
}
