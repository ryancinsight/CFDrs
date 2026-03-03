//! Shared low-level face-soup operations for arrangement seam repair.
//!
//! Centralizes boundary extraction, vertex-merge application, and unordered
//! face deduplication so arrangement stages use one canonical implementation.

use std::collections::{HashMap, HashSet};

use crate::domain::core::index::VertexId;
use crate::infrastructure::storage::face_store::FaceData;

/// Build directed boundary half-edges for a face soup.
///
/// Returns edges `(a,b)` that appear exactly once and do not have reverse
/// counterpart `(b,a)`.
pub(crate) fn boundary_half_edges(faces: &[FaceData]) -> Vec<(VertexId, VertexId)> {
    let mut half_edges: HashMap<(VertexId, VertexId), u32> = HashMap::new();
    for face in faces {
        let v = face.vertices;
        for i in 0..3 {
            let j = (i + 1) % 3;
            *half_edges.entry((v[i], v[j])).or_insert(0) += 1;
        }
    }

    half_edges
        .iter()
        .filter(|&(&(a, b), &count)| a != b && count == 1 && !half_edges.contains_key(&(b, a)))
        .map(|(&edge, _)| edge)
        .collect()
}

#[inline]
pub(crate) fn merge_root(merge: &HashMap<VertexId, VertexId>, mut v: VertexId) -> VertexId {
    while let Some(&next) = merge.get(&v) {
        if next == v {
            break;
        }
        v = next;
    }
    v
}

/// Apply a vertex merge map to faces and remove degenerate/duplicate triangles.
///
/// `merge` is interpreted as `discard -> keep` mapping. Chains are followed
/// transitively (`a->b`, `b->c` => `a->c`).
pub(crate) fn apply_vertex_merge(faces: &mut Vec<FaceData>, merge: &HashMap<VertexId, VertexId>) {
    if merge.is_empty() {
        return;
    }

    for face in faces.iter_mut() {
        for vid in &mut face.vertices {
            *vid = merge_root(merge, *vid);
        }
    }

    faces.retain(|f| {
        f.vertices[0] != f.vertices[1]
            && f.vertices[1] != f.vertices[2]
            && f.vertices[2] != f.vertices[0]
    });
    dedup_faces_unordered(faces);
}

/// Remove duplicate faces by canonicalized unordered vertex triplet.
pub(crate) fn dedup_faces_unordered(faces: &mut Vec<FaceData>) {
    let mut seen: HashSet<[VertexId; 3]> = HashSet::with_capacity(faces.len());
    faces.retain(|f| {
        let mut key = f.vertices;
        key.sort();
        seen.insert(key)
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::core::index::RegionId;

    #[test]
    fn boundary_half_edges_empty_for_closed_tetra_surface() {
        let v0 = VertexId::new(0);
        let v1 = VertexId::new(1);
        let v2 = VertexId::new(2);
        let v3 = VertexId::new(3);
        let faces = vec![
            FaceData::new(v0, v1, v2, RegionId::INVALID),
            FaceData::new(v0, v3, v1, RegionId::INVALID),
            FaceData::new(v1, v3, v2, RegionId::INVALID),
            FaceData::new(v2, v3, v0, RegionId::INVALID),
        ];
        assert!(boundary_half_edges(&faces).is_empty());
    }

    #[test]
    fn apply_vertex_merge_removes_degenerate_and_duplicate_faces() {
        let v0 = VertexId::new(0);
        let v1 = VertexId::new(1);
        let v2 = VertexId::new(2);
        let v3 = VertexId::new(3);

        let mut faces = vec![
            FaceData::untagged(v0, v1, v2),
            FaceData::untagged(v0, v1, v2), // duplicate
            FaceData::untagged(v0, v2, v3),
        ];

        // Collapse v3 onto v2 -> third face becomes degenerate.
        let mut merge = HashMap::new();
        merge.insert(v3, v2);

        apply_vertex_merge(&mut faces, &merge);
        assert_eq!(faces.len(), 1);
        assert_eq!(faces[0].vertices, [v0, v1, v2]);
    }
}
