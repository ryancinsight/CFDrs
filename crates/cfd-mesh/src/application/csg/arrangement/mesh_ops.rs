//! Shared low-level face-soup operations for arrangement seam repair.
//!
//! Centralizes boundary extraction, vertex-merge application, and unordered
//! face deduplication so arrangement stages use one canonical implementation.

use hashbrown::{HashMap, HashSet};

use crate::domain::core::index::VertexId;
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;

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
///
/// ## Degenerate face removal
///
/// After remapping, a face is removed if:
/// 1. **Collapsed vertices** — any two of the three vertex IDs are equal, OR
/// 2. **Collinear vertices** — all three vertex IDs are distinct but the
///    positions are (nearly) collinear, producing a zero-area triangle.
///
/// ## Theorem — Collinear Detection via Cross-Product Magnitude
///
/// For triangle $(p_0, p_1, p_2)$, the area is $\tfrac12 \lVert (p_1-p_0)
/// \times (p_2-p_0) \rVert$.  The face is degenerate when $\lVert (p_1-p_0)
/// \times (p_2-p_0) \rVert^2 \le \varepsilon$ for machine-zero tolerance
/// $\varepsilon = 10^{-30}$.  This catches all three failure modes: identical
/// vertices (zero edge), collinear vertices (zero cross product), and
/// near-degenerate slivers.  ∎
pub(crate) fn apply_vertex_merge(
    faces: &mut Vec<FaceData>,
    merge: &HashMap<VertexId, VertexId>,
    pool: &VertexPool,
) {
    if merge.is_empty() {
        return;
    }

    for face in faces.iter_mut() {
        for vid in &mut face.vertices {
            *vid = merge_root(merge, *vid);
        }
    }

    faces.retain(|f| {
        // Collapsed vertices (two or more equal after merge)
        if f.vertices[0] == f.vertices[1]
            || f.vertices[1] == f.vertices[2]
            || f.vertices[2] == f.vertices[0]
        {
            return false;
        }
        // Collinear but distinct vertices (zero-area triangle)
        let p0 = pool.position(f.vertices[0]);
        let p1 = pool.position(f.vertices[1]);
        let p2 = pool.position(f.vertices[2]);
        (p1 - p0).cross(&(p2 - p0)).norm_squared() > 1e-30
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
        use crate::domain::core::scalar::Real;

        let mut pool = VertexPool::new(1.0 as Real);
        let up = nalgebra::Vector3::<Real>::new(0.0, 0.0, 1.0);
        let v0 = pool.insert_unique(nalgebra::Point3::new(0.0, 0.0, 0.0), up);
        let v1 = pool.insert_unique(nalgebra::Point3::new(1.0, 0.0, 0.0), up);
        let v2 = pool.insert_unique(nalgebra::Point3::new(0.0, 1.0, 0.0), up);
        let v3 = pool.insert_unique(nalgebra::Point3::new(0.0, 2.0, 0.0), up);

        let mut faces = vec![
            FaceData::untagged(v0, v1, v2),
            FaceData::untagged(v0, v1, v2), // duplicate
            FaceData::untagged(v0, v2, v3),
        ];

        // Collapse v3 onto v2 -> third face becomes degenerate.
        let mut merge = HashMap::new();
        merge.insert(v3, v2);

        apply_vertex_merge(&mut faces, &merge, &pool);
        assert_eq!(faces.len(), 1);
        assert_eq!(faces[0].vertices, [v0, v1, v2]);
    }

    #[test]
    fn apply_vertex_merge_removes_collinear_zero_area_faces() {
        use crate::domain::core::scalar::Real;

        let mut pool = VertexPool::new(1.0 as Real);
        let up = nalgebra::Vector3::<Real>::new(0.0, 0.0, 1.0);
        let v0 = pool.insert_unique(nalgebra::Point3::new(0.0, 0.0, 0.0), up);
        let v1 = pool.insert_unique(nalgebra::Point3::new(1.0, 0.0, 0.0), up);
        let v2 = pool.insert_unique(nalgebra::Point3::new(0.5, 0.0, 0.0), up);
        // v0, v1, v2 are collinear — face area = 0
        let v3 = pool.insert_unique(nalgebra::Point3::new(2.0, 0.0, 0.0), up);

        // Merge v3 -> v1 so the merge map is non-empty (triggers filtering).
        let mut merge = HashMap::new();
        merge.insert(v3, v1);
        let mut faces = vec![
            FaceData::untagged(v0, v1, v2), // collinear → zero area
        ];
        apply_vertex_merge(&mut faces, &merge, &pool);
        assert!(faces.is_empty(), "collinear zero-area face should be culled");
    }
}
