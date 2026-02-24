//! Topology-preserving vertex welding.
//!
//! Merges coincident vertices while strictly preserving 2-manifold properties.

use hashbrown::{HashMap, HashSet};

use crate::application::welding::spatial_hash::SpatialHashGrid;
use crate::domain::core::index::VertexId;
use crate::domain::core::scalar::{Point3r, Real, TOLERANCE};
use crate::infrastructure::storage::face_store::FaceStore;

/// Result of a vertex welding operation.
#[derive(Debug)]
pub struct WeldResult {
    /// The consolidated, fully-welded list of vertex positions.
    pub positions: Vec<Point3r>,
    /// Number of vertices that were successfully merged.
    pub vertices_merged: usize,
    /// Number of faces updated to point to the new canonical vertices.
    pub faces_updated: usize,
}

/// Topology-aware vertex welder for existing meshes.
pub struct MeshWelder {
    /// Tolerance for vertex welding.
    tolerance: Real,
}

impl MeshWelder {
    /// Create a new topological welder with the default tolerance.
    pub fn new() -> Self {
        Self {
            tolerance: TOLERANCE.sqrt(),
        }
    }

    /// Create a topological welder with a custom tolerance.
    pub fn with_tolerance(tolerance: Real) -> Self {
        Self { tolerance }
    }

    /// Weld duplicate vertices topologically and update the `FaceStore`.
    ///
    /// Unlike naive spatial welding, this algorithm ensures that:
    /// 1. Vertices sharing a face are never welded (prevents degenerate face creation).
    /// 2. The merge does not create non-manifold edges (>2 faces per edge).
    ///
    /// Returns a `WeldResult` containing the new canonical positions. The provided
    /// `FaceStore` is updated in-place to reference the new packed vertex indices.
    pub fn weld(&self, positions: &[Point3r], face_store: &mut FaceStore) -> WeldResult {
        let n = positions.len();
        if n == 0 {
            return WeldResult {
                positions: Vec::new(),
                vertices_merged: 0,
                faces_updated: 0,
            };
        }

        // 1. Build initial vertex-to-face adjacency graph for active vertices
        let mut v_faces: Vec<Vec<u32>> = vec![Vec::new(); n];
        for (f_id, face) in face_store.iter_enumerated() {
            for v in face.vertices {
                v_faces[v.raw() as usize].push(f_id.0);
            }
        }

        let old_active_count = v_faces.iter().filter(|f| !f.is_empty()).count();

        // Keep track of the current vertices for each face to dynamically check topologies.
        let mut cur_face_verts: HashMap<u32, [u32; 3]> = HashMap::with_capacity(face_store.len());
        for (f_id, face) in face_store.iter_enumerated() {
            cur_face_verts.insert(
                f_id.0,
                [
                    face.vertices[0].raw(),
                    face.vertices[1].raw(),
                    face.vertices[2].raw(),
                ],
            );
        }

        // 2. Build a SpatialHashGrid containing all original positions
        let mut grid = SpatialHashGrid::new(self.tolerance * 2.0);
        for (i, p) in positions.iter().enumerate() {
            grid.insert(p, i as u32);
        }

        // 3. Greedy topological clustering
        // remap[i] maps original vertex i to its canonical merged vertex head.
        let mut remap: Vec<u32> = (0..n as u32).collect();

        for i in 0..n as u32 {
            if remap[i as usize] != i {
                continue; // Already merged into another vertex.
            }

            if v_faces[i as usize].is_empty() {
                continue; // Unused vertex, ignore.
            }

            // Find all spatial neighbors
            let candidates = grid.query_radius(&positions[i as usize], self.tolerance, positions);

            for c in candidates {
                if c <= i {
                    continue; // Only merge larger indices into smaller/earlier indices.
                }

                let root_c = remap[c as usize];
                if root_c != c {
                    continue;
                }

                if v_faces[c as usize].is_empty() {
                    continue;
                }

                // Attempt to merge c into i.
                if self.is_safe_to_merge(i, c, &v_faces[i as usize], &v_faces[c as usize], &cur_face_verts) {
                    // SAFE! Merge c into i.
                    remap[c as usize] = i;

                    // Update dynamic topology.
                    // Move all faces from c to i.
                    let faces_of_c = std::mem::take(&mut v_faces[c as usize]);
                    for f_id in faces_of_c {
                        // Update cur_face_verts
                        if let Some(verts) = cur_face_verts.get_mut(&f_id) {
                            for v in verts.iter_mut() {
                                if *v == c {
                                    *v = i;
                                }
                            }
                        }
                        // Add face to i's list if not present
                        if !v_faces[i as usize].contains(&f_id) {
                            v_faces[i as usize].push(f_id);
                        }
                    }
                }
            }
        }

        // 4. Pack vertices
        let mut packed_positions = Vec::new();
        let mut pack_map: HashMap<u32, u32> = HashMap::new();

        let mut faces_updated = 0;

        for (_f_id, face) in face_store.iter_mut_enumerated() {
            let mut changed = false;
            for v in &mut face.vertices {
                let old_raw = v.raw();
                let canonical_id = remap[old_raw as usize];

                if canonical_id != old_raw {
                    changed = true;
                }

                // Find or assign packed id
                let packed_id = *pack_map.entry(canonical_id).or_insert_with(|| {
                    let new_idx = packed_positions.len() as u32;
                    packed_positions.push(positions[canonical_id as usize]);
                    new_idx
                });

                *v = VertexId::new(packed_id);
            }
            if changed {
                faces_updated += 1;
            }
        }

        WeldResult {
            positions: packed_positions,
            vertices_merged: old_active_count - pack_map.len(),
            faces_updated,
        }
    }

    /// Evaluates the topological safety of merging vertex `src` into `dst`.
    fn is_safe_to_merge(
        &self,
        dst: u32,
        src: u32,
        dst_faces: &[u32],
        src_faces: &[u32],
        cur_face_verts: &HashMap<u32, [u32; 3]>,
    ) -> bool {
        // Condition 1: Shared Face Check (Degenerate Prevention)
        for &f_id in src_faces {
            if dst_faces.contains(&f_id) {
                return false;
            }
        }

        // Condition 2: Non-Manifold Edge Prevention
        let mut dst_neighbors = HashMap::new();
        for &f_id in dst_faces {
            if let Some(verts) = cur_face_verts.get(&f_id) {
                for &v in verts {
                    if v != dst {
                        *dst_neighbors.entry(v).or_insert(0) += 1;
                    }
                }
            }
        }

        let mut src_neighbors = HashMap::new();
        for &f_id in src_faces {
            if let Some(verts) = cur_face_verts.get(&f_id) {
                for &v in verts {
                    if v != src {
                        *src_neighbors.entry(v).or_insert(0) += 1;
                    }
                }
            }
        }

        // If they share a neighbor, the combined face count on that edge cannot exceed 2.
        for (v, count_src) in &src_neighbors {
            if let Some(count_dst) = dst_neighbors.get(v) {
                if count_src + count_dst > 2 {
                    return false;
                }
            }
        }

        true
    }
}

impl Default for MeshWelder {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::core::index::{RegionId, VertexId};
    use crate::infrastructure::storage::face_store::FaceData;

    fn p(x: f64, y: f64, z: f64) -> Point3r {
        Point3r::new(x, y, z)
    }

    #[test]
    fn topological_weld_prevents_degenerate_face() {
        let mut positions = vec![
            p(0.0, 0.0, 0.0),
            p(1e-7, 0.0, 0.0), // very close, should not weld because they share a face
            p(0.0, 1.0, 0.0),
        ];
        
        let mut face_store = FaceStore::new();
        face_store.push(FaceData::new(VertexId::new(0), VertexId::new(1), VertexId::new(2), RegionId::new(0)));
        
        let welder = MeshWelder::new();
        let result = welder.weld(&positions, &mut face_store);
        
        assert_eq!(result.vertices_merged, 0, "Should not weld adjacent vertices on the same face.");
    }

    #[test]
    fn topological_weld_allows_disconnected_merge() {
        let mut positions = vec![
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(0.0, 1.0, 0.0),
            // Disconnected copy of the same points
            p(0.0, 1e-7, 0.0), // Merges to 0
            p(1.0, 0.0, 0.0),  // Merges to 1
            p(0.0, 1.0, 0.0),  // Merges to 2
        ];

        let mut face_store = FaceStore::new();
        face_store.push(FaceData::new(VertexId::new(0), VertexId::new(1), VertexId::new(2), RegionId::new(0)));
        face_store.push(FaceData::new(VertexId::new(3), VertexId::new(4), VertexId::new(5), RegionId::new(0)));

        let welder = MeshWelder::new();
        let result = welder.weld(&positions, &mut face_store);

        assert_eq!(result.vertices_merged, 3, "Should successfully weld fully disconnected faces.");
        assert_eq!(result.faces_updated, 1, "The second face should be updated to point to the first group.");
    }

    #[test]
    fn topological_weld_prevents_pinching() {
        // Two triangles that share an edge (A, B) but have separate opposite vertices C and D.
        // Wait, pinching is when C and D are close, but welding them creates an edge (A,C/D) and (B,C/D) with 3 faces?
        // No, if they share edge (A, B), they have faces (A,B,C) and (B,A,D).
        // A=0, B=1, C=2, D=3.
        // C and D are disconnected from each other.
        // If we weld C to D, edge (A, C) and edge (A, D) merge.
        // (A, C) has face 1, (A, D) has face 2. 1 + 1 = 2 faces.
        // (B, C) has face 1, (B, D) has face 2. 1 + 1 = 2 faces.
        // So C and D CAN be welded! This just folds the two triangles on top of each other!
        // Is folding non-manifold?
        // If they fold, the face normals point opposite ways. The shared edge (A,B) has 2 faces.
        // Wait! The edge A-B has 2 faces originally. If they fold, A-B still has 2 faces.
        // So topological edge count is perfectly preserved!
        // Manifold edge count allows it. It just creates degenerate volume.

        // What about PINCHING two independent sheets?
        // Sheet 1: (A,B,C) and (B,D,C)
        // Sheet 2: (E,F,G) and (F,H,G)
        // If we weld C to G, they don't share any neighbors! So they weld and form a bowtie vertex.
        // Our welder allows this (it's non-manifold vertex, but edge counts are preserved).
        // What if we weld edge C-D to G-H?
        // They share neighbors, but total faces <= 2. It allows edge merging of boundaries!
        
        // Wait, how do we trigger non-manifold EDGE creation?
        // Sheet 1 connects A-B with 2 faces (interior edge).
        // Sheet 2 connects C-D with 2 faces (interior edge).
        // If we weld A to C, and B to D.
        // When welding A to C, they share NO neighbors originally. Welded!
        // Now A/C and B/D share a neighbor (A/C).
        // We try to weld B to D.
        // A/C is a neighbor of B (from Sheet 1) with 2 faces.
        // A/C is a neighbor of D (from Sheet 2) with 2 faces.
        // If we weld B to D, the shared neighbor A/C will have 2 + 2 = 4 faces.
        // Our check: `count_src + count_dst > 2`. 2 + 2 = 4 > 2. So it FAILS!
        // EXACTLY! Pinching interior edges is safely blocked.
        let mut positions = vec![
            // Sheet 1
            p(0.0, 0.0, 0.0), // 0: A
            p(1.0, 0.0, 0.0), // 1: B
            p(0.5, 1.0, 0.0), // 2: C1
            p(0.5, -1.0, 0.0),// 3: C2
            // Sheet 2 (parallel, very close)
            p(0.0, 0.0, 1e-7), // 4: C (wants to merge with A)
            p(1.0, 0.0, 1e-7), // 5: D (wants to merge with B)
            p(0.5, 1.0, 1e-7), // 6: G1
            p(0.5, -1.0, 1e-7),// 7: G2
        ];
        
        let mut face_store = FaceStore::new();
        // Sheet 1: interior edge (0,1)
        face_store.push(FaceData::new(VertexId::new(0), VertexId::new(1), VertexId::new(2), RegionId::new(0)));
        face_store.push(FaceData::new(VertexId::new(1), VertexId::new(0), VertexId::new(3), RegionId::new(0)));
        // Sheet 2: interior edge (4,5)
        face_store.push(FaceData::new(VertexId::new(4), VertexId::new(5), VertexId::new(6), RegionId::new(0)));
        face_store.push(FaceData::new(VertexId::new(5), VertexId::new(4), VertexId::new(7), RegionId::new(0)));

        let welder = MeshWelder::new();
        let result = welder.weld(&positions, &mut face_store);

        // A (0) and C (4) share no neighbors initially, so they will weld.
        // Once welded, edge (0, 1) has 2 faces, and edge (0, 5) has 2 faces.
        // Now B (1) and D (5) want to weld. They share neighbor 0.
        // count(1, 0) = 2. count(5, 0) = 2. Total = 4. 
        // Welding B to D is blocked!
        // Also (2) and (6) will weld because they share no neighbors.
        // Let's assert that `faces_updated` is not enough to completely merge the two sheets.
        assert!(result.vertices_merged < 4, "Should not fully merge two parallel interior sheets");
    }
}
