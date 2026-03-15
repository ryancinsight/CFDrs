//! Coplanar face dispatch for arrangement CSG.
//!
//! When two operand meshes share faces that lie in the same geometric plane,
//! the standard co-refinement + GWN classification pipeline cannot resolve
//! them: the intersection curve is degenerate (a polygon, not a curve) and
//! the generalized winding number is discontinuous on the shared plane.
//!
//! This module detects groups of coplanar faces across operands, resolves
//! each group via 2D polygon Boolean operations (delegated to
//! [`resolve_oriented_coplanar_group`]), and injects the resulting seam
//! vertices into adjacent non-coplanar ("barrel") faces so the downstream
//! CDT co-refinement can stitch everything into a watertight surface.
//!
//! ## Algorithm — Coplanar Dispatch
//!
//! 1. **Build coplanar groups** — Coplanar face pairs `(mesh_a, face_a, mesh_b, face_b)`
//!    are merged via DSU into connected components.  Each component shares
//!    one geometric plane.
//!
//! 2. **Resolve each group** — For each connected component, project faces
//!    onto the shared plane and compute 2D polygon Booleans (union /
//!    intersection / difference) via [`resolve_oriented_coplanar_group`].
//!
//! 3. **Inject seam vertices** — Boundary vertices of the resolved 2D result
//!    are injected as Steiner points into adjacent barrel faces via
//!    [`inject_cap_seam_into_barrels`], preventing T-junctions at the
//!    coplanar–non-coplanar boundary.
//!
//! ## Theorem — DSU Group Correctness
//!
//! All faces that share a geometric plane are placed in the same DSU
//! component if and only if there exists a chain of coplanar-pair entries
//! connecting them.
//!
//! *Proof.*  DSU union is transitive.  Each coplanar pair `(fa, fb)` causes
//! `union(id(fa), id(fb))`.  If `fa ∼ fb` and `fb ∼ fc` via two pairs,
//! then `find(fa) = find(fb) = find(fc)` after both unions.  Conversely,
//! faces that share no pair chain remain in disjoint components.  ∎
//!
//! ## References
//!
//! - Attene, M. (2010). "A lightweight approach to repairing digitized
//!   polygon meshes." *The Visual Computer*, 26(11), 1393–1406.
//! - Bernstein, G., & Fussell, D. (2009). "Fast, exact, linear Booleans."
//!   *Computer Graphics Forum*, 28(5), 1269–1278.

use hashbrown::{HashMap, HashSet};

use super::super::boolean::BooleanOp;
use super::super::intersect::SnapSegment;
use super::coplanar_resolution::resolve_oriented_coplanar_group;
use super::dsu::DisjointSet;
use crate::domain::core::index::VertexId;
use crate::domain::core::scalar::Point3r;
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;

/// Result of coplanar face dispatch.
///
/// Contains the faces that were resolved by 2D polygon Booleans, plus
/// bookkeeping of which original faces were consumed (so they can be
/// excluded from the standard 3D co-refinement pipeline).
pub(crate) struct CoplanarDispatchResult {
    /// Per-mesh set of face indices that were handled by coplanar resolution
    /// and should be skipped during 3D co-refinement.
    pub(crate) skipped_faces_by_mesh: Vec<HashSet<usize>>,
    /// One representative face per coplanar group, used to detect and skip
    /// fragments that lie on an already-resolved coplanar plane.
    pub(crate) representative_faces: Vec<FaceData>,
    /// Resolved 2D Boolean result faces, already in 3D coordinates.
    pub(crate) result_faces: Vec<FaceData>,
}

/// Dispatch coplanar face groups for 2D polygon Boolean resolution.
///
/// Groups coplanar face pairs via DSU, resolves each group through
/// [`resolve_oriented_coplanar_group`], and injects boundary seam
/// vertices into adjacent barrel faces.
///
/// ## Complexity
///
/// | Phase | Time | Space |
/// |-------|------|-------|
/// | DSU grouping | O(P α(P)) | O(F) |
/// | Per-group resolution | O(Σ Gᵢ log Gᵢ) | O(max Gᵢ) |
/// | Seam injection | O(F × S) per mesh | O(S) |
///
/// where P = number of coplanar pairs, F = total unique faces, Gᵢ = faces
/// in group i, S = boundary seam vertices per group.
pub(crate) fn dispatch_boolean_coplanar(
    op: BooleanOp,
    n_meshes: usize,
    meshes: &[Vec<FaceData>],
    coplanar_pairs: &[(usize, usize, usize, usize)],
    pool: &mut VertexPool,
    segs: &mut [Vec<Vec<SnapSegment>>],
) -> CoplanarDispatchResult {
    let mut skipped_faces_by_mesh = vec![HashSet::new(); n_meshes];
    if coplanar_pairs.is_empty() {
        return CoplanarDispatchResult {
            skipped_faces_by_mesh,
            representative_faces: Vec::new(),
            result_faces: Vec::new(),
        };
    }

    // 1. Gather unique indices
    let mut unique_list = Vec::new();
    let mut set = HashSet::new();
    for &(ma, fa, mb, fb) in coplanar_pairs {
        if set.insert((ma, fa)) {
            unique_list.push((ma, fa));
        }
        if set.insert((mb, fb)) {
            unique_list.push((mb, fb));
        }
    }

    let mut id_map = HashMap::new();
    for (i, &key) in unique_list.iter().enumerate() {
        id_map.insert(key, i);
    }

    let mut dsu = DisjointSet::new(unique_list.len());
    for &(ma, fa, mb, fb) in coplanar_pairs {
        dsu.union(id_map[&(ma, fa)], id_map[&(mb, fb)]);
    }

    let mut groups: HashMap<usize, Vec<(usize, usize)>> = HashMap::new();
    for (i, &key) in unique_list.iter().enumerate() {
        groups.entry(dsu.find(i)).or_default().push(key);
    }

    let mut representative_faces = Vec::new();
    let mut result = Vec::new();
    for group in groups.values() {
        // Collect all faces to find a plane basis
        let all_faces: Vec<FaceData> = group.iter().map(|&(m, f)| meshes[m][f]).collect();
        if let Some(basis) = crate::application::csg::coplanar::detect_flat_plane(&all_faces, pool)
        {
            // Valid plane group!
            let mut faces_by_mesh: Vec<Vec<FaceData>> = vec![Vec::new(); n_meshes];
            representative_faces.push(all_faces[0]);

            for &(mid, fid) in group {
                faces_by_mesh[mid].push(meshes[mid][fid]);
                skipped_faces_by_mesh[mid].insert(fid);
            }
            let res = resolve_oriented_coplanar_group(op, &faces_by_mesh, &basis, pool);

            if !res.is_empty() {
                let mut seam_vids_sorted = boundary_vertex_ids(&res);
                seam_vids_sorted.sort_unstable();
                let seam_positions: Vec<Point3r> = seam_vids_sorted
                    .iter()
                    .map(|&vid| *pool.position(vid))
                    .collect();

                if !seam_positions.is_empty() {
                    for mid in 0..n_meshes {
                        if !skipped_faces_by_mesh[mid].is_empty() {
                            let used_in_mesh: HashSet<usize> = group
                                .iter()
                                .filter(|&&(m, _)| m == mid)
                                .map(|&(_, f)| f)
                                .collect();
                            if used_in_mesh.is_empty() {
                                continue;
                            }
                            crate::application::csg::arrangement::propagate::inject_cap_seam_into_barrels(
                                &meshes[mid],
                                &used_in_mesh,
                                &basis.origin,
                                &basis.normal,
                                &seam_positions,
                                &mut segs[mid],
                                pool,
                            );
                        }
                    }
                }

                result.extend(res);
            }
        }
    }

    CoplanarDispatchResult {
        skipped_faces_by_mesh,
        representative_faces,
        result_faces: result,
    }
}

/// Extract boundary vertex IDs from a face soup.
///
/// A half-edge is on the boundary if and only if its canonical (sorted)
/// undirected edge appears exactly once across all faces.  Returns the
/// deduplicated set of vertex IDs incident to such boundary edges.
fn boundary_vertex_ids(faces: &[FaceData]) -> Vec<VertexId> {
    let mut edge_counts: HashMap<(VertexId, VertexId), usize> =
        HashMap::with_capacity(faces.len() * 3);
    for face in faces {
        for (a, b) in [
            (face.vertices[0], face.vertices[1]),
            (face.vertices[1], face.vertices[2]),
            (face.vertices[2], face.vertices[0]),
        ] {
            let edge = if a < b { (a, b) } else { (b, a) };
            *edge_counts.entry(edge).or_default() += 1;
        }
    }

    let mut boundary_vids = HashSet::new();
    for ((a, b), count) in edge_counts {
        if count == 1 {
            boundary_vids.insert(a);
            boundary_vids.insert(b);
        }
    }

    boundary_vids.into_iter().collect()
}
