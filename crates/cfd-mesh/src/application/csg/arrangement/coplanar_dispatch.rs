use std::collections::{HashMap, HashSet};

use super::super::boolean::BooleanOp;
use super::super::intersect::SnapSegment;
use super::coplanar_resolution::resolve_oriented_coplanar_group;
use super::dsu::DisjointSet;
use crate::domain::core::index::VertexId;
use crate::domain::core::scalar::Point3r;
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;

pub(crate) struct CoplanarDispatchResult {
    pub(crate) skipped_faces_by_mesh: Vec<HashSet<usize>>,
    pub(crate) representative_faces: Vec<FaceData>,
    pub(crate) result_faces: Vec<FaceData>,
}

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
