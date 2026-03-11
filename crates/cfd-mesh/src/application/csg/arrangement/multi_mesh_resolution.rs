//! Multi-mesh fragment resolution for canonical Boolean arrangement.
//!
//! Handles post-corefine fragment consolidation and classification for the
//! canonical Boolean engine across any operand count.

use std::collections::{HashMap, HashSet};

use super::boolean_csg::BooleanOp;
use super::classify::{
    centroid, classify_fragment_prepared, prepare_classification_faces, tri_normal,
};
use super::fragment_analysis::is_degenerate_sliver_with_normal;
use super::tiebreaker::FragmentClass;
use crate::domain::core::index::VertexId;
use crate::domain::core::scalar::Point3r;
use crate::domain::topology::predicates::{orient3d, Sign};
use crate::domain::geometry::aabb::Aabb;
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;

/// Face classification record across the generalized Boolean fragment set.
#[derive(Clone, Debug)]
pub struct BooleanFragmentRecord {
    pub face: FaceData,
    pub mesh_idx: usize,
    pub parent_idx: usize,
}

#[derive(Clone, Copy, Debug)]
struct CachedClassification {
    class: FragmentClass,
}

struct CoplanarPlaneInfo {
    a: Point3r,
    b: Point3r,
    c: Point3r,
}

/// Resolve Boolean fragments into result faces using one survivorship policy
/// for both binary and dense N-way inputs.
pub(crate) fn resolve_multi_mesh_fragments(
    op: BooleanOp,
    frags: &mut Vec<BooleanFragmentRecord>,
    meshes: &[Vec<FaceData>],
    mesh_aabbs: &[Aabb],
    coplanar_group_faces: &[FaceData],
    mut coplanar_result_faces: Vec<FaceData>,
    pool: &VertexPool,
) -> Vec<FaceData> {
    consolidate_cross_mesh_vertices(frags, pool);

    let coplanar_plane_infos: Vec<CoplanarPlaneInfo> = coplanar_group_faces
        .iter()
        .map(|rep_face| CoplanarPlaneInfo {
            a: *pool.position(rep_face.vertices[0]),
            b: *pool.position(rep_face.vertices[1]),
            c: *pool.position(rep_face.vertices[2]),
        })
        .collect();

    let mut prepared_meshes = Vec::with_capacity(meshes.len());
    for mesh_faces in meshes {
        prepared_meshes.push(prepare_classification_faces(mesh_faces, pool));
    }

    let mesh_count = meshes.len();
    let mut class_cache = vec![None; frags.len() * mesh_count];

    let mut result_faces = Vec::with_capacity(frags.len() + coplanar_result_faces.len());

    for (frag_index, frag) in frags.iter().enumerate() {
        let p0 = *pool.position(frag.face.vertices[0]);
        let p1 = *pool.position(frag.face.vertices[1]);
        let p2 = *pool.position(frag.face.vertices[2]);

        let tri = [p0, p1, p2];
        let normal = tri_normal(&tri);
        if is_degenerate_sliver_with_normal(&tri, &normal) {
            continue;
        }

        if lies_on_resolved_coplanar_plane(&tri, &coplanar_plane_infos) {
            continue;
        }

        let centroid = centroid(&tri);
        let mut survive = true;
        for (other_idx, prepared_faces) in prepared_meshes.iter().enumerate() {
            if other_idx == frag.mesh_idx {
                continue;
            }

            let cache_slot = frag_index * mesh_count + other_idx;
            let cached = if let Some(cached) = class_cache[cache_slot] {
                cached
            } else {
                let resolved = CachedClassification {
                    class: classify_fragment_against_mesh(
                        &centroid,
                        &normal,
                        mesh_aabbs[other_idx].contains_point(&centroid),
                        prepared_faces,
                    ),
                };
                class_cache[cache_slot] = Some(resolved);
                resolved
            };

            if !fragment_survives_against_operand(op, frag.mesh_idx, other_idx, cached.class) {
                survive = false;
                break;
            }
        }

        if survive {
            let parent_face = meshes[frag.mesh_idx][frag.parent_idx];
            let face = if op == BooleanOp::Difference && frag.mesh_idx > 0 {
                FaceData::new(
                    frag.face.vertices[0],
                    frag.face.vertices[2],
                    frag.face.vertices[1],
                    parent_face.region,
                )
            } else {
                FaceData::new(
                    frag.face.vertices[0],
                    frag.face.vertices[1],
                    frag.face.vertices[2],
                    parent_face.region,
                )
            };
            result_faces.push(face);
        }
    }

    result_faces.append(&mut coplanar_result_faces);
    result_faces
}

fn lies_on_resolved_coplanar_plane(
    tri: &[Point3r; 3],
    coplanar_plane_infos: &[CoplanarPlaneInfo],
) -> bool {
    coplanar_plane_infos.iter().any(|plane| {
        orient3d(&plane.a, &plane.b, &plane.c, &tri[0]) == Sign::Zero
            && orient3d(&plane.a, &plane.b, &plane.c, &tri[1]) == Sign::Zero
            && orient3d(&plane.a, &plane.b, &plane.c, &tri[2]) == Sign::Zero
    })
}

fn classify_fragment_against_mesh(
    centroid: &Point3r,
    frag_normal: &nalgebra::Vector3<f64>,
    aabb_contains_centroid: bool,
    prepared_faces: &[super::classify::PreparedFace],
) -> FragmentClass {
    if !aabb_contains_centroid || prepared_faces.is_empty() {
        return FragmentClass::Outside;
    }

    classify_fragment_prepared(centroid, frag_normal, prepared_faces)
}

fn fragment_survives_against_operand(
    op: BooleanOp,
    frag_mesh_idx: usize,
    other_idx: usize,
    class: FragmentClass,
) -> bool {
    match op {
        BooleanOp::Union => matches!(class, FragmentClass::Outside)
            || matches!(class, FragmentClass::CoplanarSame) && frag_mesh_idx < other_idx,
        BooleanOp::Intersection => {
            matches!(class, FragmentClass::Inside | FragmentClass::CoplanarSame)
        }
        BooleanOp::Difference => {
            if frag_mesh_idx == 0 {
                return matches!(class, FragmentClass::Outside);
            }

            if other_idx == 0 {
                return matches!(class, FragmentClass::Inside | FragmentClass::CoplanarOpposite);
            }

            matches!(class, FragmentClass::Outside)
                || matches!(class, FragmentClass::CoplanarSame) && frag_mesh_idx < other_idx
        }
    }
}

fn consolidate_cross_mesh_vertices(frags: &mut Vec<BooleanFragmentRecord>, pool: &VertexPool) {
    let tol_sq = 4.0e-8_f64; // (2e-4)^2
    let tol = 2e-4_f64;
    let inv_cell = 1.0 / tol;

    let mut all_vids = Vec::with_capacity(frags.len() * 3);
    for fragment in &*frags {
        for &vertex in &fragment.face.vertices {
            all_vids.push(vertex);
        }
    }
    all_vids.sort_unstable();
    all_vids.dedup();
    if all_vids.is_empty() {
        return;
    }

    let mut grid: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::with_capacity(all_vids.len());
    let positions: Vec<Point3r> = all_vids.iter().map(|&vid| *pool.position(vid)).collect();

    for (index, position) in positions.iter().enumerate() {
        let ix = (position.x * inv_cell).floor() as i64;
        let iy = (position.y * inv_cell).floor() as i64;
        let iz = (position.z * inv_cell).floor() as i64;
        grid.entry((ix, iy, iz)).or_default().push(index);
    }

    let mut parent: Vec<usize> = (0..all_vids.len()).collect();

    fn find_root(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            let p = parent[x];
            let gp = parent[p];
            parent[x] = gp;
            x = gp;
        }
        x
    }

    for (index, position) in positions.iter().enumerate() {
        let ix = (position.x * inv_cell).floor() as i64;
        let iy = (position.y * inv_cell).floor() as i64;
        let iz = (position.z * inv_cell).floor() as i64;

        for dx in -1_i64..=1 {
            for dy in -1_i64..=1 {
                for dz in -1_i64..=1 {
                    if let Some(candidates) = grid.get(&(ix + dx, iy + dy, iz + dz)) {
                        for &other_index in candidates {
                            if index == other_index {
                                continue;
                            }

                            let other_position = &positions[other_index];
                            let d_sq = (other_position.x - position.x).powi(2)
                                + (other_position.y - position.y).powi(2)
                                + (other_position.z - position.z).powi(2);
                            if d_sq < tol_sq {
                                let root_i = find_root(&mut parent, index);
                                let root_j = find_root(&mut parent, other_index);
                                if root_i != root_j {
                                    if root_i < root_j {
                                        parent[root_j] = root_i;
                                    } else {
                                        parent[root_i] = root_j;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    let mut merge_map: HashMap<VertexId, VertexId> = HashMap::with_capacity(all_vids.len() / 2);
    for index in 0..all_vids.len() {
        let root = find_root(&mut parent, index);
        if root != index {
            merge_map.insert(all_vids[index], all_vids[root]);
        }
    }

    if !merge_map.is_empty() {
        for fragment in frags.iter_mut() {
            for vertex in &mut fragment.face.vertices {
                if let Some(&root) = merge_map.get(vertex) {
                    *vertex = root;
                }
            }
        }
    }

    frags.retain(|fragment| {
        let v = fragment.face.vertices;
        v[0] != v[1] && v[1] != v[2] && v[2] != v[0]
    });

    let mut seen = HashSet::with_capacity(frags.len());
    frags.retain(|fragment| {
        let mut key = fragment.face.vertices;
        key.sort();
        seen.insert(key)
    });
}
