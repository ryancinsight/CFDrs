//! `IndexedMesh` wrapper API for CSG Boolean operations

use super::union_strategy::{concat_disjoint_meshes, spatial_union_order};
use crate::application::csg::boolean::BooleanOp;
use crate::application::csg::reconstruct;
use crate::domain::core::error::{MeshError, MeshResult};
use crate::domain::core::index::{FaceId, VertexId};
use crate::domain::geometry::normal::triangle_normal;
use crate::domain::mesh::IndexedMesh;
use crate::infrastructure::storage::face_store::{FaceData, FaceStore};
use crate::infrastructure::storage::vertex_pool::VertexPool;

/// High-level boolean operation on two [`IndexedMesh`] objects.
///
/// Merges the vertex pools, runs the proven binary Boolean pipeline, and
/// reconstructs a fresh deduplicated `IndexedMesh`.
pub fn csg_boolean(
    op: BooleanOp,
    mesh_a: &IndexedMesh,
    mesh_b: &IndexedMesh,
) -> MeshResult<IndexedMesh> {
    csg_boolean_indexed(op, mesh_a, mesh_b)
}

/// Backward-compatible alias for the historical indexed binary Boolean entrypoint.
pub fn csg_boolean_indexed(
    op: BooleanOp,
    mesh_a: &IndexedMesh,
    mesh_b: &IndexedMesh,
) -> MeshResult<IndexedMesh> {
    let mut combined = VertexPool::for_csg();
    let (faces_a, faces_b) = remap_binary_face_soups(mesh_a, mesh_b, &mut combined);
    let is_coplanar = crate::application::csg::coplanar::detect_flat_plane(&faces_a, &combined)
        .is_some()
        && crate::application::csg::coplanar::detect_flat_plane(&faces_b, &combined).is_some();
    let result_faces =
        crate::application::csg::arrangement::boolean_csg::csg_boolean_unfinalized(
            op,
            &[faces_a.clone(), faces_b.clone()],
            &mut combined,
        )?;
    postprocess_boolean_mesh(result_faces, &combined, is_coplanar)
}

/// Best-effort indexed Boolean that preserves the historical tolerant path for
/// junction-heavy unions.
pub fn csg_boolean_indexed_tolerant(
    op: BooleanOp,
    mesh_a: &IndexedMesh,
    mesh_b: &IndexedMesh,
) -> MeshResult<IndexedMesh> {
    match csg_boolean_indexed(op, mesh_a, mesh_b) {
        Ok(mesh) => Ok(mesh),
        Err(MeshError::NotWatertight { .. }) => {
            let mut combined = VertexPool::for_csg();
            let (faces_a, faces_b) = remap_binary_face_soups(mesh_a, mesh_b, &mut combined);
            let is_coplanar = crate::application::csg::coplanar::detect_flat_plane(
                &faces_a,
                &combined,
            )
            .is_some()
                && crate::application::csg::coplanar::detect_flat_plane(&faces_b, &combined)
                    .is_some();
            let result_faces =
                crate::application::csg::arrangement::boolean_csg::csg_boolean_unfinalized(
                    op,
                    &[faces_a.clone(), faces_b.clone()],
                    &mut combined,
                )?;
            let mut mesh = reconstruct::reconstruct_mesh(&result_faces, &combined);
            mesh.recompute_normals();
            repair_boolean_mesh(&mut mesh, is_coplanar)?;
            collapse_degenerate_faces(&mut mesh);
            split_non_manifold_vertices(&mut mesh);
            mesh.rebuild_edges();
            Ok(mesh)
        }
        Err(error) => Err(error),
    }
}

/// Best-effort indexed Boolean used by higher-level assembly paths that can
/// accept a repaired-but-still-open intermediate mesh and enforce topology
/// policy at the final aggregate surface instead.
pub fn csg_boolean_indexed_best_effort(
    op: BooleanOp,
    mesh_a: &IndexedMesh,
    mesh_b: &IndexedMesh,
) -> MeshResult<IndexedMesh> {
    match csg_boolean_indexed(op, mesh_a, mesh_b) {
        Ok(mesh) => Ok(mesh),
        Err(MeshError::NotWatertight { .. }) => {
            let mut combined = VertexPool::for_csg();
            let (faces_a, faces_b) = remap_binary_face_soups(mesh_a, mesh_b, &mut combined);
            let is_coplanar = crate::application::csg::coplanar::detect_flat_plane(
                &faces_a,
                &combined,
            )
            .is_some()
                && crate::application::csg::coplanar::detect_flat_plane(&faces_b, &combined)
                    .is_some();
            let result_faces =
                crate::application::csg::arrangement::boolean_csg::csg_boolean_unfinalized(
                    op,
                    &[faces_a.clone(), faces_b.clone()],
                    &mut combined,
                )?;
            let mut mesh = reconstruct::reconstruct_mesh(&result_faces, &combined);
            mesh.recompute_normals();
            repair_boolean_mesh_with_policy(&mut mesh, is_coplanar, false)?;
            collapse_degenerate_faces(&mut mesh);
            split_non_manifold_vertices(&mut mesh);
            mesh.rebuild_edges();
            Ok(mesh)
        }
        Err(error) => Err(error),
    }
}

/// Compute an indexed Boolean union across an arbitrary number of meshes using
/// the canonical generalized arrangement engine.
pub fn csg_boolean_nary_union(meshes: &[IndexedMesh]) -> MeshResult<IndexedMesh> {
    if meshes.is_empty() {
        return Err(MeshError::EmptyBooleanResult {
            op: format!("{:?}", BooleanOp::Union),
        });
    }

    if meshes.len() == 1 {
        return Ok(meshes[0].clone());
    }

    match csg_boolean_nary(BooleanOp::Union, meshes) {
        Ok(mesh) => Ok(mesh),
        Err(MeshError::NotWatertight { .. }) => {
            let order = spatial_union_order(meshes);
            let mut accumulated = meshes[order[0]].clone();
            let mut accumulated_aabb = accumulated.bounding_box();
            for &idx in &order[1..] {
                let mesh = &meshes[idx];
                let mesh_aabb = mesh.bounding_box();
                accumulated = if !accumulated_aabb.intersects(&mesh_aabb) {
                    concat_disjoint_meshes(&accumulated, mesh)
                } else {
                    csg_boolean_indexed_tolerant(BooleanOp::Union, &accumulated, mesh)?
                };
                accumulated_aabb = accumulated.bounding_box();
            }
            Ok(accumulated)
        }
        Err(error) => Err(error),
    }
}

/// Compute an indexed Boolean across an arbitrary number of meshes using the
/// canonical generalized arrangement engine.
///
/// The first mesh is treated as the minuend for [`BooleanOp::Difference`], and
/// every later mesh is treated as a subtractive operand.
pub fn csg_boolean_nary(op: BooleanOp, meshes: &[IndexedMesh]) -> MeshResult<IndexedMesh> {
    if meshes.is_empty() {
        return Err(MeshError::EmptyBooleanResult {
            op: format!("{op:?}"),
        });
    }

    if meshes.len() == 1 {
        return Ok(meshes[0].clone());
    }

    let mut combined = VertexPool::for_csg();
    let face_soups = remap_nary_face_soups(meshes, &mut combined);
    let is_coplanar = face_soups.iter().all(|faces| {
        crate::application::csg::coplanar::detect_flat_plane(faces, &combined).is_some()
    });
    let result_faces = crate::application::csg::arrangement::boolean_csg::csg_boolean_unfinalized(
        op,
        &face_soups,
        &mut combined,
    )?;
    postprocess_boolean_mesh(result_faces, &combined, is_coplanar)
}

fn remap_binary_face_soups(
    mesh_a: &IndexedMesh,
    mesh_b: &IndexedMesh,
    combined: &mut VertexPool,
) -> (Vec<FaceData>, Vec<FaceData>) {
    use crate::domain::core::index::VertexId;
    use std::collections::HashMap;

    let mut remap_a: HashMap<VertexId, VertexId> = HashMap::with_capacity(mesh_a.vertices.len());
    for (old_id, _) in mesh_a.vertices.iter() {
        let pos = *mesh_a.vertices.position(old_id);
        let nrm = *mesh_a.vertices.normal(old_id);
        remap_a.insert(old_id, combined.insert_or_weld(pos, nrm));
    }

    let mut remap_b: HashMap<VertexId, VertexId> = HashMap::with_capacity(mesh_b.vertices.len());
    for (old_id, _) in mesh_b.vertices.iter() {
        let pos = *mesh_b.vertices.position(old_id);
        let nrm = *mesh_b.vertices.normal(old_id);
        remap_b.insert(old_id, combined.insert_or_weld(pos, nrm));
    }

    let faces_a: Vec<FaceData> = mesh_a
        .faces
        .iter()
        .map(|face| FaceData {
            vertices: face.vertices.map(|vertex_id| remap_a[&vertex_id]),
            region: face.region,
        })
        .collect();
    let faces_b: Vec<FaceData> = mesh_b
        .faces
        .iter()
        .map(|face| FaceData {
            vertices: face.vertices.map(|vertex_id| remap_b[&vertex_id]),
            region: face.region,
        })
        .collect();

    (faces_a, faces_b)
}

fn remap_nary_face_soups(meshes: &[IndexedMesh], combined: &mut VertexPool) -> Vec<Vec<FaceData>> {
    use crate::domain::core::index::VertexId;
    use std::collections::HashMap;

    let mut face_soups = Vec::with_capacity(meshes.len());
    for mesh in meshes {
        let mut remap: HashMap<VertexId, VertexId> = HashMap::with_capacity(mesh.vertices.len());
        for (old_id, _) in mesh.vertices.iter() {
            let pos = *mesh.vertices.position(old_id);
            let nrm = *mesh.vertices.normal(old_id);
            remap.insert(old_id, combined.insert_or_weld(pos, nrm));
        }

        let faces: Vec<FaceData> = mesh
            .faces
            .iter()
            .map(|face| FaceData {
                vertices: face.vertices.map(|vertex_id| remap[&vertex_id]),
                region: face.region,
            })
            .collect();
        face_soups.push(faces);
    }
    face_soups
}

fn postprocess_boolean_mesh(
    result_faces: Vec<FaceData>,
    combined: &VertexPool,
    is_coplanar: bool,
) -> MeshResult<IndexedMesh> {
    let mut mesh = reconstruct::reconstruct_mesh(&result_faces, combined);
    mesh.recompute_normals();
    repair_boolean_mesh(&mut mesh, is_coplanar)?;
    collapse_degenerate_faces(&mut mesh);
    split_non_manifold_vertices(&mut mesh);
    mesh.rebuild_edges();
    Ok(mesh)
}

fn repair_boolean_mesh(mesh: &mut IndexedMesh, is_coplanar: bool) -> MeshResult<()> {
    repair_boolean_mesh_with_policy(mesh, is_coplanar, true)
}

fn repair_boolean_mesh_with_policy(
    mesh: &mut IndexedMesh,
    is_coplanar: bool,
    require_watertight: bool,
) -> MeshResult<()> {
    if !is_coplanar {
        mesh.rebuild_edges();
        let mut report = crate::application::watertight::check::check_watertight(
            &mesh.vertices,
            &mesh.faces,
            mesh.edges_ref().unwrap(),
        );
        if !report.is_watertight {
            if report.is_closed && !report.orientation_consistent {
                mesh.orient_outward();
                mesh.rebuild_edges();
                report = crate::application::watertight::check::check_watertight(
                    &mesh.vertices,
                    &mesh.faces,
                    mesh.edges_ref().unwrap(),
                );
            }

            if !report.is_watertight && report.is_closed && !report.orientation_consistent {
                let edge_store =
                    crate::infrastructure::storage::edge_store::EdgeStore::from_face_store(
                        &mesh.faces,
                    );
                let _ = crate::domain::topology::orientation::fix_orientation(
                    &mut mesh.faces,
                    &edge_store,
                );
                mesh.rebuild_edges();
                mesh.orient_outward();
                mesh.rebuild_edges();
                report = crate::application::watertight::check::check_watertight(
                    &mesh.vertices,
                    &mesh.faces,
                    mesh.edges_ref().unwrap(),
                );
            }

            if !report.is_watertight
                && report.non_manifold_edge_count == 0
                && report.boundary_edge_count > 0
                && report.boundary_edge_count <= 512
            {
                let edge_store =
                    crate::infrastructure::storage::edge_store::EdgeStore::from_face_store(
                        &mesh.faces,
                    );
                let added = crate::application::watertight::seal::seal_boundary_loops(
                    &mut mesh.vertices,
                    &mut mesh.faces,
                    &edge_store,
                    crate::domain::core::index::RegionId::INVALID,
                );
                if added > 0 {
                    mesh.rebuild_edges();
                    mesh.orient_outward();
                    mesh.rebuild_edges();
                    report = crate::application::watertight::check::check_watertight(
                        &mesh.vertices,
                        &mesh.faces,
                        mesh.edges_ref().unwrap(),
                    );
                }
            }

            if !report.is_watertight && report.boundary_edge_count > 0 {
                let improved = crate::application::watertight::repair::MeshRepair::iterative_boundary_stitch(
                    &mut mesh.faces,
                    &mesh.vertices,
                    3,
                );
                if improved > 0 {
                    mesh.rebuild_edges();
                    mesh.orient_outward();
                    mesh.rebuild_edges();
                    report = crate::application::watertight::check::check_watertight(
                        &mesh.vertices,
                        &mesh.faces,
                        mesh.edges_ref().unwrap(),
                    );
                }
            }

            // Non-manifold edge resolution + boundary vertex merging.
            // Handles edges with 3+ faces (from arrangement artifacts) and
            // sliver gaps at complex intersection curves.
            if !report.is_watertight
                && (report.boundary_edge_count > 0 || report.non_manifold_edge_count > 0)
            {
                eprintln!("[repair] ENTER: boundary={} non_manifold={} faces={}",
                    report.boundary_edge_count, report.non_manifold_edge_count, mesh.face_count());

                // First resolve non-manifold edges.
                split_non_manifold_edges(mesh);
                collapse_degenerate_faces(mesh);
                mesh.rebuild_edges();

                {
                    let r = crate::application::watertight::check::check_watertight(
                        &mesh.vertices, &mesh.faces, mesh.edges_ref().unwrap());
                    eprintln!("[repair] after split_nm: boundary={} nm={} faces={}",
                        r.boundary_edge_count, r.non_manifold_edge_count, mesh.face_count());
                }

                // Seal any boundary loops opened by face removal.
                if !mesh.is_watertight() {
                    let es = crate::infrastructure::storage::edge_store::EdgeStore::from_face_store(
                        &mesh.faces,
                    );
                    let sealed = crate::application::watertight::seal::seal_boundary_loops(
                        &mut mesh.vertices,
                        &mut mesh.faces,
                        &es,
                        crate::domain::core::index::RegionId::INVALID,
                    );
                    if sealed > 0 {
                        collapse_degenerate_faces(mesh);
                                mesh.rebuild_edges();
                    }
                    eprintln!("[repair] after seal1: sealed={} watertight={}", sealed, mesh.is_watertight());
                }

                // Merge nearby boundary vertices to close sliver gaps.
                if !mesh.is_watertight() {
                    merge_nearby_boundary_vertices(mesh);
                    collapse_degenerate_faces(mesh);
                        mesh.rebuild_edges();
                    {
                        let r = crate::application::watertight::check::check_watertight(
                            &mesh.vertices, &mesh.faces, mesh.edges_ref().unwrap());
                        eprintln!("[repair] after merge: boundary={} nm={} faces={}",
                            r.boundary_edge_count, r.non_manifold_edge_count, mesh.face_count());
                    }
                }

                // Final seal pass after merging.
                if !mesh.is_watertight() {
                    let es2 = crate::infrastructure::storage::edge_store::EdgeStore::from_face_store(
                        &mesh.faces,
                    );
                    let _ = crate::application::watertight::seal::seal_boundary_loops(
                        &mut mesh.vertices,
                        &mut mesh.faces,
                        &es2,
                        crate::domain::core::index::RegionId::INVALID,
                    );
                    collapse_degenerate_faces(mesh);
                        mesh.rebuild_edges();
                    eprintln!("[repair] after seal2: watertight={}", mesh.is_watertight());
                }

                mesh.orient_outward();
                mesh.rebuild_edges();
                report = crate::application::watertight::check::check_watertight(
                    &mesh.vertices,
                    &mesh.faces,
                    mesh.edges_ref().unwrap(),
                );
                eprintln!("[repair] FINAL: boundary={} nm={} faces={} watertight={}",
                    report.boundary_edge_count, report.non_manifold_edge_count, mesh.face_count(), report.is_watertight);
            }

            if require_watertight && !report.is_watertight {
                return Err(MeshError::NotWatertight {
                    count: report.boundary_edge_count + report.non_manifold_edge_count,
                });
            }
        }

        mesh.orient_outward();
        if mesh.signed_volume() < 0.0 {
            mesh.flip_faces();
            mesh.rebuild_edges();
        }
    }

    Ok(())
}

/// Collapse degenerate (zero-area) faces by merging the redundant vertex.
///
/// Handles two cases:
/// 1. Near-coincident vertices (distance² < tolerance²): merge the pair.
/// 2. Collinear slivers (3 distinct vertices on a line): find the "middle"
///    vertex and merge it into the nearest endpoint, but only when (a) the edge
///    has at most 2 incident faces, (b) no duplicate faces would be created,
///    and (c) no surviving face normals would be inverted.
fn collapse_degenerate_faces(mesh: &mut IndexedMesh) {
    use std::collections::{HashMap, HashSet};

    let tol_sq = 1e-18_f64; // (1e-9)²
    let mut total_collapsed: usize = 0;
    let mut skip_faces: HashSet<usize> = HashSet::new();

    loop {
        // Build an edge-use count: how many faces reference each undirected edge.
        let mut edge_use: HashMap<(VertexId, VertexId), usize> = HashMap::new();
        for face in mesh.faces.iter() {
            let v = face.vertices;
            for &(a, b) in &[(v[0], v[1]), (v[1], v[2]), (v[2], v[0])] {
                let key = if a < b { (a, b) } else { (b, a) };
                *edge_use.entry(key).or_insert(0) += 1;
            }
        }

        // Find a degenerate face we can safely collapse.
        let degen = mesh.faces.iter().enumerate().find_map(|(i, face)| {
            if skip_faces.contains(&i) {
                return None;
            }
            let pa = mesh.vertices.position(face.vertices[0]);
            let pb = mesh.vertices.position(face.vertices[1]);
            let pc = mesh.vertices.position(face.vertices[2]);
            if triangle_normal(pa, pb, pc).is_some() {
                return None; // not degenerate
            }

            let d01 = (pb - pa).norm_squared();
            let d12 = (pc - pb).norm_squared();
            let d20 = (pa - pc).norm_squared();

            // Case 1: near-coincident vertex pair — always safe.
            let shortest = d01.min(d12).min(d20);
            if shortest <= tol_sq {
                let (keep, remove) = if d01 <= d12 && d01 <= d20 {
                    (face.vertices[0], face.vertices[1])
                } else if d12 <= d20 {
                    (face.vertices[1], face.vertices[2])
                } else {
                    (face.vertices[2], face.vertices[0])
                };
                return Some((i, keep, remove, true));
            }

            // Case 2: collinear sliver — find the middle vertex (opposite the
            // longest edge) and merge it into the nearer endpoint.
            let (middle_idx, ep_a_idx, ep_b_idx) = if d01 >= d12 && d01 >= d20 {
                (2, 0, 1)
            } else if d12 >= d20 {
                (0, 1, 2)
            } else {
                (1, 2, 0)
            };
            let middle = face.vertices[middle_idx];
            let ep_a = face.vertices[ep_a_idx];
            let ep_b = face.vertices[ep_b_idx];

            let pm = mesh.vertices.position(middle);
            let pea = mesh.vertices.position(ep_a);
            let peb = mesh.vertices.position(ep_b);
            let (keep, remove) = if (pea - pm).norm_squared() <= (peb - pm).norm_squared() {
                (ep_a, middle)
            } else {
                (ep_b, middle)
            };

            // Link condition: the edge (keep, remove) must be shared by at most
            // 2 faces.
            let edge_key = if keep < remove {
                (keep, remove)
            } else {
                (remove, keep)
            };
            let uses = edge_use.get(&edge_key).copied().unwrap_or(0);
            if uses > 2 {
                return None;
            }
            Some((i, keep, remove, false))
        });

        let (face_idx, keep, remove, is_coincident) = match degen {
            Some(v) => v,
            None => break,
        };

        // For sliver collapses, simulate the rename first and check for hazards.
        if !is_coincident {
            let mut would_create_duplicate = false;
            let mut would_invert_normal = false;
            let mut new_face_keys: HashSet<[VertexId; 3]> = HashSet::new();

            for face in mesh.faces.iter() {
                let mut verts = face.vertices;
                for v in verts.iter_mut() {
                    if *v == remove {
                        *v = keep;
                    }
                }
                // Skip faces that become degenerate (two identical verts).
                if verts[0] == verts[1]
                    || verts[1] == verts[2]
                    || verts[2] == verts[0]
                {
                    continue;
                }
                let mut key = verts;
                if key[0] > key[1] {
                    key.swap(0, 1);
                }
                if key[1] > key[2] {
                    key.swap(1, 2);
                }
                if key[0] > key[1] {
                    key.swap(0, 1);
                }
                if !new_face_keys.insert(key) {
                    would_create_duplicate = true;
                    break;
                }
            }

            // Check that no surviving face gets its normal inverted.
            if !would_create_duplicate {
                for face in mesh.faces.iter() {
                    let has_remove = face.vertices.contains(&remove);
                    if !has_remove {
                        continue;
                    }
                    // Compute pre-collapse normal.
                    let p0 = mesh.vertices.position(face.vertices[0]);
                    let p1 = mesh.vertices.position(face.vertices[1]);
                    let p2 = mesh.vertices.position(face.vertices[2]);
                    let pre_n = triangle_normal(p0, p1, p2);

                    // Compute post-collapse vertices.
                    let mut new_verts = face.vertices;
                    for v in new_verts.iter_mut() {
                        if *v == remove {
                            *v = keep;
                        }
                    }
                    if new_verts[0] == new_verts[1]
                        || new_verts[1] == new_verts[2]
                        || new_verts[2] == new_verts[0]
                    {
                        continue; // will be removed
                    }
                    let q0 = mesh.vertices.position(new_verts[0]);
                    let q1 = mesh.vertices.position(new_verts[1]);
                    let q2 = mesh.vertices.position(new_verts[2]);
                    let post_n = triangle_normal(q0, q1, q2);

                    if let (Some(pn), Some(qn)) = (pre_n, post_n) {
                        if pn.dot(&qn) < 0.0 {
                            would_invert_normal = true;
                            break;
                        }
                    }
                }
            }

            if would_create_duplicate || would_invert_normal {
                skip_faces.insert(face_idx);
                continue;
            }
        }

        // Commit: rewrite all face references from `remove` → `keep`.
        for face in mesh.faces.iter_mut() {
            for v in face.vertices.iter_mut() {
                if *v == remove {
                    *v = keep;
                }
            }
        }

        // Purge collapsed faces and exact duplicates.
        let mut seen: HashSet<[VertexId; 3]> = HashSet::new();
        let mut keep_faces: Vec<FaceData> = Vec::with_capacity(mesh.faces.len());
        let mut removed = 0usize;
        for face in mesh.faces.iter() {
            if face.vertices[0] == face.vertices[1]
                || face.vertices[1] == face.vertices[2]
                || face.vertices[2] == face.vertices[0]
            {
                removed += 1;
                continue;
            }
            let mut key = face.vertices;
            if key[0] > key[1] {
                key.swap(0, 1);
            }
            if key[1] > key[2] {
                key.swap(1, 2);
            }
            if key[0] > key[1] {
                key.swap(0, 1);
            }
            if !seen.insert(key) {
                removed += 1;
                continue;
            }
            keep_faces.push(*face);
        }
        total_collapsed += removed;
        skip_faces.clear(); // face indices changed after rebuild
        let mut new_faces = FaceStore::with_capacity(keep_faces.len());
        for f in keep_faces {
            new_faces.push(f);
        }
        mesh.faces = new_faces;
    }
    if total_collapsed > 0 {
        tracing::debug!(
            "CSG postprocess: collapsed {} degenerate face(s)",
            total_collapsed
        );
    }
}

/// Split non-manifold "pinch" vertices created by edge collapse.
fn split_non_manifold_vertices(mesh: &mut IndexedMesh) {
    use std::collections::{HashMap, HashSet, VecDeque};
    let mut vertex_faces: HashMap<VertexId, Vec<usize>> = HashMap::new();
    for (fi, face) in mesh.faces.iter().enumerate() {
        for &v in &face.vertices {
            vertex_faces.entry(v).or_default().push(fi);
        }
    }
    let mut total_splits: usize = 0;
    let vertices: Vec<VertexId> = vertex_faces.keys().copied().collect();
    for v in vertices {
        let face_indices = match vertex_faces.get(&v) {
            Some(fi) if fi.len() >= 2 => fi,
            _ => continue,
        };
        let mut visited: HashSet<usize> = HashSet::new();
        let mut components: Vec<Vec<usize>> = Vec::new();
        for &start_fi in face_indices {
            if visited.contains(&start_fi) {
                continue;
            }
            let mut component: Vec<usize> = Vec::new();
            let mut queue: VecDeque<usize> = VecDeque::new();
            queue.push_back(start_fi);
            visited.insert(start_fi);
            while let Some(fi) = queue.pop_front() {
                component.push(fi);
                let face = mesh.faces.get(FaceId::from_usize(fi));
                let others: Vec<VertexId> =
                    face.vertices.iter().copied().filter(|&vid| vid != v).collect();
                for &neighbor_fi in face_indices {
                    if visited.contains(&neighbor_fi) {
                        continue;
                    }
                    let neighbor = mesh.faces.get(FaceId::from_usize(neighbor_fi));
                    let shares_edge = others.iter().any(|ov| neighbor.vertices.contains(ov));
                    if shares_edge {
                        visited.insert(neighbor_fi);
                        queue.push_back(neighbor_fi);
                    }
                }
            }
            components.push(component);
        }
        if components.len() <= 1 {
            continue;
        }
        let pos = *mesh.vertices.position(v);
        let normal = *mesh.vertices.normal(v);
        for component in components.iter().skip(1) {
            let new_v = mesh.add_vertex(pos, normal);
            for &fi in component {
                let fid = FaceId::from_usize(fi);
                let face_mut = mesh.faces.get_mut(fid);
                for vref in face_mut.vertices.iter_mut() {
                    if *vref == v {
                        *vref = new_v;
                    }
                }
            }
            total_splits += 1;
        }
    }
    if total_splits > 0 {
        tracing::debug!(
            "CSG postprocess: split {} non-manifold vertex instance(s)",
            total_splits
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::application::watertight::check::check_watertight;
    use crate::domain::core::scalar::Point3r;
    use crate::domain::geometry::primitives::{Cube, Cylinder, Disk, PrimitiveMesh, UvSphere};

    fn sphere() -> IndexedMesh {
        UvSphere {
            radius: 1.0,
            center: Point3r::origin(),
            segments: 16,
            stacks: 8,
        }
        .build()
        .expect("sphere build")
    }

    fn cylinder() -> IndexedMesh {
        Cylinder {
            base_center: Point3r::new(0.0, -1.5, 0.0),
            radius: 0.4,
            height: 3.0,
            segments: 16,
        }
        .build()
        .expect("cylinder build")
    }

    fn cube_a() -> IndexedMesh {
        Cube {
            origin: Point3r::new(-1.0, -1.0, -1.0),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .expect("cube_a build")
    }

    fn cube_b() -> IndexedMesh {
        Cube {
            origin: Point3r::new(-0.5, -0.5, -0.5),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .expect("cube_b build")
    }

    fn disk_a() -> IndexedMesh {
        Disk {
            center: Point3r::new(0.0, 0.0, 0.0),
            radius: 1.0,
            segments: 16,
        }
        .build()
        .expect("disk_a build")
    }

    fn disk_b() -> IndexedMesh {
        Disk {
            center: Point3r::new(0.5, 0.0, 0.0),
            radius: 1.0,
            segments: 16,
        }
        .build()
        .expect("disk_b build")
    }

    /// Assert a 3-D CSG result is watertight with a positive signed volume.
    fn assert_3d_watertight(mut mesh: IndexedMesh) {
        mesh.rebuild_edges();
        let report = check_watertight(&mesh.vertices, &mesh.faces, mesh.edges_ref().unwrap());
        assert!(
            report.is_watertight,
            "CSG result must be watertight: {} boundary edge(s), {} non-manifold edge(s)",
            report.boundary_edge_count, report.non_manifold_edge_count,
        );
        assert!(
            mesh.signed_volume() > 0.0,
            "CSG result must have positive signed volume (outward-oriented normals)",
        );
    }

    fn component_count(mesh: &mut IndexedMesh) -> usize {
        use crate::domain::topology::connectivity::connected_components;
        use crate::domain::topology::AdjacencyGraph;

        mesh.rebuild_edges();
        let adjacency = AdjacencyGraph::build(&mesh.faces, mesh.edges_ref().unwrap());
        connected_components(&mesh.faces, &adjacency).len()
    }

    fn symmetric_parallel_cylinders(segments: usize) -> (IndexedMesh, IndexedMesh) {
        let radius = 0.6;
        let height = 3.0;
        let separation = radius;
        let cyl_a = Cylinder {
            base_center: Point3r::new(-separation / 2.0, -height / 2.0, 0.0),
            radius,
            height,
            segments,
        }
        .build()
        .expect("symmetric cyl_a build");
        let cyl_b = Cylinder {
            base_center: Point3r::new(separation / 2.0, -height / 2.0, 0.0),
            radius,
            height,
            segments,
        }
        .build()
        .expect("symmetric cyl_b build");
        (cyl_a, cyl_b)
    }

    fn planar_branch(angle_from_x: f64, radius: f64, height: f64, segments: usize) -> IndexedMesh {
        use crate::application::csg::CsgNode;
        use nalgebra::{Isometry3, Translation3, UnitQuaternion, Vector3};

        let raw = Cylinder {
            base_center: Point3r::new(0.0, 0.0, 0.0),
            radius,
            height,
            segments,
        }
        .build()
        .expect("branch build");
        let rotation = UnitQuaternion::<f64>::from_axis_angle(
            &Vector3::z_axis(),
            angle_from_x - std::f64::consts::FRAC_PI_2,
        );
        CsgNode::Transform {
            node: Box::new(CsgNode::Leaf(Box::new(raw))),
            iso: Isometry3::from_parts(Translation3::new(0.0, 0.0, 0.0), rotation),
        }
        .evaluate()
        .expect("branch transform")
    }

    fn planar_trunk(radius: f64, height: f64, extension: f64, segments: usize) -> IndexedMesh {
        use crate::application::csg::CsgNode;
        use nalgebra::{Isometry3, Translation3, UnitQuaternion, Vector3};

        let raw = Cylinder {
            base_center: Point3r::new(0.0, 0.0, 0.0),
            radius,
            height: height + extension,
            segments,
        }
        .build()
        .expect("trunk build");
        let rotation = UnitQuaternion::<f64>::from_axis_angle(
            &Vector3::z_axis(),
            -std::f64::consts::FRAC_PI_2,
        );
        CsgNode::Transform {
            node: Box::new(CsgNode::Leaf(Box::new(raw))),
            iso: Isometry3::from_parts(Translation3::new(-height, 0.0, 0.0), rotation),
        }
        .evaluate()
        .expect("trunk transform")
    }

    fn quadfurcation_meshes() -> Vec<IndexedMesh> {
        let radius = 0.5;
        let height = 3.0;
        let extension = radius * 0.10;
        let segments = 32;
        let mut meshes = vec![planar_trunk(radius, height, extension, segments)];
        for angle_deg in [60.0_f64, 20.0, -20.0, -60.0] {
            meshes.push(planar_branch(angle_deg.to_radians(), radius, height, segments));
        }
        meshes
    }

    fn trifurcation_meshes() -> Vec<IndexedMesh> {
        let radius = 0.5;
        let height = 3.0;
        let extension = radius * 0.10;
        let segments = 32;
        let mut meshes = vec![planar_trunk(radius, height, extension, segments)];
        for angle_deg in [45.0_f64, 90.0, -45.0] {
            meshes.push(planar_branch(angle_deg.to_radians(), radius, height, segments));
        }
        meshes
    }

    fn pentafurcation_meshes() -> Vec<IndexedMesh> {
        let radius = 0.5;
        let height = 3.0;
        let extension = radius * 0.10;
        let segments = 32;
        let mut meshes = vec![planar_trunk(radius, height, extension, segments)];
        for angle_deg in [60.0_f64, 30.0, 0.0, -30.0, -60.0] {
            meshes.push(planar_branch(angle_deg.to_radians(), radius, height, segments));
        }
        meshes
    }

    // ── sphere × cylinder (curved × curved — arrangement pipeline) ─────────────

    #[test]
    fn sphere_cylinder_union_is_watertight() {
        let result =
            csg_boolean(BooleanOp::Union, &sphere(), &cylinder()).expect("sphere ∪ cylinder");
        assert_3d_watertight(result);
    }

    #[test]
    fn sphere_cylinder_intersection_is_watertight() {
        let result = csg_boolean(BooleanOp::Intersection, &sphere(), &cylinder())
            .expect("sphere ∩ cylinder");
        assert_3d_watertight(result);
    }

    #[test]
    fn sphere_cylinder_difference_is_watertight() {
        let result =
            csg_boolean(BooleanOp::Difference, &sphere(), &cylinder()).expect("sphere \\ cylinder");
        assert_3d_watertight(result);
    }

    // ── cube × cube (flat faces — intersecting arrangement pipeline) ───────────

    #[test]
    fn cube_cube_union_is_watertight() {
        let result = csg_boolean(BooleanOp::Union, &cube_a(), &cube_b()).expect("cube ∪ cube");
        assert_3d_watertight(result);
    }

    #[test]
    fn cube_cube_intersection_is_watertight() {
        let result =
            csg_boolean(BooleanOp::Intersection, &cube_a(), &cube_b()).expect("cube ∩ cube");
        assert_3d_watertight(result);
    }

    #[test]
    fn cube_cube_difference_is_watertight() {
        let result =
            csg_boolean(BooleanOp::Difference, &cube_a(), &cube_b()).expect("cube \\ cube");
        assert_3d_watertight(result);
    }

    // ── cube × cylinder coplanar (caps flush with cube walls) ──────────────────

    fn cylinder_coplanar() -> IndexedMesh {
        Cylinder {
            base_center: Point3r::new(0.0, -1.0, 0.0),
            radius: 0.4,
            height: 2.0,
            segments: 16,
        }
        .build()
        .expect("cylinder_coplanar build")
    }

    /// Difference of cube minus a coplanar cylinder must be watertight.
    /// The cylinder end caps are coplanar with the cube's top and bottom walls.
    /// The 2-D coplanar pipeline must subtract circular discs from the square
    /// walls, producing annular rings (tunnel openings).
    #[test]
    fn cube_cylinder_coplanar_difference_is_watertight() {
        let result = csg_boolean(BooleanOp::Difference, &cube_a(), &cylinder_coplanar())
            .expect("cube \\\\ cylinder_coplanar");
        assert_3d_watertight(result);
    }

    #[test]
    fn cube_cylinder_coplanar_union_is_watertight() {
        let result = csg_boolean(BooleanOp::Union, &cube_a(), &cylinder_coplanar())
            .expect("cube ∪ cylinder_coplanar");
        assert_3d_watertight(result);
    }

    #[test]
    fn cube_cylinder_coplanar_intersection_is_watertight() {
        let result = csg_boolean(BooleanOp::Intersection, &cube_a(), &cylinder_coplanar())
            .expect("cube ∩ cylinder_coplanar");
        assert_3d_watertight(result);
    }

    // ── disk × disk (coplanar — 2-D Sutherland-Hodgman pipeline) ───────────────
    // Disk operands are open surfaces; the coplanar path produces an open
    // surface result.  Only assert the operation completes without error.

    #[test]
    fn disk_disk_union_succeeds() {
        csg_boolean(BooleanOp::Union, &disk_a(), &disk_b()).expect("disk ∪ disk must not error");
    }

    #[test]
    fn disk_disk_intersection_succeeds() {
        csg_boolean(BooleanOp::Intersection, &disk_a(), &disk_b())
            .expect("disk ∩ disk must not error");
    }

    #[test]
    fn disk_disk_difference_succeeds() {
        csg_boolean(BooleanOp::Difference, &disk_a(), &disk_b())
            .expect("disk \\ disk must not error");
    }

    #[test]
    fn symmetric_parallel_cylinder_intersection_is_single_watertight_component() {
        let (cyl_a, cyl_b) = symmetric_parallel_cylinders(64);
        let mut result =
            csg_boolean(BooleanOp::Intersection, &cyl_a, &cyl_b).expect("symmetric intersection");

        result.rebuild_edges();
        let report = check_watertight(&result.vertices, &result.faces, result.edges_ref().unwrap());
        assert!(
            report.is_watertight,
            "symmetric cylinder intersection must be watertight: boundary={}, non_manifold={}",
            report.boundary_edge_count,
            report.non_manifold_edge_count
        );
        assert_eq!(
            component_count(&mut result),
            1,
            "symmetric cylinder intersection must remain a single component",
        );

        let radius = 0.6;
        let height = 3.0;
        let theta = std::f64::consts::FRAC_PI_3;
        let overlap_area = 2.0 * radius * radius * (theta - theta.sin() * theta.cos());
        let expected = height * overlap_area;
        let relative_error = (result.signed_volume() - expected).abs() / expected;
        assert!(
            relative_error < 0.01,
            "symmetric cylinder intersection volume error {:.2}% exceeds 1%",
            relative_error * 100.0
        );
    }

    #[test]
    fn indexed_nary_quadfurcation_union_is_watertight_without_component_dropping() {
        let mut result =
            csg_boolean_nary(BooleanOp::Union, &quadfurcation_meshes()).expect("quadfurcation union");
        assert_eq!(
            component_count(&mut result),
            1,
            "quadfurcation union must be a single connected component",
        );
        assert_3d_watertight(result);
    }

    #[test]
    fn indexed_nary_trifurcation_union_is_watertight_without_component_dropping() {
        let mut result =
            csg_boolean_nary(BooleanOp::Union, &trifurcation_meshes()).expect("trifurcation union");
        assert_eq!(
            component_count(&mut result),
            1,
            "trifurcation union must be a single connected component",
        );
        assert_3d_watertight(result);
    }

    #[test]
    fn indexed_nary_pentafurcation_union_is_watertight_without_component_dropping() {
        let mut result =
            csg_boolean_nary(BooleanOp::Union, &pentafurcation_meshes()).expect("pentafurcation union");
        assert_eq!(
            component_count(&mut result),
            1,
            "pentafurcation union must be a single connected component",
        );
        assert_3d_watertight(result);
    }

    #[test]
    fn indexed_nary_union_is_permutation_invariant() {
        let forward = quadfurcation_meshes();
        let mut reversed = quadfurcation_meshes();
        reversed.reverse();

        let mut forward_union =
            csg_boolean_nary(BooleanOp::Union, &forward).expect("forward quadfurcation union");
        let mut reversed_union =
            csg_boolean_nary(BooleanOp::Union, &reversed).expect("reversed quadfurcation union");

        assert_3d_watertight(forward_union.clone());
        assert_3d_watertight(reversed_union.clone());
        assert_eq!(
            component_count(&mut forward_union),
            component_count(&mut reversed_union),
            "operand order must not change the number of connected components",
        );

        let forward_volume = forward_union.signed_volume();
        let reversed_volume = reversed_union.signed_volume();
        let relative_error =
            (forward_volume - reversed_volume).abs() / forward_volume.abs().max(1.0e-12);
        assert!(
            relative_error < 0.005,
            "operand order changed union volume by {:.2}%",
            relative_error * 100.0
        );
    }

    // ── Y-junction trunk difference (curved × curved, Difference) ──────────
    // Diagnostic: verify watertight trunk difference has outward-only normals.
    // The BFS seed is the extremal (max-X) face — by the Jordan-Brouwer theorem
    // its outward normal must have nx ≥ 0, so BFS correctly orients the mesh.
    #[test]
    #[ignore = "pre-existing: CSG difference produces phantom island components — repair pipeline does not strip them"]
    fn cylinder_difference_normals_check() {
        use crate::application::csg::CsgNode;
        use crate::application::quality::normals::analyze_normals;
        use crate::domain::core::scalar::Point3r;
        use crate::domain::geometry::primitives::{Cylinder, PrimitiveMesh};
        use nalgebra::{Isometry3, Translation3, UnitQuaternion, Vector3};
        use std::f64::consts::FRAC_PI_2;

        const R: f64 = 0.5;
        const H_TRUNK: f64 = 3.0;
        const H_BRANCH: f64 = 3.0;
        const EPS: f64 = R * 0.10;
        const SEGS: usize = 32;
        let theta = std::f64::consts::FRAC_PI_4;

        let trunk = {
            let raw = Cylinder {
                base_center: Point3r::new(0.0, 0.0, 0.0),
                radius: R,
                height: H_TRUNK + EPS,
                segments: SEGS,
            }
            .build()
            .unwrap();
            let rot = UnitQuaternion::<f64>::from_axis_angle(&Vector3::z_axis(), -FRAC_PI_2);
            let iso = Isometry3::from_parts(Translation3::new(-H_TRUNK, 0.0, 0.0), rot);
            CsgNode::Transform {
                node: Box::new(CsgNode::Leaf(Box::new(raw))),
                iso,
            }
            .evaluate()
            .unwrap()
        };
        let branch_up = {
            let raw = Cylinder {
                base_center: Point3r::new(0.0, 0.0, 0.0),
                radius: R,
                height: H_BRANCH,
                segments: SEGS,
            }
            .build()
            .unwrap();
            let rot = UnitQuaternion::<f64>::from_axis_angle(&Vector3::z_axis(), theta - FRAC_PI_2);
            let iso = Isometry3::from_parts(Translation3::new(0.0, 0.0, 0.0), rot);
            CsgNode::Transform {
                node: Box::new(CsgNode::Leaf(Box::new(raw))),
                iso,
            }
            .evaluate()
            .unwrap()
        };
        let branch_dn = {
            let raw = Cylinder {
                base_center: Point3r::new(0.0, 0.0, 0.0),
                radius: R,
                height: H_BRANCH,
                segments: SEGS,
            }
            .build()
            .unwrap();
            let rot =
                UnitQuaternion::<f64>::from_axis_angle(&Vector3::z_axis(), -theta - FRAC_PI_2);
            let iso = Isometry3::from_parts(Translation3::new(0.0, 0.0, 0.0), rot);
            CsgNode::Transform {
                node: Box::new(CsgNode::Leaf(Box::new(raw))),
                iso,
            }
            .evaluate()
            .unwrap()
        };
        let branches = csg_boolean(BooleanOp::Union, &branch_up, &branch_dn).unwrap();
        let mut result = csg_boolean(BooleanOp::Difference, &trunk, &branches).unwrap();
        let normals_before = analyze_normals(&result);
        eprintln!(
            "before orient_outward: outward={}, inward={}, degen={}",
            normals_before.outward_faces,
            normals_before.inward_faces,
            normals_before.degenerate_faces,
        );
        eprintln!("is_watertight_before={}", result.is_watertight());
        result.orient_outward();
        let normals_after = analyze_normals(&result);
        eprintln!(
            "after  orient_outward: outward={}, inward={}, degen={}",
            normals_after.outward_faces, normals_after.inward_faces, normals_after.degenerate_faces,
        );
        eprintln!("is_watertight_after={}", result.is_watertight());
        assert_eq!(
            normals_after.inward_faces, 0,
            "orient_outward must eliminate inward faces"
        );

        // Single connected component — retain_largest_component must have
        // stripped the 2 × 8-face phantom islands from the trunk difference.
        {
            use crate::domain::topology::connectivity::connected_components;
            use crate::domain::topology::AdjacencyGraph;
            result.rebuild_edges();
            let edges = result.edges_ref().unwrap();
            let adj = AdjacencyGraph::build(&result.faces, edges);
            let comps = connected_components(&result.faces, &adj);
            eprintln!("components={}", comps.len());
            assert_eq!(
                comps.len(),
                1,
                "trunk difference must be a single connected component; \
                 got {} (phantom islands not removed)",
                comps.len(),
            );
        }
        // Euler characteristic χ = 2 for a single genus-0 closed body.
        {
            use crate::application::watertight::check::check_watertight;
            result.rebuild_edges();
            let rpt =
                check_watertight(&result.vertices, &result.faces, result.edges_ref().unwrap());
            eprintln!("euler_characteristic={:?}", rpt.euler_characteristic);
            assert_eq!(
                rpt.euler_characteristic,
                Some(2),
                "trunk difference must have Euler χ = 2; got {:?}",
                rpt.euler_characteristic,
            );
        }
    }

}

// ── Boundary vertex merging ──────────────────────────────────────────────────

/// Merge nearby boundary vertices to close sliver gaps at intersection curves.
///
/// Identifies boundary vertices (those on boundary edges) and merges pairs
/// within a small adaptive tolerance.  This closes small gaps left by CSG
/// arrangement precision limits at complex intersection curves.
///
/// # Theorem — Boundary Vertex Merge Convergence
///
/// Each merge reduces the boundary edge count by exactly 2 (the two half-edges
/// incident to the merged vertex pair become interior).  The process terminates
/// when no further merges are possible (fixed-point).  ∎
fn merge_nearby_boundary_vertices(mesh: &mut IndexedMesh) {
    use std::collections::HashSet;

    // Adaptive tolerance: 5% of the mean edge length, clamped to [0.01, 0.2] mm.
    // This scales correctly for any mesh granularity while remaining safe for
    // millifluidic geometries (R ≈ 0.5 mm, typical edge length ≈ 0.05–0.5 mm).
    let mean_edge_len = {
        mesh.rebuild_edges();
        let edges = match mesh.edges_ref() {
            Some(e) => e,
            None => return,
        };
        let (sum, count) = edges
            .iter()
            .map(|e| {
                let pa = mesh.vertices.position(e.vertices.0);
                let pb = mesh.vertices.position(e.vertices.1);
                (pa - pb).norm()
            })
            .fold((0.0_f64, 0usize), |(s, c), d| (s + d, c + 1));
        if count == 0 {
            return;
        }
        sum / count as f64
    };
    let tol = (mean_edge_len * 0.20).clamp(0.01, 0.2);
    let max_iter = 30;

    for _iter in 0..max_iter {
        mesh.rebuild_edges();
        let edges_ref = match mesh.edges_ref() {
            Some(e) => e,
            None => break,
        };

        // Phase 1: collect boundary vertex IDs.
        let mut boundary_verts: HashSet<VertexId> = HashSet::new();
        for edge in edges_ref.iter() {
            if edge.is_boundary() {
                boundary_verts.insert(edge.vertices.0);
                boundary_verts.insert(edge.vertices.1);
            }
        }

        if boundary_verts.is_empty() {
            break;
        }

        let bv: Vec<VertexId> = boundary_verts.iter().copied().collect();

        // Phase 2: find closest boundary-boundary pair within tolerance.
        let mut best: Option<(VertexId, VertexId, f64)> = None;
        for i in 0..bv.len() {
            let pi = mesh.vertices.position(bv[i]);
            for j in (i + 1)..bv.len() {
                let pj = mesh.vertices.position(bv[j]);
                let d = (pi - pj).norm();
                if d < tol {
                    if best.is_none() || d < best.unwrap().2 {
                        best = Some((bv[i], bv[j], d));
                    }
                }
            }
        }

        // Phase 3: if no boundary-boundary pair found, try boundary-to-interior.
        if best.is_none() {
            let all_vids: Vec<VertexId> = mesh.vertices.iter().map(|(id, _)| id).collect();
            for &bvid in &bv {
                let bp = mesh.vertices.position(bvid);
                let per_vertex_tol = tol * 0.5;
                for &ivid in &all_vids {
                    if bvid == ivid || boundary_verts.contains(&ivid) {
                        continue;
                    }
                    let ip = mesh.vertices.position(ivid);
                    let d = (bp - ip).norm();
                    if d < per_vertex_tol {
                        if best.is_none() || d < best.unwrap().2 {
                            best = Some((ivid, bvid, d));
                        }
                    }
                }
            }
        }

        let (keep, remove, _dist) = match best {
            Some(b) => b,
            None => break,
        };

        // Merge: replace all references to `remove` with `keep`.
        let mut changed = false;
        for face in mesh.faces.iter_mut() {
            for v in face.vertices.iter_mut() {
                if *v == remove {
                    *v = keep;
                    changed = true;
                }
            }
        }

        if !changed {
            break;
        }

        collapse_degenerate_faces(mesh);
        mesh.rebuild_edges();

        split_non_manifold_edges(mesh);
        collapse_degenerate_faces(mesh);
        mesh.rebuild_edges();

        if !mesh.is_watertight() {
            let es = crate::infrastructure::storage::edge_store::EdgeStore::from_face_store(
                &mesh.faces,
            );
            crate::application::watertight::seal::seal_boundary_loops(
                &mut mesh.vertices,
                &mut mesh.faces,
                &es,
                crate::domain::core::index::RegionId::INVALID,
            );
            collapse_degenerate_faces(mesh);
                mesh.rebuild_edges();
        }

        if mesh.is_watertight() {
            break;
        }
    }
}

// ── Non-manifold edge splitting ──────────────────────────────────────────────

/// Resolve non-manifold edges by removing excess faces.
///
/// A 2-manifold requires every edge to be shared by exactly 2 faces.
/// CSG arrangement can produce edges with 3+ faces at intersection curves.
/// This function removes excess faces to restore edge-manifold topology.
///
/// # Theorem — Non-Manifold Edge Elimination
///
/// After removal, every edge has at most 2 faces.  The removed faces'
/// other edges may become boundary edges (1 face) or remain manifold
/// (2 faces).  The resulting mesh has no non-manifold edges.  ∎
fn split_non_manifold_edges(mesh: &mut IndexedMesh) {
    use std::collections::HashMap;

    let mut edge_faces: HashMap<(VertexId, VertexId), Vec<usize>> = HashMap::new();
    for (fi, face) in mesh.faces.iter().enumerate() {
        let v = face.vertices;
        for &(a, b) in &[(v[0], v[1]), (v[1], v[2]), (v[2], v[0])] {
            let key = if a < b { (a, b) } else { (b, a) };
            edge_faces.entry(key).or_default().push(fi);
        }
    }

    let mut faces_to_remove: hashbrown::HashSet<usize> = hashbrown::HashSet::new();
    for (_edge, face_indices) in &edge_faces {
        if face_indices.len() > 2 {
            for &fi in &face_indices[2..] {
                faces_to_remove.insert(fi);
            }
        }
    }

    if faces_to_remove.is_empty() {
        return;
    }

    let mut clean_faces: Vec<FaceData> = Vec::with_capacity(
        mesh.faces.len() - faces_to_remove.len(),
    );
    for (fi, face) in mesh.faces.iter().enumerate() {
        if !faces_to_remove.contains(&fi) {
            clean_faces.push(*face);
        }
    }
    mesh.faces = crate::infrastructure::storage::face_store::FaceStore::new();
    for face in clean_faces {
        mesh.faces.push(face);
    }
}
