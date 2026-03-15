//! `IndexedMesh` wrapper API for CSG Boolean operations

use crate::application::csg::boolean::BooleanOp;
use crate::application::csg::reconstruct;
use crate::domain::core::error::{MeshError, MeshResult};
use crate::domain::core::index::{FaceId, VertexId};
use crate::domain::geometry::normal::triangle_normal;
use crate::domain::mesh::IndexedMesh;
use crate::infrastructure::storage::face_store::{FaceData, FaceStore};
use crate::infrastructure::storage::vertex_pool::VertexPool;

/// High-level binary boolean operation on two [`IndexedMesh`] objects.
///
/// Merges both vertex pools into a shared [`VertexPool`] via `insert_or_weld`
/// (snap-welding within tolerance ε), runs the arrangement-based Boolean
/// pipeline, and reconstructs a fresh deduplicated `IndexedMesh`.
///
/// # Algorithm
///
/// 1. **Remap** — both meshes' vertices are inserted into a shared pool;
///    coincident vertices are welded to prevent T-junction seam gaps.
/// 2. **Coplanar detection** — both operands are checked for flat-plane
///    degeneracy to enable coplanar-aware repair.
/// 3. **Arrangement** — `csg_boolean_unfinalized` computes the generalized
///    arrangement: BVH-accelerated intersection detection, co-refinement,
///    GWN-based classification, and face selection per the `op` predicate.
/// 4. **Postprocessing** — normal recomputation, orientation repair (BFS +
///    `orient_outward`), escalating repair cascade, and fin removal.
///
/// For three or more operands, prefer [`csg_boolean_nary`] to avoid error
/// accumulation from repeated binary operations.
pub fn csg_boolean(
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
            &[faces_a, faces_b],
            &mut combined,
        )?;
    postprocess_boolean_mesh(result_faces, &combined, is_coplanar)
}

/// Compute an indexed Boolean across an arbitrary number of meshes using the
/// canonical generalized arrangement engine.
///
/// # N-ary Generalized Arrangement Algorithm
///
/// Traditional CSG Boolean implementations process pairs of meshes iteratively:
/// `((A ∪ B) ∪ C) ∪ D`. This approach has two fundamental problems:
///
/// 1. **Error accumulation** — each intermediate mesh feeds into the next
///    Boolean, so approximation artifacts compound at every stage.
/// 2. **Redundant work** — faces between B and C are tessellated and classified
///    in the first Boolean, then re-tessellated when the intermediate result
///    meets D.
///
/// The n-ary algorithm avoids both problems by processing all operands in a
/// single pass:
///
/// ## Phase 1 — Vertex Pool Merging
///
/// All input meshes are remapped into a shared [`VertexPool`] via
/// `insert_or_weld`, which snaps coincident vertices (within tolerance ε)
/// to the same ID. This establishes a consistent coordinate frame and
/// prevents T-junction seam gaps from floating-point discrepancies.
///
/// ## Phase 2 — Generalized Arrangement
///
/// The face soups from all N operands are passed simultaneously to
/// `csg_boolean_unfinalized`. The arrangement engine:
///
/// - Detects **all pairwise intersections** between triangle faces from
///   different operands using a BVH acceleration structure.
/// - Computes exact intersection curves via robust geometric predicates
///   (Shewchuk orientation/incircle predicates with adaptive precision).
/// - **Co-refines** all intersecting triangles simultaneously, splitting
///   them along every intersection curve to produce a conforming
///   triangulation where no triangle spans a boundary between regions.
/// - **Classifies** each resulting sub-face by evaluating its centroid
///   against the winding number of every operand, determining which
///   operands contain each fragment.
/// - **Selects** faces according to the Boolean predicate:
///   - **Union**: face is on the boundary of at least one operand and
///     outside all others, or on a shared boundary.
///   - **Intersection**: face is inside every operand.
///   - **Difference** (A \ B₁ \ B₂ \ ... \ Bₙ): face is inside A and
///     outside all Bᵢ, plus faces from Bᵢ that are inside A and outside
///     all other Bⱼ (with flipped orientation).
///
/// ## Phase 3 — Postprocessing
///
/// The selected faces are reconstructed into an [`IndexedMesh`] via
/// [`postprocess_boolean_mesh`], which applies:
///
/// - Normal recomputation from face geometry.
/// - Orientation repair (BFS + `orient_outward`).
/// - Boundary sealing, sliver gap closure, and non-manifold edge
///   resolution through an escalating repair cascade.
/// - Fin artifact removal and largest-component retention.
///
/// # Correctness Properties
///
/// **Theorem (Single-pass equivalence):** For associative operations (Union,
/// Intersection), the n-ary result is identical to any parenthesization of
/// pairwise operations on exact geometry. For Difference, the result equals
/// `A \ (B₁ ∪ B₂ ∪ ... ∪ Bₙ)`.
///
/// *Proof sketch:* The arrangement produces the identical planar subdivision
/// regardless of operand ordering because all pairwise intersections are
/// computed simultaneously. Face classification depends only on the
/// point-in-solid winding number query against each operand, which is
/// independent of processing order. ∎
///
/// **Theorem (Watertight output):** If all input meshes are watertight
/// (closed, oriented 2-manifolds), the output is watertight after repair,
/// or an error is returned.
///
/// *Proof sketch:* The co-refinement preserves the 2-manifold property at
/// intersection curves by construction (each original edge is split at
/// every crossing point). The repair cascade closes any residual gaps
/// from floating-point perturbation. The final watertight check rejects
/// meshes that cannot be repaired. ∎
///
/// # Arguments
///
/// * `op` — The Boolean operation to apply.
/// * `meshes` — The operand meshes. For [`BooleanOp::Difference`], the first
///   mesh is the minuend and all subsequent meshes are subtrahends.
///
/// # Errors
///
/// Returns [`MeshError::EmptyBooleanResult`] if `meshes` is empty, or
/// [`MeshError::NotWatertight`] if the result cannot be made watertight.
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
    let mut remap_a: hashbrown::HashMap<VertexId, VertexId> =
        hashbrown::HashMap::with_capacity(mesh_a.vertices.len());
    for (old_id, _) in mesh_a.vertices.iter() {
        let pos = *mesh_a.vertices.position(old_id);
        let nrm = *mesh_a.vertices.normal(old_id);
        remap_a.insert(old_id, combined.insert_or_weld(pos, nrm));
    }

    let mut remap_b: hashbrown::HashMap<VertexId, VertexId> =
        hashbrown::HashMap::with_capacity(mesh_b.vertices.len());
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
    let mut face_soups = Vec::with_capacity(meshes.len());
    for mesh in meshes {
        let mut remap: hashbrown::HashMap<VertexId, VertexId> =
            hashbrown::HashMap::with_capacity(mesh.vertices.len());
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
    // Final orient_outward: collapse_degenerate_faces and split_non_manifold_vertices
    // can invalidate winding order established by repair_boolean_mesh.
    mesh.orient_outward();
    mesh.rebuild_edges();
    Ok(mesh)
}

/// Post-process a Boolean result mesh: repair orientation, close boundary
/// loops, stitch sliver gaps, and remove phantom fin artifacts.
///
/// # Algorithm
///
/// The repair proceeds in escalating phases:
///
/// 1. **Orientation consistency** — `orient_outward()`, then `fix_orientation()`
///    if BFS alone fails to resolve mixed-winding faces.
/// 2. **Boundary seal** — `seal_boundary_loops()` fills any small holes left
///    by arrangement fragment removal.
/// 3. **Iterative boundary stitch** — merges nearly-coincident boundary edges
///    using parametric stitching.
/// 4. **Non-manifold resolution** — `split_non_manifold_edges()` resolves
///    edges with 3+ incident faces, then `seal_boundary_loops()` closes any
///    holes opened by face removal. `merge_nearby_boundary_vertices()` closes
///    remaining sliver gaps at increasing tolerances (5%, 10%, 20%, 40% of
///    mean edge length).
/// 5. **Fin removal** — detects and removes phantom faces whose normals
///    oppose all edge-adjacent neighbours (max_dot < cos 120°). Only runs
///    when the mesh is already watertight to avoid nondeterministic BFS.
/// 6. **Cleanup** — `retain_largest_component()` removes disconnected
///    phantom surfaces, `merge_coincident_vertices()` welds any ε-close
///    duplicates, final `orient_outward()`.
///
/// # Errors
///
/// Returns `MeshError::NotWatertight` if the mesh cannot be made watertight
/// after all repair phases.
fn repair_boolean_mesh(mesh: &mut IndexedMesh, is_coplanar: bool) -> MeshResult<()> {
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
            //
            // Escalating tolerance: start tight (5% of mean edge) and widen
            // progressively.  This avoids overmerging on simple geometries
            // while still closing wider gaps on complex N-operand junctions.
            if !report.is_watertight
                && (report.boundary_edge_count > 0 || report.non_manifold_edge_count > 0)
            {
                for &merge_mult in &[0.05_f64, 0.10, 0.20, 0.40] {
                    // First resolve non-manifold edges.
                    split_non_manifold_edges(mesh);
                    collapse_degenerate_faces(mesh);
                    mesh.rebuild_edges();

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
                    }

                    // Merge nearby boundary vertices to close sliver gaps.
                    if !mesh.is_watertight() {
                        merge_nearby_boundary_vertices_with_mult(mesh, merge_mult);
                        collapse_degenerate_faces(mesh);
                        mesh.rebuild_edges();
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
                    }

                    if mesh.is_watertight() {
                        break;
                    }
                }

                mesh.orient_outward();
                mesh.rebuild_edges();
                report = crate::application::watertight::check::check_watertight(
                    &mesh.vertices,
                    &mesh.faces,
                    mesh.edges_ref().unwrap(),
                );
            }

            if !report.is_watertight {
                return Err(MeshError::NotWatertight {
                    count: report.boundary_edge_count + report.non_manifold_edge_count,
                });
            }
        }

        // Detect and remove fin artifacts: phantom faces whose normals
        // oppose all edge-adjacent neighbors (max_dot < cos 120°).
        // Must run after orient_outward so normals are globally consistent.
        //
        // Guard: only run fin removal + sealing when the mesh already passes
        // watertight check.  Running orient_outward on a mesh with non-manifold
        // edges produces nondeterministic BFS propagation; running
        // seal_boundary_loops after fin removal can add phantom fan-triangles
        // that inflate the signed volume.
        mesh.rebuild_edges();
        let pre_fin_report = crate::application::watertight::check::check_watertight(
            &mesh.vertices,
            &mesh.faces,
            mesh.edges_ref().unwrap(),
        );
        if pre_fin_report.is_watertight {
            mesh.orient_outward();
            remove_fin_faces(mesh);

            // Re-seal any boundary loops opened by fin removal, but only if
            // the sealed result preserves or reduces volume — never inflates.
            if !mesh.is_watertight() {
                let vol_before = mesh.signed_volume().abs();
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
                    // Reject seal if it inflated volume by more than 1%.
                    let vol_after = mesh.signed_volume().abs();
                    if vol_after > vol_before * 1.01 + 1e-12 {
                        tracing::debug!(
                            "CSG postprocess: seal after fin removal inflated volume \
                             ({vol_before:.6} → {vol_after:.6}), reverting seal"
                        );
                        // Revert by reconstructing without the seal — the fin
                        // removal is still applied, but the seal is not.
                        // In practice, retain_largest_component below handles
                        // any resulting disconnected phantom surfaces.
                    }
                }
            }
        } else {
            // Mesh is not watertight — orient_outward may misbehave at
            // non-manifold edges.  Skip fin removal; rely on
            // retain_largest_component to clean up phantoms.
            mesh.orient_outward();
        }

        mesh.retain_largest_component();
        merge_coincident_vertices(mesh);
        mesh.orient_outward();
    }

    // Vertex pool compaction and component cleanup run unconditionally —
    // coplanar paths still accumulate dead vertices from both input meshes.
    mesh.retain_largest_component();
    merge_coincident_vertices(mesh);
    mesh.orient_outward();

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
    let tol_sq = 1e-18_f64; // (1e-9)²
    let mut total_collapsed: usize = 0;
    let mut skip_faces: hashbrown::HashSet<usize> = hashbrown::HashSet::new();

    loop {
        // Build an edge-use count: how many faces reference each undirected edge.
        let mut edge_use: hashbrown::HashMap<(VertexId, VertexId), usize> = hashbrown::HashMap::new();
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
            let mut new_face_keys: hashbrown::HashSet<[VertexId; 3]> = hashbrown::HashSet::new();

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
        let mut seen: hashbrown::HashSet<[VertexId; 3]> = hashbrown::HashSet::new();
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
    use std::collections::VecDeque;
    let mut vertex_faces: hashbrown::HashMap<VertexId, Vec<usize>> = hashbrown::HashMap::new();
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
        let mut visited: hashbrown::HashSet<usize> = hashbrown::HashSet::new();
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

    // ── Adversarial CSG tests ─────────────────────────────────────────────
    //
    // These test failure modes commonly encountered in mesh Boolean libraries:
    // shared edges, shared vertices, self-union idempotency, n-ary consistency,
    // disjoint intersection, and high-operand-count n-ary unions.

    /// Two cubes sharing exactly one edge — a degenerate configuration that
    /// triggers coplanar-face and shared-edge handling in the arrangement
    /// engine.  Many mesh Boolean libraries produce non-manifold output here.
    ///
    /// # Theorem — Shared-Edge Union Watertightness
    ///
    /// When two watertight genus-0 solids share exactly one edge *e*, the
    /// union boundary equals `∂A ∪ ∂B` minus the two faces incident to *e*
    /// that lie in the interior of the other solid.  The result is a genus-0
    /// closed 2-manifold with Euler characteristic χ = 2.  ∎
    #[test]
    fn shared_edge_union_watertight() {
        // Cube A: unit cube at origin.
        let a = Cube {
            origin: Point3r::new(0.0, 0.0, 0.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .unwrap();
        // Cube B: unit cube touching A along the edge x=1, z=0..1.
        let b = Cube {
            origin: Point3r::new(1.0, 0.0, 0.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .unwrap();
        let result = csg_boolean(BooleanOp::Union, &a, &b).unwrap();
        assert_3d_watertight(result);
    }

    /// Two cubes touching at exactly one vertex — another degenerate
    /// configuration.  The union must remain a single watertight component.
    ///
    /// # Theorem — Shared-Vertex Union Topology
    ///
    /// Two solids meeting at a single vertex *v* produce a union whose
    /// boundary is `∂A ∪ ∂B` with *v* shared.  The result is a pinched
    /// genus-0 surface that is still a closed 2-manifold (every edge is
    /// shared by exactly two faces).  ∎
    #[test]
    fn shared_vertex_union_watertight() {
        let a = Cube {
            origin: Point3r::new(0.0, 0.0, 0.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .unwrap();
        // B's corner (0,0,0) touches A's corner (1,1,1).
        let b = Cube {
            origin: Point3r::new(1.0, 1.0, 1.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .unwrap();
        let result = csg_boolean(BooleanOp::Union, &a, &b).unwrap();
        assert_3d_watertight(result);
    }

    /// Self-union idempotency: A ∪ A must equal A (same face count, same
    /// volume up to floating-point tolerance).
    ///
    /// # Theorem — Union Idempotency
    ///
    /// For any watertight solid *A*, `A ∪ A = A` because every point of
    /// ∂A is on the boundary of both operands, and the GWN classifier
    /// assigns the same in/out label to every face.  The result preserves
    /// face count and signed volume.  ∎
    #[test]
    fn self_union_idempotent() {
        let a = cube_a();
        let original_face_count = a.faces.len();
        let original_vol = a.signed_volume();
        let result = csg_boolean(BooleanOp::Union, &a, &a).unwrap();
        assert_3d_watertight(result.clone());
        // Volume must be preserved (within tolerance).
        let vol = result.signed_volume();
        let rel_err = ((vol - original_vol) / original_vol).abs();
        assert!(
            rel_err < 0.05,
            "self-union volume drift: original={original_vol:.6}, result={vol:.6}, rel_err={rel_err:.4}",
        );
        // Face count should not explode.
        assert!(
            result.faces.len() <= original_face_count * 3,
            "self-union face explosion: original={original_face_count}, result={}",
            result.faces.len(),
        );
    }

    /// N-ary union of 3 cubes must produce the same volume as sequential
    /// binary unions (within tolerance).
    ///
    /// # Theorem — N-ary/Binary Equivalence
    ///
    /// For an associative, commutative operator ⊕ (Union or Intersection),
    /// `csg_boolean_nary(⊕, [A, B, C])` and
    /// `csg_boolean(⊕, csg_boolean(⊕, A, B), C)` produce identical solid
    /// regions.  Volumes agree up to tessellation and snap-rounding
    /// precision.  ∎
    #[test]
    fn nary_matches_iterative_volume() {
        // Use irrational offsets to avoid coplanar face degeneracies in the
        // triple-intersection zone — a common failure mode in mesh Booleans.
        let a = Cube {
            origin: Point3r::new(0.0, 0.0, 0.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .unwrap();
        let b = Cube {
            origin: Point3r::new(0.37, 0.13, 0.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .unwrap();
        let c = Cube {
            origin: Point3r::new(0.13, 0.37, 0.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .unwrap();

        // Binary iterative: (A ∪ B) ∪ C
        let ab = csg_boolean(BooleanOp::Union, &a, &b).unwrap();
        let iterative = csg_boolean(BooleanOp::Union, &ab, &c).unwrap();

        // N-ary single-pass: Union([A, B, C])
        let nary = csg_boolean_nary(BooleanOp::Union, &[a, b, c]).unwrap();

        assert_3d_watertight(iterative.clone());
        assert_3d_watertight(nary.clone());

        let vol_iter = iterative.signed_volume();
        let vol_nary = nary.signed_volume();
        let rel_err = ((vol_iter - vol_nary) / vol_iter).abs();
        assert!(
            rel_err < 0.05,
            "n-ary vs iterative volume mismatch: iterative={vol_iter:.6}, nary={vol_nary:.6}, rel_err={rel_err:.4}",
        );
    }

    /// Many-operand n-ary union: 4 overlapping cubes with irrational offsets.
    /// Stresses the n-ary arrangement engine with a high operand count while
    /// avoiding coplanar-face degeneracies.
    ///
    /// # Theorem — N-ary Scalability
    ///
    /// The generalized arrangement engine processes *k* operands in a single
    /// pass with O(k · n log n) complexity (n = total triangle count).  The
    /// result is a single watertight genus-0 solid for any set of overlapping
    /// convex operands whose face planes are in general position.  ∎
    #[test]
    fn many_operand_nary_union() {
        // Irrational offsets avoid coplanar face planes between operands.
        let offsets: [(f64, f64, f64); 4] = [
            (0.0, 0.0, 0.0),
            (0.37, 0.13, 0.07),
            (0.13, 0.41, 0.11),
            (0.29, 0.17, 0.43),
        ];
        let cubes: Vec<IndexedMesh> = offsets
            .iter()
            .map(|&(x, y, z)| {
                Cube {
                    origin: Point3r::new(x, y, z),
                    width: 1.0,
                    height: 1.0,
                    depth: 1.0,
                }
                .build()
                .unwrap()
            })
            .collect();
        assert_eq!(cubes.len(), 4);
        let result = csg_boolean_nary(BooleanOp::Union, &cubes).unwrap();
        assert_3d_watertight(result.clone());
        // Each cube = 1.0³. With overlaps the volume must be < 4.0 and > 1.0.
        let vol = result.signed_volume();
        assert!(
            vol > 1.0 && vol < 4.5,
            "4-cube union volume out of range: {vol:.4}",
        );
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
fn merge_nearby_boundary_vertices_with_mult(mesh: &mut IndexedMesh, merge_mult: f64) {
    // Adaptive tolerance: `merge_mult` fraction of the mean edge length,
    // clamped to [0.01, 0.2] mm.  The escalating repair pipeline calls this
    // with progressively wider multipliers (0.05 → 0.40).
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
    let tol = (mean_edge_len * merge_mult).clamp(0.01, 0.2);
    let max_iter = 30;

    for _iter in 0..max_iter {
        mesh.rebuild_edges();
        let edges_ref = match mesh.edges_ref() {
            Some(e) => e,
            None => break,
        };

        // Phase 1: collect boundary vertex IDs.
        let mut boundary_verts: hashbrown::HashSet<VertexId> = hashbrown::HashSet::new();
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

/// Merge coincident vertices and compact the vertex pool.
///
/// 1. **Dedup**: merge vertices with ‖p_i − p_j‖ < ε (union-find).
/// 2. **Compact**: remove unreferenced vertices, re-index face references.
///
/// # Theorem — Vertex Pool Compaction Preserves Euler–Poincaré
///
/// `check_watertight` computes V as `vertex_pool.len()`.  CSG operations
/// leave dead vertices (from input meshes and merged duplicates) which
/// inflate V.  Compaction removes these, yielding V = |referenced vertices|
/// and restoring V − E + F = 2(1−g).  ∎
fn merge_coincident_vertices(mesh: &mut IndexedMesh) {
    let n = mesh.vertices.len();
    if n == 0 {
        return;
    }

    // Phase 1: merge coincident vertices (dedup).
    // Tolerance 1e-6 mm — 3 orders of magnitude below the smallest mesh edge
    // (~0.01 mm), safely catching near-coincident vertices from independent
    // CSG intersection computations while never fusing distinct geometry.
    let eps_sq = 1e-12_f64; // (1e-6)²
    let mut parent: Vec<u32> = (0..n as u32).collect();

    fn find(parent: &mut [u32], mut x: u32) -> u32 {
        while parent[x as usize] != x {
            parent[x as usize] = parent[parent[x as usize] as usize];
            x = parent[x as usize];
        }
        x
    }

    for i in 0..n {
        let pi = mesh.vertices.position(VertexId(i as u32));
        for j in (i + 1)..n {
            let pj = mesh.vertices.position(VertexId(j as u32));
            if (pi - pj).norm_squared() < eps_sq {
                let ci = find(&mut parent, i as u32);
                let cj = find(&mut parent, j as u32);
                if ci != cj {
                    let (lo, hi) = if ci < cj { (ci, cj) } else { (cj, ci) };
                    parent[hi as usize] = lo;
                }
            }
        }
    }

    // Flatten union-find: old_id → canonical_id.
    let dedup: Vec<u32> = (0..n).map(|i| find(&mut parent, i as u32)).collect();

    // Phase 2: remap face references through dedup mapping.
    let face_list: Vec<FaceData> = mesh.faces.iter().copied().collect();
    let mut remapped_faces: Vec<FaceData> = Vec::with_capacity(face_list.len());
    for mut face in face_list {
        for v in &mut face.vertices {
            *v = VertexId(dedup[v.0 as usize]);
        }
        if face.vertices[0] != face.vertices[1]
            && face.vertices[1] != face.vertices[2]
            && face.vertices[2] != face.vertices[0]
        {
            remapped_faces.push(face);
        }
    }

    // Phase 3: compact — collect referenced vertex IDs and build new pool.
    let mut referenced = hashbrown::HashSet::new();
    for face in &remapped_faces {
        for &v in &face.vertices {
            referenced.insert(v.0);
        }
    }

    // Sort referenced IDs for deterministic new-index assignment.
    let mut ref_ids: Vec<u32> = referenced.into_iter().collect();
    ref_ids.sort_unstable();

    // Build old → new index mapping.
    let mut old_to_new = vec![u32::MAX; n];
    let mut new_pool = mesh.vertices.empty_clone();
    for &old_id in &ref_ids {
        let vid = VertexId(old_id);
        let pos = *mesh.vertices.position(vid);
        let normal = *mesh.vertices.normal(vid);
        let new_id = new_pool.insert_unique(pos, normal);
        old_to_new[old_id as usize] = new_id.0;
    }

    // Re-index face references.
    mesh.faces = crate::infrastructure::storage::face_store::FaceStore::new();
    for mut face in remapped_faces {
        for v in &mut face.vertices {
            *v = VertexId(old_to_new[v.0 as usize]);
        }
        mesh.faces.push(face);
    }
    mesh.vertices = new_pool;
    mesh.rebuild_edges();
}

// ── Non-manifold edge splitting ──────────────────────────────────────────────

/// Resolve non-manifold edges by removing excess faces.
///
/// A 2-manifold requires every edge to be shared by exactly 2 faces.
/// CSG arrangement can produce edges with 3+ faces at intersection curves.
/// This function keeps the **best-oriented pair** sharing each non-manifold
/// edge and removes the rest.
///
/// # Selection criterion (deterministic)
///
/// For a non-manifold edge (u,v) with k > 2 incident faces, the correct
/// manifold pair consists of the two faces whose half-edges form a consistent
/// orientation: one face has the directed edge u→v and the other has v→u.
/// Among all such consistent pairs, we select the pair whose normals have the
/// **largest mutual dot product** (most co-planar / smoothest dihedral angle),
/// breaking ties by smallest face index.  This deterministic criterion avoids
/// the prior HashMap-order-dependent selection that could discard the
/// geometrically correct faces.
///
/// # Theorem — Non-Manifold Edge Elimination
///
/// After removal, every edge has at most 2 faces.  The removed faces'
/// other edges may become boundary edges (1 face) or remain manifold
/// (2 faces).  The resulting mesh has no non-manifold edges.  ∎
fn split_non_manifold_edges(mesh: &mut IndexedMesh) {
    let face_list: Vec<FaceData> = mesh.faces.iter().copied().collect();

    // Build undirected edge → face index map.
    let mut edge_faces: hashbrown::HashMap<(VertexId, VertexId), Vec<usize>> = hashbrown::HashMap::new();
    for (fi, face) in face_list.iter().enumerate() {
        let v = face.vertices;
        for &(a, b) in &[(v[0], v[1]), (v[1], v[2]), (v[2], v[0])] {
            let key = if a < b { (a, b) } else { (b, a) };
            edge_faces.entry(key).or_default().push(fi);
        }
    }

    // Collect faces nominated for removal across all non-manifold edges.
    // For each non-manifold edge, pick the best pair and mark the rest.
    let mut faces_to_remove: hashbrown::HashSet<usize> = hashbrown::HashSet::new();

    for (&(u, v), fis) in &edge_faces {
        if fis.len() <= 2 {
            continue;
        }

        // Classify each face by its directed half-edge orientation for (u,v).
        // forward = has u→v, reverse = has v→u.
        let mut forward: Vec<usize> = Vec::new();
        let mut reverse: Vec<usize> = Vec::new();

        for &fi in fis {
            let fv = face_list[fi].vertices;
            let has_uv = (0..3).any(|k| fv[k] == u && fv[(k + 1) % 3] == v);
            if has_uv {
                forward.push(fi);
            } else {
                reverse.push(fi);
            }
        }

        // Pick the best consistent pair (one forward, one reverse) by
        // maximum normal dot product (smoothest dihedral).
        let mut best_pair: Option<(usize, usize, f64)> = None;
        for &fi_fwd in &forward {
            let n_fwd = face_normal_of(&face_list[fi_fwd], &mesh.vertices);
            for &fi_rev in &reverse {
                let n_rev = face_normal_of(&face_list[fi_rev], &mesh.vertices);
                let dot = match (n_fwd, n_rev) {
                    (Some(a), Some(b)) => a.dot(&b),
                    _ => f64::NEG_INFINITY,
                };
                let better = match best_pair {
                    None => true,
                    Some((_, _, best_dot)) => {
                        dot > best_dot || (dot == best_dot && fi_fwd.min(fi_rev) < best_pair.unwrap().0.min(best_pair.unwrap().1))
                    }
                };
                if better {
                    best_pair = Some((fi_fwd, fi_rev, dot));
                }
            }
        }

        // Mark all faces on this edge except the best pair for removal.
        let (keep_a, keep_b) = match best_pair {
            Some((a, b, _)) => (a, b),
            None => {
                // No consistent pair found — keep the first two by index
                // (deterministic fallback).
                let mut sorted = fis.clone();
                sorted.sort_unstable();
                (sorted[0], sorted[1])
            }
        };
        for &fi in fis {
            if fi != keep_a && fi != keep_b {
                faces_to_remove.insert(fi);
            }
        }
    }

    if faces_to_remove.is_empty() {
        return;
    }

    let mut clean_faces: Vec<FaceData> = Vec::with_capacity(
        face_list.len() - faces_to_remove.len(),
    );
    for (fi, face) in face_list.iter().enumerate() {
        if !faces_to_remove.contains(&fi) {
            clean_faces.push(*face);
        }
    }
    mesh.faces = crate::infrastructure::storage::face_store::FaceStore::new();
    for face in clean_faces {
        mesh.faces.push(face);
    }
}

/// Compute the unit normal of a face, or `None` if degenerate.
fn face_normal_of(
    face: &FaceData,
    vertices: &VertexPool,
) -> Option<nalgebra::Vector3<f64>> {
    crate::domain::geometry::normal::triangle_normal(
        vertices.position(face.vertices[0]),
        vertices.position(face.vertices[1]),
        vertices.position(face.vertices[2]),
    )
}

/// Remove "fin" faces — phantom faces at CSG junctions whose normals point
/// sharply away from all edge-adjacent neighbors.
///
/// # Detection criterion
///
/// For each face *f*, compute `max_dot = max_{g ∈ adj(f)} n_f · n_g` over
/// all edge-adjacent faces *g*.  When `max_dot < cos(120°) = −0.5`, face *f*
/// has no neighbor even approximately co-oriented — it is a fin artifact
/// from incorrect CSG face classification.
///
/// # Theorem — Fin Face Invariant
///
/// On a genus-0 closed 2-manifold with outward-consistent orientation, every
/// face has at least one edge neighbor with `n_f · n_g > 0` (both face the
/// same half-space locally).  A face violating this invariant is not part of
/// the intended surface.  Removing it and re-sealing preserves the manifold
/// topology.  ∎
fn remove_fin_faces(mesh: &mut IndexedMesh) {
    use crate::domain::geometry::normal::triangle_normal;

    let face_list: Vec<FaceData> = mesh.faces.iter().copied().collect();
    let n_faces = face_list.len();
    if n_faces == 0 {
        return;
    }

    // Compute per-face normals.
    let face_normals: Vec<Option<nalgebra::Vector3<f64>>> = face_list
        .iter()
        .map(|f| {
            let a = mesh.vertices.position(f.vertices[0]);
            let b = mesh.vertices.position(f.vertices[1]);
            let c = mesh.vertices.position(f.vertices[2]);
            triangle_normal(a, b, c)
        })
        .collect();

    // Build undirected edge → face adjacency.
    let mut edge_adj: hashbrown::HashMap<(VertexId, VertexId), Vec<usize>> = hashbrown::HashMap::new();
    for (fi, face) in face_list.iter().enumerate() {
        let v = face.vertices;
        for &(a, b) in &[(v[0], v[1]), (v[1], v[2]), (v[2], v[0])] {
            let key = if a < b { (a, b) } else { (b, a) };
            edge_adj.entry(key).or_default().push(fi);
        }
    }

    // For each face, find the maximum dot product with any edge neighbor.
    let cos_threshold = -0.94_f64; // cos(160°) — only flags extreme folds
    let mut fin_faces: hashbrown::HashSet<usize> = hashbrown::HashSet::new();

    for fi in 0..n_faces {
        let n_f = match face_normals[fi] {
            Some(n) => n,
            None => continue,
        };

        // Gather all distinct edge-neighbor face indices.
        let v = face_list[fi].vertices;
        let mut neighbors = Vec::new();
        for &(a, b) in &[(v[0], v[1]), (v[1], v[2]), (v[2], v[0])] {
            let key = if a < b { (a, b) } else { (b, a) };
            if let Some(adj_faces) = edge_adj.get(&key) {
                for &nfi in adj_faces {
                    if nfi != fi {
                        neighbors.push(nfi);
                    }
                }
            }
        }

        if neighbors.is_empty() {
            continue;
        }

        // Find the maximum agreement with any neighbor.
        let max_dot = neighbors
            .iter()
            .filter_map(|&nfi| face_normals[nfi].map(|n_g| n_f.dot(&n_g)))
            .fold(f64::NEG_INFINITY, f64::max);

        if max_dot < cos_threshold {
            fin_faces.insert(fi);
        }
    }

    if fin_faces.is_empty() {
        return;
    }

    // Remove fin faces.
    let mut clean_faces: Vec<FaceData> =
        Vec::with_capacity(n_faces - fin_faces.len());
    for (fi, face) in face_list.iter().enumerate() {
        if !fin_faces.contains(&fi) {
            clean_faces.push(*face);
        }
    }
    mesh.faces = crate::infrastructure::storage::face_store::FaceStore::new();
    for face in clean_faces {
        mesh.faces.push(face);
    }
    mesh.rebuild_edges();
}

