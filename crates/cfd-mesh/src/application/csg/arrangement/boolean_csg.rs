//! Canonical mesh-arrangement CSG pipeline for binary and multi-mesh operations.
use crate::application::csg::arrangement::coplanar_dispatch::dispatch_boolean_coplanar;
use crate::application::csg::arrangement::fragment_refinement::{
    append_corefined_fragments,
};
use crate::application::csg::arrangement::multi_mesh_resolution::resolve_multi_mesh_fragments;
use crate::application::csg::boolean::containment::{containment, Containment};
pub use crate::application::csg::arrangement::multi_mesh_resolution::BooleanFragmentRecord;
use crate::application::csg::arrangement::propagate::propagate_seam_vertices;
use crate::application::csg::arrangement::result_finalization::finalize_boolean_faces;
use crate::application::csg::broad_phase::triangle_aabb;
use crate::application::csg::intersect::{intersect_triangles, IntersectionType, SnapSegment};
use crate::domain::core::error::{MeshError, MeshResult};
use crate::domain::geometry::aabb::Aabb;
use crate::infrastructure::spatial::bvh::with_bvh;
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;

/// Broad-phase overlapping face pair between two meshes in the canonical Boolean pass.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BooleanCandidatePair {
    pub mesh_a: usize,
    pub face_a: usize,
    pub mesh_b: usize,
    pub face_b: usize,
}

/// Defines the semantic target for the generalized Boolean arrangement.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BooleanOp {
    /// Result survives if exterior to all other volumes. (A ∪ B ∪ ...)
    Union,
    /// Result survives if interior to all other volumes. (A ∩ B ∩ ...)
    Intersection,
    /// Base mesh (index 0) carving. Survives if base is exterior to all others, OR if subtractors (index > 0) are interior to base AND exterior to all other subtractors.
    Difference,
}

/// Computes an exact Boolean result for one or more meshes using a single generalized arrangement pass.
///
/// This avoids the progressive geometry degradation associated with iteratively computing
/// `(A op B) op C`, solving directly for the global multi-way topological boundaries.
pub fn csg_boolean(
    op: BooleanOp,
    meshes: &[Vec<FaceData>],
    pool: &mut VertexPool,
) -> MeshResult<Vec<FaceData>> {
    csg_boolean_with_policy(op, meshes, pool, true)
}

pub(crate) fn csg_boolean_unfinalized(
    op: BooleanOp,
    meshes: &[Vec<FaceData>],
    pool: &mut VertexPool,
) -> MeshResult<Vec<FaceData>> {
    csg_boolean_with_policy(op, meshes, pool, false)
}

fn csg_boolean_with_policy(
    op: BooleanOp,
    meshes: &[Vec<FaceData>],
    pool: &mut VertexPool,
    finalize_faces: bool,
) -> MeshResult<Vec<FaceData>> {
    match meshes {
        [] => Err(MeshError::EmptyBooleanResult {
            op: "BooleanUnion".to_string(),
        }),
        [mesh] => Ok(mesh.clone()),
        _ => execute_generalized_boolean(op, meshes, pool, finalize_faces),
    }
}

fn propagate_seam_vertices_until_stable(
    faces: &[FaceData],
    segs: &mut [Vec<SnapSegment>],
    pool: &VertexPool,
) {
    const MAX_PROPAGATION_PASSES: usize = 8;

    for _ in 0..MAX_PROPAGATION_PASSES {
        let before: usize = segs.iter().map(Vec::len).sum();
        propagate_seam_vertices(faces, segs, pool);
        let after: usize = segs.iter().map(Vec::len).sum();
        if after == before {
            break;
        }
    }
}

/// Execute the canonical arrangement pipeline after any decisive operand
/// relationships have been resolved.
fn execute_generalized_boolean(
    op: BooleanOp,
    meshes: &[Vec<FaceData>],
    pool: &mut VertexPool,
    finalize_faces: bool,
) -> MeshResult<Vec<FaceData>> {
    if let Some(short_circuit) = resolve_short_circuit_boolean(op, meshes, pool) {
        return short_circuit;
    }

    execute_arrangement_pass(op, meshes, pool, finalize_faces)
}

fn resolve_short_circuit_boolean(
    op: BooleanOp,
    meshes: &[Vec<FaceData>],
    pool: &mut VertexPool,
) -> Option<MeshResult<Vec<FaceData>>> {
    // ── Coplanar flat-mesh fast-path ─────────────────────────────────────────
    // If all N operands lie in the exact same 3D plane, dispatch directly
    // to the robust Sutherland-Hodgman 2D pipeline rather than full 3D intersections.
    let mut all_coplanar = true;
    let mut reference_basis: Option<crate::application::csg::coplanar::basis::PlaneBasis> = None;

    for mesh_faces in meshes.iter() {
        if mesh_faces.is_empty() {
            continue;
        }
        if let Some(basis) = crate::application::csg::coplanar::detect_flat_plane(mesh_faces, pool)
        {
            if let Some(ref ref_basis) = reference_basis {
                let origin_diff = basis.origin.coords - ref_basis.origin.coords;
                let cross = basis.normal.cross(&ref_basis.normal).norm();
                let dot_offset = origin_diff.dot(&ref_basis.normal).abs();
                if cross > 1e-6 || dot_offset > 1e-6 {
                    all_coplanar = false;
                    break;
                }
            } else {
                reference_basis = Some(basis);
            }
        } else {
            all_coplanar = false;
            break;
        }
    }

    if all_coplanar {
        if let Some(basis) = reference_basis {
            let mut current = meshes[0].clone();
            for next_faces in &meshes[1..] {
                current = crate::application::csg::coplanar::boolean_coplanar(
                    op, &current, next_faces, pool, &basis,
                );
                if current.is_empty() {
                    break;
                }
            }
            if current.is_empty() {
                return Some(Err(MeshError::EmptyBooleanResult {
                    op: format!("{op:?}"),
                }));
            }
            return Some(Ok(current));
        }
    }

    if let [faces_a, faces_b] = meshes {
        let result = match containment(faces_a, faces_b, pool) {
            Containment::BInsideA => match op {
                BooleanOp::Union => faces_a.clone(),
                BooleanOp::Intersection => faces_b.clone(),
                BooleanOp::Difference => return None,
            },
            Containment::AInsideB => match op {
                BooleanOp::Union => faces_b.clone(),
                BooleanOp::Intersection => faces_a.clone(),
                BooleanOp::Difference => Vec::new(),
            },
            Containment::Disjoint => match op {
                BooleanOp::Union => {
                    let mut result = Vec::with_capacity(faces_a.len() + faces_b.len());
                    result.extend_from_slice(faces_a);
                    result.extend_from_slice(faces_b);
                    result
                }
                BooleanOp::Intersection => Vec::new(),
                BooleanOp::Difference => faces_a.clone(),
            },
            Containment::Intersecting => return None,
        };

        return Some(non_empty_result(op, result));
    }

    None
}

#[inline]
fn non_empty_result(op: BooleanOp, result: Vec<FaceData>) -> MeshResult<Vec<FaceData>> {
    if result.is_empty() {
        Err(MeshError::EmptyBooleanResult {
            op: format!("{op:?}"),
        })
    } else {
        Ok(result)
    }
}

fn execute_arrangement_pass(
    op: BooleanOp,
    meshes: &[Vec<FaceData>],
    pool: &mut VertexPool,
    finalize_faces: bool,
) -> MeshResult<Vec<FaceData>> {
    let n_meshes = meshes.len();

    // ── Pre-compute global Mesh AABBs ──────────────────────────────────────────────
    let mut mesh_aabbs: Vec<Aabb> = Vec::with_capacity(n_meshes);
    for m in meshes {
        let mut bb = Aabb::empty();
        for f in m {
            bb.expand(pool.position(f.vertices[0]));
            bb.expand(pool.position(f.vertices[1]));
            bb.expand(pool.position(f.vertices[2]));
        }
        // Expand AABB slightly to defeat floating point precision misses on snapped vertices
        bb.min -= crate::domain::core::scalar::Vector3r::new(1e-6, 1e-6, 1e-6);
        bb.max += crate::domain::core::scalar::Vector3r::new(1e-6, 1e-6, 1e-6);
        mesh_aabbs.push(bb);
    }

    // ── Phase 1: generalized Boolean broad phase ─────────────────────────────────────────
    let mut aabbs: Vec<Vec<Aabb>> = Vec::with_capacity(n_meshes);
    for m in meshes {
        aabbs.push(m.iter().map(|f| triangle_aabb(f, pool)).collect());
    }

    let mut pairs = Vec::new();
    for i in 0..n_meshes {
        with_bvh(&aabbs[i], |tree, token| {
            for j in (i + 1)..n_meshes {
                if !mesh_aabbs[i].intersects(&mesh_aabbs[j]) {
                    continue;
                }
                let mut hits = Vec::new();
                for (fb, aabb_b) in aabbs[j].iter().enumerate() {
                    hits.clear();
                    tree.query_overlapping(aabb_b, &token, &mut hits);
                    for &fa in &hits {
                        pairs.push(BooleanCandidatePair {
                            mesh_a: i,
                            face_a: fa,
                            mesh_b: j,
                            face_b: fb,
                        });
                    }
                }
            }
        });
    }

    // ── Phase 2: Narrow Segment Intersect ───────────────────────────────────────────
    let mut segs: Vec<Vec<Vec<SnapSegment>>> =
        meshes.iter().map(|m| vec![Vec::new(); m.len()]).collect();
    let mut coplanar_pairs = Vec::new();

    for pair in &pairs {
        let fa = &meshes[pair.mesh_a][pair.face_a];
        let fb = &meshes[pair.mesh_b][pair.face_b];
        match intersect_triangles(fa, pool, fb, pool) {
            IntersectionType::Segment { start, end } => {
                let snap = SnapSegment { start, end };
                segs[pair.mesh_a][pair.face_a].push(snap);
                segs[pair.mesh_b][pair.face_b].push(snap);
            }
            IntersectionType::Coplanar => {
                coplanar_pairs.push((pair.mesh_a, pair.face_a, pair.mesh_b, pair.face_b));
            }
            IntersectionType::None => {}
        }
    }

    // ── Phase 2b.5: Propagate seam vertices
    for i in 0..n_meshes {
        propagate_seam_vertices_until_stable(&meshes[i], &mut segs[i], pool);
    }

    // ── Phase 2c: Coplanar Dispatches ───────────────────────────────────────────────
    let coplanar_phase = dispatch_boolean_coplanar(op, n_meshes, meshes, &coplanar_pairs, pool, &mut segs);

    // Coplanar-cap resolution can inject new seam segments onto rim triangles
    // after the initial propagation pass. Propagate once more so adjacent barrel
    // triangles sharing those edges receive matching Steiner vertices.
    for i in 0..n_meshes {
        propagate_seam_vertices_until_stable(&meshes[i], &mut segs[i], pool);
    }

    // ── Phase 3: Universal Co-refinement ────────────────────────────────────────────
    let mut frags: Vec<BooleanFragmentRecord> = Vec::new();

    for (m_idx, m_faces) in meshes.iter().enumerate() {
        append_corefined_fragments(
            &mut frags,
            m_faces,
            &coplanar_phase.skipped_faces_by_mesh[m_idx],
            &segs[m_idx],
            pool,
            |face, parent_idx| BooleanFragmentRecord {
                face,
                mesh_idx: m_idx,
                parent_idx,
            },
        );
    }

    // ── Phase 3.5/4: Generalized fragment consolidation and classification ─────────
    let mut result_faces = resolve_multi_mesh_fragments(
        op,
        &mut frags,
        meshes,
        &mesh_aabbs,
        &coplanar_phase.representative_faces,
        coplanar_phase.result_faces,
        pool,
    );

    // ── Phase 5: Global Patching ────────────────────────────────────────────────────

    if finalize_faces {
        finalize_boolean_faces(&mut result_faces, pool);
    }

    Ok(result_faces)
}
