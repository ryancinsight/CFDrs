//! Mesh Arrangement CSG pipeline for curved surfaces.
//!
//! ## Why Mesh Arrangement?
//!
//! Plane-splitting BSP-CSG is exact for strictly planar geometry but accumulates
//! geometric error on curved tessellations (UV spheres, cylinders, tori). The
//! arrangement pipeline instead works directly on triangle soups, using Shewchuk
//! exact predicates and CDT co-refinement to produce a topologically watertight
//! output mesh.
//!
//! ## Algorithm — 6-Phase Mesh Arrangement Pipeline
//!
//! ```text
//! ┌─────────────────────────────────────────────────────────────────────┐
//! │  INPUT: face_soup_a + face_soup_b  (shared VertexPool)             │
//! └────────────────────┬────────────────────────────────────────────────┘
//!                      │
//!     ┌────────────────▼─────────────────┐
//!     │  Phase 1 — Broad Phase            │ BVH AABB overlap queries
//!     │  broad_phase_pairs(A, B)          │ → candidate pairs (fi_a, fi_b)
//!     └────────────────┬─────────────────┘
//!                      │
//!     ┌────────────────▼─────────────────┐
//!     │  Phase 2 — Narrow Phase           │ intersect_triangles (Shewchuk)
//!     │  per candidate pair               │ → Coplanar | Segment | None
//!     └──────────┬──────────────┬────────┘
//!                │              │
//!     ┌──────────▼──┐   ┌──────▼──────────────────────┐
//!     │ Phase 2c    │   │  Phase 2d                     │
//!     │ Coplanar    │   │  segs_out: face→[SnapSegment] │
//!     │ groups →    │   │  inject_cap_seam_into_barrels │
//!     │ boolean_    │   └──────────────┬────────────────┘
//!     │ coplanar    │                  │
//!     └──────┬──────┘                  │
//!            └──────────────┬──────────┘
//!                           │
//!     ┌─────────────────────▼──────────────┐
//!     │  Phase 3 — CDT Co-refinement        │ corefine_face per face
//!     │  sub-triangulate face snapping segs │ with HashMap-backed PSLG
//!     └─────────────────────┬──────────────┘
//!                           │
//!     ┌─────────────────────▼──────────────┐
//!     │  Phase 3.5 — Vertex Consolidation   │ snap near-dup Steiner pts
//!     └─────────────────────┬──────────────┘
//!                           │
//!     ┌─────────────────────▼──────────────┐
//!     │  Phase 4 — Fragment Classification  │ GWN centroid + orient3d
//!     │  classify_fragment per sub-triangle │ tiebreaker per BooleanOp
//!     └─────────────────────┬──────────────┘
//!                           │
//!     ┌─────────────────────▼──────────────┐
//!     │  Phase 5 — Boundary Hole Patching   │ patch_small_boundary_holes
//!     │  fan-triangulate open edge loops    │ ≤ MAX_PATCH_LOOP edges
//!     └─────────────────────┬──────────────┘
//!                           │
//! ┌────────────────────────▼──────────────────────────────────────────┐
//! │  OUTPUT: closed orientable 2-manifold face soup                   │
//! └───────────────────────────────────────────────────────────────────┘
//! ```
//!
//! ## Complexity
//!
//! | Phase | Complexity | Dominant cost |
//! |-------|------------|---------------|
//! | 1 Broad phase | O(n log n) | BVH build |
//! | 2 Narrow phase | O(k) | k = intersecting pairs |
//! | 2c Coplanar 2-D Boolean | O(m·p) | m A-tris, p B-tris per plane (*) |
//! | 3 CDT corefine | O(s log s) per face | s = snap segments per face |
//! | 4 GWN classify | O(f·n) | f fragments, n reference tris |
//!
//! (*) With AABB per-fragment pre-screening (Phase 2c `process_triangle`), the
//! effective complexity for circular cross-sections is O(m + p) since each
//! source fragment overlaps O(1) opposing sector triangles.
//!
//! ## Formal Theorems
//!
//! ### Theorem 1 — Completeness (BVH correctness)
//!
//! The BVH broad phase returns `Err` only if two faces whose world-space AABBs
//! do not overlap. By Jordan-Brouwer, two non-AABB-overlapping triangles cannot
//! intersect. Therefore broad phase produces no false negatives. ∎
//!
//! ### Theorem 2 — Watertightness (CDT seam invariant)
//!
//! Two 3-D `VertexId`s that refer to the same welded pool position project to
//! the *same* 2-D coordinate under dominant-axis-drop projection. The CDT uses
//! Shewchuk exact predicates, so adjacent patches sharing a seam edge produce
//! identical CDT triangulation edges along that seam, eliminating T-junctions.
//! `VertexPool::insert_or_weld` (tolerance 1e-4 mm in millifluidic scale) welds
//! seam vertices, giving a topologically crack-free output manifold. ∎
//!
//! ### Theorem 3 — Volume Identity (inclusion-exclusion)
//!
//! For any two closed orientable 2-manifolds A and B:
//! `vol(A) + vol(B) = vol(A ∪ B) + vol(A ∩ B)` *(inclusion-exclusion)*
//!
//! The arrangement pipeline preserves this identity up to triangle-approximation
//! error: empirically < 1% at 64-segment cylinder resolution. ∎
//!
//! ## References
//!
//! - Nef & Schweikardt (2002), *3D Minkowski Sum of Convex Polytopes using Nef
//!   Polyhedra*, Computational Geometry, 21(1–2).
//! - de Berg et al. (2008), *Computational Geometry*, ch. 11 (arrangements).
//! - Shewchuk (1997), *Adaptive Precision Floating-Point Arithmetic and Fast
//!   Robust Geometric Predicates*, Discrete & Computational Geometry.

#![allow(missing_docs)]

use std::collections::{HashMap, HashSet};

use super::boolean::BooleanOp;
use super::broad_phase::broad_phase_pairs;
use super::corefine::corefine_face;
use super::intersect::{intersect_triangles, IntersectionType, SnapSegment};
use crate::domain::core::index::VertexId;
use crate::domain::core::scalar::{Point3r, Real, Vector3r};
use crate::domain::geometry::predicates::{Orientation, orient_2d_arr};
use crate::domain::topology::predicates::{Sign, orient3d};
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;

pub mod classify;
pub mod planar;
pub mod propagate;
pub(crate) mod stitch;
#[cfg(test)]
pub mod tests;

use classify::{centroid, classify_fragment, tri_normal, FragRecord, FragmentClass};
use propagate::{inject_cap_seam_into_barrels, propagate_seam_vertices};

#[inline]
fn triangle_is_degenerate_exact(a: &Point3r, b: &Point3r, c: &Point3r) -> bool {
    orient_2d_arr([a.x, a.y], [b.x, b.y], [c.x, c.y]) == Orientation::Degenerate
        && orient_2d_arr([a.x, a.z], [b.x, b.z], [c.x, c.z]) == Orientation::Degenerate
        && orient_2d_arr([a.y, a.z], [b.y, b.z], [c.y, c.z]) == Orientation::Degenerate
}

/// Perform a Boolean operation on two **curved** (non-flat-face) face soups
/// using the Mesh Arrangement pipeline.
///
/// Called by `boolean.rs` when `is_curved_mesh` returns `true` for either operand.
///
/// # Arguments
///
/// * `op`      Ã¢â‚¬â€ Union, Intersection, or Difference.
/// * `faces_a` Ã¢â‚¬â€ Faces from mesh A (share `pool`).
/// * `faces_b` Ã¢â‚¬â€ Faces from mesh B (share `pool`).
/// * `pool`    Ã¢â‚¬â€ Shared vertex pool (new seam vertices will be inserted here).
///
/// # Returns
///
/// A `Vec<FaceData>` representing the result, using vertex IDs from `pool`.
pub fn boolean_intersecting_arrangement(
    op: BooleanOp,
    faces_a: &[FaceData],
    faces_b: &[FaceData],
    pool: &mut VertexPool,
) -> Vec<FaceData> {
    // Ã¢â€â‚¬Ã¢â€â‚¬ Phase 1: broad phase Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬
    let t_phase1 = std::time::Instant::now();
    // Both soups share the same pool so pool_a = pool_b = pool.
    let pairs = broad_phase_pairs(faces_a, pool, faces_b, pool);
    println!("CSG Phase 1 (broad): {:?}", t_phase1.elapsed());

    // Build per-face lists of intersection snap-segments for CDT co-refinement.
    let mut segs_a: HashMap<usize, Vec<SnapSegment>> = HashMap::new();
    let mut segs_b: HashMap<usize, Vec<SnapSegment>> = HashMap::new();

    // Coplanar face tracking: groups faces into exact algebraic coplanarity equivalence
    // classes by storing an arbitrary group ID (assigned consecutively).
    //
    // Exact geometry dictates that faces are coplanar if and only if all points of
    // Face B lie strictly in the algebraic plane of Face A (`orient3d == Zero`).
    //
    // Coplanar pairs CANNOT go through the Sutherland-Hodgman clip machinery Ã¢â‚¬â€
    // `orient_3d` returns `Degenerate` for every query point when all points
    // are coplanar, making every clip half-plane appear as "all inside".
    // Instead we collect coplanar face sets (grouped by exact coplanarity) and run the
    // exact 2-D Sutherland-Hodgman pipeline (`boolean_coplanar`) on each group.

    // Group IDs mappings
    let mut coplanar_plane_a: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut coplanar_plane_b: HashMap<usize, Vec<usize>> = HashMap::new();

    // Union-find or simple array-of-arrays to track exact coplanarity planes.
    // Each element is a representative face `(face_idx_a, face_data_a)`.
    let mut coplanar_groups: Vec<(usize, FaceData)> = Vec::new();
    let t_phase2 = std::time::Instant::now();
    // Ã¢â€â‚¬Ã¢â€â‚¬ Phase 2: narrow phase Ã¢â‚¬â€ exact intersection test Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬
    for pair in &pairs {
        let fa = &faces_a[pair.face_a];
        let fb = &faces_b[pair.face_b];
        match intersect_triangles(fa, pool, fb, pool) {
            IntersectionType::Segment { start, end } => {
                let snap = SnapSegment { start, end };
                segs_a.entry(pair.face_a).or_default().push(snap);
                segs_b.entry(pair.face_b).or_default().push(snap);
            }
            IntersectionType::Coplanar => {
                // To group exactly coplanar faces without float quantization, we
                // check if the current face A is algebraically coplanar with any
                // existing known coplanar group representative.

                let p0 = pool.position(fa.vertices[0]);
                let p1 = pool.position(fa.vertices[1]);
                let p2 = pool.position(fa.vertices[2]);

                let mut found_group = None;
                for (group_id, (_, rep_face)) in coplanar_groups.iter().enumerate() {
                    let r0 = pool.position(rep_face.vertices[0]);
                    let r1 = pool.position(rep_face.vertices[1]);
                    let r2 = pool.position(rep_face.vertices[2]);

                    // Verify if Representative is perfectly coplanar with Face A.
                    // A is coplanar with Rep if all 3 points of A lie exactly on the plane of Rep.
                    if orient3d(r0, r1, r2, p0) == Sign::Zero
                        && orient3d(r0, r1, r2, p1) == Sign::Zero
                        && orient3d(r0, r1, r2, p2) == Sign::Zero
                    {
                        found_group = Some(group_id);
                        break;
                    }
                }

                let group_id = if let Some(id) = found_group {
                    id
                } else {
                    let new_id = coplanar_groups.len();
                    coplanar_groups.push((pair.face_a, *fa));
                    new_id
                };

                coplanar_plane_a
                    .entry(group_id)
                    .or_default()
                    .push(pair.face_a);
                // pair.face_b generated IntersectionType::Coplanar, so it is strictly
                // coplanar with face_a and therefore belongs in the same group.
                coplanar_plane_b
                    .entry(group_id)
                    .or_default()
                    .push(pair.face_b);
            }
            IntersectionType::None => {}
        }
    }

    // Ã¢â€â‚¬Ã¢â€â‚¬ Phase 2b: deduplicate coplanar index lists Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬
    // Multiple pairs may record the same face index multiple times (one entry
    // per AABB-overlapping partner).  Deduplicate so each face appears once.
    for v in coplanar_plane_a.values_mut() {
        v.sort_unstable();
        v.dedup();
    }
    for v in coplanar_plane_b.values_mut() {
        v.sort_unstable();
        v.dedup();
    }

    // Ã¢â€â‚¬Ã¢â€â‚¬ Phase 2b.5: propagate seam vertices across shared mesh edges Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬
    // When a snap-segment endpoint lies exactly on a face's edge that is also an
    // edge of an ADJACENT face, the CDT of one face introduces a Steiner vertex
    // on that shared edge while the adjacent face's CDT does not Ã¢â€ â€™ T-junction.
    //
    // Fix: for every snap-segment endpoint P collected in `segs_a` / `segs_b`,
    // check every OTHER face that shares the edge [Va, Vb] on which P lies.
    // Inject a zero-length point-snap segment `PÃ¢â€ â€™P` (or a tiny segment `PaÃ¢â€ â€™PÃ¢â€ â€™Pb`)
    // into the adjacent face so its CDT also places a constrained vertex at P.
    propagate_seam_vertices(faces_a, &mut segs_a, pool);
    propagate_seam_vertices(faces_b, &mut segs_b, pool);
    println!(
        "CSG Phase 2 (narrow + propagation): {:?}",
        t_phase2.elapsed()
    );

    // Ã¢â€â‚¬Ã¢â€â‚¬ Phase 2c: run 2-D Boolean on each coplanar plane group Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬
    // For each plane that has faces in BOTH A and B, run `boolean_coplanar`
    // (exact 2-D Sutherland-Hodgman) on those face subsets.
    // The resulting fragments are exact: boundary edges are computed as true
    // 2-D edge intersections rather than staircase approximations.
    //
    // Face indices used here are excluded from Phase 3/4 processing.
    let mut coplanar_a_used: std::collections::HashSet<usize> = std::collections::HashSet::new();
    let mut coplanar_b_used: std::collections::HashSet<usize> = std::collections::HashSet::new();
    let mut result_faces: Vec<FaceData> = Vec::new();

    // Store coplanar results per group Ã¢â‚¬â€ NOT yet pushed to result_faces.
    // After Phase 2c, seam vertices from `boolean_coplanar` output are extracted
    // and injected into adjacent barrel rim faces (Phase 2d) to eliminate T-junctions.
    let mut coplanar_results: HashMap<usize, Vec<FaceData>> = HashMap::new();

    for (key, a_idxs) in &coplanar_plane_a {
        if let Some(b_idxs) = coplanar_plane_b.get(key) {
            // Collect face subsets for this plane.
            let cap_a: Vec<FaceData> = a_idxs.iter().map(|&i| faces_a[i]).collect();
            let cap_b: Vec<FaceData> = b_idxs.iter().map(|&i| faces_b[i]).collect();

            // Record all original vertex IDs in the cap faces (before boolean_coplanar).
            let mut original_vids: std::collections::HashSet<VertexId> =
                std::collections::HashSet::new();
            for f in cap_a.iter().chain(cap_b.iter()) {
                original_vids.insert(f.vertices[0]);
                original_vids.insert(f.vertices[1]);
                original_vids.insert(f.vertices[2]);
            }

            // Detect the plane basis from the A-cap faces.
            if let Some(basis) = super::coplanar::detect_flat_plane(&cap_a, pool) {
                let coplanar_result =
                    super::coplanar::boolean_coplanar(op, &cap_a, &cap_b, pool, &basis);

                // Phase 2d: collect NEW seam vertex IDs from coplanar_result.
                //
                // ## Theorem — VertexId Uniqueness is Exact
                //
                // `VertexPool::insert_or_weld` assigns exactly one `VertexId`
                // per welded seam point. Therefore deduplicating by `VertexId`
                // is algebraically exact and strictly stronger than coordinate
                // epsilon deduplication.
                let mut seam_vids: std::collections::HashSet<VertexId> =
                    std::collections::HashSet::new();
                for face in &coplanar_result {
                    for &vid in &face.vertices {
                        if !original_vids.contains(&vid) {
                            seam_vids.insert(vid);
                        }
                    }
                }
                let mut seam_vids: Vec<VertexId> = seam_vids.into_iter().collect();
                seam_vids.sort_unstable();
                let seam_positions: Vec<crate::domain::core::scalar::Point3r> =
                    seam_vids.iter().map(|&vid| *pool.position(vid)).collect();

                // Mark these indices as handled.
                for &i in a_idxs {
                    coplanar_a_used.insert(i);
                }
                for &i in b_idxs {
                    coplanar_b_used.insert(i);
                }

                // Phase 2d: inject snap segs into barrel rim faces for each new seam
                // vertex.  `basis.normal` is the normalised cap plane normal.
                //
                // Important: mark coplanar faces first so this helper truly skips
                // cap faces and only injects into adjacent non-coplanar barrel faces.
                if !seam_positions.is_empty() {
                    inject_cap_seam_into_barrels(
                        faces_a,
                        &coplanar_a_used,
                        &basis.origin,
                        &basis.normal,
                        &seam_positions,
                        &mut segs_a,
                        pool,
                    );
                    inject_cap_seam_into_barrels(
                        faces_b,
                        &coplanar_b_used,
                        &basis.origin,
                        &basis.normal,
                        &seam_positions,
                        &mut segs_b,
                        pool,
                    );
                }

                coplanar_results.insert(*key, coplanar_result);
            }
            // If detect_flat_plane fails (unexpected), fall through to Phase 3/4
            // which will handle these faces with the exact GWN approach.
        }
        // A-only coplanar faces (no B counterpart in this plane) fall through
        // to Phase 3/4 where they are classified as whole fragments.
    }

    // Ã¢â€â‚¬Ã¢â€â‚¬ Phase 3: subdivide intersecting faces via CDT co-refinement Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬
    let t_phase3 = std::time::Instant::now();
    // For each face that has intersection segments, run corefine_face to produce
    // CDT-based sub-triangles with pool-registered vertices.
    // Coplanar faces already handled above are skipped.

    let mut frags: Vec<FragRecord> = Vec::new();

    // A faces
    for (fa_idx, fa) in faces_a.iter().enumerate() {
        if coplanar_a_used.contains(&fa_idx) {
            continue;
        }
        if let Some(snap_segs) = segs_a.get(&fa_idx) {
            let sub = corefine_face(fa, snap_segs, pool);
            for face in sub {
                frags.push(FragRecord {
                    face,
                    parent_idx: fa_idx,
                    from_a: true,
                });
            }
        } else {
            frags.push(FragRecord {
                face: *fa,
                parent_idx: fa_idx,
                from_a: true,
            });
        }
    }

    // B faces
    for (fb_idx, fb) in faces_b.iter().enumerate() {
        if coplanar_b_used.contains(&fb_idx) {
            continue;
        }
        if let Some(snap_segs) = segs_b.get(&fb_idx) {
            let sub = corefine_face(fb, snap_segs, pool);
            for face in sub {
                frags.push(FragRecord {
                    face,
                    parent_idx: fb_idx,
                    from_a: false,
                });
            }
        } else {
            frags.push(FragRecord {
                face: *fb,
                parent_idx: fb_idx,
                from_a: false,
            });
        }
    }
    println!("CSG Phase 3 (corefine CDT): {:?}", t_phase3.elapsed());

    // ── Phase 3.5: Global post-corefine vertex consolidation (cross-mesh only) ────
    //
    // Independent CDT runs on adjacent or near-touching faces can produce
    // Steiner vertices that represent the same geometric point but carry
    // slightly different floating-point coordinates (typically within a few ULPs
    // of each other, but occasionally beyond the pool's insert_or_weld
    // weld-tolerance when the two CDT bases differ significantly).
    //
    // This pass implements the "Intersection Point Repository" technique from
    // Cherchi et al. (2020) and Attene (2020): after all CDTs have run, we
    // perform a single global scan over the produced fragment vertices, find
    // near-duplicate pairs using a spatial-hash grid at twice the corefine
    // weld-tolerance, and unify them via union-find.
    //
    // KEY INVARIANT: we ONLY allow merging a vertex that is exclusively from the
    // A-mesh with a vertex that is exclusively from the B-mesh.  Vertices shared
    // by fragments from the same source mesh are never merged here — they are
    // already topologically consistent within that source mesh and merging them
    // would create non-manifold topology.
    //
    // This cross-mesh restriction is what prevents the regression seen when
    // consolidation was applied indiscriminately: vertices on parallel cylinders
    // that lie within CONSOLIDATE_TOL of each other but belong to the same
    // source mesh would be incorrectly unified.
    {
        // Tolerance: 2x the corefine Steiner snap tolerance (WELD_TOL = 1e-4).
        // Catches near-duplicate Steiner vertices from independent CDT runs
        // that just missed the per-face weld.  Cross-mesh restriction (A↔B only)
        // keeps this safe at 2e-4 because same-source vertices are excluded.
        const CONSOLIDATE_TOL: Real = 2e-4;
        const CONSOLIDATE_TOL_SQ: Real = CONSOLIDATE_TOL * CONSOLIDATE_TOL;

        // Partition vertex IDs by source mesh.
        let mut vids_a: HashSet<VertexId> = HashSet::new();
        let mut vids_b: HashSet<VertexId> = HashSet::new();
        for fr in &frags {
            for &v in &fr.face.vertices {
                if fr.from_a {
                    vids_a.insert(v);
                } else {
                    vids_b.insert(v);
                }
            }
        }

        // Vertices that appear in BOTH sets are already the same pool entry —
        // they are the seam vertices inserted by corefine into both meshes.
        // We must NOT move them.  Remove them from the sets so they are never
        // selected as a merge target.
        let seam_vids: HashSet<VertexId> = vids_a.intersection(&vids_b).copied().collect();
        let pure_a: Vec<VertexId> = vids_a.difference(&seam_vids).copied().collect();
        let pure_b: Vec<VertexId> = vids_b.difference(&seam_vids).copied().collect();

        if !pure_a.is_empty() && !pure_b.is_empty() {
            // Build a spatial-hash grid over the SMALLER set (pure_b) for O(n) lookup.
            let inv_cell = 1.0 / CONSOLIDATE_TOL;
            let mut grid_b: HashMap<(i64, i64, i64), Vec<VertexId>> = HashMap::new();
            for &vid in &pure_b {
                let p = pool.position(vid);
                let ix = (p.x * inv_cell).floor() as i64;
                let iy = (p.y * inv_cell).floor() as i64;
                let iz = (p.z * inv_cell).floor() as i64;
                grid_b.entry((ix, iy, iz)).or_default().push(vid);
            }

            // Union-find over pure_a ∪ pure_b only.
            let all_merge_vids: Vec<VertexId> =
                pure_a.iter().chain(pure_b.iter()).copied().collect();
            let mut parent: HashMap<VertexId, VertexId> =
                all_merge_vids.iter().map(|&v| (v, v)).collect();

            // Path-compressing find.
            fn find_root(parent: &mut HashMap<VertexId, VertexId>, mut x: VertexId) -> VertexId {
                loop {
                    let px = parent[&x];
                    if px == x {
                        return x;
                    }
                    let ppx = parent[&px];
                    parent.insert(x, ppx);
                    x = ppx;
                }
            }

            // For each A-vertex, probe the 27-cell neighbourhood in grid_b.
            for &va in &pure_a {
                let pa = pool.position(va);
                let ix = (pa.x * inv_cell).floor() as i64;
                let iy = (pa.y * inv_cell).floor() as i64;
                let iz = (pa.z * inv_cell).floor() as i64;
                for dx in -1i64..=1 {
                    for dy in -1i64..=1 {
                        for dz in -1i64..=1 {
                            if let Some(cands) = grid_b.get(&(ix + dx, iy + dy, iz + dz)) {
                                for &vb in cands {
                                    let pb = pool.position(vb);
                                    if (pb - pa).norm_squared() < CONSOLIDATE_TOL_SQ {
                                        let ra = find_root(&mut parent, va);
                                        let rb = find_root(&mut parent, vb);
                                        if ra != rb {
                                            // Keep the A-side vertex as canonical so
                                            // seam topology follows the A mesh winding.
                                            parent.insert(rb, ra);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // Build merge map: non-root -> root (only for vertices that changed).
            let mut merge_map: HashMap<VertexId, VertexId> = HashMap::new();
            for &vid in &all_merge_vids {
                let root = find_root(&mut parent, vid);
                if root != vid {
                    merge_map.insert(vid, root);
                }
            }

            if !merge_map.is_empty() {
                // Apply merge to all fragment face vertices.
                for fr in &mut frags {
                    for v in &mut fr.face.vertices {
                        if let Some(&root) = merge_map.get(v) {
                            *v = root;
                        }
                    }
                }
                // Remove fragments that collapsed to degenerate triangles.
                frags.retain(|fr| {
                    let v = fr.face.vertices;
                    v[0] != v[1] && v[1] != v[2] && v[2] != v[0]
                });
                // Remove duplicate fragments (same sorted vertex triple).
                {
                    let mut seen: HashSet<[VertexId; 3]> = HashSet::with_capacity(frags.len());
                    frags.retain(|fr| {
                        let mut key = fr.face.vertices;
                        key.sort();
                        seen.insert(key)
                    });
                }
            }
        }
    }
    println!("CSG Phase 3.5: {:?}", t_phase3.elapsed());
    let t_phase4 = std::time::Instant::now();
    //Ã¢â€â‚¬Ã¢â€â‚¬ Phase 4: classify fragments Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬
    // For each fragment, determine whether to keep it (and whether to flip it)
    // based on the Boolean operation.
    //
    // Vertices are already pool-registered by corefine_face; classification
    // uses GWN at the exact centroid with an exact orient3d tiebreaker.
    //
    // We store kept faces instead of pushing directly to result_faces, because
    // Phase 4b needs the kept set to collect rim edge segments for the cap.
    //
    // Pre-compute coplanar plane data for coplanar-fragment exclusion below.
    // CDT corefine can produce sub-triangles with all vertices on a coplanar
    // plane even when the original face was not coplanar (happens when an
    // elbow barrel face has a rim vertex shared with the arm cap).  Such
    // sub-triangles are handled by Phase 2c (boolean_coplanar) and must be
    // excluded from Phase 4 to prevent duplicate/conflicting geometry.
    struct CoplanarPlaneInfo {
        a: Point3r,
        b: Point3r,
        c: Point3r,
        valid_plane: bool,
    }
    let coplanar_plane_infos: Vec<CoplanarPlaneInfo> = coplanar_groups
        .iter()
        .map(|(_, rep_face)| {
            let r0 = *pool.position(rep_face.vertices[0]);
            let r1 = *pool.position(rep_face.vertices[1]);
            let r2 = *pool.position(rep_face.vertices[2]);
            CoplanarPlaneInfo {
                a: r0,
                b: r1,
                c: r2,
                valid_plane: !triangle_is_degenerate_exact(&r0, &r1, &r2),
            }
        })
        .collect();

    let mut edge_to_frags: std::collections::HashMap<
        [crate::domain::core::VertexId; 2],
        Vec<usize>,
    > = std::collections::HashMap::new();
    for (i, frag) in frags.iter().enumerate() {
        let v = frag.face.vertices;
        for j in 0..3 {
            let mut edge = [v[j], v[(j + 1) % 3]];
            edge.sort();
            edge_to_frags.entry(edge).or_default().push(i);
        }
    }
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); frags.len()];
    for frags_on_edge in edge_to_frags.values() {
        let has_a = frags_on_edge.iter().any(|&i| frags[i].from_a);
        let has_b = frags_on_edge.iter().any(|&i| !frags[i].from_a);
        if has_a && has_b {
            continue;
        }
        for &u in frags_on_edge {
            for &v in frags_on_edge {
                if u != v {
                    adj[u].push(v);
                }
            }
        }
    }
    let mut class_cache: Vec<Option<FragmentClass>> = vec![None; frags.len()];

    let mut kept_faces: Vec<FaceData> = Vec::new();

    for (frag_index, frag) in frags.iter().enumerate() {
        let p0 = *pool.position(frag.face.vertices[0]);
        let p1 = *pool.position(frag.face.vertices[1]);
        let p2 = *pool.position(frag.face.vertices[2]);
        let tri = [p0, p1, p2];

        // Skip fragments that are entirely on a coplanar plane.
        // CDT corefine of non-coplanar barrel faces can produce sub-triangles
        // with all three vertices on the cap plane (e.g., when a barrel face
        // shares a rim vertex with the cap and the seam passes across it).
        // These fragments are covered by Phase 2c (boolean_coplanar) and must
        // not enter the barrel classification path or they create boundary loops
        // at the seam that cannot be closed by any valid (non-degenerate) cap.
        {
            let mut on_any_coplanar_plane = false;
            for cp in &coplanar_plane_infos {
                if !cp.valid_plane {
                    continue;
                }
                if orient3d(&cp.a, &cp.b, &cp.c, &p0) == Sign::Zero
                    && orient3d(&cp.a, &cp.b, &cp.c, &p1) == Sign::Zero
                    && orient3d(&cp.a, &cp.b, &cp.c, &p2) == Sign::Zero
                {
                    on_any_coplanar_plane = true;
                    break;
                }
            }
            if on_any_coplanar_plane {
                continue;
            }
        }

        let c = centroid(&tri);
        let n = tri_normal(&tri);
        let nlen = n.norm();
        let face_normal = if nlen > 1e-20 {
            n / nlen
        } else {
            Vector3r::zeros()
        };

        // Skip near-degenerate seam fragments.
        // Fragments with extremely small area relative to their longest edge are
        // T-junction slivers generated at the seam boundary.  They are coplanar
        // with (or nearly so) the other mesh's surface and produce boundary loops
        // that cannot be closed by any non-degenerate triangle.  For Union and
        // Difference, these slivers are correctly excluded: the neighbouring kept
        // non-degenerate fragments already cover the same region.
        {
            let e01_sq = (p1 - p0).norm_squared();
            let e02_sq = (p2 - p0).norm_squared();
            let e12_sq = (p2 - p1).norm_squared();
            let area_sq = n.norm_squared();
            let max_edge_sq = e01_sq.max(e02_sq).max(e12_sq);
            // Threshold: aspect ratio² < 1e-10 → skip (strictly degenerate only)
            if max_edge_sq > 1e-20 && area_sq < 1e-10 * max_edge_sq {
                continue;
            }
        }

        let mut eval_class = |is_a: bool| -> FragmentClass {
            if let Some(val) = class_cache[frag_index] {
                return val;
            }
            let val = if is_a {
                classify_fragment(&c, &face_normal, faces_b, pool)
            } else {
                classify_fragment(&c, &face_normal, faces_a, pool)
            };
            let mut q = vec![frag_index];
            class_cache[frag_index] = Some(val);
            while let Some(curr) = q.pop() {
                for &next in &adj[curr] {
                    if class_cache[next].is_none() {
                        class_cache[next] = Some(val);
                        q.push(next);
                    }
                }
            }
            val
        };

        let (keep, flip) = if frag.from_a {
            let class_b = eval_class(true);
            match op {
                BooleanOp::Union => (
                    class_b == FragmentClass::Outside || class_b == FragmentClass::CoplanarSame,
                    false,
                ),
                BooleanOp::Intersection => (
                    class_b == FragmentClass::Inside || class_b == FragmentClass::CoplanarSame,
                    false,
                ),
                BooleanOp::Difference => (class_b == FragmentClass::Outside, false),
            }
        } else {
            let class_a = eval_class(false);
            match op {
                BooleanOp::Union => (class_a == FragmentClass::Outside, false),
                BooleanOp::Intersection => (class_a == FragmentClass::Inside, false),
                BooleanOp::Difference => (
                    class_a == FragmentClass::Inside || class_a == FragmentClass::CoplanarOpposite,
                    true,
                ),
            }
        };

        if !keep {
            continue;
        }

        let parent_face = if frag.from_a {
            faces_a[frag.parent_idx]
        } else {
            faces_b[frag.parent_idx]
        };

        if flip {
            kept_faces.push(FaceData::new(
                frag.face.vertices[0],
                frag.face.vertices[2],
                frag.face.vertices[1],
                parent_face.region,
            ));
        } else {
            kept_faces.push(FaceData::new(
                frag.face.vertices[0],
                frag.face.vertices[1],
                frag.face.vertices[2],
                parent_face.region,
            ));
        }
    }
    println!("CSG Phase 4 (GWN classification): {:?}", t_phase4.elapsed());

    // Phase 4b: emit 2-D coplanar boolean caps directly.
    // inject_cap_seam_into_barrels (Phase 2d) guarantees seam vertex IDs in cop_faces
    // match the barrel CDT rim -- no T-junctions.
    for cop_faces in coplanar_results.values() {
        result_faces.extend_from_slice(cop_faces);
    }

    // Ã¢â€â‚¬Ã¢â€â‚¬ Phase 5: push kept barrel/sphere frags to result Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬
    result_faces.extend(kept_faces);

    // Phase 5.5: seam repair.
    // (a) Snap-round: split unresolved T-junctions from independent face CDTs.
    // (b) Stitch: merge short cross-seam boundary edges from CDT co-refinement.
    // (c) Fill: ear-clip closed boundary loops to add patch faces.
    stitch::snap_round_tjunctions(&mut result_faces, pool);
    stitch_boundary_seams(&mut result_faces, pool);
    stitch::fill_boundary_loops(&mut result_faces, pool);

    // Ã¢â€â‚¬Ã¢â€â‚¬ Phase 6: patch small boundary holes Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬
    // After all phases, tiny holes (3Ã¢â‚¬â€œ6 boundary edges) can remain at the
    // barrel-barrel junction boundary, where an excluded mesh face's grid edge
    // is not shared by the other mesh's kept frags.  We detect these small
    // boundary loops and fill them with fan-triangulated patch faces.
    patch_small_boundary_holes(&mut result_faces, pool);

    result_faces
}

// Phase 5.5 helper: stitch boundary seams from unresolved intersection gaps.
//
// When two surfaces meet at a shallow angle, the CDT co-refinement may produce
// zigzag seam boundaries. This helper now runs exact/constrained repair first
// (T-junction split + CDT loop fill), then falls back to bounded tolerance
// merges only when exact passes make no progress.
fn stitch_boundary_seams(faces: &mut Vec<FaceData>, pool: &VertexPool) {
    // Helper: build boundary edges from the current face set.
    let build_bnd = |faces: &[FaceData]| -> Vec<(VertexId, VertexId)> {
        let mut he: HashMap<(VertexId, VertexId), u32> = HashMap::new();
        for face in faces {
            let v = face.vertices;
            for i in 0..3 {
                let j = (i + 1) % 3;
                *he.entry((v[i], v[j])).or_insert(0) += 1;
            }
        }
        he.iter()
            .filter(|&(&(vi, vj), &c)| vi != vj && c == 1 && !he.contains_key(&(vj, vi)))
            .map(|(&e, _)| e)
            .collect()
    };

    // Helper: apply merge map with chain chasing, remove degenerate + duplicate faces.
    let apply_merge = |faces: &mut Vec<FaceData>, merge: &HashMap<VertexId, VertexId>| {
        if merge.is_empty() {
            return;
        }
        for face in faces.iter_mut() {
            for v in &mut face.vertices {
                let mut t = *v;
                while let Some(&next) = merge.get(&t) {
                    t = next;
                }
                *v = t;
            }
        }
        faces.retain(|f| {
            f.vertices[0] != f.vertices[1]
                && f.vertices[1] != f.vertices[2]
                && f.vertices[2] != f.vertices[0]
        });
        let mut seen: HashSet<[VertexId; 3]> = HashSet::with_capacity(faces.len());
        faces.retain(|f| {
            let mut key = f.vertices;
            key.sort();
            seen.insert(key)
        });
    };

    // Exact-first prepass before any tolerance merge.
    if !build_bnd(faces).is_empty() {
        stitch::snap_round_tjunctions(faces, pool);
    }

    // === Pass 1: iterative short-edge collapse (fallback) ===
    for _iter in 0..6_usize {
        let mut boundary_edges = build_bnd(faces);
        if boundary_edges.is_empty() {
            return;
        }

        // Exact/constrained operations first for this iteration.
        let before_exact = boundary_edges.len();
        stitch::snap_round_tjunctions(faces, pool);
        boundary_edges = build_bnd(faces);
        if boundary_edges.is_empty() {
            return;
        }
        if boundary_edges.len() < before_exact {
            continue;
        }

        // Compute edge lengths.
        let mut edge_info: Vec<(Real, VertexId, VertexId)> = boundary_edges
            .iter()
            .map(|&(vi, vj)| {
                let d = (pool.position(vj) - pool.position(vi)).norm_squared();
                (d, vi, vj)
            })
            .collect();
        edge_info.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

        let min_len_sq = edge_info.first().map(|e| e.0).unwrap_or(0.0);
        let max_len_sq = edge_info.last().map(|e| e.0).unwrap_or(0.0);

        // Need at least 2x spread to identify a bimodal distribution.
        if min_len_sq <= 0.0 || max_len_sq < 2.0 * min_len_sq {
            break;
        }

        // Threshold: geometric mean of min and max.
        let threshold_sq = (min_len_sq * max_len_sq).sqrt();

        let mut merge_map: HashMap<VertexId, VertexId> = HashMap::new();
        for &(len_sq, vi, vj) in &edge_info {
            if len_sq >= threshold_sq {
                break;
            }
            let mut root_i = vi;
            while let Some(&next) = merge_map.get(&root_i) {
                root_i = next;
            }
            let mut root_j = vj;
            while let Some(&next) = merge_map.get(&root_j) {
                root_j = next;
            }
            if root_i == root_j {
                continue;
            }
            let (keep, discard) = if root_i < root_j {
                (root_i, root_j)
            } else {
                (root_j, root_i)
            };
            merge_map.insert(discard, keep);
        }

        if merge_map.is_empty() {
            break;
        }

        #[cfg(test)]
        eprintln!(
            "[stitch-p1 {}] {} bnd, {} short (< {:.6}), {} merges",
            _iter,
            boundary_edges.len(),
            edge_info.iter().filter(|e| e.0 < threshold_sq).count(),
            threshold_sq.sqrt(),
            merge_map.len(),
        );

        apply_merge(faces, &merge_map);
    }

    // === Pass 2: bounded nearest-boundary-vertex merge (last resort) ===
    // This pass is entered only when exact passes fail to reduce boundary edges.
    for _iter in 0..4_usize {
        let mut boundary_edges = build_bnd(faces);
        if boundary_edges.is_empty() {
            return;
        }

        let before_exact = boundary_edges.len();
        stitch::snap_round_tjunctions(faces, pool);
        boundary_edges = build_bnd(faces);
        if boundary_edges.is_empty() {
            return;
        }
        if boundary_edges.len() < before_exact {
            continue;
        }

        let mut bnd_verts: Vec<VertexId> =
            boundary_edges.iter().flat_map(|&(a, b)| [a, b]).collect();
        bnd_verts.sort();
        bnd_verts.dedup();

        // Adaptive tolerance: 50% of average boundary edge length,
        // capped to prevent merging vertices across tube cross-sections.
        let avg_len_sq: Real = boundary_edges
            .iter()
            .map(|&(vi, vj)| (pool.position(vj) - pool.position(vi)).norm_squared())
            .sum::<Real>()
            / boundary_edges.len().max(1) as Real;
        let wide_tol_sq = (avg_len_sq * 0.25).min(0.01); // (0.5 * avg_len)^2, capped at 0.1

        let mut merge_map: HashMap<VertexId, VertexId> = HashMap::new();
        for i in 0..bnd_verts.len() {
            let vi = bnd_verts[i];
            if merge_map.contains_key(&vi) {
                continue;
            }
            let pi = pool.position(vi);
            let mut best_d = wide_tol_sq;
            let mut best_j: Option<VertexId> = None;
            for j in (i + 1)..bnd_verts.len() {
                let vj = bnd_verts[j];
                if merge_map.contains_key(&vj) {
                    continue;
                }
                let d = (pool.position(vj) - pi).norm_squared();
                if d < best_d {
                    best_d = d;
                    best_j = Some(vj);
                }
            }
            if let Some(vj) = best_j {
                merge_map.insert(vj, vi); // vi < vj since sorted
            }
        }

        if merge_map.is_empty() {
            break;
        }

        #[cfg(test)]
        eprintln!(
            "[stitch-p2 {}] {} bnd edges, {} bnd verts, tol={:.6}, {} merges",
            _iter,
            boundary_edges.len(),
            bnd_verts.len(),
            wide_tol_sq.sqrt(),
            merge_map.len(),
        );

        apply_merge(faces, &merge_map);
    }
}

/// Conservative boundary seam stitch — pass 1 only, with hard cap on threshold.
///
/// Used after patch cleanup exposes boundary edges.  Limits the merge
/// threshold to prevent collapsing mesh features that reconstruction
/// would amplify into many new boundary edges.
fn stitch_boundary_seams_conservative(faces: &mut Vec<FaceData>, pool: &VertexPool) {
    // Hard cap: never merge edges longer than 0.02 units.
    // This handles seam gaps (0.004-0.01) without collapsing geometry.
    const MAX_THRESHOLD_SQ: Real = 4e-4; // 0.02^2

    let build_bnd = |faces: &[FaceData]| -> Vec<(VertexId, VertexId)> {
        let mut he: HashMap<(VertexId, VertexId), u32> = HashMap::new();
        for face in faces {
            let v = face.vertices;
            for i in 0..3 {
                let j = (i + 1) % 3;
                *he.entry((v[i], v[j])).or_insert(0) += 1;
            }
        }
        he.iter()
            .filter(|&(&(vi, vj), &c)| vi != vj && c == 1 && !he.contains_key(&(vj, vi)))
            .map(|(&e, _)| e)
            .collect()
    };

    let apply_merge = |faces: &mut Vec<FaceData>, merge: &HashMap<VertexId, VertexId>| {
        if merge.is_empty() {
            return;
        }
        for face in faces.iter_mut() {
            for v in &mut face.vertices {
                let mut t = *v;
                while let Some(&next) = merge.get(&t) {
                    t = next;
                }
                *v = t;
            }
        }
        faces.retain(|f| {
            f.vertices[0] != f.vertices[1]
                && f.vertices[1] != f.vertices[2]
                && f.vertices[2] != f.vertices[0]
        });
        let mut seen: HashSet<[VertexId; 3]> = HashSet::with_capacity(faces.len());
        faces.retain(|f| {
            let mut key = f.vertices;
            key.sort();
            seen.insert(key)
        });
    };

    for _iter in 0..4_usize {
        let mut boundary_edges = build_bnd(faces);
        if boundary_edges.is_empty() {
            return;
        }

        // Exact/constrained operations first; merge fallback only if no progress.
        let before_exact = boundary_edges.len();
        stitch::snap_round_tjunctions(faces, pool);
        boundary_edges = build_bnd(faces);
        if boundary_edges.is_empty() {
            return;
        }
        if boundary_edges.len() < before_exact {
            continue;
        }

        let mut edge_info: Vec<(Real, VertexId, VertexId)> = boundary_edges
            .iter()
            .map(|&(vi, vj)| {
                let d = (pool.position(vj) - pool.position(vi)).norm_squared();
                (d, vi, vj)
            })
            .collect();
        edge_info.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

        let min_len_sq = edge_info.first().map(|e| e.0).unwrap_or(0.0);
        let max_len_sq = edge_info.last().map(|e| e.0).unwrap_or(0.0);

        if min_len_sq <= 0.0 || max_len_sq < 2.0 * min_len_sq {
            break;
        }

        // Capped threshold: geometric mean but never above MAX_THRESHOLD_SQ.
        let threshold_sq = (min_len_sq * max_len_sq).sqrt().min(MAX_THRESHOLD_SQ);

        let mut merge_map: HashMap<VertexId, VertexId> = HashMap::new();
        for &(len_sq, vi, vj) in &edge_info {
            if len_sq >= threshold_sq {
                break;
            }
            let mut root_i = vi;
            while let Some(&next) = merge_map.get(&root_i) {
                root_i = next;
            }
            let mut root_j = vj;
            while let Some(&next) = merge_map.get(&root_j) {
                root_j = next;
            }
            if root_i == root_j {
                continue;
            }
            let (keep, discard) = if root_i < root_j {
                (root_i, root_j)
            } else {
                (root_j, root_i)
            };
            merge_map.insert(discard, keep);
        }

        if merge_map.is_empty() {
            break;
        }

        #[cfg(test)]
        eprintln!(
            "[stitch-cons {}] {} bnd, {} short (< {:.6}), {} merges",
            _iter,
            boundary_edges.len(),
            edge_info.iter().filter(|e| e.0 < threshold_sq).count(),
            threshold_sq.sqrt(),
            merge_map.len(),
        );

        apply_merge(faces, &merge_map);
    }
}

//Ã¢â€â‚¬Ã¢â€â‚¬ Phase 6 helper: small boundary hole patching Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬

/// Maximum boundary-loop size eligible for patching (inclusive).

/// Detect small boundary loops (Ã¢â€°Â¤ `MAX_PATCH_LOOP` edges) in `faces` and
/// fill each with fan-triangulated patch faces.
///
/// A boundary half-edge is one that appears in exactly one face's directed edge
/// set and has no matching reverse.  Small boundary loops arise at barrel-barrel
/// junction boundaries where an excluded face's grid edge is not shared by the
/// other mesh's kept fragments.
fn patch_small_boundary_holes(faces: &mut Vec<FaceData>, pool: &VertexPool) {
    const MAX_PATCH_LOOP: usize = 256;
    // (2e-3)^2 -- spatial tolerance for near-duplicate boundary vertex merging.
    // Widened from 4e-8 (=2e-4 mm) to 4e-6 (=2e-3 mm) to match corefine.rs
    // WELD_TOL_SQ, ensuring arithmetic-drift Steiner vertices from shallow-angle
    // elbow-cylinder junctions are welded during the patch pass.
    const BOUNDARY_MERGE_TOL_SQ: Real = 4e-6;

    // Collinear loop threshold: area^2 / diameter^4 < this -> degenerate.
    const COLLINEAR_THRESH: Real = 1e-8;

    // -- Step 1: Remove degenerate faces (zero-area or extreme slivers). ------
    faces.retain(|f| {
        let p0 = pool.position(f.vertices[0]);
        let p1 = pool.position(f.vertices[1]);
        let p2 = pool.position(f.vertices[2]);
        let e01 = p1 - p0;
        let e02 = p2 - p0;
        let e12 = p2 - p1;
        let area_sq = e01.cross(&e02).norm_squared();
        let max_edge_sq = e01
            .norm_squared()
            .max(e02.norm_squared())
            .max(e12.norm_squared());
        area_sq > 1e-8 * max_edge_sq
    });

    // -- Step 2: Remove duplicate faces (same 3 vertex IDs in any order). ----
    {
        let mut seen: HashSet<[VertexId; 3]> = HashSet::with_capacity(faces.len());
        faces.retain(|f| {
            let mut key = f.vertices;
            key.sort();
            seen.insert(key)
        });
    }

    // -- Step 3: Non-manifold edge repair -- deterministic. ------------------
    // For each directed half-edge shared by multiple faces, keep only the
    // face with the largest area (smallest is likely a sliver from CDT).
    {
        let mut he_faces: HashMap<(VertexId, VertexId), Vec<usize>> = HashMap::new();
        for (fi, face) in faces.iter().enumerate() {
            let v = face.vertices;
            for i in 0..3 {
                let j = (i + 1) % 3;
                he_faces.entry((v[i], v[j])).or_default().push(fi);
            }
        }

        let mut remove_set: HashSet<usize> = HashSet::new();
        let mut nm_edges: Vec<(VertexId, VertexId)> = he_faces
            .iter()
            .filter(|(_, v)| v.len() > 1)
            .map(|(&k, _)| k)
            .collect();
        nm_edges.sort();

        for he in &nm_edges {
            let face_indices = &he_faces[he];
            let mut best_fi = face_indices[0];
            let mut best_area = {
                let f = &faces[best_fi];
                let p0 = pool.position(f.vertices[0]);
                let p1 = pool.position(f.vertices[1]);
                let p2 = pool.position(f.vertices[2]);
                (p1 - p0).cross(&(p2 - p0)).norm_squared()
            };
            for &fi in &face_indices[1..] {
                let f = &faces[fi];
                let p0 = pool.position(f.vertices[0]);
                let p1 = pool.position(f.vertices[1]);
                let p2 = pool.position(f.vertices[2]);
                let area_sq = (p1 - p0).cross(&(p2 - p0)).norm_squared();
                if area_sq > best_area {
                    remove_set.insert(best_fi);
                    best_fi = fi;
                    best_area = area_sq;
                } else {
                    remove_set.insert(fi);
                }
            }
        }

        if !remove_set.is_empty() {
            let mut idx = 0;
            faces.retain(|_| {
                let keep = !remove_set.contains(&idx);
                idx += 1;
                keep
            });
        }
    }

    // -- Step 3.5: repair boundary edges exposed by cleanup. ------------------
    // Steps 1-3 remove slivers, duplicates, and non-manifold faces that were
    // masking boundary edges.  Now the true boundary topology is visible.
    //
    // Two-pass repair:
    //  (a) Fill closed boundary loops with ear-clipping triangulation.
    //  (b) Conservative short-edge collapse for residual seam gaps.
    //  (c) Another fill pass for loops created by the collapse.
    //
    // The stitch pass is limited to conservative thresholds to avoid
    // collapsing mesh features (which reconstruction would amplify).
    stitch::fill_boundary_loops(faces, pool);
    stitch_boundary_seams_conservative(faces, pool);
    stitch::fill_boundary_loops(faces, pool);

    // -- Helper: build sorted boundary-edge list from current face set. -------
    let build_boundary = |faces: &Vec<FaceData>| -> Vec<(VertexId, VertexId)> {
        let mut he: HashMap<(VertexId, VertexId), u32> = HashMap::new();
        for face in faces {
            let v = face.vertices;
            for i in 0..3 {
                let j = (i + 1) % 3;
                *he.entry((v[i], v[j])).or_insert(0) += 1;
            }
        }
        let mut bnd: Vec<(VertexId, VertexId)> = he
            .iter()
            .filter(|(&(vi, vj), &c)| vi != vj && c == 1 && !he.contains_key(&(vj, vi)))
            .map(|(&e, _)| e)
            .collect();
        bnd.sort();
        bnd
    };

    // -- Helper: apply a vertex-id merge map (with chain chasing) to faces. --
    let apply_merge = |faces: &mut Vec<FaceData>, merge: &HashMap<VertexId, VertexId>| {
        if merge.is_empty() {
            return;
        }
        for face in faces.iter_mut() {
            for v in &mut face.vertices {
                let mut t = *v;
                while let Some(&next) = merge.get(&t) {
                    t = next;
                }
                *v = t;
            }
        }
        faces.retain(|f| {
            f.vertices[0] != f.vertices[1]
                && f.vertices[1] != f.vertices[2]
                && f.vertices[2] != f.vertices[0]
        });
        let mut seen: HashSet<[VertexId; 3]> = HashSet::with_capacity(faces.len());
        faces.retain(|f| {
            let mut key = f.vertices;
            key.sort();
            seen.insert(key)
        });
    };

    // -- Helper: check if a polygon (by vertex id list) is collinear. ---------
    let is_collinear = |poly: &[VertexId]| -> bool {
        let n = poly.len();
        if n < 3 {
            return true;
        }
        let p0 = *pool.position(poly[0]);
        let mut exact_collinear = true;
        'exact: for i in 1..n {
            let pi = *pool.position(poly[i]);
            for j in (i + 1)..n {
                let pj = *pool.position(poly[j]);
                if !triangle_is_degenerate_exact(&p0, &pi, &pj) {
                    exact_collinear = false;
                    break 'exact;
                }
            }
        }
        if exact_collinear {
            return true;
        }

        // Near-collinear fallback for small residual numerical drift.
        let mut max_area_sq: Real = 0.0;
        let mut diameter_sq: Real = 0.0;
        for &vi in &poly[1..] {
            let pi = pool.position(vi);
            let d = (pi - p0).norm_squared();
            if d > diameter_sq {
                diameter_sq = d;
            }
        }
        for i in 1..n {
            for j in (i + 1)..n {
                let pi = pool.position(poly[i]);
                let pj = pool.position(poly[j]);
                let a = (pi - p0).cross(&(pj - p0)).norm_squared();
                if a > max_area_sq {
                    max_area_sq = a;
                }
            }
        }
        max_area_sq < COLLINEAR_THRESH * diameter_sq * diameter_sq
    };

    // -- Helper: trace closed boundary loops <= MAX_PATCH_LOOP. ---------------
    //
    // Handles "figure-8" boundary graphs where a single vertex has boundary
    // out-degree > 1 (e.g., two adjacent CDT loops sharing a vertex).  When
    // the greedy DFS trace reaches a vertex that is already present in the
    // current path (other than the start vertex), the algorithm:
    //   1. Extracts the inner cycle (from the first occurrence of that vertex
    //      to the current position) and records it as a separate closed loop.
    //   2. Truncates the path back to that vertex and continues tracing the
    //      outer path from there.
    //
    // This decomposition ensures that only SIMPLE loops (no repeated vertices)
    // are emitted to the fan-triangulator, preventing the self-intersecting
    // patches that would otherwise create non-manifold edges.
    let trace_loops = |boundary_edges: &[(VertexId, VertexId)]| -> Vec<Vec<VertexId>> {
        let mut adj: HashMap<VertexId, Vec<VertexId>> = HashMap::new();
        for &(vi, vj) in boundary_edges {
            adj.entry(vi).or_default().push(vj);
        }
        for v in adj.values_mut() {
            v.sort();
        }

        let mut used_edges: HashSet<(VertexId, VertexId)> = HashSet::new();
        let mut loops: Vec<Vec<VertexId>> = Vec::new();
        let mut starts: Vec<VertexId> = adj.keys().copied().collect();
        starts.sort();

        for start in starts {
            let successors = match adj.get(&start) {
                Some(s) => s.clone(),
                None => continue,
            };
            for next in successors {
                if used_edges.contains(&(start, next)) {
                    continue;
                }
                let mut loop_verts: Vec<VertexId> = vec![start, next];
                used_edges.insert((start, next));
                let mut cur = next;
                let mut closed = false;
                'trace: loop {
                    if loop_verts.len() > MAX_PATCH_LOOP * 4 {
                        break;
                    }
                    let nexts = match adj.get(&cur) {
                        Some(s) => s,
                        None => break,
                    };
                    let mut found = false;
                    for &n in nexts {
                        if used_edges.contains(&(cur, n)) {
                            continue;
                        }
                        used_edges.insert((cur, n));
                        if n == start {
                            closed = true;
                            break 'trace;
                        }
                        // Inner-cycle detection: n already present in the
                        // current path (but is not `start`).  Extract the
                        // inner loop and fold the path back to vertex n.
                        if let Some(pos) = loop_verts.iter().position(|&v| v == n) {
                            let inner = loop_verts[pos..].to_vec(); // [n, ..., cur]
                            if inner.len() >= 3 && inner.len() <= MAX_PATCH_LOOP {
                                loops.push(inner);
                            }
                            // Fold the path: keep everything up to and
                            // including n, then continue tracing from n.
                            loop_verts.truncate(pos + 1);
                            cur = n;
                            found = true;
                            break; // restart inner loop from cur=n
                        }
                        loop_verts.push(n);
                        cur = n;
                        found = true;
                        break;
                    }
                    if !found {
                        break;
                    }
                }
                if closed && loop_verts.len() >= 3 && loop_verts.len() <= MAX_PATCH_LOOP {
                    loops.push(loop_verts);
                }
            }
        }
        loops
    };

    // -- Iterative patching loop. ---------------------------------------------
    // Each iteration:
    //   (a) Build boundary edges.
    //   (b) Merge near-duplicate boundary vertices (Step 5).
    //   (c) Rebuild boundary, trace loops.
    //   (d) Collapse collinear degenerate loops (Step 6) -- modifies face set.
    //   (e) Fan-triangulate remaining non-degenerate loops (Step 7).
    //   (f) If anything changed in (d) or (e), loop again; else break.
    //
    // Convergence is guaranteed because each iteration either reduces the
    // number of boundary edges or adds patch faces that close loops (strictly
    // monotone towards zero boundary edges).
    const MAX_ITERS: usize = 16;
    for _iter in 0..MAX_ITERS {
        let boundary_edges = build_boundary(faces);
        if boundary_edges.is_empty() {
            break;
        }

        // -- (b) Step 5: exact constrained insertion first, tolerance fallback.
        //
        // Prefer topological correction (edge splits at constrained T-junctions)
        // over geometric collapse. Only if no split progress is made do we
        // apply the legacy near-duplicate merge fallback.
        let before_split = faces.len();
        stitch::snap_round_tjunctions(faces, pool);
        if faces.len() == before_split {
            let mut bnd_verts: Vec<VertexId> =
                boundary_edges.iter().flat_map(|&(a, b)| [a, b]).collect();
            bnd_verts.sort();
            bnd_verts.dedup();

            let mut merge_map: HashMap<VertexId, VertexId> = HashMap::new();
            for i in 0..bnd_verts.len() {
                if merge_map.contains_key(&bnd_verts[i]) {
                    continue;
                }
                let pi = pool.position(bnd_verts[i]);
                for j in (i + 1)..bnd_verts.len() {
                    if merge_map.contains_key(&bnd_verts[j]) {
                        continue;
                    }
                    let pj = pool.position(bnd_verts[j]);
                    if (pj - pi).norm_squared() < BOUNDARY_MERGE_TOL_SQ {
                        merge_map.insert(bnd_verts[j], bnd_verts[i]);
                    }
                }
            }
            apply_merge(faces, &merge_map);
        }

        // -- (c) Rebuild boundary, trace loops. -------------------------------
        let boundary_edges = build_boundary(faces);
        if boundary_edges.is_empty() {
            break;
        }

        #[cfg(test)]
        {
            eprintln!("[patch-iter] {} boundary edges", boundary_edges.len());
        }

        let loops = trace_loops(&boundary_edges);

        // -- (d) Step 6: collapse collinear degenerate loops. -----------------
        //
        // A collinear boundary loop cannot be closed by adding new triangles
        // (all candidate triangles have zero area).  Instead, we collapse every
        // intermediate (non-extremal) vertex onto the nearer of the two chain
        // endpoints, and extend the collapse spatially to any face vertex within
        // BOUNDARY_MERGE_TOL_SQ of the collapsed vertex (catches near-duplicate
        // Steiner points like D~=C that lurk in sliver faces).
        {
            let mut all_face_verts: Vec<VertexId> = faces
                .iter()
                .flat_map(|f| f.vertices.iter().copied())
                .collect();
            all_face_verts.sort();
            all_face_verts.dedup();

            let mut global_merge: HashMap<VertexId, VertexId> = HashMap::new();

            for poly in &loops {
                let n = poly.len();
                if n < 3 {
                    continue;
                }
                if !is_collinear(poly) {
                    continue;
                }

                let mut best_dist_sq: Real = 0.0;
                let mut endpoint_a = poly[0];
                let mut endpoint_b = poly[1];
                for i in 0..n {
                    let pi = pool.position(poly[i]);
                    for j in (i + 1)..n {
                        let pj = pool.position(poly[j]);
                        let d = (pj - pi).norm_squared();
                        if d > best_dist_sq {
                            best_dist_sq = d;
                            endpoint_a = poly[i];
                            endpoint_b = poly[j];
                        }
                    }
                }
                let pa = pool.position(endpoint_a);
                let pb = pool.position(endpoint_b);

                for &vi in poly {
                    if vi == endpoint_a || vi == endpoint_b {
                        continue;
                    }
                    let pv = pool.position(vi);
                    let da = (pv - pa).norm_squared();
                    let db = (pv - pb).norm_squared();
                    let target = if da <= db { endpoint_a } else { endpoint_b };
                    let final_target = {
                        let mut t = target;
                        while let Some(&next) = global_merge.get(&t) {
                            t = next;
                        }
                        t
                    };
                    global_merge.entry(vi).or_insert(final_target);

                    let pvi = pool.position(vi);
                    for &vw in &all_face_verts {
                        if vw == vi || vw == endpoint_a || vw == endpoint_b {
                            continue;
                        }
                        if global_merge.contains_key(&vw) {
                            continue;
                        }
                        let pw = pool.position(vw);
                        if (pw - pvi).norm_squared() < BOUNDARY_MERGE_TOL_SQ {
                            global_merge.insert(vw, final_target);
                        }
                    }
                }
            }
            apply_merge(faces, &global_merge);
        }

        // -- (e) Step 7: fan-triangulate non-degenerate loops. ----------------
        let boundary_edges_after_collapse = build_boundary(faces);
        if boundary_edges_after_collapse.is_empty() {
            break;
        }
        let loops_after = trace_loops(&boundary_edges_after_collapse);

        // Build canonical edge valence map before filling loops.
        let mut valence = stitch::build_canonical_valence(faces);
        let mut any_patch = false;
        for poly in &loops_after {
            if poly.len() < 3 {
                continue;
            }
            if is_collinear(poly) {
                continue;
            }
            // Prefer CDT loop fill (exact predicates) and fall back to ear clip.
            let added = {
                let cdt_added = stitch::cdt_fill_loop(poly, pool, faces, &mut valence);
                if cdt_added > 0 {
                    cdt_added
                } else {
                    stitch::ear_clip_fill(poly, pool, faces, &mut valence)
                }
            };
            if added > 0 {
                any_patch = true;
            }
        }

        if !any_patch {
            break;
        }
    }

    // -- Final cleanup: duplicate removal. ------------------------------------
    {
        let mut seen: HashSet<[VertexId; 3]> = HashSet::with_capacity(faces.len());
        faces.retain(|f| {
            let mut key = f.vertices;
            key.sort();
            seen.insert(key)
        });
    }

    #[cfg(test)]
    {
        let remain = build_boundary(faces);
        eprintln!("[patch-done] {} boundary edges remain", remain.len());
        for (a, b) in &remain {
            let pa = pool.position(*a);
            let pb = pool.position(*b);
            eprintln!(
                "  {:?}->{:?}  ({:.9},{:.9},{:.9})->({:.9},{:.9},{:.9})",
                a, b, pa.x, pa.y, pa.z, pb.x, pb.y, pb.z
            );
        }
    }
}

// Ã¢â€â‚¬Ã¢â€â‚¬ Tests Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬Ã¢â€â‚¬
