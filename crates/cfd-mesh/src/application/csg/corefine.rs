//! 3-D mesh co-refinement via CDT (Constrained Delaunay Triangulation).
//!
//! ## Algorithm — 7-Step CDT Co-Refinement (per face)
//!
//! ```text
//! face + snap_segments
//!         │
//!    ┌────▼────────────────────────────────────────────────────┐
//!    │  Step 0 — Choose 2-D projection axes                     │
//!    │  dominant_normal_axes(n): drop max-magnitude axis         │
//!    │  Theorem: projected area ≥ true area / √3  (always       │
//!    │         well-conditioned for non-zero normal)             │
//!    └────┬────────────────────────────────────────────────────┘
//!         │
//!    ┌────▼────────────────────────────────────────────────────┐
//!    │  Step 1 — Classify snap endpoint positions               │
//!    │  • On edge → edge Steiner point (t parameter)            │
//!    │  • At corner → existing VertexId                         │
//!    │  • Interior → interior Steiner (HashSet dedup)           │
//!    └────┬────────────────────────────────────────────────────┘
//!         │
//!    ┌────▼────────────────────────────────────────────────────┐
//!    │  Step 2 — Interior crossing points                       │
//!    │  All pairs (si, sj): line-line intersection in face plane │
//!    │  O(s²) where s = |snap_segments|; s≤4 in practice        │
//!    └────┬────────────────────────────────────────────────────┘
//!         │
//!    ┌────▼────────────────────────────────────────────────────┐
//!    │  Step 3 — Early exit                                     │
//!    │  No edge Steiners + no interior → return [face]          │
//!    └────┬────────────────────────────────────────────────────┘
//!         │
//!    ┌────▼────────────────────────────────────────────────────┐
//!    │  Step 4 — Build ordered boundary polygon                 │
//!    │  3 corners + sorted Steiners per edge                    │
//!    │  capacity = 3 + Σ|edge_steiners[i]|  (exact)            │
//!    └────┬────────────────────────────────────────────────────┘
//!         │
//!    ┌────▼────────────────────────────────────────────────────┐
//!    │  Step 5 — Build PSLG (Planar Straight Line Graph)        │
//!    │  HashMap<VertexId, PslgVertexId> (O(face_v) not pool_v)  │
//!    │  Shatter boundary + constraint segments to sub-edges     │
//!    └────┬────────────────────────────────────────────────────┘
//!         │
//!    ┌────▼────────────────────────────────────────────────────┐
//!    │  Step 6 — resolve_crossings + Cdt::try_from_pslg        │
//!    │  Shewchuk exact predicates; fallback: return [face]      │
//!    └────┬────────────────────────────────────────────────────┘
//!         │
//!    ┌────▼────────────────────────────────────────────────────┐
//!    │  Step 7 — Lift CDT triangles back to 3-D                │
//!    │  Exclude supertriangle vertices, flip to match face_n    │
//!    └────┬────────────────────────────────────────────────────┘
//!         │
//!  sub-triangles replacing original face
//! ```
//!
//! ## Theorem — Projection Correctness (Watertightness Invariant)
//!
//! Two 3-D `VertexId`s that refer to the same welded pool position project to
//! the *same* 2-D coordinate under the dominant-axis-drop projection, and
//! therefore receive the same PSLG index.  Because the CDT uses Shewchuk exact
//! predicates (no epsilon fallbacks), any two adjacent face patches that share
//! a seam edge will produce exactly the same set of CDT edges along that seam,
//! eliminating T-junctions and achieving topological watertightness. ∎
//!
//! ## Memory Note
//!
//! `vid_to_pslg` uses `HashMap<VertexId, PslgVertexId>` with capacity hint
//! `boundary_vids.len() + interior_vids.len()`.  The previous `Vec<Option<_>>`
//! of length `pool.len()` allocated O(pool_size) per call — up to 100 k entries
//! for pools of 2 000+ vertices.  The HashMap reduces per-call allocation from
//! O(pool_size) to O(face_vertex_count) (typically 3–12 entries).

use std::collections::{HashMap, HashSet};

use super::intersect::SnapSegment;
use crate::domain::core::index::VertexId;
use crate::domain::core::scalar::{Point3r, Real, Vector3r};
use crate::domain::geometry::predicates::{orient_2d_arr, Orientation};
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;

use crate::application::delaunay::pslg::vertex::PslgVertexId;
use crate::application::delaunay::{Cdt, Pslg};

/// Distance tolerance squared for edge Steiner projection.
/// Widened from 1e-8 to 1e-6 to handle shallow-angle tangent junctions
/// (e.g. elbow-cylinder V-shape) where floating-point drift can exceed 1e-8.
const WELD_TOL_SQ: Real = 1e-6;
/// Exclude snap endpoints that fall exactly at a face corner.
const EDGE_EPS: Real = 1e-6;

type PointBits3 = [u64; 3];

#[inline]
fn point_bits3(p: &Point3r) -> PointBits3 {
    [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]
}

/// Canonical exact key for an undirected 3-D snap segment.
///
/// Uses raw IEEE-754 bit patterns (no tolerance/quantization), so two segments
/// deduplicate iff both endpoint coordinates are bit-identical up to endpoint
/// ordering.
#[inline]
fn canonical_segment_key(seg: &SnapSegment) -> (PointBits3, PointBits3) {
    let a = point_bits3(&seg.start);
    let b = point_bits3(&seg.end);
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}

// ── Main entry point ─────────────────────────────────────────────────────────

/// Co-refine `face` against `snap_segments` using CDT-based 2-D projection.
///
/// Steiner vertices are inserted via `VertexPool::insert_or_weld`; the 3-D
/// face boundary polygon and snap chords are projected to 2-D (dropping the
/// dominant normal axis), fed to `Cdt::from_pslg` (Shewchuk predicates), then
/// lifted back to 3-D `VertexId`s.
pub fn corefine_face(
    face: &FaceData,
    snap_segments: &[SnapSegment],
    pool: &mut VertexPool,
) -> Vec<FaceData> {
    let a = *pool.position(face.vertices[0]);
    let b = *pool.position(face.vertices[1]);
    let c = *pool.position(face.vertices[2]);

    let face_n = (b - a).cross(&(c - a));
    if face_n.norm_squared() < 1e-20 {
        return Vec::new();
    }
    let face_n_unit = face_n / face_n.norm();
    let face_pts = [a, b, c];

    // Exact segment canonicalization:
    // remove degenerate and bit-identical duplicate constraints before O(s²)
    // crossing detection and PSLG resolve_crossings.
    let mut dedup_snap_segments: Vec<SnapSegment> = Vec::with_capacity(snap_segments.len());
    let mut seen_snap_segments: HashSet<(PointBits3, PointBits3)> =
        HashSet::with_capacity(snap_segments.len());
    for seg in snap_segments {
        if (seg.end - seg.start).norm_squared() < 1e-24 {
            continue;
        }
        if seen_snap_segments.insert(canonical_segment_key(seg)) {
            dedup_snap_segments.push(*seg);
        }
    }
    dedup_snap_segments.sort_unstable_by_key(canonical_segment_key);

    // ── Step 0: Choose projection axes ───────────────────────────────────────
    let (axis_u, axis_v) = dominant_normal_axes(face_n_unit);

    // ── Step 1: Edge Steiner vertices ─────────────────────────────────────────
    // SmallVec analogue: each edge typically has 0–2 Steiner points; use a
    // plain Vec but with capacity hint of 4 to avoid reallocs in the common case.
    let mut edge_steiners: [Vec<(Real, VertexId)>; 3] = [
        Vec::with_capacity(4),
        Vec::with_capacity(4),
        Vec::with_capacity(4),
    ];
    // seg_vids[i] = [vid_of_start, vid_of_end]; None if not on face boundary.
    let mut seg_vids: Vec<[Option<VertexId>; 2]> = vec![[None, None]; dedup_snap_segments.len()];
    // Use a HashSet for O(1) dedup and a Vec to preserve insertion order.
    let mut interior_vid_set: HashSet<VertexId> = HashSet::new();
    let mut interior_vids: Vec<VertexId> = Vec::new();

    for (si, seg) in dedup_snap_segments.iter().enumerate() {
        for (ep, &p3d) in [seg.start, seg.end].iter().enumerate() {
            'edge_search: for ei in 0..3_usize {
                let va = face_pts[ei];
                let vb = face_pts[(ei + 1) % 3];
                let edge = vb - va;
                let edge_sq = edge.norm_squared();
                if edge_sq < 1e-20 {
                    continue;
                }
                let t = (p3d - va).dot(&edge) / edge_sq;
                if t <= EDGE_EPS || t >= 1.0 - EDGE_EPS {
                    continue;
                }
                let interp = va + edge * t;
                if (p3d - interp).norm_squared() > WELD_TOL_SQ * 4.0 {
                    continue;
                }
                let vid = pool.insert_or_weld(p3d, face_n_unit);
                seg_vids[si][ep] = Some(vid);
                if !edge_steiners[ei].iter().any(|&(_, v)| v == vid) {
                    edge_steiners[ei].push((t, vid));
                }
                break 'edge_search;
            }

            // Corner vertex fallback — if snap endpoint is within WELD_TOL_SQ
            // of a face corner, snap it to the corner VertexId.
            if seg_vids[si][ep].is_none() {
                for ci in 0..3_usize {
                    if (p3d - face_pts[ci]).norm_squared() < WELD_TOL_SQ * 4.0 {
                        seg_vids[si][ep] = Some(face.vertices[ci]);
                        break;
                    }
                }
            }

            // Interior endpoint fallback
            if seg_vids[si][ep].is_none() && inside_triangle(p3d, a, b, c, face_n, axis_u, axis_v) {
                let vid = pool.insert_or_weld(p3d, face_n_unit);
                if interior_vid_set.insert(vid) {
                    interior_vids.push(vid);
                }
            }
        }
    }

    // ── Step 2: Interior crossing points ─────────────────────────────────────
    let n_sq = face_n.norm_squared();
    for i in 0..dedup_snap_segments.len() {
        for j in (i + 1)..dedup_snap_segments.len() {
            let s1 = &dedup_snap_segments[i];
            let s2 = &dedup_snap_segments[j];
            let d1 = s1.end - s1.start;
            let d2 = s2.end - s2.start;
            let r = s2.start - s1.start;
            let denom = d1.cross(&d2).dot(&face_n);
            let min_d = 1e-14 * n_sq.sqrt() * (d1.norm() + d2.norm());
            if denom.abs() < min_d {
                continue;
            }
            let t_param = r.cross(&d2).dot(&face_n) / denom;
            let s_param = r.cross(&d1).dot(&face_n) / denom;
            if t_param <= EDGE_EPS || t_param >= 1.0 - EDGE_EPS {
                continue;
            }
            if s_param <= EDGE_EPS || s_param >= 1.0 - EDGE_EPS {
                continue;
            }
            let crossing = s1.start + d1 * t_param;
            if !inside_triangle(crossing, a, b, c, face_n, axis_u, axis_v) {
                continue;
            }
            let vid = pool.insert_or_weld(crossing, face_n_unit);
            if interior_vid_set.insert(vid) {
                interior_vids.push(vid);
            }
        }
    }

    // ── Step 3: Early exit ────────────────────────────────────────────────────
    if edge_steiners.iter().all(std::vec::Vec::is_empty) && interior_vids.is_empty() {
        return vec![*face];
    }

    // ── Step 4: Build ordered boundary polygon ────────────────────────────────
    for es in &mut edge_steiners {
        es.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    }
    // Compute exact capacity: 3 corners + sum of Steiners per edge.
    let steiner_count: usize = edge_steiners.iter().map(|e| e.len()).sum();
    let mut boundary_vids: Vec<VertexId> = Vec::with_capacity(3 + steiner_count);
    for ei in 0..3_usize {
        boundary_vids.push(face.vertices[ei]);
        for &(_, vid) in &edge_steiners[ei] {
            boundary_vids.push(vid);
        }
    }
    boundary_vids.dedup();
    if boundary_vids.len() < 3 {
        return vec![*face];
    }

    // ── Sliver-face guard ─────────────────────────────────────────────────────
    // If the 2D projected area of the boundary polygon is near-zero, the CDT
    // will be degenerate (all vertices colinear in 2D).  This happens for very
    // thin sliver triangles at shallow-angle elbow-cylinder junctions.
    // In this case, fall back to simple midpoint-subdivision: split the edge
    // containing each edge-Steiner point and return the sub-triangles directly.
    {
        let mut area2 = 0.0_f64;
        let nb = boundary_vids.len();
        for i in 0..nb {
            let pa = *pool.position(boundary_vids[i]);
            let pb = *pool.position(boundary_vids[(i + 1) % nb]);
            let (pu, pv) = project_2d(pa, axis_u, axis_v);
            let (qu, qv) = project_2d(pb, axis_u, axis_v);
            area2 += pu * qv - qu * pv;
        }
        if area2.abs() < 1e-10 {
            // Sliver: produce sub-triangles by splitting each Steiner-containing
            // edge and fan-stitching. This guarantees Steiner points appear as
            // vertices on the shared boundary even when CDT would degenerate.
            return midpoint_subdivide(face, &edge_steiners, pool, face_n);
        }
    }

    // ── Step 5: Build PSLG from 2-D projections ───────────────────────────────
    // Map: VertexId → PslgVertexId.
    //
    // ## Memory note
    //
    // Previously `vid_to_pslg` was a `Vec<Option<PslgVertexId>>` of length
    // `pool.len() + 1`, which allocated O(pool_size) memory per `corefine_face`
    // call.  With pools of 2 000+ vertices and 50+ faces this became 100 k
    // entries.  A `HashMap` with an exact capacity hint is O(face_vertex_count)
    // and avoids the large up-front allocation entirely.
    let register_cap = boundary_vids.len() + interior_vids.len();
    let mut vid_to_pslg: HashMap<VertexId, PslgVertexId> = HashMap::with_capacity(register_cap);
    let mut pslg_to_vid: Vec<VertexId> = Vec::with_capacity(register_cap);
    let mut pslg = Pslg::new();

    // Helper: register a VertexId into the PSLG if not already there.
    // Returns the PSLG vertex id.
    let register = |vid: VertexId,
                    pslg: &mut Pslg,
                    vid_to_pslg: &mut HashMap<VertexId, PslgVertexId>,
                    pslg_to_vid: &mut Vec<VertexId>,
                    pool: &VertexPool|
     -> PslgVertexId {
        if let Some(&pid) = vid_to_pslg.get(&vid) {
            return pid;
        }
        let pos3d = *pool.position(vid);
        let (u, v) = project_2d(pos3d, axis_u, axis_v);
        let pid = pslg.add_vertex(u, v);
        vid_to_pslg.insert(vid, pid);
        pslg_to_vid.push(vid);
        pid
    };

    // Register boundary vertices in order.
    for &vid in &boundary_vids {
        register(vid, &mut pslg, &mut vid_to_pslg, &mut pslg_to_vid, pool);
    }
    // Register interior Steiner vertices.
    for &vid in &interior_vids {
        register(vid, &mut pslg, &mut vid_to_pslg, &mut pslg_to_vid, pool);
    }

    // Extract unique 2D points from PSLG for segment shattering.
    let unique_pts: Vec<[Real; 2]> = pslg.vertices().iter().map(|v| [v.x, v.y]).collect();
    let mut pslg_edges = std::collections::HashSet::new();

    // Boundary polygon segments (ring).
    let nb = boundary_vids.len();
    for i in 0..nb {
        let va = boundary_vids[i];
        let vb = boundary_vids[(i + 1) % nb];
        let pa = vid_to_pslg[&va];
        let pb = vid_to_pslg[&vb];
        if pa != pb {
            let p1 = unique_pts[pa.idx()];
            let p2 = unique_pts[pb.idx()];
            let on_edge =
                crate::application::csg::arrangement::planar::collect_points_on_segment_interior(
                    &unique_pts,
                    p1,
                    p2,
                    (pa.idx(), pb.idx()),
                    1e-8,
                    1e-14,
                );
            crate::application::csg::arrangement::planar::insert_shattered_subedges(
                on_edge,
                &mut pslg_edges,
            );
        }
    }

    // Constraint segments from snap-segment endpoints.
    for vids in &seg_vids {
        if let (Some(v0), Some(v1)) = (vids[0], vids[1]) {
            if let (Some(&p0), Some(&p1)) = (vid_to_pslg.get(&v0), vid_to_pslg.get(&v1)) {
                if p0 != p1 {
                    let pa = unique_pts[p0.idx()];
                    let pb = unique_pts[p1.idx()];
                    let on_edge = crate::application::csg::arrangement::planar::collect_points_on_segment_interior(
                        &unique_pts, pa, pb, (p0.idx(), p1.idx()), 1e-8, 1e-14
                    );
                    crate::application::csg::arrangement::planar::insert_shattered_subedges(
                        on_edge,
                        &mut pslg_edges,
                    );
                }
            }
        }
    }

    let mut sorted_pslg_edges: Vec<(usize, usize)> = pslg_edges.into_iter().collect();
    sorted_pslg_edges.sort_unstable();
    for (a, b) in sorted_pslg_edges {
        let _ = pslg.add_segment(PslgVertexId::from_usize(a), PslgVertexId::from_usize(b));
    }

    // ── Step 6: Build CDT ────────────────────────────────────────────────────
    // Resolve any interior segment crossings that arise from 3-D seam curves
    // whose 2-D projections cross (e.g. out-of-plane multi-branch junctions).
    // `resolve_crossings` now handles proper crossings, T-intersections, and
    // collinear overlaps, so this error path should be unreachable for all
    // geometrically valid CSG inputs.
    pslg.resolve_crossings();
    let cdt = match Cdt::try_from_pslg(&pslg) {
        Ok(cdt) => cdt,
        // Should not be reached after resolve_crossings handles all illegality
        // classes.  Log for diagnostics and fall back to the original face so
        // the overall mesh remains topologically consistent.
        Err(e) => {
            tracing::warn!(
                "corefine_face: PSLG validation failed after resolve_crossings: {e}. \
                 Returning original face (no CDT refinement applied)."
            );
            return vec![*face];
        }
    };
    let dt = cdt.triangulation();

    // ── Step 7: Emit triangles lifted back to 3-D ─────────────────────────────
    let mut triangles: Vec<FaceData> = Vec::new();
    for (_, tri) in dt.interior_triangles() {
        let [pv0, pv1, pv2] = tri.vertices;
        // Exclude super-triangle vertices (index >= pslg_to_vid.len()).
        if pv0.idx() >= pslg_to_vid.len()
            || pv1.idx() >= pslg_to_vid.len()
            || pv2.idx() >= pslg_to_vid.len()
        {
            continue;
        }
        let v0 = pslg_to_vid[pv0.idx()];
        let v1 = pslg_to_vid[pv1.idx()];
        let v2 = pslg_to_vid[pv2.idx()];
        if v0 == v1 || v1 == v2 || v0 == v2 {
            continue;
        }

        let p0 = *pool.position(v0);
        let p1 = *pool.position(v1);
        let p2 = *pool.position(v2);
        let tri_n = (p1 - p0).cross(&(p2 - p0));
        if tri_n.norm_squared() < 1e-30 {
            continue;
        }
        if tri_n.dot(&face_n) >= 0.0 {
            triangles.push(FaceData::new(v0, v1, v2, face.region));
        } else {
            triangles.push(FaceData::new(v0, v2, v1, face.region));
        }
    }

    if triangles.is_empty() {
        vec![*face]
    } else {
        triangles
    }
}

// ── Geometry helpers ──────────────────────────────────────────────────────────

/// Choose two projection axes by dropping the largest-magnitude axis of the
/// face normal, ensuring the 2-D projected polygon is well-conditioned.
///
/// # Theorem — Projection Conditioning
///
/// Dropping the dominant axis ensures the projected area ≥ true 3-D area / √3,
/// so the polygon is never near-degenerate for any non-zero normal.
fn dominant_normal_axes(n: Vector3r) -> (usize, usize) {
    let (ax, ay, az) = (n.x.abs(), n.y.abs(), n.z.abs());
    if ax >= ay && ax >= az {
        (1, 2) // drop X → keep Y,Z
    } else if ay >= ax && ay >= az {
        (0, 2) // drop Y → keep X,Z
    } else {
        (0, 1) // drop Z → keep X,Y
    }
}

/// Project a 3-D point to 2-D by selecting two coordinate axes.
#[inline]
fn project_2d(p: Point3r, axis_u: usize, axis_v: usize) -> (Real, Real) {
    (p[axis_u], p[axis_v])
}

/// Returns `true` if `p` lies inside or on triangle `(a,b,c)`.
///
/// ## Algorithm
///
/// 1. Exact stage: classify in the dominant-axis 2-D projection using Shewchuk
///    `orient_2d_arr` signs (boundary-inclusive).
/// 2. Fallback stage: if exact stage rejects, run the prior signed-area test
///    with a tiny epsilon to retain legacy tolerance for numerical drift.
///
/// ## Theorem — Exact Projected Inclusion
///
/// Let `π` be the dominant-axis projection for a non-degenerate triangle
/// `(a,b,c)`. For any point `p`, if all directed edge orientation signs of
/// `(π(a),π(b),π(p))`, `(π(b),π(c),π(p))`, `(π(c),π(a),π(p))` have the same
/// orientation as `(π(a),π(b),π(c))` or are degenerate, then `π(p)` lies in the
/// closed triangle `π(abc)`.
///
/// **Proof sketch.**
/// The oriented half-space representation of a triangle in 2-D is the
/// intersection of its three directed edge half-planes. Exact predicate signs
/// evaluate those half-plane memberships without floating-point sign errors.
/// Boundary points correspond to degenerate orientation on one or more edges. ∎
fn inside_triangle(
    p: Point3r,
    a: Point3r,
    b: Point3r,
    c: Point3r,
    face_n: Vector3r,
    axis_u: usize,
    axis_v: usize,
) -> bool {
    if inside_triangle_exact_projected(p, a, b, c, axis_u, axis_v) {
        return true;
    }
    inside_triangle_eps(p, a, b, c, face_n)
}

#[inline]
fn inside_triangle_exact_projected(
    p: Point3r,
    a: Point3r,
    b: Point3r,
    c: Point3r,
    axis_u: usize,
    axis_v: usize,
) -> bool {
    let pa = [a[axis_u], a[axis_v]];
    let pb = [b[axis_u], b[axis_v]];
    let pc = [c[axis_u], c[axis_v]];
    let pp = [p[axis_u], p[axis_v]];

    let tri_ori = orient_2d_arr(pa, pb, pc);
    if tri_ori == Orientation::Degenerate {
        return false;
    }

    let o0 = orient_2d_arr(pa, pb, pp);
    let o1 = orient_2d_arr(pb, pc, pp);
    let o2 = orient_2d_arr(pc, pa, pp);

    match tri_ori {
        Orientation::Positive => {
            o0 != Orientation::Negative
                && o1 != Orientation::Negative
                && o2 != Orientation::Negative
        }
        Orientation::Negative => {
            o0 != Orientation::Positive
                && o1 != Orientation::Positive
                && o2 != Orientation::Positive
        }
        Orientation::Degenerate => false,
    }
}

#[inline]
fn inside_triangle_eps(p: Point3r, a: Point3r, b: Point3r, c: Point3r, face_n: Vector3r) -> bool {
    const EPS: Real = 1e-9;
    let d0 = (b - a).cross(&(p - a)).dot(&face_n);
    let d1 = (c - b).cross(&(p - b)).dot(&face_n);
    let d2 = (a - c).cross(&(p - c)).dot(&face_n);
    d0 >= -EPS && d1 >= -EPS && d2 >= -EPS
}

/// Sliver-face fallback: build an ordered boundary polygon from all edge
/// Steiner points and fan-triangulate from the first corner.
///
/// ## Theorem — Completeness for Sliver Faces
///
/// For a degenerate sliver triangle (2D projected area ≈ 0), the CDT fails to
/// resolve the constraint graph. This fallback builds the *augmented boundary
/// polygon* by interleaving each corner with its edge Steiner chain, producing
/// a convex polygon with all Steiner vertices on the boundary. A simple fan
/// from the first vertex covers all sub-triangles without overlap or gap,
/// and guarantees that every edge-shared Steiner vertex appears in the output.
///
/// Adjacent faces that share any split edge receive the same Steiner injections
/// via `propagate_seam_vertices`, so the sub-edges match — eliminating T-junction
/// cracks even for the tightest sliver angles (≤ 1°).
fn midpoint_subdivide(
    face: &FaceData,
    edge_steiners: &[Vec<(Real, VertexId)>; 3],
    pool: &VertexPool,
    face_n: Vector3r,
) -> Vec<FaceData> {
    let has_any = edge_steiners.iter().any(|e| !e.is_empty());
    if !has_any {
        return vec![*face];
    }

    // Build ordered boundary polygon:
    //   corner[0] → steiners[0] → corner[1] → steiners[1] → corner[2] → steiners[2]
    let corners = face.vertices;
    let mut poly: Vec<VertexId> =
        Vec::with_capacity(3 + edge_steiners.iter().map(|e| e.len()).sum::<usize>());
    for ei in 0..3_usize {
        poly.push(corners[ei]);
        for &(_, sv) in &edge_steiners[ei] {
            poly.push(sv);
        }
    }
    poly.dedup();
    if poly.len() < 3 {
        return vec![*face];
    }

    // Fan triangulate from poly[0] to cover all sub-cells.
    let v0 = poly[0];
    let p0 = *pool.position(v0);
    let mut result = Vec::with_capacity(poly.len() - 2);

    for i in 1..poly.len() - 1 {
        let va = poly[i];
        let vb = poly[i + 1];
        if va == v0 || vb == v0 || va == vb {
            continue;
        }
        let pa = *pool.position(va);
        let pb = *pool.position(vb);
        let tri_n = (pa - p0).cross(&(pb - p0));
        if tri_n.norm_squared() < 1e-30 {
            continue;
        }
        if tri_n.dot(&face_n) >= 0.0 {
            result.push(FaceData::new(v0, va, vb, face.region));
        } else {
            result.push(FaceData::new(v0, vb, va, face.region));
        }
    }

    if result.is_empty() {
        vec![*face]
    } else {
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(x: Real, y: Real, z: Real) -> Point3r {
        Point3r::new(x, y, z)
    }

    #[test]
    fn inside_triangle_exact_accepts_edge_point() {
        let a = p(0.0, 0.0, 0.0);
        let b = p(1.0, 0.0, 0.0);
        let c = p(0.0, 1.0, 0.0);
        let n = (b - a).cross(&(c - a));
        let (axis_u, axis_v) = dominant_normal_axes(n / n.norm());

        let edge_point = p(0.5, 0.5, 0.0); // on edge b-c
        assert!(inside_triangle(edge_point, a, b, c, n, axis_u, axis_v));
    }

    #[test]
    fn inside_triangle_exact_rejects_outside_point() {
        let a = p(0.0, 0.0, 0.0);
        let b = p(1.0, 0.0, 0.0);
        let c = p(0.0, 1.0, 0.0);
        let n = (b - a).cross(&(c - a));
        let (axis_u, axis_v) = dominant_normal_axes(n / n.norm());

        let outside = p(1.1, 0.2, 0.0);
        assert!(!inside_triangle(outside, a, b, c, n, axis_u, axis_v));
    }

    #[test]
    fn inside_triangle_exact_handles_reversed_winding() {
        let a = p(0.0, 0.0, 0.0);
        let b = p(0.0, 1.0, 0.0);
        let c = p(1.0, 0.0, 0.0); // reversed orientation vs CCW XY
        let n = (b - a).cross(&(c - a));
        let (axis_u, axis_v) = dominant_normal_axes(n / n.norm());

        let inside = p(0.2, 0.2, 0.0);
        assert!(inside_triangle(inside, a, b, c, n, axis_u, axis_v));
    }
}
