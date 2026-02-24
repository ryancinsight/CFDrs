//! 3-D mesh co-refinement via CDT (Constrained Delaunay Triangulation).
//!
//! ## Algorithm — Topological 2-D Projection + CDT
//!
//! ### Theorem — Projection Correctness (Watertightness Invariant)
//!
//! Two 3-D `VertexId`s that refer to the same welded pool position project to
//! the *same* 2-D coordinate under the dominant-axis-drop projection, and
//! therefore receive the same PSLG index.  Because the CDT uses Shewchuk exact
//! predicates (no epsilon fallbacks), any two adjacent face patches that share
//! a seam edge will produce exactly the same set of CDT edges along that seam,
//! eliminating T-junctions and achieving topological watertightness. ∎

use super::intersect::SnapSegment;
use crate::domain::core::index::VertexId;
use crate::domain::core::scalar::{Point3r, Real, Vector3r};
use crate::infrastructure::storage::face_store::FaceData;
use crate::infrastructure::storage::vertex_pool::VertexPool;

use crate::application::delaunay::pslg::vertex::PslgVertexId;
use crate::application::delaunay::{Cdt, Pslg};

/// Distance tolerance squared for edge Steiner projection.
const WELD_TOL_SQ: Real = 1e-8;
/// Exclude snap endpoints that fall exactly at a face corner.
const EDGE_EPS: Real = 1e-7;

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

    // ── Step 0: Choose projection axes ───────────────────────────────────────
    let (axis_u, axis_v) = dominant_normal_axes(face_n_unit);

    // ── Step 1: Edge Steiner vertices ─────────────────────────────────────────
    let mut edge_steiners: [Vec<(Real, VertexId)>; 3] = [Vec::new(), Vec::new(), Vec::new()];
    // seg_vids[i] = [vid_of_start, vid_of_end]; None if not on face boundary.
    let mut seg_vids: Vec<[Option<VertexId>; 2]> = vec![[None, None]; snap_segments.len()];
    let mut interior_vids: Vec<VertexId> = Vec::new();

    for (si, seg) in snap_segments.iter().enumerate() {
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

            // Corner vertex fallback
            if seg_vids[si][ep].is_none() {
                for ci in 0..3_usize {
                    if (p3d - face_pts[ci]).norm_squared() < WELD_TOL_SQ * 4.0 {
                        seg_vids[si][ep] = Some(face.vertices[ci]);
                        break;
                    }
                }
            }

            // Interior endpoint fallback
            if seg_vids[si][ep].is_none() && inside_triangle(p3d, a, b, c, face_n) {
                let vid = pool.insert_or_weld(p3d, face_n_unit);
                if !interior_vids.contains(&vid) {
                    interior_vids.push(vid);
                }
            }
        }
    }

    // ── Step 2: Interior crossing points ─────────────────────────────────────
    let n_sq = face_n.norm_squared();
    for i in 0..snap_segments.len() {
        for j in (i + 1)..snap_segments.len() {
            let s1 = &snap_segments[i];
            let s2 = &snap_segments[j];
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
            if !inside_triangle(crossing, a, b, c, face_n) {
                continue;
            }
            let vid = pool.insert_or_weld(crossing, face_n_unit);
            if !interior_vids.contains(&vid) {
                interior_vids.push(vid);
            }
        }
    }

    // ── Step 3: Early exit ────────────────────────────────────────────────────
    if edge_steiners.iter().all(|v| v.is_empty()) && interior_vids.is_empty() {
        return vec![*face];
    }

    // ── Step 4: Build ordered boundary polygon ────────────────────────────────
    for es in &mut edge_steiners {
        es.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    }
    let mut boundary_vids: Vec<VertexId> = Vec::with_capacity(9);
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

    // ── Step 5: Build PSLG from 2-D projections ───────────────────────────────
    // Map: VertexId (u32) → PslgVertexId index.
    let pool_len = pool.len();
    let mut vid_to_pslg: Vec<Option<PslgVertexId>> = vec![None; pool_len + 1];
    let mut pslg_to_vid: Vec<VertexId> = Vec::new();
    let mut pslg = Pslg::new();

    // Helper: register a VertexId into the PSLG if not already there.
    // Returns the PSLG vertex id.
    let register = |vid: VertexId,
                        pslg: &mut Pslg,
                        vid_to_pslg: &mut Vec<Option<PslgVertexId>>,
                        pslg_to_vid: &mut Vec<VertexId>,
                        pool: &VertexPool|
     -> PslgVertexId {
        let ui = vid.as_usize();
        if ui >= vid_to_pslg.len() {
            vid_to_pslg.resize(ui + 1, None);
        }
        if let Some(pid) = vid_to_pslg[ui] {
            return pid;
        }
        let pos3d = *pool.position(vid);
        let (u, v) = project_2d(pos3d, axis_u, axis_v);
        let pid = pslg.add_vertex(u, v);
        vid_to_pslg[ui] = Some(pid);
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
        let pa = vid_to_pslg[va.as_usize()].unwrap();
        let pb = vid_to_pslg[vb.as_usize()].unwrap();
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
            if let (Some(p0), Some(p1)) = (
                vid_to_pslg.get(v0.as_usize()).copied().flatten(),
                vid_to_pslg.get(v1.as_usize()).copied().flatten(),
            ) {
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

    for (a, b) in pslg_edges {
        let _ = pslg.add_segment(PslgVertexId::from_usize(a), PslgVertexId::from_usize(b));
    }

    // ── Step 6: Build CDT ────────────────────────────────────────────────────
    let cdt = Cdt::from_pslg(&pslg);
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

/// Returns `true` if `p` lies inside or on the boundary of triangle `(a,b,c)`.
fn inside_triangle(p: Point3r, a: Point3r, b: Point3r, c: Point3r, face_n: Vector3r) -> bool {
    const EPS: Real = 1e-9;
    let d0 = (b - a).cross(&(p - a)).dot(&face_n);
    let d1 = (c - b).cross(&(p - b)).dot(&face_n);
    let d2 = (a - c).cross(&(p - c)).dot(&face_n);
    d0 >= -EPS && d1 >= -EPS && d2 >= -EPS
}

/// Extract the region from FaceData (re-exported from the domain).
#[allow(dead_code)]
pub use crate::domain::core::index::RegionId as _RegionId;
