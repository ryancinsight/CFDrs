//! CDT-based polygon clipping backend (canonical production path).
//!
//! Triangulates both polygons and their intersections using constrained
//! Delaunay triangulation, then classifies triangles by centroid
//! winding-number tests.
//!
//! This is the production backend behind `clip2d::boolean_clip`.

use super::geometry::{ensure_ccw, point_in_polygon, seg_intersect};
use super::ClipOp;
use crate::application::csg::arrangement::planar::{
    build_pslg_from_points_and_edges, collect_points_on_segment_interior,
    insert_shattered_subedges, PlanarEdgeKey,
};
use crate::domain::core::scalar::Real;

/// CDT-based polygon clipping using constrained Delaunay triangulation.
///
/// # Algorithm
///
/// 1. Collect all vertices from both polygons and all edge-edge intersections.
/// 2. Build a PSLG with constraint edges from both polygons.
/// 3. Compute CDT.
/// 4. Classify each output triangle by testing its centroid against both
///    input polygons (winding number test).
/// 5. Emit triangles matching the requested Boolean operation.
pub fn cdt_clip(subject: &[[Real; 2]], clip: &[[Real; 2]], op: ClipOp) -> Vec<Vec<[Real; 2]>> {
    use crate::application::delaunay::Cdt;

    if subject.len() < 3 || clip.len() < 3 {
        return Vec::new();
    }

    let mut subj = subject.to_vec();
    let mut clp = clip.to_vec();
    ensure_ccw(&mut subj);
    ensure_ccw(&mut clp);

    // ── 1. Collect all vertices ──────────────────────────────────────────────
    let mut points: Vec<[Real; 2]> = Vec::new();
    let mut subj_indices: Vec<usize> = Vec::new();
    let mut clip_indices: Vec<usize> = Vec::new();

    for &p in &subj {
        subj_indices.push(points.len());
        points.push(p);
    }
    for &p in &clp {
        clip_indices.push(points.len());
        points.push(p);
    }

    // ── 2. Compute edge-edge intersections ───────────────────────────────────
    let sn = subj.len();
    let cn = clp.len();
    for i in 0..sn {
        let j = (i + 1) % sn;
        for k in 0..cn {
            let l = (k + 1) % cn;
            if let Some((t, s)) = seg_intersect(subj[i], subj[j], clp[k], clp[l]) {
                if t > 1e-10 && t < 1.0 - 1e-10 && s > 1e-10 && s < 1.0 - 1e-10 {
                    let px = subj[i][0] + t * (subj[j][0] - subj[i][0]);
                    let py = subj[i][1] + t * (subj[j][1] - subj[i][1]);
                    points.push([px, py]);
                }
            }
        }
    }

    // ── 3. Deduplicate ───────────────────────────────────────────────────────
    const WELD_TOL: Real = 1e-8;
    let mut canonical: Vec<usize> = (0..points.len()).collect();
    for i in 0..points.len() {
        for j in 0..i {
            let dx = points[i][0] - points[j][0];
            let dy = points[i][1] - points[j][1];
            if dx * dx + dy * dy < WELD_TOL * WELD_TOL {
                canonical[i] = canonical[j];
                break;
            }
        }
    }

    // Build unique point list.
    let mut unique: Vec<[Real; 2]> = Vec::new();
    let mut remap: Vec<usize> = vec![0; points.len()];
    let mut seen = vec![false; points.len()];
    for i in 0..points.len() {
        let c = canonical[i];
        if !seen[c] {
            seen[c] = true;
            remap[c] = unique.len();
            unique.push(points[c]);
        }
        remap[i] = remap[c];
    }

    if unique.len() < 3 {
        return Vec::new();
    }

    // ── 4. Build PSLG ───────────────────────────────────────────────────────
    let mut pslg_edges = std::collections::HashSet::<PlanarEdgeKey>::new();
    // Subject edges (shattered at intersection points).
    add_shattered_edges(&subj, &subj_indices, &remap, &unique, &mut pslg_edges);
    // Clip edges (shattered at intersection points).
    add_shattered_edges(&clp, &clip_indices, &remap, &unique, &mut pslg_edges);

    let pslg = build_pslg_from_points_and_edges(&unique, &pslg_edges);

    // ── 5. Build CDT ─────────────────────────────────────────────────────────
    let cdt = Cdt::from_pslg(&pslg);
    let dt = cdt.triangulation();

    // ── 6. Classify and emit triangles ───────────────────────────────────────
    let verts = dt.vertices();
    let mut result: Vec<Vec<[Real; 2]>> = Vec::new();

    for (_, tri) in dt.interior_triangles() {
        let [v0, v1, v2] = tri.vertices;
        let p0 = [verts[v0.idx()].x, verts[v0.idx()].y];
        let p1 = [verts[v1.idx()].x, verts[v1.idx()].y];
        let p2 = [verts[v2.idx()].x, verts[v2.idx()].y];

        let mx = (p0[0] + p1[0] + p2[0]) / 3.0;
        let my = (p0[1] + p1[1] + p2[1]) / 3.0;

        let in_subj = point_in_polygon(mx, my, &subj);
        let in_clip = point_in_polygon(mx, my, &clp);

        let keep = match op {
            ClipOp::Intersection => in_subj && in_clip,
            ClipOp::Union => in_subj || in_clip,
            ClipOp::Difference => in_subj && !in_clip,
        };

        if keep {
            result.push(vec![p0, p1, p2]);
        }
    }

    result
}

fn add_shattered_edges(
    poly: &[[Real; 2]],
    indices: &[usize],
    remap: &[usize],
    unique: &[[Real; 2]],
    pslg_edges: &mut std::collections::HashSet<PlanarEdgeKey>,
) {
    let n = poly.len();
    for i in 0..n {
        let j = (i + 1) % n;
        let ri = remap[indices[i]];
        let rj = remap[indices[j]];
        if ri == rj {
            continue;
        }
        let on_edge =
            collect_points_on_segment_interior(unique, poly[i], poly[j], (ri, rj), 1e-8, 1e-12);
        insert_shattered_subedges(on_edge, pslg_edges);
    }
}
