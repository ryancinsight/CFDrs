//! Shared planar arrangement helpers for CDT-backed CSG subroutines.
//!
//! These utilities consolidate the common "shatter constraints into PSLG
//! sub-edges" workflow used by:
//! - `clip::polygon2d::cdt` (native 2-D polygon Boolean)
//! - `corefine` (3-D face projection -> planar CDT subdivision)

use std::collections::HashSet;

use crate::application::delaunay::{Pslg, PslgVertexId};
use crate::domain::core::scalar::Real;

/// Canonical undirected edge key between two point slots.
pub(crate) type PlanarEdgeKey = (usize, usize);

/// Collect candidate `(t, slot)` pairs that lie on a segment interior (plus
/// explicit endpoints) using raw line projection.
///
/// This variant is used by the native 2-D polygon clipping path where polygon
/// loop endpoints are known explicitly and interior intersection points are
/// added from the merged point set.
pub(crate) fn collect_points_on_segment_interior(
    unique_pts: &[[Real; 2]],
    p1: [Real; 2],
    p2: [Real; 2],
    endpoint_slots: (usize, usize),
    t_eps: Real,
    dist_sq_tol: Real,
) -> Vec<(Real, usize)> {
    let (ri, rj) = endpoint_slots;
    let dx = p2[0] - p1[0];
    let dy = p2[1] - p1[1];
    let l2 = dx * dx + dy * dy;

    let mut out = vec![(0.0, ri), (1.0, rj)];
    if l2 < 1e-24 {
        return out;
    }

    for (slot, upt) in unique_pts.iter().enumerate() {
        if slot == ri || slot == rj {
            continue;
        }

        let t = ((upt[0] - p1[0]) * dx + (upt[1] - p1[1]) * dy) / l2;
        if t < t_eps || t > 1.0 - t_eps {
            continue;
        }

        let fx = p1[0] + t * dx;
        let fy = p1[1] + t * dy;
        let ex = upt[0] - fx;
        let ey = upt[1] - fy;
        let d2 = ex * ex + ey * ey;
        if d2 < dist_sq_tol {
            out.push((t, slot));
        }
    }

    out
}

/// Sort a `(t, slot)` point list, deduplicate slots, and emit shattered
/// consecutive sub-edges into `edges`.
pub(crate) fn insert_shattered_subedges(
    mut on_seg: Vec<(Real, usize)>,
    edges: &mut HashSet<PlanarEdgeKey>,
) {
    on_seg.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    on_seg.dedup_by_key(|x| x.1);

    for w in on_seg.windows(2) {
        let (a, b) = (w[0].1, w[1].1);
        if a == b {
            continue;
        }
        let key = if a < b { (a, b) } else { (b, a) };
        edges.insert(key);
    }
}

/// Build a PSLG from planar points and deduplicated undirected edge keys.
pub(crate) fn build_pslg_from_points_and_edges(
    points: &[[Real; 2]],
    edges: &HashSet<PlanarEdgeKey>,
) -> Pslg {
    let mut pslg = Pslg::with_capacity(points.len(), edges.len());
    for &p in points {
        pslg.add_vertex(p[0], p[1]);
    }
    for &(a, b) in edges {
        if a == b {
            continue;
        }
        let id_a = PslgVertexId::from_usize(a);
        let id_b = PslgVertexId::from_usize(b);
        pslg.add_segment(id_a, id_b);
    }
    pslg
}
