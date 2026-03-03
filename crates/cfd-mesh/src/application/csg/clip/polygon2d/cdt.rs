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
use std::collections::{HashMap, HashSet};

const WELD_TOL: Real = 1e-8;
const INTERIOR_TOL: Real = 1e-10;

#[derive(Clone, Copy, Debug)]
struct EdgeAabb2d {
    a: [Real; 2],
    b: [Real; 2],
    min_x: Real,
    min_y: Real,
    max_x: Real,
    max_y: Real,
}

impl EdgeAabb2d {
    fn new(a: [Real; 2], b: [Real; 2]) -> Self {
        Self {
            a,
            b,
            min_x: a[0].min(b[0]),
            min_y: a[1].min(b[1]),
            max_x: a[0].max(b[0]),
            max_y: a[1].max(b[1]),
        }
    }

    fn from_polygon(poly: &[[Real; 2]]) -> Vec<Self> {
        let n = poly.len();
        let mut edges = Vec::with_capacity(n);
        for i in 0..n {
            let j = (i + 1) % n;
            let a = poly[i];
            let b = poly[j];
            if a != b {
                edges.push(Self::new(a, b));
            }
        }
        edges
    }

    #[inline]
    fn aabb_overlaps(self, other: Self) -> bool {
        self.min_x <= other.max_x
            && other.min_x <= self.max_x
            && self.min_y <= other.max_y
            && other.min_y <= self.max_y
    }
}

/// Uniform-grid spatial welding for 2-D points.
///
/// # Theorem — 3x3 neighborhood suffices
///
/// Let the grid cell width be `h = tol`. For any two points `p` and `q` with
/// `||p-q||_2 < tol`, their cell coordinates differ by at most 1 along each
/// axis, so `q` lies in one of the 9 neighboring cells of `p`.
///
/// **Proof sketch**
///
/// For each axis, `|p_i-q_i| <= ||p-q||_2 < h`. Dividing by `h` and taking
/// floor indices implies cell index difference in `{-1,0,1}`. Cartesian product
/// over two axes gives exactly the 3x3 neighborhood. ∎
struct SpatialHashWeld2d {
    inv_cell: Real,
    tol_sq: Real,
    bins: HashMap<(i64, i64), Vec<usize>>,
}

impl SpatialHashWeld2d {
    fn new(tol: Real) -> Self {
        Self {
            inv_cell: 1.0 / tol,
            tol_sq: tol * tol,
            bins: HashMap::new(),
        }
    }

    #[inline]
    fn cell_of(&self, p: [Real; 2]) -> (i64, i64) {
        (
            (p[0] * self.inv_cell).floor() as i64,
            (p[1] * self.inv_cell).floor() as i64,
        )
    }

    fn insert_or_weld(&mut self, p: [Real; 2], unique: &mut Vec<[Real; 2]>) -> usize {
        let (cx, cy) = self.cell_of(p);

        for dx in -1_i64..=1 {
            for dy in -1_i64..=1 {
                if let Some(candidates) = self.bins.get(&(cx + dx, cy + dy)) {
                    for &idx in candidates {
                        let q = unique[idx];
                        let ex = p[0] - q[0];
                        let ey = p[1] - q[1];
                        if ex * ex + ey * ey <= self.tol_sq {
                            return idx;
                        }
                    }
                }
            }
        }

        let idx = unique.len();
        unique.push(p);
        self.bins.entry((cx, cy)).or_default().push(idx);
        idx
    }
}

/// Sweep index for clip edges sorted by `min_x`.
struct ClipEdgeSweepIndex2d {
    by_min_x: Vec<usize>,
}

impl ClipEdgeSweepIndex2d {
    fn build(clip_edges: &[EdgeAabb2d]) -> Self {
        let mut by_min_x: Vec<usize> = (0..clip_edges.len()).collect();
        by_min_x.sort_unstable_by(|&i, &j| clip_edges[i].min_x.total_cmp(&clip_edges[j].min_x));
        Self { by_min_x }
    }

    /// Collect all clip-edge indices whose AABBs overlap `subject_edge`.
    ///
    /// # Theorem — completeness
    ///
    /// Let `L = {c | c.min_x <= s.max_x}`. Any clip edge whose AABB overlaps
    /// `s` must be in `L` (x-overlap necessity), so scanning only `L` cannot
    /// miss true overlaps. Filtering by `c.max_x >= s.min_x` and y-overlap then
    /// yields exactly the AABB-overlap set.
    fn query_overlaps(
        &self,
        subject_edge: EdgeAabb2d,
        clip_edges: &[EdgeAabb2d],
        out: &mut Vec<usize>,
    ) {
        out.clear();
        let limit = self
            .by_min_x
            .partition_point(|&idx| clip_edges[idx].min_x <= subject_edge.max_x);
        for &cj in &self.by_min_x[..limit] {
            let ce = clip_edges[cj];
            if ce.max_x < subject_edge.min_x {
                continue;
            }
            if subject_edge.aabb_overlaps(ce) {
                out.push(cj);
            }
        }
        // Preserve original clip-edge order for deterministic welding behavior.
        out.sort_unstable();
    }
}

/// Emit all subject-vs-clip edge pairs whose AABBs overlap using sweep queries.
fn collect_overlapping_edge_pairs_sweep(
    subject_edges: &[EdgeAabb2d],
    clip_edges: &[EdgeAabb2d],
    mut emit: impl FnMut(usize, usize),
) {
    let index = ClipEdgeSweepIndex2d::build(clip_edges);
    let mut candidates: Vec<usize> = Vec::new();
    for (si, &se) in subject_edges.iter().enumerate() {
        index.query_overlaps(se, clip_edges, &mut candidates);
        for &cj in &candidates {
            emit(si, cj);
        }
    }
}

/// CDT-based polygon clipping using constrained Delaunay triangulation.
///
/// # Algorithm
///
/// 1. Weld subject and clip vertices into a shared unique set (`O(n)` expected).
/// 2. Find edge-edge intersection points using sweep broad-phase + exact segment test.
/// 3. Build a PSLG with shattered constraint edges from both polygons.
/// 4. Compute CDT.
/// 5. Classify each output triangle by testing its centroid against both
///    input polygons (winding number test).
pub fn cdt_clip(subject: &[[Real; 2]], clip: &[[Real; 2]], op: ClipOp) -> Vec<Vec<[Real; 2]>> {
    use crate::application::delaunay::Cdt;

    if subject.len() < 3 || clip.len() < 3 {
        return Vec::new();
    }

    let mut subj = subject.to_vec();
    let mut clp = clip.to_vec();
    ensure_ccw(&mut subj);
    ensure_ccw(&mut clp);

    let mut unique: Vec<[Real; 2]> = Vec::with_capacity(subj.len() + clp.len());
    let mut welder = SpatialHashWeld2d::new(WELD_TOL);

    let mut subj_indices: Vec<usize> = Vec::with_capacity(subj.len());
    let mut clip_indices: Vec<usize> = Vec::with_capacity(clp.len());

    for &p in &subj {
        subj_indices.push(welder.insert_or_weld(p, &mut unique));
    }
    for &p in &clp {
        clip_indices.push(welder.insert_or_weld(p, &mut unique));
    }

    let subj_edges = EdgeAabb2d::from_polygon(&subj);
    let clip_edges = EdgeAabb2d::from_polygon(&clp);

    collect_overlapping_edge_pairs_sweep(&subj_edges, &clip_edges, |si, cj| {
        let se = subj_edges[si];
        let ce = clip_edges[cj];
        if let Some((t, s)) = seg_intersect(se.a, se.b, ce.a, ce.b) {
            if t > INTERIOR_TOL
                && t < 1.0 - INTERIOR_TOL
                && s > INTERIOR_TOL
                && s < 1.0 - INTERIOR_TOL
            {
                let px = se.a[0] + t * (se.b[0] - se.a[0]);
                let py = se.a[1] + t * (se.b[1] - se.a[1]);
                let _ = welder.insert_or_weld([px, py], &mut unique);
            }
        }
    });

    if unique.len() < 3 {
        return Vec::new();
    }

    let mut pslg_edges = HashSet::<PlanarEdgeKey>::new();
    add_shattered_edges(&subj, &subj_indices, &unique, &mut pslg_edges);
    add_shattered_edges(&clp, &clip_indices, &unique, &mut pslg_edges);

    let pslg = build_pslg_from_points_and_edges(&unique, &pslg_edges);

    let cdt = Cdt::from_pslg(&pslg);
    let dt = cdt.triangulation();

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
    unique: &[[Real; 2]],
    pslg_edges: &mut HashSet<PlanarEdgeKey>,
) {
    let n = poly.len();
    for i in 0..n {
        let j = (i + 1) % n;
        let ri = indices[i];
        let rj = indices[j];
        if ri == rj {
            continue;
        }
        let on_edge =
            collect_points_on_segment_interior(unique, poly[i], poly[j], (ri, rj), 1e-8, 1e-12);
        insert_shattered_subedges(on_edge, pslg_edges);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;
    use std::collections::BTreeSet;

    fn brute_overlapping_pairs(
        subject: &[EdgeAabb2d],
        clip: &[EdgeAabb2d],
    ) -> BTreeSet<(usize, usize)> {
        let mut set = BTreeSet::new();
        for (i, &se) in subject.iter().enumerate() {
            for (j, &ce) in clip.iter().enumerate() {
                if se.aabb_overlaps(ce) {
                    let _ = set.insert((i, j));
                }
            }
        }
        set
    }

    fn sweep_overlapping_pairs(
        subject: &[EdgeAabb2d],
        clip: &[EdgeAabb2d],
    ) -> BTreeSet<(usize, usize)> {
        let mut set = BTreeSet::new();
        collect_overlapping_edge_pairs_sweep(subject, clip, |i, j| {
            let _ = set.insert((i, j));
        });
        set
    }

    proptest! {
        #[test]
        fn sweep_pairs_match_bruteforce(
            subj in prop::collection::vec((-40_i16..40_i16, -40_i16..40_i16), 3..10),
            clip in prop::collection::vec((-40_i16..40_i16, -40_i16..40_i16), 3..10)
        ) {
            let subj_poly: Vec<[Real;2]> = subj.into_iter()
                .map(|(x,y)| [x as Real * 0.25, y as Real * 0.25])
                .collect();
            let clip_poly: Vec<[Real;2]> = clip.into_iter()
                .map(|(x,y)| [x as Real * 0.25, y as Real * 0.25])
                .collect();

            let subj_edges = EdgeAabb2d::from_polygon(&subj_poly);
            let clip_edges = EdgeAabb2d::from_polygon(&clip_poly);

            let brute = brute_overlapping_pairs(&subj_edges, &clip_edges);
            let sweep = sweep_overlapping_pairs(&subj_edges, &clip_edges);
            prop_assert_eq!(sweep, brute);
        }

        #[test]
        fn spatial_hash_weld_is_idempotent_and_merges_close_points(
            x in -10_000.0_f64..10_000.0_f64,
            y in -10_000.0_f64..10_000.0_f64,
            dx in -2.0e-9_f64..2.0e-9_f64,
            dy in -2.0e-9_f64..2.0e-9_f64,
        ) {
            let mut unique = Vec::new();
            let mut welder = SpatialHashWeld2d::new(WELD_TOL);

            let p = [x, y];
            let q = [x + dx, y + dy];

            let i0 = welder.insert_or_weld(p, &mut unique);
            let i1 = welder.insert_or_weld(p, &mut unique);
            let i2 = welder.insert_or_weld(q, &mut unique);

            prop_assert_eq!(i0, i1);
            prop_assert_eq!(i0, i2);
            prop_assert_eq!(unique.len(), 1);
        }
    }

    #[test]
    fn sweep_finds_crossing_pair() {
        let subject = vec![[0.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
        let clip = vec![[0.0, 2.0], [2.0, 0.0], [2.0, 2.0]];

        let subject_edges = EdgeAabb2d::from_polygon(&subject);
        let clip_edges = EdgeAabb2d::from_polygon(&clip);

        let pairs = sweep_overlapping_pairs(&subject_edges, &clip_edges);
        assert!(
            !pairs.is_empty(),
            "crossing polygons must produce overlap candidates"
        );
    }
}
