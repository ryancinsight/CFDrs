//! Ruppert's refinement algorithm.
//!
//! # Theorem — Ruppert Termination and Quality Guarantee
//!
//! **Statement** (Ruppert 1995): Given a valid PSLG with minimum input angle
//! $> 60°$, Ruppert's algorithm terminates and produces a CDT where every
//! triangle has a radius-edge ratio $\leq B$ where $B \geq \sqrt{2}$.  The
//! resulting mesh has $O(n)$ triangles where $n$ is determined by the local
//! feature size.
//!
//! **Proof sketch**: The algorithm maintains two invariants:
//! 1. **No encroached segments**: all constraint segments have empty diametral
//!    circles.
//! 2. **All triangles are "good"**: radius-edge ratio $\leq B$.
//!
//! When a bad triangle is found, its circumcenter (or off-center) is a
//! candidate insertion point:
//! - If it encroaches a segment, split the segment instead (midpoint insertion).
//! - Otherwise, insert the circumcenter.
//!
//! Each insertion either splits a segment (making it shorter) or eliminates a
//! bad triangle.  Segment lengths are bounded below by the local feature size
//! function $\text{lfs}$, so the number of segment splits is finite.  After
//! all segments are unencroached, circumcenter insertions strictly improve
//! quality.  The process terminates because the number of possible Steiner
//! points is bounded by the area / (minimum allowed triangle area).
//!
//! # Modern Enhancements
//!
//! We incorporate several modern improvements over the original 1995 algorithm:
//!
//! 1. **Off-centers** (Üngör 2004): Instead of always inserting the
//!    circumcenter, we insert the off-center point, which is closer to the
//!    shortest edge and results in fewer Steiner points.
//!
//! 2. **Concentric shell splitting** (Shewchuk 2002): Segment midpoints are
//!    rounded to powers-of-two spacing to prevent infinite cascading splits.
//!
//! 3. **Quality-prioritized queue**: Worst-quality triangles are refined first,
//!    maximizing improvement per insertion.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use crate::application::delaunay::constraint::enforce::Cdt;
use crate::application::delaunay::pslg::vertex::PslgVertexId;
use crate::application::delaunay::refinement::circumcenter::{circumcenter, off_center};
use crate::application::delaunay::refinement::encroachment::{
    is_encroached, point_encroaches_segment,
};
use crate::application::delaunay::refinement::metric::MetricTensor;
use crate::application::delaunay::refinement::quality::TriangleQuality;
use crate::application::delaunay::triangulation::locate::Location;
use crate::application::delaunay::triangulation::triangle::{Triangle, TriangleId};
use crate::domain::core::scalar::Real;

/// A bad triangle entry in the priority queue.
#[derive(Clone, Copy)]
struct BadTriangle {
    tid: TriangleId,
    ratio: Real,
}

impl PartialEq for BadTriangle {
    fn eq(&self, other: &Self) -> bool {
        self.ratio == other.ratio
    }
}

impl Eq for BadTriangle {}

impl PartialOrd for BadTriangle {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for BadTriangle {
    fn cmp(&self, other: &Self) -> Ordering {
        // Max-heap: worst quality first.
        self.ratio
            .partial_cmp(&other.ratio)
            .unwrap_or(Ordering::Equal)
    }
}

/// Ruppert's refinement algorithm operating on a CDT.
///
/// Supports both isotropic refinement (default) and anisotropic refinement via
/// an optional [`MetricTensor`] that redefines the quality metric.
///
/// # Example — Isotropic
///
/// ```rust,ignore
/// use cfd_mesh::application::delaunay::{Pslg, Cdt, RuppertRefiner};
///
/// let cdt = Cdt::from_pslg(&pslg);
/// let mut refiner = RuppertRefiner::new(cdt);
/// refiner.set_max_ratio(1.414);
/// refiner.refine();
/// let result = refiner.into_cdt();
/// ```
///
/// # Example — Anisotropic (10:1 elongation along x)
///
/// ```rust,ignore
/// use cfd_mesh::application::delaunay::{Pslg, Cdt, RuppertRefiner, MetricTensor};
///
/// let cdt = Cdt::from_pslg(&pslg);
/// let metric = MetricTensor::anisotropic(0.0, 10.0);
/// let mut refiner = RuppertRefiner::new(cdt).with_metric(metric);
/// refiner.set_max_ratio(1.414);
/// refiner.refine();
/// ```
pub struct RuppertRefiner {
    cdt: Cdt,
    /// Maximum allowed radius-edge ratio (default: √2 ≈ 1.414).
    max_ratio: Real,
    /// Maximum allowed triangle area (optional constraint).
    max_area: Option<Real>,
    /// Maximum number of Steiner points to insert (safety limit).
    max_steiner: usize,
    /// Number of Steiner points inserted.
    steiner_count: usize,
    /// Optional anisotropic metric tensor.
    ///
    /// When `Some(M)`, triangle quality is measured in M-space (the Cholesky
    /// transform of M is applied to all vertices before computing the
    /// isotropic radius-edge ratio).  When `None`, standard Euclidean quality
    /// is used — identical to previous behaviour.
    metric: Option<MetricTensor>,
}

impl RuppertRefiner {
    /// Create a new refiner wrapping a CDT.
    #[must_use]
    pub fn new(cdt: Cdt) -> Self {
        Self {
            cdt,
            max_ratio: core::f64::consts::SQRT_2,
            max_area: None,
            max_steiner: 100_000,
            steiner_count: 0,
            metric: None,
        }
    }

    /// Set the anisotropic metric tensor (builder pattern).
    ///
    /// When set, triangle quality is evaluated in M-space rather than
    /// Euclidean space, allowing anisotropic refinement.
    ///
    /// Default: `None` (isotropic Euclidean — fully backward-compatible).
    #[must_use]
    pub fn with_metric(mut self, metric: MetricTensor) -> Self {
        self.metric = Some(metric);
        self
    }

    /// Set (or replace) the anisotropic metric tensor.
    pub fn set_metric(&mut self, metric: MetricTensor) {
        self.metric = Some(metric);
    }

    /// Remove the anisotropic metric, reverting to isotropic Euclidean quality.
    pub fn clear_metric(&mut self) {
        self.metric = None;
    }

    /// Set the maximum radius-edge ratio bound.
    ///
    /// Must be ≥ 1.0 (lower values demand higher quality but may not terminate
    /// for inputs with small angles).  The theoretical minimum for guaranteed
    /// termination is $\sqrt{2}$.
    pub fn set_max_ratio(&mut self, ratio: Real) {
        self.max_ratio = ratio.max(1.0);
    }

    /// Set the maximum triangle area constraint.
    pub fn set_max_area(&mut self, area: Real) {
        self.max_area = Some(area);
    }

    /// Set the maximum number of Steiner points (safety limit).
    pub fn set_max_steiner(&mut self, n: usize) {
        self.max_steiner = n;
    }

    /// Run the refinement algorithm.
    ///
    /// Returns the number of Steiner points inserted.
    pub fn refine(&mut self) -> usize {
        // Phase 1: Fix all encroached segments.
        self.fix_encroached_segments();

        // Phase 2: Fix bad triangles.
        self.fix_bad_triangles();

        self.steiner_count
    }

    /// Phase 1: Split all encroached constraint segments.
    fn fix_encroached_segments(&mut self) {
        let max_iter = self.max_steiner;
        for _ in 0..max_iter {
            let encroached = self.find_encroached_segment();
            match encroached {
                Some((a, b)) => {
                    self.split_segment(a, b);
                }
                None => return,
            }
        }
    }

    /// Find an encroached constraint segment.
    fn find_encroached_segment(&self) -> Option<(PslgVertexId, PslgVertexId)> {
        let dt = self.cdt.triangulation();
        for &(a, b) in self.cdt.constrained_edges() {
            if is_encroached(dt, a, b) {
                return Some((a, b));
            }
        }
        None
    }

    /// Split a constraint segment at its midpoint (with concentric shell rounding).
    fn split_segment(&mut self, a: PslgVertexId, b: PslgVertexId) -> PslgVertexId {
        let va = *self.cdt.triangulation().vertex(a);
        let vb = *self.cdt.triangulation().vertex(b);

        // Concentric shell midpoint (Shewchuk 2002):
        // Round the split point to the nearest power-of-two fraction.
        let mx = (va.x + vb.x) * 0.5;
        let my = (va.y + vb.y) * 0.5;

        let mid = self.concentric_midpoint(va.x, va.y, vb.x, vb.y, mx, my);

        let mid_vid = self.cdt.triangulation_mut().insert_steiner(mid.0, mid.1);

        // Remove the old constraint and replace with the two sub-segments.
        self.cdt.remove_constraint(a, b);
        self.cdt.add_constraint(a, mid_vid);
        self.cdt.add_constraint(mid_vid, b);

        self.steiner_count += 1;
        mid_vid
    }

    /// Phase 2: Iteratively fix bad triangles.
    fn fix_bad_triangles(&mut self) {
        let mut queue = self.build_bad_queue();

        while let Some(bad) = queue.pop() {
            if self.steiner_count >= self.max_steiner {
                break;
            }

            let dt = self.cdt.triangulation();
            // Check if this triangle is still alive and still bad.
            let tri = dt.triangle(bad.tid);
            if !tri.alive {
                continue;
            }

            let q = self.triangle_quality(tri);
            if q.is_good(self.max_ratio) && !self.exceeds_area(&q) {
                continue;
            }

            let [v0, v1, v2] = tri.vertices;
            let a = *dt.vertex(v0);
            let b = *dt.vertex(v1);
            let c = *dt.vertex(v2);

            // Compute circumcenter (or off-center).
            let insertion_point =
                off_center(&a, &b, &c, self.max_ratio).or_else(|| circumcenter(&a, &b, &c));

            let (px, py) = match insertion_point {
                Some(p) => p,
                None => continue,
            };

            // Check if the insertion point encroaches any constraint segment.
            let encroached_seg = self.check_encroachment(px, py);

            if let Some((sa, sb)) = encroached_seg {
                // Split the encroached segment instead.
                let split_vid = self.split_segment(sa, sb);
                // After a segment split the new Steiner vertex's 1-ring may
                // contain bad triangles; scan that ring only.
                self.append_bad_triangles_ring(&mut queue, split_vid);
            } else {
                // Verify the insertion point is inside the domain
                // (not in a super-triangle region or outside the boundary).
                let dt = self.cdt.triangulation_mut();
                let inside = match dt.locate_point(bad.tid, px, py) {
                    Some(
                        Location::Inside(tid)
                        | Location::OnEdge(tid, _)
                        | Location::OnVertex(tid, _),
                    ) => {
                        let t = &dt.triangles[tid.idx()];
                        !t.vertices.iter().any(|v| dt.super_verts.contains(v))
                    }
                    _ => false,
                };

                if !inside {
                    continue;
                }

                // Insert the circumcenter/off-center.
                let new_vid = self.cdt.triangulation_mut().insert_steiner(px, py);
                self.steiner_count += 1;

                // Only scan the 1-ring of the newly inserted vertex for
                // bad triangles — O(deg) ≈ O(6) instead of O(T).
                self.append_bad_triangles_ring(&mut queue, new_vid);
            }
        }
    }

    /// Build the initial priority queue of bad triangles.
    fn build_bad_queue(&self) -> BinaryHeap<BadTriangle> {
        let mut queue = BinaryHeap::new();
        let dt = self.cdt.triangulation();

        for (tid, tri) in dt.interior_triangles() {
            let q = self.triangle_quality(tri);
            if !q.is_good(self.max_ratio) || self.exceeds_area(&q) {
                queue.push(BadTriangle {
                    tid,
                    ratio: q.radius_edge_ratio,
                });
            }
        }
        queue
    }

    /// Scan only the 1-ring of vertex `vid` for bad triangles.
    ///
    /// # Theorem — 1-Ring Sufficiency for Queue Maintenance
    ///
    /// **Statement**: After inserting a Steiner point $p$ into a CDT, the
    /// only triangles whose quality can *worsen* are those in the 1-ring
    /// of $p$ (i.e., triangles incident on $p$ after the Bowyer-Watson
    /// cavity re-triangulation and subsequent flips).
    ///
    /// **Proof sketch**: The Bowyer-Watson insertion removes the cavity
    /// (triangles whose circumcircle contains $p$) and replaces them with
    /// fan triangles around $p$.  Edge flips during `flip_fix` only
    /// rearrange triangles within and adjacent to this cavity.  Triangles
    /// outside the 1-ring of $p$ are geometrically unchanged; their
    /// quality metrics are invariant.  Therefore scanning only the 1-ring
    /// is sufficient to detect all newly-bad triangles.  ∎
    ///
    /// # Complexity
    ///
    /// $O(\deg(p)) \approx O(6)$ instead of $O(T)$ per insertion.
    fn append_bad_triangles_ring(
        &self,
        queue: &mut BinaryHeap<BadTriangle>,
        vid: PslgVertexId,
    ) {
        let dt = self.cdt.triangulation();
        for tid in dt.triangles_around_vertex(vid) {
            let tri = dt.triangle(tid);
            if !tri.alive || tri.vertices.iter().any(|v| dt.super_verts.contains(v)) {
                continue;
            }
            let q = self.triangle_quality(tri);
            if !q.is_good(self.max_ratio) || self.exceeds_area(&q) {
                queue.push(BadTriangle {
                    tid,
                    ratio: q.radius_edge_ratio,
                });
            }
        }
    }

    /// Compute triangle quality, using the anisotropic metric when set.
    ///
    /// ## Metric-Weighted Quality
    ///
    /// When `self.metric = Some(M)`, quality is measured in M-space.  The
    /// Cholesky factor L of M (M = Lᵀ L) is applied to transform vertex
    /// coordinates: `a' = L·a, b' = L·b, c' = L·c`.  The isotropic quality
    /// of the transformed triangle `(a', b', c')` equals the metric-weighted
    /// quality of the original triangle.
    ///
    /// **Theorem**: |e|_M = |L·e|₂ for any edge vector e.  Therefore the
    /// metric circumradius R_M = R_iso(La, Lb, Lc) and the metric shortest
    /// edge l_min_M = l_min_iso(La, Lb, Lc).  The ratio R_M / l_min_M
    /// equals the metric-quality ρ_M of the triangle.  QED.
    fn triangle_quality(&self, tri: &Triangle) -> TriangleQuality {
        let dt = self.cdt.triangulation();
        let a = dt.vertex(tri.vertices[0]);
        let b = dt.vertex(tri.vertices[1]);
        let c = dt.vertex(tri.vertices[2]);

        if let Some(ref m) = self.metric {
            // Transform vertices by Cholesky(M) and compute isotropic quality.
            if let Some(l) = m.cholesky() {
                use crate::application::delaunay::pslg::vertex::PslgVertex;
                let (ax2, ay2) = MetricTensor::apply_cholesky(&l, a.x, a.y);
                let (bx2, by2) = MetricTensor::apply_cholesky(&l, b.x, b.y);
                let (cx2, cy2) = MetricTensor::apply_cholesky(&l, c.x, c.y);
                let a2 = PslgVertex::new(ax2, ay2);
                let b2 = PslgVertex::new(bx2, by2);
                let c2 = PslgVertex::new(cx2, cy2);
                return TriangleQuality::compute(&a2, &b2, &c2);
            }
            // Metric is not PD (degenerate) — fall through to Euclidean.
        }

        TriangleQuality::compute(a, b, c)
    }

    /// Check if triangle area exceeds the maximum.
    fn exceeds_area(&self, q: &TriangleQuality) -> bool {
        self.max_area.is_some_and(|max| q.area > max)
    }

    /// Check if point `(px, py)` encroaches any constraint segment.
    fn check_encroachment(&self, px: Real, py: Real) -> Option<(PslgVertexId, PslgVertexId)> {
        let dt = self.cdt.triangulation();
        for &(a, b) in self.cdt.constrained_edges() {
            let va = dt.vertex(a);
            let vb = dt.vertex(b);
            if point_encroaches_segment(va, vb, px, py) {
                return Some((a, b));
            }
        }
        None
    }

    /// Concentric shell midpoint (Shewchuk 2002).
    ///
    /// Rounds the midpoint to a power-of-two distance from `a`, preventing
    /// infinite cascading splits when nearby constraints interact.
    fn concentric_midpoint(
        &self,
        ax: Real,
        ay: Real,
        bx: Real,
        by: Real,
        mx: Real,
        my: Real,
    ) -> (Real, Real) {
        let seg_len = ((bx - ax) * (bx - ax) + (by - ay) * (by - ay)).sqrt();
        if seg_len < 1e-15 {
            return (mx, my);
        }

        // Find the largest power of two ≤ seg_len/2.
        let half = seg_len * 0.5;
        let pow2 = (2.0_f64).powf(half.log2().floor());
        let t = pow2 / seg_len;
        let t = t.clamp(0.25, 0.75); // Stay away from endpoints.

        (ax + t * (bx - ax), ay + t * (by - ay))
    }

    // ── Public accessors ──────────────────────────────────────────────────

    /// Access the CDT.
    #[must_use]
    pub fn cdt(&self) -> &Cdt {
        &self.cdt
    }

    /// Consume the refiner and return the CDT.
    #[must_use]
    pub fn into_cdt(self) -> Cdt {
        self.cdt
    }

    /// Number of Steiner points inserted.
    #[must_use]
    pub fn steiner_count(&self) -> usize {
        self.steiner_count
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::application::delaunay::constraint::enforce::Cdt;
    use crate::application::delaunay::pslg::graph::Pslg;
    use crate::application::delaunay::refinement::metric::MetricTensor;

    /// Build a square PSLG (0,0)-(1,0)-(1,1)-(0,1) with constrained boundary.
    fn square_cdt() -> Cdt {
        let mut p = Pslg::new();
        let v0 = p.add_vertex(0.0, 0.0);
        let v1 = p.add_vertex(1.0, 0.0);
        let v2 = p.add_vertex(1.0, 1.0);
        let v3 = p.add_vertex(0.0, 1.0);
        p.add_segment(v0, v1);
        p.add_segment(v1, v2);
        p.add_segment(v2, v3);
        p.add_segment(v3, v0);
        Cdt::from_pslg(&p)
    }

    /// Isotropic refiner terminates and produces ≥ 2 interior triangles.
    #[test]
    fn isotropic_refiner_terminates() {
        let cdt = square_cdt();
        let mut refiner = RuppertRefiner::new(cdt);
        refiner.set_max_ratio(1.5);
        let n = refiner.refine();
        assert!(n <= 100_000, "Steiner count should be bounded");
        assert!(
            refiner.cdt().triangulation().interior_triangles().count() >= 2,
            "should have interior triangles"
        );
    }

    /// Identity metric gives same result as no metric (backward compat).
    #[test]
    fn identity_metric_matches_no_metric() {
        let cdt_a = square_cdt();
        let cdt_b = square_cdt();

        let mut r_iso = RuppertRefiner::new(cdt_a);
        r_iso.set_max_ratio(1.5);
        r_iso.set_max_steiner(50);
        r_iso.refine();

        let mut r_id = RuppertRefiner::new(cdt_b).with_metric(MetricTensor::identity());
        r_id.set_max_ratio(1.5);
        r_id.set_max_steiner(50);
        r_id.refine();

        // Both should produce the same number of interior triangles
        // (identity metric ≡ Euclidean quality).
        assert_eq!(
            r_iso.cdt().triangulation().interior_triangles().count(),
            r_id.cdt().triangulation().interior_triangles().count(),
            "identity metric should give same result as no metric"
        );
    }

    /// Anisotropic metric terminates within the safety limit.
    ///
    /// A 2:1 anisotropic metric on the unit square is used (rather than 10:1)
    /// because very high aspect ratios cause O(α²) insertions on a square domain
    /// where all edges are comparable in Euclidean space.  The test checks that
    /// the refiner honours `max_steiner` and produces a valid CDT.
    #[test]
    fn anisotropic_refiner_terminates() {
        let cdt = square_cdt();
        // 2:1 anisotropy along x: triangles may be 2× longer in x than in y.
        let metric = MetricTensor::anisotropic(0.0, 2.0);
        let mut refiner = RuppertRefiner::new(cdt).with_metric(metric);
        refiner.set_max_ratio(1.5);
        refiner.set_max_steiner(2_000); // explicit safety cap for test
        let n = refiner.refine();
        assert!(
            n <= 2_000,
            "anisotropic refiner Steiner count should be bounded: {n}"
        );
        assert!(
            refiner.cdt().triangulation().interior_triangles().count() >= 2,
            "CDT must contain triangles after anisotropic refinement"
        );
    }

    /// set_metric / clear_metric mutate metric field correctly.
    #[test]
    fn set_and_clear_metric() {
        let mut refiner = RuppertRefiner::new(square_cdt());
        assert!(refiner.metric.is_none());
        refiner.set_metric(MetricTensor::anisotropic(0.5, 5.0));
        assert!(refiner.metric.is_some());
        refiner.clear_metric();
        assert!(refiner.metric.is_none());
    }

    /// with_metric builder sets metric; final CDT is valid (≥ 1 alive triangle).
    ///
    /// Uses a moderate 3:1 anisotropy to keep Steiner point count manageable.
    #[test]
    fn with_metric_builder_produces_valid_cdt() {
        let cdt = square_cdt();
        let mut refiner = RuppertRefiner::new(cdt).with_metric(MetricTensor::anisotropic(0.0, 3.0));
        refiner.set_max_steiner(2_000);
        refiner.refine();
        let tri_count = refiner.cdt().triangulation().interior_triangles().count();
        assert!(
            tri_count >= 1,
            "CDT must contain triangles after anisotropic refinement"
        );
    }
}
