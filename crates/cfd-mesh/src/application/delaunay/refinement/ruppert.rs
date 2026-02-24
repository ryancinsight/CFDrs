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
use crate::application::delaunay::refinement::quality::TriangleQuality;
use crate::application::delaunay::triangulation::locate::{locate, Location};
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
/// # Example
///
/// ```rust,ignore
/// use cfd_mesh::application::delaunay::{Pslg, Cdt, RuppertRefiner};
///
/// let mut pslg = Pslg::new();
/// pslg.add_vertex(0.0, 0.0);
/// pslg.add_vertex(10.0, 0.0);
/// pslg.add_vertex(10.0, 10.0);
/// pslg.add_vertex(0.0, 10.0);
/// // Add boundary segments...
///
/// let cdt = Cdt::from_pslg(&pslg);
/// let mut refiner = RuppertRefiner::new(cdt);
/// refiner.set_max_ratio(1.414);
/// refiner.refine();
///
/// let result = refiner.into_cdt();
/// // All triangles now satisfy the quality bound.
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
}

impl RuppertRefiner {
    /// Create a new refiner wrapping a CDT.
    pub fn new(cdt: Cdt) -> Self {
        Self {
            cdt,
            max_ratio: core::f64::consts::SQRT_2,
            max_area: None,
            max_steiner: 100_000,
            steiner_count: 0,
        }
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
    fn split_segment(&mut self, a: PslgVertexId, b: PslgVertexId) {
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

            match encroached_seg {
                Some((sa, sb)) => {
                    // Split the encroached segment instead.
                    self.split_segment(sa, sb);
                }
                None => {
                    // Verify the insertion point is inside the domain
                    // (not in a super-triangle region or outside the boundary).
                    let dt = self.cdt.triangulation();
                    let inside = match locate(dt.vertices(), dt.triangles_slice(), bad.tid, px, py)
                    {
                        Some(Location::Inside(tid))
                        | Some(Location::OnEdge(tid, _))
                        | Some(Location::OnVertex(tid, _)) => {
                            let t = dt.triangle(tid);
                            !t.vertices.iter().any(|v| dt.super_verts.contains(v))
                        }
                        _ => false,
                    };

                    if !inside {
                        continue;
                    }

                    // Insert the circumcenter/off-center.
                    self.cdt.triangulation_mut().insert_steiner(px, py);
                    self.steiner_count += 1;
                }
            }

            // Re-scan for new bad triangles near the insertion.
            self.append_bad_triangles(&mut queue);
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

    /// Scan all interior triangles and add bad ones to the queue.
    fn append_bad_triangles(&self, queue: &mut BinaryHeap<BadTriangle>) {
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
    }

    /// Compute triangle quality.
    fn triangle_quality(&self, tri: &Triangle) -> TriangleQuality {
        let dt = self.cdt.triangulation();
        let a = dt.vertex(tri.vertices[0]);
        let b = dt.vertex(tri.vertices[1]);
        let c = dt.vertex(tri.vertices[2]);
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
    pub fn cdt(&self) -> &Cdt {
        &self.cdt
    }

    /// Consume the refiner and return the CDT.
    pub fn into_cdt(self) -> Cdt {
        self.cdt
    }

    /// Number of Steiner points inserted.
    pub fn steiner_count(&self) -> usize {
        self.steiner_count
    }
}
