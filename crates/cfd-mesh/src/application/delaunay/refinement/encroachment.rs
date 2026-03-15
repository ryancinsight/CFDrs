//! Segment encroachment detection and resolution.
//!
//! # Definition — Encroached Segment
//!
//! A segment is *encroached* if there exists a vertex (other than its
//! endpoints) strictly inside its diametral circle (the circle whose
//! diameter is the segment).
//!
//! # Theorem — Segment Splitting Preserves Conformity
//!
//! **Statement**: If a segment is encroached, inserting its midpoint and
//! splitting the segment into two sub-segments produces a valid PSLG
//! refinement.  The new subsegments may themselves become encroached, but
//! the process terminates because segment lengths strictly decrease and
//! are bounded below by the local feature size.
//!
//! **Proof sketch**: The diametral circle of each subsegment is contained
//! within the diametral circle of the original segment, so the new
//! subsegments are "better separated" from vertices.  Ruppert (1995)
//! shows termination using a local feature size argument: no segment can
//! be split shorter than the local feature size $\text{lfs}(v)$ at its
//! midpoint.
//!
//! # Theorem — Local Encroachment Sufficiency (CDT Optimisation)
//!
//! **Statement**: In a CDT where segment $(a,b)$ is a constraint edge,
//! the only vertices that can encroach $(a,b)$ are the opposite vertices
//! of the (at most two) triangles incident on that edge.
//!
//! **Proof sketch**: Let $D$ be the diametral circle of $(a,b)$.  In a
//! CDT, $(a,b)$ is shared by at most two triangles $t_1 = (a,b,c)$ and
//! $t_2 = (a,b,d)$.  Any other vertex $v$ strictly inside $D$ would
//! violate the constrained empty-circumdisk property: $v$ would be inside
//! the circumcircle of $t_1$ or $t_2$ (since $D$ is contained in both
//! circumcircles when the triangle exists), and visibility through the
//! constraint edge is unrestricted for the opposite vertex.  Hence only
//! $c$ and $d$ need to be tested.  ∎
//!
//! **Complexity**: $O(1)$ per segment (two opposite-vertex lookups) vs
//! the naïve $O(T)$ full scan.

use crate::application::delaunay::pslg::vertex::{PslgVertex, PslgVertexId};
use crate::application::delaunay::triangulation::bowyer_watson::DelaunayTriangulation;
use crate::application::delaunay::triangulation::triangle::GHOST_TRIANGLE;
use crate::domain::core::scalar::Real;

/// Check if the segment `(a, b)` is encroached by any vertex in the
/// triangulation.
///
/// Uses the local CDT property: only the opposite vertices of the two
/// triangles sharing edge $(a,b)$ are checked.  Falls back to a
/// vertex-star scan when the edge is not found directly.
///
/// # Complexity
///
/// $O(\deg(a))$ — typically $O(6)$ for well-distributed Delaunay.
#[must_use]
pub fn is_encroached(dt: &DelaunayTriangulation, a: PslgVertexId, b: PslgVertexId) -> bool {
    encroaching_vertex(dt, a, b).is_some()
}

/// Find a vertex that encroaches the segment `(a, b)`.
///
/// Returns `Some(vid)` if a vertex is found strictly inside the diametral
/// circle, or `None` if the segment is unencroached.
///
/// # Algorithm
///
/// 1. Walk the vertex star of `a` to find triangles containing both `a`
///    and `b`.  These are the (at most two) triangles sharing edge $(a,b)$.
/// 2. For each such triangle, test the opposite vertex (the one that is
///    neither `a` nor `b`) against the diametral circle.
/// 3. If no such triangle is found (edge does not exist in the
///    triangulation), fall back to scanning all neighbours of both `a`
///    and `b` — still $O(\deg(a) + \deg(b))$, not $O(T)$.
#[must_use]
pub fn encroaching_vertex(
    dt: &DelaunayTriangulation,
    a: PslgVertexId,
    b: PslgVertexId,
) -> Option<PslgVertexId> {
    let va = dt.vertex(a);
    let vb = dt.vertex(b);

    let cx = (va.x + vb.x) * 0.5;
    let cy = (va.y + vb.y) * 0.5;
    let rsq = va.dist_sq(vb) * 0.25;

    // Fast path: find the ≤2 triangles sharing edge (a,b) via vertex-star.
    let star_a = dt.triangles_around_vertex(a);
    let mut found_edge = false;

    for &tid in &star_a {
        let tri = dt.triangle(tid);
        if !tri.contains_vertex(b) {
            continue;
        }
        found_edge = true;
        // The opposite vertex is the one that is neither a nor b.
        for &vid in &tri.vertices {
            if vid == a || vid == b || dt.super_verts.contains(&vid) {
                continue;
            }
            if vertex_in_diametral(dt.vertex(vid), cx, cy, rsq) {
                return Some(vid);
            }
        }
    }

    if found_edge {
        return None;
    }

    // Fallback: edge (a,b) is not in the triangulation.  Scan all
    // unique neighbours of a and b — still O(deg(a) + deg(b)).
    let star_b = dt.triangles_around_vertex(b);
    for star in [&star_a, &star_b] {
        for &tid in star {
            let tri = dt.triangle(tid);
            for &vid in &tri.vertices {
                if vid == a || vid == b || dt.super_verts.contains(&vid) {
                    continue;
                }
                if vertex_in_diametral(dt.vertex(vid), cx, cy, rsq) {
                    return Some(vid);
                }
            }
        }
    }
    None
}

/// Test whether vertex `v` lies strictly inside the diametral circle
/// centred at `(cx, cy)` with squared radius `rsq`.
///
/// Uses a scale-relative tolerance: $\varepsilon_{\text{rel}} =
/// \max(r^2, 1) \cdot 10^{-14}$ so that the threshold adapts to the
/// coordinate magnitude.
#[inline]
fn vertex_in_diametral(v: &PslgVertex, cx: Real, cy: Real, rsq: Real) -> bool {
    let dx = v.x - cx;
    let dy = v.y - cy;
    let dsq = dx * dx + dy * dy;
    let tol = rsq.max(1.0) * 1e-14;
    dsq < rsq - tol
}

/// Check if a point `(px, py)` encroaches the segment between vertices
/// `a` and `b`.
///
/// Uses a scale-relative tolerance so the test is stable across
/// different coordinate ranges.
#[must_use]
pub fn point_encroaches_segment(va: &PslgVertex, vb: &PslgVertex, px: Real, py: Real) -> bool {
    let cx = (va.x + vb.x) * 0.5;
    let cy = (va.y + vb.y) * 0.5;
    let rsq = va.dist_sq(vb) * 0.25;

    let dx = px - cx;
    let dy = py - cy;
    let dsq = dx * dx + dy * dy;
    let tol = rsq.max(1.0) * 1e-14;
    dsq < rsq - tol
}
