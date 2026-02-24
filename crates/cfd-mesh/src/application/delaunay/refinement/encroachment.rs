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

use crate::application::delaunay::pslg::vertex::{PslgVertex, PslgVertexId};
use crate::application::delaunay::triangulation::bowyer_watson::DelaunayTriangulation;
use crate::domain::core::scalar::Real;

/// Check if the segment `(a, b)` is encroached by any vertex in the
/// triangulation.
///
/// A segment is encroached if a vertex (other than `a` or `b`) lies
/// strictly inside the diametral circle of `(a, b)`.
pub fn is_encroached(dt: &DelaunayTriangulation, a: PslgVertexId, b: PslgVertexId) -> bool {
    encroaching_vertex(dt, a, b).is_some()
}

/// Find a vertex that encroaches the segment `(a, b)`.
///
/// Returns `Some(vid)` if a vertex is found strictly inside the diametral
/// circle, or `None` if the segment is unencroached.
pub fn encroaching_vertex(
    dt: &DelaunayTriangulation,
    a: PslgVertexId,
    b: PslgVertexId,
) -> Option<PslgVertexId> {
    let va = dt.vertex(a);
    let vb = dt.vertex(b);

    // Diametral circle center and radius².
    let cx = (va.x + vb.x) * 0.5;
    let cy = (va.y + vb.y) * 0.5;
    let rsq = va.dist_sq(vb) * 0.25;

    for (_tid, tri) in dt.all_alive_triangles() {
        for &vid in &tri.vertices {
            if vid == a || vid == b {
                continue;
            }
            // Skip super-triangle vertices.
            if dt.super_verts.contains(&vid) {
                continue;
            }

            let v = dt.vertex(vid);
            let dx = v.x - cx;
            let dy = v.y - cy;
            let dsq = dx * dx + dy * dy;

            if dsq < rsq - 1e-14 {
                return Some(vid);
            }
        }
    }
    None
}

/// Check if a point `(px, py)` encroaches the segment between vertices `a` and `b`.
pub fn point_encroaches_segment(va: &PslgVertex, vb: &PslgVertex, px: Real, py: Real) -> bool {
    let cx = (va.x + vb.x) * 0.5;
    let cy = (va.y + vb.y) * 0.5;
    let rsq = va.dist_sq(vb) * 0.25;

    let dx = px - cx;
    let dy = py - cy;
    let dsq = dx * dx + dy * dy;

    dsq < rsq - 1e-14
}
