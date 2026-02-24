//! Circumcenter computation for triangles.
//!
//! # Theorem — Circumcenter as Voronoi Vertex
//!
//! **Statement**: The circumcenter of a Delaunay triangle is the Voronoi
//! vertex dual to the three generators.  It is equidistant from all three
//! vertices of the triangle.
//!
//! **Proof sketch**: The circumcenter is defined as the intersection of
//! the perpendicular bisectors of any two edges.  By construction, it is
//! equidistant from all three vertices, hence it is the center of the
//! unique circumscribed circle.

use crate::application::delaunay::pslg::vertex::PslgVertex;
use crate::domain::core::scalar::Real;

/// Compute the circumcenter of triangle `(a, b, c)`.
///
/// Returns `(cx, cy)`.
///
/// # Panics
///
/// Returns `None` if the three points are collinear (degenerate triangle).
pub fn circumcenter(a: &PslgVertex, b: &PslgVertex, c: &PslgVertex) -> Option<(Real, Real)> {
    // Using the formula from Shewchuk's "Triangle" library.
    let ax = a.x - c.x;
    let ay = a.y - c.y;
    let bx = b.x - c.x;
    let by = b.y - c.y;

    let d = 2.0 * (ax * by - ay * bx);
    if d.abs() < 1e-30 {
        return None; // Degenerate.
    }

    let a_sq = ax * ax + ay * ay;
    let b_sq = bx * bx + by * by;

    let cx = (a_sq * by - b_sq * ay) / d + c.x;
    let cy = (b_sq * ax - a_sq * bx) / d + c.y;

    Some((cx, cy))
}

/// Compute the off-center for a skinny triangle.
///
/// The off-center lies on the line from the shortest-edge midpoint through
/// the circumcenter, at a distance equal to the quality bound times the
/// shortest edge length.  This reduces the number of Steiner points needed
/// compared to always inserting at the circumcenter.
///
/// # Reference
///
/// Üngör (2004), "Off-centers: A new type of Steiner points for computing
/// size-optimal quality-guaranteed Delaunay triangulations."
pub fn off_center(
    a: &PslgVertex,
    b: &PslgVertex,
    c: &PslgVertex,
    ratio_bound: Real,
) -> Option<(Real, Real)> {
    let (ccx, ccy) = circumcenter(a, b, c)?;

    // Find the shortest edge.
    let edges = [
        ((a, b), a.dist_sq(b)),
        ((b, c), b.dist_sq(c)),
        ((c, a), c.dist_sq(a)),
    ];

    let (shortest, shortest_sq) = edges
        .iter()
        .min_by(|(_, d1), (_, d2)| d1.partial_cmp(d2).unwrap())
        .unwrap();

    let shortest_len = shortest_sq.sqrt();
    let mid_x = (shortest.0.x + shortest.1.x) * 0.5;
    let mid_y = (shortest.0.y + shortest.1.y) * 0.5;

    // Vector from midpoint to circumcenter.
    let dx = ccx - mid_x;
    let dy = ccy - mid_y;
    let dist_to_cc = (dx * dx + dy * dy).sqrt();

    // Desired distance: ratio_bound * shortest_len.
    let desired = ratio_bound * shortest_len;

    if dist_to_cc <= desired {
        // Circumcenter is close enough; use it directly.
        Some((ccx, ccy))
    } else {
        // Place the off-center at distance `desired` from the midpoint.
        let scale = desired / dist_to_cc;
        Some((mid_x + dx * scale, mid_y + dy * scale))
    }
}
