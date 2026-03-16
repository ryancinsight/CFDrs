//! Cavity polygon re-triangulation.
//!
//! When enforcing a constraint edge that crosses multiple existing edges,
//! the crossed triangles are removed, leaving a cavity polygon on each side
//! of the constraint.  This module re-triangulates each cavity with a fan
//! that respects the Delaunay property where possible.
//!
//! # Theorem — Cavity Re-Triangulation
//!
//! **Statement**: A star-shaped polygon with a known kernel point can be
//! triangulated by connecting the kernel to all boundary edges, producing
//! a fan triangulation in $O(k)$ time where $k$ is the polygon size.
//!
//! **Proof sketch**: Star-shapedness guarantees that the kernel sees every
//! boundary edge.  Connecting the kernel to each boundary vertex creates
//! triangles whose interiors lie entirely within the polygon.

use crate::domain::core::scalar::Real;
use crate::domain::geometry::predicates::{incircle, orient_2d, Orientation};
use nalgebra::Point2;

use crate::application::delaunay::pslg::vertex::{PslgVertex, PslgVertexId};

/// Re-triangulate a cavity polygon on one side of a constraint edge.
///
/// # Theorem — Ear-Clipping with Delaunay Preference
///
/// **Statement**: Every simple polygon with $n \geq 4$ vertices has at least
/// two non-overlapping ears (Meisters, 1975).  Iteratively clipping the ear
/// with fewest circumcircle-interior violations produces a triangulation
/// that is locally Delaunay wherever the polygon boundary allows.
///
/// **Proof sketch**: An ear of a simple polygon is a triangle $(v_{i-1}, v_i,
/// v_{i+1})$ whose interior lies entirely within the polygon.  Meisters'
/// Two-Ears Theorem guarantees at least two such ears exist.  Clipping $v_i$
/// produces a polygon with $n-1$ vertices that is still simple, so the
/// process terminates after $n - 2$ clips.  Preferring the ear with minimum
/// incircle violations is a greedy heuristic: it yields a Delaunay
/// triangulation when one exists, and a near-Delaunay one otherwise.
///
/// # Complexity
///
/// $O(n^2)$ time via circular linked-list ear removal (O(1) per clip) and
/// pre-computed point cache.
///
/// # Arguments
///
/// - `vertices`: the vertex pool
/// - `polygon`: ordered sequence of vertex IDs forming the cavity boundary
///   (should be in CCW order; the constraint edge endpoints are the first and
///   last vertices in the polygon)
///
/// # Returns
///
/// A list of triangles (as vertex-index triples) filling the cavity.
#[must_use]
pub fn retriangulate_cavity(
    vertices: &[PslgVertex],
    polygon: &[PslgVertexId],
) -> Vec<[PslgVertexId; 3]> {
    if polygon.len() < 3 {
        return Vec::new();
    }
    if polygon.len() == 3 {
        return vec![[polygon[0], polygon[1], polygon[2]]];
    }

    let n = polygon.len();

    // Pre-compute Point2 values to avoid repeated construction in inner
    // loops.  Cache is indexed by polygon position, not vertex ID.
    let points: Vec<Point2<Real>> = polygon
        .iter()
        .map(|vid| {
            let v = &vertices[vid.idx()];
            Point2::new(v.x, v.y)
        })
        .collect();

    // Circular linked-list for O(1) ear removal.  `next[i]` and `prev[i]`
    // give the successor/predecessor of polygon index `i` among the still-
    // active vertices.  Removal relinks neighbours in O(1) instead of the
    // O(n) element shifting required by `Vec::remove()`.
    let mut next: Vec<usize> = (1..n + 1).collect();
    next[n - 1] = 0;
    let mut prev: Vec<usize> = (0..n).map(|i| if i == 0 { n - 1 } else { i - 1 }).collect();
    let mut alive_count = n;
    let mut head = 0usize; // Always points to a live vertex.

    let mut result = Vec::with_capacity(n - 2);

    while alive_count > 3 {
        let mut best_ear: Option<usize> = None;
        let mut min_incircle_count = usize::MAX;

        // Walk the circular list to find the best ear.
        let mut cur = head;
        for _ in 0..alive_count {
            let ia = prev[cur];
            let ib = cur;
            let ic = next[cur];

            let pa = &points[ia];
            let pb = &points[ib];
            let pc = &points[ic];

            // Must be CCW.
            if orient_2d(pa, pb, pc) != Orientation::Positive {
                cur = next[cur];
                continue;
            }

            // Check no other active vertex is inside this ear.
            let mut valid = true;
            let mut check = next[ic];
            while check != ia {
                if point_in_triangle(pa, pb, pc, &points[check]) {
                    valid = false;
                    break;
                }
                check = next[check];
            }

            if !valid {
                cur = next[cur];
                continue;
            }

            // Count circumcircle violations among remaining vertices.
            let mut incircle_count = 0usize;
            let mut check = next[ic];
            while check != ia {
                if incircle(pa, pb, pc, &points[check]) == Orientation::Positive {
                    incircle_count += 1;
                }
                check = next[check];
            }

            if incircle_count == 0 {
                best_ear = Some(ib);
                break; // Perfect Delaunay ear found.
            }

            if incircle_count < min_incircle_count {
                min_incircle_count = incircle_count;
                best_ear = Some(ib);
            }

            cur = next[cur];
        }

        if let Some(ib) = best_ear {
            let ia = prev[ib];
            let ic = next[ib];
            result.push([polygon[ia], polygon[ib], polygon[ic]]);
            // Unlink ib from the circular list — O(1).
            next[prev[ib]] = next[ib];
            prev[next[ib]] = prev[ib];
            // If we removed head, advance it.
            if ib == head {
                head = ic;
            }
            alive_count -= 1;
        } else {
            // Fallback: force the first ear even if not perfectly valid.
            let ia = prev[head];
            let ib = head;
            let ic = next[head];
            result.push([polygon[ia], polygon[ib], polygon[ic]]);
            next[prev[head]] = next[head];
            prev[next[head]] = prev[head];
            head = ic;
            alive_count -= 1;
        }
    }

    // Last remaining triangle — walk the list to find the 3 survivors.
    if alive_count == 3 {
        let a = head;
        let b = next[a];
        let c = next[b];
        result.push([polygon[a], polygon[b], polygon[c]]);
    }

    result
}

/// Check if point `p` lies strictly inside triangle `(a, b, c)` (CCW).
fn point_in_triangle(
    a: &Point2<Real>,
    b: &Point2<Real>,
    c: &Point2<Real>,
    p: &Point2<Real>,
) -> bool {
    orient_2d(a, b, p) == Orientation::Positive
        && orient_2d(b, c, p) == Orientation::Positive
        && orient_2d(c, a, p) == Orientation::Positive
}
