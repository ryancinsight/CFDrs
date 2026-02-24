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

    // Ear-clipping with Delaunay preference.
    // For each candidate ear (a, b, c) where b is the ear tip:
    //   1. Check that (a, b, c) is CCW.
    //   2. Check that no other polygon vertex lies inside the triangle.
    //   3. Among valid ears, prefer the one whose circumcircle is most empty.
    let mut remaining: Vec<usize> = (0..polygon.len()).collect();
    let mut result = Vec::with_capacity(polygon.len() - 2);

    while remaining.len() > 3 {
        let n = remaining.len();
        let mut best_ear: Option<usize> = None;
        let mut min_incircle_count = usize::MAX;

        for i in 0..n {
            let ia = remaining[(i + n - 1) % n];
            let ib = remaining[i];
            let ic = remaining[(i + 1) % n];

            let a = &vertices[polygon[ia].idx()];
            let b = &vertices[polygon[ib].idx()];
            let c = &vertices[polygon[ic].idx()];

            let pa = Point2::new(a.x, a.y);
            let pb = Point2::new(b.x, b.y);
            let pc = Point2::new(c.x, c.y);

            // Must be CCW.
            if orient_2d(&pa, &pb, &pc) != Orientation::Positive {
                continue;
            }

            // Check no other vertex is inside this ear.
            let mut valid = true;
            for j in 0..n {
                if j == (i + n - 1) % n || j == i || j == (i + 1) % n {
                    continue;
                }
                let d = &vertices[polygon[remaining[j]].idx()];
                let pd = Point2::new(d.x, d.y);
                if point_in_triangle(&pa, &pb, &pc, &pd) {
                    valid = false;
                    break;
                }
            }

            if !valid {
                continue;
            }

            // Count how many vertices are strictly inside the circumcircle
            let mut incircle_count = 0;
            for j in 0..n {
                if j == (i + n - 1) % n || j == i || j == (i + 1) % n {
                    continue;
                }
                let d = &vertices[polygon[remaining[j]].idx()];
                let pd = Point2::new(d.x, d.y);
                if incircle(&pa, &pb, &pc, &pd) == Orientation::Positive {
                    incircle_count += 1;
                }
            }

            if incircle_count == 0 {
                best_ear = Some(i);
                break; // Perfect Delaunay ear found
            }

            if incircle_count < min_incircle_count {
                min_incircle_count = incircle_count;
                best_ear = Some(i);
            }
        }

        match best_ear {
            Some(i) => {
                let ia = remaining[(i + remaining.len() - 1) % remaining.len()];
                let ib = remaining[i];
                let ic = remaining[(i + 1) % remaining.len()];
                result.push([polygon[ia], polygon[ib], polygon[ic]]);
                remaining.remove(i);
            }
            None => {
                // Fallback: force the first ear even if not perfectly valid.
                let ia = remaining[0];
                let ib = remaining[1];
                let ic = remaining[2];
                result.push([polygon[ia], polygon[ib], polygon[ic]]);
                remaining.remove(1);
            }
        }
    }

    // Last remaining triangle.
    if remaining.len() == 3 {
        result.push([
            polygon[remaining[0]],
            polygon[remaining[1]],
            polygon[remaining[2]],
        ]);
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
