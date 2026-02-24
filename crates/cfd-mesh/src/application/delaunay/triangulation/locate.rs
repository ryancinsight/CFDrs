//! Point location via Lawson's oriented walk.
//!
//! # Theorem — Lawson Walk Correctness
//!
//! **Statement**: Starting from any triangle in a Delaunay triangulation and
//! walking toward the query point by crossing edges where `orient_2d` indicates
//! the point lies on the opposite side, the walk reaches the containing
//! triangle in $O(\sqrt{n})$ expected steps for uniformly distributed points.
//!
//! **Proof sketch**: At each step, the signed area of the walk triangle with
//! respect to the query point strictly decreases.  Since the triangulation has
//! a finite number of triangles and the walk never revisits a triangle, it
//! terminates.  The $O(\sqrt{n})$ bound follows from a geometric argument
//! showing that the walk crosses at most $O(\sqrt{n})$ Delaunay edges on
//! average for random point distributions (Devroye, Lemaire, Moreau 2004).
//!
//! # Algorithm
//!
//! ```text
//! LOCATE(q, start_triangle):
//!   t ← start_triangle
//!   LOOP:
//!     for each edge i of t:
//!       let (a, b) = edge vertices of t.edge[i]
//!       if orient_2d(a, b, q) < 0:  // q is on the wrong side
//!         t ← t.adj[i]
//!         CONTINUE LOOP
//!     RETURN t  // q is inside t (or on its boundary)
//! ```

use crate::domain::core::scalar::Real;
use crate::domain::geometry::predicates::{orient_2d, Orientation};
use nalgebra::Point2;

use super::triangle::{TriangleId, GHOST_TRIANGLE};
use crate::application::delaunay::pslg::vertex::PslgVertex;

/// Result of a point location query.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Location {
    /// The query point lies strictly inside triangle `tid`.
    Inside(TriangleId),
    /// The query point lies on edge `edge` of triangle `tid`.
    OnEdge(TriangleId, usize),
    /// The query point coincides with vertex at local index `vi` of triangle `tid`.
    OnVertex(TriangleId, usize),
}

/// Locate the triangle containing point `(qx, qy)` using Lawson's oriented walk.
///
/// # Arguments
///
/// - `vertices`: the vertex pool
/// - `triangles`: the triangle pool (may contain tombstoned entries)
/// - `start`: the triangle to begin the walk from
/// - `qx`, `qy`: query point coordinates
///
/// # Returns
///
/// A [`Location`] describing where the query point lies relative to the
/// triangulation.  Returns `None` if the point is outside the convex hull
/// (walked off the boundary).
pub fn locate(
    vertices: &[PslgVertex],
    triangles: &[super::triangle::Triangle],
    start: TriangleId,
    qx: Real,
    qy: Real,
) -> Option<Location> {
    let q = Point2::new(qx, qy);
    let mut tid = start;
    let max_steps = triangles.len() * 2;

    for _ in 0..max_steps {
        if tid == GHOST_TRIANGLE {
            return None;
        }
        let tri = &triangles[tid.idx()];
        if !tri.alive {
            return None;
        }

        let v0 = &vertices[tri.vertices[0].idx()];
        let v1 = &vertices[tri.vertices[1].idx()];
        let v2 = &vertices[tri.vertices[2].idx()];

        let p0 = Point2::new(v0.x, v0.y);
        let p1 = Point2::new(v1.x, v1.y);
        let p2 = Point2::new(v2.x, v2.y);

        // Check each edge (opposite vertex i = edge from v[(i+1)%3] to v[(i+2)%3])
        let o0 = orient_2d(&p1, &p2, &q); // edge opposite v0
        if o0 == Orientation::Negative {
            tid = tri.adj[0];
            continue;
        }

        let o1 = orient_2d(&p2, &p0, &q); // edge opposite v1
        if o1 == Orientation::Negative {
            tid = tri.adj[1];
            continue;
        }

        let o2 = orient_2d(&p0, &p1, &q); // edge opposite v2
        if o2 == Orientation::Negative {
            tid = tri.adj[2];
            continue;
        }

        // q is inside or on the boundary of this triangle.
        // Check for vertex-coincidence first.
        if o1 == Orientation::Degenerate && o2 == Orientation::Degenerate {
            return Some(Location::OnVertex(tid, 0));
        }
        if o0 == Orientation::Degenerate && o2 == Orientation::Degenerate {
            return Some(Location::OnVertex(tid, 1));
        }
        if o0 == Orientation::Degenerate && o1 == Orientation::Degenerate {
            return Some(Location::OnVertex(tid, 2));
        }

        // Check for edge-coincidence.
        if o0 == Orientation::Degenerate {
            return Some(Location::OnEdge(tid, 0));
        }
        if o1 == Orientation::Degenerate {
            return Some(Location::OnEdge(tid, 1));
        }
        if o2 == Orientation::Degenerate {
            return Some(Location::OnEdge(tid, 2));
        }

        return Some(Location::Inside(tid));
    }

    None
}
