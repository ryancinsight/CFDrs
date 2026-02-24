//! Triangle-triangle adjacency bookkeeping.
//!
//! Provides helper functions for maintaining the symmetric adjacency invariant:
//!
//! # Invariant
//!
//! For every triangle `t` and edge index `i`:
//! - If `t.adj[i] == GHOST_TRIANGLE`, the edge is on the convex hull boundary.
//! - Otherwise, `t.adj[i]` is the triangle sharing that edge, and there exists
//!   an edge index `j` in `t.adj[i]` such that `t.adj[i].adj[j] == t`.
//!
//! This symmetry mirrors the twin-involution property of the half-edge mesh.

use super::triangle::{Triangle, TriangleId, GHOST_TRIANGLE};
use crate::application::delaunay::pslg::vertex::PslgVertexId;

/// Adjacency helper for the triangle pool.
///
/// Wraps a `&mut Vec<Triangle>` and exposes adjacency maintenance operations.
pub struct Adjacency;

impl Adjacency {
    /// Set the adjacency between two triangles across a shared edge.
    ///
    /// `t1.adj[e1] = t2` and `t2.adj[e2] = t1`.
    #[inline]
    pub fn link(triangles: &mut [Triangle], t1: TriangleId, e1: usize, t2: TriangleId, e2: usize) {
        triangles[t1.idx()].adj[e1] = t2;
        triangles[t2.idx()].adj[e2] = t1;
    }

    /// Set `t1.adj[e1] = t2` (one-sided).  Used when `t2` is `GHOST_TRIANGLE`.
    #[inline]
    pub fn link_one(triangles: &mut [Triangle], t1: TriangleId, e1: usize, t2: TriangleId) {
        triangles[t1.idx()].adj[e1] = t2;
    }

    /// Find which edge index of triangle `t` is adjacent to `other`.
    ///
    /// Returns `None` if `other` is not adjacent.
    #[inline]
    pub fn find_edge(triangles: &[Triangle], t: TriangleId, other: TriangleId) -> Option<usize> {
        triangles[t.idx()].adj.iter().position(|&a| a == other)
    }

    /// Find the edge index in triangle `t` that is opposite vertex `v`.
    #[inline]
    pub fn edge_opposite_vertex(
        triangles: &[Triangle],
        t: TriangleId,
        v: PslgVertexId,
    ) -> Option<usize> {
        triangles[t.idx()].opposite_edge(v)
    }

    /// Return the triangle across edge `e` of triangle `t`.
    #[inline]
    pub fn neighbor(triangles: &[Triangle], t: TriangleId, e: usize) -> TriangleId {
        triangles[t.idx()].adj[e]
    }

    /// Check whether triangle `t` is on the convex hull boundary
    /// (has at least one `GHOST_TRIANGLE` neighbor).
    pub fn is_hull_triangle(triangles: &[Triangle], t: TriangleId) -> bool {
        triangles[t.idx()].adj.iter().any(|&a| a == GHOST_TRIANGLE)
    }

    /// Verify the symmetric adjacency invariant for all alive triangles.
    ///
    /// Returns `true` if every adjacency link is symmetric.
    #[cfg(any(test, debug_assertions))]
    pub fn verify_symmetry(triangles: &[Triangle]) -> bool {
        for (i, tri) in triangles.iter().enumerate() {
            if !tri.alive {
                continue;
            }
            let tid = TriangleId::from_usize(i);
            for e in 0..3 {
                let nbr = tri.adj[e];
                if nbr == GHOST_TRIANGLE {
                    continue;
                }
                if nbr.idx() >= triangles.len() || !triangles[nbr.idx()].alive {
                    return false;
                }
                if Self::find_edge(triangles, nbr, tid).is_none() {
                    return false;
                }
            }
        }
        true
    }
}
