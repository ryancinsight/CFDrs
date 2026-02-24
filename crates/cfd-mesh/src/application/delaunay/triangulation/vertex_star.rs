//! Vertex-star traversal for the Delaunay triangulation.
//!
//! Provides efficient O(deg(v)) iteration over all triangles incident on a
//! given vertex by fan-walking through the triangle adjacency structure.
//!
//! # Theorem — Vertex-Star Degree
//!
//! **Statement**: In a Delaunay triangulation of $n$ well-distributed points,
//! the expected degree of any vertex is $\le 6$ (exactly 6 for interior
//! vertices of a triangulation of a convex point set by Euler's formula).
//!
//! **Proof sketch**: By the Euler formula $V - E + F = 2$ and the relation
//! $2E = 3F$ for a triangulation, the average vertex degree is
//! $\bar{d} = 2E/V = 6 - 12/V \to 6$ as $V \to \infty$.

use super::bowyer_watson::DelaunayTriangulation;
use super::triangle::{TriangleId, GHOST_TRIANGLE};
use crate::application::delaunay::pslg::vertex::PslgVertexId;
use crate::domain::geometry::predicates::{incircle, orient_2d, Orientation};

impl DelaunayTriangulation {
    /// Collect all alive triangles incident on vertex `v` by fan-walking
    /// around the vertex star using adjacency links.
    ///
    /// Returns an empty `Vec` if the vertex has no recorded incident triangle
    /// (e.g. before insertion).
    ///
    /// # Complexity
    ///
    /// O(deg(v)), typically ≈ 6 for a well-distributed Delaunay triangulation.
    pub(crate) fn triangles_around_vertex(&self, v: PslgVertexId) -> Vec<TriangleId> {
        let start = self.vert_to_tri[v.idx()];
        if start == GHOST_TRIANGLE || !self.triangles[start.idx()].alive {
            return self.triangles_around_vertex_linear(v);
        }
        let mut result = Vec::with_capacity(8);
        let mut cur = start;
        loop {
            if !self.triangles[cur.idx()].alive {
                return self.triangles_around_vertex_linear(v);
            }
            result.push(cur);
            let tri = &self.triangles[cur.idx()];
            let li = match tri.vertex_index(v) {
                Some(i) => i,
                None => return self.triangles_around_vertex_linear(v),
            };
            let next = tri.adj[(li + 1) % 3];
            if next == GHOST_TRIANGLE || next == start {
                break;
            }
            cur = next;
        }
        // If we stopped at GHOST (hull vertex), walk the other direction.
        if self.triangles[start.idx()].alive {
            let tri = &self.triangles[start.idx()];
            if let Some(li) = tri.vertex_index(v) {
                let mut cur2 = tri.adj[(li + 2) % 3];
                while cur2 != GHOST_TRIANGLE && cur2 != start {
                    if !self.triangles[cur2.idx()].alive {
                        break;
                    }
                    result.push(cur2);
                    let t2 = &self.triangles[cur2.idx()];
                    let li2 = match t2.vertex_index(v) {
                        Some(i) => i,
                        None => break,
                    };
                    cur2 = t2.adj[(li2 + 2) % 3];
                }
            }
        }
        result
    }

    /// Linear-scan fallback when the `vert_to_tri` hint is stale.
    fn triangles_around_vertex_linear(&self, v: PslgVertexId) -> Vec<TriangleId> {
        self.triangles
            .iter()
            .enumerate()
            .filter(|(_, t)| t.alive && t.contains_vertex(v))
            .map(|(i, _)| TriangleId::from_usize(i))
            .collect()
    }

    /// Check whether edge (a, b) exists in the triangulation.
    ///
    /// Uses the vertex-star walk: O(deg(a)) ≈ O(6) instead of O(T).
    pub(crate) fn edge_exists_fast(&self, a: PslgVertexId, b: PslgVertexId) -> bool {
        for tid in self.triangles_around_vertex(a) {
            let tri = &self.triangles[tid.idx()];
            if tri.contains_vertex(b) {
                return true;
            }
        }
        false
    }

    /// Find the triangle containing edge (a, b) and the local edge index.
    ///
    /// Uses the vertex-star walk: O(deg(a)) instead of O(T).
    pub(crate) fn find_edge_fast(
        &self,
        a: PslgVertexId,
        b: PslgVertexId,
    ) -> Option<(TriangleId, usize)> {
        for tid in self.triangles_around_vertex(a) {
            let tri = &self.triangles[tid.idx()];
            for edge_idx in 0..3 {
                let (ea, eb) = tri.edge_vertices(edge_idx);
                if (ea == a && eb == b) || (ea == b && eb == a) {
                    return Some((tid, edge_idx));
                }
            }
        }
        None
    }

    /// Verify the Delaunay property for all interior triangles.
    ///
    /// # Theorem — Delaunay Verification
    ///
    /// A triangulation is Delaunay iff for every interior edge shared by
    /// triangles `(a, b, c)` and `(a, c, d)`, point `d` does not lie strictly
    /// inside the circumcircle of `(a, b, c)`.
    pub fn is_delaunay(&self) -> bool {
        for (tid, tri) in self.all_alive_triangles() {
            for edge in 0..3 {
                if tri.constrained[edge] {
                    continue;
                }
                let nbr = tri.adj[edge];
                if nbr == GHOST_TRIANGLE {
                    continue;
                }
                let nbr_tri = &self.triangles[nbr.idx()];
                if !nbr_tri.alive {
                    continue;
                }
                let nbr_edge = match nbr_tri.shared_edge(tid) {
                    Some(e) => e,
                    None => return false,
                };
                let v_opp = nbr_tri.vertices[nbr_edge];
                let [a, b, c] = tri.vertices;
                let pa = self.vertices[a.idx()].to_point2();
                let pb = self.vertices[b.idx()].to_point2();
                let pc = self.vertices[c.idx()].to_point2();
                let pd = self.vertices[v_opp.idx()].to_point2();
                let ort = orient_2d(&pa, &pb, &pc);
                let inside = if ort == Orientation::Positive {
                    incircle(&pa, &pb, &pc, &pd) == Orientation::Positive
                } else if ort == Orientation::Negative {
                    incircle(&pb, &pa, &pc, &pd) == Orientation::Positive
                } else {
                    false
                };
                if inside {
                    return false;
                }
            }
        }
        true
    }
}
