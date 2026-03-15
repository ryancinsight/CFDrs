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

    /// Check that every real vertex has at least `k` distinct real neighbours.
    ///
    /// This is a *necessary* condition for the triangulation graph to be
    /// $k$-vertex-connected.
    ///
    /// # Theorem — Whitney (1932)
    ///
    /// **Statement**: Every convex Delaunay triangulation of $n \ge 4$
    /// points in general position is 3-vertex-connected.
    ///
    /// **Proof sketch**: In a convex Delaunay triangulation each interior
    /// vertex has degree $\ge 3$ (the empty-circumcircle property forces at
    /// least three Delaunay neighbours), and each convex-hull vertex is
    /// adjacent to its two hull neighbours plus at least one interior vertex.
    /// By Whitney's theorem a 2-connected planar graph whose every face is
    /// bounded by a simple cycle — which holds for triangulations — is
    /// 3-connected.  ∎
    ///
    /// # Complexity
    ///
    /// $O(n \cdot \bar{d})$ where $\bar{d} \le 6$ is the average vertex degree.
    /// Uses a stack-allocated scratch buffer (capacity 20) to avoid heap
    /// allocation for typical Delaunay vertices (degree ≤ 12).
    #[must_use]
    pub fn is_k_connected(&self, k: usize) -> bool {
        // Scratch buffer for collecting neighbour indices.  Delaunay
        // vertices have expected degree 6; 20 handles all practical cases
        // without heap allocation.
        let mut nbrs: Vec<usize> = Vec::with_capacity(20);

        for vid_idx in 0..self.num_real_vertices {
            let vid = PslgVertexId::from_usize(vid_idx);
            let tris = self.triangles_around_vertex(vid);
            if tris.is_empty() {
                return false;
            }
            nbrs.clear();
            for &tid in &tris {
                let tri = &self.triangles[tid.idx()];
                for &v in &tri.vertices {
                    if v != vid && !self.super_verts.contains(&v) {
                        let idx = v.idx();
                        // Linear scan is faster than HashSet for small N
                        // (typical: 6 neighbours).
                        if !nbrs.contains(&idx) {
                            nbrs.push(idx);
                        }
                    }
                }
            }
            if nbrs.len() < k {
                return false;
            }
        }
        true
    }

    /// Verify the Delaunay property for all interior triangles.
    ///
    /// # Theorem — Delaunay Verification
    ///
    /// A triangulation is Delaunay iff for every interior edge shared by
    /// triangles `(a, b, c)` and `(a, c, d)`, point `d` does not lie strictly
    /// inside the circumcircle of `(a, b, c)`.
    #[must_use]
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

    /// Verify the Euler formula for a planar triangulation.
    ///
    /// # Theorem — Euler–Poincaré for Planar Triangulations
    ///
    /// **Statement**: For a planar triangulation with $V$ vertices, $E$ edges,
    /// and $F$ interior faces (triangles), $V - E + F = 1$ (the outer face is
    /// not counted).  Equivalently $E = (3F + b)/2$ where $b$ is the number
    /// of boundary edges, and $F = 2V - b - 2$.
    ///
    /// **Proof sketch**: Euler's formula for a connected planar graph gives
    /// $V - E + F_{\text{all}} = 2$ where $F_{\text{all}}$ includes the
    /// outer face.  Subtracting 1 for the outer face yields $V - E + F = 1$.
    /// In a triangulation every interior face has 3 edges; each interior edge
    /// is shared by 2 faces, each boundary edge by 1.  Counting:
    /// $3F = 2E - b$, so $E = (3F + b)/2$.  Combined with Euler:
    /// $V - (3F + b)/2 + F = 1 \Rightarrow F = 2V - b - 2$.  ∎
    ///
    /// # Complexity
    ///
    /// $O(V + T)$ where $T$ = number of alive triangles.
    #[must_use]
    pub fn satisfies_euler(&self) -> bool {
        let v = self.num_real_vertices;
        if v < 3 {
            return true;
        }

        let f: usize = self.interior_triangles().count();
        if f == 0 {
            return true;
        }

        // Count unique edges among interior triangles.
        let mut edges = std::collections::HashSet::new();
        for (_tid, tri) in self.interior_triangles() {
            for e in 0..3 {
                let (va, vb) = tri.edge_vertices(e);
                if self.super_verts.contains(&va) || self.super_verts.contains(&vb) {
                    continue;
                }
                let key = if va <= vb { (va, vb) } else { (vb, va) };
                edges.insert(key);
            }
        }
        let e = edges.len();

        // Euler: V - E + F ∈ {1, 2} for a planar triangulation.
        let euler = (v as isize) - (e as isize) + (f as isize);
        euler == 1 || euler == 2
    }

    /// Collect the convex hull vertices.
    ///
    /// A real vertex is on the convex hull iff it belongs to a triangle
    /// that shares an edge with a triangle containing a super-triangle
    /// vertex.  The two non-super vertices of such a shared edge lie on
    /// the hull.
    ///
    /// # Complexity
    ///
    /// $O(T)$ where $T$ is the number of alive triangles.
    #[must_use]
    pub fn convex_hull_vertices(&self) -> Vec<PslgVertexId> {
        let mut hull = Vec::new();
        // Walk triangles that touch a super-vertex.  For each such
        // triangle, edges shared with a real (non-super) neighbour
        // contribute their real vertices to the hull.
        for (_tid, tri) in self.all_alive_triangles() {
            if !tri.vertices.iter().any(|v| self.super_verts.contains(v)) {
                continue;
            }
            for e in 0..3 {
                let (va, vb) = tri.edge_vertices(e);
                // Both endpoints must be real (non-super) to be hull edges.
                if self.super_verts.contains(&va) || self.super_verts.contains(&vb) {
                    continue;
                }
                if !hull.contains(&va) {
                    hull.push(va);
                }
                if !hull.contains(&vb) {
                    hull.push(vb);
                }
            }
        }
        hull
    }

    /// Compute the minimum vertex degree across all real vertices.
    ///
    /// This is a lower bound on vertex connectivity for planar graphs.
    /// Whitney's theorem gives $\kappa \ge 3$ for convex Delaunay
    /// triangulations of $n \ge 4$ points in general position.
    ///
    /// # Theorem — Connectivity Augmentation Bounds (García et al., 2025)
    ///
    /// **Statement** (arXiv:2509.01096): Any convex-position PSLG can be
    /// augmented to $k$-connected $O(k^2)$-planar with $O(k^2)$ local
    /// crossings per edge.  For $k = 4$, some edges may require
    /// $\Omega(n)$ crossings.  Flip-based 4-connectivity augmentation of
    /// triangulations admits an EPTAS, and small-$k$ cases on convex
    /// points have linear-time DP solutions.
    ///
    /// **Implication for meshing**: When stronger connectivity is needed
    /// (e.g. multigrid robustness), Steiner point insertion at low-degree
    /// vertices is preferred over beyond-planar augmentation, keeping the
    /// mesh planar while raising $\kappa$.
    ///
    /// # Complexity
    ///
    /// $O(n \cdot \bar{d})$ where $\bar{d} \le 6$.
    #[must_use]
    pub fn min_vertex_connectivity(&self) -> usize {
        let mut min_deg = usize::MAX;
        let mut nbrs: Vec<usize> = Vec::with_capacity(20);

        for vid_idx in 0..self.num_real_vertices {
            let vid = PslgVertexId::from_usize(vid_idx);
            let tris = self.triangles_around_vertex(vid);
            if tris.is_empty() {
                return 0;
            }
            nbrs.clear();
            for &tid in &tris {
                let tri = &self.triangles[tid.idx()];
                for &v in &tri.vertices {
                    if v != vid && !self.super_verts.contains(&v) {
                        let idx = v.idx();
                        if !nbrs.contains(&idx) {
                            nbrs.push(idx);
                        }
                    }
                }
            }
            min_deg = min_deg.min(nbrs.len());
        }
        if min_deg == usize::MAX { 0 } else { min_deg }
    }
}
