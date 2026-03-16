//! Incremental Bowyer-Watson Delaunay triangulation.
//!
//! # Theorem — Bowyer-Watson Correctness
//!
//! **Statement** (Bowyer 1981, Watson 1981): Given $DT(P)$ and a new point
//! $p$, the algorithm identifies and removes the *bad cavity* (triangles
//! whose circumcircle contains $p$), then fan-connects $p$ to the cavity
//! boundary.  The result is $DT(P \cup \{p\})$.
//!
//! **Proof sketch**: By the empty-circumcircle characterisation, exactly
//! those triangles whose circumcircle contains $p$ are invalidated.
//! Fan-connecting $p$ to the cavity boundary creates triangles whose
//! circumcircles contain only $p$ and the two edge endpoints — hence
//! Delaunay.
//!
//! # Complexity
//!
//! Expected $O(n \log n)$ (uniformly distributed); worst-case $O(n^2)$.

use super::locate::Location;
use super::ordering::hilbert_order;
use super::triangle::{Triangle, TriangleId, GHOST_TRIANGLE};
use crate::application::delaunay::pslg::vertex::{PslgVertex, PslgVertexId};
use crate::domain::core::scalar::Real;
use crate::domain::geometry::predicates::{incircle, orient_2d, Orientation};

/// The Delaunay triangulation data structure.
///
/// Stores vertices and triangles as flat arrays with explicit adjacency.
/// Supports incremental point insertion via the Bowyer-Watson algorithm.
///
/// # Vertex-to-Triangle Adjacency
///
/// A `vert_to_tri` map is maintained: for each vertex `v`, `vert_to_tri[v]`
/// contains one incident alive triangle.  This enables O(deg(v)) edge
/// lookups instead of O(T) full scans — a **critical** performance
/// improvement for constraint enforcement where `edge_exists()`,
/// `collect_crossing_edges()`, and `mark_edge_constrained()` are called
/// per constraint segment.
///
/// # Example
///
/// ```rust,ignore
/// use cfd_mesh::application::delaunay::triangulation::DelaunayTriangulation;
///
/// let points = vec![(0.0, 0.0), (1.0, 0.0), (0.5, 1.0), (0.5, 0.3)];
/// let dt = DelaunayTriangulation::from_points(&points);
/// assert!(dt.is_delaunay());
/// ```
pub struct DelaunayTriangulation {
    /// Vertex pool (includes 3 super-triangle vertices at the end).
    pub(crate) vertices: Vec<PslgVertex>,
    /// Triangle pool (tombstoned entries have `alive == false`).
    pub(crate) triangles: Vec<Triangle>,
    /// Index of the last triangle used for point-location (walk hint).
    last_triangle: TriangleId,
    /// Number of real (non-super-triangle) vertices.
    pub(crate) num_real_vertices: usize,
    /// IDs of the three super-triangle vertices.
    pub(crate) super_verts: [PslgVertexId; 3],
    /// Vertex → one incident alive triangle.  Enables O(deg(v)) local
    /// traversal for edge queries instead of O(T) global scans.
    pub(crate) vert_to_tri: Vec<TriangleId>,
    /// Epoch stamps for O(1) amortised cycle detection in the Lawson walk.
    /// `locate_epoch[tid]` holds the generation when triangle `tid` was last
    /// visited.
    locate_epoch: Vec<u32>,
    /// Generation counter, incremented once per `locate_point` call.
    locate_gen: u32,
}

impl DelaunayTriangulation {
    /// Build a Delaunay triangulation from a set of 2-D points.
    ///
    /// The points are inserted incrementally via Bowyer-Watson.
    #[must_use]
    pub fn from_points(points: &[(Real, Real)]) -> Self {
        let vertices: Vec<PslgVertex> =
            points.iter().map(|&(x, y)| PslgVertex::new(x, y)).collect();
        Self::from_vertices(vertices)
    }

    /// Build from an existing vertex array.
    ///
    /// Uses Hilbert-curve ordering (BRIO-inspired) for optimal O(n log n)
    /// insertion performance.  By sorting vertices along a space-filling
    /// curve, consecutive insertions are spatially adjacent, reducing the
    /// Lawson walk distance from O(√n) to near-O(1) per insertion.
    #[must_use]
    pub fn from_vertices(real_vertices: Vec<PslgVertex>) -> Self {
        let n = real_vertices.len();
        let mut dt = Self::init_with_super_triangle(&real_vertices);

        if n == 0 {
            return dt;
        }

        // Compute Hilbert-curve insertion order.
        let order = hilbert_order(&real_vertices);

        for &orig_idx in &order {
            let vid = PslgVertexId::from_usize(orig_idx);
            dt.insert_vertex(vid);
        }

        dt.num_real_vertices = n;
        dt
    }

    /// Initialize the triangulation with a super-triangle that encloses all
    /// input vertices.
    fn init_with_super_triangle(real_vertices: &[PslgVertex]) -> Self {
        // Compute bounding box and inflate it significantly.
        let (mut min_x, mut min_y) = (Real::INFINITY, Real::INFINITY);
        let (mut max_x, mut max_y) = (Real::NEG_INFINITY, Real::NEG_INFINITY);

        for v in real_vertices {
            if v.x < min_x {
                min_x = v.x;
            }
            if v.y < min_y {
                min_y = v.y;
            }
            if v.x > max_x {
                max_x = v.x;
            }
            if v.y > max_y {
                max_y = v.y;
            }
        }

        // Handle degenerate case (0 or 1 point).
        if !min_x.is_finite() {
            min_x = -1.0;
            min_y = -1.0;
            max_x = 1.0;
            max_y = 1.0;
        }

        let dx = (max_x - min_x).max(1e-10);
        let dy = (max_y - min_y).max(1e-10);
        let dmax = dx.max(dy);
        let cx = (min_x + max_x) * 0.5;
        let cy = (min_y + max_y) * 0.5;

        // Super-triangle vertices far outside the bounding box.
        let margin = 20.0 * dmax;
        let sv0 = PslgVertex::new(cx - margin, cy - margin);
        let sv1 = PslgVertex::new(cx + margin, cy - margin);
        let sv2 = PslgVertex::new(cx, cy + margin);

        let n = real_vertices.len();
        let sv0_id = PslgVertexId::from_usize(n);
        let sv1_id = PslgVertexId::from_usize(n + 1);
        let sv2_id = PslgVertexId::from_usize(n + 2);

        let mut vertices: Vec<PslgVertex> = Vec::with_capacity(n + 3);
        vertices.extend_from_slice(real_vertices);
        vertices.push(sv0);
        vertices.push(sv1);
        vertices.push(sv2);

        let super_tri = Triangle::new(sv0_id, sv1_id, sv2_id);
        let triangles = vec![super_tri];
        let super_tid = TriangleId::new(0);

        // Initialize vert_to_tri: super-triangle vertices point to t0,
        // real vertices start as GHOST (updated on insertion).
        let mut vert_to_tri = vec![GHOST_TRIANGLE; n + 3];
        vert_to_tri[sv0_id.idx()] = super_tid;
        vert_to_tri[sv1_id.idx()] = super_tid;
        vert_to_tri[sv2_id.idx()] = super_tid;

        Self {
            vertices,
            triangles,
            last_triangle: TriangleId::new(0),
            num_real_vertices: 0,
            super_verts: [sv0_id, sv1_id, sv2_id],
            vert_to_tri,
            locate_epoch: vec![0; 1],
            locate_gen: 0,
        }
    }

    /// Locate the triangle containing `(qx, qy)` via Lawson's oriented walk
    /// with epoch-based cycle detection.
    ///
    /// # Theorem — Epoch-Stamped Visited Tracking
    ///
    /// **Statement**: A per-triangle epoch vector `E[0..T]` and a global
    /// generation counter `g` provide amortised $O(1)$ visited-set
    /// operations with $O(T)$ memory allocated once, eliminating the
    /// per-call `HashSet` allocation incurred by the cold-path
    /// [`locate()`](super::locate::locate) free function.
    ///
    /// **Proof**: Each call increments `g` and marks visited triangles by
    /// setting `E[tid] = g`.  A triangle is "visited in this call" iff
    /// `E[tid] == g`.  Since `g` is unique per call, all previous marks
    /// are implicitly invalidated without clearing `E`.  The `Vec` grows
    /// monotonically (amortised $O(1)$ push); no other allocation occurs.
    /// When `g` wraps around $2^{32}$ we reset `E` to zero — this happens
    /// once per ~4 billion locate calls and costs $O(T)$.  ∎
    pub(crate) fn locate_point(
        &mut self,
        start: TriangleId,
        qx: Real,
        qy: Real,
    ) -> Option<Location> {
        use nalgebra::Point2;

        self.locate_gen = self.locate_gen.wrapping_add(1);
        if self.locate_gen == 0 {
            // Epoch wrapped — clear all stamps.
            self.locate_epoch.fill(0);
            self.locate_gen = 1;
        }
        let gen = self.locate_gen;

        // Ensure epoch vec covers all triangles.
        if self.locate_epoch.len() < self.triangles.len() {
            self.locate_epoch.resize(self.triangles.len(), 0);
        }

        let q = Point2::new(qx, qy);
        let mut tid = start;
        let max_steps = self.triangles.len() * 3;

        for _ in 0..max_steps {
            if tid == GHOST_TRIANGLE {
                return None;
            }
            let idx = tid.idx();
            if idx >= self.triangles.len() {
                return None;
            }
            if self.locate_epoch[idx] == gen {
                return None; // Cycle detected.
            }
            self.locate_epoch[idx] = gen;

            let tri = &self.triangles[idx];
            if !tri.alive {
                return None;
            }

            let v0 = &self.vertices[tri.vertices[0].idx()];
            let v1 = &self.vertices[tri.vertices[1].idx()];
            let v2 = &self.vertices[tri.vertices[2].idx()];

            let p0 = Point2::new(v0.x, v0.y);
            let p1 = Point2::new(v1.x, v1.y);
            let p2 = Point2::new(v2.x, v2.y);

            let o0 = orient_2d(&p1, &p2, &q);
            if o0 == Orientation::Negative {
                tid = tri.adj[0];
                continue;
            }

            let o1 = orient_2d(&p2, &p0, &q);
            if o1 == Orientation::Negative {
                tid = tri.adj[1];
                continue;
            }

            let o2 = orient_2d(&p0, &p1, &q);
            if o2 == Orientation::Negative {
                tid = tri.adj[2];
                continue;
            }

            // Vertex coincidence.
            if o1 == Orientation::Degenerate && o2 == Orientation::Degenerate {
                return Some(Location::OnVertex(tid, 0));
            }
            if o0 == Orientation::Degenerate && o2 == Orientation::Degenerate {
                return Some(Location::OnVertex(tid, 1));
            }
            if o0 == Orientation::Degenerate && o1 == Orientation::Degenerate {
                return Some(Location::OnVertex(tid, 2));
            }

            // Edge coincidence.
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

    /// Insert a vertex into the triangulation via Bowyer-Watson.
    pub(crate) fn insert_vertex(&mut self, vid: PslgVertexId) {
        let v = self.vertices[vid.idx()];

        // 1. Locate the containing triangle (epoch-based, zero-alloc).
        let loc = if let Some(l) = self.locate_point(self.last_triangle, v.x, v.y) {
            l
        } else {
            // Fallback: linear scan for a valid starting triangle.
            let start = self
                .triangles
                .iter()
                .position(|t| t.alive)
                .map_or(TriangleId::new(0), TriangleId::from_usize);

            match self.locate_point(start, v.x, v.y) {
                Some(l) => l,
                None => return, // Point outside super-triangle — skip.
            }
        };

        match loc {
            Location::OnVertex(_, _) => {
                // Duplicate point — skip insertion.
            }
            Location::Inside(tid) => {
                self.insert_in_triangle(vid, tid);
            }
            Location::OnEdge(tid, edge) => {
                self.insert_on_edge(vid, tid, edge);
            }
        }
    }

    /// Insert vertex `vid` strictly inside triangle `tid`.
    ///
    /// Splits `tid` into 3 sub-triangles and restores Delaunay via edge flips.
    fn insert_in_triangle(&mut self, vid: PslgVertexId, tid: TriangleId) {
        let old = self.triangles[tid.idx()].clone();
        let [v0, v1, v2] = old.vertices;
        let [a0, a1, a2] = old.adj;

        // Kill the old triangle.
        self.triangles[tid.idx()].alive = false;

        // Create 3 new triangles.
        //   t0: (vid, v0, v1)  — replaces tid's slot is reused for t0
        //   t1: (vid, v1, v2)
        //   t2: (vid, v2, v0)
        let t0 = self.push_triangle(Triangle::new(vid, v0, v1));
        let t1 = self.push_triangle(Triangle::new(vid, v1, v2));
        let t2 = self.push_triangle(Triangle::new(vid, v2, v0));

        // Internal adjacency: each new triangle is adjacent to the other two.
        //  t0: edge opposite vid (edge 0) is between v0,v1 → adj[0] = old adj[2]
        //      edge opposite v0 (edge 1) is vid,v1 → adj[1] = t1
        //      edge opposite v1 (edge 2) is vid,v0 → adj[2] = t2
        self.triangles[t0.idx()].adj = [a2, t1, t2];
        self.triangles[t1.idx()].adj = [a0, t2, t0];
        self.triangles[t2.idx()].adj = [a1, t0, t1];

        // Fix external adjacency: the old neighbors pointed to `tid` — update them.
        self.fix_adjacency(a2, tid, t0);
        self.fix_adjacency(a0, tid, t1);
        self.fix_adjacency(a1, tid, t2);

        // Update vert_to_tri for all affected vertices.
        self.vert_to_tri[vid.idx()] = t0;
        self.vert_to_tri[v0.idx()] = t0;
        self.vert_to_tri[v1.idx()] = t1;
        self.vert_to_tri[v2.idx()] = t2;

        self.last_triangle = t0;

        // Restore Delaunay property via edge flips.
        self.flip_fix(t0, 0);
        self.flip_fix(t1, 0);
        self.flip_fix(t2, 0);
    }

    /// Insert vertex `vid` on edge `edge` of triangle `tid`.
    ///
    /// Splits the two triangles sharing the edge into 4 sub-triangles.
    fn insert_on_edge(&mut self, vid: PslgVertexId, tid: TriangleId, edge: usize) {
        let nbr_tid = self.triangles[tid.idx()].adj[edge];

        if nbr_tid == GHOST_TRIANGLE {
            // Edge is on the convex hull — only split the single triangle.
            self.insert_on_hull_edge(vid, tid, edge);
            return;
        }

        let old_t = self.triangles[tid.idx()].clone();
        let old_n = self.triangles[nbr_tid.idx()].clone();

        // Vertices of tid: v_opp is opposite the shared edge.
        let v_opp_t = old_t.vertices[edge];
        let (va, vb) = old_t.edge_vertices(edge); // shared edge vertices

        // Find the opposite vertex in the neighbor.
        let nbr_edge = old_n.shared_edge(tid).expect("adjacency broken");
        let v_opp_n = old_n.vertices[nbr_edge];

        // Kill old triangles.
        self.triangles[tid.idx()].alive = false;
        self.triangles[nbr_tid.idx()].alive = false;

        // Create 4 new triangles:
        //   t0: (vid, v_opp_t, va)
        //   t1: (vid, vb, v_opp_t)
        //   t2: (vid, va, v_opp_n)
        //   t3: (vid, v_opp_n, vb)
        let t0 = self.push_triangle(Triangle::new(vid, v_opp_t, va));
        let t1 = self.push_triangle(Triangle::new(vid, vb, v_opp_t));
        let t2 = self.push_triangle(Triangle::new(vid, va, v_opp_n));
        let t3 = self.push_triangle(Triangle::new(vid, v_opp_n, vb));

        // Adjacency of old_t edges (excluding the shared edge):
        let adj_t_va = old_t.adj[(edge + 2) % 3]; // edge opposite vb in old_t = v_opp_t → va
        let adj_t_vb = old_t.adj[(edge + 1) % 3]; // edge opposite va in old_t = vb → v_opp_t

        let adj_n_va = old_n.adj[(nbr_edge + 1) % 3]; // edge opposite vb in old_n
        let adj_n_vb = old_n.adj[(nbr_edge + 2) % 3]; // edge opposite va in old_n

        // t0: (vid, v_opp_t, va)
        //  edge 0 (opp vid) = v_opp_t→va → adj_t_va
        //  edge 1 (opp v_opp_t) = vid→va → t2
        //  edge 2 (opp va) = vid→v_opp_t → t1
        self.triangles[t0.idx()].adj = [adj_t_va, t2, t1];

        // t1: (vid, vb, v_opp_t)
        //  edge 0 (opp vid) = vb→v_opp_t → adj_t_vb
        //  edge 1 (opp vb) = vid→v_opp_t → t0
        //  edge 2 (opp v_opp_t) = vid→vb → t3
        self.triangles[t1.idx()].adj = [adj_t_vb, t0, t3];

        // t2: (vid, va, v_opp_n)
        //  edge 0 (opp vid) = va→v_opp_n → adj_n_va
        //  edge 1 (opp va) = vid→v_opp_n → t3
        //  edge 2 (opp v_opp_n) = vid→va → t0
        self.triangles[t2.idx()].adj = [adj_n_va, t3, t0];

        // t3: (vid, v_opp_n, vb)
        //  edge 0 (opp vid) = v_opp_n→vb → adj_n_vb
        //  edge 1 (opp v_opp_n) = vid→vb → t1
        //  edge 2 (opp vb) = vid→v_opp_n → t2
        self.triangles[t3.idx()].adj = [adj_n_vb, t1, t2];

        // Fix external neighbors.
        self.fix_adjacency(adj_t_va, tid, t0);
        self.fix_adjacency(adj_t_vb, tid, t1);
        self.fix_adjacency(adj_n_va, nbr_tid, t2);
        self.fix_adjacency(adj_n_vb, nbr_tid, t3);

        // Update vert_to_tri for all affected vertices.
        self.vert_to_tri[vid.idx()] = t0;
        self.vert_to_tri[v_opp_t.idx()] = t0;
        self.vert_to_tri[va.idx()] = t2;
        self.vert_to_tri[vb.idx()] = t1;
        self.vert_to_tri[v_opp_n.idx()] = t2;

        self.last_triangle = t0;

        // Restore Delaunay.
        self.flip_fix(t0, 0);
        self.flip_fix(t1, 0);
        self.flip_fix(t2, 0);
        self.flip_fix(t3, 0);
    }

    /// Handle the special case where the edge is on the hull boundary.
    fn insert_on_hull_edge(&mut self, vid: PslgVertexId, tid: TriangleId, edge: usize) {
        let old = self.triangles[tid.idx()].clone();
        let v_opp = old.vertices[edge];
        let (va, vb) = old.edge_vertices(edge);

        self.triangles[tid.idx()].alive = false;

        // Split into 2 triangles.
        let t0 = self.push_triangle(Triangle::new(vid, v_opp, va));
        let t1 = self.push_triangle(Triangle::new(vid, vb, v_opp));

        let adj_va = old.adj[(edge + 2) % 3];
        let adj_vb = old.adj[(edge + 1) % 3];

        self.triangles[t0.idx()].adj = [adj_va, GHOST_TRIANGLE, t1];
        self.triangles[t1.idx()].adj = [adj_vb, t0, GHOST_TRIANGLE];

        self.fix_adjacency(adj_va, tid, t0);
        self.fix_adjacency(adj_vb, tid, t1);

        // Update vert_to_tri for all affected vertices.
        self.vert_to_tri[vid.idx()] = t0;
        self.vert_to_tri[v_opp.idx()] = t0;
        self.vert_to_tri[va.idx()] = t0;
        self.vert_to_tri[vb.idx()] = t1;

        self.last_triangle = t0;

        self.flip_fix(t0, 0);
        self.flip_fix(t1, 0);
    }

    /// Push a new triangle and return its ID.
    fn push_triangle(&mut self, tri: Triangle) -> TriangleId {
        let id = TriangleId::from_usize(self.triangles.len());
        self.triangles.push(tri);
        id
    }

    /// Fix external adjacency: replace `old_tid` with `new_tid` in the
    /// neighbor's adjacency list.
    fn fix_adjacency(&mut self, nbr: TriangleId, old_tid: TriangleId, new_tid: TriangleId) {
        if nbr == GHOST_TRIANGLE {
            return;
        }
        for a in &mut self.triangles[nbr.idx()].adj {
            if *a == old_tid {
                *a = new_tid;
                return;
            }
        }
    }

    /// Iterative Delaunay edge-flip restoration.
    ///
    /// If the edge `edge` of triangle `tid` violates the Delaunay criterion
    /// (the opposite vertex is inside the circumcircle), flip the edge and
    /// push the two new external edges onto the work stack.
    ///
    /// # Theorem — Flip Termination
    ///
    /// Each flip strictly increases the minimum angle of the affected
    /// quadrilateral.  Since the angle space is bounded and the triangulation
    /// is finite, the flip sequence terminates.
    ///
    /// # Implementation Note
    ///
    /// Uses an explicit stack instead of recursion to avoid stack overflow on
    /// large meshes where deep flip cascades can occur.
    fn flip_fix(&mut self, start_tid: TriangleId, start_edge: usize) {
        let mut stack: Vec<(TriangleId, usize)> = vec![(start_tid, start_edge)];

        while let Some((tid, edge)) = stack.pop() {
            if !self.triangles[tid.idx()].alive {
                continue;
            }

            let nbr = self.triangles[tid.idx()].adj[edge];
            if nbr == GHOST_TRIANGLE {
                continue;
            }

            // Check the constraint flag — never flip a constrained edge.
            if self.triangles[tid.idx()].constrained[edge] {
                continue;
            }

            let tri = &self.triangles[tid.idx()];
            let nbr_tri = &self.triangles[nbr.idx()];

            // Preserve existing constrained flags on non-flipped boundary edges.
            let tri_cons = tri.constrained;
            let nbr_cons = nbr_tri.constrained;

            let v_opp_t = tri.vertices[edge]; // vertex opposite the shared edge in tid
            let nbr_edge = nbr_tri.shared_edge(tid).expect("adjacency broken");
            let v_opp_n = nbr_tri.vertices[nbr_edge]; // vertex opposite in nbr

            // Shared edge vertices.
            let (va, vb) = tri.edge_vertices(edge);

            // In-circle test: is v_opp_n inside circumcircle of (v_opp_t, va, vb)?
            let pa = self.vertices[va.idx()].to_point2();
            let pb = self.vertices[vb.idx()].to_point2();
            let pc = self.vertices[v_opp_t.idx()].to_point2();
            let pd = self.vertices[v_opp_n.idx()].to_point2();

            // orient_2d(a,b,c) must be positive for the incircle test to be correct.
            let ort = orient_2d(&pa, &pb, &pc);
            let inside = if ort == Orientation::Positive {
                incircle(&pa, &pb, &pc, &pd) == Orientation::Positive
            } else if ort == Orientation::Negative {
                incircle(&pb, &pa, &pc, &pd) == Orientation::Positive
            } else {
                // Degenerate triangle — skip.
                continue;
            };

            if !inside {
                continue;
            }

            // Perform the edge flip.
            let adj_tid_va = tri.adj[(edge + 1) % 3];
            let adj_tid_vb = tri.adj[(edge + 2) % 3];
            let adj_nbr_va = nbr_tri.adj[(nbr_edge + 2) % 3];
            let adj_nbr_vb = nbr_tri.adj[(nbr_edge + 1) % 3];

            // Rewrite tid → (v_opp_t, v_opp_n, vb).
            self.triangles[tid.idx()].vertices = [v_opp_t, v_opp_n, vb];
            self.triangles[tid.idx()].adj = [adj_nbr_va, adj_tid_va, nbr];
            self.triangles[tid.idx()].constrained = [
                nbr_cons[(nbr_edge + 2) % 3],
                tri_cons[(edge + 1) % 3],
                false, // New diagonal is never constrained.
            ];

            // Rewrite nbr → (v_opp_n, v_opp_t, va).
            self.triangles[nbr.idx()].vertices = [v_opp_n, v_opp_t, va];
            self.triangles[nbr.idx()].adj = [adj_tid_vb, adj_nbr_vb, tid];
            self.triangles[nbr.idx()].constrained = [
                tri_cons[(edge + 2) % 3],
                nbr_cons[(nbr_edge + 1) % 3],
                false, // New diagonal is never constrained.
            ];

            // Fix external adjacency.
            self.fix_adjacency(adj_nbr_va, nbr, tid);
            self.fix_adjacency(adj_tid_vb, tid, nbr);

            // Update vert_to_tri: vertices may have moved between triangles.
            self.vert_to_tri[v_opp_t.idx()] = tid;
            self.vert_to_tri[v_opp_n.idx()] = nbr;
            self.vert_to_tri[va.idx()] = nbr;
            self.vert_to_tri[vb.idx()] = tid;

            // Push the two new external edges for further checking.
            stack.push((tid, 0)); // (v_opp_n, vb) from old nbr
            stack.push((nbr, 1)); // (va, v_opp_n) from old nbr
        }
    }

    // ── Public query API ──────────────────────────────────────────────────

    /// Number of real (non-super-triangle) vertices.
    #[must_use]
    pub fn vertex_count(&self) -> usize {
        self.num_real_vertices
    }

    /// Number of alive triangles (including those incident to super-triangle vertices).
    #[must_use]
    pub fn triangle_count_raw(&self) -> usize {
        self.triangles.iter().filter(|t| t.alive).count()
    }

    /// Number of interior triangles (excluding those touching super-triangle vertices).
    #[must_use]
    pub fn triangle_count(&self) -> usize {
        self.triangles
            .iter()
            .filter(|t| t.alive && !self.is_super_triangle(t))
            .count()
    }

    /// Check if a triangle is incident to a super-triangle vertex.
    fn is_super_triangle(&self, tri: &Triangle) -> bool {
        tri.vertices.iter().any(|v| self.super_verts.contains(v))
    }

    /// Iterate over all interior (non-super) alive triangles.
    pub fn interior_triangles(&self) -> impl Iterator<Item = (TriangleId, &Triangle)> {
        self.triangles
            .iter()
            .enumerate()
            .filter(|(_, t)| t.alive && !self.is_super_triangle(t))
            .map(|(i, t)| (TriangleId::from_usize(i), t))
    }

    /// Iterate over all alive triangles (including super-triangle ones).
    pub fn all_alive_triangles(&self) -> impl Iterator<Item = (TriangleId, &Triangle)> {
        self.triangles
            .iter()
            .enumerate()
            .filter(|(_, t)| t.alive)
            .map(|(i, t)| (TriangleId::from_usize(i), t))
    }

    /// Access a vertex by ID.
    #[inline]
    #[must_use]
    pub fn vertex(&self, id: PslgVertexId) -> &PslgVertex {
        &self.vertices[id.idx()]
    }

    /// Access a triangle by ID.
    #[inline]
    #[must_use]
    pub fn triangle(&self, id: TriangleId) -> &Triangle {
        &self.triangles[id.idx()]
    }

    /// Access a mutable triangle by ID.
    #[inline]
    pub(crate) fn triangle_mut(&mut self, id: TriangleId) -> &mut Triangle {
        &mut self.triangles[id.idx()]
    }

    /// Access the full vertex slice.
    #[must_use]
    pub fn vertices(&self) -> &[PslgVertex] {
        &self.vertices
    }

    /// Access the full triangle slice.
    #[must_use]
    pub fn triangles_slice(&self) -> &[Triangle] {
        &self.triangles
    }

    /// Mutable access to the triangle slice.
    pub(crate) fn triangles_mut(&mut self) -> &mut Vec<Triangle> {
        &mut self.triangles
    }

    /// Insert a new vertex into the vertex pool and return its ID.
    ///
    /// Used by the refinement algorithm.
    pub(crate) fn add_vertex(&mut self, v: PslgVertex) -> PslgVertexId {
        let id = PslgVertexId::from_usize(self.vertices.len());
        self.vertices.push(v);
        self.vert_to_tri.push(GHOST_TRIANGLE);
        id
    }

    /// Insert a Steiner point into the triangulation.
    ///
    /// Returns the new vertex ID.
    pub(crate) fn insert_steiner(&mut self, x: Real, y: Real) -> PslgVertexId {
        let vid = self.add_vertex(PslgVertex::new(x, y));
        self.insert_vertex(vid);
        self.num_real_vertices += 1;
        vid
    }

    /// Read-only access to the vertex→triangle hint array.
    ///
    /// For each vertex `v`, `vert_to_tri_slice()[v.idx()]` holds one alive
    /// triangle incident to `v` (or [`GHOST_TRIANGLE`] if the vertex has not
    /// been inserted yet).
    ///
    /// # Invariant — Vertex-to-Triangle Validity
    ///
    /// **Statement**: After every insertion, for every inserted real vertex
    /// $v_i$ with $i < n$, $\text{vert\_to\_tri}[i]$ refers to an alive
    /// triangle whose vertex list contains $v_i$.
    ///
    /// **Proof**: Each of `insert_in_triangle`, `insert_on_edge`, and
    /// `insert_on_hull_edge` explicitly sets `vert_to_tri[v]` for every
    /// vertex of every newly created triangle.  The subsequent `flip_fix`
    /// calls update the hints for all four affected vertices after each
    /// flip.  Since flips only rearrange existing alive triangles (marking
    /// old ones dead and rewriting in-place), the hint is always updated
    /// to an alive triangle containing the vertex.  ∎
    #[must_use]
    pub fn vert_to_tri_slice(&self) -> &[TriangleId] {
        &self.vert_to_tri
    }

    /// Remove dead (tombstoned) triangles and remap all adjacency links.
    ///
    /// # Theorem — Compaction Invariant Preservation
    ///
    /// **Statement**: Let $T$ be a triangulation with $D$ dead and $A$ alive
    /// triangles.  `compact()` produces a triangulation $T'$ with exactly
    /// $A$ triangles such that:
    ///
    /// 1. Every alive triangle in $T$ maps bijectively to a triangle in $T'$.
    /// 2. Adjacency is preserved: $\text{adj}_{T'}(f(t), e) = f(\text{adj}_T(t, e))$
    ///    for all alive $t$ and edges $e$, where $f$ is the remapping.
    /// 3. Vertex→triangle hints remain valid.
    /// 4. Constrained-edge flags are preserved.
    ///
    /// **Proof sketch**: The remapping $f : [0 \dots A{+}D) \to [0 \dots A)$
    /// is a monotone injection on alive indices.  Because dead triangles are
    /// never referenced by alive triangles' adjacency lists (each insertion
    /// and flip only writes alive triangle IDs), every adjacency entry is
    /// either `GHOST_TRIANGLE` or an alive triangle ID.  Replacing each
    /// alive ID by its $f$-image preserves the bijection.  The `vert_to_tri`
    /// hints point to alive triangles by the vertex-to-triangle invariant
    /// above, so remapping them through $f$ keeps them valid.  ∎
    ///
    /// # Complexity
    ///
    /// $O(A + D)$ time and $O(A + D)$ auxiliary space for the remap table.
    pub fn compact(&mut self) {
        let n = self.triangles.len();
        let mut remap = vec![GHOST_TRIANGLE; n];
        let mut new_idx = 0u32;
        for i in 0..n {
            if self.triangles[i].alive {
                remap[i] = TriangleId::new(new_idx);
                new_idx += 1;
            }
        }

        // Compact the triangle array in-place.
        let mut write = 0;
        for read in 0..n {
            if self.triangles[read].alive {
                self.triangles[write] = self.triangles[read].clone();
                for a in &mut self.triangles[write].adj {
                    if *a != GHOST_TRIANGLE {
                        *a = remap[a.idx()];
                    }
                }
                write += 1;
            }
        }
        self.triangles.truncate(write);

        // Remap vertex→triangle hints.
        for hint in &mut self.vert_to_tri {
            if *hint != GHOST_TRIANGLE && hint.idx() < n {
                *hint = remap[hint.idx()];
            }
        }

        // Remap last_triangle.
        if self.last_triangle != GHOST_TRIANGLE && self.last_triangle.idx() < n {
            self.last_triangle = remap[self.last_triangle.idx()];
        } else {
            self.last_triangle = if self.triangles.is_empty() {
                GHOST_TRIANGLE
            } else {
                TriangleId::new(0)
            };
        }
    }
}
