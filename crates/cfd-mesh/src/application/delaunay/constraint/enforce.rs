//! Constraint edge enforcement — turning a Delaunay triangulation into a CDT.
//!
//! # Algorithm — Edge-Flip Constraint Recovery
//!
//! Given a Delaunay triangulation and a constraint segment `(a, b)`:
//!
//! 1. If the edge `(a, b)` already exists in the triangulation, mark it as
//!    constrained and return.
//! 2. Otherwise, collect all edges that cross the segment `(a, b)`.
//! 3. Iteratively flip crossing edges.  After each flip, the new diagonal
//!    either no longer crosses `(a, b)` (progress) or it crosses but with a
//!    shorter intersection (bounded by the geometry).
//! 4. When the flipping sequence terminates, the edge `(a, b)` exists in the
//!    triangulation.  Mark it as constrained.
//!
//! # Theorem — Constraint Recovery Termination
//!
//! **Statement**: For a set of non-crossing constraint segments, the iterative
//! flip-based constraint recovery terminates in $O(n^2)$ flips in the worst
//! case, and $O(n)$ flips on average for well-distributed inputs.
//!
//! **Proof sketch**: Each flip replaces a crossing edge with a non-crossing
//! one relative to the constraint segment.  Since the number of possible
//! triangulations is finite and we make monotone progress (the set of
//! crossing edges shrinks), the process terminates.
//!
//! # Hole Removal
//!
//! After all constraints are enforced, triangles inside hole regions (identified
//! by flood-fill from hole seed points) are removed.

use std::collections::{HashSet, VecDeque};

use crate::domain::core::scalar::Real;
use crate::domain::geometry::predicates::{orient_2d, Orientation};
use nalgebra::Point2;

use crate::application::delaunay::pslg::graph::Pslg;
use crate::application::delaunay::pslg::graph::PslgValidationError;
use crate::application::delaunay::pslg::vertex::PslgVertexId;
use crate::application::delaunay::triangulation::bowyer_watson::DelaunayTriangulation;
use crate::application::delaunay::triangulation::triangle::{TriangleId, GHOST_TRIANGLE};

/// Constrained Delaunay Triangulation.
///
/// Wraps a [`DelaunayTriangulation`] with enforced constraint edges from a PSLG.
///
/// # Theorem — CDT Uniqueness
///
/// **Statement**: Given a valid PSLG (no crossing segments), the CDT is unique
/// up to co-circular degeneracies.  That is, if no four input points are
/// co-circular, the CDT is the unique triangulation that:
/// 1. Contains all constraint segments as edges.
/// 2. Maximises the minimum angle among all such triangulations.
///
/// **Proof sketch**: Among all triangulations containing the constraints, the
/// CDT is obtained by making every non-constrained edge locally Delaunay.
/// The Delaunay criterion (maximize minimum angle) has a unique optimum when
/// no four points are co-circular.
///
/// # Example
///
/// ```rust,ignore
/// use cfd_mesh::application::delaunay::{Pslg, Cdt};
///
/// let mut pslg = Pslg::new();
/// let a = pslg.add_vertex(0.0, 0.0);
/// let b = pslg.add_vertex(1.0, 0.0);
/// let c = pslg.add_vertex(0.5, 1.0);
/// pslg.add_segment(a, b);
/// pslg.add_segment(b, c);
/// pslg.add_segment(c, a);
///
/// let cdt = Cdt::from_pslg(&pslg);
/// assert!(cdt.triangulation().is_delaunay());
/// ```
pub struct Cdt {
    /// The underlying triangulation.
    dt: DelaunayTriangulation,
    /// The set of constrained edges (canonical vertex pairs).
    constrained_edges: HashSet<(PslgVertexId, PslgVertexId)>,
}

impl Cdt {
    /// Build a CDT from a PSLG.
    ///
    /// 1. Inserts all PSLG vertices into a Delaunay triangulation.
    /// 2. Enforces all constraint segments.
    /// 3. Removes triangles in hole regions.
    pub fn from_pslg(pslg: &Pslg) -> Self {
        Self::try_from_pslg(pslg)
            .unwrap_or_else(|err| panic!("Invalid PSLG for CDT construction: {err}"))
    }

    /// Build a CDT from a PSLG with validation.
    ///
    /// Returns an error if the PSLG contains degenerate/duplicate segments,
    /// out-of-range segment endpoints, or segment intersections.
    pub fn try_from_pslg(pslg: &Pslg) -> Result<Self, PslgValidationError> {
        pslg.validate()?;

        // Build the Delaunay triangulation from the PSLG vertices.
        let points: Vec<(Real, Real)> = pslg.vertices().iter().map(|v| (v.x, v.y)).collect();

        let dt = DelaunayTriangulation::from_points(&points);

        let mut cdt = Self {
            dt,
            constrained_edges: HashSet::new(),
        };

        // Enforce each constraint segment.
        for seg in pslg.segments() {
            cdt.enforce_constraint(seg.start, seg.end);
        }

        // Remove hole triangles.
        for hole in pslg.holes() {
            cdt.remove_hole_triangles(hole.x, hole.y);
        }

        Ok(cdt)
    }

    /// Enforce a constraint edge between vertices `a` and `b`.
    fn enforce_constraint(&mut self, a: PslgVertexId, b: PslgVertexId) {
        // Record the constraint.
        let canonical = if a <= b { (a, b) } else { (b, a) };
        self.constrained_edges.insert(canonical);

        // Check if the edge already exists.
        if self.edge_exists(a, b) {
            self.mark_edge_constrained(a, b);
            return;
        }

        // Flip edges crossing the segment (a, b).
        self.flip_to_enforce(a, b);

        // Mark the edge as constrained.
        self.mark_edge_constrained(a, b);
    }

    /// Check if an edge between vertices `a` and `b` exists in the triangulation.
    ///
    /// Uses vertex-star walk: O(deg(a)) ≈ O(6) instead of O(T).
    fn edge_exists(&self, a: PslgVertexId, b: PslgVertexId) -> bool {
        self.dt.edge_exists_fast(a, b)
    }

    /// Flip crossing edges until `(a, b)` exists in the triangulation.
    ///
    /// Uses the iterative approach: collect edges that intersect the segment
    /// `(a, b)`, then flip them one by one.
    /// Flip edges to recover the constraint `(a, b)`.
    ///
    /// Uses the Sloan (1993) queue-based algorithm:
    /// 1. Collect **all** edges crossing the segment `(a, b)`.
    /// 2. Process the queue round-robin: for each crossing edge, attempt a flip.
    ///    - If the quad is convex and the flip succeeds, check whether the new
    ///      diagonal still crosses `(a, b)`.  If so, re-enqueue it.
    ///    - If the quad is NOT convex, defer the edge to the back of the queue.
    /// 3. When all crossings are resolved the edge `(a, b)` exists.
    ///
    /// Steiner point insertion is only used as a last-resort fallback when a
    /// full pass through the queue makes zero progress.
    fn flip_to_enforce(&mut self, a: PslgVertexId, b: PslgVertexId) {
        // Collect all crossing edges into a work queue.
        let mut queue: VecDeque<(TriangleId, usize)> = VecDeque::new();
        self.collect_crossing_edges(a, b, &mut queue);
        if queue.is_empty() {
            return;
        }

        // Safety bound: O(n²) iterations for n crossing edges.
        let max_rounds = (queue.len() + 1) * (queue.len() + 1) * 2;
        let mut iter = 0usize;

        while let Some((tid, edge)) = queue.pop_front() {
            iter += 1;
            if iter > max_rounds {
                // Exceeded iteration budget — fallback: insert Steiner point
                // for first remaining crossing and restart.
                self.steiner_fallback(a, b, tid, edge);
                // Re-collect crossings and try again.
                queue.clear();
                self.collect_crossing_edges(a, b, &mut queue);
                if queue.is_empty() {
                    return;
                }
                iter = 0;
                continue;
            }

            // Validate: the triangle might have been killed by a prior flip.
            if !self.dt.triangle(tid).alive {
                continue;
            }

            // Verify the edge still crosses (a, b).
            let (va, vb) = self.dt.triangle(tid).edge_vertices(edge);
            {
                let ea = self.dt.vertex(va).to_point2();
                let eb = self.dt.vertex(vb).to_point2();
                let sa = self.dt.vertex(a).to_point2();
                let sb = self.dt.vertex(b).to_point2();
                if !segments_cross(&sa, &sb, &ea, &eb) {
                    continue; // No longer crossing — already resolved.
                }
            }

            if self.can_flip(tid, edge) {
                // Record the new diagonal endpoints BEFORE the flip mutates state.
                let v_opp_t = self.dt.triangle(tid).vertices[edge];
                let nbr_tid = self.dt.triangle(tid).adj[edge];
                let nbr_edge = self
                    .dt
                    .triangle(nbr_tid)
                    .shared_edge(tid)
                    .expect("adjacency broken");
                let v_opp_n = self.dt.triangle(nbr_tid).vertices[nbr_edge];

                self.perform_flip(tid, edge);

                // Check if the new diagonal (v_opp_t — v_opp_n) crosses (a, b).
                let new_a = self.dt.vertex(v_opp_t).to_point2();
                let new_b = self.dt.vertex(v_opp_n).to_point2();
                let sa = self.dt.vertex(a).to_point2();
                let sb = self.dt.vertex(b).to_point2();
                if segments_cross(&sa, &sb, &new_a, &new_b) {
                    // New diagonal still crosses — find which triangle now
                    // contains that edge and re-enqueue.
                    if let Some(new_loc) = self.find_edge_in_triangles(v_opp_t, v_opp_n) {
                        queue.push_back(new_loc);
                    }
                }
            } else {
                // Not convex — defer to the back of the queue.
                queue.push_back((tid, edge));
            }
        }
    }

    /// Collect all edges currently crossing segment `(a, b)`.
    ///
    /// Walks from triangle at `a` towards `b`, collecting crossing edges
    /// along the way.  Falls back to a full scan only when the directed
    /// walk fails (e.g. tombstoned triangles).
    fn collect_crossing_edges(
        &self,
        a: PslgVertexId,
        b: PslgVertexId,
        queue: &mut VecDeque<(TriangleId, usize)>,
    ) {
        let sa = self.dt.vertex(a).to_point2();
        let sb = self.dt.vertex(b).to_point2();

        // Directed walk: start from triangles around `a`, walk towards `b`.
        let star_a = self.dt.triangles_around_vertex(a);
        let mut visited = HashSet::new();

        // Find the starting triangle whose edge opposite `a` crosses (a,b).
        let mut walk_queue: VecDeque<TriangleId> = VecDeque::new();
        for &tid in &star_a {
            let tri = self.dt.triangle(tid);
            for edge in 0..3 {
                let (va, vb) = tri.edge_vertices(edge);
                if va == a || va == b || vb == a || vb == b {
                    continue;
                }
                if tri.constrained[edge] {
                    continue;
                }
                let ea = self.dt.vertex(va).to_point2();
                let eb = self.dt.vertex(vb).to_point2();
                if segments_cross(&sa, &sb, &ea, &eb) {
                    queue.push_back((tid, edge));
                    visited.insert(tid);
                    let nbr = tri.adj[edge];
                    if nbr != GHOST_TRIANGLE && !visited.contains(&nbr) {
                        walk_queue.push_back(nbr);
                    }
                }
            }
        }

        // Continue walking through neighbors of found crossings.
        while let Some(tid) = walk_queue.pop_front() {
            if visited.contains(&tid) {
                continue;
            }
            visited.insert(tid);
            let tri = self.dt.triangle(tid);
            if !tri.alive {
                continue;
            }
            for edge in 0..3 {
                let (va, vb) = tri.edge_vertices(edge);
                if va == a || va == b || vb == a || vb == b {
                    // Check if we reached b — stop walking this direction.
                    continue;
                }
                if tri.constrained[edge] {
                    continue;
                }
                let ea = self.dt.vertex(va).to_point2();
                let eb = self.dt.vertex(vb).to_point2();
                if segments_cross(&sa, &sb, &ea, &eb) {
                    queue.push_back((tid, edge));
                    let nbr = tri.adj[edge];
                    if nbr != GHOST_TRIANGLE && !visited.contains(&nbr) {
                        walk_queue.push_back(nbr);
                    }
                }
            }
        }
    }

    /// Find the triangle and local edge index containing edge `(u, v)`.
    ///
    /// Uses vertex-star walk: O(deg(u)) ≈ O(6) instead of O(T).
    fn find_edge_in_triangles(
        &self,
        u: PslgVertexId,
        v: PslgVertexId,
    ) -> Option<(TriangleId, usize)> {
        self.dt.find_edge_fast(u, v)
    }

    /// Last-resort: insert a Steiner point at the intersection of the
    /// constraint segment and a crossing edge.
    fn steiner_fallback(&mut self, a: PslgVertexId, b: PslgVertexId, tid: TriangleId, edge: usize) {
        let (va, vb) = self.dt.triangle(tid).edge_vertices(edge);
        let pa = self.dt.vertex(va);
        let pb = self.dt.vertex(vb);
        let pc = self.dt.vertex(a);
        let pd = self.dt.vertex(b);
        if let Some((ix, iy)) = segment_intersection(pa.x, pa.y, pb.x, pb.y, pc.x, pc.y, pd.x, pd.y)
        {
            self.dt.insert_steiner(ix, iy);
        }
    }

    /// Check if the quadrilateral formed by flipping edge `edge` of triangle
    /// `tid` is convex (required for a valid flip).
    fn can_flip(&self, tid: TriangleId, edge: usize) -> bool {
        let tri = self.dt.triangle(tid);
        let nbr_tid = tri.adj[edge];
        if nbr_tid == GHOST_TRIANGLE {
            return false;
        }

        let nbr = self.dt.triangle(nbr_tid);
        let nbr_edge = match nbr.shared_edge(tid) {
            Some(e) => e,
            None => return false,
        };

        let v_opp_t = tri.vertices[edge];
        let v_opp_n = nbr.vertices[nbr_edge];
        let (va, vb) = tri.edge_vertices(edge);

        let pa = self.dt.vertex(va).to_point2();
        let pb = self.dt.vertex(vb).to_point2();
        let pc = self.dt.vertex(v_opp_t).to_point2();
        let pd = self.dt.vertex(v_opp_n).to_point2();

        // The quad (v_opp_t, va, v_opp_n, vb) must be convex.
        // Check that the diagonals cross.
        let o1 = orient_2d(&pc, &pd, &pa);
        let o2 = orient_2d(&pc, &pd, &pb);
        let o3 = orient_2d(&pa, &pb, &pc);
        let o4 = orient_2d(&pa, &pb, &pd);

        (o1 != o2 && o1 != Orientation::Degenerate && o2 != Orientation::Degenerate)
            && (o3 != o4 && o3 != Orientation::Degenerate && o4 != Orientation::Degenerate)
    }

    /// Perform an edge flip on edge `edge` of triangle `tid`.
    fn perform_flip(&mut self, tid: TriangleId, edge: usize) {
        let tri = self.dt.triangle(tid).clone();
        let nbr_tid = tri.adj[edge];
        let nbr = self.dt.triangle(nbr_tid).clone();

        let nbr_edge = nbr.shared_edge(tid).expect("adjacency broken");

        let v_opp_t = tri.vertices[edge];
        let v_opp_n = nbr.vertices[nbr_edge];
        let (va, vb) = tri.edge_vertices(edge);

        let adj_tid_1 = tri.adj[(edge + 1) % 3];
        let adj_tid_2 = tri.adj[(edge + 2) % 3];
        // In the neighbor (CCW), shared edge traverses vb→va:
        //   adj[(nbr_edge+1)%3] = opp vb → edge (v_opp_n, va)
        //   adj[(nbr_edge+2)%3] = opp va → edge (v_opp_n, vb)
        let adj_nbr_opp_vb = nbr.adj[(nbr_edge + 1) % 3]; // across (v_opp_n, va)
        let adj_nbr_opp_va = nbr.adj[(nbr_edge + 2) % 3]; // across (v_opp_n, vb)

        let tri_cons = tri.constrained;
        let nbr_cons = nbr.constrained;

        let ts = self.dt.triangles_mut();

        // Rewrite tid → (v_opp_t, v_opp_n, vb)
        //   edge 0 opp v_opp_t = (v_opp_n, vb) → adj_nbr_opp_va
        //   edge 1 opp v_opp_n = (vb, v_opp_t) → adj_tid_1
        //   edge 2 opp vb      = shared → nbr_tid
        ts[tid.idx()].vertices = [v_opp_t, v_opp_n, vb];
        ts[tid.idx()].adj = [adj_nbr_opp_va, adj_tid_1, nbr_tid];
        ts[tid.idx()].constrained = [
            nbr_cons[(nbr_edge + 2) % 3],
            tri_cons[(edge + 1) % 3],
            false, // New diagonal is unconstrained.
        ];

        // Rewrite nbr → (v_opp_n, v_opp_t, va)
        //   edge 0 opp v_opp_n = (v_opp_t, va) → adj_tid_2
        //   edge 1 opp v_opp_t = (va, v_opp_n) → adj_nbr_opp_vb
        //   edge 2 opp va      = shared → tid
        ts[nbr_tid.idx()].vertices = [v_opp_n, v_opp_t, va];
        ts[nbr_tid.idx()].adj = [adj_tid_2, adj_nbr_opp_vb, tid];
        ts[nbr_tid.idx()].constrained = [
            tri_cons[(edge + 2) % 3],
            nbr_cons[(nbr_edge + 1) % 3],
            false, // New diagonal is unconstrained.
        ];

        // Fix external adjacency: edges that moved between triangle IDs.
        if adj_nbr_opp_va != GHOST_TRIANGLE {
            for a in &mut ts[adj_nbr_opp_va.idx()].adj {
                if *a == nbr_tid {
                    *a = tid;
                    break;
                }
            }
        }
        if adj_tid_2 != GHOST_TRIANGLE {
            for a in &mut ts[adj_tid_2.idx()].adj {
                if *a == tid {
                    *a = nbr_tid;
                    break;
                }
            }
        }

        // Update vert_to_tri for affected vertices.
        self.dt.vert_to_tri[v_opp_t.idx()] = tid;
        self.dt.vert_to_tri[v_opp_n.idx()] = nbr_tid;
        self.dt.vert_to_tri[va.idx()] = nbr_tid;
        self.dt.vert_to_tri[vb.idx()] = tid;
    }

    /// Mark an edge between `a` and `b` as constrained in both triangles
    /// that share it.
    ///
    /// Uses vertex-star walk: O(deg(a)) ≈ O(6) instead of O(T).
    fn mark_edge_constrained(&mut self, a: PslgVertexId, b: PslgVertexId) {
        let tris = self.dt.triangles_around_vertex(a);
        for tid in tris {
            let tri = self.dt.triangle_mut(tid);
            for edge in 0..3 {
                let (va, vb) = tri.edge_vertices(edge);
                if (va == a && vb == b) || (va == b && vb == a) {
                    tri.constrained[edge] = true;
                }
            }
        }
    }

    /// Remove triangles that are inside a hole region via flood-fill.
    fn remove_hole_triangles(&mut self, hx: Real, hy: Real) {
        // Find the triangle containing the hole seed.
        let start = {
            let first_alive = self
                .dt
                .triangles_slice()
                .iter()
                .position(|t| t.alive)
                .map(TriangleId::from_usize);

            match first_alive {
                Some(s) => s,
                None => return,
            }
        };

        use crate::application::delaunay::triangulation::locate::{locate, Location};

        let loc = locate(self.dt.vertices(), self.dt.triangles_slice(), start, hx, hy);
        let seed_tid = match loc {
            Some(Location::Inside(tid))
            | Some(Location::OnEdge(tid, _))
            | Some(Location::OnVertex(tid, _)) => tid,
            None => return,
        };

        // Flood-fill from the seed, stopping at constrained edges.
        let mut queue = VecDeque::new();
        let mut visited = HashSet::new();
        queue.push_back(seed_tid);
        visited.insert(seed_tid);

        while let Some(tid) = queue.pop_front() {
            let tri = self.dt.triangle(tid).clone();
            self.dt.triangle_mut(tid).alive = false;

            for edge in 0..3 {
                if tri.constrained[edge] {
                    continue; // Don't cross constraint edges.
                }
                let nbr = tri.adj[edge];
                if nbr != GHOST_TRIANGLE && !visited.contains(&nbr) {
                    if self.dt.triangle(nbr).alive {
                        visited.insert(nbr);
                        queue.push_back(nbr);
                    }
                }
            }
        }
    }

    // ── Public API ────────────────────────────────────────────────────────

    /// Access the underlying triangulation.
    pub fn triangulation(&self) -> &DelaunayTriangulation {
        &self.dt
    }

    /// Consume the CDT, returning the underlying triangulation.
    pub fn into_triangulation(self) -> DelaunayTriangulation {
        self.dt
    }

    /// Mutable access to the underlying triangulation.
    pub(crate) fn triangulation_mut(&mut self) -> &mut DelaunayTriangulation {
        &mut self.dt
    }

    /// Check if an edge is constrained.
    pub fn is_constrained(&self, a: PslgVertexId, b: PslgVertexId) -> bool {
        let canonical = if a <= b { (a, b) } else { (b, a) };
        self.constrained_edges.contains(&canonical)
    }

    /// The set of all constrained edges.
    pub fn constrained_edges(&self) -> &HashSet<(PslgVertexId, PslgVertexId)> {
        &self.constrained_edges
    }

    /// Insert a new constraint edge (used by refinement after segment splits).
    pub(crate) fn add_constraint(&mut self, a: PslgVertexId, b: PslgVertexId) {
        self.enforce_constraint(a, b);
    }

    /// Remove a constraint edge (used when splitting a segment).
    ///
    /// Uses vertex-star walk: O(deg(a)) instead of O(T).
    pub(crate) fn remove_constraint(&mut self, a: PslgVertexId, b: PslgVertexId) {
        let canonical = if a <= b { (a, b) } else { (b, a) };
        self.constrained_edges.remove(&canonical);
        // Unmark the edge in the triangulation (it may no longer exist).
        let tris = self.dt.triangles_around_vertex(a);
        for tid in tris {
            let tri = self.dt.triangle_mut(tid);
            for edge in 0..3 {
                let (va, vb) = tri.edge_vertices(edge);
                if (va == a && vb == b) || (va == b && vb == a) {
                    tri.constrained[edge] = false;
                }
            }
        }
    }
}

// ── Geometric helpers ─────────────────────────────────────────────────────────

/// Test if two line segments properly cross (share an interior point).
fn segments_cross(
    a1: &Point2<Real>,
    a2: &Point2<Real>,
    b1: &Point2<Real>,
    b2: &Point2<Real>,
) -> bool {
    let o1 = orient_2d(a1, a2, b1);
    let o2 = orient_2d(a1, a2, b2);
    let o3 = orient_2d(b1, b2, a1);
    let o4 = orient_2d(b1, b2, a2);

    // Proper crossing: endpoints of each segment lie on opposite sides of the other.
    o1 != o2
        && o3 != o4
        && o1 != Orientation::Degenerate
        && o2 != Orientation::Degenerate
        && o3 != Orientation::Degenerate
        && o4 != Orientation::Degenerate
}

/// Compute the intersection point of two line segments.
///
/// Returns `None` if they are parallel or do not intersect.
fn segment_intersection(
    ax1: Real,
    ay1: Real,
    ax2: Real,
    ay2: Real,
    bx1: Real,
    by1: Real,
    bx2: Real,
    by2: Real,
) -> Option<(Real, Real)> {
    let dx_a = ax2 - ax1;
    let dy_a = ay2 - ay1;
    let dx_b = bx2 - bx1;
    let dy_b = by2 - by1;

    let denom = dx_a * dy_b - dy_a * dx_b;
    if denom.abs() < 1e-30 {
        return None; // Parallel.
    }

    let t = ((bx1 - ax1) * dy_b - (by1 - ay1) * dx_b) / denom;

    Some((ax1 + t * dx_a, ay1 + t * dy_a))
}
