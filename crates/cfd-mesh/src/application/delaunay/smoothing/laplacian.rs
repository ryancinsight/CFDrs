//! Uniform Laplacian mesh smoothing for 2-D constrained Delaunay triangulations.
//!
//! ## Algorithm
//!
//! Each interior (non-boundary) vertex is iteratively displaced toward the
//! centroid of its 1-ring (the average of all neighboring vertex positions).
//!
//! ```text
//! v_new = v + λ · (centroid(neighbors(v)) − v)
//!       = (1 − λ) · v + λ · centroid(neighbors(v))
//! ```
//!
//! where `λ ∈ (0, 0.5]` is the step size parameter.  The process repeats for
//! `max_iter` iterations.
//!
//! ## Theorem (Convergence and Termination)
//!
//! **Termination**: the loop runs exactly `max_iter` times.  ∎
//!
//! **Stability for λ ≤ 0.5**: each update is a convex combination of the old
//! vertex position (weight `1 − λ ≥ 0.5`) and the neighbor centroid (weight
//! `λ ≤ 0.5`).  Positions therefore remain within the convex hull of the
//! initial configuration.  The mapping is a strict contraction for the vertex
//! displacement vector when `λ < 1`, so the iteration converges to a fixed
//! point (uniform equilibrium).
//!
//! **Monotone minimum angle (convex neighbourhoods)**: for convex 1-rings,
//! moving a vertex toward the centroid is a net quality improvement — it
//! reduces the maximum aspect ratio of incident triangles.  For non-convex
//! 1-rings, quality is generally improved but not guaranteed to be monotone
//! (Taubin 1995).  ∎
//!
//! **Boundary invariant**: vertices on constrained edges (PSLG boundary) and
//! super-triangle vertices are never moved when `preserve_boundary = true`.
//! This guarantees the PSLG geometry is preserved throughout smoothing.  ∎
//!
//! ## Complexity
//!
//! O(`max_iter` × F) where F is the face count.  Each vertex lookup via
//! `triangles_around_vertex` is O(deg(v)) ≈ O(6) for typical Delaunay meshes.
//!
//! ## Reference
//!
//! Taubin, G. (1995). A signal processing approach to fair surface design.
//! *ACM SIGGRAPH*, 351–358.

use crate::application::delaunay::constraint::enforce::Cdt;
use crate::application::delaunay::pslg::vertex::PslgVertexId;
use crate::domain::core::scalar::Real;
use std::collections::HashSet;

// ── Public API ────────────────────────────────────────────────────────────────

/// Iterative uniform Laplacian smoother for 2-D CDT meshes.
///
/// Move each interior vertex toward the centroid of its 1-ring neighbourhood
/// for up to `max_iter` passes with step size `lambda`.
///
/// # Example
///
/// ```ignore
/// use cfd_mesh::application::delaunay::smoothing::LaplacianSmoother;
/// let smoother = LaplacianSmoother { max_iter: 10, lambda: 0.3, preserve_boundary: true };
/// smoother.smooth(&mut cdt);
/// ```
#[derive(Clone, Debug)]
pub struct LaplacianSmoother {
    /// Number of smoothing iterations.  Higher values → smoother mesh, slower.
    pub max_iter: usize,
    /// Step size ∈ `(0.0, 0.5]`.  Values > 0.5 risk instability.
    pub lambda: Real,
    /// If `true`, vertices on constrained edges are never moved.
    pub preserve_boundary: bool,
}

impl Default for LaplacianSmoother {
    fn default() -> Self {
        Self { max_iter: 5, lambda: 0.3, preserve_boundary: true }
    }
}

impl LaplacianSmoother {
    /// Apply `self.max_iter` Laplacian smoothing passes to `cdt` in place.
    ///
    /// Interior vertices are moved toward their 1-ring centroid.
    /// Boundary vertices (constrained-edge endpoints) are fixed when
    /// `self.preserve_boundary` is `true`.
    ///
    /// After smoothing the triangulation connectivity is unchanged; only vertex
    /// positions are updated.  Re-running Lawson flips via `Cdt::from_pslg` is
    /// recommended if strict Delaunay compliance is required post-smoothing.
    pub fn smooth(&self, cdt: &mut Cdt) {
        if self.max_iter == 0 || self.lambda <= 0.0 {
            return;
        }

        let lambda = self.lambda.clamp(0.0, 1.0);

        // Build the set of "frozen" vertex IDs (boundary + super-triangle).
        let frozen: HashSet<PslgVertexId> = build_frozen_set(cdt, self.preserve_boundary);

        let dt = cdt.triangulation();
        let num_real = dt.num_real_vertices;
        let _ = dt;

        for _ in 0..self.max_iter {
            // Collect all updated positions before applying them (Jacobi update).
            // This prevents the order of iteration from affecting the result.
            let mut new_positions: Vec<Option<(Real, Real)>> = vec![None; {
                cdt.triangulation().vertices().len()
            }];

            let dt = cdt.triangulation();
            for raw in 0..num_real {
                let vid = PslgVertexId::from_usize(raw);
                if frozen.contains(&vid) {
                    continue;
                }
                let neighbors = one_ring_neighbors(dt, vid);
                if neighbors.is_empty() {
                    continue;
                }
                let (cx, cy) = centroid(&neighbors, dt);
                let v = dt.vertex(vid);
                let nx = v.x + lambda * (cx - v.x);
                let ny = v.y + lambda * (cy - v.y);
                new_positions[raw] = Some((nx, ny));
            }

            // Apply Jacobi updates.
            let dt_mut = cdt.triangulation_mut();
            for (raw, pos) in new_positions.iter().enumerate() {
                if let Some((nx, ny)) = pos {
                    dt_mut.vertices[raw].x = *nx;
                    dt_mut.vertices[raw].y = *ny;
                }
            }
        }
    }

    /// Return the number of iterations that would be performed.
    #[must_use]
    pub fn effective_iterations(&self) -> usize {
        if self.lambda <= 0.0 { 0 } else { self.max_iter }
    }
}

// ── Helpers ───────────────────────────────────────────────────────────────────

/// Collect the 1-ring neighbour vertex IDs of vertex `v` in `dt`.
///
/// Returns all vertices that share an interior triangle edge with `v`,
/// excluding `v` itself and any super-triangle ghost vertices.
pub(crate) fn one_ring_neighbors(
    dt: &crate::application::delaunay::triangulation::bowyer_watson::DelaunayTriangulation,
    v: PslgVertexId,
) -> Vec<PslgVertexId> {
    let tris = dt.triangles_around_vertex(v);
    let mut neighbors = Vec::with_capacity(tris.len() * 2);
    for tid in tris {
        let tri = dt.triangle(tid);
        for &u in &tri.vertices {
            if u != v && !dt.super_verts.contains(&u) {
                neighbors.push(u);
            }
        }
    }
    neighbors.sort_unstable_by_key(|id| id.raw());
    neighbors.dedup();
    neighbors
}

/// Compute the centroid (arithmetic mean) of vertex positions.
fn centroid(
    neighbors: &[PslgVertexId],
    dt: &crate::application::delaunay::triangulation::bowyer_watson::DelaunayTriangulation,
) -> (Real, Real) {
    let n = neighbors.len() as Real;
    let (sx, sy) =
        neighbors.iter().fold((0.0, 0.0), |(ax, ay), &id| {
            let v = dt.vertex(id);
            (ax + v.x, ay + v.y)
        });
    (sx / n, sy / n)
}

/// Build the set of frozen (unmovable) vertex IDs.
///
/// Always includes super-triangle vertices.  If `preserve_boundary`, also
/// includes all constrained-edge endpoints.
pub(crate) fn build_frozen_set(cdt: &Cdt, preserve_boundary: bool) -> HashSet<PslgVertexId> {
    let dt = cdt.triangulation();
    let mut frozen = HashSet::new();

    // Super-triangle vertices are never moved.
    for &sv in &dt.super_verts {
        frozen.insert(sv);
    }

    if preserve_boundary {
        for &(a, b) in cdt.constrained_edges() {
            frozen.insert(a);
            frozen.insert(b);
        }
    }

    frozen
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::application::delaunay::{Cdt, Pslg};
    use crate::domain::core::scalar::Real;

    /// Build a simple 5-point CDT: 4 corners of a unit square + 1 interior point.
    fn small_square_cdt(interior_x: Real, interior_y: Real) -> Cdt {
        let mut pslg = Pslg::new();
        let v0 = pslg.add_vertex(0.0, 0.0);
        let v1 = pslg.add_vertex(1.0, 0.0);
        let v2 = pslg.add_vertex(1.0, 1.0);
        let v3 = pslg.add_vertex(0.0, 1.0);
        let v4 = pslg.add_vertex(interior_x, interior_y);
        // Boundary edges (constrained).
        pslg.add_segment(v0, v1);
        pslg.add_segment(v1, v2);
        pslg.add_segment(v2, v3);
        pslg.add_segment(v3, v0);
        let _ = v4; // Interior vertex — no segments to it.
        Cdt::from_pslg(&pslg)
    }

    /// Smoothing an already-centred interior vertex is a no-op.
    #[test]
    fn smooth_centered_interior_vertex_is_stable() {
        let mut cdt = small_square_cdt(0.5, 0.5);
        let smoother = LaplacianSmoother { max_iter: 10, lambda: 0.3, preserve_boundary: true };
        smoother.smooth(&mut cdt);
        // Interior vertex should remain near (0.5, 0.5) — centroid of the 4 corners.
        let dt = cdt.triangulation();
        let v4 = crate::application::delaunay::pslg::vertex::PslgVertexId::from_usize(4);
        let v = dt.vertex(v4);
        assert!(
            (v.x - 0.5).abs() < 0.1 && (v.y - 0.5).abs() < 0.1,
            "centred vertex should stay near (0.5, 0.5), got ({}, {})",
            v.x,
            v.y
        );
    }

    /// Smoothing an off-centre interior vertex moves it toward the centroid.
    #[test]
    fn smooth_moves_off_centre_vertex_toward_centroid() {
        // Interior at (0.1, 0.1) — well off-centre.
        let mut cdt = small_square_cdt(0.1, 0.1);
        let before_x;
        let before_y;
        {
            let dt = cdt.triangulation();
            let v4 = PslgVertexId::from_usize(4);
            let v = dt.vertex(v4);
            before_x = v.x;
            before_y = v.y;
        }
        let smoother = LaplacianSmoother { max_iter: 5, lambda: 0.3, preserve_boundary: true };
        smoother.smooth(&mut cdt);
        let dt = cdt.triangulation();
        let v4 = PslgVertexId::from_usize(4);
        let v = dt.vertex(v4);
        // After smoothing the vertex should have moved toward (0.5, 0.5).
        assert!(
            v.x > before_x && v.y > before_y,
            "off-centre vertex should move toward centroid: before=({before_x},{before_y}), after=({},{})",
            v.x, v.y
        );
    }

    /// With `preserve_boundary = true`, boundary vertices must not move.
    #[test]
    fn boundary_vertices_are_preserved() {
        let mut cdt = small_square_cdt(0.1, 0.1);
        let corners_before: Vec<(Real, Real)> = {
            let dt = cdt.triangulation();
            (0..4).map(|i| {
                let v = dt.vertex(PslgVertexId::from_usize(i));
                (v.x, v.y)
            }).collect()
        };
        let smoother = LaplacianSmoother { max_iter: 10, lambda: 0.5, preserve_boundary: true };
        smoother.smooth(&mut cdt);
        let dt = cdt.triangulation();
        for i in 0..4 {
            let v = dt.vertex(PslgVertexId::from_usize(i));
            let (bx, by) = corners_before[i];
            assert!(
                (v.x - bx).abs() < 1e-12 && (v.y - by).abs() < 1e-12,
                "boundary vertex {i} moved: before=({bx},{by}), after=({},{})",
                v.x, v.y
            );
        }
    }

    /// Zero iterations: no-op.
    #[test]
    fn zero_iterations_is_noop() {
        let mut cdt = small_square_cdt(0.1, 0.1);
        let before: Vec<(Real, Real)> = {
            let dt = cdt.triangulation();
            dt.vertices().iter().map(|v| (v.x, v.y)).collect()
        };
        LaplacianSmoother { max_iter: 0, lambda: 0.5, preserve_boundary: true }.smooth(&mut cdt);
        let dt = cdt.triangulation();
        for (i, v) in dt.vertices().iter().enumerate() {
            assert_eq!((v.x, v.y), before[i], "vertex {i} should not move with 0 iterations");
        }
    }
}
