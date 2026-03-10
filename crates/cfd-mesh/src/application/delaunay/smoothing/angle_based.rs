//! Angle-based Laplacian mesh smoothing for 2-D CDT meshes.
//!
//! ## Algorithm — Angle-Defect–Weighted Laplacian
//!
//! This smoother augments uniform Laplacian smoothing with an angle-quality
//! guard: after each vertex move, the proposed new position is only accepted
//! if it does **not** decrease the minimum interior angle of any incident
//! triangle below the current minimum.
//!
//! Formally, for vertex `v` with 1-ring neighbours `u_i`:
//!
//! ```text
//! v_candidate = v + λ · (centroid(u_i) − v)
//! v_new = v_candidate   if min_angle(ring(v_candidate)) ≥ min_angle(ring(v))
//!       = v             otherwise  (revert — preserve existing quality)
//! ```
//!
//! This conservative acceptance criterion ensures that the minimum angle is
//! monotone non-decreasing across iterations (unlike the bare Laplacian,
//! which can decrease angles in non-convex 1-rings).
//!
//! ## Theorem (Monotone Minimum Angle)
//!
//! Let `θ_k = min_angle(mesh)` after iteration `k`.  By the acceptance
//! criterion, `θ_{k+1} ≥ θ_k` for every vertex update.  Therefore the
//! global minimum angle is monotone non-decreasing across all iterations
//! of `AngleBasedSmoother::smooth`.  QED.
//!
//! ## Theorem (Termination)
//!
//! The loop performs exactly `max_iter` iterations regardless of acceptance;
//! each iteration does O(F) work.  Total complexity: O(`max_iter` × F).  ∎
//!
//! ## Reference
//!
//! Canann, S. A., Tristano, J. R., & Staten, M. L. (1998). An approach to
//! combined Laplacian and optimisation-based smoothing for triangular,
//! quadrilateral, and quad-dominant meshes. *IMECE/ASME*, 95–104.

use crate::application::delaunay::constraint::enforce::Cdt;
use crate::application::delaunay::pslg::vertex::PslgVertexId;
use crate::application::delaunay::smoothing::laplacian::{build_frozen_set, one_ring_neighbors};
use crate::domain::core::scalar::Real;

// ── Public API ────────────────────────────────────────────────────────────────

/// Conservative angle-quality–preserving Laplacian smoother for 2-D CDT meshes.
///
/// Applies the Laplacian step only when the resulting minimum incident angle is
/// at least as large as before the move.  This guarantees monotone minimum-angle
/// improvement across all iterations.
#[derive(Clone, Debug)]
pub struct AngleBasedSmoother {
    /// Number of smoothing passes.
    pub max_iter: usize,
    /// Laplacian step size ∈ `(0.0, 0.5]`.
    pub lambda: Real,
    /// If `true`, vertices on constrained edges are never moved.
    pub preserve_boundary: bool,
}

impl Default for AngleBasedSmoother {
    fn default() -> Self {
        Self {
            max_iter: 5,
            lambda: 0.3,
            preserve_boundary: true,
        }
    }
}

impl AngleBasedSmoother {
    /// Apply angle-quality–guarded Laplacian smoothing to `cdt`.
    ///
    /// Each vertex move is accepted only if it does not decrease the minimum
    /// angle of any incident triangle.  This is a conservative approximation
    /// of ODT (Optimal Delaunay Triangulation) smoothing.
    pub fn smooth(&self, cdt: &mut Cdt) {
        if self.max_iter == 0 || self.lambda <= 0.0 {
            return;
        }
        let lambda = self.lambda.clamp(0.0, 1.0);
        let frozen = build_frozen_set(cdt, self.preserve_boundary);
        let num_real = cdt.triangulation().num_real_vertices;

        for _ in 0..self.max_iter {
            // Serial update (Gauss-Seidel): apply each accepted move immediately.
            // This can improve convergence rate vs. Jacobi at the cost of
            // order-dependence — acceptable for quality improvement heuristics.
            for raw in 0..num_real {
                let vid = PslgVertexId::from_usize(raw);
                if frozen.contains(&vid) {
                    continue;
                }
                let dt = cdt.triangulation();
                let neighbors = one_ring_neighbors(dt, vid);
                if neighbors.is_empty() {
                    continue;
                }

                // Compute centroid of 1-ring neighbours.
                let n = neighbors.len() as Real;
                let (sx, sy) = neighbors.iter().fold((0.0, 0.0), |(ax, ay), &id| {
                    let u = dt.vertex(id);
                    (ax + u.x, ay + u.y)
                });
                let (cx, cy) = (sx / n, sy / n);

                let vx = dt.vertex(vid).x;
                let vy = dt.vertex(vid).y;
                let nx = vx + lambda * (cx - vx);
                let ny = vy + lambda * (cy - vy);

                // Compute current minimum angle for incident triangles.
                let tris = dt.triangles_around_vertex(vid);
                let current_min = min_angle_in_ring(&tris, dt, vid, vx, vy);
                let proposed_min = min_angle_in_ring(&tris, dt, vid, nx, ny);

                // Accept move only if minimum angle is non-decreasing.
                if proposed_min >= current_min {
                    let dt_mut = cdt.triangulation_mut();
                    dt_mut.vertices[raw].x = nx;
                    dt_mut.vertices[raw].y = ny;
                }
            }
        }
    }
}

// ── Helpers ───────────────────────────────────────────────────────────────────

/// Minimum interior angle (in radians) across all triangles incident on `vid`,
/// with `vid` placed at `(vx, vy)` (may differ from current stored position).
fn min_angle_in_ring(
    tris: &[crate::application::delaunay::triangulation::triangle::TriangleId],
    dt: &crate::application::delaunay::triangulation::bowyer_watson::DelaunayTriangulation,
    vid: PslgVertexId,
    vx: Real,
    vy: Real,
) -> Real {
    use crate::application::delaunay::triangulation::triangle::GHOST_TRIANGLE;

    let mut min_ang = Real::INFINITY;
    for &tid in tris {
        if tid == GHOST_TRIANGLE {
            continue;
        }
        let tri = dt.triangle(tid);
        // Get the three vertex positions, substituting (vx,vy) for vid.
        let pts: [(Real, Real); 3] = tri.vertices.map(|u| {
            if u == vid {
                (vx, vy)
            } else {
                let w = dt.vertex(u);
                (w.x, w.y)
            }
        });
        let ang = min_angle_triangle(pts[0], pts[1], pts[2]);
        if ang < min_ang {
            min_ang = ang;
        }
    }
    min_ang
}

/// Minimum interior angle (radians) of a triangle with vertices `a`, `b`, `c`.
///
/// Uses the law of cosines: `cos θ = (e1² + e2² − e3²) / (2 e1 e2)` for each
/// vertex angle, then returns the minimum.
fn min_angle_triangle(a: (Real, Real), b: (Real, Real), c: (Real, Real)) -> Real {
    let edge_sq = |p: (Real, Real), q: (Real, Real)| {
        let dx = p.0 - q.0;
        let dy = p.1 - q.1;
        dx * dx + dy * dy
    };
    let ab = edge_sq(a, b);
    let bc = edge_sq(b, c);
    let ca = edge_sq(c, a);

    if ab < 1e-30 || bc < 1e-30 || ca < 1e-30 {
        return 0.0; // Degenerate triangle — treat as zero minimum angle.
    }

    // Angle at a = acos((ab + ca − bc) / (2 * sqrt(ab) * sqrt(ca)))
    // Using the clamped cos to handle rounding near ±1.
    let angle_at = |adj1: Real, adj2: Real, opp: Real| -> Real {
        let cos = (adj1 + adj2 - opp) / (2.0 * adj1.sqrt() * adj2.sqrt());
        cos.clamp(-1.0, 1.0).acos()
    };

    let ang_a = angle_at(ab, ca, bc);
    let ang_b = angle_at(ab, bc, ca);
    let ang_c = angle_at(bc, ca, ab);

    ang_a.min(ang_b).min(ang_c)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::application::delaunay::{Cdt, Pslg};

    /// Build a simple unit-square CDT with one interior point.
    fn square_cdt(ix: f64, iy: f64) -> Cdt {
        let mut p = Pslg::new();
        let v0 = p.add_vertex(0.0, 0.0);
        let v1 = p.add_vertex(1.0, 0.0);
        let v2 = p.add_vertex(1.0, 1.0);
        let v3 = p.add_vertex(0.0, 1.0);
        p.add_segment(v0, v1);
        p.add_segment(v1, v2);
        p.add_segment(v2, v3);
        p.add_segment(v3, v0);
        let _vi = p.add_vertex(ix, iy);
        Cdt::from_pslg(&p)
    }

    /// min_angle_triangle: equilateral triangle has all angles = 60° = π/3.
    #[test]
    fn equilateral_triangle_min_angle_is_60_deg() {
        use std::f64::consts::FRAC_PI_3;
        let ang = min_angle_triangle((0.0, 0.0), (1.0, 0.0), (0.5, 3.0_f64.sqrt() / 2.0));
        assert!(
            (ang - FRAC_PI_3).abs() < 1e-6,
            "equilateral min angle should be π/3, got {ang}"
        );
    }

    /// Angle-based smoother terminates and does not increase vertex count.
    #[test]
    fn angle_based_smooth_terminates_and_preserves_vertex_count() {
        let mut cdt = square_cdt(0.2, 0.2);
        let before_count = cdt.triangulation().vertex_count();
        AngleBasedSmoother::default().smooth(&mut cdt);
        assert_eq!(
            cdt.triangulation().vertex_count(),
            before_count,
            "smoothing must not add or remove vertices"
        );
    }

    /// Boundary vertices are unchanged after angle-based smoothing.
    #[test]
    fn angle_based_preserves_boundary_vertices() {
        let mut cdt = square_cdt(0.2, 0.2);
        let corners_before: Vec<(f64, f64)> = {
            let dt = cdt.triangulation();
            (0..4)
                .map(|i| {
                    let v = dt.vertex(PslgVertexId::from_usize(i));
                    (v.x, v.y)
                })
                .collect()
        };
        AngleBasedSmoother {
            max_iter: 10,
            lambda: 0.5,
            preserve_boundary: true,
        }
        .smooth(&mut cdt);
        let dt = cdt.triangulation();
        for i in 0..4 {
            let v = dt.vertex(PslgVertexId::from_usize(i));
            let (bx, by) = corners_before[i];
            assert!(
                (v.x - bx).abs() < 1e-12 && (v.y - by).abs() < 1e-12,
                "boundary vertex {i} moved from ({bx},{by}) to ({},{})",
                v.x,
                v.y
            );
        }
    }
}
