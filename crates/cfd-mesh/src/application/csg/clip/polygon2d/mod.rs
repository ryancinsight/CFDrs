//! 2-D polygon Boolean operations (CDT-backed production path + test oracles).
//!
//! Canonical SPI for coplanar polygon Boolean algebra.  All callers outside
//! this sub-module should import from [`crate::application::csg::clip`] rather
//! than addressing `polygon2d` directly.
//!
//! ## Algorithms
//!
//! | Algorithm | Strengths | Weaknesses |
//! |-----------|-----------|------------|
//! | Sutherland-Hodgman | O(n) per half-plane; exact for convex clips | Convex clip region only |
//! | CDT-based | Handles all degeneracies via triangulation | O(n log n) overhead |
//!
//! ## Public SPI
//!
//! - [`boolean_clip`]: canonical polygon Boolean entrypoint (CDT-backed)
//! - [`clip_polygon_to_triangle`]: specialized convex-triangle clip helper used
//!   by the coplanar CSG pipeline
//! - [`split_polygon_outside_triangle`]: specialized complement decomposition
//!   helper used by the coplanar CSG pipeline
//!
//! Classical algorithms (e.g., Sutherland-Hodgman)
//! are retained as internal utilities (half-plane clipping).
//!
//! ## References
//!
//! - Sutherland & Hodgman (1974), "Reentrant polygon clipping"

pub(crate) mod cdt;
pub(crate) mod geometry;
pub(crate) mod sutherland_hodgman;
pub use geometry::polygon_area;

use crate::domain::core::scalar::Real;
use crate::domain::geometry::predicates::{orient_2d_arr, Orientation};
use cdt::cdt_clip;
use sutherland_hodgman::{sh_clip_convex, sh_clip_halfplane};

/// 2-D Boolean operation type for polygon clipping.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ClipOp {
    /// Subject ∩ Clip
    Intersection,
    /// Subject ∪ Clip
    Union,
    /// Subject \ Clip
    Difference,
}

/// Canonical polygon Boolean clipping entrypoint.
///
/// This is the only production SPI for general polygon-vs-polygon Boolean
/// clipping. Internally it uses the CDT-based constrained arrangement backend,
/// which is robust under shared vertices/edges and other degeneracies.
pub fn boolean_clip(subject: &[[Real; 2]], clip: &[[Real; 2]], op: ClipOp) -> Vec<Vec<[Real; 2]>> {
    cdt_clip(subject, clip, op)
}

// ── Convenience wrappers ──────────────────────────────────────────────────────

/// Clip a 2-D polygon to the inside of a CCW triangle.
///
/// Used by the coplanar CSG pipeline.
pub fn clip_polygon_to_triangle(
    poly: &[[Real; 2]],
    dx: Real,
    dy: Real,
    ex: Real,
    ey: Real,
    fx: Real,
    fy: Real,
) -> Vec<[Real; 2]> {
    if poly.len() < 3 {
        return Vec::new();
    }

    // Ensure CCW triangle winding.
    let d = [dx, dy];
    let e = [ex, ey];
    let f = [fx, fy];
    let ori = orient_2d_arr(d, e, f);

    let tri = match ori {
        Orientation::Negative => vec![d, f, e],
        Orientation::Degenerate => return Vec::new(),
        Orientation::Positive => vec![d, e, f],
    };

    sh_clip_convex(poly, &tri)
}

/// Decompose `poly \ triangle(d,e,f)` into disjoint pieces.
///
/// Used by the coplanar CSG pipeline.
pub fn split_polygon_outside_triangle(
    poly: &[[Real; 2]],
    dx: Real,
    dy: Real,
    ex: Real,
    ey: Real,
    fx: Real,
    fy: Real,
) -> Vec<Vec<[Real; 2]>> {
    if poly.len() < 3 {
        return Vec::new();
    }

    // Ensure CCW triangle winding.
    let d = [dx, dy];
    let e = [ex, ey];
    let f = [fx, fy];
    let ori = orient_2d_arr(d, e, f);

    let (dx, dy, ex, ey, fx, fy) = match ori {
        Orientation::Negative => (dx, dy, fx, fy, ex, ey),
        Orientation::Degenerate => return vec![poly.to_vec()],
        Orientation::Positive => (dx, dy, ex, ey, fx, fy),
    };

    // 3-piece complement decomposition via S-H half-plane clips.
    //
    // For a CCW triangle (d,e,f), the complement of its interior is:
    //   R² \ tri = right(d→e) ∪ (left(d→e) ∩ right(e→f)) ∪ (left(d→e) ∩ left(e→f) ∩ right(f→d))
    //
    // These three regions are disjoint and partition the complement.
    let mut out = Vec::with_capacity(3);

    // piece_0: right of d→e (= outside edge d→e for CCW triangle)
    let p0 = sh_clip_halfplane(poly, ex, ey, dx, dy); // reversed = right side
    if p0.len() >= 3 {
        out.push(p0);
    }

    // piece_1 & 2: first restrict to left of d→e (inside that edge)
    let in_de = sh_clip_halfplane(poly, dx, dy, ex, ey);
    if in_de.len() >= 3 {
        // piece_1: left(d→e) ∩ right(e→f)
        let p1 = sh_clip_halfplane(&in_de, fx, fy, ex, ey);
        if p1.len() >= 3 {
            out.push(p1);
        }

        // piece_2: left(d→e) ∩ left(e→f) ∩ right(f→d)
        let in_ef = sh_clip_halfplane(&in_de, ex, ey, fx, fy);
        if in_ef.len() >= 3 {
            let p2 = sh_clip_halfplane(&in_ef, dx, dy, fx, fy);
            if p2.len() >= 3 {
                out.push(p2);
            }
        }
    }
    out
}

// ══════════════════════════════════════════════════════════════════════════════
// Tests
// ══════════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(a: Real, b: Real, tol: Real) -> bool {
        (a - b).abs() < tol
    }

    fn total_area(polys: &[Vec<[Real; 2]>]) -> Real {
        polys.iter().map(|p| polygon_area(p)).sum()
    }

    // ── clip_polygon_to_triangle ─────────────────────────────────────────────

    #[test]
    fn clip_to_triangle_full_inside() {
        let small = vec![[0.2, 0.2], [0.4, 0.2], [0.3, 0.4]];
        let result = clip_polygon_to_triangle(&small, 0.0, 0.0, 1.0, 0.0, 0.5, 1.0);
        let area = polygon_area(&result);
        let expected = polygon_area(&small);
        assert!(
            approx_eq(area, expected, 1e-10),
            "small triangle inside big triangle"
        );
    }

    #[test]
    fn clip_to_triangle_partial_overlap() {
        let sq = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let result = clip_polygon_to_triangle(&sq, -0.5, -0.5, 1.5, -0.5, 0.5, 0.5);
        assert!(result.len() >= 3, "should produce a valid polygon");
        let area = polygon_area(&result);
        assert!(area > 0.0 && area < 1.0, "partial overlap area={area}");
    }

    // ── split_polygon_outside_triangle ───────────────────────────────────────

    #[test]
    fn split_outside_area_conservation() {
        let poly = vec![[0.0, 0.0], [2.0, 0.0], [1.0, 2.0]];
        let (dx, dy, ex, ey, fx, fy) = (0.5, 0.0, 1.5, 0.0, 1.0, 1.0);
        let inside = clip_polygon_to_triangle(&poly, dx, dy, ex, ey, fx, fy);
        let outside = split_polygon_outside_triangle(&poly, dx, dy, ex, ey, fx, fy);

        let area_poly = polygon_area(&poly);
        let area_in = polygon_area(&inside);
        let area_out: Real = outside.iter().map(|p| polygon_area(p)).sum();
        let err = ((area_in + area_out) - area_poly).abs();
        assert!(err < 0.01, "area conservation: err={err:.2e}");
    }

    // ── Inclusion-exclusion identity ─────────────────────────────────────────

    #[test]
    fn canonical_boolean_clip_inclusion_exclusion_identity() {
        let a = vec![[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
        let b = vec![[1.0, 1.0], [3.0, 1.0], [3.0, 3.0], [1.0, 3.0]];

        let area_a = polygon_area(&a);
        let area_b = polygon_area(&b);

        let union = boolean_clip(&a, &b, ClipOp::Union);
        let inter = boolean_clip(&a, &b, ClipOp::Intersection);

        let area_union: Real = union.iter().map(|p| polygon_area(p)).sum();
        let area_inter: Real = inter.iter().map(|p| polygon_area(p)).sum();

        let lhs = area_a + area_b;
        let rhs = area_union + area_inter;
        let err = (lhs - rhs).abs() / lhs;
        assert!(
            err < 0.05,
            "inclusion-exclusion: |A|+|B|={lhs:.4} |A∪B|+|A∩B|={rhs:.4} err={:.1}%",
            err * 100.0
        );
    }


    // ── Oracle-comparison regressions for known degeneracies ────────────────

    #[test]
    fn canonical_shared_edge_contact_regression() {
        let a = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let b = vec![[1.0, 0.0], [2.0, 0.0], [2.0, 1.0], [1.0, 1.0]];

        for (op, expected_area) in [
            (ClipOp::Union, 2.0),
            (ClipOp::Intersection, 0.0),
            (ClipOp::Difference, 1.0),
        ] {
            let canonical = boolean_clip(&a, &b, op);
            let a_can = total_area(&canonical);
            assert!(
                approx_eq(a_can, expected_area, 0.05),
                "shared-edge {op:?}: canonical area={a_can}, expected={expected_area}",
            );
        }
    }

    #[test]
    fn canonical_touching_vertex_contact_regression() {
        let a = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let b = vec![[1.0, 1.0], [2.0, 1.0], [2.0, 2.0], [1.0, 2.0]];

        for (op, expected_area) in [
            (ClipOp::Union, 2.0),
            (ClipOp::Intersection, 0.0),
            (ClipOp::Difference, 1.0),
        ] {
            let canonical = boolean_clip(&a, &b, op);
            let a_can = total_area(&canonical);

            assert!(
                approx_eq(a_can, expected_area, 0.05),
                "touching-vertex {op:?}: canonical area={a_can}, expected={expected_area}",
            );
        }
    }

    #[test]
    fn canonical_hole_like_nested_difference_regression() {
        let outer = vec![[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0]];
        let inner = vec![[1.0, 1.0], [3.0, 1.0], [3.0, 3.0], [1.0, 3.0]];
        let expected_area = 12.0;

        let canonical = boolean_clip(&outer, &inner, ClipOp::Difference);
        let a_can = total_area(&canonical);

        assert!(
            approx_eq(a_can, expected_area, 0.1),
            "hole-like difference: canonical area={a_can}, expected={expected_area}"
        );
    }
}
