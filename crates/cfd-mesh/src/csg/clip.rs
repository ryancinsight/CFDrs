//! Sutherland–Hodgman polygon clipping with exact orientation predicates.
//!
//! Clips a convex polygon against a half-space.  The half-space is defined by
//! three CCW-oriented points `(pa, pb, pc)` on the boundary plane: a point P
//! is considered **inside** iff
//! `orient_3d(pa, pb, pc, P) ≥ 0` (i.e. `Positive` or `Degenerate`).
//!
//! Using exact [`orient_3d`] for the inside/outside decision means that no
//! false clipping occurs due to floating-point cancellation exactly at the
//! boundary.  The intersection-point computation itself uses ordinary `f64`
//! arithmetic, which is acceptable because position errors do not affect the
//! topological correctness of the Boolean result.
//!
//! ## Algorithm
//!
//! Sutherland–Hodgman for a single clipping half-space (one pass):
//!
//! ```text
//! output = []
//! for each directed edge (S → E) in the polygon:
//!   S inside,  E inside   →  append E
//!   S inside,  E outside  →  append clip(S, E)
//!   S outside, E inside   →  append clip(S, E), append E
//!   S outside, E outside  →  (skip)
//! ```
//!
//! To clip against a full convex polytope, call [`clip_polygon_to_halfplane`]
//! once for each face of the polytope.
//!
//! ## References
//!
//! Sutherland & Hodgman (1974), "Reentrant polygon clipping",
//! *Communications of the ACM*, 17(1), 32–42.

use crate::core::scalar::{Real, Point3r};
use crate::geometry::predicates::{orient_3d, Orientation};

// ── Polygon clipping ──────────────────────────────────────────────────────────

/// Clip a convex polygon against a half-space using Sutherland–Hodgman.
///
/// The half-space is defined by three points `pa`, `pb`, `pc` in
/// **counter-clockwise** order.  A point P is inside iff
/// `orient_3d(pa, pb, pc, P) ≥ 0`.
///
/// # Returns
///
/// The clipped polygon vertices, or an empty `Vec` when the entire polygon
/// lies outside.  Returns empty for inputs with fewer than 3 vertices.
///
/// # Example
///
/// ```rust
/// # use cfd_mesh::csg::clip::clip_polygon_to_halfplane;
/// # use nalgebra::Point3;
/// let square: Vec<Point3<f64>> = vec![
///     Point3::new(0.0, 0.0, 0.0),
///     Point3::new(2.0, 0.0, 0.0),
///     Point3::new(2.0, 2.0, 0.0),
///     Point3::new(0.0, 2.0, 0.0),
/// ];
/// // Clip plane z=0 with inside = z≥0 (normal +z: pa→pb→pc CCW viewed from +z).
/// let pa = Point3::new(0.0_f64, 0.0, 0.0);
/// let pb = Point3::new(1.0, 0.0, 0.0);
/// let pc = Point3::new(0.0, 1.0, 0.0);
/// let clipped = clip_polygon_to_halfplane(&square, &pa, &pb, &pc);
/// // All four square vertices lie in z=0 (Degenerate = inside) → fully kept.
/// assert_eq!(clipped.len(), 4);
/// ```
pub fn clip_polygon_to_halfplane(
    polygon: &[Point3r],
    pa: &Point3r,
    pb: &Point3r,
    pc: &Point3r,
) -> Vec<Point3r> {
    if polygon.len() < 3 {
        return Vec::new();
    }

    let plane_normal = (pb - pa).cross(&(pc - pa));

    let is_inside = |p: &Point3r| -> bool {
        let arr = |q: &Point3r| [q.x, q.y, q.z];
        orient_3d(arr(pa), arr(pb), arr(pc), arr(p)) != Orientation::Negative
    };

    // Plane-edge intersection parameter t ∈ [0,1] along S→E.
    let clip_point = |s: &Point3r, e: &Point3r| -> Point3r {
        let ds = plane_normal.dot(&(s - pa));
        let de = plane_normal.dot(&(e - pa));
        let denom = ds - de;
        if denom.abs() < 1e-20 {
            return *s; // Edge is parallel to the plane.
        }
        let t = ds / denom;
        Point3r::new(
            s.x + (e.x - s.x) * t,
            s.y + (e.y - s.y) * t,
            s.z + (e.z - s.z) * t,
        )
    };

    let n_verts = polygon.len();
    let mut output = Vec::with_capacity(n_verts);

    for i in 0..n_verts {
        let s = &polygon[i];
        let e = &polygon[(i + 1) % n_verts];
        let s_in = is_inside(s);
        let e_in = is_inside(e);

        match (s_in, e_in) {
            (true, true)   => output.push(*e),
            (true, false)  => output.push(clip_point(s, e)),
            (false, true)  => { output.push(clip_point(s, e)); output.push(*e); }
            (false, false) => {}
        }
    }

    output
}

// ── Triangle convenience wrappers ─────────────────────────────────────────────

/// Clip a triangle against a half-space and fan-triangulate the result.
///
/// Returns the sub-triangles that lie inside the half-space, or an empty
/// `Vec` if the triangle is fully clipped.
pub fn clip_triangle_to_halfplane(
    a: &Point3r, b: &Point3r, c: &Point3r,
    pa: &Point3r, pb: &Point3r, pc: &Point3r,
) -> Vec<[Point3r; 3]> {
    let polygon = [*a, *b, *c];
    let clipped = clip_polygon_to_halfplane(&polygon, pa, pb, pc);
    fan_triangulate(&clipped)
}

// ── Fan triangulation ─────────────────────────────────────────────────────────

/// Fan-triangulate a convex polygon represented as an ordered vertex list.
///
/// Returns `n − 2` triangles for an n-gon (n ≥ 3).
/// Returns an empty `Vec` when the polygon has fewer than 3 vertices.
///
/// # Theorem: Fan Triangulation is valid for Convex Polygons
/// For a convex n-gon with vertices V0 … V_{n-1}, the triangles
/// (V0, V_i, V_{i+1}) for i = 1 … n−2 are non-overlapping and cover the
/// polygon exactly.  All edges V0–V_i lie strictly inside (or on the boundary
/// of) the polygon by convexity. ∎
pub fn fan_triangulate(polygon: &[Point3r]) -> Vec<[Point3r; 3]> {
    if polygon.len() < 3 {
        return Vec::new();
    }
    let root = polygon[0];
    (1..polygon.len() - 1)
        .map(|i| [root, polygon[i], polygon[i + 1]])
        .collect()
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // CCW plane normal points in +z: pa=(0,0,0), pb=(1,0,0), pc=(0,1,0)
    // → inside = z ≥ 0.
    fn z_plane() -> (Point3r, Point3r, Point3r) {
        (
            Point3r::new(0.0, 0.0, 0.0),
            Point3r::new(1.0, 0.0, 0.0),
            Point3r::new(0.0, 1.0, 0.0),
        )
    }

    #[test]
    fn square_on_boundary_fully_kept() {
        let poly: Vec<Point3r> = vec![
            Point3r::new(0.0, 0.0, 0.0),
            Point3r::new(1.0, 0.0, 0.0),
            Point3r::new(1.0, 1.0, 0.0),
            Point3r::new(0.0, 1.0, 0.0),
        ];
        let (pa, pb, pc) = z_plane();
        let result = clip_polygon_to_halfplane(&poly, &pa, &pb, &pc);
        // All vertices are on the plane (Degenerate = inside).
        assert_eq!(result.len(), 4);
    }

    #[test]
    fn triangle_below_plane_fully_clipped() {
        let poly: Vec<Point3r> = vec![
            Point3r::new(0.0, 0.0, -1.0),
            Point3r::new(1.0, 0.0, -1.0),
            Point3r::new(0.0, 1.0, -1.0),
        ];
        let (pa, pb, pc) = z_plane();
        let result = clip_polygon_to_halfplane(&poly, &pa, &pb, &pc);
        assert!(result.is_empty(),
            "triangle fully below z=0 should be fully clipped, got {} vertices", result.len());
    }

    #[test]
    fn triangle_above_plane_fully_kept() {
        let poly: Vec<Point3r> = vec![
            Point3r::new(0.0, 0.0, 1.0),
            Point3r::new(1.0, 0.0, 1.0),
            Point3r::new(0.0, 1.0, 1.0),
        ];
        let (pa, pb, pc) = z_plane();
        let result = clip_polygon_to_halfplane(&poly, &pa, &pb, &pc);
        assert_eq!(result.len(), 3, "triangle above z=0 should be fully kept");
    }

    #[test]
    fn straddling_triangle_clipped_correctly() {
        // Triangle straddles z=0: one vertex below, two above.
        let poly: Vec<Point3r> = vec![
            Point3r::new(0.0, 0.0, -1.0), // below
            Point3r::new(2.0, 0.0,  1.0), // above
            Point3r::new(0.0, 2.0,  1.0), // above
        ];
        let (pa, pb, pc) = z_plane();
        let result = clip_polygon_to_halfplane(&poly, &pa, &pb, &pc);
        // Result should be a quadrilateral (4 vertices): the two intersection
        // points plus the two above-plane vertices.
        assert_eq!(result.len(), 4,
            "straddling triangle should produce quadrilateral, got {result:?}");
        for p in &result {
            assert!(p.z >= -1e-9, "all result vertices should be ≥ z=0, got z={}", p.z);
        }
    }

    #[test]
    fn fan_triangulate_pentagon() {
        let poly: Vec<Point3r> = (0..5)
            .map(|i| {
                let angle = i as Real * std::f64::consts::TAU / 5.0;
                Point3r::new(angle.cos(), angle.sin(), 0.0)
            })
            .collect();
        let tris = fan_triangulate(&poly);
        assert_eq!(tris.len(), 3, "pentagon → 3 triangles");
    }

    #[test]
    fn fan_triangulate_degenerate_returns_empty() {
        let poly: Vec<Point3r> = vec![
            Point3r::new(0.0, 0.0, 0.0),
            Point3r::new(1.0, 0.0, 0.0),
        ];
        let tris = fan_triangulate(&poly);
        assert!(tris.is_empty());
    }
}
