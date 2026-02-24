//! Sutherland-Hodgman convex polygon clipping.
//!
//! Optimal for convex clip regions: O(n) per clipping edge, O(n·m) total.
//! Only supports convex clip polygons — for concave clips, use the hybrid
//! dispatcher which selects an appropriate general algorithm.
//!
//! # References
//!
//! Sutherland & Hodgman (1974), "Reentrant polygon clipping"

use crate::domain::core::scalar::Real;

/// Evaluates the 2-D cross product (unscaled signed distance).
/// MUST be evaluated in standard floats to construct precise `t` interpolations.
#[inline]
fn edge_distance(ax: Real, ay: Real, bx: Real, by: Real, px: Real, py: Real) -> Real {
    (bx - ax) * (py - ay) - (by - ay) * (px - ax)
}

/// Clip a polygon against the left half-plane of directed edge (ax,ay)→(bx,by).
///
/// Retained from the original implementation for backward compatibility.
/// Optimal for convex clip regions (one pass per edge, O(n) total).
pub fn sh_clip_halfplane(
    poly: &[[Real; 2]],
    ax: Real,
    ay: Real,
    bx: Real,
    by: Real,
) -> Vec<[Real; 2]> {
    if poly.len() < 2 {
        return Vec::new();
    }
    let n = poly.len();
    let mut out = Vec::with_capacity(n + 1);
    for i in 0..n {
        let s = poly[i];
        let e = poly[(i + 1) % n];
        let sc = edge_distance(ax, ay, bx, by, s[0], s[1]);
        let ec = edge_distance(ax, ay, bx, by, e[0], e[1]);

        // Use the same float values for inside checking to perfectly synchronize
        // with the numeric interpolation branching.
        let s_in = sc >= 0.0;
        let e_in = ec >= 0.0;
        match (s_in, e_in) {
            (true, true) => out.push(e),
            (true, false) => {
                let denom = sc - ec;
                if denom.abs() > 1e-30 {
                    let t = sc / denom;
                    out.push([s[0] + (e[0] - s[0]) * t, s[1] + (e[1] - s[1]) * t]);
                }
            }
            (false, true) => {
                let denom = sc - ec;
                if denom.abs() > 1e-30 {
                    let t = sc / denom;
                    out.push([s[0] + (e[0] - s[0]) * t, s[1] + (e[1] - s[1]) * t]);
                }
                out.push(e);
            }
            (false, false) => {}
        }
    }
    out
}

/// Clip subject polygon to the inside of a CCW convex clip polygon
/// using iterated Sutherland-Hodgman half-plane clips.
pub fn sh_clip_convex(subject: &[[Real; 2]], clip: &[[Real; 2]]) -> Vec<[Real; 2]> {
    let n = clip.len();
    if n < 3 || subject.len() < 3 {
        return Vec::new();
    }
    let mut result = subject.to_vec();
    for i in 0..n {
        let j = (i + 1) % n;
        result = sh_clip_halfplane(&result, clip[i][0], clip[i][1], clip[j][0], clip[j][1]);
        if result.len() < 3 {
            return Vec::new();
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::super::geometry::polygon_area;
    use super::*;

    fn approx_eq(a: Real, b: Real, tol: Real) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn sh_clip_two_squares_intersection() {
        let subject = vec![[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
        let clip = vec![[1.0, 1.0], [3.0, 1.0], [3.0, 3.0], [1.0, 3.0]];
        let result = sh_clip_convex(&subject, &clip);
        assert!(result.len() >= 3);
        let area = polygon_area(&result);
        assert!(
            approx_eq(area, 1.0, 0.01),
            "intersection area should be 1.0, got {area}"
        );
    }

    #[test]
    fn sh_clip_triangle_inside_square() {
        let tri = vec![[0.5, 0.5], [1.5, 0.5], [1.0, 1.5]];
        let sq = vec![[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
        let result = sh_clip_convex(&tri, &sq);
        assert!(result.len() >= 3);
        let area = polygon_area(&result);
        let tri_area = polygon_area(&tri);
        assert!(
            approx_eq(area, tri_area, 1e-10),
            "triangle fully inside square"
        );
    }
}
