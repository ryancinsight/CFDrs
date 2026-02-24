//! Geometry helpers for 2-D polygon clipping.
//!
//! Provides fundamental geometric predicates and utilities shared across
//! all clipping algorithms: signed area, convexity test, winding-number
//! point-in-polygon, and segment-segment intersection.

use crate::domain::core::scalar::Real;
use crate::domain::geometry::predicates::{orient_2d_arr, Orientation};

/// Signed area of a 2-D polygon (positive = CCW, negative = CW).
pub(crate) fn signed_area(poly: &[[Real; 2]]) -> Real {
    let n = poly.len();
    if n < 3 {
        return 0.0;
    }
    let mut sum = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        sum += poly[i][0] * poly[j][1] - poly[j][0] * poly[i][1];
    }
    sum * 0.5
}

/// Unsigned area of a 2-D polygon.
#[inline]
pub fn polygon_area(poly: &[[Real; 2]]) -> Real {
    signed_area(poly).abs()
}

/// Test if a simple polygon is convex.
#[allow(dead_code)]
pub(crate) fn is_convex(poly: &[[Real; 2]]) -> bool {
    let n = poly.len();
    if n < 3 {
        return true;
    }
    let mut sign = Orientation::Degenerate;
    for i in 0..n {
        let j = (i + 1) % n;
        let k = (i + 2) % n;
        let ori = orient_2d_arr(poly[i], poly[j], poly[k]);
        if ori == Orientation::Degenerate {
            continue;
        }
        if sign == Orientation::Degenerate {
            sign = ori;
        } else if sign != ori {
            return false;
        }
    }
    true
}

/// Ensure a polygon is in CCW winding order.
pub(crate) fn ensure_ccw(poly: &mut Vec<[Real; 2]>) {
    // We can use signed area for CCW check. Alternatively, any sequence
    // of points gives area orientation. For pure boolean predicates,
    // Shewchuk orient_2d is better but for macroscopic ordering area is fine.
    if signed_area(poly) < 0.0 {
        poly.reverse();
    }
}

/// Segment-segment intersection parameter.
/// Returns `(t, s)` where `t` is the parameter along `(p1→p2)` and `s` along `(p3→p4)`.
pub(crate) fn seg_intersect(
    p1: [Real; 2],
    p2: [Real; 2],
    p3: [Real; 2],
    p4: [Real; 2],
) -> Option<(Real, Real)> {
    let d1x = p2[0] - p1[0];
    let d1y = p2[1] - p1[1];
    let d2x = p4[0] - p3[0];
    let d2y = p4[1] - p3[1];
    let denom = d1x * d2y - d1y * d2x;
    if denom.abs() < 1e-20 {
        return None;
    }
    let dx = p3[0] - p1[0];
    let dy = p3[1] - p1[1];
    let t = (dx * d2y - dy * d2x) / denom;
    let s = (dx * d1y - dy * d1x) / denom;
    Some((t, s))
}

/// Point-in-polygon test using winding number (robust for concave polygons).
pub(crate) fn point_in_polygon(px: Real, py: Real, poly: &[[Real; 2]]) -> bool {
    let n = poly.len();
    if n < 3 {
        return false;
    }
    let p = [px, py];
    let mut winding = 0i32;
    for i in 0..n {
        let j = (i + 1) % n;
        let yi = poly[i][1];
        let yj = poly[j][1];
        if yi <= py {
            if yj > py {
                if orient_2d_arr(poly[i], poly[j], p) == Orientation::Positive {
                    winding += 1;
                }
            }
        } else if yj <= py {
            if orient_2d_arr(poly[i], poly[j], p) == Orientation::Negative {
                winding -= 1;
            }
        }
    }
    winding != 0
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(a: Real, b: Real, tol: Real) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn test_signed_area_ccw_triangle() {
        let tri = vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let area = signed_area(&tri);
        assert!(area > 0.0, "CCW triangle should have positive signed area");
        assert!(approx_eq(area, 0.5, 1e-12));
    }

    #[test]
    fn test_signed_area_cw_triangle() {
        let tri = vec![[0.0, 0.0], [0.0, 1.0], [1.0, 0.0]];
        let area = signed_area(&tri);
        assert!(area < 0.0, "CW triangle should have negative signed area");
    }

    #[test]
    fn test_is_convex_square() {
        let sq = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        assert!(is_convex(&sq));
    }

    #[test]
    fn test_is_convex_l_shape() {
        let l = vec![
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 1.0],
            [1.0, 1.0],
            [1.0, 2.0],
            [0.0, 2.0],
        ];
        assert!(!is_convex(&l));
    }

    #[test]
    fn test_point_in_polygon_inside() {
        let sq = vec![[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
        assert!(point_in_polygon(1.0, 1.0, &sq));
    }

    #[test]
    fn test_point_in_polygon_outside() {
        let sq = vec![[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
        assert!(!point_in_polygon(3.0, 1.0, &sq));
    }
}
