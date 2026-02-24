//! 2-D point and AABB helpers.

use crate::domain::core::scalar::Real;
use crate::domain::geometry::predicates::{orient_2d_arr, Orientation};

/// Test whether 2-D point `(px,py)` lies inside or on the boundary of the
/// CCW-wound triangle `(ax,ay)→(bx,by)→(cx,cy)` using exact arithmetic.
#[inline]
pub(crate) fn point_in_tri_2d_exact(
    px: Real,
    py: Real,
    ax: Real,
    ay: Real,
    bx: Real,
    by: Real,
    cx: Real,
    cy: Real,
) -> bool {
    let p = [px, py];
    let a = [ax, ay];
    let b = [bx, by];
    let c = [cx, cy];

    let d0 = orient_2d_arr(a, b, p);
    let d1 = orient_2d_arr(b, c, p);
    let d2 = orient_2d_arr(c, a, p);

    let neg =
        d0 == Orientation::Negative || d1 == Orientation::Negative || d2 == Orientation::Negative;
    let pos =
        d0 == Orientation::Positive || d1 == Orientation::Positive || d2 == Orientation::Positive;

    // Inside or on boundary iff it's not strictly outside on some edges and strictly inside on others.
    !(neg && pos)
}

/// Test whether 2-D point is inside the union of a set of triangles
/// (given as flat `[ax,ay,bx,by,cx,cy]` arrays).
#[inline]
pub(crate) fn point_in_union_2d_exact(px: Real, py: Real, tris: &[[Real; 6]]) -> bool {
    tris.iter()
        .any(|t| point_in_tri_2d_exact(px, py, t[0], t[1], t[2], t[3], t[4], t[5]))
}

/// 2-D AABB of a triangle: `[min_u, min_v, max_u, max_v]`.
#[inline]
pub(crate) fn aabb2(ax: Real, ay: Real, bx: Real, by: Real, cx: Real, cy: Real) -> [Real; 4] {
    [
        ax.min(bx).min(cx),
        ay.min(by).min(cy),
        ax.max(bx).max(cx),
        ay.max(by).max(cy),
    ]
}

/// True if two 2-D AABBs intersect (inclusive boundary).
#[inline]
pub(crate) fn aabb_overlaps(a: &[Real; 4], b: &[Real; 4]) -> bool {
    a[0] <= b[2] && b[0] <= a[2] && a[1] <= b[3] && b[1] <= a[3]
}

/// Unsigned area of a 2-D polygon via the shoelace formula.
#[cfg(test)]
#[inline]
pub(crate) fn polygon_area_2d(poly: &[[Real; 2]]) -> Real {
    let n = poly.len();
    if n < 3 {
        return 0.0;
    }
    let mut sum = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        sum += poly[i][0] * poly[j][1] - poly[j][0] * poly[i][1];
    }
    sum.abs() * 0.5
}
