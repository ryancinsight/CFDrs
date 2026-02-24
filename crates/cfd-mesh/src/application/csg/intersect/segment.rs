//! Floating-point line intersection math.

use crate::domain::core::scalar::{Point3r, Real};
use super::types::IntersectionType;

/// Compute the 3-D intersection segment for two triangles known to straddle
/// each other's planes.
pub(crate) fn compute_segment(
    a: &Point3r,
    b: &Point3r,
    c: &Point3r,
    d: &Point3r,
    e: &Point3r,
    f: &Point3r,
) -> IntersectionType {
    let n1 = (b - a).cross(&(c - a));
    let n2 = (e - d).cross(&(f - d));

    // Direction of the intersection line (L = n1 × n2).
    let dir = n1.cross(&n2);
    if dir.norm_squared() < 1e-20 {
        // Planes are nearly parallel despite the straddling test passing.
        // Treat as coplanar.
        return IntersectionType::Coplanar;
    }

    // Choose the axis with the largest |dir| component to maximise numerical
    // stability of the 1-D projection (Möller 1997 §3).
    let abs_dir = [dir.x.abs(), dir.y.abs(), dir.z.abs()];
    let max_axis = if abs_dir[0] >= abs_dir[1] && abs_dir[0] >= abs_dir[2] {
        0
    } else if abs_dir[1] >= abs_dir[2] {
        1
    } else {
        2
    };
    let coord = |p: &Point3r| match max_axis {
        0 => p.x,
        1 => p.y,
        _ => p.z,
    };

    // Signed distances of T1 vertices to T2's plane (un-normalised).
    let da = n2.dot(&(a - d));
    let db = n2.dot(&(b - d));
    let dc = n2.dot(&(c - d));
    // Signed distances of T2 vertices to T1's plane.
    let dd = n1.dot(&(d - a));
    let de = n1.dot(&(e - a));
    let df = n1.dot(&(f - a));

    // 1-D projections onto the intersection line (max-axis coordinate).
    let (pa, pb, pc) = (coord(a), coord(b), coord(c));
    let (pd, pe, pf) = (coord(d), coord(e), coord(f));

    let seg1 = edge_crossings_interval([a, b, c], [pa, pb, pc], [da, db, dc]);
    let seg2 = edge_crossings_interval([d, e, f], [pd, pe, pf], [dd, de, df]);

    let (t1_min, t1_max, pt1_enter, pt1_leave) = match seg1 {
        Some(s) => s,
        None => return IntersectionType::None,
    };
    let (t2_min, t2_max, pt2_enter, pt2_leave) = match seg2 {
        Some(s) => s,
        None => return IntersectionType::None,
    };

    // Overlap of the two intervals.
    let t_enter = t1_min.max(t2_min);
    let t_leave = t1_max.min(t2_max);

    if t_enter > t_leave + 1e-12 {
        return IntersectionType::None;
    }

    // Select the 3-D endpoints from whichever triangle provides each bound.
    let start = if t1_min >= t2_min {
        pt1_enter
    } else {
        pt2_enter
    };
    let end = if t1_max <= t2_max {
        pt1_leave
    } else {
        pt2_leave
    };

    if (end - start).norm_squared() < 1e-20 {
        return IntersectionType::None; // Touching at a single point only.
    }

    IntersectionType::Segment { start, end }
}

/// Compute the interval `[t_min, t_max]` where a triangle's edges cross the
/// other triangle's plane, together with the 3-D crossing points.
///
/// Returns `None` for degenerate configurations where fewer than two distinct
/// crossings are found.
fn edge_crossings_interval(
    verts: [&Point3r; 3],
    projs: [Real; 3],
    dists: [Real; 3],
) -> Option<(Real, Real, Point3r, Point3r)> {
    let mut crossings: Vec<(Real, Point3r)> = Vec::with_capacity(2);

    for i in 0..3 {
        let j = (i + 1) % 3;
        let di = dists[i];
        let dj = dists[j];

        if di * dj < 0.0 {
            // Edge crosses: di and dj have strictly opposite signs.
            let denom = di - dj;
            if denom.abs() < 1e-30 {
                continue;
            }
            let t = di / denom; // parameter ∈ (0,1)
            let tp = projs[i] + (projs[j] - projs[i]) * t;
            crossings.push((tp, lerp(verts[i], verts[j], t)));
        } else if di.abs() < 1e-12 && dj.abs() > 1e-12 {
            // Vertex i lies exactly on the plane; it is a crossing point.
            crossings.push((projs[i], *verts[i]));
        }
    }

    // Deduplicate by 1-D projection (vertex-on-plane may appear from two edges).
    crossings.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap_or(std::cmp::Ordering::Equal));
    crossings.dedup_by(|x, y| (x.0 - y.0).abs() < 1e-12);

    if crossings.len() < 2 {
        return None;
    }

    let (t_min, pt_min) = crossings[0];
    let (t_max, pt_max) = crossings[crossings.len() - 1];
    Some((t_min, t_max, pt_min, pt_max))
}

/// Linear interpolation between two 3-D points.
#[inline]
fn lerp(a: &Point3r, b: &Point3r, t: Real) -> Point3r {
    Point3r::new(
        a.x + (b.x - a.x) * t,
        a.y + (b.y - a.y) * t,
        a.z + (b.z - a.z) * t,
    )
}
