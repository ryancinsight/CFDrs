//! Geometry helpers: surface area, axis selection, AABB extents.

use crate::domain::geometry::aabb::Aabb;

/// Half surface area of an AABB, proportional to full surface area.
///
/// Avoids the ×2 factor throughout — the SAH cost ratio only requires relative
/// surface areas, so the constant factor cancels.
#[inline]
pub(super) fn surface_area(aabb: &Aabb) -> f64 {
    let (w, h, d) = aabb.extents();
    w * h + h * d + d * w
}

/// Extent of an AABB along one axis (0 = X, 1 = Y, 2 = Z).
#[inline]
pub(super) fn axis_extent(aabb: &Aabb, axis: usize) -> f64 {
    let (w, h, d) = aabb.extents();
    match axis {
        0 => w,
        1 => h,
        _ => d,
    }
}

/// Minimum coordinate of an AABB along one axis.
#[inline]
pub(super) fn axis_min(aabb: &Aabb, axis: usize) -> f64 {
    match axis {
        0 => aabb.min.x,
        1 => aabb.min.y,
        _ => aabb.min.z,
    }
}

/// Coordinate of a point along one axis.
#[inline]
pub(super) fn axis_value(p: &nalgebra::Point3<f64>, axis: usize) -> f64 {
    match axis {
        0 => p.x,
        1 => p.y,
        _ => p.z,
    }
}

/// Return `(axis, extent)` for the longest axis of an AABB.
#[inline]
pub(super) fn longest_axis(aabb: &Aabb) -> (usize, f64) {
    let (w, h, d) = aabb.extents();
    if w >= h && w >= d {
        (0, w)
    } else if h >= d {
        (1, h)
    } else {
        (2, d)
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
pub(super) mod tests {
    use super::*;
    use nalgebra::Point3;

    fn pt(x: f64, y: f64, z: f64) -> Point3<f64> {
        Point3::new(x, y, z)
    }

    #[test]
    fn surface_area_unit_cube() {
        let a = Aabb::new(pt(0.0, 0.0, 0.0), pt(1.0, 1.0, 1.0));
        // Half SA of a unit cube: 3 × (1×1) = 3.
        assert!((surface_area(&a) - 3.0).abs() < 1e-12);
    }

    #[test]
    fn longest_axis_selects_x() {
        let a = Aabb::new(pt(0.0, 0.0, 0.0), pt(10.0, 1.0, 1.0));
        let (axis, extent) = longest_axis(&a);
        assert_eq!(axis, 0);
        assert!((extent - 10.0).abs() < 1e-12);
    }
}
