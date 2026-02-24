//! Per-triangle quality measurements.

use crate::domain::core::constants;
use crate::domain::core::scalar::{Point3r, Real};

/// Triangle aspect ratio: longest_edge / shortest_altitude.
///
/// A perfect equilateral triangle has aspect ratio ~1.155.
/// CFD meshes typically want aspect ratio < 5.
pub fn aspect_ratio(a: &Point3r, b: &Point3r, c: &Point3r) -> Real {
    let edges = [(b - a).norm(), (c - b).norm(), (a - c).norm()];
    let longest = edges.iter().copied().fold(0.0, Real::max);
    let area = 0.5 * (b - a).cross(&(c - a)).norm();

    if area < Real::EPSILON {
        return Real::INFINITY; // Degenerate
    }

    let shortest_alt = 2.0 * area / longest;
    longest / shortest_alt
}

/// Minimum interior angle of a triangle (in radians).
pub fn min_angle(a: &Point3r, b: &Point3r, c: &Point3r) -> Real {
    let ab = (b - a).normalize();
    let ac = (c - a).normalize();
    let ba = (a - b).normalize();
    let bc = (c - b).normalize();
    let ca = (a - c).normalize();
    let cb = (b - c).normalize();

    let angle_a = ab.dot(&ac).clamp(-1.0, 1.0).acos();
    let angle_b = ba.dot(&bc).clamp(-1.0, 1.0).acos();
    let angle_c = ca.dot(&cb).clamp(-1.0, 1.0).acos();

    angle_a.min(angle_b).min(angle_c)
}

/// Maximum interior angle of a triangle (in radians).
pub fn max_angle(a: &Point3r, b: &Point3r, c: &Point3r) -> Real {
    let ab = (b - a).normalize();
    let ac = (c - a).normalize();
    let ba = (a - b).normalize();
    let bc = (c - b).normalize();
    let ca = (a - c).normalize();
    let cb = (b - c).normalize();

    let angle_a = ab.dot(&ac).clamp(-1.0, 1.0).acos();
    let angle_b = ba.dot(&bc).clamp(-1.0, 1.0).acos();
    let angle_c = ca.dot(&cb).clamp(-1.0, 1.0).acos();

    angle_a.max(angle_b).max(angle_c)
}

/// Equiangle skewness: measures deviation from equilateral.
///
/// 0 = equilateral, 1 = degenerate.
pub fn equiangle_skewness(a: &Point3r, b: &Point3r, c: &Point3r) -> Real {
    let ideal = constants::PI / 3.0; // 60° for equilateral triangle
    let max_a = max_angle(a, b, c);
    let min_a = min_angle(a, b, c);

    let skew_max = (max_a - ideal) / (constants::PI - ideal);
    let skew_min = (ideal - min_a) / ideal;

    skew_max.max(skew_min)
}

/// Edge length ratio: shortest_edge / longest_edge.
///
/// 1.0 = all edges equal, approaching 0 = bad quality.
pub fn edge_length_ratio(a: &Point3r, b: &Point3r, c: &Point3r) -> Real {
    let edges = [(b - a).norm(), (c - b).norm(), (a - c).norm()];
    let shortest = edges.iter().copied().fold(Real::INFINITY, Real::min);
    let longest = edges.iter().copied().fold(0.0, Real::max);

    if longest < Real::EPSILON {
        return 0.0;
    }

    shortest / longest
}
