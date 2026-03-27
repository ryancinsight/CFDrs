//! Ray casting — CPU-based picking fallback using ray-geometry intersection.
//!
//! # Theorem — Moller-Trumbore Intersection
//!
//! The Moller-Trumbore algorithm solves the ray-triangle intersection by
//! expressing the intersection point as `O + tD = (1-u-v)V0 + uV1 + vV2`
//! and solving for `(t, u, v)` via Cramer's rule. The intersection is valid
//! when `t > 0`, `u >= 0`, `v >= 0`, and `u + v <= 1`.  ∎

use nalgebra::{Point3, Vector3};

/// A ray in 3D space defined by an origin and direction.
#[derive(Clone, Debug)]
pub struct Ray3 {
    pub origin: Point3<f64>,
    pub direction: Vector3<f64>,
}

impl Ray3 {
    /// Create a new ray. The direction should be normalized for correct `t` values.
    #[must_use]
    pub fn new(origin: Point3<f64>, direction: Vector3<f64>) -> Self {
        Self { origin, direction }
    }

    /// Compute a point along the ray at parameter `t`.
    #[must_use]
    pub fn at(&self, t: f64) -> Point3<f64> {
        self.origin + self.direction * t
    }
}

/// Moller-Trumbore ray-triangle intersection.
///
/// Returns the parameter `t` along the ray if the ray hits the triangle.
/// The triangle vertices must be in counter-clockwise winding order for
/// a front-face hit (positive `t`).
#[must_use]
pub fn ray_triangle_intersection(
    ray: &Ray3,
    v0: &Point3<f64>,
    v1: &Point3<f64>,
    v2: &Point3<f64>,
) -> Option<f64> {
    let edge1 = v1 - v0;
    let edge2 = v2 - v0;
    let h = ray.direction.cross(&edge2);
    let a = edge1.dot(&h);

    if a.abs() < EPSILON {
        return None; // Parallel to the triangle.
    }

    let f = 1.0 / a;
    let s = ray.origin - v0;
    let u = f * s.dot(&h);
    if !(0.0..=1.0).contains(&u) {
        return None;
    }

    let q = s.cross(&edge1);
    let v = f * ray.direction.dot(&q);
    if v < 0.0 || u + v > 1.0 {
        return None;
    }

    let t = f * edge2.dot(&q);
    if t > EPSILON {
        Some(t)
    } else {
        None // Behind the ray origin.
    }
}

/// Ray-AABB intersection using the slab method.
///
/// Returns `(t_near, t_far)` if the ray intersects the box, or `None`
/// if it misses entirely.
#[must_use]
pub fn ray_aabb_intersection(
    ray: &Ray3,
    aabb_min: &Point3<f64>,
    aabb_max: &Point3<f64>,
) -> Option<(f64, f64)> {
    let mut t_near = f64::NEG_INFINITY;
    let mut t_far = f64::INFINITY;

    for i in 0..3 {
        let inv_d = 1.0 / ray.direction[i];
        let mut t0 = (aabb_min[i] - ray.origin[i]) * inv_d;
        let mut t1 = (aabb_max[i] - ray.origin[i]) * inv_d;
        if inv_d < 0.0 {
            std::mem::swap(&mut t0, &mut t1);
        }
        t_near = t_near.max(t0);
        t_far = t_far.min(t1);
        if t_near > t_far {
            return None;
        }
    }

    if t_far < 0.0 {
        return None; // Box is entirely behind the ray.
    }

    Some((t_near, t_far))
}

const EPSILON: f64 = 1e-12;

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn ray_hits_triangle() {
        let ray = Ray3::new(Point3::new(0.25, 0.25, -1.0), Vector3::z());
        let v0 = Point3::new(0.0, 0.0, 0.0);
        let v1 = Point3::new(1.0, 0.0, 0.0);
        let v2 = Point3::new(0.0, 1.0, 0.0);
        let t = ray_triangle_intersection(&ray, &v0, &v1, &v2).expect("should hit");
        assert_relative_eq!(t, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn ray_misses_triangle() {
        let ray = Ray3::new(Point3::new(2.0, 2.0, -1.0), Vector3::z());
        let v0 = Point3::new(0.0, 0.0, 0.0);
        let v1 = Point3::new(1.0, 0.0, 0.0);
        let v2 = Point3::new(0.0, 1.0, 0.0);
        assert!(ray_triangle_intersection(&ray, &v0, &v1, &v2).is_none());
    }

    #[test]
    fn ray_hits_aabb() {
        let ray = Ray3::new(Point3::new(0.5, 0.5, -5.0), Vector3::z());
        let min = Point3::new(0.0, 0.0, 0.0);
        let max = Point3::new(1.0, 1.0, 1.0);
        let (t_near, t_far) = ray_aabb_intersection(&ray, &min, &max).expect("should hit");
        assert_relative_eq!(t_near, 5.0, epsilon = 1e-10);
        assert_relative_eq!(t_far, 6.0, epsilon = 1e-10);
    }

    #[test]
    fn ray_misses_aabb() {
        let ray = Ray3::new(Point3::new(5.0, 5.0, -5.0), Vector3::z());
        let min = Point3::new(0.0, 0.0, 0.0);
        let max = Point3::new(1.0, 1.0, 1.0);
        assert!(ray_aabb_intersection(&ray, &min, &max).is_none());
    }
}
