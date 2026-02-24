//! Axis-Aligned Bounding Box — generic over scalar precision `T: Scalar`.

use crate::domain::core::scalar::Scalar;
use nalgebra::Point3;

/// Axis-aligned bounding box generic over scalar precision.
///
/// The default precision is `f64` so that all existing code that writes
/// `Aabb::new(...)` or `Aabb::empty()` without type annotations continues
/// to compile unchanged.  New code may instantiate `Aabb<f32>` for
/// GPU-staging geometry without any feature flag.
#[derive(Clone, Copy, Debug)]
pub struct Aabb<T: Scalar = f64> {
    /// Minimum corner.
    pub min: Point3<T>,
    /// Maximum corner.
    pub max: Point3<T>,
}

impl<T: Scalar> Aabb<T> {
    /// Create an AABB from explicit min and max corners.
    pub fn new(min: Point3<T>, max: Point3<T>) -> Self {
        Self { min, max }
    }

    /// Create an inverted (empty) AABB — grows to contain the first [`expand`]ed point.
    pub fn empty() -> Self {
        Self {
            min: Point3::new(T::infinity(), T::infinity(), T::infinity()),
            max: Point3::new(T::neg_infinity(), T::neg_infinity(), T::neg_infinity()),
        }
    }

    /// Expand to include `p`.
    pub fn expand(&mut self, p: &Point3<T>) {
        if p.x < self.min.x {
            self.min.x = p.x;
        }
        if p.y < self.min.y {
            self.min.y = p.y;
        }
        if p.z < self.min.z {
            self.min.z = p.z;
        }
        if p.x > self.max.x {
            self.max.x = p.x;
        }
        if p.y > self.max.y {
            self.max.y = p.y;
        }
        if p.z > self.max.z {
            self.max.z = p.z;
        }
    }

    /// Union with another AABB.
    pub fn union(&self, other: &Aabb<T>) -> Aabb<T> {
        Aabb {
            min: Point3::new(
                if self.min.x < other.min.x {
                    self.min.x
                } else {
                    other.min.x
                },
                if self.min.y < other.min.y {
                    self.min.y
                } else {
                    other.min.y
                },
                if self.min.z < other.min.z {
                    self.min.z
                } else {
                    other.min.z
                },
            ),
            max: Point3::new(
                if self.max.x > other.max.x {
                    self.max.x
                } else {
                    other.max.x
                },
                if self.max.y > other.max.y {
                    self.max.y
                } else {
                    other.max.y
                },
                if self.max.z > other.max.z {
                    self.max.z
                } else {
                    other.max.z
                },
            ),
        }
    }

    /// Check if two AABBs overlap (inclusive of boundary).
    pub fn intersects(&self, other: &Aabb<T>) -> bool {
        self.min.x <= other.max.x
            && self.max.x >= other.min.x
            && self.min.y <= other.max.y
            && self.max.y >= other.min.y
            && self.min.z <= other.max.z
            && self.max.z >= other.min.z
    }

    /// Check if this AABB completely contains another AABB.
    pub fn contains_aabb(&self, other: &Aabb<T>) -> bool {
        self.min.x <= other.min.x
            && self.max.x >= other.max.x
            && self.min.y <= other.min.y
            && self.max.y >= other.max.y
            && self.min.z <= other.min.z
            && self.max.z >= other.max.z
    }

    /// Check if `p` lies inside or on the boundary of this AABB.
    pub fn contains_point(&self, p: &Point3<T>) -> bool {
        p.x >= self.min.x
            && p.x <= self.max.x
            && p.y >= self.min.y
            && p.y <= self.max.y
            && p.z >= self.min.z
            && p.z <= self.max.z
    }

    /// (Width, height, depth) of the box.
    pub fn extents(&self) -> (T, T, T) {
        (
            self.max.x - self.min.x,
            self.max.y - self.min.y,
            self.max.z - self.min.z,
        )
    }

    /// Centre point.
    pub fn center(&self) -> Point3<T> {
        let half = <T as Scalar>::from_f64(0.5);
        Point3::new(
            (self.min.x + self.max.x) * half,
            (self.min.y + self.max.y) * half,
            (self.min.z + self.max.z) * half,
        )
    }

    /// Enclosed volume.
    pub fn volume(&self) -> T {
        let (w, h, d) = self.extents();
        w * h * d
    }

    /// Build an AABB from an iterator of points.
    pub fn from_points<'a>(points: impl Iterator<Item = &'a Point3<T>>) -> Self
    where
        T: 'a,
    {
        let mut aabb = Self::empty();
        for p in points {
            aabb.expand(p);
        }
        aabb
    }
}

impl<T: Scalar> Default for Aabb<T> {
    fn default() -> Self {
        Self::empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Point3;

    #[test]
    fn empty_aabb_grows_on_expand_f64() {
        let mut b = Aabb::<f64>::empty();
        b.expand(&Point3::new(1.0, 2.0, 3.0));
        b.expand(&Point3::new(-1.0, 0.0, 5.0));
        assert_eq!(b.min, Point3::new(-1.0, 0.0, 3.0));
        assert_eq!(b.max, Point3::new(1.0, 2.0, 5.0));
    }

    #[test]
    fn empty_aabb_grows_on_expand_f32() {
        let mut b = Aabb::<f32>::empty();
        b.expand(&Point3::new(1.0_f32, 2.0, 3.0));
        b.expand(&Point3::new(-1.0_f32, 0.0, 5.0));
        assert_eq!(b.min, Point3::new(-1.0_f32, 0.0, 3.0));
        assert_eq!(b.max, Point3::new(1.0_f32, 2.0, 5.0));
    }

    #[test]
    fn center_f64() {
        let b = Aabb::<f64>::new(Point3::new(0.0, 0.0, 0.0), Point3::new(2.0, 4.0, 6.0));
        let c = b.center();
        assert!((c.x - 1.0).abs() < 1e-12);
        assert!((c.y - 2.0).abs() < 1e-12);
        assert!((c.z - 3.0).abs() < 1e-12);
    }

    #[test]
    fn volume_f32() {
        let b = Aabb::<f32>::new(
            Point3::new(0.0_f32, 0.0, 0.0),
            Point3::new(2.0_f32, 3.0, 4.0),
        );
        assert!((b.volume() - 24.0_f32).abs() < 1e-5);
    }

    #[test]
    fn intersects_and_disjoint() {
        let a = Aabb::<f64>::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0));
        let b = Aabb::<f64>::new(Point3::new(0.5, 0.5, 0.5), Point3::new(1.5, 1.5, 1.5));
        let c = Aabb::<f64>::new(Point3::new(2.0, 2.0, 2.0), Point3::new(3.0, 3.0, 3.0));
        assert!(a.intersects(&b));
        assert!(!a.intersects(&c));
    }
}
