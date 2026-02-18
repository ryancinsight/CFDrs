//! Axis-Aligned Bounding Box.

use crate::core::scalar::{Real, Point3r};

/// Axis-aligned bounding box.
#[derive(Clone, Copy, Debug)]
pub struct Aabb {
    /// Minimum corner.
    pub min: Point3r,
    /// Maximum corner.
    pub max: Point3r,
}

impl Aabb {
    /// Create an AABB from min and max corners.
    pub fn new(min: Point3r, max: Point3r) -> Self {
        Self { min, max }
    }

    /// Create an empty (inverted) AABB that will grow on first `expand`.
    pub fn empty() -> Self {
        Self {
            min: Point3r::new(Real::INFINITY, Real::INFINITY, Real::INFINITY),
            max: Point3r::new(Real::NEG_INFINITY, Real::NEG_INFINITY, Real::NEG_INFINITY),
        }
    }

    /// Expand to include a point.
    pub fn expand(&mut self, p: &Point3r) {
        self.min.x = self.min.x.min(p.x);
        self.min.y = self.min.y.min(p.y);
        self.min.z = self.min.z.min(p.z);
        self.max.x = self.max.x.max(p.x);
        self.max.y = self.max.y.max(p.y);
        self.max.z = self.max.z.max(p.z);
    }

    /// Expand to include another AABB.
    pub fn union(&self, other: &Aabb) -> Aabb {
        Aabb {
            min: Point3r::new(
                self.min.x.min(other.min.x),
                self.min.y.min(other.min.y),
                self.min.z.min(other.min.z),
            ),
            max: Point3r::new(
                self.max.x.max(other.max.x),
                self.max.y.max(other.max.y),
                self.max.z.max(other.max.z),
            ),
        }
    }

    /// Check if two AABBs overlap.
    pub fn intersects(&self, other: &Aabb) -> bool {
        self.min.x <= other.max.x
            && self.max.x >= other.min.x
            && self.min.y <= other.max.y
            && self.max.y >= other.min.y
            && self.min.z <= other.max.z
            && self.max.z >= other.min.z
    }

    /// Check if a point is inside (or on the boundary of) this AABB.
    pub fn contains_point(&self, p: &Point3r) -> bool {
        p.x >= self.min.x
            && p.x <= self.max.x
            && p.y >= self.min.y
            && p.y <= self.max.y
            && p.z >= self.min.z
            && p.z <= self.max.z
    }

    /// Dimensions (width, height, depth).
    pub fn extents(&self) -> (Real, Real, Real) {
        (
            self.max.x - self.min.x,
            self.max.y - self.min.y,
            self.max.z - self.min.z,
        )
    }

    /// Center point.
    pub fn center(&self) -> Point3r {
        Point3r::new(
            (self.min.x + self.max.x) * 0.5,
            (self.min.y + self.max.y) * 0.5,
            (self.min.z + self.max.z) * 0.5,
        )
    }

    /// Volume.
    pub fn volume(&self) -> Real {
        let (w, h, d) = self.extents();
        w * h * d
    }

    /// Build AABB from an iterator of points.
    pub fn from_points<'a>(points: impl Iterator<Item = &'a Point3r>) -> Self {
        let mut aabb = Self::empty();
        for p in points {
            aabb.expand(p);
        }
        aabb
    }
}

impl Default for Aabb {
    fn default() -> Self {
        Self::empty()
    }
}
