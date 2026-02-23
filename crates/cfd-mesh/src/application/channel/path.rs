//! Channel centerline paths.

use crate::core::scalar::{Real, Point3r, Vector3r};

/// A channel centerline path defined by ordered waypoints.
#[derive(Clone, Debug)]
pub struct ChannelPath {
    /// Ordered waypoints (3D positions of the channel centerline).
    points: Vec<Point3r>,
}

impl ChannelPath {
    /// Create a path from a set of waypoints.
    pub fn new(points: Vec<Point3r>) -> Self {
        assert!(points.len() >= 2, "path must have at least 2 points");
        Self { points }
    }

    /// Create a straight-line path between two points.
    pub fn straight(start: Point3r, end: Point3r) -> Self {
        Self {
            points: vec![start, end],
        }
    }

    /// Get the waypoints.
    pub fn points(&self) -> &[Point3r] {
        &self.points
    }

    /// Number of segments.
    pub fn segment_count(&self) -> usize {
        self.points.len() - 1
    }

    /// Total path length.
    pub fn length(&self) -> Real {
        self.points
            .windows(2)
            .map(|w| (w[1] - w[0]).norm())
            .sum()
    }

    /// Direction at a given segment index (normalized).
    pub fn segment_direction(&self, segment: usize) -> Vector3r {
        let dir = self.points[segment + 1] - self.points[segment];
        dir.normalize()
    }

    /// Compute a local frame (tangent, normal, binormal) at each waypoint
    /// using the Frenet–Serret formulas with fallback for straight segments.
    pub fn compute_frames(&self) -> Vec<FrenetFrame> {
        let n = self.points.len();
        let mut frames = Vec::with_capacity(n);

        for i in 0..n {
            let tangent = if i == 0 {
                (self.points[1] - self.points[0]).normalize()
            } else if i == n - 1 {
                (self.points[n - 1] - self.points[n - 2]).normalize()
            } else {
                ((self.points[i + 1] - self.points[i]).normalize()
                    + (self.points[i] - self.points[i - 1]).normalize())
                .normalize()
            };

            // Find a normal perpendicular to tangent
            let normal = find_perpendicular(&tangent);
            let binormal = tangent.cross(&normal).normalize();

            frames.push(FrenetFrame {
                position: self.points[i],
                tangent,
                normal,
                binormal,
            });
        }

        frames
    }
}

/// Local coordinate frame at a point along the path.
#[derive(Clone, Debug)]
pub struct FrenetFrame {
    /// Position on the path.
    pub position: Point3r,
    /// Tangent direction (forward along path).
    pub tangent: Vector3r,
    /// Normal direction (perpendicular to tangent).
    pub normal: Vector3r,
    /// Binormal direction (tangent × normal).
    pub binormal: Vector3r,
}

/// Find a vector perpendicular to `v`.
fn find_perpendicular(v: &Vector3r) -> Vector3r {
    let candidate = if v.x.abs() < 0.9 {
        Vector3r::x()
    } else {
        Vector3r::y()
    };
    let perp = v.cross(&candidate);
    perp.normalize()
}
