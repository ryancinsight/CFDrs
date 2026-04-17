//! Hidden-line removal — depth-based classification of projected edges.
//!
//! This module provides a software depth buffer approach for classifying
//! projected 2D edges as visible or hidden based on the 3D mesh geometry.

use nalgebra::{Point2, Point3, Vector3};

/// Simple software depth buffer for hidden-line classification.
pub struct DepthClassifier {
    width: usize,
    height: usize,
    depth_buffer: Vec<f64>,
    scale: f64,
    offset_x: f64,
    offset_y: f64,
}

impl DepthClassifier {
    /// Create a depth buffer at the given resolution.
    ///
    /// The `bounds` define the 2D bounding box of the projected geometry
    /// as (`min_x`, `min_y`, `max_x`, `max_y`).
    #[must_use]
    pub fn new(width: usize, height: usize, bounds: (f64, f64, f64, f64)) -> Self {
        let (min_x, min_y, max_x, max_y) = bounds;
        let range_x = max_x - min_x;
        let range_y = max_y - min_y;
        let scale = if range_x > range_y {
            (width as f64) / range_x
        } else {
            (height as f64) / range_y
        };
        Self {
            width,
            height,
            depth_buffer: vec![f64::INFINITY; width * height],
            scale,
            offset_x: -min_x,
            offset_y: -min_y,
        }
    }

    /// Rasterize a triangle into the depth buffer.
    pub fn rasterize_triangle(
        &mut self,
        p0: &Point3<f64>,
        p1: &Point3<f64>,
        p2: &Point3<f64>,
        view_dir: &Vector3<f64>,
    ) {
        // Project to 2D + depth (depth = dot with view direction).
        let d0 = p0.coords.dot(view_dir);
        let d1 = p1.coords.dot(view_dir);
        let d2 = p2.coords.dot(view_dir);

        let s0 = self.to_screen(p0.x, p0.y);
        let s1 = self.to_screen(p1.x, p1.y);
        let s2 = self.to_screen(p2.x, p2.y);

        // Compute bounding box in screen space.
        let min_x = s0.0.min(s1.0).min(s2.0).max(0) as usize;
        let max_x = (s0.0.max(s1.0).max(s2.0) as usize).min(self.width - 1);
        let min_y = s0.1.min(s1.1).min(s2.1).max(0) as usize;
        let max_y = (s0.1.max(s1.1).max(s2.1) as usize).min(self.height - 1);

        for y in min_y..=max_y {
            for x in min_x..=max_x {
                let px = x as f64 + 0.5;
                let py = y as f64 + 0.5;
                if let Some((u, v, w)) = barycentric(
                    &TriPoint { px, py },
                    &Triangle {
                        ax: f64::from(s0.0), ay: f64::from(s0.1),
                        bx: f64::from(s1.0), by: f64::from(s1.1),
                        cx: f64::from(s2.0), cy: f64::from(s2.1),
                    },
                ) {
                    let depth = u * d0 + v * d1 + w * d2;
                    let idx = y * self.width + x;
                    if depth < self.depth_buffer[idx] {
                        self.depth_buffer[idx] = depth;
                    }
                }
            }
        }
    }

    /// Test if a 2D point at the given depth is visible (in front of the depth buffer).
    #[must_use]
    pub fn is_visible(&self, p: &Point2<f64>, depth: f64) -> bool {
        let (sx, sy) = self.to_screen(p.x, p.y);
        let x = sx as usize;
        let y = sy as usize;
        if x >= self.width || y >= self.height {
            return true;
        }
        let idx = y * self.width + x;
        depth <= self.depth_buffer[idx] + 1e-6
    }

    fn to_screen(&self, x: f64, y: f64) -> (i32, i32) {
        (
            ((x + self.offset_x) * self.scale) as i32,
            ((y + self.offset_y) * self.scale) as i32,
        )
    }
}

/// Barycentric coordinates of a point relative to a triangle.
struct TriPoint {
    px: f64,
    py: f64,
}

struct Triangle {
    ax: f64,
    ay: f64,
    bx: f64,
    by: f64,
    cx: f64,
    cy: f64,
}

/// Compute barycentric coordinates of point in triangle.
/// Returns None if point is outside the triangle.
fn barycentric(pt: &TriPoint, tri: &Triangle) -> Option<(f64, f64, f64)> {
    let v0x = tri.bx - tri.ax;
    let v0y = tri.by - tri.ay;
    let v1x = tri.cx - tri.ax;
    let v1y = tri.cy - tri.ay;
    let v2x = pt.px - tri.ax;
    let v2y = pt.py - tri.ay;

    let d00 = v0x * v0x + v0y * v0y;
    let d01 = v0x * v1x + v0y * v1y;
    let d11 = v1x * v1x + v1y * v1y;
    let d20 = v2x * v0x + v2y * v0y;
    let d21 = v2x * v1x + v2y * v1y;

    let denom = d00 * d11 - d01 * d01;
    if denom.abs() < 1e-12 {
        return None;
    }

    let v = (d11 * d20 - d01 * d21) / denom;
    let w = (d00 * d21 - d01 * d20) / denom;
    let u = 1.0 - v - w;

    if u >= 0.0 && v >= 0.0 && w >= 0.0 {
        Some((u, v, w))
    } else {
        None
    }
}
