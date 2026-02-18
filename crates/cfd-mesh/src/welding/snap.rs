//! Snap-to-grid utilities.

use crate::core::scalar::{Real, Point3r};

/// Configuration for snapping.
#[derive(Clone, Debug)]
pub struct SnapConfig {
    /// Grid spacing for snap-to-grid.
    pub grid_spacing: Real,
}

impl SnapConfig {
    /// Snap a coordinate to the nearest grid point.
    #[inline]
    pub fn snap_value(&self, v: Real) -> Real {
        (v / self.grid_spacing).round() * self.grid_spacing
    }

    /// Snap a point to the nearest grid point.
    pub fn snap_point(&self, p: &Point3r) -> Point3r {
        Point3r::new(
            self.snap_value(p.x),
            self.snap_value(p.y),
            self.snap_value(p.z),
        )
    }

    /// Snap all positions in a slice to the grid.
    pub fn snap_all(&self, positions: &mut [Point3r]) {
        for p in positions.iter_mut() {
            *p = self.snap_point(p);
        }
    }
}

impl Default for SnapConfig {
    fn default() -> Self {
        Self {
            grid_spacing: 1e-6, // 1 μm precision
        }
    }
}
