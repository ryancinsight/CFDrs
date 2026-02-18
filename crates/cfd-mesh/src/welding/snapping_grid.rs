//! Unified snapping grid combining coordinate snapping and spatial deduplication.
//!
//! `SnappingGrid` replaces the split between `SnapConfig` (coordinate-only
//! quantisation) and `SpatialHashGrid` (point-set deduplication).  It supports
//! both usage patterns:
//!
//! - **Snap-to-grid**: round coordinates to the nearest grid point via
//!   [`snap_point`](SnappingGrid::snap_point).
//! - **Insert + nearest-neighbor dedup**: accumulate a set of 3-D points and
//!   find the nearest existing point within a tolerance via
//!   [`insert_or_find`](SnappingGrid::insert_or_find).
//!
//! Both operations share the same grid spacing, ensuring that the snapping
//! radius and the dedup lookup radius are automatically consistent.
//!
//! ## 26-neighbor search
//!
//! A query inspects the center cell plus all 26 face-, edge- and
//! corner-adjacent cells (the 3×3×3 = 27-cell Moore neighborhood).  This
//! guarantees that no point within the tolerance can fall in an un-searched
//! cell, provided the cell size ≥ tolerance.
//!
//! ```text
//!  ┌─────┬─────┬─────┐
//!  │  ·  │  ·  │  ·  │   ← slice z-1
//!  ├─────┼─────┼─────┤
//!  │  ·  │  ·  │  ·  │
//!  └─────┴─────┴─────┘
//!  ┌─────┬─────┬─────┐
//!  │  ·  │  ·  │  ·  │   ← slice z (query cell at centre)
//!  ├─────┼─────┼─────┤
//!  │  ·  │  Q  │  ·  │   Q = query point
//!  └─────┴─────┴─────┘
//!  ┌─────┬─────┬─────┐
//!  │  ·  │  ·  │  ·  │   ← slice z+1
//!  ├─────┼─────┼─────┤   27 cells total searched
//!  │  ·  │  ·  │  ·  │
//!  └─────┴─────┴─────┘
//! ```

use hashbrown::HashMap;
use crate::core::scalar::{Real, Point3r};

// ── Grid cell key ─────────────────────────────────────────────────────────────

/// Integer grid cell coordinates.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
struct CellKey {
    x: i64,
    y: i64,
    z: i64,
}

impl CellKey {
    #[inline]
    fn from_point(p: &Point3r, inv_cell: Real) -> Self {
        Self {
            x: (p.x * inv_cell).floor() as i64,
            y: (p.y * inv_cell).floor() as i64,
            z: (p.z * inv_cell).floor() as i64,
        }
    }

    /// Iterate over the 3×3×3 = 27 Moore-neighborhood cells.
    #[inline]
    fn neighborhood(self) -> [CellKey; 27] {
        let mut cells = [CellKey { x: 0, y: 0, z: 0 }; 27];
        let mut i = 0;
        let mut dz = -1i64;
        while dz <= 1 {
            let mut dy = -1i64;
            while dy <= 1 {
                let mut dx = -1i64;
                while dx <= 1 {
                    cells[i] = CellKey { x: self.x + dx, y: self.y + dy, z: self.z + dz };
                    i += 1;
                    dx += 1;
                }
                dy += 1;
            }
            dz += 1;
        }
        cells
    }
}

// ── SnappingGrid ──────────────────────────────────────────────────────────────

/// Unified spatial-hash grid for coordinate snapping and point deduplication.
///
/// # Configuration
///
/// - `cell_size`: the spatial hash grid spacing.  **Must be ≥ `tolerance`** so
///   that the 27-cell neighborhood search is guaranteed to find all points
///   within the tolerance radius.
/// - `tolerance`: the maximum distance at which two points are considered
///   coincident during deduplication (squared comparison for speed).
///
/// # Typical usage
///
/// ```rust,ignore
/// use cfd_mesh::welding::SnappingGrid;
///
/// let mut grid = SnappingGrid::new(0.01, 1e-4); // 10 µm cell, 100 nm weld
///
/// let a = grid.insert_or_find([0.0, 0.0, 0.0].into(), &mut positions);
/// let b = grid.insert_or_find([0.0, 0.0, 0.0].into(), &mut positions); // deduped
/// assert_eq!(a, b);
/// ```
pub struct SnappingGrid {
    /// Spatial hash: cell → list of point indices into the external buffer.
    buckets:      HashMap<CellKey, Vec<u32>>,
    /// `1 / cell_size` for fast quantisation.
    inv_cell:     Real,
    /// Cell size (used for coordinate snapping).
    cell_size:    Real,
    /// Squared weld tolerance.
    tolerance_sq: Real,
}

impl SnappingGrid {
    /// Create a new grid.
    ///
    /// # Panics
    /// Panics if `cell_size <= 0` or `tolerance > cell_size`.
    pub fn new(cell_size: Real, tolerance: Real) -> Self {
        assert!(cell_size > 0.0, "cell_size must be positive");
        assert!(
            tolerance <= cell_size,
            "tolerance ({:.2e}) must be ≤ cell_size ({:.2e}) for safe 27-cell search",
            tolerance, cell_size,
        );
        Self {
            buckets:      HashMap::new(),
            inv_cell:     1.0 / cell_size,
            cell_size,
            tolerance_sq: tolerance * tolerance,
        }
    }

    /// Default configuration suited for millifluidic mm-scale meshes.
    ///
    /// - cell size = 0.01 mm (10 µm)
    /// - tolerance = 1e-4 mm (100 nm)
    pub fn millifluidic() -> Self {
        Self::new(0.01, 1e-4)
    }

    // ── Snap-to-grid ──────────────────────────────────────────────────────

    /// Snap a single coordinate to the nearest grid point.
    #[inline]
    pub fn snap_value(&self, v: Real) -> Real {
        (v * self.inv_cell).round() * self.cell_size
    }

    /// Snap a point to the nearest spatial-hash grid point.
    #[inline]
    pub fn snap_point(&self, p: &Point3r) -> Point3r {
        Point3r::new(
            self.snap_value(p.x),
            self.snap_value(p.y),
            self.snap_value(p.z),
        )
    }

    /// Snap all points in `positions` to the grid in-place.
    pub fn snap_all(&self, positions: &mut [Point3r]) {
        for p in positions.iter_mut() {
            *p = self.snap_point(p);
        }
    }

    // ── Insert + dedup ────────────────────────────────────────────────────

    /// Insert `point` into the grid at index `idx` (no dedup check).
    pub fn insert(&mut self, point: &Point3r, idx: u32) {
        let cell = CellKey::from_point(point, self.inv_cell);
        self.buckets
            .entry(cell)
            .or_insert_with(|| Vec::with_capacity(4))
            .push(idx);
    }

    /// Find the existing index nearest to `point` within `tolerance`, or
    /// `None` if no point is near enough.
    ///
    /// `positions` is the external slice from which actual coordinates are read.
    pub fn find_nearest(&self, point: &Point3r, positions: &[Point3r]) -> Option<u32> {
        let cell = CellKey::from_point(point, self.inv_cell);
        let mut best: Option<(u32, Real)> = None;

        for &neighbor in &cell.neighborhood() {
            if let Some(indices) = self.buckets.get(&neighbor) {
                for &idx in indices {
                    let d2 = (positions[idx as usize] - point).norm_squared();
                    if d2 <= self.tolerance_sq {
                        if best.map_or(true, |(_, bd)| d2 < bd) {
                            best = Some((idx, d2));
                        }
                    }
                }
            }
        }

        best.map(|(idx, _)| idx)
    }

    /// Try to find an existing point within tolerance; if none found, push
    /// `point` onto `positions` and register it in the grid.
    ///
    /// Returns the index (existing or new) of the representative point.
    pub fn insert_or_find(&mut self, point: Point3r, positions: &mut Vec<Point3r>) -> u32 {
        if let Some(existing) = self.find_nearest(&point, positions) {
            return existing;
        }
        let idx = positions.len() as u32;
        positions.push(point);
        self.insert(&point, idx);
        idx
    }

    // ── Query utilities ───────────────────────────────────────────────────

    /// Return all indices within `radius` of `point`.
    pub fn query_radius(&self, point: &Point3r, radius: Real, positions: &[Point3r]) -> Vec<u32> {
        let cell = CellKey::from_point(point, self.inv_cell);
        let r2 = radius * radius;
        let mut out = Vec::new();
        for &neighbor in &cell.neighborhood() {
            if let Some(indices) = self.buckets.get(&neighbor) {
                for &idx in indices {
                    if (positions[idx as usize] - point).norm_squared() <= r2 {
                        out.push(idx);
                    }
                }
            }
        }
        out
    }

    /// Clear the grid (does not affect any external positions buffer).
    pub fn clear(&mut self) {
        self.buckets.clear();
    }

    /// Number of occupied cells.
    pub fn cell_count(&self) -> usize {
        self.buckets.len()
    }
}

impl Default for SnappingGrid {
    fn default() -> Self {
        Self::millifluidic()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Point3;

    fn pt(x: f64, y: f64, z: f64) -> Point3r { Point3::new(x, y, z) }

    #[test]
    fn snap_value_rounds_to_nearest_cell() {
        let g = SnappingGrid::new(0.1, 1e-4);
        assert!((g.snap_value(0.05) - 0.1).abs() < 1e-10); // rounds up to 0.1
        assert!((g.snap_value(0.04) - 0.0).abs() < 1e-10); // rounds down to 0.0
    }

    #[test]
    fn insert_or_find_deduplicates_near_points() {
        let mut g = SnappingGrid::new(0.01, 1e-4);
        let mut positions: Vec<Point3r> = Vec::new();
        let a = g.insert_or_find(pt(0.0, 0.0, 0.0), &mut positions);
        // second point is within 1e-5 mm — should be welded
        let b = g.insert_or_find(pt(1e-5, 0.0, 0.0), &mut positions);
        assert_eq!(a, b, "nearby points should be deduped");
        assert_eq!(positions.len(), 1);
    }

    #[test]
    fn insert_or_find_keeps_distant_points() {
        let mut g = SnappingGrid::new(0.01, 1e-4);
        let mut positions: Vec<Point3r> = Vec::new();
        let a = g.insert_or_find(pt(0.0, 0.0, 0.0), &mut positions);
        let b = g.insert_or_find(pt(1.0, 0.0, 0.0), &mut positions);
        assert_ne!(a, b);
        assert_eq!(positions.len(), 2);
    }

    #[test]
    fn neighborhood_covers_27_cells() {
        let k = CellKey { x: 0, y: 0, z: 0 };
        let nbrs = k.neighborhood();
        assert_eq!(nbrs.len(), 27);
        // Contains self
        assert!(nbrs.contains(&k));
    }
}
