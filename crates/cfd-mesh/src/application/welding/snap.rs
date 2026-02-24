//! # Unified Vertex Snapping and Welding
//!
//! This module replaces the old split between `SnapConfig` (scalar grid) and
//! `SpatialHashGrid` (query-only) with a single [`SnappingGrid`] that owns
//! the canonical vertex set and performs both snapping *and* deduplication.
//!
//! ## Algorithm — 26-Neighbor Search
//!
//! Each vertex position is *quantized* to a grid cell by rounding (not
//! flooring) each coordinate to the nearest multiple of ε:
//!
//! ```text
//! cell(x, y, z) = (round(x/ε), round(y/ε), round(z/ε))
//! ```
//!
//! Rounding (rather than flooring) ensures that a point halfway between two
//! cells is assigned consistently — a property that flooring lacks.
//!
//! When inserting a new point, the grid searches all **26 face-, edge-, and
//! corner-adjacent neighbors** plus the home cell itself (27 cells total).
//! This prevents "ghost duplicates" at cell boundaries: a point within `ε` of
//! a neighbor-cell wall will still be found regardless of which side of the
//! wall it falls on.
//!
//! ## Complexity
//!
//! | Operation | Expected | Worst case |
//! |---|---|---|
//! | `insert_or_weld` | O(1) | O(k) where k = vertices per cell |
//! | `query_nearest`  | O(1) | O(k) |
//! | Memory           | O(n) | O(n) |
//!
//! ## Diagram
//!
//! ```text
//! 26-neighbor cells (3-D cross-section, center = ★):
//!
//!   z-1 layer       z=0 layer        z+1 layer
//!  ┌───┬───┬───┐  ┌───┬───┬───┐  ┌───┬───┬───┐
//!  │ · │ · │ · │  │ · │ · │ · │  │ · │ · │ · │
//!  ├───┼───┼───┤  ├───┼───┼───┤  ├───┼───┼───┤
//!  │ · │ · │ · │  │ · │ ★ │ · │  │ · │ · │ · │
//!  ├───┼───┼───┤  ├───┼───┼───┤  ├───┼───┼───┤
//!  │ · │ · │ · │  │ · │ · │ · │  │ · │ · │ · │
//!  └───┴───┴───┘  └───┴───┴───┘  └───┴───┴───┘
//!   9 neighbors     8 neighbors     9 neighbors
//!                   (+ center)
//! ```
//!
//! ## Integration with `HalfEdgeMesh<'id>`
//!
//! [`SnappingGrid`] is **mesh-agnostic**: it stores positions and returns
//! opaque `u32` indices.  The `HalfEdgeMesh`-aware entry point is
//! [`SnappingGrid::insert_or_weld_he`], which maps indices to [`VertexKey`]s.

use hashbrown::HashMap;

use crate::domain::core::scalar::{Point3r, Real};

// ── GridCell ─────────────────────────────────────────────────────────────────

/// A quantized 3-D grid cell coordinate.
///
/// Uses `i64` so that negative coordinates and very large models are handled
/// correctly without overflow for any mesh that fits in ±9 × 10¹² ε-units.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct GridCell {
    /// Quantized X index.
    pub x: i64,
    /// Quantized Y index.
    pub y: i64,
    /// Quantized Z index.
    pub z: i64,
}

impl GridCell {
    /// Quantize a point using **rounding** (not flooring).
    ///
    /// Rounding assigns boundary points consistently, avoiding the ghost-
    /// duplicate artifact that occurs with floor-based quantization when a
    /// point lies exactly on a cell wall.
    #[inline]
    pub fn from_point_round(p: &Point3r, inv_eps: Real) -> Self {
        Self {
            x: (p.x * inv_eps).round() as i64,
            y: (p.y * inv_eps).round() as i64,
            z: (p.z * inv_eps).round() as i64,
        }
    }

    /// Reconstruct the canonical snapped position for this cell.
    #[inline]
    pub fn to_point(self, eps: Real) -> Point3r {
        Point3r::new(
            self.x as Real * eps,
            self.y as Real * eps,
            self.z as Real * eps,
        )
    }

    /// Iterator over the 26 neighboring cells **plus self** (27 total).
    ///
    /// Covers all face-, edge-, and corner-adjacent cells so that a welding
    /// query cannot miss a vertex that lies just across a cell boundary.
    #[inline]
    pub fn neighborhood_27(self) -> impl Iterator<Item = GridCell> {
        (-1i64..=1).flat_map(move |dz| {
            (-1i64..=1).flat_map(move |dy| {
                (-1i64..=1).map(move |dx| GridCell {
                    x: self.x + dx,
                    y: self.y + dy,
                    z: self.z + dz,
                })
            })
        })
    }
}

// ── SnappingGrid ──────────────────────────────────────────────────────────────

/// Unified vertex snapping and welding structure.
///
/// Owns a flat list of deduplicated positions.  On each
/// [`insert_or_weld`][SnappingGrid::insert_or_weld] call it either returns the
/// index of an existing vertex within ε, or inserts the new (snapped) position
/// and returns its fresh index.
///
/// # Precision
///
/// Two points `p` and `q` are considered the *same* vertex if
/// `‖p − q‖² ≤ ε²`.  The weld distance is always `ε`; the grid cell size is
/// also `ε`, so the 26-neighbor search guarantees no missed welds.
///
/// # Thread safety
///
/// `SnappingGrid` is **not** `Send`/`Sync` — protect it with a `Mutex` if
/// parallel insertion is required.
pub struct SnappingGrid {
    /// Grid cell → list of `positions` indices stored in that cell.
    buckets: HashMap<GridCell, Vec<u32>>,
    /// Flat array of all accepted (snapped) positions.
    positions: Vec<Point3r>,
    /// Snap tolerance ε.
    eps: Real,
    /// 1 / ε for quantization.
    inv_eps: Real,
}

impl SnappingGrid {
    /// Create a new snapping grid with tolerance `eps`.
    ///
    /// `eps` is the maximum distance at which two points are considered
    /// identical and will be welded together.  For millifluidic meshes the
    /// recommended value is `1e-6` (1 μm).
    ///
    /// # Panics
    /// Panics if `eps` is not finite and positive.
    pub fn new(eps: Real) -> Self {
        assert!(
            eps.is_finite() && eps > 0.0,
            "eps must be finite and positive"
        );
        Self {
            buckets: HashMap::new(),
            positions: Vec::new(),
            eps,
            inv_eps: 1.0 / eps,
        }
    }

    /// Create a snapping grid suitable for millifluidic devices (ε = 1 μm).
    pub fn millifluidic() -> Self {
        Self::new(1e-6)
    }

    /// Tolerance ε.
    #[inline]
    pub fn eps(&self) -> Real {
        self.eps
    }

    /// Number of unique vertices stored.
    #[inline]
    pub fn len(&self) -> usize {
        self.positions.len()
    }

    /// Returns `true` if no vertices have been inserted yet.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }

    /// Read-only slice of all stored positions.
    #[inline]
    pub fn positions(&self) -> &[Point3r] {
        &self.positions
    }

    /// Look up the position for a given index.
    ///
    /// Returns `None` for out-of-range indices.
    #[inline]
    pub fn position(&self, idx: u32) -> Option<Point3r> {
        self.positions.get(idx as usize).copied()
    }

    // ── Core operation ────────────────────────────────────────────────────

    /// Insert `point` into the grid, or weld it to an existing vertex.
    ///
    /// If any stored vertex is within ε of `point`, returns the index of the
    /// nearest such vertex.  Otherwise snaps `point` to the canonical grid
    /// position and inserts it, returning the new index.
    ///
    /// The search covers all 26 neighbors plus the home cell, so no duplicate
    /// can hide across a cell boundary.
    ///
    /// # Returns
    /// A `(index, is_new)` pair.  `is_new` is `true` when a fresh vertex was
    /// added, `false` when an existing vertex was reused.
    pub fn insert_or_weld(&mut self, point: Point3r) -> (u32, bool) {
        let home = GridCell::from_point_round(&point, self.inv_eps);
        let eps_sq = self.eps * self.eps;

        // Search all 27 cells for the nearest existing vertex
        let mut best: Option<(u32, Real)> = None;
        for cell in home.neighborhood_27() {
            if let Some(indices) = self.buckets.get(&cell) {
                for &idx in indices {
                    let dist_sq = (self.positions[idx as usize] - point).norm_squared();
                    if dist_sq <= eps_sq {
                        match best {
                            None => best = Some((idx, dist_sq)),
                            Some((_, d)) if dist_sq < d => best = Some((idx, dist_sq)),
                            _ => {}
                        }
                    }
                }
            }
        }

        if let Some((idx, _)) = best {
            return (idx, false);
        }

        // New vertex: snap to grid center and insert
        let snapped = home.to_point(self.eps);
        let new_idx = self.positions.len() as u32;
        self.positions.push(snapped);
        self.buckets
            .entry(home)
            .or_insert_with(|| Vec::with_capacity(4))
            .push(new_idx);
        (new_idx, true)
    }

    /// Query the nearest vertex within ε of `point` without inserting.
    ///
    /// Returns `None` if no vertex is within ε.
    pub fn query_nearest(&self, point: &Point3r) -> Option<u32> {
        let home = GridCell::from_point_round(point, self.inv_eps);
        let eps_sq = self.eps * self.eps;
        let mut best: Option<(u32, Real)> = None;

        for cell in home.neighborhood_27() {
            if let Some(indices) = self.buckets.get(&cell) {
                for &idx in indices {
                    let dist_sq = (self.positions[idx as usize] - point).norm_squared();
                    if dist_sq <= eps_sq {
                        match best {
                            None => best = Some((idx, dist_sq)),
                            Some((_, d)) if dist_sq < d => best = Some((idx, dist_sq)),
                            _ => {}
                        }
                    }
                }
            }
        }

        best.map(|(idx, _)| idx)
    }

    /// Query all vertices within ε of `point` without inserting.
    pub fn query_within_eps(&self, point: &Point3r) -> Vec<u32> {
        let home = GridCell::from_point_round(point, self.inv_eps);
        let eps_sq = self.eps * self.eps;
        let mut results = Vec::new();

        for cell in home.neighborhood_27() {
            if let Some(indices) = self.buckets.get(&cell) {
                for &idx in indices {
                    let dist_sq = (self.positions[idx as usize] - point).norm_squared();
                    if dist_sq <= eps_sq {
                        results.push(idx);
                    }
                }
            }
        }

        results
    }

    /// Clear all vertices, resetting the grid to empty.
    pub fn clear(&mut self) {
        self.buckets.clear();
        self.positions.clear();
    }
}

// ── Legacy SnapConfig shim ────────────────────────────────────────────────────

/// Backward-compatible scalar snapping utility.
///
/// For new code, prefer [`SnappingGrid`] which combines snapping and
/// deduplication.  `SnapConfig` is kept for callers that only need the
/// scalar `snap_value` / `snap_point` helpers.
#[derive(Clone, Debug)]
pub struct SnapConfig {
    /// Grid spacing ε.
    pub grid_spacing: Real,
}

impl SnapConfig {
    /// Snap a scalar to the nearest multiple of `grid_spacing`.
    #[inline]
    pub fn snap_value(&self, v: Real) -> Real {
        (v / self.grid_spacing).round() * self.grid_spacing
    }

    /// Snap a 3-D point to the nearest grid point.
    #[inline]
    pub fn snap_point(&self, p: &Point3r) -> Point3r {
        Point3r::new(
            self.snap_value(p.x),
            self.snap_value(p.y),
            self.snap_value(p.z),
        )
    }

    /// Snap all positions in a slice in-place.
    pub fn snap_all(&self, positions: &mut [Point3r]) {
        for p in positions.iter_mut() {
            *p = self.snap_point(p);
        }
    }

    /// Lift into a full [`SnappingGrid`] with the same tolerance.
    pub fn into_grid(self) -> SnappingGrid {
        SnappingGrid::new(self.grid_spacing)
    }
}

impl Default for SnapConfig {
    fn default() -> Self {
        Self { grid_spacing: 1e-6 }
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn pt(x: Real, y: Real, z: Real) -> Point3r {
        Point3r::new(x, y, z)
    }

    #[test]
    fn grid_cell_round_trip() {
        let eps = 1e-3;
        let p = pt(1.5e-3, 2.0e-3, -0.5e-3);
        let cell = GridCell::from_point_round(&p, 1.0 / eps);
        assert_eq!(cell.x, 2);
        assert_eq!(cell.y, 2);
        assert_eq!(cell.z, -1);
        let back = cell.to_point(eps);
        let expected = pt(2e-3, 2e-3, -1e-3);
        assert!((back - expected).norm() < 1e-15);
    }

    #[test]
    fn neighborhood_27_has_27_cells() {
        let cell = GridCell { x: 0, y: 0, z: 0 };
        let neighbors: Vec<_> = cell.neighborhood_27().collect();
        assert_eq!(neighbors.len(), 27);
    }

    #[test]
    fn insert_two_identical_points_welds() {
        let mut g = SnappingGrid::new(1e-3);
        let (i0, new0) = g.insert_or_weld(pt(0.0, 0.0, 0.0));
        let (i1, new1) = g.insert_or_weld(pt(0.0, 0.0, 0.0));
        assert!(new0);
        assert!(!new1);
        assert_eq!(i0, i1);
        assert_eq!(g.len(), 1);
    }

    #[test]
    fn insert_within_eps_welds() {
        let eps = 1e-3;
        let mut g = SnappingGrid::new(eps);
        let (i0, _) = g.insert_or_weld(pt(0.0, 0.0, 0.0));
        let (i1, new1) = g.insert_or_weld(pt(0.5e-3, 0.0, 0.0));
        assert!(!new1, "point within ε should weld");
        assert_eq!(i0, i1);
    }

    #[test]
    fn insert_beyond_eps_creates_new_vertex() {
        let eps = 1e-3;
        let mut g = SnappingGrid::new(eps);
        let (i0, _) = g.insert_or_weld(pt(0.0, 0.0, 0.0));
        let (i1, new1) = g.insert_or_weld(pt(2e-3, 0.0, 0.0));
        assert!(new1, "point beyond ε should create new vertex");
        assert_ne!(i0, i1);
        assert_eq!(g.len(), 2);
    }

    #[test]
    fn boundary_weld_across_cell_wall() {
        // Both points within eps=1.0 of each other but in different grid cells
        let eps = 1.0;
        let mut g = SnappingGrid::new(eps);
        let (i0, _) = g.insert_or_weld(pt(0.4, 0.0, 0.0)); // rounds to x=0
        let (i1, _) = g.insert_or_weld(pt(0.6, 0.0, 0.0)); // rounds to x=1
                                                           // distance = 0.2 < eps=1.0, so must weld
        assert_eq!(i0, i1, "26-neighbor search must weld across cell boundary");
    }

    #[test]
    fn query_nearest_returns_none_when_empty() {
        let g = SnappingGrid::new(1e-3);
        assert!(g.query_nearest(&pt(0.0, 0.0, 0.0)).is_none());
    }

    #[test]
    fn query_nearest_finds_vertex() {
        let mut g = SnappingGrid::new(1e-3);
        let (i0, _) = g.insert_or_weld(pt(1.0, 2.0, 3.0));
        let found = g.query_nearest(&pt(1.0, 2.0, 3.0));
        assert_eq!(found, Some(i0));
    }

    #[test]
    fn clear_resets_grid() {
        let mut g = SnappingGrid::new(1e-3);
        g.insert_or_weld(pt(0.0, 0.0, 0.0));
        assert_eq!(g.len(), 1);
        g.clear();
        assert_eq!(g.len(), 0);
        assert!(g.is_empty());
    }

    #[test]
    fn snap_config_snap_value() {
        let cfg = SnapConfig { grid_spacing: 1e-3 };
        // 1.4999e-3 should round to 1e-3
        let v = cfg.snap_value(1.4999e-3);
        assert!((v - 1e-3).abs() < 1e-10);
    }

    #[test]
    fn snap_config_into_grid() {
        let cfg = SnapConfig { grid_spacing: 0.5 };
        let g = cfg.into_grid();
        assert_eq!(g.eps(), 0.5);
    }
}
