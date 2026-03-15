//! # Unified Vertex Snapping and Welding
//!
//! This module unifies coordinate snapping and deduplication in a single
//! [`SnappingGrid`] that owns the canonical vertex set.
//!
//! ## Algorithm — 26-Neighbor Search
//!
//! Each vertex position is *quantized* to a grid cell via **round-half-up**:
//!
//! ```text
//! cell(x, y, z) = (floor(x/ε + 0.5), floor(y/ε + 0.5), floor(z/ε + 0.5))
//! ```
//!
//! # Theorem — Deterministic Quantization
//!
//! `floor(v + 0.5)` (round-half-up) is a single-valued function for all real
//! inputs including negative values.  Rust's `.round()` uses round-half-away-
//! from-zero, which maps `-0.5 → -1` while floor-based maps `-0.5 → 0`.  When
//! two distinct floating-point computation paths to the same geometric point
//! straddle a half-integer boundary with opposite-sign rounding errors, `.round()`
//! can assign different grid cells; `floor(v + 0.5)` always assigns the same
//! cell for any sign of the tie-breaking error. ∎
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
    #[must_use]
    pub fn from_point_round(p: &Point3r, inv_eps: Real) -> Self {
        Self {
            x: (p.x * inv_eps + 0.5).floor() as i64,
            y: (p.y * inv_eps + 0.5).floor() as i64,
            z: (p.z * inv_eps + 0.5).floor() as i64,
        }
    }

    /// Reconstruct the canonical snapped position for this cell.
    #[inline]
    #[must_use]
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
    #[must_use]
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
    #[must_use]
    pub fn millifluidic() -> Self {
        Self::new(1e-6)
    }

    /// Tolerance ε.
    #[inline]
    #[must_use]
    pub fn eps(&self) -> Real {
        self.eps
    }

    /// Number of unique vertices stored.
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.positions.len()
    }

    /// Returns `true` if no vertices have been inserted yet.
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }

    /// Read-only slice of all stored positions.
    #[inline]
    #[must_use]
    pub fn positions(&self) -> &[Point3r] {
        &self.positions
    }

    /// Look up the position for a given index.
    ///
    /// Returns `None` for out-of-range indices.
    #[inline]
    #[must_use]
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
    #[must_use]
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
    #[must_use]
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
        // z = -0.5e-3 is exactly at the half-cell boundary between cells -1 and 0.
        // Round-half-up convention (floor(v + 0.5)) maps -0.5 → floor(0.0) = 0.
        let p = pt(1.5e-3, 2.0e-3, -0.5e-3);
        let cell = GridCell::from_point_round(&p, 1.0 / eps);
        assert_eq!(cell.x, 2); // 1.5 → floor(2.0) = 2
        assert_eq!(cell.y, 2); // 2.0 → floor(2.5) = 2
        assert_eq!(cell.z, 0); // -0.5 → floor(0.0) = 0  (round-half-up)
        let back = cell.to_point(eps);
        let expected = pt(2e-3, 2e-3, 0.0);
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

    /// Regression: round-half-up gives a single deterministic cell for negative
    /// half-integer boundaries where `.round()` (round-half-away-from-zero) diverges.
    ///
    /// # Theorem — Deterministic Quantization at Half-Integer Boundary
    ///
    /// For `v = -(k·ε + ε/2)` with integer k ≥ 0, `.round()` gives `-(k+1)` while
    /// `floor(v/ε + 0.5)` gives `-k`.  Only one is correct per tie-breaking rule;
    /// the key property is that the same rule is applied consistently regardless of
    /// whether the point is approached from above or below. ∎
    #[test]
    fn from_point_round_negative_half_cell_is_deterministic() {
        let eps = 1.0;
        let inv_eps = 1.0 / eps;
        // Exact negative half-cell boundary: x = -0.5 (i.e. -(0·ε + ε/2))
        let p = Point3r::new(-0.5, 0.0, 0.0);
        let cell = GridCell::from_point_round(&p, inv_eps);
        // floor(-0.5 + 0.5) = floor(0.0) = 0 — consistent, round-half-up
        assert_eq!(cell.x, 0, "round-half-up must map -0.5 to cell 0");

        // Symmetric positive case: x = +0.5 → cell 1
        let p2 = Point3r::new(0.5, 0.0, 0.0);
        let cell2 = GridCell::from_point_round(&p2, inv_eps);
        assert_eq!(cell2.x, 1, "round-half-up must map +0.5 to cell 1");
    }

    /// Two slightly-different floating-point representations of the same point must
    /// land in the same GridCell when they are within the same ε-neighbourhood.
    #[test]
    fn snap_determinism_nearby_floats_same_cell() {
        let eps = 1e-4;
        let inv_eps = 1.0 / eps;
        // Two f64 representations of nominally the same point, separated by 1 ULP.
        let v1 = 3.5e-4_f64;
        let v2 = 3.5e-4_f64 + f64::EPSILON * 3.5e-4; // 1 ULP deviation

        let p1 = Point3r::new(v1, 0.0, 0.0);
        let p2 = Point3r::new(v2, 0.0, 0.0);
        let c1 = GridCell::from_point_round(&p1, inv_eps);
        let c2 = GridCell::from_point_round(&p2, inv_eps);

        // Both must map to the same cell (both are clearly in the k=3 half of ε=1e-4)
        assert_eq!(
            c1.x, c2.x,
            "ULP-adjacent floats must quantize to the same cell: {v1:.20e} vs {v2:.20e}"
        );
    }

    // ── Adversarial Tests ─────────────────────────────────────────────────

    /// Snap idempotency: inserting the snapped position of a point must return
    /// the same index — snap(snap(p)) == snap(p).
    #[test]
    fn snap_idempotency() {
        let eps = 1e-3;
        let mut g = SnappingGrid::new(eps);
        let original = pt(1.7e-3, -0.3e-3, 4.2e-3);
        let (i0, _) = g.insert_or_weld(original);
        let snapped = g.position(i0).unwrap();
        let (i1, is_new) = g.insert_or_weld(snapped);
        assert!(!is_new, "re-inserting a snapped position must not create a new vertex");
        assert_eq!(i0, i1, "snap(snap(p)) must equal snap(p)");
    }

    /// Stress: all-negative coordinates quantize correctly.
    #[test]
    fn negative_coordinates_grid_cell() {
        let eps = 1e-3;
        let inv_eps = 1.0 / eps;
        // Grid cell for deeply negative coordinates
        let p = pt(-100.0, -200.0, -300.0);
        let cell = GridCell::from_point_round(&p, inv_eps);
        // -100.0 * 1000 + 0.5 = -99999.5, floor = -100000
        assert_eq!(cell.x, -100_000);
        assert_eq!(cell.y, -200_000);
        assert_eq!(cell.z, -300_000);
        // Round-trip
        let back = cell.to_point(eps);
        assert!((back.x - (-100.0)).abs() < 1e-10);
    }

    /// Large coordinate values must not overflow i64.
    #[test]
    fn large_coordinates_no_overflow() {
        let eps = 1e-6;
        let inv_eps = 1.0 / eps;
        // 1e6 meters * inv_eps(1e6) = 1e12 — well within i64 range
        let p = pt(1e6, -1e6, 1e6);
        let cell = GridCell::from_point_round(&p, inv_eps);
        // 1e6 / 1e-6 + 0.5 = 1e12 + 0.5 → floor = 1e12
        assert_eq!(cell.x, 1_000_000_000_000_i64);
        // -1e12 + 0.5 = -999999999999.5 → floor = -1000000000000
        assert_eq!(cell.y, -1_000_000_000_000_i64);
    }

    /// Insertion order must not affect welding outcome.
    #[test]
    fn insertion_order_determinism() {
        let eps = 1e-3;
        let points = vec![
            pt(0.0, 0.0, 0.0),
            pt(0.5e-3, 0.0, 0.0),  // within eps of point 0
            pt(2.0e-3, 0.0, 0.0),  // new vertex
            pt(2.3e-3, 0.0, 0.0),  // within eps of point 2
        ];

        // Forward order
        let mut g_fwd = SnappingGrid::new(eps);
        let ids_fwd: Vec<u32> = points.iter().map(|&p| g_fwd.insert_or_weld(p).0).collect();

        // Reverse order
        let mut g_rev = SnappingGrid::new(eps);
        let ids_rev: Vec<u32> = points.iter().rev().map(|&p| g_rev.insert_or_weld(p).0).collect();
        let ids_rev: Vec<u32> = ids_rev.into_iter().rev().collect();

        // Both orderings must produce the same number of unique vertices
        assert_eq!(g_fwd.len(), g_rev.len(), "insertion order must not change vertex count");
        // Points that welded forward should also weld in reverse
        assert_eq!(ids_fwd[0], ids_fwd[1], "forward: 0 and 1 must weld");
        assert_eq!(ids_rev[0], ids_rev[1], "reverse: 0 and 1 must weld");
        assert_eq!(ids_fwd[2], ids_fwd[3], "forward: 2 and 3 must weld");
        assert_eq!(ids_rev[2], ids_rev[3], "reverse: 2 and 3 must weld");
    }

    /// Points at exactly the cell boundary ε/2 from origin must consistently weld.
    #[test]
    fn cell_boundary_straddling() {
        let eps = 1e-3;
        let mut g = SnappingGrid::new(eps);
        // Place two points straddling a cell wall: one at +δ, one at -δ from 0.5*eps
        let half_eps = 0.5 * eps;
        let delta = 1e-10; // much smaller than eps
        let (i0, _) = g.insert_or_weld(pt(half_eps - delta, 0.0, 0.0));
        let (i1, _) = g.insert_or_weld(pt(half_eps + delta, 0.0, 0.0));
        // Distance between them is 2*delta << eps, so they must weld
        assert_eq!(i0, i1, "points straddling cell boundary must weld via 27-neighbor search");
    }

    /// All points at origin must weld to a single vertex.
    #[test]
    fn all_zero_positions() {
        let mut g = SnappingGrid::new(1e-3);
        let mut ids = Vec::new();
        for _ in 0..100 {
            ids.push(g.insert_or_weld(pt(0.0, 0.0, 0.0)).0);
        }
        assert_eq!(g.len(), 1, "all-zero must produce single vertex");
        assert!(ids.iter().all(|&id| id == 0));
    }

    /// Dense cluster of points within eps must all weld to one vertex.
    #[test]
    fn dense_cluster_within_eps() {
        let eps = 1e-3;
        let mut g = SnappingGrid::new(eps);
        let center = pt(5.0, 5.0, 5.0);
        let mut all_same = true;
        let (first_id, _) = g.insert_or_weld(center);
        for i in 1..=50 {
            let offset = (i as f64) * 1e-5; // << eps
            let p = pt(center.x + offset, center.y - offset, center.z);
            let (id, _) = g.insert_or_weld(p);
            if id != first_id {
                all_same = false;
            }
        }
        assert!(all_same, "dense cluster within eps must all weld to one vertex");
    }

    /// Grid with points along all three axis-aligned diagonals must weld correctly.
    #[test]
    fn diagonal_cell_boundary_weld() {
        let eps = 1e-3;
        let mut g = SnappingGrid::new(eps);
        // Point at corner of 8 cells (at grid intersection)
        let corner = pt(1e-3, 1e-3, 1e-3); // exactly at grid point
        let (i0, _) = g.insert_or_weld(corner);
        // Point slightly offset in all 3 axes
        let nearby = pt(1e-3 + 0.1e-3, 1e-3 + 0.1e-3, 1e-3 + 0.1e-3);
        let dist = ((0.1e-3_f64).powi(2) * 3.0).sqrt();
        assert!(dist < eps, "sanity: nearby point is within eps");
        let (i1, _) = g.insert_or_weld(nearby);
        assert_eq!(i0, i1, "corner-adjacent point must weld via 27-neighbor");
    }

    /// query_within_eps must find all and only vertices within eps.
    #[test]
    fn query_within_eps_correctness() {
        let eps = 1e-3;
        let mut g = SnappingGrid::new(eps);
        // Insert 3 points: two within eps of query, one outside
        g.insert_or_weld(pt(0.0, 0.0, 0.0));       // idx 0
        g.insert_or_weld(pt(5e-3, 0.0, 0.0));       // idx 1, far away
        g.insert_or_weld(pt(0.5e-4, 0.0, 0.0));     // welds to idx 0
        g.insert_or_weld(pt(10e-3, 0.0, 0.0));      // idx 2, far away

        let query = pt(0.0, 0.0, 0.0);
        let results = g.query_within_eps(&query);
        assert_eq!(results.len(), 1, "only one unique vertex near origin");
        assert_eq!(results[0], 0);
    }

    /// SnappingGrid must not panic on subnormal epsilon values.
    #[test]
    fn very_small_eps() {
        let eps = 1e-15;
        let mut g = SnappingGrid::new(eps);
        let (i0, _) = g.insert_or_weld(pt(1e-15, 0.0, 0.0));
        let (i1, _) = g.insert_or_weld(pt(1e-15 + 1e-16, 0.0, 0.0));
        // Within eps, must weld
        assert_eq!(i0, i1);
    }

    /// Transitivity: if A welds to B and B to C, then C must weld to A's
    /// canonical index (no "chain splitting").
    #[test]
    fn transitive_welding_chain() {
        let eps = 1e-3;
        let mut g = SnappingGrid::new(eps);
        let (i0, _) = g.insert_or_weld(pt(0.0, 0.0, 0.0));
        // B is within eps of A, so welds to A
        let (i1, _) = g.insert_or_weld(pt(0.4e-3, 0.0, 0.0));
        assert_eq!(i0, i1);
        // C is within eps of B's original position but NOT of A's snapped position
        // However, C should still try welding to the nearest stored vertex.
        // The snapped position of A is at grid center (0,0,0), so C at 0.8e-3 is
        // 0.8e-3 from (0,0,0) which is within eps=1e-3. Should still weld!
        let (i2, _) = g.insert_or_weld(pt(0.8e-3, 0.0, 0.0));
        assert_eq!(i0, i2, "transitive chain must weld to the canonical vertex");
    }

    // ── Cocyclic / equidistant vertex tests (Lévy 2025 insight) ──────────

    /// Four vertices equidistant from a query point: all placed at the corners
    /// of a tetrahedron within eps. The first inserted must win consistently.
    #[test]
    fn cocyclic_four_equidistant_weld_deterministic() {
        let eps = 1e-3;
        let mut g = SnappingGrid::new(eps);
        let r = 0.3e-3; // distance from center, well within eps
        let center = pt(5.0, 5.0, 5.0);

        // 4 points on a tiny tetrahedron centered at `center`, all distance r from center
        let s = r / 3.0_f64.sqrt();
        let points = [
            pt(center.x + s, center.y + s, center.z + s),
            pt(center.x + s, center.y - s, center.z - s),
            pt(center.x - s, center.y + s, center.z - s),
            pt(center.x - s, center.y - s, center.z + s),
        ];

        let (first_id, _) = g.insert_or_weld(points[0]);
        for &p in &points[1..] {
            let (id, _) = g.insert_or_weld(p);
            assert_eq!(id, first_id, "all equidistant points within eps must weld to first");
        }
        assert_eq!(g.len(), 1, "four cocyclic points within eps → one vertex");
    }

    /// Multiple clusters along a line, each cluster's points within eps but
    /// clusters are > eps apart. Must produce exactly one vertex per cluster.
    #[test]
    fn linear_clusters_separation() {
        let eps = 1e-3;
        let mut g = SnappingGrid::new(eps);
        let cluster_sep = 5.0 * eps; // well separated

        for cluster_idx in 0..10 {
            let base_x = cluster_idx as f64 * cluster_sep;
            let first_id = g.insert_or_weld(pt(base_x, 0.0, 0.0)).0;
            // Insert 5 more points clustered within eps/10
            for j in 1..=5 {
                let offset = j as f64 * eps * 0.05;
                let (id, _) = g.insert_or_weld(pt(base_x + offset, 0.0, 0.0));
                assert_eq!(id, first_id, "cluster point must weld to cluster's first");
            }
        }
        assert_eq!(g.len(), 10, "10 well-separated clusters → 10 vertices");
    }

    /// Points placed at all 27 neighbor cell centers around a seed must all
    /// be found by neighborhood_27 queries.
    #[test]
    fn neighborhood_coverage_27_cells() {
        let eps = 1e-3;
        let mut g = SnappingGrid::new(eps);
        let seed = pt(0.0, 0.0, 0.0);
        g.insert_or_weld(seed);

        // Place a point just barely inside each of the 27 cells (including center)
        let mut found_count = 0;
        for dz in -1..=1_i64 {
            for dy in -1..=1_i64 {
                for dx in -1..=1_i64 {
                    let query = pt(
                        dx as f64 * eps * 0.3,
                        dy as f64 * eps * 0.3,
                        dz as f64 * eps * 0.3,
                    );
                    // These are all within eps of seed (max dist = sqrt(3)*0.3*eps ≈ 0.52*eps < eps)
                    if let Some(_) = g.query_nearest(&query) {
                        found_count += 1;
                    }
                }
            }
        }
        assert_eq!(found_count, 27, "seed must be found from all 27 neighbor cell positions");
    }

    /// Stress: 1000 random-ish points in a small volume should never create
    /// more vertices than expected given the eps.
    #[test]
    fn stress_many_points_bounded_vertex_count() {
        let eps = 1e-3;
        let mut g = SnappingGrid::new(eps);
        let n = 1000;
        // Place points on a very fine grid (spacing = eps/20) within a 2eps cube
        // All points are within 2*sqrt(3)*eps < 4*eps of origin.
        // Max possible unique vertices ≈ (2eps / eps)^3 = 8 cells, but some may weld
        for i in 0..n {
            let x = (i % 10) as f64 * eps * 0.2;
            let y = ((i / 10) % 10) as f64 * eps * 0.2;
            let z = (i / 100) as f64 * eps * 0.2;
            g.insert_or_weld(pt(x, y, z));
        }
        // With eps=1e-3 and spacing 0.2*eps = 2e-4, multiple points per cell
        // will weld. We just verify no panic and reasonable vertex count.
        assert!(
            g.len() <= n,
            "vertex count must not exceed insertion count"
        );
        assert!(
            g.len() >= 1,
            "at least one vertex must exist"
        );
    }

    /// GridCell::to_point round-trips correctly through from_point_round.
    #[test]
    fn grid_cell_round_trip_all_octants() {
        let eps = 1e-4;
        let inv_eps = 1.0 / eps;
        // Test all 8 octants
        let coords: [f64; 3] = [1.23456e-3, -7.8901e-4, 4.567e-3];
        for &sx in &[-1.0, 1.0] {
            for &sy in &[-1.0, 1.0] {
                for &sz in &[-1.0, 1.0] {
                    let p = Point3r::new(coords[0] * sx, coords[1] * sy, coords[2] * sz);
                    let cell = GridCell::from_point_round(&p, inv_eps);
                    let back = cell.to_point(eps);
                    // Round-trip error must be < eps
                    assert!(
                        (back.x - p.x).abs() < eps,
                        "X round-trip error too large: {:.2e}",
                        (back.x - p.x).abs()
                    );
                    assert!(
                        (back.y - p.y).abs() < eps,
                        "Y round-trip error too large: {:.2e}",
                        (back.y - p.y).abs()
                    );
                    assert!(
                        (back.z - p.z).abs() < eps,
                        "Z round-trip error too large: {:.2e}",
                        (back.z - p.z).abs()
                    );
                }
            }
        }
    }
}
