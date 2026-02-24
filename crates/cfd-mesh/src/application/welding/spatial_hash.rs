//! Standalone spatial hash grid for point queries.
//!
//! This is the same algorithm used inside `VertexPool`, but exposed as a
//! generic utility for any point set.

use hashbrown::HashMap;

use crate::domain::core::scalar::{Point3r, Real};

/// Key for a spatial hash grid cell.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct GridCell {
    /// X grid coordinate.
    pub x: i64,
    /// Y grid coordinate.
    pub y: i64,
    /// Z grid coordinate.
    pub z: i64,
}

impl GridCell {
    /// Quantize a point to a grid cell.
    #[inline]
    pub fn from_point(p: &Point3r, inv_cell_size: Real) -> Self {
        Self {
            x: (p.x * inv_cell_size).floor() as i64,
            y: (p.y * inv_cell_size).floor() as i64,
            z: (p.z * inv_cell_size).floor() as i64,
        }
    }

    /// Iterate over the 3×3×3 neighborhood (27 cells including self).
    pub fn neighborhood(self) -> impl Iterator<Item = GridCell> {
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

/// A spatial hash grid mapping 3D points to bucket indices.
pub struct SpatialHashGrid {
    /// Buckets: grid cell → list of indices.
    buckets: HashMap<GridCell, Vec<u32>>,
    /// Reciprocal of cell size.
    inv_cell_size: Real,
}

impl SpatialHashGrid {
    /// Create a new spatial hash grid.
    ///
    /// `cell_size` should be roughly the search radius for typical queries.
    pub fn new(cell_size: Real) -> Self {
        assert!(cell_size > 0.0, "cell_size must be positive");
        Self {
            buckets: HashMap::new(),
            inv_cell_size: 1.0 / cell_size,
        }
    }

    /// Insert a point with a given index.
    pub fn insert(&mut self, point: &Point3r, index: u32) {
        let cell = GridCell::from_point(point, self.inv_cell_size);
        self.buckets
            .entry(cell)
            .or_insert_with(|| Vec::with_capacity(4))
            .push(index);
    }

    /// Query: find all indices within `radius` of `point`.
    ///
    /// Searches the 3×3×3 neighborhood and filters by actual distance.
    pub fn query_radius(&self, point: &Point3r, radius: Real, positions: &[Point3r]) -> Vec<u32> {
        let cell = GridCell::from_point(point, self.inv_cell_size);
        let radius_sq = radius * radius;
        let mut results = Vec::new();

        for neighbor in cell.neighborhood() {
            if let Some(indices) = self.buckets.get(&neighbor) {
                for &idx in indices {
                    let dist_sq = (positions[idx as usize] - point).norm_squared();
                    if dist_sq <= radius_sq {
                        results.push(idx);
                    }
                }
            }
        }

        results
    }

    /// Find the nearest point within `radius`, if any.
    pub fn query_nearest(
        &self,
        point: &Point3r,
        radius: Real,
        positions: &[Point3r],
    ) -> Option<u32> {
        let cell = GridCell::from_point(point, self.inv_cell_size);
        let radius_sq = radius * radius;
        let mut best: Option<(u32, Real)> = None;

        for neighbor in cell.neighborhood() {
            if let Some(indices) = self.buckets.get(&neighbor) {
                for &idx in indices {
                    let dist_sq = (positions[idx as usize] - point).norm_squared();
                    if dist_sq <= radius_sq {
                        if best.map_or(true, |(_, d)| dist_sq < d) {
                            best = Some((idx, dist_sq));
                        }
                    }
                }
            }
        }

        best.map(|(idx, _)| idx)
    }

    /// Clear the grid.
    pub fn clear(&mut self) {
        self.buckets.clear();
    }
}
