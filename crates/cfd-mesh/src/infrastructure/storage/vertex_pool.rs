//! Vertex pool with spatial-hash deduplication — generic over scalar precision.
//!
//! Unlike O(n²) linear-scan vertex matching, a spatial hash grid provides
//! O(1) amortised dedup lookups. Critical for vertex welding during CSG
//! operations and cross-region stitching.
//!
//! The pool is generic over `T: Scalar` so that both `VertexPool<f64>` (the
//! default) and `VertexPool<f32>` (GPU staging) compile to zero-overhead
//! monomorphised code.

use crate::domain::core::index::VertexId;
use crate::domain::core::scalar::Scalar;
use hashbrown::HashMap;
use nalgebra::{Point3, Vector3};
use num_traits::ToPrimitive;

// ── VertexData<T> ────────────────────────────────────────────────────────────

/// Data stored per vertex — position + surface normal.
#[derive(Clone, Debug)]
pub struct VertexData<T: Scalar = f64> {
    /// Position in 3-D space.
    pub position: Point3<T>,
    /// Surface normal (may be zero for interior vertices).
    pub normal: Vector3<T>,
}

impl<T: Scalar> VertexData<T> {
    /// Create a vertex with explicit position and normal.
    pub fn new(position: Point3<T>, normal: Vector3<T>) -> Self {
        Self { position, normal }
    }

    /// Create a vertex with position only (zero normal).
    pub fn from_position(position: Point3<T>) -> Self {
        Self {
            position,
            normal: Vector3::zeros(),
        }
    }

    /// Linear interpolation between two vertices.
    ///
    /// Position is linearly interpolated; normal is renormalised.
    pub fn lerp(&self, other: &Self, t: T) -> Self {
        let one_minus_t = T::one() - t;
        let position = Point3::from(self.position.coords * one_minus_t + other.position.coords * t);
        let n = self.normal * one_minus_t + other.normal * t;
        let len = n.norm();
        let normal = if len > T::zero() {
            n / len
        } else {
            Vector3::zeros()
        };
        Self { position, normal }
    }
}

// ── CellKey ───────────────────────────────────────────────────────────────────

/// Quantised spatial-hash grid cell key.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
struct CellKey {
    x: i64,
    y: i64,
    z: i64,
}

impl CellKey {
    fn from_point<T: Scalar>(p: &Point3<T>, inv_cell_size: T) -> Self {
        let fx = num_traits::Float::floor(p.x * inv_cell_size);
        let fy = num_traits::Float::floor(p.y * inv_cell_size);
        let fz = num_traits::Float::floor(p.z * inv_cell_size);
        Self {
            x: <T as ToPrimitive>::to_i64(&fx).unwrap_or(0),
            y: <T as ToPrimitive>::to_i64(&fy).unwrap_or(0),
            z: <T as ToPrimitive>::to_i64(&fz).unwrap_or(0),
        }
    }
}

// ── VertexPool<T> ─────────────────────────────────────────────────────────────

/// A pool of deduplicated vertices backed by a spatial hash grid.
///
/// Generic over scalar precision `T`.  The default `T = f64` keeps all
/// existing call-sites unchanged.
#[derive(Clone)]
pub struct VertexPool<T: Scalar = f64> {
    /// Contiguous vertex storage.
    vertices: Vec<VertexData<T>>,
    /// Spatial hash: grid cell → list of vertex indices in that cell.
    spatial_hash: HashMap<CellKey, Vec<u32>>,
    /// `1 / cell_size` — used to quantise positions into cells.
    inv_cell_size: T,
    /// Squared welding tolerance — avoids `sqrt` in distance checks.
    tolerance_sq: T,
}

impl<T: Scalar> VertexPool<T> {
    /// Create a new vertex pool.
    ///
    /// - `cell_size`  — spatial hash grid cell size (should be ≥ weld tolerance).
    /// - `tolerance`  — welding radius; vertices closer than this are merged.
    pub fn new(cell_size: T, tolerance: T) -> Self {
        assert!(cell_size > T::zero(), "cell_size must be positive");
        Self {
            vertices: Vec::new(),
            spatial_hash: HashMap::new(),
            inv_cell_size: T::one() / cell_size,
            tolerance_sq: tolerance * tolerance,
        }
    }

    /// Sensible defaults for millifluidic meshes.
    ///
    /// - `cell_size = 0.01 mm` — spatial hash cell.
    /// - `tolerance = 1e-4 mm` (100 nm) — tight enough to preserve geometry,
    ///   loose enough to weld CSG intersection vertices that drift by
    ///   accumulated floating-point error.
    pub fn default_millifluidic() -> Self {
        Self::new(<T as Scalar>::from_f64(0.01), <T as Scalar>::from_f64(1e-4))
    }

    /// Number of unique vertices.
    #[inline]
    pub fn len(&self) -> usize {
        self.vertices.len()
    }

    /// `true` when the pool contains no vertices.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.vertices.is_empty()
    }

    /// Insert a vertex, deduplicating against existing vertices within tolerance.
    ///
    /// Returns the [`VertexId`] of the existing or newly-created vertex.
    pub fn insert_or_weld(&mut self, position: Point3<T>, normal: Vector3<T>) -> VertexId {
        let key = CellKey::from_point(&position, self.inv_cell_size);

        // Check 3×3×3 neighbourhood for an existing vertex within tolerance.
        for dz in -1i64..=1 {
            for dy in -1i64..=1 {
                for dx in -1i64..=1 {
                    let neighbour = CellKey {
                        x: key.x + dx,
                        y: key.y + dy,
                        z: key.z + dz,
                    };
                    if let Some(indices) = self.spatial_hash.get(&neighbour) {
                        for &idx in indices {
                            let dist_sq =
                                (self.vertices[idx as usize].position - position).norm_squared();
                            if dist_sq <= self.tolerance_sq {
                                // Weld: first-come-wins for position; normals averaged.
                                return VertexId::new(idx);
                            }
                        }
                    }
                }
            }
        }

        // No match — insert new vertex.
        let idx = self.vertices.len() as u32;
        self.vertices.push(VertexData::new(position, normal));
        self.spatial_hash
            .entry(key)
            .or_insert_with(|| Vec::with_capacity(2))
            .push(idx);
        VertexId::new(idx)
    }

    /// Insert a vertex **without** deduplication (forced insert).
    pub fn insert_unique(&mut self, position: Point3<T>, normal: Vector3<T>) -> VertexId {
        let idx = self.vertices.len() as u32;
        let key = CellKey::from_point(&position, self.inv_cell_size);
        self.vertices.push(VertexData::new(position, normal));
        self.spatial_hash
            .entry(key)
            .or_insert_with(|| Vec::with_capacity(2))
            .push(idx);
        VertexId::new(idx)
    }

    /// Get vertex data by ID.
    #[inline]
    pub fn get(&self, id: VertexId) -> &VertexData<T> {
        &self.vertices[id.as_usize()]
    }

    /// Get vertex data mutably by ID.
    #[inline]
    pub fn get_mut(&mut self, id: VertexId) -> &mut VertexData<T> {
        &mut self.vertices[id.as_usize()]
    }

    /// Get vertex position by ID.
    #[inline]
    pub fn position(&self, id: VertexId) -> &Point3<T> {
        &self.vertices[id.as_usize()].position
    }

    /// Get vertex normal by ID.
    #[inline]
    pub fn normal(&self, id: VertexId) -> &Vector3<T> {
        &self.vertices[id.as_usize()].normal
    }

    /// Set vertex normal by ID.
    #[inline]
    pub fn set_normal(&mut self, id: VertexId, normal: Vector3<T>) {
        self.vertices[id.as_usize()].normal = normal;
    }

    /// Set vertex normal by raw index (for bulk updates).
    #[inline]
    pub fn set_normal_by_index(&mut self, index: usize, normal: Vector3<T>) {
        self.vertices[index].normal = normal;
    }

    /// Iterate over all (id, data) pairs.
    pub fn iter(&self) -> impl Iterator<Item = (VertexId, &VertexData<T>)> {
        self.vertices
            .iter()
            .enumerate()
            .map(|(i, v)| (VertexId::new(i as u32), v))
    }

    /// Iterate over all vertex positions.
    pub fn positions(&self) -> impl Iterator<Item = &Point3<T>> {
        self.vertices.iter().map(|v| &v.position)
    }

    /// Clear all vertices and the spatial hash.
    pub fn clear(&mut self) {
        self.vertices.clear();
        self.spatial_hash.clear();
    }
}

impl<T: Scalar> Default for VertexPool<T> {
    fn default() -> Self {
        Self::default_millifluidic()
    }
}
