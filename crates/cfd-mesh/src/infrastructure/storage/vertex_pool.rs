//! Vertex pool with spatial-hash deduplication — generic over scalar precision.
//!
//! ## Architecture
//!
//! ```text
//! insert_or_weld(point p)
//!         │
//!  ┌──────▼─────────────────────────────────────────┐
//!  │  Snap-rounding mode (tolerance_sq = None):       │
//!  │  key = floor(p / cell_size)                      │
//!  │  if spatial_hash[key] exists → return first id   │
//!  │  else → insert new vertex                        │
//!  └──────┬─────────────────────────────────────────┘
//!         │ OR
//!  ┌──────▼─────────────────────────────────────────┐
//!  │  Tolerance-based mode (tolerance_sq = Some(ε²)): │
//!  │  search 3×3×3 = 27 neighbouring cells            │
//!  │  return closest existing vertex within ε         │
//!  │  else → insert new vertex                        │
//!  └────────────────────────────────────────────────┘
//! ```
//!
//! ## Theorem — 27-Neighbour Search Completeness
//!
//! For `tolerance ≤ cell_size`, any point p' within distance `tolerance` of p
//! must lie in one of the 27 cells surrounding p's cell.  Proof: the maximum
//! distance from p to any corner of its cell is `cell_size · √3`.  A sphere of
//! radius `tolerance ≤ cell_size` around p is contained within the 27-cell
//! neighbourhood since each adjacent cell is at most one `cell_size` away
//! (Chebyshev distance ≤ 1).  Therefore the 27-cell scan never misses a valid
//! weld candidate.  ∎
//!
//! Unlike O(n²) linear-scan vertex matching, a spatial hash grid provides
//! O(1) amortised dedup lookups. Critical for vertex welding during CSG
//! operations and cross-region stitching.
//!
//! The pool is generic over `T: Scalar` so that both `VertexPool<f64>` (the
//! default) and `VertexPool<f32>` (GPU staging) compile to zero-overhead
//! monomorphised code.
//!
//! ## Exactness Mandate
//!
//! **Snap-rounding mode** uses pure integer grid-cell comparison — no
//! floating-point distance tolerances or epsilon spheres.  Points are
//! deterministically quantized to the nearest grid intersection via
//! `floor(p · inv_cell_size)`.  Two points in the same cell are always welded;
//! two points in different cells are never welded (regardless of geometric
//! proximity).  This gives reproducible, order-independent results.

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
///
/// # Welding Modes
///
/// - **Snap-rounding mode** (`tolerance_sq = None`): pure grid cell hash lookup.
///   Vertices in the same `cell_size`-grid cell are welded unconditionally.  Used
///   for primitive building where grid placement is exact and no distance drift
///   occurs between adjacent vertices.
///
/// - **Tolerance-based mode** (`tolerance_sq = Some(tol_sq)`): checks all 27
///   neighbouring grid cells for any existing vertex within distance `√tol_sq`.
///   Used for CSG operations where Steiner points computed by adjacent faces via
///   different numerical paths may differ by a small floating-point epsilon.
#[derive(Clone)]
pub struct VertexPool<T: Scalar = f64> {
    /// Contiguous vertex storage.
    vertices: Vec<VertexData<T>>,
    /// Spatial hash: grid cell → list of vertex indices in that cell.
    spatial_hash: HashMap<CellKey, Vec<u32>>,
    /// Inverse of the quantization cell size (1 / `cell_size`).
    inv_cell_size: T,
    /// When `Some(tol_sq)`, `insert_or_weld` performs a distance check in the
    /// 27-neighbour cell neighbourhood instead of a single-cell hash lookup.
    tolerance_sq: Option<T>,
}

impl<T: Scalar> VertexPool<T> {
    /// Create a new vertex pool in pure snap-rounding mode.
    pub fn new(cell_size: T) -> Self {
        Self {
            vertices: Vec::new(),
            spatial_hash: HashMap::new(),
            inv_cell_size: T::one() / cell_size,
            tolerance_sq: None,
        }
    }

    /// Create a new vertex pool in tolerance-based welding mode.
    ///
    /// Points within `tolerance` mm of an existing vertex are welded to it.
    /// Internally uses a spatial hash grid with `cell_size = tolerance` and
    /// checks all 27 neighbouring cells for a distance match.
    pub fn with_tolerance(cell_size: T, tolerance: T) -> Self {
        Self {
            vertices: Vec::new(),
            spatial_hash: HashMap::new(),
            inv_cell_size: T::one() / cell_size,
            tolerance_sq: Some(tolerance * tolerance),
        }
    }

    /// Sensible defaults for millifluidic meshes (pure snap-rounding).
    ///
    /// - `cell_size = 1e-4 mm` (100 nm) — points in the same cell are welded.
    #[must_use] 
    pub fn default_millifluidic() -> Self {
        Self::new(<T as Scalar>::from_f64(1e-4))
    }

    /// Pool for CSG boolean operations (tolerance-based welding).
    ///
    /// Uses a `1e-4 mm` grid cell with `1e-4 mm` distance-squared tolerance.
    /// This welds Steiner points whose floating-point positions differ by up to
    /// `1e-4 mm` (100 nm) due to different numerical computation paths across
    /// adjacent faces, preventing T-junction seam gaps in the CSG result.
    #[must_use]
    pub fn for_csg() -> Self {
        Self::with_tolerance(
            <T as Scalar>::from_f64(1e-4),
            <T as Scalar>::from_f64(1e-4),
        )
    }

    /// Create an empty clone with the exact same cell_size and tolerance_sq
    /// settings as this pool. Useful for reconstructing a mesh without regressing
    /// to default scalar tolerances.
    #[must_use]
    pub fn empty_clone(&self) -> Self {
        Self {
            vertices: Vec::new(),
            spatial_hash: HashMap::new(),
            inv_cell_size: self.inv_cell_size,
            tolerance_sq: self.tolerance_sq,
        }
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

    /// Insert a vertex, welding if an existing vertex is nearby.
    ///
    /// In **snap-rounding mode** (`tolerance_sq = None`): welds unconditionally
    /// to the first vertex in the same grid cell.
    ///
    /// In **tolerance-based mode** (`tolerance_sq = Some(tol_sq)`): searches all
    /// 27 neighbouring cells and returns the closest existing vertex whose
    /// distance-squared is ≤ `tol_sq`. If none found, inserts a new vertex.
    pub fn insert_or_weld(&mut self, position: Point3<T>, normal: Vector3<T>) -> VertexId {
        let key = CellKey::from_point(&position, self.inv_cell_size);

        if let Some(tol_sq) = self.tolerance_sq {
            // Tolerance-based mode: search 27-neighbour cells.
            let mut best_id: Option<u32> = None;
            let mut best_dist_sq = tol_sq;
            for dz in -1i64..=1 {
                for dy in -1i64..=1 {
                    for dx in -1i64..=1 {
                        let nk = CellKey { x: key.x + dx, y: key.y + dy, z: key.z + dz };
                        if let Some(indices) = self.spatial_hash.get(&nk) {
                            for &idx in indices {
                                let v = &self.vertices[idx as usize];
                                let d = (v.position - position).norm_squared();
                                if d <= best_dist_sq {
                                    best_dist_sq = d;
                                    best_id = Some(idx);
                                }
                            }
                        }
                    }
                }
            }
            if let Some(idx) = best_id {
                return VertexId::new(idx);
            }
        } else {
            // Snap-rounding mode: single cell lookup.
            if let Some(indices) = self.spatial_hash.get(&key) {
                return VertexId::new(indices[0]);
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
