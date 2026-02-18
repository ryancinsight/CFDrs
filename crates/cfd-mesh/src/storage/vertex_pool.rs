//! Vertex pool with spatial-hash deduplication.
//!
//! Unlike csgrs's O(n²) linear-scan vertex matching, we use a spatial hash
//! grid for O(1) amortized dedup lookups. This is critical for vertex
//! welding during CSG operations and cross-region stitching.

use hashbrown::HashMap;

use crate::core::index::VertexId;
use crate::core::scalar::{Real, Point3r, Vector3r, TOLERANCE};

/// Data stored per vertex.
#[derive(Clone, Debug)]
pub struct VertexData {
    /// Position in 3D space.
    pub position: Point3r,
    /// Surface normal (may be zero for interior vertices).
    pub normal: Vector3r,
}

impl VertexData {
    /// Create a new vertex with position and normal.
    pub fn new(position: Point3r, normal: Vector3r) -> Self {
        Self { position, normal }
    }

    /// Create a vertex with position only (zero normal).
    pub fn from_position(position: Point3r) -> Self {
        Self {
            position,
            normal: Vector3r::zeros(),
        }
    }

    /// Linear interpolation between two vertices.
    pub fn lerp(&self, other: &Self, t: Real) -> Self {
        let position = Point3r::from(
            self.position.coords * (1.0 - t) + other.position.coords * t,
        );
        let normal = (self.normal * (1.0 - t) + other.normal * t).normalize();
        Self { position, normal }
    }
}

/// Quantized grid cell key for spatial hashing.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
struct CellKey {
    x: i64,
    y: i64,
    z: i64,
}

impl CellKey {
    fn from_point(p: &Point3r, inv_cell_size: Real) -> Self {
        Self {
            x: (p.x * inv_cell_size).floor() as i64,
            y: (p.y * inv_cell_size).floor() as i64,
            z: (p.z * inv_cell_size).floor() as i64,
        }
    }
}

/// A pool of deduplicated vertices backed by a spatial hash grid.
pub struct VertexPool {
    /// Contiguous vertex storage.
    vertices: Vec<VertexData>,
    /// Spatial hash: grid cell → list of vertex indices in that cell.
    spatial_hash: HashMap<CellKey, Vec<u32>>,
    /// 1.0 / cell_size for quantization.
    inv_cell_size: Real,
    /// Welding tolerance (squared, for distance checks).
    tolerance_sq: Real,
}

impl VertexPool {
    /// Create a new vertex pool.
    ///
    /// - `cell_size`: spatial hash grid cell size (should be ≥ welding tolerance).
    /// - `tolerance`: welding tolerance — vertices closer than this are merged.
    pub fn new(cell_size: Real, tolerance: Real) -> Self {
        assert!(cell_size > 0.0, "cell_size must be positive");
        Self {
            vertices: Vec::new(),
            spatial_hash: HashMap::new(),
            inv_cell_size: 1.0 / cell_size,
            tolerance_sq: tolerance * tolerance,
        }
    }

    /// Create with sensible defaults for millifluidic meshes (cell_size = 0.01 mm).
    pub fn default_millifluidic() -> Self {
        Self::new(0.01, TOLERANCE)
    }

    /// Number of unique vertices.
    #[inline]
    pub fn len(&self) -> usize {
        self.vertices.len()
    }

    /// Is the pool empty?
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.vertices.is_empty()
    }

    /// Insert a vertex, deduplicating against existing vertices within tolerance.
    ///
    /// Returns the `VertexId` of the existing or newly-created vertex.
    pub fn insert_or_weld(&mut self, position: Point3r, normal: Vector3r) -> VertexId {
        let key = CellKey::from_point(&position, self.inv_cell_size);

        // Check the 3×3×3 neighborhood of grid cells for matches.
        for dz in -1i64..=1 {
            for dy in -1i64..=1 {
                for dx in -1i64..=1 {
                    let neighbor = CellKey {
                        x: key.x + dx,
                        y: key.y + dy,
                        z: key.z + dz,
                    };
                    if let Some(indices) = self.spatial_hash.get(&neighbor) {
                        for &idx in indices {
                            let existing = &self.vertices[idx as usize];
                            let dist_sq = (existing.position - position).norm_squared();
                            if dist_sq <= self.tolerance_sq {
                                // Weld: average the normals for smoother shading.
                                // We don't move the position — first-come wins.
                                return VertexId::new(idx);
                            }
                        }
                    }
                }
            }
        }

        // No match found — insert new vertex.
        let idx = self.vertices.len() as u32;
        self.vertices.push(VertexData::new(position, normal));
        self.spatial_hash
            .entry(key)
            .or_insert_with(|| Vec::with_capacity(2))
            .push(idx);
        VertexId::new(idx)
    }

    /// Insert a vertex **without** deduplication (forced insert).
    pub fn insert_unique(&mut self, position: Point3r, normal: Vector3r) -> VertexId {
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
    pub fn get(&self, id: VertexId) -> &VertexData {
        &self.vertices[id.as_usize()]
    }

    /// Get vertex data mutably by ID.
    #[inline]
    pub fn get_mut(&mut self, id: VertexId) -> &mut VertexData {
        &mut self.vertices[id.as_usize()]
    }

    /// Get vertex position by ID.
    #[inline]
    pub fn position(&self, id: VertexId) -> &Point3r {
        &self.vertices[id.as_usize()].position
    }

    /// Get vertex normal by ID.
    #[inline]
    pub fn normal(&self, id: VertexId) -> &Vector3r {
        &self.vertices[id.as_usize()].normal
    }

    /// Iterate over all vertices with IDs.
    pub fn iter(&self) -> impl Iterator<Item = (VertexId, &VertexData)> {
        self.vertices
            .iter()
            .enumerate()
            .map(|(i, v)| (VertexId::new(i as u32), v))
    }

    /// Access all positions as a slice (for bulk operations).
    pub fn positions(&self) -> impl Iterator<Item = &Point3r> {
        self.vertices.iter().map(|v| &v.position)
    }

    /// Clear all vertices and the spatial hash.
    pub fn clear(&mut self) {
        self.vertices.clear();
        self.spatial_hash.clear();
    }
}

impl Default for VertexPool {
    fn default() -> Self {
        Self::default_millifluidic()
    }
}
