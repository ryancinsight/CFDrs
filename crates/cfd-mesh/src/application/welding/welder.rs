//! Batch vertex welding for existing meshes.
//!
//! After CSG operations or mesh import, we may need to re-weld vertices that
//! should be coincident but aren't (e.g., from separate polygon soups).

use hashbrown::HashMap;

use crate::application::welding::spatial_hash::SpatialHashGrid;
use crate::domain::core::index::VertexId;
use crate::domain::core::scalar::{Point3r, Real, TOLERANCE};
use crate::infrastructure::storage::face_store::FaceStore;

/// Weld result.
#[derive(Debug)]
pub struct WeldResult {
    /// Number of vertices that were merged.
    pub vertices_merged: usize,
    /// Number of faces updated.
    pub faces_updated: usize,
}

/// Performs vertex welding on a mesh.
pub struct MeshWelder {
    /// Tolerance for vertex welding.
    tolerance: Real,
}

impl MeshWelder {
    /// Create a new welder with the default tolerance.
    pub fn new() -> Self {
        Self {
            tolerance: TOLERANCE.sqrt(),
        }
    }

    /// Create a welder with a custom tolerance.
    pub fn with_tolerance(tolerance: Real) -> Self {
        Self { tolerance }
    }

    /// Weld duplicate vertices in a position list and remap face references.
    ///
    /// Returns the deduplicated positions and the remap table.
    pub fn weld(&self, positions: &[Point3r], face_store: &mut FaceStore) -> WeldResult {
        let mut grid = SpatialHashGrid::new(self.tolerance * 2.0);
        let mut canonical: Vec<Point3r> = Vec::with_capacity(positions.len());
        let mut remap: HashMap<u32, u32> = HashMap::with_capacity(positions.len());

        for (i, pos) in positions.iter().enumerate() {
            if let Some(existing) = grid.query_nearest(pos, self.tolerance, &canonical) {
                remap.insert(i as u32, existing);
            } else {
                let new_idx = canonical.len() as u32;
                grid.insert(pos, new_idx);
                canonical.push(*pos);
                remap.insert(i as u32, new_idx);
            }
        }

        let vertices_merged = positions.len() - canonical.len();

        // Remap face vertex references
        let mut faces_updated = 0usize;
        for (_, face) in face_store.iter_mut_enumerated() {
            let mut changed = false;
            for v in &mut face.vertices {
                let old_raw = v.raw();
                if let Some(&new_raw) = remap.get(&old_raw) {
                    if new_raw != old_raw {
                        *v = VertexId::new(new_raw);
                        changed = true;
                    }
                }
            }
            if changed {
                faces_updated += 1;
            }
        }

        WeldResult {
            vertices_merged,
            faces_updated,
        }
    }
}

impl Default for MeshWelder {
    fn default() -> Self {
        Self::new()
    }
}
