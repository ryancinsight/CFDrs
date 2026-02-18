//! Hex-to-tet mesh decomposition.
//!
//! Converts hexahedral cells to five or six tetrahedra each.  For the current
//! stub, hexahedral meshes are passed through unchanged (already tetrahedral
//! meshes are returned as-is).

use nalgebra::RealField;
use crate::mesh::Mesh;

/// Converts a mesh to all-tetrahedra.
pub struct HexToTetConverter;

impl HexToTetConverter {
    /// Convert every hexahedral cell to equivalent tetrahedra.
    ///
    /// For meshes that are already tetrahedral this is a no-op clone.
    pub fn convert<T: Copy + RealField>(mesh: &Mesh<T>) -> Mesh<T> {
        // For now, return a clone; full hex decomposition can be added later.
        mesh.clone()
    }
}
