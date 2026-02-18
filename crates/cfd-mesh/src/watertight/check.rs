//! Watertight checking.

use crate::core::error::{MeshError, MeshResult};
use crate::storage::edge_store::EdgeStore;
use crate::storage::face_store::FaceStore;
use crate::storage::vertex_pool::VertexPool;
use crate::topology::manifold;
use crate::topology::orientation;
use crate::geometry::measure;

/// Comprehensive watertight status report.
#[derive(Clone, Debug)]
pub struct WatertightReport {
    /// Is the mesh a closed 2-manifold (no boundary edges)?
    pub is_closed: bool,
    /// Number of boundary edges.
    pub boundary_edge_count: usize,
    /// Number of non-manifold edges.
    pub non_manifold_edge_count: usize,
    /// Is orientation consistent?
    pub orientation_consistent: bool,
    /// Signed volume (should be positive for outward-oriented mesh).
    pub signed_volume: f64,
    /// Is the mesh watertight (all checks pass)?
    pub is_watertight: bool,
}

/// Check if a mesh is watertight.
pub fn check_watertight(
    vertex_pool: &VertexPool,
    face_store: &FaceStore,
    edge_store: &EdgeStore,
) -> WatertightReport {
    let manifold_report = manifold::check_manifold(edge_store);
    let orientation_ok = orientation::check_orientation(face_store, edge_store).is_ok();

    // Compute signed volume
    let signed_vol = measure::total_signed_volume(
        face_store.iter_enumerated().map(|(_, face)| {
            (
                vertex_pool.position(face.vertices[0]),
                vertex_pool.position(face.vertices[1]),
                vertex_pool.position(face.vertices[2]),
            )
        }),
    );

    let is_closed = manifold_report.is_closed_manifold;

    WatertightReport {
        is_closed,
        boundary_edge_count: manifold_report.boundary_edges,
        non_manifold_edge_count: manifold_report.non_manifold_edges,
        orientation_consistent: orientation_ok,
        signed_volume: signed_vol as f64,
        is_watertight: is_closed && orientation_ok,
    }
}

/// Assert the mesh is watertight, returning an error if not.
pub fn assert_watertight(
    vertex_pool: &VertexPool,
    face_store: &FaceStore,
    edge_store: &EdgeStore,
) -> MeshResult<WatertightReport> {
    let report = check_watertight(vertex_pool, face_store, edge_store);
    if !report.is_watertight {
        return Err(MeshError::NotWatertight {
            count: report.boundary_edge_count,
        });
    }
    Ok(report)
}
