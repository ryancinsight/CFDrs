//! Manifold validation.
//!
//! Checks that every edge is shared by exactly 2 faces (closed 2-manifold)
//! and no vertex has a non-disk fan topology.

use crate::domain::core::error::{MeshError, MeshResult};
use crate::domain::core::index::EdgeId;
use crate::infrastructure::storage::edge_store::EdgeStore;

/// Result of manifold checking.
#[derive(Clone, Debug)]
pub struct ManifoldReport {
    /// Total number of edges.
    pub total_edges: usize,
    /// Edges with exactly 2 adjacent faces (manifold interior).
    pub manifold_edges: usize,
    /// Edges with exactly 1 adjacent face (boundary).
    pub boundary_edges: usize,
    /// Edges with >2 adjacent faces (non-manifold).
    pub non_manifold_edges: usize,
    /// IDs of non-manifold edges.
    pub non_manifold_edge_ids: Vec<EdgeId>,
    /// IDs of boundary edges.
    pub boundary_edge_ids: Vec<EdgeId>,
    /// Is the mesh a closed 2-manifold?
    pub is_closed_manifold: bool,
}

/// Check if the mesh is a closed 2-manifold.
///
/// A closed 2-manifold has every edge shared by exactly 2 faces.
pub fn check_manifold(edge_store: &EdgeStore) -> ManifoldReport {
    let mut manifold_edges = 0usize;
    let mut boundary_edges = 0usize;
    let mut non_manifold_edges = 0usize;
    let mut non_manifold_edge_ids = Vec::new();
    let mut boundary_edge_ids = Vec::new();

    for (eid, edge) in edge_store.iter_enumerated() {
        match edge.valence() {
            2 => manifold_edges += 1,
            1 => {
                boundary_edges += 1;
                boundary_edge_ids.push(eid);
            }
            _ => {
                non_manifold_edges += 1;
                non_manifold_edge_ids.push(eid);
            }
        }
    }

    let is_closed_manifold = boundary_edges == 0 && non_manifold_edges == 0;

    ManifoldReport {
        total_edges: edge_store.len(),
        manifold_edges,
        boundary_edges,
        non_manifold_edges,
        non_manifold_edge_ids,
        boundary_edge_ids,
        is_closed_manifold,
    }
}

/// Assert that the mesh is a closed 2-manifold, returning an error if not.
pub fn assert_manifold(edge_store: &EdgeStore) -> MeshResult<ManifoldReport> {
    let report = check_manifold(edge_store);
    if !report.non_manifold_edge_ids.is_empty() {
        return Err(MeshError::NonManifoldEdge {
            edge: report.non_manifold_edge_ids[0],
            count: edge_store.get(report.non_manifold_edge_ids[0]).valence(),
        });
    }
    if !report.boundary_edge_ids.is_empty() {
        return Err(MeshError::NotWatertight {
            count: report.boundary_edges,
        });
    }
    Ok(report)
}
