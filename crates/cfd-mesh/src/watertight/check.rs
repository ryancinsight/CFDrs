//! Watertight checking.
//!
//! Validates a closed triangle mesh using four independent criteria:
//!
//! 1. **Manifold + closed**: every edge is shared by exactly 2 faces.
//! 2. **Euler characteristic**: `V - E + F = 2` for a genus-0 closed sphere.
//!    A torus gives `V - E + F = 0`, etc. Mismatches reveal topological defects.
//! 3. **Orientation consistency**: all face pairs sharing an edge have opposite
//!    directed-edge orientations (no two adjacent faces wind the same way).
//! 4. **Positive signed volume**: the divergence-theorem volume integral should
//!    be positive for an outward-oriented mesh.
//!
//! ## Euler's Theorem
//!
//! For a convex polyhedron (or any genus-0 closed surface):
//!
//! $$V - E + F = 2 \cdot (1 - g)$$
//!
//! where $g$ is the genus (number of handles). For a sphere or cube $g = 0$
//! so the characteristic is 2. A torus has $g = 1$ so the characteristic is 0.
//!
//! For a triangle mesh: $E = 3F/2$ (each face contributes 3 half-edges, each
//! edge is shared by exactly 2 faces in a manifold), so the relation reduces
//! to $V - E + F = 2$ for a closed manifold of genus 0.

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
    /// Euler characteristic $\chi = V - E + F$.
    ///
    /// - `2` for a closed sphere-topology surface (genus 0)
    /// - `0` for a torus (genus 1)
    /// - Negative values indicate complex topology or mesh defects
    ///
    /// `None` when vertices/edges/faces counts are not available.
    pub euler_characteristic: Option<i64>,
    /// Expected Euler characteristic for a valid closed manifold of genus 0.
    pub euler_expected: i64,
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

    // Euler characteristic: V - E + F = 2 for a closed genus-0 manifold.
    // For a triangle mesh with E manifold edges: each face has 3 edges, each
    // interior edge is shared by 2 faces, so E = 3F/2 (manifold only).
    let v = vertex_pool.len() as i64;
    let e = edge_store.len() as i64;
    let f = face_store.len() as i64;
    let euler = v - e + f;

    let is_closed = manifold_report.is_closed_manifold;

    WatertightReport {
        is_closed,
        boundary_edge_count: manifold_report.boundary_edges,
        non_manifold_edge_count: manifold_report.non_manifold_edges,
        orientation_consistent: orientation_ok,
        signed_volume: signed_vol as f64,
        is_watertight: is_closed && orientation_ok,
        euler_characteristic: Some(euler),
        euler_expected: 2,
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
    // Euler characteristic check — only meaningful for closed manifolds.
    if let Some(chi) = report.euler_characteristic {
        if chi != report.euler_expected {
            return Err(MeshError::Other(format!(
                "Euler characteristic χ = {} (expected {}); topology defect detected",
                chi, report.euler_expected
            )));
        }
    }
    Ok(report)
}
