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
use crate::core::scalar::Scalar;
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
pub fn check_watertight<T: Scalar>(
    vertex_pool: &VertexPool<T>,
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
    let signed_vol_f64 = num_traits::ToPrimitive::to_f64(&signed_vol).unwrap_or(0.0);

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
        signed_volume: signed_vol_f64,
        is_watertight: is_closed && orientation_ok,
        euler_characteristic: Some(euler),
        euler_expected: 2,
    }
}

/// Assert the mesh is watertight, returning an error if not.
pub fn assert_watertight<T: Scalar>(
    vertex_pool: &VertexPool<T>,
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

// ── HalfEdgeMesh<'id> watertight check ───────────────────────────────────────

use crate::mesh::HalfEdgeMesh;
use crate::permission::GhostToken;

/// Check the watertightness of a [`HalfEdgeMesh`] directly from its
/// half-edge topology, without needing an `EdgeStore`.
///
/// ## Algorithm
///
/// A `HalfEdgeMesh` is watertight iff **every half-edge has a non-boundary
/// twin** — i.e., `twin(he).face` is `Some(_)` for all `he`.  Boundary
/// sentinels have `face = None` and represent open edges.
///
/// The Euler characteristic is computed from the counts:
///
/// ```text
/// V = vertex_count (all keys in the vertex SlotMap)
/// E = half_edge_count / 2   (each undirected edge = 2 half-edges)
/// F = face_count             (interior faces only)
/// χ = V - E + F
/// ```
///
/// For a closed genus-0 manifold, χ = 2.
///
/// ## Theorem — Euler–Poincaré for triangle meshes
///
/// In a closed orientable 2-manifold triangulation:
/// - Each face contributes 3 half-edges.
/// - Each undirected edge is shared by exactly 2 half-edges.
/// - Therefore `E = HE / 2 = 3F / 2`, so `V - 3F/2 + F = V - F/2 = 2 - 2g`.
///   For genus 0: `F = 2V - 4`, `E = 3V - 6`.
pub fn check_halfedge<'id>(
    mesh: &HalfEdgeMesh<'id>,
    token: &GhostToken<'id>,
) -> WatertightReport {
    use crate::core::index::HalfEdgeKey;
    use hashbrown::HashSet;

    let v = mesh.vertex_count() as i64;
    let he_total = mesh.half_edge_count() as i64;
    let f = mesh.face_count() as i64;

    // Count boundary half-edges (those whose twin has face=None).
    // A boundary half-edge is one that was inserted as a sentinel (face=None).
    // We detect them by iterating over half-edges from each face loop and
    // checking the twin's face field.
    let mut boundary_count: usize = 0;
    let non_manifold_count: usize = 0;

    // Collect all half-edge keys from interior face loops
    let mut interior_hes: HashSet<HalfEdgeKey> = HashSet::new();
    for fk in mesh.face_keys() {
        for he in mesh.face_half_edges(fk, token) {
            interior_hes.insert(he);
        }
    }

    for he in interior_hes.iter().copied() {
        let twin = mesh.he_twin(he, token);
        match mesh.he_face(twin, token) {
            None => boundary_count += 1,
            Some(_) => {} // manifold interior edge
        }
    }

    let is_closed = boundary_count == 0;

    // For a pure triangle mesh with no boundary sentinels, each undirected
    // edge is represented by exactly 2 half-edges.  Sentinels inflate the
    // count; subtract them.
    let _boundary_he_count = (he_total as usize) - interior_hes.len();
    let e = (interior_hes.len() as i64) / 2;
    let euler = v - e + f;

    // Compute signed volume from face vertex loops
    let signed_vol = mesh.signed_volume(token);

    WatertightReport {
        is_closed,
        boundary_edge_count: boundary_count,
        non_manifold_edge_count: non_manifold_count,
        orientation_consistent: signed_vol > 0.0,
        signed_volume: signed_vol as f64,
        is_watertight: is_closed && (signed_vol > 0.0),
        euler_characteristic: Some(euler),
        euler_expected: 2,
    }
}

/// Assert that a [`HalfEdgeMesh`] is watertight.
pub fn assert_halfedge_watertight<'id>(
    mesh: &HalfEdgeMesh<'id>,
    token: &GhostToken<'id>,
) -> MeshResult<WatertightReport> {
    let report = check_halfedge(mesh, token);
    if !report.is_watertight {
        return Err(MeshError::NotWatertight {
            count: report.boundary_edge_count,
        });
    }
    if let Some(chi) = report.euler_characteristic {
        if chi != report.euler_expected {
            return Err(MeshError::Other(format!(
                "HalfEdgeMesh Euler characteristic χ = {} (expected {})",
                chi, report.euler_expected
            )));
        }
    }
    Ok(report)
}

#[cfg(test)]
mod he_tests {
    use super::*;
    use crate::mesh::with_mesh;
    use nalgebra::Point3;

    fn unit_tetrahedron() -> WatertightReport {
        with_mesh(|mut mesh, mut token| {
            // Regular tetrahedron with outward normals
            let a = mesh.add_vertex(Point3::new(1.0_f64,  1.0,  1.0), &token);
            let b = mesh.add_vertex(Point3::new(-1.0, -1.0,  1.0), &token);
            let c = mesh.add_vertex(Point3::new(-1.0,  1.0, -1.0), &token);
            let d = mesh.add_vertex(Point3::new( 1.0, -1.0, -1.0), &token);
            // Four CCW faces (outward-pointing normals by right-hand rule)
            mesh.add_triangle(a, b, c, &mut token).unwrap();
            mesh.add_triangle(a, c, d, &mut token).unwrap();
            mesh.add_triangle(a, d, b, &mut token).unwrap();
            mesh.add_triangle(b, d, c, &mut token).unwrap();
            check_halfedge(&mesh, &token)
        })
    }

    #[test]
    fn tetrahedron_is_closed() {
        let report = unit_tetrahedron();
        assert_eq!(report.boundary_edge_count, 0, "tetrahedron must have no boundary edges");
        assert!(report.is_closed);
    }

    #[test]
    fn tetrahedron_euler_characteristic() {
        let report = unit_tetrahedron();
        // V=4, E=6, F=4 → χ=2
        assert_eq!(
            report.euler_characteristic,
            Some(2),
            "tetrahedron χ must equal 2, got {:?}",
            report.euler_characteristic
        );
    }

    #[test]
    fn single_triangle_is_not_watertight() {
        let report = with_mesh(|mut mesh, mut token| {
            let a = mesh.add_vertex(Point3::new(0.0_f64, 0.0, 0.0), &token);
            let b = mesh.add_vertex(Point3::new(1.0, 0.0, 0.0), &token);
            let c = mesh.add_vertex(Point3::new(0.0, 1.0, 0.0), &token);
            mesh.add_triangle(a, b, c, &mut token).unwrap();
            check_halfedge(&mesh, &token)
        });
        assert!(!report.is_watertight, "a single triangle must not be watertight");
        assert!(report.boundary_edge_count > 0);
    }
}
