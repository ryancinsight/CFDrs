//! Tetrahedral mesh utility functions for FEM element assembly.
//!
//! Provides vertex-index extraction (P1 and P2 variants), corner ordering
//! for positive Jacobian determinant, and bounding-box mesh scale computation.

use cfd_core::error::Result;
use cfd_mesh::domain::topology::Cell;
use cfd_mesh::IndexedMesh;
use nalgebra::{RealField, Vector3};
use num_traits::{Float, FromPrimitive};
use std::collections::HashSet;

/// Extract vertex indices from a cell for element assembly.
///
/// # Theorem — Canonical Tet10 Ordering
///
/// For a conforming P2 tetrahedral mesh, the returned 10-node index vector
/// `[c0, c1, c2, c3, m01, m12, m20, m03, m13, m23]` satisfies:
///
/// 1. The first 4 entries are corner nodes ordered so that
///    `det(J) = (p1−p0)×(p2−p0)·(p3−p0) > 0` (positive Jacobian).
/// 2. Each mid-edge node `m_ij` is the vertex closest to `(p_i + p_j) / 2`
///    among all non-corner vertices, guaranteeing geometric uniqueness on
///    conforming meshes.
///
/// **Proof sketch**: Positive-Jacobian ordering is enforced by `order_tet_corners`
/// which evaluates all 24 permutations. Mid-edge assignment by nearest-midpoint
/// is unique when the mesh is P2-conforming (each edge has exactly one mid-node
/// at its geometric centre).
pub fn extract_vertex_indices<T: cfd_mesh::domain::core::Scalar + RealField + Copy + Float>(
    cell: &Cell,
    mesh: &IndexedMesh<T>,
    n_corner_nodes: usize,
) -> Result<Vec<usize>> {
    // Pre-allocate for tet element: at most 10 unique vertices (4 corners + 6 mid-edges).
    let mut counts = std::collections::HashMap::with_capacity(10);
    for &f_idx in &cell.faces {
        if f_idx < mesh.face_count() {
            let f = mesh
                .faces
                .get(cfd_mesh::domain::core::index::FaceId::from_usize(f_idx));
            for &v_id in &f.vertices {
                *counts.entry(v_id.as_usize()).or_insert(0) += 1;
            }
        }
    }

    // In a tetrahedron, corners are shared by 3 faces, mid-edges by 2.
    let mut corners = Vec::with_capacity(4);
    let mut mid_edges = Vec::with_capacity(6);
    for (&v_idx, &count) in &counts {
        if count == 3 {
            corners.push(v_idx);
        } else if count == 2 {
            mid_edges.push(v_idx);
        }
    }

    if corners.len() == 4 && mid_edges.is_empty() && counts.len() == 4 {
        let ordered = order_tet_corners(&corners, mesh);

        if mesh.vertex_count() > n_corner_nodes {
            // P2 mesh: Face geometry only has corners, missing mid-edges.
            // Recover mid-edges via geometric search over extra nodes.
            let mut final_nodes = ordered.clone();
            let edges = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)];
            for &(i, j) in &edges {
                let v_i = ordered[i];
                let v_j = ordered[j];

                let selected = nearest_mid_node(
                    mesh,
                    n_corner_nodes,
                    v_i,
                    v_j,
                    n_corner_nodes..mesh.vertex_count(),
                    None,
                );
                final_nodes.push(selected);
            }
            return Ok(final_nodes);
        }
        // P1 Tet
        return Ok(ordered);
    }

    if corners.len() == 4 && mid_edges.len() == 6 {
        // P2 Tet
        // We need them in canonical Tet10 order for LagrangeTet10:
        // Corners: 0, 1, 2, 3
        // Mid-edges: 4:(0,1), 5:(1,2), 6:(2,0), 7:(0,3), 8:(1,3), 9:(2,3)
        let ordered = order_tet_corners(&corners, mesh);
        let mut final_nodes = ordered.clone();
        let mut used_mid_edges = std::collections::HashSet::new();

        let edges = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)];
        for &(i, j) in &edges {
            let v_i = ordered[i];
            let v_j = ordered[j];

            let selected = nearest_mid_node(
                mesh,
                n_corner_nodes,
                v_i,
                v_j,
                mid_edges.iter().copied(),
                Some(&used_mid_edges),
            );

            used_mid_edges.insert(selected);
            final_nodes.push(selected);
        }
        return Ok(final_nodes);
    }

    // Fallback for other elements (e.g. Hex)
    let mut all: Vec<usize> = counts.keys().copied().collect();
    all.sort_unstable();
    Ok(all)
}

/// Cache-accelerated variant of [`extract_vertex_indices`] for P2 FEM assembly.
///
/// # Performance (GAP-PERF-001)
///
/// When `mid_cache` is non-empty (P2 mesh), mid-node lookup for each of the 6
/// edges is O(1) amortised via the pre-built `HashMap<(usize,usize), usize>`,
/// replacing the O(n_mid) brute-force nearest-midpoint scan in `extract_vertex_indices`.
///
/// For P1 meshes (cache empty), falls back to the corner-only path identical to
/// the uncached version.
///
/// # Theorem — Correctness Equivalence
///
/// For any conforming P2 tetrahedral mesh where `MidNodeCache::build` was invoked
/// with the same mesh and `n_corner_nodes`, the output of `extract_vertex_indices_cached`
/// is element-wise identical to the output of `extract_vertex_indices`.
///
/// **Proof**: `MidNodeCache::build` maps each canonical edge `(min,max)` to the
/// unique mid-node closest to the geometric midpoint (geometric uniqueness from P2
/// conformity). `extract_vertex_indices` finds the same node via exhaustive nearest-
/// midpoint scan. Both converge to the same index by definition of the minimum.
pub fn extract_vertex_indices_cached<
    T: cfd_mesh::domain::core::Scalar + RealField + Copy + Float,
>(
    cell: &Cell,
    mesh: &IndexedMesh<T>,
    n_corner_nodes: usize,
    mid_cache: &crate::fem::mid_node_cache::MidNodeCache,
) -> Result<Vec<usize>> {
    // Pre-allocate for tet element: at most 10 unique vertices (4 corners + 6 mid-edges).
    let mut counts = std::collections::HashMap::with_capacity(10);
    for &f_idx in &cell.faces {
        if f_idx < mesh.face_count() {
            let f = mesh
                .faces
                .get(cfd_mesh::domain::core::index::FaceId::from_usize(f_idx));
            for &v_id in &f.vertices {
                *counts.entry(v_id.as_usize()).or_insert(0) += 1;
            }
        }
    }

    let mut corners = Vec::with_capacity(4);
    let mut mid_edges = Vec::with_capacity(6);
    for (&v_idx, &count) in &counts {
        if count == 3 {
            corners.push(v_idx);
        } else if count == 2 {
            mid_edges.push(v_idx);
        }
    }

    if corners.len() == 4 && mid_edges.is_empty() && counts.len() == 4 {
        let ordered = order_tet_corners(&corners, mesh);

        if mesh.vertex_count() > n_corner_nodes {
            // P2 mesh: use cache for O(1) mid-node lookup (GAP-PERF-001)
            let mut final_nodes = ordered.clone();
            let mut used_mid_edges = HashSet::new();
            let edges = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)];
            for &(i, j) in &edges {
                let v_i = ordered[i];
                let v_j = ordered[j];

                if let Some(m_idx) = mid_cache.get(v_i, v_j) {
                    used_mid_edges.insert(m_idx);
                    final_nodes.push(m_idx);
                } else {
                    // Cache miss: fall back to geometric search over the full
                    // mid-node range [n_corner_nodes, vertex_count).  In this
                    // branch mid_edges is EMPTY (the cell's P1-structured face
                    // data hasn't been enriched yet), so we must scan the
                    // global mid-node set rather than the empty local vec.
                    let selected = nearest_mid_node(
                        mesh,
                        n_corner_nodes,
                        v_i,
                        v_j,
                        n_corner_nodes..mesh.vertex_count(),
                        Some(&used_mid_edges),
                    );
                    used_mid_edges.insert(selected);
                    final_nodes.push(selected);
                }
            }
            return Ok(final_nodes);
        }
        // P1 Tet
        return Ok(ordered);
    }

    if corners.len() == 4 && mid_edges.len() == 6 {
        // P2 Tet: corners + mid-edges both recovered from face data
        let ordered = order_tet_corners(&corners, mesh);
        let mut final_nodes = ordered.clone();
        let mut used_mid_edges = std::collections::HashSet::new();

        let edges = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)];
        for &(i, j) in &edges {
            let v_i = ordered[i];
            let v_j = ordered[j];

            // Use cache for O(1) lookup if available
            if let Some(m_idx) = mid_cache.get(v_i, v_j) {
                used_mid_edges.insert(m_idx);
                final_nodes.push(m_idx);
                continue;
            }

            // Fallback: geometric search among mid_edges
            let p_i = mesh
                .vertices
                .position(cfd_mesh::domain::core::index::VertexId::from_usize(v_i))
                .coords;
            let p_j = mesh
                .vertices
                .position(cfd_mesh::domain::core::index::VertexId::from_usize(v_j))
                .coords;
            let target = (p_i + p_j)
                * <T as FromPrimitive>::from_f64(0.5)
                    .expect("0.5 is exactly representable in IEEE 754");

            let mut best_node = None;
            let mut min_dist = T::infinity();
            for &m_idx in &mid_edges {
                if used_mid_edges.contains(&m_idx) {
                    continue;
                }
                let dist = (mesh
                    .vertices
                    .position(cfd_mesh::domain::core::index::VertexId::from_usize(m_idx))
                    .coords
                    - target)
                    .norm();
                if dist < min_dist {
                    min_dist = dist;
                    best_node = Some(m_idx);
                }
            }
            if let Some(m_idx) = best_node {
                used_mid_edges.insert(m_idx);
                final_nodes.push(m_idx);
            }
        }
        return Ok(final_nodes);
    }

    // Fallback for other elements (e.g. Hex)
    let mut all: Vec<usize> = counts.keys().copied().collect();
    all.sort_unstable();
    Ok(all)
}

/// Order tetrahedron corners for positive Jacobian determinant.
///
/// # Theorem — Positive-Orientation Ordering
///
/// Given 4 corner vertices of a non-degenerate tetrahedron, at least one of the
/// 24 permutations `(v0, v1, v2, v3)` satisfies:
///
/// ```text
/// det(J) = (p1 − p0) × (p2 − p0) · (p3 − p0) > 0
/// ```
///
/// **Proof**: For a non-degenerate tetrahedron `det(J) ≠ 0`. Since swapping any
/// two vertices negates the determinant, exactly 12 of the 24 permutations yield
/// `det > 0`. This function selects the lexicographically smallest such permutation.
fn order_tet_corners<T: cfd_mesh::domain::core::Scalar + RealField + Copy + Float>(
    corners: &[usize],
    mesh: &IndexedMesh<T>,
) -> Vec<usize> {
    let perms: [[usize; 4]; 24] = [
        [0, 1, 2, 3],
        [0, 1, 3, 2],
        [0, 2, 1, 3],
        [0, 2, 3, 1],
        [0, 3, 1, 2],
        [0, 3, 2, 1],
        [1, 0, 2, 3],
        [1, 0, 3, 2],
        [1, 2, 0, 3],
        [1, 2, 3, 0],
        [1, 3, 0, 2],
        [1, 3, 2, 0],
        [2, 0, 1, 3],
        [2, 0, 3, 1],
        [2, 1, 0, 3],
        [2, 1, 3, 0],
        [2, 3, 0, 1],
        [2, 3, 1, 0],
        [3, 0, 1, 2],
        [3, 0, 2, 1],
        [3, 1, 0, 2],
        [3, 1, 2, 0],
        [3, 2, 0, 1],
        [3, 2, 1, 0],
    ];

    let mut best: Option<Vec<usize>> = None;
    let mut best_det = T::neg_infinity();

    for perm in &perms {
        let v0 = corners[perm[0]];
        let v1 = corners[perm[1]];
        let v2 = corners[perm[2]];
        let v3 = corners[perm[3]];

        let p0 = mesh
            .vertices
            .position(cfd_mesh::domain::core::index::VertexId::from_usize(v0))
            .coords;
        let p1 = mesh
            .vertices
            .position(cfd_mesh::domain::core::index::VertexId::from_usize(v1))
            .coords;
        let p2 = mesh
            .vertices
            .position(cfd_mesh::domain::core::index::VertexId::from_usize(v2))
            .coords;
        let p3 = mesh
            .vertices
            .position(cfd_mesh::domain::core::index::VertexId::from_usize(v3))
            .coords;

        let det = (p1 - p0).cross(&(p2 - p0)).dot(&(p3 - p0));
        if det > T::zero() {
            let candidate = vec![v0, v1, v2, v3];
            let take = match &best {
                None => true,
                Some(existing) => candidate < *existing,
            };
            if take {
                best = Some(candidate);
                best_det = det;
            }
        } else if best.is_none() && det > best_det {
            best_det = det;
            best = Some(vec![v0, v1, v2, v3]);
        }
    }

    best.unwrap_or_else(|| corners.to_vec())
}

/// Compute the bounding-box diagonal of a mesh (characteristic length scale).
///
/// Returns `‖max − min‖₂` where `min` and `max` are the component-wise
/// extrema of all vertex positions.
pub(crate) fn compute_mesh_scale<T: cfd_mesh::domain::core::Scalar + RealField + Copy + Float>(
    mesh: &IndexedMesh<T>,
) -> T {
    let mut min = Vector3::new(T::infinity(), T::infinity(), T::infinity());
    let mut max = Vector3::new(T::neg_infinity(), T::neg_infinity(), T::neg_infinity());
    for v in mesh.vertices.iter() {
        let p = v.1.position.coords;
        min.x = Float::min(min.x, p.x);
        min.y = Float::min(min.y, p.y);
        min.z = Float::min(min.z, p.z);
        max.x = Float::max(max.x, p.x);
        max.y = Float::max(max.y, p.y);
        max.z = Float::max(max.z, p.z);
    }
    (max - min).norm()
}

fn nearest_mid_node<T, I>(
    mesh: &IndexedMesh<T>,
    n_corner_nodes: usize,
    v_i: usize,
    v_j: usize,
    candidates: I,
    excluded: Option<&HashSet<usize>>,
) -> usize
where
    T: cfd_mesh::domain::core::Scalar + RealField + Copy + Float,
    I: IntoIterator<Item = usize> + Clone,
{
    use cfd_mesh::domain::core::index::VertexId;

    let p_i = mesh.vertices.position(VertexId::from_usize(v_i)).coords;
    let p_j = mesh.vertices.position(VertexId::from_usize(v_j)).coords;
    let target = (p_i + p_j)
        * <T as FromPrimitive>::from_f64(0.5).expect("0.5 is exactly representable in IEEE 754");

    let search = |skip_used: bool| -> Option<usize> {
        let mut best_node = None;
        let mut min_dist_sq = T::infinity();

        for m_idx in candidates.clone() {
            if skip_used && excluded.is_some_and(|used| used.contains(&m_idx)) {
                continue;
            }

            let pm = mesh.vertices.position(VertexId::from_usize(m_idx)).coords;
            let dist_sq = (pm - target).norm_squared();
            if dist_sq < min_dist_sq {
                min_dist_sq = dist_sq;
                best_node = Some(m_idx);
            }
        }

        best_node
    };

    search(true)
        .or_else(|| search(false))
        .unwrap_or(n_corner_nodes)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fem::mid_node_cache::MidNodeCache;
    use cfd_mesh::domain::topology::Cell;
    use cfd_mesh::IndexedMesh;
    use nalgebra::Point3;

    /// Build a minimal P1 tet mesh: 4 vertices, 4 triangular faces, 1 cell.
    fn build_single_tet_mesh() -> IndexedMesh<f64> {
        let mut mesh = IndexedMesh::<f64>::new();
        let v0 = mesh.add_vertex_pos(Point3::new(0.0, 0.0, 0.0));
        let v1 = mesh.add_vertex_pos(Point3::new(1.0, 0.0, 0.0));
        let v2 = mesh.add_vertex_pos(Point3::new(0.0, 1.0, 0.0));
        let v3 = mesh.add_vertex_pos(Point3::new(0.0, 0.0, 1.0));

        let f0 = mesh.add_face(v0, v1, v2);
        let f1 = mesh.add_face(v0, v1, v3);
        let f2 = mesh.add_face(v0, v2, v3);
        let f3 = mesh.add_face(v1, v2, v3);

        let cell = Cell::tetrahedron(f0.as_usize(), f1.as_usize(), f2.as_usize(), f3.as_usize());
        mesh.cells.push(cell);
        mesh
    }

    /// Build a minimal conforming P2 tet mesh with the same cell as the P1 fixture.
    fn build_single_p2_tet_mesh() -> IndexedMesh<f64> {
        let mut mesh = build_single_tet_mesh();

        let p01 = Point3::new(0.5, 0.0, 0.0);
        let p12 = Point3::new(0.5, 0.5, 0.0);
        let p20 = Point3::new(0.0, 0.5, 0.0);
        let p03 = Point3::new(0.0, 0.0, 0.5);
        let p13 = Point3::new(0.5, 0.0, 0.5);
        let p23 = Point3::new(0.0, 0.5, 0.5);

        mesh.add_vertex_pos(p01);
        mesh.add_vertex_pos(p12);
        mesh.add_vertex_pos(p20);
        mesh.add_vertex_pos(p03);
        mesh.add_vertex_pos(p13);
        mesh.add_vertex_pos(p23);

        mesh
    }

    #[test]
    fn extract_p1_tet_returns_four_corner_indices() {
        let mesh = build_single_tet_mesh();
        let cell = &mesh.cells[0];
        let n_corner_nodes = mesh.vertex_count(); // P1: all vertices are corners

        let indices = extract_vertex_indices(cell, &mesh, n_corner_nodes)
            .expect("extract_vertex_indices should succeed for a valid P1 tet");

        assert_eq!(
            indices.len(),
            4,
            "P1 tet must return exactly 4 corner indices"
        );
    }

    #[test]
    fn extracted_indices_are_within_mesh_bounds() {
        let mesh = build_single_tet_mesh();
        let cell = &mesh.cells[0];
        let n_corner_nodes = mesh.vertex_count();

        let indices = extract_vertex_indices(cell, &mesh, n_corner_nodes)
            .expect("extract_vertex_indices should succeed");

        for &idx in &indices {
            assert!(
                idx < mesh.vertex_count(),
                "vertex index {} out of bounds (mesh has {} vertices)",
                idx,
                mesh.vertex_count()
            );
        }
    }

    #[test]
    fn degenerate_cell_with_no_valid_faces_returns_empty() {
        let mesh = build_single_tet_mesh();
        // Create a cell whose face indices are all out-of-range.
        let bogus_cell = Cell::tetrahedron(999, 1000, 1001, 1002);

        let indices = extract_vertex_indices(&bogus_cell, &mesh, mesh.vertex_count())
            .expect("should not error, just return an empty fallback");

        assert!(
            indices.is_empty(),
            "cell with out-of-range faces should yield no vertex indices"
        );
    }

    #[test]
    fn compute_mesh_scale_nonzero_for_unit_tet() {
        let mesh = build_single_tet_mesh();
        let scale = compute_mesh_scale(&mesh);
        assert!(
            scale > 0.0,
            "mesh scale of a non-degenerate tet must be positive"
        );
        // The bounding box diagonal of the unit tet is sqrt(1^2+1^2+1^2) = sqrt(3).
        let expected = 3.0_f64.sqrt();
        assert!(
            (scale - expected).abs() < 1e-10,
            "expected mesh scale ~{}, got {}",
            expected,
            scale
        );
    }

    #[test]
    fn cached_and_uncached_p2_extraction_match() {
        let mesh = build_single_p2_tet_mesh();
        let cell = &mesh.cells[0];
        let n_corner_nodes = 4;
        let mid_cache = MidNodeCache::build(&mesh, n_corner_nodes);

        let uncached = extract_vertex_indices(cell, &mesh, n_corner_nodes)
            .expect("uncached extraction should succeed");
        let cached = extract_vertex_indices_cached(cell, &mesh, n_corner_nodes, &mid_cache)
            .expect("cached extraction should succeed");

        assert_eq!(uncached.len(), 10);
        assert_eq!(uncached, cached);
    }
}
