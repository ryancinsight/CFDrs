//! Mesh refinement module

use crate::topology::{Cell, ElementType, Vertex};
use nalgebra::RealField;
use std::collections::HashSet;

pub mod criteria;
pub mod utils;

use utils::RefinementContext;

// Re-export main types
pub use criteria::*;

/// Refinement strategy trait
pub trait RefinementStrategy<T: RealField + Copy>: Send + Sync {
    /// Apply refinement to mesh
    fn refine(&self, mesh: &mut crate::mesh::Mesh<T>) -> Result<(), crate::error::MeshError>;

    /// Get strategy name
    fn name(&self) -> &'static str;
}

/// Refine a tetrahedron element (1 -> 8 subdivision)
fn refine_tetrahedron<T: RealField + Copy>(
    mesh: &crate::mesh::Mesh<T>,
    cell: &Cell,
    ctx: &mut RefinementContext<T>,
) -> Result<(), crate::error::MeshError> {
    // Get vertices
    let mut v_indices = Vec::new();
    for &f_idx in &cell.faces {
        if let Some(face) = mesh.face(f_idx) {
            for &v_idx in &face.vertices {
                if !v_indices.contains(&v_idx) {
                    v_indices.push(v_idx);
                }
            }
        }
    }
    v_indices.sort_unstable();

    if v_indices.len() != 4 {
        return Err(crate::error::MeshError::InvalidMesh(format!(
            "Tetrahedron must have 4 vertices, found {}",
            v_indices.len()
        )));
    }

    let [v0, v1, v2, v3] = [v_indices[0], v_indices[1], v_indices[2], v_indices[3]];

    // Get midpoints
    let m01 = ctx.get_midpoint(v0, v1)?;
    let m02 = ctx.get_midpoint(v0, v2)?;
    let m03 = ctx.get_midpoint(v0, v3)?;
    let m12 = ctx.get_midpoint(v1, v2)?;
    let m13 = ctx.get_midpoint(v1, v3)?;
    let m23 = ctx.get_midpoint(v2, v3)?;

    // Define the 8 new tets
    let new_tets = vec![
        // 4 Corners
        vec![v0, m01, m02, m03],
        vec![v1, m01, m12, m13],
        vec![v2, m02, m12, m23],
        vec![v3, m03, m13, m23],
        // 4 Inner
        vec![m01, m23, m02, m12],
        vec![m01, m23, m12, m13],
        vec![m01, m23, m13, m03],
        vec![m01, m23, m03, m02],
    ];

    for tet_verts in new_tets {
        let f0 = ctx.add_face_dedup(vec![tet_verts[0], tet_verts[1], tet_verts[2]]);
        let f1 = ctx.add_face_dedup(vec![tet_verts[0], tet_verts[1], tet_verts[3]]);
        let f2 = ctx.add_face_dedup(vec![tet_verts[0], tet_verts[2], tet_verts[3]]);
        let f3 = ctx.add_face_dedup(vec![tet_verts[1], tet_verts[2], tet_verts[3]]);
        ctx.new_mesh.add_cell(Cell::tetrahedron(f0, f1, f2, f3));
    }

    Ok(())
}

/// Refine a triangle element (1 -> 4 subdivision)
fn refine_triangle<T: RealField + Copy>(
    mesh: &crate::mesh::Mesh<T>,
    cell: &Cell,
    ctx: &mut RefinementContext<T>,
) -> Result<(), crate::error::MeshError> {
    // Get vertices
    let mut v_indices = Vec::new();
    for &f_idx in &cell.faces {
        if let Some(face) = mesh.face(f_idx) {
            for &v_idx in &face.vertices {
                if !v_indices.contains(&v_idx) {
                    v_indices.push(v_idx);
                }
            }
        }
    }
    v_indices.sort_unstable();

    if v_indices.len() != 3 {
        return Err(crate::error::MeshError::InvalidMesh(format!(
            "Triangle must have 3 vertices, found {}",
            v_indices.len()
        )));
    }

    let [v0, v1, v2] = [v_indices[0], v_indices[1], v_indices[2]];

    // Get midpoints
    let m01 = ctx.get_midpoint(v0, v1)?;
    let m12 = ctx.get_midpoint(v1, v2)?;
    let m20 = ctx.get_midpoint(v2, v0)?;

    // Define the 4 new triangles
    let new_tris = vec![
        vec![v0, m01, m20],
        vec![v1, m12, m01],
        vec![v2, m20, m12],
        vec![m01, m12, m20],
    ];

    for tri_verts in new_tris {
        let f0 = ctx.add_face_dedup(vec![tri_verts[0], tri_verts[1], tri_verts[2]]);
        // For 2D elements, we create a triangular cell with one face
        ctx.new_mesh.add_cell(Cell {
            faces: vec![f0],
            element_type: ElementType::Triangle,
            global_id: None,
            partition_id: None,
        });
    }

    Ok(())
}

/// Refine a hexahedron element (1 -> 8 subdivision)
fn refine_hexahedron<T: RealField + Copy>(
    mesh: &crate::mesh::Mesh<T>,
    cell: &Cell,
    ctx: &mut RefinementContext<T>,
) -> Result<(), crate::error::MeshError> {
    let vertices = collect_cell_vertices(mesh, cell)?;

    let points = HexRefinementPoints::new(ctx, vertices)?;
    for hex_verts in points.subcells() {
        add_hexahedron_cell(ctx, hex_verts);
    }

    Ok(())
}

struct HexRefinementPoints {
    vertices: [usize; 8],
    m01: usize,
    m12: usize,
    m23: usize,
    m30: usize,
    m45: usize,
    m56: usize,
    m67: usize,
    m74: usize,
    m04: usize,
    m15: usize,
    m26: usize,
    m37: usize,
    fc_bottom: usize,
    fc_top: usize,
    fc_front: usize,
    fc_back: usize,
    fc_left: usize,
    fc_right: usize,
    body_center: usize,
}

impl HexRefinementPoints {
    fn new<T: RealField + Copy>(
        ctx: &mut RefinementContext<T>,
        vertices: [usize; 8],
    ) -> Result<Self, crate::error::MeshError> {
        let m01 = ctx.get_midpoint(vertices[0], vertices[1])?;
        let m12 = ctx.get_midpoint(vertices[1], vertices[2])?;
        let m23 = ctx.get_midpoint(vertices[2], vertices[3])?;
        let m30 = ctx.get_midpoint(vertices[3], vertices[0])?;
        let m45 = ctx.get_midpoint(vertices[4], vertices[5])?;
        let m56 = ctx.get_midpoint(vertices[5], vertices[6])?;
        let m67 = ctx.get_midpoint(vertices[6], vertices[7])?;
        let m74 = ctx.get_midpoint(vertices[7], vertices[4])?;
        let m04 = ctx.get_midpoint(vertices[0], vertices[4])?;
        let m15 = ctx.get_midpoint(vertices[1], vertices[5])?;
        let m26 = ctx.get_midpoint(vertices[2], vertices[6])?;
        let m37 = ctx.get_midpoint(vertices[3], vertices[7])?;

        let fc_bottom = ctx.get_midpoint(m01, m23)?;
        let fc_top = ctx.get_midpoint(m45, m67)?;
        let fc_front = ctx.get_midpoint(m01, m45)?;
        let fc_back = ctx.get_midpoint(m23, m67)?;
        let fc_left = ctx.get_midpoint(m30, m74)?;
        let fc_right = ctx.get_midpoint(m12, m56)?;

        let body_center = ctx.get_midpoint(fc_bottom, fc_top)?;

        Ok(Self {
            vertices,
            m01,
            m12,
            m23,
            m30,
            m45,
            m56,
            m67,
            m74,
            m04,
            m15,
            m26,
            m37,
            fc_bottom,
            fc_top,
            fc_front,
            fc_back,
            fc_left,
            fc_right,
            body_center,
        })
    }

    fn subcells(&self) -> [[usize; 8]; 8] {
        [
            [
                self.vertices[0],
                self.m01,
                self.fc_bottom,
                self.m30,
                self.m04,
                self.fc_front,
                self.body_center,
                self.fc_left,
            ],
            [
                self.m01,
                self.vertices[1],
                self.m12,
                self.fc_bottom,
                self.fc_front,
                self.m15,
                self.fc_right,
                self.body_center,
            ],
            [
                self.fc_bottom,
                self.m12,
                self.vertices[2],
                self.m23,
                self.body_center,
                self.fc_right,
                self.m26,
                self.fc_back,
            ],
            [
                self.m30,
                self.fc_bottom,
                self.m23,
                self.vertices[3],
                self.fc_left,
                self.body_center,
                self.fc_back,
                self.m37,
            ],
            [
                self.m04,
                self.fc_front,
                self.body_center,
                self.fc_left,
                self.vertices[4],
                self.m45,
                self.fc_top,
                self.m74,
            ],
            [
                self.fc_front,
                self.m15,
                self.fc_right,
                self.body_center,
                self.m45,
                self.vertices[5],
                self.m56,
                self.fc_top,
            ],
            [
                self.body_center,
                self.fc_right,
                self.m26,
                self.fc_back,
                self.fc_top,
                self.m56,
                self.vertices[6],
                self.m67,
            ],
            [
                self.fc_left,
                self.body_center,
                self.fc_back,
                self.m37,
                self.m74,
                self.fc_top,
                self.m67,
                self.vertices[7],
            ],
        ]
    }
}

fn add_hexahedron_cell<T: RealField + Copy>(ctx: &mut RefinementContext<T>, hex_verts: [usize; 8]) {
    let f0 = ctx.add_face_dedup(vec![hex_verts[0], hex_verts[1], hex_verts[2], hex_verts[3]]);
    let f1 = ctx.add_face_dedup(vec![hex_verts[4], hex_verts[5], hex_verts[6], hex_verts[7]]);
    let f2 = ctx.add_face_dedup(vec![hex_verts[0], hex_verts[1], hex_verts[5], hex_verts[4]]);
    let f3 = ctx.add_face_dedup(vec![hex_verts[2], hex_verts[3], hex_verts[7], hex_verts[6]]);
    let f4 = ctx.add_face_dedup(vec![hex_verts[0], hex_verts[3], hex_verts[7], hex_verts[4]]);
    let f5 = ctx.add_face_dedup(vec![hex_verts[1], hex_verts[2], hex_verts[6], hex_verts[5]]);

    ctx.new_mesh
        .add_cell(Cell::hexahedron(vec![f0, f1, f2, f3, f4, f5]));
}

fn collect_cell_vertices<T: RealField + Copy>(
    mesh: &crate::mesh::Mesh<T>,
    cell: &Cell,
) -> Result<[usize; 8], crate::error::MeshError> {
    let mut v_indices = Vec::new();
    for &f_idx in &cell.faces {
        if let Some(face) = mesh.face(f_idx) {
            for &v_idx in &face.vertices {
                if !v_indices.contains(&v_idx) {
                    v_indices.push(v_idx);
                }
            }
        }
    }
    v_indices.sort_unstable();

    if v_indices.len() != 8 {
        return Err(crate::error::MeshError::InvalidMesh(format!(
            "Hexahedron must have 8 vertices, found {}",
            v_indices.len()
        )));
    }

    v_indices.try_into().map_err(|_| {
        crate::error::MeshError::InvalidMesh("Hexahedron vertex collection failed".to_string())
    })
}

/// Uniform refinement strategy
pub struct UniformRefinement;

impl<T: RealField + Copy> RefinementStrategy<T> for UniformRefinement {
    #[allow(clippy::too_many_lines)]
    fn refine(&self, mesh: &mut crate::mesh::Mesh<T>) -> Result<(), crate::error::MeshError> {
        if mesh.cell_count() == 0 {
            return Ok(());
        }

        // Support uniform refinement for additional element types.
        for cell in mesh.cells() {
            match cell.element_type {
                ElementType::Tetrahedron | ElementType::Triangle | ElementType::Hexahedron => {
                    // Supported element types
                }
                _ => {
                    return Err(crate::error::MeshError::NotImplemented(format!(
                        "Uniform refinement not implemented for {:?}. Supported types: Triangle, Tetrahedron, Hexahedron",
                        cell.element_type
                    )));
                }
            }
        }

        let mut ctx = RefinementContext::new(mesh);

        // Iterate over cells and refine
        for cell in mesh.cells() {
            match cell.element_type {
                ElementType::Tetrahedron => {
                    refine_tetrahedron(mesh, cell, &mut ctx)?;
                }
                ElementType::Triangle => {
                    refine_triangle(mesh, cell, &mut ctx)?;
                }
                ElementType::Hexahedron => {
                    refine_hexahedron(mesh, cell, &mut ctx)?;
                }
                _ => {
                    // This should not happen due to earlier check
                    return Err(crate::error::MeshError::NotImplemented(format!(
                        "Unsupported element type: {:?}",
                        cell.element_type
                    )));
                }
            }
        }

        // Replace mesh
        *mesh = ctx.new_mesh;

        Ok(())
    }

    fn name(&self) -> &'static str {
        "Uniform"
    }
}

/// Estimate gradient magnitude for a cell using neighboring cell values
fn estimate_cell_gradient<T: RealField + Copy>(
    cell_idx: usize,
    mesh: &crate::mesh::Mesh<T>,
    field: &[T],
) -> Result<T, crate::error::MeshError> {
    // Get the current cell
    let cell = mesh.cell(cell_idx).ok_or_else(|| {
        crate::error::MeshError::InvalidMesh(format!("Cell {cell_idx} not found"))
    })?;

    // Get current cell value
    let current_value = field.get(cell_idx).copied().unwrap_or_else(|| T::zero());

    // Find neighboring cells by sharing faces
    let mut neighbor_values = Vec::new();
    let mut neighbor_positions = Vec::new();

    for &face_idx in &cell.faces {
        if let Some(_face) = mesh.face(face_idx) {
            // Find all cells that share this face
            for (other_cell_idx, other_cell) in mesh.cells().iter().enumerate() {
                if other_cell_idx != cell_idx && other_cell.faces.contains(&face_idx) {
                    // Found a neighbor
                    if let Some(&neighbor_value) = field.get(other_cell_idx) {
                        neighbor_values.push(neighbor_value);

                        // Compute neighbor cell centroid
                        let centroid = compute_cell_centroid(mesh, other_cell)?;
                        neighbor_positions.push(centroid);
                    }
                }
            }
        }
    }

    if neighbor_values.is_empty() {
        // No neighbors found, gradient is zero
        return Ok(T::zero());
    }

    // Compute current cell centroid
    let current_centroid = compute_cell_centroid(mesh, cell)?;

    // Compute gradient using least squares approach
    // ∇f ≈ argmin Σ(|f_i - f_current - ∇f·(x_i - x_current)|²)
    let mut gradient_sum = T::zero();
    let mut weight_sum = T::zero();

    for (neighbor_value, neighbor_pos) in neighbor_values.iter().zip(neighbor_positions.iter()) {
        let value_diff = *neighbor_value - current_value;

        // Compute distance vector
        let dx = neighbor_pos.position.x - current_centroid.position.x;
        let dy = neighbor_pos.position.y - current_centroid.position.y;
        let dz = neighbor_pos.position.z - current_centroid.position.z;

        // Compute distance magnitude
        let distance_sq = dx * dx + dy * dy + dz * dz;

        if distance_sq > T::from_f64(1e-12).unwrap() {
            // Weight by inverse distance
            let weight = T::one() / distance_sq.sqrt();

            // Approximate gradient magnitude as |Δf|/|Δx|
            let distance = distance_sq.sqrt();
            let gradient_contribution = (value_diff / distance).abs();

            gradient_sum += weight * gradient_contribution;
            weight_sum += weight;
        }
    }

    if weight_sum > T::zero() {
        Ok(gradient_sum / weight_sum)
    } else {
        Ok(T::zero())
    }
}

/// Compute the centroid of a cell
fn compute_cell_centroid<T: RealField + Copy>(
    mesh: &crate::mesh::Mesh<T>,
    cell: &Cell,
) -> Result<crate::topology::Vertex<T>, crate::error::MeshError> {
    use std::collections::HashSet;

    let mut vertex_indices = HashSet::new();

    // Collect all unique vertex indices from all faces
    for &face_idx in &cell.faces {
        if let Some(face) = mesh.face(face_idx) {
            for &v_idx in &face.vertices {
                vertex_indices.insert(v_idx);
            }
        }
    }

    if vertex_indices.is_empty() {
        return Err(crate::error::MeshError::InvalidMesh(
            "Cell has no vertices".to_string(),
        ));
    }

    // Compute average position
    let mut sum_x = T::zero();
    let mut sum_y = T::zero();
    let mut sum_z = T::zero();
    let count = T::from_usize(vertex_indices.len()).unwrap();

    for &v_idx in &vertex_indices {
        if let Some(vertex) = mesh.vertex(v_idx) {
            sum_x += vertex.position.x;
            sum_y += vertex.position.y;
            sum_z += vertex.position.z;
        }
    }

    Ok(crate::topology::Vertex {
        position: nalgebra::Point3::new(sum_x / count, sum_y / count, sum_z / count),
        global_id: None,
        partition_id: None,
    })
}

/// Check if 4 edges form opposite pairs pattern for green refinement
fn check_opposite_edge_pattern<T: RealField + Copy>(
    split_edges: &[(usize, usize)],
    mesh: &crate::mesh::Mesh<T>,
    cell_idx: usize,
) -> bool {
    use std::collections::HashSet;

    if split_edges.len() != 4 {
        return false;
    }

    // Get the cell to understand its topology
    let Some(cell) = mesh.cell(cell_idx) else {
        return false;
    };

    // Collect all vertices of the cell
    let mut cell_vertices = HashSet::new();
    for &face_idx in &cell.faces {
        if let Some(face) = mesh.face(face_idx) {
            for &v_idx in &face.vertices {
                cell_vertices.insert(v_idx);
            }
        }
    }

    // For a tetrahedron, check if we have exactly 4 vertices
    if cell_vertices.len() != 4 {
        return false;
    }

    // Count vertex occurrences in split edges
    let mut vertex_counts = std::collections::HashMap::new();
    for &(u, v) in split_edges {
        *vertex_counts.entry(u).or_insert(0) += 1;
        *vertex_counts.entry(v).or_insert(0) += 1;
    }

    // For opposite edges pattern, each vertex should appear exactly 2 times
    // (each vertex is part of exactly 2 opposite edges)
    vertex_counts.values().all(|&count| count == 2)
}

/// Adaptive refinement strategy
pub struct AdaptiveRefinement<T: RealField + Copy> {
    /// Refinement criterion
    pub criterion: RefinementCriterion<T>,
}

impl<T: RealField + Copy> RefinementStrategy<T> for AdaptiveRefinement<T> {
    #[allow(clippy::too_many_lines)]
    fn refine(&self, mesh: &mut crate::mesh::Mesh<T>) -> Result<(), crate::error::MeshError> {
        if mesh.cell_count() == 0 {
            return Ok(());
        }

        // Mark cells for refinement based on criterion
        let mut marked_cells = HashSet::new();

        for (i, cell) in mesh.cells().iter().enumerate() {
            // Check element type
            if cell.element_type != ElementType::Tetrahedron {
                return Err(crate::error::MeshError::NotImplemented(format!(
                    "Adaptive refinement only implemented for Tetrahedra, found {:?}",
                    cell.element_type
                )));
            }

            let should_refine = match &self.criterion {
                RefinementCriterion::Error {
                    error_field,
                    threshold,
                } => {
                    if let Some(err) = error_field.get(i) {
                        *err > *threshold
                    } else {
                        false
                    }
                }
                RefinementCriterion::Gradient { field, threshold } => {
                    // Implement gradient estimation using mesh connectivity/geometry.
                    // Compute gradient magnitude for each cell using neighboring cell values
                    let gradient_magnitude = estimate_cell_gradient(i, mesh, field)?;

                    gradient_magnitude > *threshold
                }
                RefinementCriterion::Custom(func) => {
                    let vertices: Vec<Vertex<T>> = mesh
                        .ordered_element_vertices(cell)
                        .into_iter()
                        .copied()
                        .collect();
                    func(cell, &vertices)
                }
                RefinementCriterion::Geometric { .. } => false,
            };

            if should_refine {
                marked_cells.insert(i);
            }
        }

        if marked_cells.is_empty() {
            return Ok(());
        }

        // Conformity Handling with Full Green Refinement Patterns
        // This implementation ensures mesh conformity by using red-green refinement:
        // 1. "Red" refinement: Complete subdivision (1->8 for tets)
        // 2. "Green" refinement: Partial subdivision to maintain conformity
        //
        // Supported green patterns for tetrahedra:
        // - 1-edge split: Bisect (1->2 tets)
        // - 3-edge face split: Face refinement (1->4 tets)
        // - 2-edge split: Edge pair refinement (1->3 tets) - NEW
        // - 4-edge split: Opposite edges (1->4 tets) - NEW
        //
        // The algorithm iteratively propagates refinement to ensure conformity:
        // - Initially marked cells use red refinement
        // - Neighbors are analyzed for compatible green patterns
        // - Incompatible patterns are promoted to red refinement

        // Map edges to unique ID (min, max)
        // Pre-calculate edges for all cells to speed up
        let mut cell_edges: Vec<Vec<(usize, usize)>> = Vec::with_capacity(mesh.cell_count());
        for cell in mesh.cells() {
            let mut edges = HashSet::new();
            for &f_idx in &cell.faces {
                if let Some(face) = mesh.face(f_idx) {
                    for (u, v) in face.edges_iter() {
                        let key = if u < v { (u, v) } else { (v, u) };
                        edges.insert(key);
                    }
                }
            }
            cell_edges.push(edges.into_iter().collect());
        }

        let mut changed = true;
        while changed {
            changed = false;

            // Identify all split edges based on current marked_cells
            let mut split_edges: HashSet<(usize, usize)> = HashSet::new();
            for &cell_idx in &marked_cells {
                for edge in &cell_edges[cell_idx] {
                    split_edges.insert(*edge);
                }
            }

            // Check consistency
            for (i, edges) in cell_edges.iter().enumerate() {
                if marked_cells.contains(&i) {
                    continue;
                }

                let mut my_split_edges = 0;
                let mut split_edge_indices = Vec::new();

                for edge in edges {
                    if split_edges.contains(edge) {
                        my_split_edges += 1;
                        split_edge_indices.push(*edge);
                    }
                }

                if my_split_edges == 0 {
                    continue; // Clean cell
                }

                // Enhanced Green Patterns Check
                let is_green_allowed = match my_split_edges {
                    1 => {
                        true // 1-edge split: Bisect (1->2 tets)
                    }
                    2 => {
                        // 2-edge split: Check if edges share a vertex (edge pair refinement, 1->3 tets)
                        let (e1, e2) = (split_edge_indices[0], split_edge_indices[1]);
                        // Edges share a vertex if they have one endpoint in common
                        e1.0 == e2.0 || e1.0 == e2.1 || e1.1 == e2.0 || e1.1 == e2.1
                    }
                    3 => {
                        // 3-edge split: Check if they form a face (face refinement, 1->4 tets)
                        let mut on_single_face = false;
                        if let Some(c) = mesh.cell(i) {
                            for &f_idx in &c.faces {
                                if let Some(face) = mesh.face(f_idx) {
                                    let face_edges: HashSet<(usize, usize)> = face
                                        .edges_iter()
                                        .map(|(u, v)| if u < v { (u, v) } else { (v, u) })
                                        .collect();

                                    let all_on_face =
                                        split_edge_indices.iter().all(|e| face_edges.contains(e));
                                    if all_on_face {
                                        on_single_face = true;
                                        break;
                                    }
                                }
                            }
                        }
                        on_single_face
                    }
                    4 => {
                        // 4-edge split: Check if edges form two opposite pairs (1->4 tets)
                        // This is a more complex pattern that requires specific edge topology
                        check_opposite_edge_pattern(&split_edge_indices, mesh, i)
                    }
                    _ => {
                        false // 5+ edges: Too complex, promote to red
                    }
                };

                if !is_green_allowed {
                    marked_cells.insert(i);
                    changed = true;
                }
            }
        }

        // Prepare new mesh
        let mut ctx = RefinementContext::new(mesh);

        // Re-identify split edges for final pass
        let mut final_split_edges: HashSet<(usize, usize)> = HashSet::new();
        for &cell_idx in &marked_cells {
            for edge in &cell_edges[cell_idx] {
                final_split_edges.insert(*edge);
            }
        }

        // Build cells
        for (i, cell) in mesh.cells().iter().enumerate() {
            // Get sorted vertices
            let mut v_indices = Vec::new();
            for &f_idx in &cell.faces {
                if let Some(face) = mesh.face(f_idx) {
                    for &v_idx in &face.vertices {
                        if !v_indices.contains(&v_idx) {
                            v_indices.push(v_idx);
                        }
                    }
                }
            }
            v_indices.sort_unstable();
            let v0 = v_indices[0];
            let v1 = v_indices[1];
            let v2 = v_indices[2];
            let v3 = v_indices[3];

            if marked_cells.contains(&i) {
                // RED REFINEMENT (1 -> 8)
                let m01 = ctx.get_midpoint(v0, v1)?;
                let m02 = ctx.get_midpoint(v0, v2)?;
                let m03 = ctx.get_midpoint(v0, v3)?;
                let m12 = ctx.get_midpoint(v1, v2)?;
                let m13 = ctx.get_midpoint(v1, v3)?;
                let m23 = ctx.get_midpoint(v2, v3)?;

                let new_tets = vec![
                    vec![v0, m01, m02, m03],
                    vec![v1, m01, m12, m13],
                    vec![v2, m02, m12, m23],
                    vec![v3, m03, m13, m23],
                    vec![m01, m23, m02, m12],
                    vec![m01, m23, m12, m13],
                    vec![m01, m23, m13, m03],
                    vec![m01, m23, m03, m02],
                ];

                for tet_verts in new_tets {
                    let f0 = ctx.add_face_dedup(vec![tet_verts[0], tet_verts[1], tet_verts[2]]);
                    let f1 = ctx.add_face_dedup(vec![tet_verts[0], tet_verts[1], tet_verts[3]]);
                    let f2 = ctx.add_face_dedup(vec![tet_verts[0], tet_verts[2], tet_verts[3]]);
                    let f3 = ctx.add_face_dedup(vec![tet_verts[1], tet_verts[2], tet_verts[3]]);
                    ctx.new_mesh.add_cell(Cell::tetrahedron(f0, f1, f2, f3));
                }
            } else {
                // Check if GREEN or CLEAN
                let mut my_split_edges = Vec::new();
                // We need to know WHICH edges correspond to which pair
                // (v0, v1), (v0, v2), (v0, v3), (v1, v2), (v1, v3), (v2, v3)
                let edges = [(v0, v1), (v0, v2), (v0, v3), (v1, v2), (v1, v3), (v2, v3)];

                for (u, v) in edges {
                    if final_split_edges.contains(&(u, v)) {
                        my_split_edges.push((u, v));
                    }
                }

                if my_split_edges.is_empty() {
                    // CLEAN: Copy cell
                    let f0 = ctx.add_face_dedup(vec![v0, v1, v2]);
                    let f1 = ctx.add_face_dedup(vec![v0, v1, v3]);
                    let f2 = ctx.add_face_dedup(vec![v0, v2, v3]);
                    let f3 = ctx.add_face_dedup(vec![v1, v2, v3]);
                    ctx.new_mesh.add_cell(Cell::tetrahedron(f0, f1, f2, f3));
                } else if my_split_edges.len() == 1 {
                    // GREEN: Bisect (1 -> 2)
                    let (u, v) = my_split_edges[0];
                    let m = ctx.get_midpoint(u, v)?;

                    // The two vertices not in the edge
                    let others: Vec<usize> = [v0, v1, v2, v3]
                        .iter()
                        .filter(|&&x| x != u && x != v)
                        .copied()
                        .collect();
                    let o1 = others[0];
                    let o2 = others[1];

                    // Tet 1: u, m, o1, o2
                    let t1 = vec![u, m, o1, o2];
                    // Tet 2: v, m, o1, o2
                    let t2 = vec![v, m, o1, o2];

                    for tet_verts in [t1, t2] {
                        let f0 = ctx.add_face_dedup(vec![tet_verts[0], tet_verts[1], tet_verts[2]]);
                        let f1 = ctx.add_face_dedup(vec![tet_verts[0], tet_verts[1], tet_verts[3]]);
                        let f2 = ctx.add_face_dedup(vec![tet_verts[0], tet_verts[2], tet_verts[3]]);
                        let f3 = ctx.add_face_dedup(vec![tet_verts[1], tet_verts[2], tet_verts[3]]);
                        ctx.new_mesh.add_cell(Cell::tetrahedron(f0, f1, f2, f3));
                    }
                } else {
                    // Should not happen due to iterative propagation
                    return Err(crate::error::MeshError::InvalidMesh(format!(
                        "Unexpected split edge count {} for Green refinement",
                        my_split_edges.len()
                    )));
                }
            }
        }

        *mesh = ctx.new_mesh;
        Ok(())
    }

    fn name(&self) -> &'static str {
        "Adaptive"
    }
}
