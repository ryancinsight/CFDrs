//! Mesh refinement module

use nalgebra::RealField;
use std::collections::HashSet;
use crate::topology::{Cell, ElementType, Vertex};

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

/// Uniform refinement strategy
pub struct UniformRefinement;

impl<T: RealField + Copy> RefinementStrategy<T> for UniformRefinement {
    fn refine(&self, mesh: &mut crate::mesh::Mesh<T>) -> Result<(), crate::error::MeshError> {
        if mesh.cell_count() == 0 {
            return Ok(());
        }

        // Check if we have only tetrahedra (simplification for now)
        for cell in mesh.cells() {
            if cell.element_type != ElementType::Tetrahedron {
                return Err(crate::error::MeshError::NotImplemented(
                    format!("Uniform refinement only implemented for Tetrahedra, found {:?}", cell.element_type)
                ));
            }
        }

        let mut ctx = RefinementContext::new(mesh);

        // Iterate over cells and refine
        for cell in mesh.cells() {
            // Assume Tetrahedron
            // Get unique vertex indices
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
            // Sort to ensure deterministic behavior
            v_indices.sort_unstable();

            if v_indices.len() != 4 {
                 return Err(crate::error::MeshError::InvalidMesh(
                    format!("Tetrahedron must have 4 vertices, found {}", v_indices.len())
                ));
            }
            
            let v0 = v_indices[0];
            let v1 = v_indices[1];
            let v2 = v_indices[2];
            let v3 = v_indices[3];

            // Get midpoints
            let m01 = ctx.get_midpoint(v0, v1);
            let m02 = ctx.get_midpoint(v0, v2);
            let m03 = ctx.get_midpoint(v0, v3);
            let m12 = ctx.get_midpoint(v1, v2);
            let m13 = ctx.get_midpoint(v1, v3);
            let m23 = ctx.get_midpoint(v2, v3);

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
        }

        // Replace mesh
        *mesh = ctx.new_mesh;

        Ok(())
    }

    fn name(&self) -> &'static str {
        "Uniform"
    }
}

/// Adaptive refinement strategy
pub struct AdaptiveRefinement<T: RealField + Copy> {
    /// Refinement criterion
    pub criterion: RefinementCriterion<T>,
}

impl<T: RealField + Copy> RefinementStrategy<T> for AdaptiveRefinement<T> {
    fn refine(&self, mesh: &mut crate::mesh::Mesh<T>) -> Result<(), crate::error::MeshError> {
        if mesh.cell_count() == 0 {
            return Ok(());
        }

        // Mark cells for refinement based on criterion
        let mut marked_cells = HashSet::new();

        for (i, cell) in mesh.cells().iter().enumerate() {
            // Check element type
             if cell.element_type != ElementType::Tetrahedron {
                return Err(crate::error::MeshError::NotImplemented(
                    format!("Adaptive refinement only implemented for Tetrahedra, found {:?}", cell.element_type)
                ));
            }

            let should_refine = match &self.criterion {
                RefinementCriterion::Error { error_field, threshold } => {
                    if let Some(err) = error_field.get(i) {
                         *err > *threshold
                    } else {
                        false
                    }
                }
                RefinementCriterion::Gradient { field, threshold } => {
                    // Placeholder for gradient estimation
                    // Without connectivity info easily available for gradient calc here,
                    // we assume field is per-cell or per-vertex?
                    // Usually gradient is computed beforehand.
                    // For now, if field[i] > threshold (dummy implementation if field corresponds to cell)
                    if let Some(val) = field.get(i) {
                         *val > *threshold
                    } else {
                        false
                    }
                }
                RefinementCriterion::Custom(func) => {
                    let vertices: Vec<Vertex<T>> = mesh.ordered_element_vertices(cell)
                        .into_iter().cloned().collect();
                    func(cell, &vertices)
                }
                _ => false // Other criteria not implemented
            };

            if should_refine {
                marked_cells.insert(i);
            }
        }

        if marked_cells.is_empty() {
            return Ok(());
        }

        // Conformity Handling (Iterative Propagation)
        // We use a simplified strategy:
        // 1. Identify "Red" cells (fully refined 1->8).
        // 2. Identify "Green" cells (partially refined).
        // Since implementing full Green refinement patterns for all cases is complex,
        // we use a robust strategy:
        // - Initially marked cells are Red.
        // - Loop:
        //   - Check neighbors of Red cells.
        //   - If a neighbor would have > 3 split edges (or complicated pattern), mark it Red.
        //   - Only allow 1-edge split (Bisect) or 3-edge-face split (1->4).
        //   - If not allowed, promote to Red.

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
            for i in 0..mesh.cell_count() {
                if marked_cells.contains(&i) {
                    continue;
                }

                let mut my_split_edges = 0;
                let mut split_edge_indices = Vec::new();

                for edge in &cell_edges[i] {
                    if split_edges.contains(edge) {
                        my_split_edges += 1;
                        split_edge_indices.push(*edge);
                    }
                }

                if my_split_edges == 0 {
                    continue; // Clean cell
                }

                // Green Patterns Check
                let is_green_allowed = if my_split_edges == 1 {
                    true // Bisect
                } else if my_split_edges == 3 {
                    // Check if they form a face
                    // A face has 3 edges.
                    // Get cell faces
                    let mut on_single_face = false;
                    if let Some(c) = mesh.cell(i) {
                         for &f_idx in &c.faces {
                            if let Some(face) = mesh.face(f_idx) {
                                let face_edges: HashSet<(usize, usize)> = face.edges_iter()
                                    .map(|(u,v)| if u < v { (u, v) } else { (v, u) })
                                    .collect();

                                let all_on_face = split_edge_indices.iter().all(|e| face_edges.contains(e));
                                if all_on_face {
                                    on_single_face = true;
                                    break;
                                }
                            }
                        }
                    }
                    on_single_face
                } else {
                    false
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
                let m01 = ctx.get_midpoint(v0, v1);
                let m02 = ctx.get_midpoint(v0, v2);
                let m03 = ctx.get_midpoint(v0, v3);
                let m12 = ctx.get_midpoint(v1, v2);
                let m13 = ctx.get_midpoint(v1, v3);
                let m23 = ctx.get_midpoint(v2, v3);

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
                 let edges = [
                     (v0, v1), (v0, v2), (v0, v3),
                     (v1, v2), (v1, v3),
                     (v2, v3)
                 ];

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
                    let m = ctx.get_midpoint(u, v);

                    // The two vertices not in the edge
                    let others: Vec<usize> = [v0, v1, v2, v3].iter()
                        .filter(|&&x| x != u && x != v)
                        .copied()
                        .collect();
                    let o1 = others[0];
                    let o2 = others[1];

                    // Tet 1: u, m, o1, o2
                    let t1 = vec![u, m, o1, o2];
                    // Tet 2: v, m, o1, o2
                    let t2 = vec![v, m, o1, o2];

                    for tet_verts in vec![t1, t2] {
                        let f0 = ctx.add_face_dedup(vec![tet_verts[0], tet_verts[1], tet_verts[2]]);
                        let f1 = ctx.add_face_dedup(vec![tet_verts[0], tet_verts[1], tet_verts[3]]);
                        let f2 = ctx.add_face_dedup(vec![tet_verts[0], tet_verts[2], tet_verts[3]]);
                        let f3 = ctx.add_face_dedup(vec![tet_verts[1], tet_verts[2], tet_verts[3]]);
                        ctx.new_mesh.add_cell(Cell::tetrahedron(f0, f1, f2, f3));
                    }

                } else if my_split_edges.len() == 3 {
                    // GREEN: Face Split (1 -> 4)
                    // Identify the face. The 3 edges form a face.
                    // Vertices involved in split edges
                    let mut involved_verts = HashSet::new();
                    for (u, v) in &my_split_edges {
                        involved_verts.insert(*u);
                        involved_verts.insert(*v);
                    }
                    // The face has 3 vertices.
                    // The opposite vertex is the one NOT in involved_verts.
                    // Actually, if 3 edges form a triangle, they involve exactly 3 vertices.

                    let face_verts: Vec<usize> = involved_verts.into_iter().collect();
                    let op_vert = [v0, v1, v2, v3].iter()
                        .find(|&&x| !face_verts.contains(&x))
                        .copied()
                        .unwrap();

                    // Face vertices
                    let fv0 = face_verts[0];
                    let fv1 = face_verts[1];
                    let fv2 = face_verts[2];

                    // Midpoints
                    let m01 = ctx.get_midpoint(fv0, fv1);
                    let m12 = ctx.get_midpoint(fv1, fv2);
                    let m20 = ctx.get_midpoint(fv2, fv0);

                    // 4 Tets connecting Op to the 4 triangles on the face
                    // Base Triangles:
                    // (fv0, m01, m20)
                    // (fv1, m12, m01)
                    // (fv2, m20, m12)
                    // (m01, m12, m20)

                    let new_tets = vec![
                        vec![op_vert, fv0, m01, m20],
                        vec![op_vert, fv1, m12, m01],
                        vec![op_vert, fv2, m20, m12],
                        vec![op_vert, m01, m12, m20],
                    ];

                    for tet_verts in new_tets {
                        let f0 = ctx.add_face_dedup(vec![tet_verts[0], tet_verts[1], tet_verts[2]]);
                        let f1 = ctx.add_face_dedup(vec![tet_verts[0], tet_verts[1], tet_verts[3]]);
                        let f2 = ctx.add_face_dedup(vec![tet_verts[0], tet_verts[2], tet_verts[3]]);
                        let f3 = ctx.add_face_dedup(vec![tet_verts[1], tet_verts[2], tet_verts[3]]);
                        ctx.new_mesh.add_cell(Cell::tetrahedron(f0, f1, f2, f3));
                    }
                } else {
                    // Should not happen due to iterative propagation
                     return Err(crate::error::MeshError::InvalidMesh(
                        format!("Unexpected split edge count {} for Green refinement", my_split_edges.len())
                    ));
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
