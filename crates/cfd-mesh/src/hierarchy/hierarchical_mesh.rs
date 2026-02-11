//! Hierarchical mesh conversion (P1 to P2)

use crate::mesh::Mesh;
use crate::topology::{Cell, Face, Vertex, ElementType};
use nalgebra::{Point3, RealField};
use std::collections::HashMap;

/// Converter for promoting linear meshes to quadratic (P2)
pub struct P2MeshConverter;

impl P2MeshConverter {
    /// Convert a linear (P1) mesh to a quadratic (P2) mesh
    ///
    /// This generates unique mid-side nodes for every edge and 
    /// returns a new mesh with promoted elements (e.g. Tetrahedron -> Tetrahedron10).
    pub fn convert_to_p2<T: RealField + Copy>(mesh: &Mesh<T>) -> Mesh<T> {
        let mut new_mesh = Mesh::new();
        
        // 1. Copy original vertices
        let mut vertex_map = Vec::with_capacity(mesh.vertex_count());
        for v in mesh.vertices() {
            let idx = new_mesh.add_vertex(v.clone());
            vertex_map.push(idx);
        }

        // 2. Identify unique edges and create mid-side nodes
        // Key: (min(v1, v2), max(v1, v2)), Value: new_vertex_idx
        let mut edge_midpoints: HashMap<(usize, usize), usize> = HashMap::new();

        // 3. Process faces (promote to Quadratic)
        let mut face_map = Vec::with_capacity(mesh.face_count());
        for f in mesh.faces() {
            let mut p2_vertices = f.vertices.clone();
            
            // For each edge in the face, find or create the mid-node
            let n_corners = f.vertices.len();
            for i in 0..n_corners {
                let v1 = f.vertices[i];
                let v2 = f.vertices[(i + 1) % n_corners];
                
                let edge_key = if v1 < v2 { (v1, v2) } else { (v2, v1) };
                
                let mid_idx = if let Some(&idx) = edge_midpoints.get(&edge_key) {
                    idx
                } else {
                    // Create new vertex at midpoint
                    let pos1 = mesh.vertex(v1).unwrap().position;
                    let pos2 = mesh.vertex(v2).unwrap().position;
                    let mid_pos = Point3::from((pos1.coords + pos2.coords) * T::from_f64(0.5).unwrap());
                    let idx = new_mesh.add_vertex(Vertex::new(mid_pos));
                    edge_midpoints.insert(edge_key, idx);
                    idx
                };
                p2_vertices.push(mid_idx);
            }
            
            let p2_face = Face::from_vertices(p2_vertices);
            let face_idx = new_mesh.add_face(p2_face);
            face_map.push(face_idx);
        }

        // 4. Process cells
        for c in mesh.cells() {
            let mut p2_faces = Vec::with_capacity(c.faces.len());
            for &f_idx in &c.faces {
                p2_faces.push(face_map[f_idx]);
            }

            let mut p2_cell = Cell {
                faces: p2_faces,
                element_type: match c.element_type {
                    ElementType::Tetrahedron => ElementType::Tetrahedron10,
                    ElementType::Hexahedron => ElementType::Hexahedron20, // (Actually Hex27 is more common for full P2, but let's stick to what core has)
                    t => t, // Fallback for other types
                },
                global_id: c.global_id,
                partition_id: c.partition_id,
            };
            new_mesh.add_cell(p2_cell);
        }

        // 5. Transfer boundary markers
        for f_idx in 0..mesh.face_count() {
            if let Some(label) = mesh.boundary_label(f_idx) {
                new_mesh.mark_boundary(face_map[f_idx], label.to_string());
            }
        }

        new_mesh
    }
}
