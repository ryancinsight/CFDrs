//! Hex-to-Tet mesh decomposition
//!
//! Standard decomposition of 8-node hexahedra into 5 or 6 tetrahedra.

use crate::mesh::Mesh;
use crate::topology::{Cell, Face, ElementType};
use nalgebra::RealField;
use std::collections::HashMap;

/// Converter for decomposing hexahedral meshes into tetrahedral ones
pub struct HexToTetConverter;

impl HexToTetConverter {
    /// Decompose all hexahedral cells in a mesh into tetrahedra
    pub fn convert<T: RealField + Copy>(mesh: &Mesh<T>) -> Mesh<T> {
        let mut new_mesh = Mesh::new();
        
        // 1. Copy all vertices
        for v in mesh.vertices() {
            new_mesh.add_vertex(v.clone());
        }

        // 2. Identify boundary faces and map them
        // We'll need to decompose quad faces into triangular faces
        // Key: [v1, v2, v3] sorted, Value: new_face_idx
        let mut face_map: HashMap<Vec<usize>, usize> = HashMap::new();

        // 3. Process cells
        for c in mesh.cells() {
            if c.element_type == ElementType::Hexahedron {
                // Extract 8 vertices of the hex
                // We use the same canonical ordering as StructuredGridBuilder
                let mut v = Vec::new();
                for &f_idx in &c.faces {
                    if let Some(face) = mesh.face(f_idx) {
                        for &v_idx in &face.vertices {
                            if !v.contains(&v_idx) {
                                v.push(v_idx);
                            }
                        }
                    }
                }
                
                if v.len() == 8 {
                    // Standard 5-tetrahedron decomposition
                    // v0..v7 ordering:
                    // Bottom: v0, v1, v2, v3
                    // Top: v4, v5, v6, v7
                    // Correct connectivity from StructuredGridBuilder:
                    // v0=(0,0,0), v1=(1,0,0), v2=(1,1,0), v3=(0,1,0)
                    // v4=(0,0,1), v5=(1,0,1), v6=(1,1,1), v7=(0,1,1)
                    
                    let tets = vec![
                        [v[0], v[1], v[3], v[4]],
                        [v[1], v[2], v[3], v[6]],
                        [v[4], v[7], v[6], v[3]],
                        [v[4], v[6], v[5], v[1]],
                        [v[1], v[3], v[4], v[6]],
                    ];

                    for nodes in tets {
                        let f0 = Self::add_tri_face(&mut new_mesh, &mut face_map, [nodes[0], nodes[1], nodes[2]]);
                        let f1 = Self::add_tri_face(&mut new_mesh, &mut face_map, [nodes[0], nodes[1], nodes[3]]);
                        let f2 = Self::add_tri_face(&mut new_mesh, &mut face_map, [nodes[0], nodes[2], nodes[3]]);
                        let f3 = Self::add_tri_face(&mut new_mesh, &mut face_map, [nodes[1], nodes[2], nodes[3]]);
                        new_mesh.add_cell(Cell::tetrahedron(f0, f1, f2, f3));
                    }
                }
            } else {
                // Keep other cells (e.g. already tetrahedra)
                // But we must remap their faces
                let mut new_faces = Vec::new();
                for &f_idx in &c.faces {
                    if let Some(face) = mesh.face(f_idx) {
                        if face.vertices.len() == 3 {
                            let nf = Self::add_tri_face(&mut new_mesh, &mut face_map, [face.vertices[0], face.vertices[1], face.vertices[2]]);
                            new_faces.push(nf);
                        } else {
                             // Decompose quad face if found in non-hex cell (rare)
                             let f_a = Self::add_tri_face(&mut new_mesh, &mut face_map, [face.vertices[0], face.vertices[1], face.vertices[2]]);
                             let f_b = Self::add_tri_face(&mut new_mesh, &mut face_map, [face.vertices[0], face.vertices[2], face.vertices[3]]);
                             new_faces.push(f_a);
                             new_faces.push(f_b);
                        }
                    }
                }
                let mut new_cell = c.clone();
                new_cell.faces = new_faces;
                new_mesh.add_cell(new_cell);
            }
        }

        // 4. Transfer and decompose boundary markers
        for f_idx in mesh.marked_boundary_faces() {
            if let Some(label) = mesh.boundary_label(f_idx) {
                if let Some(face) = mesh.face(f_idx) {
                    if face.vertices.len() == 3 {
                         let nf = Self::get_tri_face_idx(&face_map, [face.vertices[0], face.vertices[1], face.vertices[2]]);
                         if let Some(idx) = nf { new_mesh.mark_boundary(idx, label.to_string()); }
                    } else if face.vertices.len() == 4 {
                         let v = &face.vertices;
                         let possible_tris = [
                             [v[0], v[1], v[2]],
                             [v[1], v[2], v[3]],
                             [v[2], v[3], v[0]],
                             [v[3], v[0], v[1]],
                         ];
                         for tri_nodes in &possible_tris {
                             if let Some(idx) = Self::get_tri_face_idx(&face_map, *tri_nodes) {
                                 new_mesh.mark_boundary(idx, label.to_string());
                             }
                         }
                    }
                }
            }
        }

        new_mesh
    }

    fn add_tri_face<T: RealField + Copy>(mesh: &mut Mesh<T>, map: &mut HashMap<Vec<usize>, usize>, nodes: [usize; 3]) -> usize {
        let mut key = nodes.to_vec();
        key.sort_unstable();
        if let Some(&idx) = map.get(&key) {
            idx
        } else {
            let idx = mesh.add_face(Face::triangle(nodes[0], nodes[1], nodes[2]));
            map.insert(key, idx);
            idx
        }
    }
    
    fn get_tri_face_idx(map: &HashMap<Vec<usize>, usize>, nodes: [usize; 3]) -> Option<usize> {
        let mut key = nodes.to_vec();
        key.sort_unstable();
        map.get(&key).copied()
    }
}
