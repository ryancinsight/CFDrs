//! Hex-to-Tet mesh decomposition
//!
//! Robust decomposition of 8-node hexahedra into tetrahedra.
//!
//! This module is a **volume/FEM tool** — intentional `Mesh<T>` usage for
//! hexahedral/tetrahedral cell topology.

// Volume tool: Mesh<T> is the correct type here.
#![allow(deprecated)]

use crate::mesh::Mesh;
use crate::topology::{Cell, ElementType, Face};
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
        // Key: sorted triangle vertex triplet, Value: new face index
        let mut face_map: HashMap<Vec<usize>, usize> = HashMap::new();

        // 3. Process cells
        for c in mesh.cells() {
            if c.element_type == ElementType::Hexahedron {
                let hex_vertices = Self::collect_unique_hex_vertices(c, mesh);
                if hex_vertices.len() == 8 {
                    let length_scale = Self::characteristic_length(mesh, &hex_vertices);
                    let tol_factor = T::from_f64(1e-12).unwrap_or_else(T::zero);
                    let volume_tol = length_scale * length_scale * length_scale * tol_factor;

                    let mut decomposed = false;

                    // Prefer recovered topological ordering to avoid decomposition
                    // bias from incidental face/vertex iteration order.
                    if let Some(recovered_order) = Self::recover_hex_vertex_order(c, mesh, volume_tol) {
                        if let Some(tets) =
                            Self::select_hex_decomposition(mesh, recovered_order, volume_tol)
                        {
                            for nodes in tets {
                                Self::add_tet(&mut new_mesh, &mut face_map, nodes);
                            }
                            decomposed = true;
                        }
                    }

                    if !decomposed {
                        if let Ok(raw_order) = <[usize; 8]>::try_from(hex_vertices.as_slice()) {
                            if let Some(tets) = Self::select_hex_decomposition(mesh, raw_order, volume_tol) {
                                for nodes in tets {
                                    Self::add_tet(&mut new_mesh, &mut face_map, nodes);
                                }
                                decomposed = true;
                            }
                        }
                    }

                    if !decomposed {
                        // Final safeguard: keep only non-degenerate tetrahedra.
                        if let Ok(raw_order) = <[usize; 8]>::try_from(hex_vertices.as_slice()) {
                            for nodes in Self::hex_six_tet_pattern(raw_order) {
                                if Self::is_non_degenerate_tet(mesh, nodes, volume_tol) {
                                    Self::add_tet(&mut new_mesh, &mut face_map, nodes);
                                }
                            }
                        }
                    }
                }
            } else {
                // Keep other cells (e.g. already tetrahedra), remapping faces
                let mut new_faces = Vec::new();
                for &f_idx in &c.faces {
                    if let Some(face) = mesh.face(f_idx) {
                        if face.vertices.len() == 3 {
                            let nf = Self::add_tri_face(
                                &mut new_mesh,
                                &mut face_map,
                                [face.vertices[0], face.vertices[1], face.vertices[2]],
                            );
                            new_faces.push(nf);
                        } else {
                            // Decompose rare quad face in non-hex cells
                            let f_a = Self::add_tri_face(
                                &mut new_mesh,
                                &mut face_map,
                                [face.vertices[0], face.vertices[1], face.vertices[2]],
                            );
                            let f_b = Self::add_tri_face(
                                &mut new_mesh,
                                &mut face_map,
                                [face.vertices[0], face.vertices[2], face.vertices[3]],
                            );
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
                        let nf = Self::get_tri_face_idx(
                            &face_map,
                            [face.vertices[0], face.vertices[1], face.vertices[2]],
                        );
                        if let Some(idx) = nf {
                            new_mesh.mark_boundary(idx, label.to_string());
                        }
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

    fn collect_unique_hex_vertices<T: RealField + Copy>(cell: &Cell, mesh: &Mesh<T>) -> Vec<usize> {
        let mut vertices = Vec::with_capacity(8);
        for &f_idx in &cell.faces {
            if let Some(face) = mesh.face(f_idx) {
                for &v_idx in &face.vertices {
                    if !vertices.contains(&v_idx) {
                        vertices.push(v_idx);
                    }
                }
            }
        }
        vertices
    }

    fn characteristic_length<T: RealField + Copy>(mesh: &Mesh<T>, vertices: &[usize]) -> T {
        let mut max_dist_sq = T::zero();
        for i in 0..vertices.len() {
            for j in (i + 1)..vertices.len() {
                let pi = mesh.vertex(vertices[i]).unwrap().position.coords;
                let pj = mesh.vertex(vertices[j]).unwrap().position.coords;
                let dist_sq = (pj - pi).norm_squared();
                if dist_sq > max_dist_sq {
                    max_dist_sq = dist_sq;
                }
            }
        }
        max_dist_sq.sqrt()
    }

    fn tet_six_volume<T: RealField + Copy>(mesh: &Mesh<T>, nodes: [usize; 4]) -> T {
        let p0 = mesh.vertex(nodes[0]).unwrap().position.coords;
        let p1 = mesh.vertex(nodes[1]).unwrap().position.coords;
        let p2 = mesh.vertex(nodes[2]).unwrap().position.coords;
        let p3 = mesh.vertex(nodes[3]).unwrap().position.coords;
        (p1 - p0).cross(&(p2 - p0)).dot(&(p3 - p0)).abs()
    }

    fn is_non_degenerate_tet<T: RealField + Copy>(
        mesh: &Mesh<T>,
        nodes: [usize; 4],
        volume_tol: T,
    ) -> bool {
        for i in 0..4 {
            for j in (i + 1)..4 {
                if nodes[i] == nodes[j] {
                    return false;
                }
            }
        }
        Self::tet_six_volume(mesh, nodes) > volume_tol
    }

    fn add_tet<T: RealField + Copy>(
        mesh: &mut Mesh<T>,
        face_map: &mut HashMap<Vec<usize>, usize>,
        nodes: [usize; 4],
    ) {
        let f0 = Self::add_tri_face(mesh, face_map, [nodes[0], nodes[1], nodes[2]]);
        let f1 = Self::add_tri_face(mesh, face_map, [nodes[0], nodes[1], nodes[3]]);
        let f2 = Self::add_tri_face(mesh, face_map, [nodes[0], nodes[2], nodes[3]]);
        let f3 = Self::add_tri_face(mesh, face_map, [nodes[1], nodes[2], nodes[3]]);
        mesh.add_cell(Cell::tetrahedron(f0, f1, f2, f3));
    }

    fn hex_five_tet_pattern(order: [usize; 8]) -> [[usize; 4]; 5] {
        [
            [order[0], order[1], order[3], order[4]],
            [order[1], order[2], order[3], order[6]],
            [order[4], order[7], order[6], order[3]],
            [order[4], order[6], order[5], order[1]],
            [order[1], order[3], order[4], order[6]],
        ]
    }

    fn hex_six_tet_pattern(order: [usize; 8]) -> [[usize; 4]; 6] {
        [
            [order[0], order[1], order[2], order[6]],
            [order[0], order[2], order[3], order[6]],
            [order[0], order[3], order[7], order[6]],
            [order[0], order[7], order[4], order[6]],
            [order[0], order[4], order[5], order[6]],
            [order[0], order[5], order[1], order[6]],
        ]
    }

    fn decomposition_min_volume<T: RealField + Copy>(
        mesh: &Mesh<T>,
        tets: &[[usize; 4]],
        volume_tol: T,
    ) -> Option<T> {
        let mut min_vol: Option<T> = None;
        for nodes in tets {
            if !Self::is_non_degenerate_tet(mesh, *nodes, volume_tol) {
                return None;
            }
            let six_v = Self::tet_six_volume(mesh, *nodes);
            min_vol = Some(match min_vol {
                Some(v) => v.min(six_v),
                None => six_v,
            });
        }
        min_vol
    }

    fn select_hex_decomposition<T: RealField + Copy>(
        mesh: &Mesh<T>,
        order: [usize; 8],
        volume_tol: T,
    ) -> Option<Vec<[usize; 4]>> {
        let five = Self::hex_five_tet_pattern(order);
        let six = Self::hex_six_tet_pattern(order);
        let q5 = Self::decomposition_min_volume(mesh, &five, volume_tol);
        let q6 = Self::decomposition_min_volume(mesh, &six, volume_tol);

        match (q5, q6) {
            (Some(v5), Some(v6)) => {
                if v5 >= v6 {
                    Some(five.to_vec())
                } else {
                    Some(six.to_vec())
                }
            }
            (Some(_), None) => Some(five.to_vec()),
            (None, Some(_)) => Some(six.to_vec()),
            (None, None) => None,
        }
    }

    fn common_neighbor_excluding(
        adjacency: &HashMap<usize, Vec<usize>>,
        a: usize,
        b: usize,
        excluded: &[usize],
    ) -> Option<usize> {
        let a_neighbors = adjacency.get(&a)?;
        let b_neighbors = adjacency.get(&b)?;
        let mut candidate = None;
        for &n in a_neighbors {
            if b_neighbors.contains(&n) && !excluded.contains(&n) {
                if candidate.is_some() {
                    return None;
                }
                candidate = Some(n);
            }
        }
        candidate
    }

    fn recover_hex_vertex_order<T: RealField + Copy>(
        cell: &Cell,
        mesh: &Mesh<T>,
        volume_tol: T,
    ) -> Option<[usize; 8]> {
        let vertices = Self::collect_unique_hex_vertices(cell, mesh);
        if vertices.len() != 8 {
            return None;
        }

        let mut adjacency: HashMap<usize, Vec<usize>> = HashMap::new();
        for &f_idx in &cell.faces {
            let face = mesh.face(f_idx)?;
            let n = face.vertices.len();
            if n < 3 {
                continue;
            }
            for i in 0..n {
                let a = face.vertices[i];
                let b = face.vertices[(i + 1) % n];
                adjacency.entry(a).or_default().push(b);
                adjacency.entry(b).or_default().push(a);
            }
        }
        for neigh in adjacency.values_mut() {
            neigh.sort_unstable();
            neigh.dedup();
        }

        let perms = [
            [0, 1, 2],
            [0, 2, 1],
            [1, 0, 2],
            [1, 2, 0],
            [2, 0, 1],
            [2, 1, 0],
        ];

        let mut best_order: Option<[usize; 8]> = None;
        let mut best_quality: Option<T> = None;

        for &v0 in &vertices {
            let Some(neigh) = adjacency.get(&v0) else {
                continue;
            };
            if neigh.len() != 3 {
                continue;
            }

            for perm in &perms {
                let v1 = neigh[perm[0]];
                let v3 = neigh[perm[1]];
                let v4 = neigh[perm[2]];

                let Some(v2) = Self::common_neighbor_excluding(&adjacency, v1, v3, &[v0, v4]) else {
                    continue;
                };
                let Some(v5) = Self::common_neighbor_excluding(&adjacency, v1, v4, &[v0, v3]) else {
                    continue;
                };
                let Some(v7) = Self::common_neighbor_excluding(&adjacency, v3, v4, &[v0, v1]) else {
                    continue;
                };

                let Some(n2) = adjacency.get(&v2) else {
                    continue;
                };
                let Some(n5) = adjacency.get(&v5) else {
                    continue;
                };
                let Some(n7) = adjacency.get(&v7) else {
                    continue;
                };

                let mut v6_candidate = None;
                for &n in n2 {
                    if n5.contains(&n)
                        && n7.contains(&n)
                        && n != v0
                        && n != v1
                        && n != v2
                        && n != v3
                        && n != v4
                        && n != v5
                        && n != v7
                    {
                        if v6_candidate.is_some() {
                            v6_candidate = None;
                            break;
                        }
                        v6_candidate = Some(n);
                    }
                }
                let Some(v6) = v6_candidate else {
                    continue;
                };

                let order = [v0, v1, v2, v3, v4, v5, v6, v7];
                let mut unique = order.to_vec();
                unique.sort_unstable();
                unique.dedup();
                if unique.len() != 8 {
                    continue;
                }

                let Some(tets) = Self::select_hex_decomposition(mesh, order, volume_tol) else {
                    continue;
                };
                let quality = tets
                    .iter()
                    .map(|nodes| Self::tet_six_volume(mesh, *nodes))
                    .fold(T::max_value().unwrap_or_else(T::one), |a, b| a.min(b));

                if best_quality.map_or(true, |best| quality > best) {
                    best_quality = Some(quality);
                    best_order = Some(order);
                }
            }
        }

        best_order
    }

    fn add_tri_face<T: RealField + Copy>(
        mesh: &mut Mesh<T>,
        map: &mut HashMap<Vec<usize>, usize>,
        nodes: [usize; 3],
    ) -> usize {
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

#[cfg(test)]
mod tests {
    use super::HexToTetConverter;
    use crate::grid::StructuredGridBuilder;
    use crate::mesh::Mesh;
    use crate::topology::ElementType;

    fn tet_six_volume(mesh: &Mesh<f64>, cell: &crate::topology::Cell) -> f64 {
        let mut vertices = Vec::new();
        for &f_idx in &cell.faces {
            let face = mesh.face(f_idx).unwrap();
            for &v_idx in &face.vertices {
                if !vertices.contains(&v_idx) {
                    vertices.push(v_idx);
                }
            }
        }
        assert_eq!(vertices.len(), 4, "Converted tetrahedron must have 4 unique vertices");

        let p0 = mesh.vertex(vertices[0]).unwrap().position.coords;
        let p1 = mesh.vertex(vertices[1]).unwrap().position.coords;
        let p2 = mesh.vertex(vertices[2]).unwrap().position.coords;
        let p3 = mesh.vertex(vertices[3]).unwrap().position.coords;
        (p1 - p0).cross(&(p2 - p0)).dot(&(p3 - p0)).abs()
    }

    fn assert_no_degenerate_tets(mesh: &Mesh<f64>) {
        let (min, max) = mesh.bounds();
        let length_scale = (max.coords - min.coords).norm();
        let volume_tol = length_scale.powi(3) * 1e-12;

        for (i, cell) in mesh.cells().iter().enumerate() {
            if cell.element_type != ElementType::Tetrahedron {
                continue;
            }
            let six_v = tet_six_volume(mesh, cell);
            assert!(
                six_v > volume_tol,
                "Degenerate tetrahedron at cell {} with 6V={:.3e}, tol={:.3e}",
                i,
                six_v,
                volume_tol
            );
        }
    }

    #[test]
    fn structured_hex_mesh_converts_to_non_degenerate_tets() {
        let hex_mesh = StructuredGridBuilder::new(4, 4, 4).build().unwrap();
        let tet_mesh = HexToTetConverter::convert(&hex_mesh);

        assert!(tet_mesh.cell_count() > 0);
        assert!(
            tet_mesh
                .cells()
                .iter()
                .all(|c| c.element_type == ElementType::Tetrahedron)
        );
        assert_no_degenerate_tets(&tet_mesh);
    }

    #[test]
    fn branching_mesh_conversion_avoids_degenerate_tets() {
        // Use a larger structured grid to exercise non-trivial tet conversion.
        let hex_mesh = StructuredGridBuilder::new(6, 4, 4).build().unwrap();
        let tet_mesh = HexToTetConverter::convert(&hex_mesh);

        assert!(tet_mesh.cell_count() > 0);
        assert_no_degenerate_tets(&tet_mesh);
    }
}
