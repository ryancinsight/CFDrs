//! Mesh connectivity and topological relationships
//!
//! This module provides structures for traversing the mesh topology,
//! such as finding neighbors, edges, and boundary faces.

use super::Mesh;
use nalgebra::RealField;
use std::collections::{HashMap, HashSet};

/// Connectivity information for mesh topology
#[derive(Debug, Clone)]
pub struct Connectivity {
    /// Node to elements mapping
    pub node_to_elements: HashMap<usize, Vec<usize>>,
    /// Element neighbors (elements sharing at least one node)
    pub element_neighbors: HashMap<usize, Vec<usize>>,
    /// Boundary faces
    pub boundary_faces: Vec<Face>,
}

/// Face representation
#[derive(Debug, Clone)]
pub struct Face {
    /// Node indices defining the face
    pub nodes: Vec<usize>,
    /// Parent element index
    pub element: usize,
    /// Face local index in element
    pub local_index: usize,
}

/// Edge connectivity for 2D/3D meshes
#[derive(Debug, Clone)]
pub struct EdgeConnectivity {
    /// Edges in the mesh
    pub edges: Vec<Edge>,
    /// Node to edges mapping
    pub node_to_edges: HashMap<usize, Vec<usize>>,
}

/// Edge representation
#[derive(Debug, Clone)]
pub struct Edge {
    /// Start node index
    pub start: usize,
    /// End node index
    pub end: usize,
    /// Adjacent elements
    pub elements: Vec<usize>,
}

impl<T: RealField + Copy> Mesh<T> {
    /// Build connectivity information
    #[must_use]
    pub fn build_connectivity(&self) -> Connectivity {
        let mut node_to_elements: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut element_neighbors: HashMap<usize, Vec<usize>> = HashMap::new();

        // Build node to elements mapping
        for (elem_idx, element) in self.elements.iter().enumerate() {
            for &node_idx in &element.nodes {
                node_to_elements.entry(node_idx).or_default().push(elem_idx);
            }
        }

        // Find element neighbors through shared nodes
        for (elem_idx, element) in self.elements.iter().enumerate() {
            let mut neighbor_set = HashSet::new();
            for &node_idx in &element.nodes {
                if let Some(connected_elements) = node_to_elements.get(&node_idx) {
                    for &other_elem in connected_elements {
                        if other_elem != elem_idx {
                            neighbor_set.insert(other_elem);
                        }
                    }
                }
            }
            element_neighbors.insert(elem_idx, neighbor_set.into_iter().collect());
        }

        Connectivity {
            node_to_elements,
            element_neighbors,
            boundary_faces: Vec::new(),
        }
    }

    /// Build edge connectivity
    pub fn build_edge_connectivity(&self) -> EdgeConnectivity {
        let mut edges_list: Vec<(usize, usize)> = Vec::new();
        let mut node_to_edges = HashMap::new();

        for element in &self.elements {
            let n = element.nodes.len();
            match element.element_type {
                crate::geometry::mesh::ElementType::Line | crate::geometry::mesh::ElementType::Line3 => {
                    if n >= 2 {
                        edges_list.push(order_edge(element.nodes[0], element.nodes[1]));
                    }
                }
                crate::geometry::mesh::ElementType::Triangle
                | crate::geometry::mesh::ElementType::Triangle6 => {
                    if n >= 3 {
                        edges_list.push(order_edge(element.nodes[0], element.nodes[1]));
                        edges_list.push(order_edge(element.nodes[1], element.nodes[2]));
                        edges_list.push(order_edge(element.nodes[2], element.nodes[0]));
                    }
                }
                crate::geometry::mesh::ElementType::Quadrilateral
                | crate::geometry::mesh::ElementType::Quadrilateral9 => {
                    if n >= 4 {
                        edges_list.push(order_edge(element.nodes[0], element.nodes[1]));
                        edges_list.push(order_edge(element.nodes[1], element.nodes[2]));
                        edges_list.push(order_edge(element.nodes[2], element.nodes[3]));
                        edges_list.push(order_edge(element.nodes[3], element.nodes[0]));
                    }
                }
                crate::geometry::mesh::ElementType::Tetrahedron
                | crate::geometry::mesh::ElementType::Tetrahedron10 => {
                    if n >= 4 {
                        // Base
                        edges_list.push(order_edge(element.nodes[0], element.nodes[1]));
                        edges_list.push(order_edge(element.nodes[1], element.nodes[2]));
                        edges_list.push(order_edge(element.nodes[2], element.nodes[0]));
                        // Sides to apex
                        edges_list.push(order_edge(element.nodes[0], element.nodes[3]));
                        edges_list.push(order_edge(element.nodes[1], element.nodes[3]));
                        edges_list.push(order_edge(element.nodes[2], element.nodes[3]));
                    }
                }
                crate::geometry::mesh::ElementType::Hexahedron
                | crate::geometry::mesh::ElementType::Hexahedron20 => {
                    if n >= 8 {
                        // Bottom face
                        edges_list.push(order_edge(element.nodes[0], element.nodes[1]));
                        edges_list.push(order_edge(element.nodes[1], element.nodes[2]));
                        edges_list.push(order_edge(element.nodes[2], element.nodes[3]));
                        edges_list.push(order_edge(element.nodes[3], element.nodes[0]));
                        // Top face
                        edges_list.push(order_edge(element.nodes[4], element.nodes[5]));
                        edges_list.push(order_edge(element.nodes[5], element.nodes[6]));
                        edges_list.push(order_edge(element.nodes[6], element.nodes[7]));
                        edges_list.push(order_edge(element.nodes[7], element.nodes[4]));
                        // Vertical edges
                        edges_list.push(order_edge(element.nodes[0], element.nodes[4]));
                        edges_list.push(order_edge(element.nodes[1], element.nodes[5]));
                        edges_list.push(order_edge(element.nodes[2], element.nodes[6]));
                        edges_list.push(order_edge(element.nodes[3], element.nodes[7]));
                    }
                }
                crate::geometry::mesh::ElementType::Pyramid => {
                    if n >= 5 {
                        // Base
                        edges_list.push(order_edge(element.nodes[0], element.nodes[1]));
                        edges_list.push(order_edge(element.nodes[1], element.nodes[2]));
                        edges_list.push(order_edge(element.nodes[2], element.nodes[3]));
                        edges_list.push(order_edge(element.nodes[3], element.nodes[0]));
                        // Sides to apex
                        edges_list.push(order_edge(element.nodes[0], element.nodes[4]));
                        edges_list.push(order_edge(element.nodes[1], element.nodes[4]));
                        edges_list.push(order_edge(element.nodes[2], element.nodes[4]));
                        edges_list.push(order_edge(element.nodes[3], element.nodes[4]));
                    }
                }
                crate::geometry::mesh::ElementType::Prism => {
                    if n >= 6 {
                        // Bottom triangle
                        edges_list.push(order_edge(element.nodes[0], element.nodes[1]));
                        edges_list.push(order_edge(element.nodes[1], element.nodes[2]));
                        edges_list.push(order_edge(element.nodes[2], element.nodes[0]));
                        // Top triangle
                        edges_list.push(order_edge(element.nodes[3], element.nodes[4]));
                        edges_list.push(order_edge(element.nodes[4], element.nodes[5]));
                        edges_list.push(order_edge(element.nodes[5], element.nodes[3]));
                        // Vertical edges
                        edges_list.push(order_edge(element.nodes[0], element.nodes[3]));
                        edges_list.push(order_edge(element.nodes[1], element.nodes[4]));
                        edges_list.push(order_edge(element.nodes[2], element.nodes[5]));
                    }
                }
            }
        }

        // De-duplicate and build mapping
        let mut unique_edges = HashSet::new();
        let mut final_edges = Vec::new();
        for (u, v) in edges_list {
            if unique_edges.insert((u, v)) {
                let idx = final_edges.len();
                final_edges.push(Edge {
                    start: u,
                    end: v,
                    elements: Vec::new(),
                });
                node_to_edges.entry(u).or_insert_with(Vec::new).push(idx);
                node_to_edges.entry(v).or_insert_with(Vec::new).push(idx);
            }
        }

        EdgeConnectivity {
            edges: final_edges,
            node_to_edges,
        }
    }
}

fn order_edge(u: usize, v: usize) -> (usize, usize) {
    if u < v {
        (u, v)
    } else {
        (v, u)
    }
}

#[cfg(test)]
mod tests {
    use crate::geometry::mesh::{ElementType, Mesh};
    use nalgebra::Point3;

    #[test]
    fn test_tetrahedron_edges() {
        let mut mesh = Mesh::new("test_tet".to_string(), 3);
        // Add 4 nodes for tetrahedron
        for i in 0..4 {
            mesh.add_node(Point3::new(i as f64, 0.0, 0.0));
        }

        mesh.add_element(ElementType::Tetrahedron, vec![0, 1, 2, 3], 0);

        let edge_conn = mesh.build_edge_connectivity();
        assert_eq!(edge_conn.edges.len(), 6);

        // Check edges existence
        let edges: Vec<(usize, usize)> = edge_conn.edges.iter()
            .map(|e| (e.start, e.end))
            .collect();

        assert!(edges.contains(&(0, 1)));
        assert!(edges.contains(&(1, 2)));
        assert!(edges.contains(&(0, 2)));
        assert!(edges.contains(&(0, 3)));
        assert!(edges.contains(&(1, 3)));
        assert!(edges.contains(&(2, 3)));
    }

    #[test]
    fn test_hexahedron_edges() {
        let mut mesh = Mesh::new("test_hex".to_string(), 3);
        // Add 8 nodes
        for i in 0..8 {
            mesh.add_node(Point3::new(i as f64, 0.0, 0.0));
        }

        mesh.add_element(ElementType::Hexahedron, vec![0, 1, 2, 3, 4, 5, 6, 7], 0);

        let edge_conn = mesh.build_edge_connectivity();
        assert_eq!(edge_conn.edges.len(), 12);
    }

    #[test]
    fn test_pyramid_edges() {
        let mut mesh = Mesh::new("test_pyr".to_string(), 3);
        for i in 0..5 {
            mesh.add_node(Point3::new(i as f64, 0.0, 0.0));
        }

        mesh.add_element(ElementType::Pyramid, vec![0, 1, 2, 3, 4], 0);

        let edge_conn = mesh.build_edge_connectivity();
        assert_eq!(edge_conn.edges.len(), 8);
    }

    #[test]
    fn test_prism_edges() {
        let mut mesh = Mesh::new("test_prism".to_string(), 3);
        for i in 0..6 {
            mesh.add_node(Point3::new(i as f64, 0.0, 0.0));
        }

        mesh.add_element(ElementType::Prism, vec![0, 1, 2, 3, 4, 5], 0);

        let edge_conn = mesh.build_edge_connectivity();
        assert_eq!(edge_conn.edges.len(), 9);
    }
}
