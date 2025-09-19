//! Mesh connectivity and topology

use super::types::{Element, Mesh};
use nalgebra::RealField;
use std::collections::{HashMap, HashSet};

/// Connectivity information for mesh topology
#[derive(Debug, Clone)]
pub struct Connectivity {
    /// Node to elements mapping
    pub node_to_elements: HashMap<usize, Vec<usize>>,
    /// Element neighbors
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

/// Face connectivity for 3D meshes
#[derive(Debug, Clone)]
pub struct FaceConnectivity {
    /// Faces in the mesh
    pub faces: Vec<Face>,
    /// Element to faces mapping
    pub element_to_faces: HashMap<usize, Vec<usize>>,
}

impl<T: RealField + Copy> Mesh<T> {
    /// Build connectivity information
    #[must_use]
    pub fn build_connectivity(&self) -> Connectivity {
        let mut node_to_elements: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut element_neighbors: HashMap<usize, Vec<usize>> = HashMap::new();
        let boundary_faces = Vec::new();

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
            boundary_faces,
        }
    }

    /// Build edge connectivity
    pub fn build_edge_connectivity(&self) -> EdgeConnectivity {
        let mut edges: Vec<Edge> = Vec::new();
        let mut node_to_edges = HashMap::new();
        let mut edge_map: HashMap<(usize, usize), usize> = HashMap::new();

        for (elem_idx, element) in self.elements.iter().enumerate() {
            // Extract edges from element
            let element_edges = self.extract_element_edges(element);

            for (n1, n2) in element_edges {
                let edge_key = if n1 < n2 { (n1, n2) } else { (n2, n1) };

                if let Some(&edge_idx) = edge_map.get(&edge_key) {
                    edges[edge_idx].elements.push(elem_idx);
                } else {
                    let edge_idx = edges.len();
                    edges.push(Edge {
                        start: edge_key.0,
                        end: edge_key.1,
                        elements: vec![elem_idx],
                    });
                    edge_map.insert(edge_key, edge_idx);

                    node_to_edges
                        .entry(edge_key.0)
                        .or_insert_with(Vec::new)
                        .push(edge_idx);
                    node_to_edges
                        .entry(edge_key.1)
                        .or_insert_with(Vec::new)
                        .push(edge_idx);
                }
            }
        }

        EdgeConnectivity {
            edges,
            node_to_edges,
        }
    }

    fn extract_element_edges(&self, element: &Element) -> Vec<(usize, usize)> {
        let mut edges = Vec::new();
        let n = element.nodes.len();

        match element.element_type {
            crate::domains::mesh_operations::element::ElementType::Triangle if n == 3 => {
                edges.push((element.nodes[0], element.nodes[1]));
                edges.push((element.nodes[1], element.nodes[2]));
                edges.push((element.nodes[2], element.nodes[0]));
            }
            crate::domains::mesh_operations::element::ElementType::Quadrilateral if n == 4 => {
                edges.push((element.nodes[0], element.nodes[1]));
                edges.push((element.nodes[1], element.nodes[2]));
                edges.push((element.nodes[2], element.nodes[3]));
                edges.push((element.nodes[3], element.nodes[0]));
            }
            crate::domains::mesh_operations::element::ElementType::Tetrahedron if n == 4 => {
                for i in 0..n {
                    for j in (i + 1)..n {
                        edges.push((element.nodes[i], element.nodes[j]));
                    }
                }
            }
            _ => {}
        }

        edges
    }
}
