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
                crate::geometry::mesh::ElementType::Triangle if n == 3 => {
                    edges_list.push(order_edge(element.nodes[0], element.nodes[1]));
                    edges_list.push(order_edge(element.nodes[1], element.nodes[2]));
                    edges_list.push(order_edge(element.nodes[2], element.nodes[0]));
                }
                crate::geometry::mesh::ElementType::Quadrilateral if n == 4 => {
                    edges_list.push(order_edge(element.nodes[0], element.nodes[1]));
                    edges_list.push(order_edge(element.nodes[1], element.nodes[2]));
                    edges_list.push(order_edge(element.nodes[2], element.nodes[3]));
                    edges_list.push(order_edge(element.nodes[3], element.nodes[0]));
                }
                _ => {} // Simplified for now
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
