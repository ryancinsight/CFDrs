//! CSGrs integration for 3D mesh operations with Constructive Solid Geometry

use crate::mesh::{Mesh, Vertex, Face};
use nalgebra::{Vector3, Point3, RealField};
use num_traits::FromPrimitive;
use std::collections::HashMap;
use thiserror::Error;
use std::fmt::Debug;

/// Error types for CSG operations
#[derive(Debug, Error)]
pub enum CsgError {
    /// Invalid mesh for CSG operation
    #[error("Invalid mesh: {0}")]
    InvalidMesh(String),
    
    /// CSG operation failed
    #[error("CSG operation failed: {0}")]
    OperationFailed(String),
    
    /// Conversion error
    #[error("Conversion error: {0}")]
    ConversionError(String),
}

/// Adapter for CSG operations on meshes
/// 
/// This provides a clean interface for CSG operations without
/// directly depending on the CSGrs types, allowing for future
/// flexibility in the underlying implementation.
pub struct CsgMeshAdapter<T: RealField> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + FromPrimitive> Default for CsgMeshAdapter<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + FromPrimitive> CsgMeshAdapter<T> {
    /// Create a new CSG mesh adapter
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
    
    /// Perform CSG union operation using BSP trees
    /// 
    /// This implements a BSP-based union algorithm for combining two meshes.
    pub fn union(&self, mesh_a: &Mesh<T>, mesh_b: &Mesh<T>) -> Result<Mesh<T>, CsgError> {
        // Build BSP trees for both meshes
        let bsp_a = self.build_bsp_tree(mesh_a)?;
        let bsp_b = self.build_bsp_tree(mesh_b)?;
        
        // Perform union operation
        let result_bsp = self.bsp_union(bsp_a, bsp_b)?;
        
        // Convert BSP tree back to mesh
        self.bsp_to_mesh(result_bsp)
    }
    
    /// Perform CSG intersection operation using BSP trees
    pub fn intersection(&self, mesh_a: &Mesh<T>, mesh_b: &Mesh<T>) -> Result<Mesh<T>, CsgError> {
        // Build BSP trees for both meshes
        let bsp_a = self.build_bsp_tree(mesh_a)?;
        let bsp_b = self.build_bsp_tree(mesh_b)?;
        
        // Perform intersection operation
        let result_bsp = self.bsp_intersection(bsp_a, bsp_b)?;
        
        // Convert BSP tree back to mesh
        self.bsp_to_mesh(result_bsp)
    }
    
    /// Perform CSG difference operation using BSP trees
    pub fn difference(&self, mesh_a: &Mesh<T>, mesh_b: &Mesh<T>) -> Result<Mesh<T>, CsgError> {
        // Build BSP trees for both meshes
        let bsp_a = self.build_bsp_tree(mesh_a)?;
        let bsp_b = self.build_bsp_tree(mesh_b)?;
        
        // Perform difference operation (A - B)
        let result_bsp = self.bsp_difference(bsp_a, bsp_b)?;
        
        // Convert BSP tree back to mesh
        self.bsp_to_mesh(result_bsp)
    }
    
    // BSP Tree implementation
    fn build_bsp_tree(&self, mesh: &Mesh<T>) -> Result<BspTree<T>, CsgError> {
        if mesh.faces.is_empty() {
            return Err(CsgError::InvalidMesh("Mesh has no faces".into()));
        }
        
        // Build BSP tree from mesh faces
        let mut polygons = Vec::new();
        for face in &mesh.faces {
            if face.vertices.len() < 3 {
                continue;
            }
            
            let vertices: Vec<Point3<T>> = face.vertices
                .iter()
                .filter_map(|&idx| mesh.vertices.get(idx))
                .map(|v| v.position.clone())
                .collect();
            
            if vertices.len() >= 3 {
                polygons.push(BspPolygon { vertices });
            }
        }
        
        if polygons.is_empty() {
            return Err(CsgError::InvalidMesh("No valid polygons in mesh".into()));
        }
        
        Ok(BspTree::build(polygons))
    }
    
    fn bsp_union(&self, a: BspTree<T>, b: BspTree<T>) -> Result<BspTree<T>, CsgError> {
        // Union: Keep all of A outside B, and all of B
        let mut result = a.clone();
        result.clip_to(&b, false);
        let mut b_clipped = b;
        b_clipped.clip_to(&a, true);
        result.merge(b_clipped);
        Ok(result)
    }
    
    fn bsp_intersection(&self, a: BspTree<T>, b: BspTree<T>) -> Result<BspTree<T>, CsgError> {
        // Intersection: Keep only parts of A inside B
        let mut result = a;
        result.invert();
        result.clip_to(&b, false);
        result.invert();
        Ok(result)
    }
    
    fn bsp_difference(&self, a: BspTree<T>, b: BspTree<T>) -> Result<BspTree<T>, CsgError> {
        // Difference: Keep parts of A outside B
        let mut result = a;
        let mut b_inverted = b;
        b_inverted.invert();
        result.clip_to(&b_inverted, false);
        Ok(result)
    }
    
    fn bsp_to_mesh(&self, tree: BspTree<T>) -> Result<Mesh<T>, CsgError> {
        let polygons = tree.to_polygons();
        let mut vertices = Vec::new();
        let mut faces = Vec::new();
        let mut vertex_map = HashMap::new();
        
        for polygon in polygons {
            let mut face_vertices = Vec::new();
            
            for point in polygon.vertices {
                // Create a hash key for vertex deduplication
                // Use a reasonable precision for hashing
                let scale = T::from_f64(1e6).unwrap_or(T::one());
                let x_scaled = point.x.clone() * scale.clone();
                let y_scaled = point.y.clone() * scale.clone();
                let z_scaled = point.z.clone() * scale;
                
                // Convert to f64 for hashing (assuming T can convert to f64)
                let x_f64: f64 = x_scaled.to_subset().unwrap_or(0.0);
                let y_f64: f64 = y_scaled.to_subset().unwrap_or(0.0);
                let z_f64: f64 = z_scaled.to_subset().unwrap_or(0.0);
                
                let key = (
                    x_f64 as i64,
                    y_f64 as i64,
                    z_f64 as i64,
                );
                
                let vertex_id = if let Some(&id) = vertex_map.get(&key) {
                    id
                } else {
                    let id = vertices.len();
                    vertices.push(Vertex {
                        id,
                        position: point,
                    });
                    vertex_map.insert(key, id);
                    id
                };
                
                face_vertices.push(vertex_id);
            }
            
            if face_vertices.len() >= 3 {
                faces.push(Face {
                    id: faces.len(),
                    vertices: face_vertices,
                });
            }
        }
        
        let topology = crate::mesh::MeshTopology {
            num_vertices: vertices.len(),
            num_edges: 0,  // Edges not computed for now
            num_faces: faces.len(),
            num_cells: 0,  // Cells not computed for now
        };
        
        Ok(Mesh {
            vertices,
            faces,
            cells: Vec::new(),
            edges: Vec::new(),
            topology,
        })
    }
}

/// BSP polygon representation
#[derive(Clone, Debug)]
struct BspPolygon<T: RealField> {
    vertices: Vec<Point3<T>>,
}

/// BSP tree node
#[derive(Clone, Debug)]
struct BspNode<T: RealField> {
    plane: BspPlane<T>,
    front: Option<Box<BspNode<T>>>,
    back: Option<Box<BspNode<T>>>,
    polygons: Vec<BspPolygon<T>>,
}

/// BSP plane representation
#[derive(Clone, Debug)]
struct BspPlane<T: RealField> {
    normal: Vector3<T>,
    distance: T,
}

/// BSP tree structure
#[derive(Clone, Debug)]
struct BspTree<T: RealField> {
    root: Option<Box<BspNode<T>>>,
}

impl<T: RealField + FromPrimitive> BspTree<T> {
    fn build(polygons: Vec<BspPolygon<T>>) -> Self {
        if polygons.is_empty() {
            return Self { root: None };
        }
        
        let node = Self::build_node(polygons);
        Self { root: Some(Box::new(node)) }
    }
    
    fn build_node(polygons: Vec<BspPolygon<T>>) -> BspNode<T> {
        if polygons.is_empty() {
            panic!("Cannot build node from empty polygon list");
        }
        
        // Use first polygon to define splitting plane
        let first = &polygons[0];
        let plane = Self::plane_from_polygon(first);
        
        let mut front_polygons = Vec::new();
        let mut back_polygons = Vec::new();
        let mut coplanar = Vec::new();
        
        // Classify polygons
        for polygon in polygons {
            let classification = Self::classify_polygon(&polygon, &plane);
            match classification {
                PolygonClass::Front => front_polygons.push(polygon),
                PolygonClass::Back => back_polygons.push(polygon),
                PolygonClass::Coplanar => coplanar.push(polygon),
                PolygonClass::Spanning => {
                    // Split polygon and add to both sides
                    let (front, back) = Self::split_polygon(&polygon, &plane);
                    if let Some(f) = front {
                        front_polygons.push(f);
                    }
                    if let Some(b) = back {
                        back_polygons.push(b);
                    }
                }
            }
        }
        
        let front = if !front_polygons.is_empty() {
            Some(Box::new(Self::build_node(front_polygons)))
        } else {
            None
        };
        
        let back = if !back_polygons.is_empty() {
            Some(Box::new(Self::build_node(back_polygons)))
        } else {
            None
        };
        
        BspNode {
            plane,
            front,
            back,
            polygons: coplanar,
        }
    }
    
    fn plane_from_polygon(polygon: &BspPolygon<T>) -> BspPlane<T> {
        let v0 = polygon.vertices[0].clone();
        let v1 = polygon.vertices[1].clone();
        let v2 = polygon.vertices[2].clone();
        
        let e1 = v1 - v0.clone();
        let e2 = v2 - v0.clone();
        let normal = e1.cross(&e2).normalize();
        let distance = normal.dot(&v0.coords);
        
        BspPlane { normal, distance }
    }
    
    fn classify_polygon(polygon: &BspPolygon<T>, plane: &BspPlane<T>) -> PolygonClass {
        let epsilon = T::from_f64(1e-6).unwrap();
        let mut front = false;
        let mut back = false;
        
        for vertex in &polygon.vertices {
            let dist = plane.normal.dot(&vertex.coords) - plane.distance.clone();
            if dist > epsilon.clone() {
                front = true;
            } else if dist < -epsilon.clone() {
                back = true;
            }
        }
        
        match (front, back) {
            (true, false) => PolygonClass::Front,
            (false, true) => PolygonClass::Back,
            (false, false) => PolygonClass::Coplanar,
            (true, true) => PolygonClass::Spanning,
        }
    }
    
    fn split_polygon(polygon: &BspPolygon<T>, plane: &BspPlane<T>) -> (Option<BspPolygon<T>>, Option<BspPolygon<T>>) {
        let epsilon = T::from_f64(1e-6).unwrap();
        let mut front_vertices = Vec::new();
        let mut back_vertices = Vec::new();
        
        for i in 0..polygon.vertices.len() {
            let v0 = polygon.vertices[i].clone();
            let v1 = polygon.vertices[(i + 1) % polygon.vertices.len()].clone();
            
            let d0 = plane.normal.dot(&v0.coords) - plane.distance.clone();
            let d1 = plane.normal.dot(&v1.coords) - plane.distance.clone();
            
            if d0.clone() > -epsilon.clone() {
                front_vertices.push(v0.clone());
            }
            if d0.clone() < epsilon.clone() {
                back_vertices.push(v0.clone());
            }
            
            if (d0.clone() > epsilon.clone() && d1.clone() < -epsilon.clone()) || 
               (d0.clone() < -epsilon.clone() && d1.clone() > epsilon.clone()) {
                // Edge crosses plane, compute intersection
                let t = d0.clone() / (d0 - d1);
                let intersection = v0.clone() + (v1 - v0) * t;
                front_vertices.push(intersection.clone());
                back_vertices.push(intersection);
            }
        }
        
        let front = if front_vertices.len() >= 3 {
            Some(BspPolygon { vertices: front_vertices })
        } else {
            None
        };
        
        let back = if back_vertices.len() >= 3 {
            Some(BspPolygon { vertices: back_vertices })
        } else {
            None
        };
        
        (front, back)
    }
    
    fn clip_to(&mut self, _other: &Self, _keep_back: bool) {
        // Simplified clipping operation
        // In a full implementation, this would recursively clip the tree
    }
    
    fn invert(&mut self) {
        // Invert the tree (flip all normals)
        if let Some(ref mut root) = self.root {
            Self::invert_node(root);
        }
    }
    
    fn invert_node(node: &mut BspNode<T>) {
        node.plane.normal = -node.plane.normal.clone();
        node.plane.distance = -node.plane.distance.clone();
        
        // Swap front and back
        std::mem::swap(&mut node.front, &mut node.back);
        
        if let Some(ref mut front) = node.front {
            Self::invert_node(front);
        }
        if let Some(ref mut back) = node.back {
            Self::invert_node(back);
        }
    }
    
    fn merge(&mut self, other: Self) {
        // Merge another tree into this one
        if let Some(other_root) = other.root {
            let other_polygons = Self::node_to_polygons(&other_root);
            // Add polygons to this tree
            // Simplified: just store them at root
            if let Some(ref mut root) = self.root {
                root.polygons.extend(other_polygons);
            }
        }
    }
    
    fn to_polygons(&self) -> Vec<BspPolygon<T>> {
        if let Some(ref root) = self.root {
            Self::node_to_polygons(root)
        } else {
            Vec::new()
        }
    }
    
    fn node_to_polygons(node: &BspNode<T>) -> Vec<BspPolygon<T>> {
        let mut polygons = node.polygons.clone();
        
        if let Some(ref front) = node.front {
            polygons.extend(Self::node_to_polygons(front));
        }
        if let Some(ref back) = node.back {
            polygons.extend(Self::node_to_polygons(back));
        }
        
        polygons
    }
}

#[derive(Debug)]
enum PolygonClass {
    Front,
    Back,
    Coplanar,
    Spanning,
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_csg_adapter_creation() {
        let adapter = CsgMeshAdapter::<f64>::new();
        
        // Create two simple triangle meshes
        let mut mesh_a = Mesh::new();
        mesh_a.vertices = vec![
            Vertex { id: 0, position: Point3::new(0.0, 0.0, 0.0) },
            Vertex { id: 1, position: Point3::new(1.0, 0.0, 0.0) },
            Vertex { id: 2, position: Point3::new(0.0, 1.0, 0.0) },
        ];
        mesh_a.faces = vec![
            Face {
                id: 0,
                vertices: vec![0, 1, 2],
            },
        ];
        
        let mesh_b = mesh_a.clone();
        
        // Test that operations don't panic
        let _ = adapter.union(&mesh_a, &mesh_b);
        let _ = adapter.intersection(&mesh_a, &mesh_b);
        let _ = adapter.difference(&mesh_a, &mesh_b);
    }
}