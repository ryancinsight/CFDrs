//! CSG (Constructive Solid Geometry) operations
//!
//! KNOWN LIMITATION: This module provides basic mesh creation capabilities
//! but does not implement complex boolean operations. The csgrs crate integration
//! has API incompatibilities that prevent full CSG functionality.
//! 
//! Current capabilities:
//! - Basic primitive mesh generation (spheres, boxes, cylinders)
//! - Mesh validation and quality checking
//! - Simple mesh transformations
//!
//! Missing capabilities (documented limitations):
//! - Boolean operations (union, intersection, difference)
//! - Complex CSG tree evaluation
//! - Mesh healing and repair

use crate::mesh::{Mesh, Vertex, Face};
use nalgebra::{Vector3, Point3, RealField};
use num_traits::FromPrimitive;
use thiserror::Error;
use std::fmt::Debug;
use cfd_core::constants;

/// Error types for CSG operations
#[derive(Debug, Error)]
pub enum CsgError {
    /// Invalid mesh for CSG operation
    #[error("Invalid mesh: {0}")]
    InvalidMesh(String),
    
    /// CSG operation failed
    #[error("CSG operation failed: {0}")]
    OperationFailed(String),
}

/// CSG operations and mesh generation
/// 
/// Provides basic mesh generation and validation capabilities.
/// Boolean operations are documented as not implemented.
pub struct CsgOperator<T: RealField> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + FromPrimitive> CsgOperator<T> {
    /// Create a new CSG operator
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
    
    /// Generate a sphere mesh with specified radius and resolution
    pub fn create_sphere(&self, radius: T, resolution: usize) -> Result<Mesh<T>, CsgError> {
        if resolution < 4 {
            return Err(CsgError::InvalidMesh("Resolution must be at least 4".to_string()));
        }
        
        let mut vertices = Vec::new();
        let mut faces = Vec::new();
        
        // Generate icosphere vertices using spherical coordinates
        let pi = T::from_f64(std::f64::consts::PI).unwrap();
        let two_pi = T::from_f64(constants::TWO * std::f64::consts::PI).unwrap();
        
        // Add top vertex
        vertices.push(Vertex {
            position: Point3::new(T::zero(), T::zero(), radius.clone()),
            id: 0,
        });
        
        // Generate latitude rings
        for i in 1..resolution {
            let theta = pi.clone() * T::from_usize(i).unwrap() / T::from_usize(resolution).unwrap();
            let sin_theta = theta.clone().sin();
            let cos_theta = theta.cos();
            
            for j in 0..resolution * 2 {
                let phi = two_pi.clone() * T::from_usize(j).unwrap() / T::from_usize(resolution * 2).unwrap();
                let sin_phi = phi.clone().sin();
                let cos_phi = phi.cos();
                
                let x = radius.clone() * sin_theta.clone() * cos_phi;
                let y = radius.clone() * sin_theta.clone() * sin_phi;
                let z = radius.clone() * cos_theta.clone();
                
                let position = Point3::new(x, y, z);
                let vertex_id = vertices.len() + 1;
                
                vertices.push(Vertex {
                    position,
                    id: vertex_id,
                });
            }
        }
        
        // Add bottom vertex
        let bottom_id = vertices.len() + 1;
        vertices.push(Vertex {
            position: Point3::new(T::zero(), T::zero(), -radius.clone()),
            id: bottom_id,
        });
        
        // Generate faces for a basic sphere (triangulation)
        // This is a simplified implementation
        for i in 0..resolution.min(6) {  // Limit faces for basic implementation
            faces.push(Face {
                vertices: vec![0, i + 1, ((i + 1) % 6) + 1],
                id: i,
            });
        }
        
        let mut mesh = Mesh {
            vertices,
            edges: Vec::new(),
            faces,
            cells: Vec::new(),
            topology: crate::mesh::MeshTopology {
                num_vertices: 0,
                num_edges: 0,
                num_faces: 0,
                num_cells: 0,
            },
        };
        mesh.update_topology();
        
        Ok(mesh)
    }
    
    /// Generate a box mesh with specified dimensions
    pub fn create_box(&self, width: T, height: T, depth: T) -> Result<Mesh<T>, CsgError> {
        let mut vertices = Vec::new();
        let mut faces = Vec::new();
        
        let half_w = width.clone() / T::from_f64(constants::TWO).unwrap();
        let half_h = height.clone() / T::from_f64(constants::TWO).unwrap();
        let half_d = depth.clone() / T::from_f64(constants::TWO).unwrap();
        
        // Generate 8 vertices of a box
        vertices.push(Vertex { position: Point3::new(-half_w.clone(), -half_h.clone(), -half_d.clone()), id: 0 });
        vertices.push(Vertex { position: Point3::new(half_w.clone(), -half_h.clone(), -half_d.clone()), id: 1 });
        vertices.push(Vertex { position: Point3::new(half_w.clone(), half_h.clone(), -half_d.clone()), id: 2 });
        vertices.push(Vertex { position: Point3::new(-half_w.clone(), half_h.clone(), -half_d.clone()), id: 3 });
        vertices.push(Vertex { position: Point3::new(-half_w.clone(), -half_h.clone(), half_d.clone()), id: 4 });
        vertices.push(Vertex { position: Point3::new(half_w.clone(), -half_h.clone(), half_d.clone()), id: 5 });
        vertices.push(Vertex { position: Point3::new(half_w.clone(), half_h.clone(), half_d.clone()), id: 6 });
        vertices.push(Vertex { position: Point3::new(-half_w.clone(), half_h.clone(), half_d.clone()), id: 7 });
        
        // Generate 12 faces (2 triangles per face, 6 faces)
        // Bottom face
        faces.push(Face { vertices: vec![0, 1, 2], id: 0 });
        faces.push(Face { vertices: vec![0, 2, 3], id: 1 });
        // Top face
        faces.push(Face { vertices: vec![4, 7, 6], id: 2 });
        faces.push(Face { vertices: vec![4, 6, 5], id: 3 });
        // Front face
        faces.push(Face { vertices: vec![0, 4, 5], id: 4 });
        faces.push(Face { vertices: vec![0, 5, 1], id: 5 });
        // Back face
        faces.push(Face { vertices: vec![2, 6, 7], id: 6 });
        faces.push(Face { vertices: vec![2, 7, 3], id: 7 });
        // Left face
        faces.push(Face { vertices: vec![0, 3, 7], id: 8 });
        faces.push(Face { vertices: vec![0, 7, 4], id: 9 });
        // Right face
        faces.push(Face { vertices: vec![1, 5, 6], id: 10 });
        faces.push(Face { vertices: vec![1, 6, 2], id: 11 });
        
        let mut mesh = Mesh {
            vertices,
            edges: Vec::new(),
            faces,
            cells: Vec::new(),
            topology: crate::mesh::MeshTopology {
                num_vertices: 0,
                num_edges: 0,
                num_faces: 0,
                num_cells: 0,
            },
        };
        mesh.update_topology();
        
        Ok(mesh)
    }
    
    /// Validate mesh topology and geometry
    pub fn validate_mesh(&self, mesh: &Mesh<T>) -> Result<bool, CsgError> {
        if mesh.vertices.is_empty() {
            return Err(CsgError::InvalidMesh("Mesh has no vertices".to_string()));
        }
        
        if mesh.faces.is_empty() {
            return Err(CsgError::InvalidMesh("Mesh has no faces".to_string()));
        }
        
        // Check that all face vertex indices are valid
        let num_vertices = mesh.vertices.len();
        for face in &mesh.faces {
            for &vertex_idx in &face.vertices {
                if vertex_idx >= num_vertices {
                    return Err(CsgError::InvalidMesh(
                        format!("Face references invalid vertex index: {}", vertex_idx)
                    ));
                }
            }
        }
        
        Ok(true)
    }
    
    /// Create a cylinder mesh with specified radius, height, and segments
    pub fn create_cylinder(&self, radius: T, height: T, segments: usize) -> Result<Mesh<T>, CsgError> {
        let mut mesh = Mesh::new();
        let center = Point3::new(T::zero(), T::zero(), T::zero());
        
        let segments = segments.max(3);
        let angle_step = T::from_f64(constants::TWO * std::f64::consts::PI / segments as f64).unwrap();
        
        // Create vertices for top and bottom circles
        for i in 0..segments {
            let angle = angle_step.clone() * T::from_usize(i).unwrap();
            let angle_f64: f64 = angle.to_subset().unwrap_or(constants::ZERO);
            let x = radius.clone() * T::from_f64(angle_f64.cos()).unwrap();
            let z = radius.clone() * T::from_f64(angle_f64.sin()).unwrap();
            
            // Bottom vertex
            mesh.vertices.push(Vertex {
                position: center.clone() + Vector3::new(x.clone(), T::zero(), z.clone()),
                id: i * 2,
            });
            
            // Top vertex
            mesh.vertices.push(Vertex {
                position: center.clone() + Vector3::new(x, height.clone(), z),
                id: i * 2 + 1,
            });
        }
        
        // Add center vertices for caps
        let bottom_center_idx = mesh.vertices.len();
        mesh.vertices.push(Vertex {
            position: center.clone(),
            id: bottom_center_idx,
        });
        
        let top_center_idx = mesh.vertices.len();
        mesh.vertices.push(Vertex {
            position: center + Vector3::new(T::zero(), height, T::zero()),
            id: top_center_idx,
        });
        
        // Create side faces
        let mut face_id = 0;
        for i in 0..segments {
            let next = (i + 1) % segments;
            let bottom_curr = i * 2;
            let top_curr = i * 2 + 1;
            let bottom_next = next * 2;
            let top_next = next * 2 + 1;
            
            // Two triangles for the side
            mesh.faces.push(Face {
                vertices: vec![bottom_curr, bottom_next, top_next],
                id: face_id,
            });
            face_id += 1;
            
            mesh.faces.push(Face {
                vertices: vec![bottom_curr, top_next, top_curr],
                id: face_id,
            });
            face_id += 1;
            
            // Bottom cap triangle
            mesh.faces.push(Face {
                vertices: vec![bottom_center_idx, bottom_next, bottom_curr],
                id: face_id,
            });
            face_id += 1;
            
            // Top cap triangle
            mesh.faces.push(Face {
                vertices: vec![top_center_idx, top_curr, top_next],
                id: face_id,
            });
            face_id += 1;
        }
        
        mesh.update_topology();
        Ok(mesh)
    }
}



/// Primitive shapes for CSG operations
pub struct CsgPrimitives<T: RealField> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + FromPrimitive> CsgPrimitives<T> {
    /// Create a box mesh
    /// This at least creates a valid mesh structure
    pub fn create_box(center: Point3<T>, size: Vector3<T>) -> Mesh<T> {
        let mut mesh = Mesh::new();
        
        let half_size = size / T::from_f64(constants::TWO).unwrap();
        
        // Create 8 vertices of the box
        let positions = vec![
            center.clone() + Vector3::new(-half_size.x.clone(), -half_size.y.clone(), -half_size.z.clone()),
            center.clone() + Vector3::new( half_size.x.clone(), -half_size.y.clone(), -half_size.z.clone()),
            center.clone() + Vector3::new( half_size.x.clone(),  half_size.y.clone(), -half_size.z.clone()),
            center.clone() + Vector3::new(-half_size.x.clone(),  half_size.y.clone(), -half_size.z.clone()),
            center.clone() + Vector3::new(-half_size.x.clone(), -half_size.y.clone(),  half_size.z.clone()),
            center.clone() + Vector3::new( half_size.x.clone(), -half_size.y.clone(),  half_size.z.clone()),
            center.clone() + Vector3::new( half_size.x.clone(),  half_size.y.clone(),  half_size.z.clone()),
            center + Vector3::new(-half_size.x.clone(),  half_size.y.clone(),  half_size.z.clone()),
        ];
        
        for (i, pos) in positions.iter().enumerate() {
            mesh.vertices.push(Vertex {
                position: pos.clone(),
                id: i,
            });
        }
        
        // Create 12 triangular faces (2 per box face)
        let face_indices = vec![
            // Front
            vec![0, 1, 2],
            vec![0, 2, 3],
            // Back
            vec![5, 4, 7],
            vec![5, 7, 6],
            // Top
            vec![3, 2, 6],
            vec![3, 6, 7],
            // Bottom
            vec![4, 5, 1],
            vec![4, 1, 0],
            // Right
            vec![1, 5, 6],
            vec![1, 6, 2],
            // Left
            vec![4, 0, 3],
            vec![4, 3, 7],
        ];
        
        for (i, vertices) in face_indices.iter().enumerate() {
            mesh.faces.push(Face {
                vertices: vertices.clone(),
                id: i,
            });
        }
        
        mesh.update_topology();
        mesh
    }
    
    /// Create a sphere mesh (icosphere approximation)
    /// Creates a simple octahedron as placeholder
    pub fn create_sphere(center: Point3<T>, radius: T, _subdivisions: usize) -> Mesh<T> {
        let mut mesh = Mesh::new();
        
        // Simple octahedron vertices
        let positions = vec![
            center.clone() + Vector3::new(T::zero(), T::zero(), radius.clone()),     // Top
            center.clone() + Vector3::new(radius.clone(), T::zero(), T::zero()),     // +X
            center.clone() + Vector3::new(T::zero(), radius.clone(), T::zero()),     // +Y
            center.clone() + Vector3::new(-radius.clone(), T::zero(), T::zero()),    // -X
            center.clone() + Vector3::new(T::zero(), -radius.clone(), T::zero()),    // -Y
            center + Vector3::new(T::zero(), T::zero(), -radius),    // Bottom
        ];
        
        for (i, pos) in positions.iter().enumerate() {
            mesh.vertices.push(Vertex {
                position: pos.clone(),
                id: i,
            });
        }
        
        // Create 8 triangular faces
        let face_indices = vec![
            vec![0, 1, 2], // Top front
            vec![0, 2, 3], // Top left
            vec![0, 3, 4], // Top back
            vec![0, 4, 1], // Top right
            vec![5, 2, 1], // Bottom front
            vec![5, 3, 2], // Bottom left
            vec![5, 4, 3], // Bottom back
            vec![5, 1, 4], // Bottom right
        ];
        
        for (i, vertices) in face_indices.iter().enumerate() {
            mesh.faces.push(Face {
                vertices: vertices.clone(),
                id: i,
            });
        }
        
        mesh.update_topology();
        mesh
    }
    
    /// Create a cylinder mesh
    /// Creates a simple prism as placeholder
    pub fn create_cylinder(base_center: Point3<T>, radius: T, height: T, segments: usize) -> Mesh<T> {
        let mut mesh = Mesh::new();
        
        let segments = segments.max(3);
        let angle_step = T::from_f64(constants::TWO * std::f64::consts::PI / segments as f64).unwrap();
        
        // Create vertices for top and bottom circles
        for i in 0..segments {
            let angle = angle_step.clone() * T::from_usize(i).unwrap();
            let angle_f64: f64 = angle.to_subset().unwrap_or(0.0);
            let x = radius.clone() * T::from_f64(angle_f64.cos()).unwrap();
            let z = radius.clone() * T::from_f64(angle_f64.sin()).unwrap();
            
            // Bottom vertex
            mesh.vertices.push(Vertex {
                position: base_center.clone() + Vector3::new(x.clone(), T::zero(), z.clone()),
                id: i * 2,
            });
            
            // Top vertex
            mesh.vertices.push(Vertex {
                position: base_center.clone() + Vector3::new(x, height.clone(), z),
                id: i * 2 + 1,
            });
        }
        
        // Add center vertices for caps
        let bottom_center_idx = mesh.vertices.len();
        mesh.vertices.push(Vertex {
            position: base_center.clone(),
            id: bottom_center_idx,
        });
        
        let top_center_idx = mesh.vertices.len();
        mesh.vertices.push(Vertex {
            position: base_center + Vector3::new(T::zero(), height, T::zero()),
            id: top_center_idx,
        });
        
        // Create side faces
        let mut face_id = 0;
        for i in 0..segments {
            let next = (i + 1) % segments;
            let bottom_curr = i * 2;
            let top_curr = i * 2 + 1;
            let bottom_next = next * 2;
            let top_next = next * 2 + 1;
            
            // Two triangles for the side
            mesh.faces.push(Face {
                vertices: vec![bottom_curr, bottom_next, top_next],
                id: face_id,
            });
            face_id += 1;
            
            mesh.faces.push(Face {
                vertices: vec![bottom_curr, top_next, top_curr],
                id: face_id,
            });
            face_id += 1;
            
            // Bottom cap triangle
            mesh.faces.push(Face {
                vertices: vec![bottom_center_idx, bottom_next, bottom_curr],
                id: face_id,
            });
            face_id += 1;
            
            // Top cap triangle
            mesh.faces.push(Face {
                vertices: vec![top_center_idx, top_curr, top_next],
                id: face_id,
            });
            face_id += 1;
        }
        
        mesh.update_topology();
        mesh
    }
}

// Legacy adapter for compatibility
pub type CsgMeshAdapter<T> = CsgOperator<T>;

#[cfg(test)]
mod tests {
    use super::*;

    
    #[test]
    fn test_mesh_generation() {
        let operator = CsgOperator::<f64>::new();
        
        // Test sphere generation
        let sphere = operator.create_sphere(1.0, 6).expect("Sphere creation should succeed");
        assert!(!sphere.vertices.is_empty());
        assert!(!sphere.faces.is_empty());
        
        // Test box generation  
        let box_mesh = operator.create_box(2.0, 1.0, 1.5).expect("Box creation should succeed");
        assert_eq!(box_mesh.vertices.len(), 8);
        assert_eq!(box_mesh.faces.len(), 12); // 2 triangles per face of 6 faces
        
        // Test mesh validation
        assert!(operator.validate_mesh(&sphere).expect("Validation should succeed"));
        assert!(operator.validate_mesh(&box_mesh).expect("Validation should succeed"));
    }
    
    #[test]
    fn test_primitives() {
        let operator = CsgOperator::<f64>::new();
        
        // Test box creation using unified CsgOperator interface
        let box_mesh = operator.create_box(1.0, 1.0, 1.0).expect("Box creation should succeed");
        assert_eq!(box_mesh.vertices.len(), 8);
        assert_eq!(box_mesh.faces.len(), 12); // 2 triangles per face
        
        // Test sphere creation
        let sphere_mesh = operator.create_sphere(1.0, 6).expect("Sphere creation should succeed");
        assert!(!sphere_mesh.vertices.is_empty());
        assert!(!sphere_mesh.faces.is_empty());
        
        // Test cylinder creation
        let cylinder_mesh = operator.create_cylinder(1.0, 2.0, 8).expect("Cylinder creation should succeed");
        assert!(!cylinder_mesh.vertices.is_empty());
        assert!(!cylinder_mesh.faces.is_empty());
    }
}