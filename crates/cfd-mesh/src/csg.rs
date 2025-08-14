//! CSG (Constructive Solid Geometry) operations using the csgrs library
//!
//! This module provides comprehensive CSG mesh generation and boolean operations
//! using the csgrs crate, enabling creation of complex geometries through
//! constructive solid geometry techniques.
//! 
//! Features:
//! - Full CSG boolean operations (union, intersection, difference, XOR)
//! - Primitive geometry generation:
//!   - Boxes/Cuboids
//!   - Spheres
//!   - Cylinders
//!   - Frustums (truncated cones)
//!   - Cones (frustums with top radius = 0)
//! - Mesh transformations (translate, rotate, scale)
//! - STL export functionality
//! - Integration with CFD mesh structures

use crate::mesh::{Mesh as CfdMesh, Vertex, Face};
use nalgebra::{Vector3, Point3, RealField};
use num_traits::{FromPrimitive, ToPrimitive};
use thiserror::Error;
use std::fmt::Debug;
// Unused import removed
use csgrs::mesh::Mesh as CsgMesh;
use csgrs::traits::CSG;

/// Error types for CSG operations
#[derive(Debug, Error)]
pub enum CsgError {
    /// Invalid mesh for CSG operation
    #[error("Invalid mesh: {0}")]
    InvalidMesh(String),
    
    /// CSG operation failed
    #[error("CSG operation failed: {0}")]
    OperationFailed(String),
    
    /// Invalid geometry parameters
    #[error("Invalid geometry parameters: {0}")]
    InvalidParameters(String),
    
    /// STL export error
    #[error("STL export failed: {0}")]
    ExportError(String),
}

/// CSG operations and mesh generation using csgrs
/// 
/// Provides comprehensive CSG functionality including boolean operations,
/// primitive generation, and mesh transformations.
pub struct CsgOperator<T: RealField> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + FromPrimitive + ToPrimitive> CsgOperator<T> {
    /// Create a new CSG operator
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
    
    /// Create a cube with specified dimensions
    pub fn create_cube(&self, width: T, height: T, depth: T) -> Result<CsgGeometry<T>, CsgError> {
        let w = width.to_f64().ok_or_else(|| CsgError::InvalidParameters("Invalid width".to_string()))?;
        let h = height.to_f64().ok_or_else(|| CsgError::InvalidParameters("Invalid height".to_string()))?;
        let d = depth.to_f64().ok_or_else(|| CsgError::InvalidParameters("Invalid depth".to_string()))?;
        
        let zero = T::from_f64(0.0).unwrap_or_else(|| T::zero());
        if width <= zero || height <= zero || depth <= zero {
            return Err(CsgError::InvalidParameters("Dimensions must be positive".to_string()));
        }
        
        // Use csgrs Mesh::cuboid for non-uniform dimensions
        let csg = CsgMesh::<()>::cuboid(w, h, d, None);
        Ok(CsgGeometry::new(csg))
    }
    
    /// Create a sphere with specified radius and resolution
    pub fn create_sphere(&self, radius: T, horizontal_segments: usize, vertical_segments: usize) -> Result<CsgGeometry<T>, CsgError> {
        let r = radius.to_f64().ok_or_else(|| CsgError::InvalidParameters("Invalid radius".to_string()))?;
        
        let zero = T::from_f64(0.0).unwrap_or_else(|| T::zero());
        if radius <= zero {
            return Err(CsgError::InvalidParameters("Radius must be positive".to_string()));
        }
        
        if horizontal_segments < 3 || vertical_segments < 2 {
            return Err(CsgError::InvalidParameters("Insufficient resolution for sphere".to_string()));
        }
        
        let csg = CsgMesh::<()>::sphere(r, horizontal_segments, vertical_segments, None);
        Ok(CsgGeometry::new(csg))
    }
    
    /// Create a cylinder with specified radius, height, and resolution
    pub fn create_cylinder(&self, radius: T, height: T, segments: usize) -> Result<CsgGeometry<T>, CsgError> {
        let r = radius.to_f64().ok_or_else(|| CsgError::InvalidParameters("Invalid radius".to_string()))?;
        let h = height.to_f64().ok_or_else(|| CsgError::InvalidParameters("Invalid height".to_string()))?;
        
        let zero = T::from_f64(0.0).unwrap_or_else(|| T::zero());
        if radius <= zero || height <= zero {
            return Err(CsgError::InvalidParameters("Radius and height must be positive".to_string()));
        }
        
        if segments < 3 {
            return Err(CsgError::InvalidParameters("Insufficient segments for cylinder".to_string()));
        }
        
        let csg = CsgMesh::<()>::cylinder(r, h, segments, None);
        Ok(CsgGeometry::new(csg))
    }
    
    /// Create a frustum (truncated cone) with specified base radius, top radius, height, and resolution
    /// For a cone (top_radius = 0), this creates a proper cone shape
    pub fn create_frustum(&self, base_radius: T, top_radius: T, height: T, segments: usize) -> Result<CsgGeometry<T>, CsgError> {
        let r1 = base_radius.to_f64().ok_or_else(|| CsgError::InvalidParameters("Invalid base radius".to_string()))?;
        let r2 = top_radius.to_f64().ok_or_else(|| CsgError::InvalidParameters("Invalid top radius".to_string()))?;
        let h = height.to_f64().ok_or_else(|| CsgError::InvalidParameters("Invalid height".to_string()))?;
        
        let zero = T::from_f64(0.0).unwrap_or_else(|| T::zero());
        if base_radius < zero || top_radius < zero || height <= zero {
            return Err(CsgError::InvalidParameters("Radii must be non-negative and height positive".to_string()));
        }
        
        if segments < 3 {
            return Err(CsgError::InvalidParameters("Insufficient segments for frustum".to_string()));
        }
        
        // Use csgrs frustum function which properly handles both cones and truncated cones
        // frustum_ptp creates a frustum from point to point
        let bottom = nalgebra::Point3::new(0.0, 0.0, 0.0);
        let top = nalgebra::Point3::new(0.0, 0.0, h);
        
        let csg = CsgMesh::<()>::frustum_ptp(bottom, top, r1, r2, segments, None);
        
        Ok(CsgGeometry::new(csg))
    }
    
    /// Create a cone (special case of frustum with top radius = 0)
    pub fn create_cone(&self, base_radius: T, height: T, segments: usize) -> Result<CsgGeometry<T>, CsgError> {
        self.create_frustum(base_radius, T::zero(), height, segments)
    }
}

impl<T: RealField + FromPrimitive + ToPrimitive> Default for CsgOperator<T> {
    fn default() -> Self {
        Self::new()
    }
}

/// Wrapper for CSG geometry with type safety and CFD integration
pub struct CsgGeometry<T: RealField> {
    csg: CsgMesh<()>,
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + FromPrimitive + ToPrimitive> CsgGeometry<T> {
    /// Create a new CSG geometry wrapper
    pub fn new(csg: CsgMesh<()>) -> Self {
        Self {
            csg,
            _phantom: std::marker::PhantomData,
        }
    }
    
    /// Get the underlying CSG mesh
    pub fn inner(&self) -> &CsgMesh<()> {
        &self.csg
    }
    
    /// Perform union operation with another geometry
    pub fn union(&self, other: &CsgGeometry<T>) -> CsgGeometry<T> {
        let result = self.csg.union(&other.csg);
        CsgGeometry::new(result)
    }
    
    /// Perform difference operation (subtract other from self)
    pub fn difference(&self, other: &CsgGeometry<T>) -> CsgGeometry<T> {
        let result = self.csg.difference(&other.csg);
        CsgGeometry::new(result)
    }
    
    /// Perform intersection operation with another geometry
    pub fn intersection(&self, other: &CsgGeometry<T>) -> CsgGeometry<T> {
        let result = self.csg.intersection(&other.csg);
        CsgGeometry::new(result)
    }
    
    /// Perform XOR operation with another geometry
    pub fn xor(&self, other: &CsgGeometry<T>) -> CsgGeometry<T> {
        let result = self.csg.xor(&other.csg);
        CsgGeometry::new(result)
    }
    
    /// Translate the geometry by the given vector
    pub fn translate(&mut self, translation: &Vector3<T>) -> Result<(), CsgError> {
        let tx = translation.x.to_f64().ok_or_else(|| CsgError::InvalidParameters("Invalid translation x".to_string()))?;
        let ty = translation.y.to_f64().ok_or_else(|| CsgError::InvalidParameters("Invalid translation y".to_string()))?;
        let tz = translation.z.to_f64().ok_or_else(|| CsgError::InvalidParameters("Invalid translation z".to_string()))?;
        
        self.csg = self.csg.translate(tx, ty, tz);
        Ok(())
    }
    
    /// Rotate the geometry around axes by given angles in degrees
    pub fn rotate(&mut self, x_deg: T, y_deg: T, z_deg: T) -> Result<(), CsgError> {
        let x = x_deg.to_f64().ok_or_else(|| CsgError::InvalidParameters("Invalid x rotation".to_string()))?;
        let y = y_deg.to_f64().ok_or_else(|| CsgError::InvalidParameters("Invalid y rotation".to_string()))?;
        let z = z_deg.to_f64().ok_or_else(|| CsgError::InvalidParameters("Invalid z rotation".to_string()))?;
        
        self.csg = self.csg.rotate(x, y, z);
        Ok(())
    }
    
    /// Scale the geometry by the given factors
    pub fn scale(&mut self, scale_x: T, scale_y: T, scale_z: T) -> Result<(), CsgError> {
        let sx = scale_x.to_f64().ok_or_else(|| CsgError::InvalidParameters("Invalid scale x".to_string()))?;
        let sy = scale_y.to_f64().ok_or_else(|| CsgError::InvalidParameters("Invalid scale y".to_string()))?;
        let sz = scale_z.to_f64().ok_or_else(|| CsgError::InvalidParameters("Invalid scale z".to_string()))?;
        
        let zero = T::from_f64(0.0).unwrap_or_else(|| T::zero());
        if scale_x <= zero || scale_y <= zero || scale_z <= zero {
            return Err(CsgError::InvalidParameters("Scale factors must be positive".to_string()));
        }
        
        self.csg = self.csg.scale(sx, sy, sz);
        Ok(())
    }
    
    /// Mirror the geometry across a plane defined by a normal vector
    pub fn mirror(&mut self, plane: csgrs::mesh::plane::Plane) -> Result<(), CsgError> {
        self.csg = self.csg.mirror(plane);
        Ok(())
    }
    
    /// Convert to CFD mesh format
    pub fn to_mesh(&self) -> Result<CfdMesh<T>, CsgError> {
        let polygons = &self.csg.polygons;
        let mut vertices = Vec::new();
        let mut faces = Vec::new();
        let mut vertex_map = std::collections::HashMap::new();
        let mut vertex_id = 1;
        
        for polygon in polygons {
            let mut face_vertices = Vec::new();
            
            for vertex in &polygon.vertices {
                let pos = vertex.pos;
                let point = Point3::new(
                    T::from_f64(pos.x as f64).ok_or_else(|| CsgError::OperationFailed("Invalid vertex x coordinate".to_string()))?,
                    T::from_f64(pos.y as f64).ok_or_else(|| CsgError::OperationFailed("Invalid vertex y coordinate".to_string()))?,
                    T::from_f64(pos.z as f64).ok_or_else(|| CsgError::OperationFailed("Invalid vertex z coordinate".to_string()))?,
                );
                
                // Use a key based on position with small tolerance for deduplication
                let scale = 1e6;
                let key = (
                    (pos.x as f64 * scale).round() as i64,
                    (pos.y as f64 * scale).round() as i64,
                    (pos.z as f64 * scale).round() as i64,
                );
                
                let id = if let Some(&existing_id) = vertex_map.get(&key) {
                    existing_id
                } else {
                    let id = vertex_id;
                    vertices.push(Vertex {
                        position: point,
                        id,
                    });
                    vertex_map.insert(key, id);
                    vertex_id += 1;
                    id
                };
                
                face_vertices.push(id);
            }
            
            if face_vertices.len() >= 3 {
                faces.push(Face {
                    vertices: face_vertices,
                    id: faces.len() + 1,
                });
            }
        }
        
        if vertices.is_empty() {
            return Err(CsgError::OperationFailed("No vertices generated".to_string()));
        }
        
        // Create mesh manually since Mesh::new() doesn't take parameters
        let mut mesh = CfdMesh::new();
        mesh.vertices = vertices;
        mesh.faces = faces;
        mesh.update_topology();
        
        Ok(mesh)
    }
    
    /// Export to STL format
    pub fn to_stl(&self, _name: &str) -> Result<String, CsgError> {
        // Use csgrs built-in STL export
        #[cfg(feature = "stl-io")]
        {
            use csgrs::io::stl::write_stl_ascii;
            let stl_string = write_stl_ascii(&self.csg.triangulate(), _name);
            Ok(stl_string)
        }
        
        #[cfg(not(feature = "stl-io"))]
        {
            // Manual STL generation if feature not enabled
            let mut stl = format!("solid {}\n", _name);
            
            // Triangulate the mesh first
            let triangulated = self.csg.triangulate();
            
            for polygon in &triangulated.polygons {
                if polygon.vertices.len() == 3 {
                    // Calculate normal (assuming CCW winding)
                    let v0 = &polygon.vertices[0].pos;
                    let v1 = &polygon.vertices[1].pos;
                    let v2 = &polygon.vertices[2].pos;
                    
                    let edge1 = nalgebra::Vector3::new(
                        (v1.x - v0.x) as f32,
                        (v1.y - v0.y) as f32,
                        (v1.z - v0.z) as f32,
                    );
                    let edge2 = nalgebra::Vector3::new(
                        (v2.x - v0.x) as f32,
                        (v2.y - v0.y) as f32,
                        (v2.z - v0.z) as f32,
                    );
                    let normal = edge1.cross(&edge2).normalize();
                    
                    stl.push_str(&format!("  facet normal {} {} {}\n", normal.x, normal.y, normal.z));
                    stl.push_str("    outer loop\n");
                    
                    for vertex in &polygon.vertices {
                        stl.push_str(&format!("      vertex {} {} {}\n", 
                            vertex.pos.x, vertex.pos.y, vertex.pos.z));
                    }
                    
                    stl.push_str("    endloop\n");
                    stl.push_str("  endfacet\n");
                }
            }
            
            stl.push_str(&format!("endsolid {}\n", _name));
            Ok(stl)
        }
    }
    
    /// Export to binary STL format
    pub fn to_stl_binary(&self, name: &str) -> Result<Vec<u8>, CsgError> {
        // For binary STL, we need to implement the format manually
        // Binary STL format:
        // - 80 byte header
        // - 4 byte unsigned integer (number of triangles)
        // - For each triangle:
        //   - 12 bytes (3 floats) for normal
        //   - 36 bytes (9 floats) for 3 vertices
        //   - 2 bytes attribute byte count
        
        let triangulated = self.csg.triangulate();
        let num_triangles = triangulated.polygons.len() as u32;
        
        let mut buffer = Vec::new();
        
        // Write header (80 bytes)
        let header = format!("Binary STL from csgrs: {}", name);
        let mut header_bytes = header.as_bytes().to_vec();
        header_bytes.resize(80, 0);
        buffer.extend_from_slice(&header_bytes);
        
        // Write number of triangles
        buffer.extend_from_slice(&num_triangles.to_le_bytes());
        
        // Write each triangle
        for polygon in &triangulated.polygons {
            if polygon.vertices.len() == 3 {
                // Calculate normal
                let v0 = &polygon.vertices[0].pos;
                let v1 = &polygon.vertices[1].pos;
                let v2 = &polygon.vertices[2].pos;
                
                let edge1 = nalgebra::Vector3::new(
                    (v1.x - v0.x) as f32,
                    (v1.y - v0.y) as f32,
                    (v1.z - v0.z) as f32,
                );
                let edge2 = nalgebra::Vector3::new(
                    (v2.x - v0.x) as f32,
                    (v2.y - v0.y) as f32,
                    (v2.z - v0.z) as f32,
                );
                let normal = edge1.cross(&edge2).normalize();
                
                // Write normal (3 floats)
                buffer.extend_from_slice(&(normal.x as f32).to_le_bytes());
                buffer.extend_from_slice(&(normal.y as f32).to_le_bytes());
                buffer.extend_from_slice(&(normal.z as f32).to_le_bytes());
                
                // Write vertices (9 floats)
                for vertex in &polygon.vertices {
                    buffer.extend_from_slice(&(vertex.pos.x as f32).to_le_bytes());
                    buffer.extend_from_slice(&(vertex.pos.y as f32).to_le_bytes());
                    buffer.extend_from_slice(&(vertex.pos.z as f32).to_le_bytes());
                }
                
                // Write attribute byte count (2 bytes, usually 0)
                buffer.extend_from_slice(&0u16.to_le_bytes());
            }
        }
        
        Ok(buffer)
    }
    
    /// Get the bounding box of the geometry
    pub fn bounding_box(&self) -> Result<(Point3<T>, Point3<T>), CsgError> {
        let aabb = self.csg.bounding_box();
        
        let min_point = Point3::new(
            T::from_f64(aabb.mins.x as f64).ok_or_else(|| CsgError::OperationFailed("Invalid min coordinate".to_string()))?,
            T::from_f64(aabb.mins.y as f64).ok_or_else(|| CsgError::OperationFailed("Invalid min coordinate".to_string()))?,
            T::from_f64(aabb.mins.z as f64).ok_or_else(|| CsgError::OperationFailed("Invalid min coordinate".to_string()))?,
        );
        
        let max_point = Point3::new(
            T::from_f64(aabb.maxs.x as f64).ok_or_else(|| CsgError::OperationFailed("Invalid max coordinate".to_string()))?,
            T::from_f64(aabb.maxs.y as f64).ok_or_else(|| CsgError::OperationFailed("Invalid max coordinate".to_string()))?,
            T::from_f64(aabb.maxs.z as f64).ok_or_else(|| CsgError::OperationFailed("Invalid max coordinate".to_string()))?,
        );
        
        Ok((min_point, max_point))
    }
    
    /// Get the number of vertices in the geometry
    pub fn vertex_count(&self) -> usize {
        self.csg.vertices().len()
    }
    
    /// Get the number of faces in the geometry
    pub fn face_count(&self) -> usize {
        self.csg.polygons.len()
    }
}

impl<T: RealField + FromPrimitive + ToPrimitive> Clone for CsgGeometry<T> {
    fn clone(&self) -> Self {
        Self {
            csg: self.csg.clone(),
            _phantom: std::marker::PhantomData,
        }
    }
}

/// CSG mesh adapter for backward compatibility
pub type CsgMeshAdapter<T> = CsgOperator<T>;

/// Builder pattern for complex CSG operations
pub struct CsgBuilder<T: RealField> {
    geometry: Option<CsgGeometry<T>>,
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + FromPrimitive + ToPrimitive> CsgBuilder<T> {
    /// Create a new CSG builder
    pub fn new() -> Self {
        Self {
            geometry: None,
            _phantom: std::marker::PhantomData,
        }
    }
    
    /// Start with a cube
    pub fn cube(mut self, width: T, height: T, depth: T) -> Result<Self, CsgError> {
        let operator = CsgOperator::new();
        self.geometry = Some(operator.create_cube(width, height, depth)?);
        Ok(self)
    }
    
    /// Start with a sphere
    pub fn sphere(mut self, radius: T, h_segments: usize, v_segments: usize) -> Result<Self, CsgError> {
        let operator = CsgOperator::new();
        self.geometry = Some(operator.create_sphere(radius, h_segments, v_segments)?);
        Ok(self)
    }
    
    /// Start with a cylinder
    pub fn cylinder(mut self, radius: T, height: T, segments: usize) -> Result<Self, CsgError> {
        let operator = CsgOperator::new();
        self.geometry = Some(operator.create_cylinder(radius, height, segments)?);
        Ok(self)
    }
    
    /// Add another geometry using union
    pub fn add(mut self, other: CsgGeometry<T>) -> Result<Self, CsgError> {
        if let Some(geometry) = self.geometry.take() {
            self.geometry = Some(geometry.union(&other));
        } else {
            self.geometry = Some(other);
        }
        Ok(self)
    }
    
    /// Subtract another geometry
    pub fn subtract(mut self, other: CsgGeometry<T>) -> Result<Self, CsgError> {
        if let Some(geometry) = self.geometry.take() {
            self.geometry = Some(geometry.difference(&other));
            Ok(self)
        } else {
            Err(CsgError::OperationFailed("No base geometry for subtraction".to_string()))
        }
    }
    
    /// Intersect with another geometry
    pub fn intersect(mut self, other: CsgGeometry<T>) -> Result<Self, CsgError> {
        if let Some(geometry) = self.geometry.take() {
            self.geometry = Some(geometry.intersection(&other));
            Ok(self)
        } else {
            Err(CsgError::OperationFailed("No base geometry for intersection".to_string()))
        }
    }
    
    /// Apply translation
    pub fn translate(mut self, translation: Vector3<T>) -> Result<Self, CsgError> {
        if let Some(mut geometry) = self.geometry.take() {
            geometry.translate(&translation)?;
            self.geometry = Some(geometry);
            Ok(self)
        } else {
            Err(CsgError::OperationFailed("No geometry to translate".to_string()))
        }
    }
    
    /// Apply rotation (in degrees)
    pub fn rotate(mut self, x_deg: T, y_deg: T, z_deg: T) -> Result<Self, CsgError> {
        if let Some(mut geometry) = self.geometry.take() {
            geometry.rotate(x_deg, y_deg, z_deg)?;
            self.geometry = Some(geometry);
            Ok(self)
        } else {
            Err(CsgError::OperationFailed("No geometry to rotate".to_string()))
        }
    }
    
    /// Apply scaling
    pub fn scale(mut self, scale_x: T, scale_y: T, scale_z: T) -> Result<Self, CsgError> {
        if let Some(mut geometry) = self.geometry.take() {
            geometry.scale(scale_x, scale_y, scale_z)?;
            self.geometry = Some(geometry);
            Ok(self)
        } else {
            Err(CsgError::OperationFailed("No geometry to scale".to_string()))
        }
    }
    
    /// Build the final geometry
    pub fn build(self) -> Result<CsgGeometry<T>, CsgError> {
        self.geometry.ok_or_else(|| CsgError::OperationFailed("No geometry built".to_string()))
    }
}

impl<T: RealField + FromPrimitive + ToPrimitive> Default for CsgBuilder<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_create_cube() {
        let operator = CsgOperator::<f64>::new();
        let cube = operator.create_cube(2.0, 2.0, 2.0).unwrap();
        
        assert!(cube.vertex_count() > 0);
        assert!(cube.face_count() > 0);
    }

    #[test]
    fn test_create_sphere() {
        let operator = CsgOperator::<f64>::new();
        let sphere = operator.create_sphere(1.0, 16, 8).unwrap();
        
        assert!(sphere.vertex_count() > 0);
        assert!(sphere.face_count() > 0);
    }

    #[test]
    fn test_boolean_operations() {
        let operator = CsgOperator::<f64>::new();
        let cube = operator.create_cube(2.0, 2.0, 2.0).unwrap();
        let sphere = operator.create_sphere(1.0, 16, 8).unwrap();
        
        // Test union
        let union_result = cube.union(&sphere);
        assert!(union_result.vertex_count() > 0);
        
        // Test difference
        let diff_result = cube.difference(&sphere);
        assert!(diff_result.vertex_count() > 0);
        
        // Test intersection
        let intersect_result = cube.intersection(&sphere);
        assert!(intersect_result.vertex_count() > 0);
    }

    #[test]
    fn test_transformations() {
        let operator = CsgOperator::<f64>::new();
        let mut cube = operator.create_cube(1.0, 1.0, 1.0).unwrap();
        
        // Test translation
        let translation = Vector3::new(1.0, 2.0, 3.0);
        cube.translate(&translation).unwrap();
        
        // Test rotation (in degrees)
        cube.rotate(45.0, 0.0, 0.0).unwrap();
        
        // Test scaling
        cube.scale(2.0, 2.0, 2.0).unwrap();
        
        assert!(cube.vertex_count() > 0);
    }

    #[test]
    fn test_builder_pattern() {
        let result = CsgBuilder::<f64>::new()
            .cube(2.0, 2.0, 2.0).unwrap()
            .translate(Vector3::new(1.0, 0.0, 0.0)).unwrap()
            .build().unwrap();
        
        assert!(result.vertex_count() > 0);
    }

    #[test]
    fn test_mesh_conversion() {
        let operator = CsgOperator::<f64>::new();
        let cube = operator.create_cube(1.0, 1.0, 1.0).unwrap();
        
        let mesh = cube.to_mesh().unwrap();
        assert!(!mesh.vertices.is_empty());
        assert!(!mesh.faces.is_empty());
    }

    #[test]
    fn test_stl_export() {
        let operator = CsgOperator::<f64>::new();
        let cube = operator.create_cube(1.0, 1.0, 1.0).unwrap();
        
        let stl_content = cube.to_stl("test_cube").unwrap();
        assert!(stl_content.contains("solid test_cube"));
        assert!(stl_content.contains("facet normal"));
        assert!(stl_content.contains("vertex"));
        assert!(stl_content.contains("endsolid"));
    }

    #[test]
    fn test_bounding_box() {
        let operator = CsgOperator::<f64>::new();
        let cube = operator.create_cube(2.0, 2.0, 2.0).unwrap();
        
        let (min_point, max_point) = cube.bounding_box().unwrap();
        
        // Cube should be centered at origin, so bounds should be symmetric
        assert_relative_eq!(min_point.x, -1.0, epsilon = 1e-6);
        assert_relative_eq!(max_point.x, 1.0, epsilon = 1e-6);
        assert_relative_eq!(min_point.y, -1.0, epsilon = 1e-6);
        assert_relative_eq!(max_point.y, 1.0, epsilon = 1e-6);
        assert_relative_eq!(min_point.z, -1.0, epsilon = 1e-6);
        assert_relative_eq!(max_point.z, 1.0, epsilon = 1e-6);
    }

    #[test]
    fn test_frustum_creation() {
        let operator = CsgOperator::<f64>::new();
        
        // Test frustum (truncated cone)
        let frustum = operator.create_frustum(2.0, 1.0, 3.0, 16).unwrap();
        assert!(frustum.vertex_count() > 0);
        assert!(frustum.face_count() > 0);
        
        // Test cone (frustum with top radius = 0)
        let cone = operator.create_cone(2.0, 3.0, 16).unwrap();
        assert!(cone.vertex_count() > 0);
        assert!(cone.face_count() > 0);
        
        // Verify cone has fewer vertices than frustum (no top circle)
        assert!(cone.vertex_count() <= frustum.vertex_count());
    }
}