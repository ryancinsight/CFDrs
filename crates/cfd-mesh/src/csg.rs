//! CSG (Constructive Solid Geometry) operations using the csgrs library
//!
//! This module provides comprehensive CSG mesh generation and boolean operations
//! using the csgrs crate, enabling creation of complex geometries through
//! constructive solid geometry techniques.
//! 
//! Features:
//! - Full CSG boolean operations (union, intersection, difference, XOR)
//! - Primitive geometry generation (spheres, boxes, cylinders, etc.)
//! - Mesh transformations (translate, rotate, scale, mirror)
//! - STL export functionality
//! - Integration with CFD mesh structures

use crate::mesh::{Mesh, Vertex, Face};
use nalgebra::{Vector3, Point3, RealField};
use num_traits::{FromPrimitive, ToPrimitive};
use thiserror::Error;
use std::fmt::Debug;
use cfd_core::constants;
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
        
        let csg = csgrs::Cube::new(w, h, d);
        Ok(CsgGeometry::new(Box::new(csg)))
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
        
        let csg = csgrs::Sphere::new(r, horizontal_segments, vertical_segments);
        Ok(CsgGeometry::new(Box::new(csg)))
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
        
        let csg = csgrs::Cylinder::new(r, h, segments);
        Ok(CsgGeometry::new(Box::new(csg)))
    }
    
    /// Create a cone with specified base radius, top radius, height, and resolution
    pub fn create_cone(&self, base_radius: T, top_radius: T, height: T, segments: usize) -> Result<CsgGeometry<T>, CsgError> {
        let r1 = base_radius.to_f64().ok_or_else(|| CsgError::InvalidParameters("Invalid base radius".to_string()))?;
        let r2 = top_radius.to_f64().ok_or_else(|| CsgError::InvalidParameters("Invalid top radius".to_string()))?;
        let h = height.to_f64().ok_or_else(|| CsgError::InvalidParameters("Invalid height".to_string()))?;
        
        let zero = T::from_f64(0.0).unwrap_or_else(|| T::zero());
        if base_radius < zero || top_radius < zero || height <= zero {
            return Err(CsgError::InvalidParameters("Radii must be non-negative and height positive".to_string()));
        }
        
        if segments < 3 {
            return Err(CsgError::InvalidParameters("Insufficient segments for cone".to_string()));
        }
        
        let csg = csgrs::Cone::new(r1, r2, h, segments);
        Ok(CsgGeometry::new(Box::new(csg)))
    }
}

impl<T: RealField + FromPrimitive + ToPrimitive> Default for CsgOperator<T> {
    fn default() -> Self {
        Self::new()
    }
}

/// Wrapper for CSG geometry with type safety and CFD integration
pub struct CsgGeometry<T: RealField> {
    csg: Box<dyn CSG>,
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + FromPrimitive + ToPrimitive> CsgGeometry<T> {
    /// Create a new CSG geometry wrapper
    pub fn new(csg: Box<dyn CSG>) -> Self {
        Self {
            csg,
            _phantom: std::marker::PhantomData,
        }
    }
    
    /// Get the underlying CSG object
    pub fn inner(&self) -> &dyn CSG {
        &*self.csg
    }
    
    /// Perform union operation with another geometry
    pub fn union(&self, other: &CsgGeometry<T>) -> CsgGeometry<T> {
        let result = self.csg.union(&*other.csg);
        CsgGeometry::new(result)
    }
    
    /// Perform difference operation (subtract other from self)
    pub fn difference(&self, other: &CsgGeometry<T>) -> CsgGeometry<T> {
        let result = self.csg.difference(&*other.csg);
        CsgGeometry::new(result)
    }
    
    /// Perform intersection operation with another geometry
    pub fn intersection(&self, other: &CsgGeometry<T>) -> CsgGeometry<T> {
        let result = self.csg.intersection(&*other.csg);
        CsgGeometry::new(result)
    }
    
    /// Perform XOR operation with another geometry
    pub fn xor(&self, other: &CsgGeometry<T>) -> CsgGeometry<T> {
        let result = self.csg.xor(&*other.csg);
        CsgGeometry::new(result)
    }
    
    /// Translate the geometry by the given vector
    pub fn translate(&mut self, translation: &Vector3<T>) -> Result<(), CsgError> {
        let tx = translation.x.to_f64().ok_or_else(|| CsgError::InvalidParameters("Invalid translation x".to_string()))?;
        let ty = translation.y.to_f64().ok_or_else(|| CsgError::InvalidParameters("Invalid translation y".to_string()))?;
        let tz = translation.z.to_f64().ok_or_else(|| CsgError::InvalidParameters("Invalid translation z".to_string()))?;
        
        self.csg.translate(tx, ty, tz);
        Ok(())
    }
    
    /// Rotate the geometry around the given axis by the given angle (in radians)
    pub fn rotate(&mut self, axis: &Vector3<T>, angle: T) -> Result<(), CsgError> {
        let ax = axis.x.to_f64().ok_or_else(|| CsgError::InvalidParameters("Invalid axis x".to_string()))?;
        let ay = axis.y.to_f64().ok_or_else(|| CsgError::InvalidParameters("Invalid axis y".to_string()))?;
        let az = axis.z.to_f64().ok_or_else(|| CsgError::InvalidParameters("Invalid axis z".to_string()))?;
        let angle_rad = angle.to_f64().ok_or_else(|| CsgError::InvalidParameters("Invalid angle".to_string()))?;
        
        self.csg.rotate(ax, ay, az, angle_rad);
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
        
        self.csg.scale(sx, sy, sz);
        Ok(())
    }
    
    /// Mirror the geometry across a plane defined by a normal vector
    pub fn mirror(&mut self, normal: &Vector3<T>) -> Result<(), CsgError> {
        let nx = normal.x.to_f64().ok_or_else(|| CsgError::InvalidParameters("Invalid normal x".to_string()))?;
        let ny = normal.y.to_f64().ok_or_else(|| CsgError::InvalidParameters("Invalid normal y".to_string()))?;
        let nz = normal.z.to_f64().ok_or_else(|| CsgError::InvalidParameters("Invalid normal z".to_string()))?;
        
        // Normalize the normal vector
        let length = (nx * nx + ny * ny + nz * nz).sqrt();
        let small_number = T::from_f64(constants::SMALL_NUMBER).unwrap_or_else(|| T::from_f64(1e-12).unwrap());
        if T::from_f64(length).unwrap() < small_number {
            return Err(CsgError::InvalidParameters("Normal vector too small".to_string()));
        }
        
        self.csg.mirror(nx / length, ny / length, nz / length);
        Ok(())
    }
    
    /// Convert to CFD mesh format
    pub fn to_mesh(&self) -> Result<Mesh<T>, CsgError> {
        let polygons = self.csg.to_polygons();
        let mut vertices = Vec::new();
        let mut faces = Vec::new();
        let mut vertex_map = std::collections::HashMap::new();
        let mut vertex_id = 1;
        
        for polygon in polygons {
            let mut face_vertices = Vec::new();
            
            for vertex in polygon.vertices {
                let pos = vertex.pos;
                let point = Point3::new(
                    T::from_f64(pos.x).ok_or_else(|| CsgError::OperationFailed("Invalid vertex x coordinate".to_string()))?,
                    T::from_f64(pos.y).ok_or_else(|| CsgError::OperationFailed("Invalid vertex y coordinate".to_string()))?,
                    T::from_f64(pos.z).ok_or_else(|| CsgError::OperationFailed("Invalid vertex z coordinate".to_string()))?,
                );
                
                // Use a key based on position with small tolerance for deduplication
                let scale = 1e6;
                let key = (
                    (pos.x * scale).round() as i64,
                    (pos.y * scale).round() as i64,
                    (pos.z * scale).round() as i64,
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
        let mut mesh = Mesh::new();
        mesh.vertices = vertices;
        mesh.faces = faces;
        mesh.update_topology();
        
        Ok(mesh)
    }
    
    /// Export to STL format
    pub fn to_stl(&self, name: &str) -> Result<String, CsgError> {
        Ok(self.csg.to_stl_ascii(name))
    }
    
    /// Export to binary STL format
    pub fn to_stl_binary(&self, name: &str) -> Result<Vec<u8>, CsgError> {
        Ok(self.csg.to_stl_binary(name))
    }
    
    /// Get the bounding box of the geometry
    pub fn bounding_box(&self) -> Result<(Point3<T>, Point3<T>), CsgError> {
        let polygons = self.csg.to_polygons();
        
        if polygons.is_empty() {
            return Err(CsgError::OperationFailed("Empty geometry".to_string()));
        }
        
        let mut min_x = f64::INFINITY;
        let mut min_y = f64::INFINITY;
        let mut min_z = f64::INFINITY;
        let mut max_x = f64::NEG_INFINITY;
        let mut max_y = f64::NEG_INFINITY;
        let mut max_z = f64::NEG_INFINITY;
        
        for polygon in polygons {
            for vertex in polygon.vertices {
                let pos = vertex.pos;
                min_x = min_x.min(pos.x);
                min_y = min_y.min(pos.y);
                min_z = min_z.min(pos.z);
                max_x = max_x.max(pos.x);
                max_y = max_y.max(pos.y);
                max_z = max_z.max(pos.z);
            }
        }
        
        let min_point = Point3::new(
            T::from_f64(min_x).ok_or_else(|| CsgError::OperationFailed("Invalid min coordinate".to_string()))?,
            T::from_f64(min_y).ok_or_else(|| CsgError::OperationFailed("Invalid min coordinate".to_string()))?,
            T::from_f64(min_z).ok_or_else(|| CsgError::OperationFailed("Invalid min coordinate".to_string()))?,
        );
        
        let max_point = Point3::new(
            T::from_f64(max_x).ok_or_else(|| CsgError::OperationFailed("Invalid max coordinate".to_string()))?,
            T::from_f64(max_y).ok_or_else(|| CsgError::OperationFailed("Invalid max coordinate".to_string()))?,
            T::from_f64(max_z).ok_or_else(|| CsgError::OperationFailed("Invalid max coordinate".to_string()))?,
        );
        
        Ok((min_point, max_point))
    }
    
    /// Get the number of vertices in the geometry
    pub fn vertex_count(&self) -> usize {
        self.csg.to_polygons().iter()
            .map(|p| p.vertices.len())
            .sum()
    }
    
    /// Get the number of faces in the geometry
    pub fn face_count(&self) -> usize {
        self.csg.to_polygons().len()
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
    
    /// Apply rotation
    pub fn rotate(mut self, axis: Vector3<T>, angle: T) -> Result<Self, CsgError> {
        if let Some(mut geometry) = self.geometry.take() {
            geometry.rotate(&axis, angle)?;
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
        
        // Test rotation
        let axis = Vector3::new(0.0, 0.0, 1.0);
        cube.rotate(&axis, std::f64::consts::PI / 4.0).unwrap();
        
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
        assert!(!mesh.vertices().is_empty());
        assert!(!mesh.faces().is_empty());
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
}