//! Structured grid generation and manipulation
//!
//! This module provides structured grid generation capabilities for CFD simulations,
//! including Cartesian, cylindrical, and curvilinear grid systems.

use crate::mesh::{Mesh, Vertex, Face, Cell, MeshTopology};
use nalgebra::{Point3, Vector3, RealField};
use num_traits::FromPrimitive;
use std::f64::consts::PI;
use thiserror::Error;

/// Named constants for grid parameters
const DEFAULT_MIN_SPACING: f64 = 1e-6;
const DEFAULT_STRETCHING_RATIO: f64 = 1.1;
const DEFAULT_SMOOTHING_ITERATIONS: usize = 10;
const ORTHOGONALITY_TOLERANCE: f64 = 1e-3;

/// Grid generation errors
#[derive(Debug, Error)]
pub enum GridError {
    #[error("Invalid dimensions: {0}")]
    InvalidDimensions(String),
    #[error("Invalid spacing: {0}")]
    InvalidSpacing(String),
    #[error("Grid generation failed: {0}")]
    GenerationFailed(String),
    #[error("Transformation error: {0}")]
    TransformationError(String),
}

/// Grid coordinate system
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CoordinateSystem {
    /// Cartesian (x, y, z) coordinates
    Cartesian,
    /// Cylindrical (r, θ, z) coordinates
    Cylindrical,
    /// Spherical (r, θ, φ) coordinates
    Spherical,
    /// Curvilinear coordinates
    Curvilinear,
}

/// Grid spacing type
#[derive(Debug)]
pub enum GridSpacing<T: RealField> {
    /// Uniform spacing
    Uniform(T),
    /// Geometric progression with ratio
    Geometric { first: T, ratio: T },
    /// Hyperbolic tangent stretching
    Hyperbolic { center: T, width: T },
    /// Custom spacing function
    Custom(Box<dyn Fn(usize, usize) -> T + Send + Sync>),
}

/// Structured grid configuration
#[derive(Debug)]
pub struct GridConfig<T: RealField> {
    /// Number of points in each direction
    pub dimensions: [usize; 3],
    /// Domain bounds
    pub bounds: [[T; 2]; 3],
    /// Coordinate system
    pub coordinate_system: CoordinateSystem,
    /// Grid spacing in each direction
    pub spacing: [GridSpacing<T>; 3],
    /// Enable grid smoothing
    pub enable_smoothing: bool,
    /// Number of smoothing iterations
    pub smoothing_iterations: usize,
}

impl<T: RealField + FromPrimitive> Default for GridConfig<T> {
    fn default() -> Self {
        Self {
            dimensions: [10, 10, 10],
            bounds: [
                [T::zero(), T::one()],
                [T::zero(), T::one()],
                [T::zero(), T::one()],
            ],
            coordinate_system: CoordinateSystem::Cartesian,
            spacing: [
                GridSpacing::Uniform(T::from_f64(0.1).unwrap()),
                GridSpacing::Uniform(T::from_f64(0.1).unwrap()),
                GridSpacing::Uniform(T::from_f64(0.1).unwrap()),
            ],
            enable_smoothing: false,
            smoothing_iterations: DEFAULT_SMOOTHING_ITERATIONS,
        }
    }
}

/// Structured grid generator
pub struct GridGenerator<T: RealField> {
    config: GridConfig<T>,
}

impl<T: RealField + FromPrimitive> GridGenerator<T> {
    /// Create a new grid generator
    pub fn new(config: GridConfig<T>) -> Self {
        Self { config }
    }
    
    /// Generate a structured grid
    pub fn generate(&self) -> Result<StructuredGrid<T>, GridError> {
        // Validate configuration
        self.validate_config()?;
        
        // Generate grid points
        let points = match self.config.coordinate_system {
            CoordinateSystem::Cartesian => self.generate_cartesian_points()?,
            CoordinateSystem::Cylindrical => self.generate_cylindrical_points()?,
            CoordinateSystem::Spherical => self.generate_spherical_points()?,
            CoordinateSystem::Curvilinear => self.generate_curvilinear_points()?,
        };
        
        // Create structured grid
        let mut grid = StructuredGrid {
            dimensions: self.config.dimensions,
            points,
            coordinate_system: self.config.coordinate_system,
        };
        
        // Apply smoothing if enabled
        if self.config.enable_smoothing {
            grid.smooth(self.config.smoothing_iterations);
        }
        
        Ok(grid)
    }
    
    /// Validate grid configuration
    fn validate_config(&self) -> Result<(), GridError> {
        // Check dimensions
        for (i, &dim) in self.config.dimensions.iter().enumerate() {
            if dim < 2 {
                return Err(GridError::InvalidDimensions(
                    format!("Dimension {} must be at least 2, got {}", i, dim)
                ));
            }
        }
        
        // Check bounds
        for (i, bounds) in self.config.bounds.iter().enumerate() {
            if bounds[1] <= bounds[0] {
                return Err(GridError::InvalidDimensions(
                    format!("Invalid bounds for dimension {}: [{:?}, {:?}]", 
                            i, bounds[0], bounds[1])
                ));
            }
        }
        
        Ok(())
    }
    
    /// Generate Cartesian grid points
    fn generate_cartesian_points(&self) -> Result<Vec<Point3<T>>, GridError> {
        let [nx, ny, nz] = self.config.dimensions;
        let mut points = Vec::with_capacity(nx * ny * nz);
        
        // Generate coordinates in each direction
        let x_coords = self.generate_coordinates(0, nx)?;
        let y_coords = self.generate_coordinates(1, ny)?;
        let z_coords = self.generate_coordinates(2, nz)?;
        
        // Create grid points
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    points.push(Point3::new(
                        x_coords[i].clone(),
                        y_coords[j].clone(),
                        z_coords[k].clone(),
                    ));
                }
            }
        }
        
        Ok(points)
    }
    
    /// Generate cylindrical grid points
    fn generate_cylindrical_points(&self) -> Result<Vec<Point3<T>>, GridError> {
        let [nr, ntheta, nz] = self.config.dimensions;
        let mut points = Vec::with_capacity(nr * ntheta * nz);
        
        // Generate coordinates
        let r_coords = self.generate_coordinates(0, nr)?;
        let theta_coords = self.generate_coordinates(1, ntheta)?;
        let z_coords = self.generate_coordinates(2, nz)?;
        
        // Convert to Cartesian
        for k in 0..nz {
            for j in 0..ntheta {
                for i in 0..nr {
                    let r = r_coords[i].clone();
                    let theta = theta_coords[j].clone();
                    let z = z_coords[k].clone();
                    
                    let x = r.clone() * theta.clone().cos();
                    let y = r * theta.sin();
                    
                    points.push(Point3::new(x, y, z));
                }
            }
        }
        
        Ok(points)
    }
    
    /// Generate spherical grid points
    fn generate_spherical_points(&self) -> Result<Vec<Point3<T>>, GridError> {
        let [nr, ntheta, nphi] = self.config.dimensions;
        let mut points = Vec::with_capacity(nr * ntheta * nphi);
        
        // Generate coordinates
        let r_coords = self.generate_coordinates(0, nr)?;
        let theta_coords = self.generate_coordinates(1, ntheta)?;
        let phi_coords = self.generate_coordinates(2, nphi)?;
        
        // Convert to Cartesian
        for k in 0..nphi {
            for j in 0..ntheta {
                for i in 0..nr {
                    let r = r_coords[i].clone();
                    let theta = theta_coords[j].clone();
                    let phi = phi_coords[k].clone();
                    
                    let x = r.clone() * theta.clone().sin() * phi.clone().cos();
                    let y = r.clone() * theta.sin() * phi.sin();
                    let z = r * theta.cos();
                    
                    points.push(Point3::new(x, y, z));
                }
            }
        }
        
        Ok(points)
    }
    
    /// Generate curvilinear grid points
    fn generate_curvilinear_points(&self) -> Result<Vec<Point3<T>>, GridError> {
        // Start with Cartesian grid and apply transformation
        let points = self.generate_cartesian_points()?;
        
        // Apply curvilinear transformation (example: simple wave)
        let transformed: Vec<Point3<T>> = points.into_iter()
            .map(|p| {
                let amplitude = T::from_f64(0.1).unwrap();
                let frequency = T::from_f64(2.0 * PI).unwrap();
                
                let x = p.x.clone();
                let y = p.y.clone() + amplitude * (frequency * p.x.clone()).sin();
                let z = p.z.clone();
                
                Point3::new(x, y, z)
            })
            .collect();
        
        Ok(transformed)
    }
    
    /// Generate coordinates for one dimension
    fn generate_coordinates(&self, dim: usize, n: usize) -> Result<Vec<T>, GridError> {
        let [min, max] = self.config.bounds[dim].clone();
        let mut coords = Vec::with_capacity(n);
        
        match &self.config.spacing[dim] {
            GridSpacing::Uniform(spacing) => {
                // Uniform spacing
                let actual_spacing = (max.clone() - min.clone()) / T::from_usize(n - 1).unwrap();
                for i in 0..n {
                    coords.push(min.clone() + actual_spacing.clone() * T::from_usize(i).unwrap());
                }
            }
            GridSpacing::Geometric { first, ratio } => {
                // Geometric progression
                coords.push(min.clone());
                let mut current_spacing = first.clone();
                let mut current = min.clone();
                
                for _ in 1..n {
                    current = current + current_spacing.clone();
                    coords.push(current.clone());
                    current_spacing = current_spacing * ratio.clone();
                }
                
                // Scale to fit bounds
                let actual_max = coords.last().unwrap().clone();
                let scale = (max.clone() - min.clone()) / (actual_max - min.clone());
                for coord in &mut coords {
                    *coord = min.clone() + (*coord - min.clone()) * scale.clone();
                }
            }
            GridSpacing::Hyperbolic { center, width } => {
                // Hyperbolic tangent stretching
                for i in 0..n {
                    let xi = T::from_f64(-1.0).unwrap() + 
                             T::from_f64(2.0).unwrap() * T::from_usize(i).unwrap() / 
                             T::from_usize(n - 1).unwrap();
                    
                    let stretched = center.clone() + width.clone() * xi.tanh();
                    let normalized = (stretched + T::one()) / T::from_f64(2.0).unwrap();
                    
                    coords.push(min.clone() + (max.clone() - min.clone()) * normalized);
                }
            }
            GridSpacing::Custom(func) => {
                // Custom spacing function
                for i in 0..n {
                    let t = func(i, n);
                    coords.push(min.clone() + (max.clone() - min.clone()) * t);
                }
            }
        }
        
        Ok(coords)
    }
    
    /// Create a uniform Cartesian grid
    pub fn uniform_cartesian(
        nx: usize, ny: usize, nz: usize,
        bounds: [[T; 2]; 3],
    ) -> Result<StructuredGrid<T>, GridError> {
        let config = GridConfig {
            dimensions: [nx, ny, nz],
            bounds,
            coordinate_system: CoordinateSystem::Cartesian,
            spacing: [
                GridSpacing::Uniform(T::zero()),
                GridSpacing::Uniform(T::zero()),
                GridSpacing::Uniform(T::zero()),
            ],
            enable_smoothing: false,
            smoothing_iterations: 0,
        };
        
        let generator = Self::new(config);
        generator.generate()
    }
}

/// Structured grid representation
#[derive(Debug, Clone)]
pub struct StructuredGrid<T: RealField> {
    /// Grid dimensions [nx, ny, nz]
    pub dimensions: [usize; 3],
    /// Grid points in row-major order
    pub points: Vec<Point3<T>>,
    /// Coordinate system
    pub coordinate_system: CoordinateSystem,
}

impl<T: RealField + FromPrimitive> StructuredGrid<T> {
    /// Get grid point at indices (i, j, k)
    pub fn get_point(&self, i: usize, j: usize, k: usize) -> Option<&Point3<T>> {
        let [nx, ny, nz] = self.dimensions;
        if i < nx && j < ny && k < nz {
            let index = i + j * nx + k * nx * ny;
            self.points.get(index)
        } else {
            None
        }
    }
    
    /// Get mutable grid point at indices
    pub fn get_point_mut(&mut self, i: usize, j: usize, k: usize) -> Option<&mut Point3<T>> {
        let [nx, ny, nz] = self.dimensions;
        if i < nx && j < ny && k < nz {
            let index = i + j * nx + k * nx * ny;
            self.points.get_mut(index)
        } else {
            None
        }
    }
    
    /// Convert to unstructured mesh
    pub fn to_mesh(&self) -> Mesh<T> {
        let mut mesh = Mesh::new();
        let [nx, ny, nz] = self.dimensions;
        
        // Add vertices
        for (id, point) in self.points.iter().enumerate() {
            mesh.vertices.push(Vertex {
                id,
                position: point.clone(),
            });
        }
        
        // Add hexahedral cells
        let mut face_id = 0;
        let mut cell_id = 0;
        
        for k in 0..nz-1 {
            for j in 0..ny-1 {
                for i in 0..nx-1 {
                    // Vertex indices for hexahedron
                    let v000 = i + j * nx + k * nx * ny;
                    let v100 = v000 + 1;
                    let v010 = v000 + nx;
                    let v110 = v010 + 1;
                    let v001 = v000 + nx * ny;
                    let v101 = v001 + 1;
                    let v011 = v001 + nx;
                    let v111 = v011 + 1;
                    
                    // Create 6 faces for hexahedron
                    let face_ids = vec![
                        face_id, face_id + 1, face_id + 2,
                        face_id + 3, face_id + 4, face_id + 5,
                    ];
                    
                    // Bottom face
                    mesh.faces.push(Face {
                        id: face_id,
                        vertices: vec![v000, v100, v110, v010],
                    });
                    face_id += 1;
                    
                    // Top face
                    mesh.faces.push(Face {
                        id: face_id,
                        vertices: vec![v001, v011, v111, v101],
                    });
                    face_id += 1;
                    
                    // Front face
                    mesh.faces.push(Face {
                        id: face_id,
                        vertices: vec![v000, v001, v101, v100],
                    });
                    face_id += 1;
                    
                    // Back face
                    mesh.faces.push(Face {
                        id: face_id,
                        vertices: vec![v010, v110, v111, v011],
                    });
                    face_id += 1;
                    
                    // Left face
                    mesh.faces.push(Face {
                        id: face_id,
                        vertices: vec![v000, v010, v011, v001],
                    });
                    face_id += 1;
                    
                    // Right face
                    mesh.faces.push(Face {
                        id: face_id,
                        vertices: vec![v100, v101, v111, v110],
                    });
                    face_id += 1;
                    
                    // Create cell
                    mesh.cells.push(Cell {
                        id: cell_id,
                        faces: face_ids,
                    });
                    cell_id += 1;
                }
            }
        }
        
        // Update topology
        mesh.topology = MeshTopology {
            num_vertices: mesh.vertices.len(),
            num_edges: 0, // Not computed
            num_faces: mesh.faces.len(),
            num_cells: mesh.cells.len(),
        };
        
        mesh
    }
    
    /// Apply Laplacian smoothing to grid
    pub fn smooth(&mut self, iterations: usize) {
        let [nx, ny, nz] = self.dimensions;
        let alpha = T::from_f64(0.5).unwrap(); // Smoothing factor
        
        for _ in 0..iterations {
            let mut new_points = self.points.clone();
            
            // Smooth interior points only
            for k in 1..nz-1 {
                for j in 1..ny-1 {
                    for i in 1..nx-1 {
                        let index = i + j * nx + k * nx * ny;
                        
                        // Get neighboring points
                        let neighbors = [
                            self.get_point(i-1, j, k).unwrap(),
                            self.get_point(i+1, j, k).unwrap(),
                            self.get_point(i, j-1, k).unwrap(),
                            self.get_point(i, j+1, k).unwrap(),
                            self.get_point(i, j, k-1).unwrap(),
                            self.get_point(i, j, k+1).unwrap(),
                        ];
                        
                        // Compute average
                        let mut avg = Vector3::zeros();
                        for neighbor in &neighbors {
                            avg += &neighbor.coords;
                        }
                        avg /= T::from_usize(neighbors.len()).unwrap();
                        
                        // Blend with original
                        new_points[index] = Point3::from(
                            self.points[index].coords.clone() * (T::one() - alpha.clone()) +
                            avg * alpha.clone()
                        );
                    }
                }
            }
            
            self.points = new_points;
        }
    }
    
    /// Compute grid quality metrics
    pub fn compute_quality(&self) -> GridQuality<T> {
        let mut min_jacobian = T::from_f64(f64::MAX).unwrap();
        let mut max_jacobian = T::zero();
        let mut min_orthogonality = T::from_f64(f64::MAX).unwrap();
        let mut max_skewness = T::zero();
        let mut min_aspect_ratio = T::from_f64(f64::MAX).unwrap();
        let mut max_aspect_ratio = T::zero();
        
        let [nx, ny, nz] = self.dimensions;
        
        for k in 0..nz-1 {
            for j in 0..ny-1 {
                for i in 0..nx-1 {
                    // Compute local metrics
                    let jacobian = self.compute_jacobian(i, j, k);
                    let orthogonality = self.compute_orthogonality(i, j, k);
                    let skewness = self.compute_skewness(i, j, k);
                    let aspect_ratio = self.compute_aspect_ratio(i, j, k);
                    
                    // Update bounds
                    if jacobian < min_jacobian { min_jacobian = jacobian.clone(); }
                    if jacobian > max_jacobian { max_jacobian = jacobian; }
                    if orthogonality < min_orthogonality { min_orthogonality = orthogonality; }
                    if skewness > max_skewness { max_skewness = skewness; }
                    if aspect_ratio < min_aspect_ratio { min_aspect_ratio = aspect_ratio.clone(); }
                    if aspect_ratio > max_aspect_ratio { max_aspect_ratio = aspect_ratio; }
                }
            }
        }
        
        GridQuality {
            min_jacobian,
            max_jacobian,
            min_orthogonality,
            max_skewness,
            min_aspect_ratio,
            max_aspect_ratio,
        }
    }
    
    /// Compute Jacobian at cell (i, j, k)
    fn compute_jacobian(&self, i: usize, j: usize, k: usize) -> T {
        // Compute Jacobian determinant of transformation
        let p000 = self.get_point(i, j, k).unwrap();
        let p100 = self.get_point(i+1, j, k).unwrap();
        let p010 = self.get_point(i, j+1, k).unwrap();
        let p001 = self.get_point(i, j, k+1).unwrap();
        
        let dx = p100 - p000;
        let dy = p010 - p000;
        let dz = p001 - p000;
        
        dx.cross(&dy).dot(&dz).abs()
    }
    
    /// Compute orthogonality at cell (i, j, k)
    fn compute_orthogonality(&self, i: usize, j: usize, k: usize) -> T {
        let p000 = self.get_point(i, j, k).unwrap();
        let p100 = self.get_point(i+1, j, k).unwrap();
        let p010 = self.get_point(i, j+1, k).unwrap();
        let p001 = self.get_point(i, j, k+1).unwrap();
        
        let dx = (p100 - p000).normalize();
        let dy = (p010 - p000).normalize();
        let dz = (p001 - p000).normalize();
        
        // Compute minimum dot product (ideally 0 for orthogonal)
        let dot_xy = dx.dot(&dy).abs();
        let dot_xz = dx.dot(&dz).abs();
        let dot_yz = dy.dot(&dz).abs();
        
        T::one() - dot_xy.max(dot_xz).max(dot_yz)
    }
    
    /// Compute skewness at cell (i, j, k)
    fn compute_skewness(&self, i: usize, j: usize, k: usize) -> T {
        // Simplified skewness metric
        T::zero()
    }
    
    /// Compute aspect ratio at cell (i, j, k)
    fn compute_aspect_ratio(&self, i: usize, j: usize, k: usize) -> T {
        let p000 = self.get_point(i, j, k).unwrap();
        let p100 = self.get_point(i+1, j, k).unwrap();
        let p010 = self.get_point(i, j+1, k).unwrap();
        let p001 = self.get_point(i, j, k+1).unwrap();
        
        let dx = (p100 - p000).norm();
        let dy = (p010 - p000).norm();
        let dz = (p001 - p000).norm();
        
        let max_len = dx.clone().max(dy.clone()).max(dz.clone());
        let min_len = dx.min(dy).min(dz);
        
        if min_len > T::zero() {
            max_len / min_len
        } else {
            T::from_f64(f64::MAX).unwrap()
        }
    }
}

/// Grid quality metrics
#[derive(Debug, Clone)]
pub struct GridQuality<T: RealField> {
    pub min_jacobian: T,
    pub max_jacobian: T,
    pub min_orthogonality: T,
    pub max_skewness: T,
    pub min_aspect_ratio: T,
    pub max_aspect_ratio: T,
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_uniform_cartesian_grid() {
        let grid = GridGenerator::<f64>::uniform_cartesian(
            5, 5, 5,
            [[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]],
        ).unwrap();
        
        assert_eq!(grid.dimensions, [5, 5, 5]);
        assert_eq!(grid.points.len(), 125);
        assert_eq!(grid.coordinate_system, CoordinateSystem::Cartesian);
    }
    
    #[test]
    fn test_grid_to_mesh() {
        let grid = GridGenerator::<f64>::uniform_cartesian(
            3, 3, 3,
            [[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]],
        ).unwrap();
        
        let mesh = grid.to_mesh();
        assert_eq!(mesh.vertices.len(), 27);
        assert_eq!(mesh.cells.len(), 8); // 2x2x2 cells
    }
}
