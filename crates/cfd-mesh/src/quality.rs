//! Mesh quality analysis and metrics.
//!
//! This module provides comprehensive tools for analyzing mesh quality using various 
//! metrics such as aspect ratio, skewness, and orthogonality. It supports all standard
//! element types and uses efficient single-pass algorithms for large meshes.
//!
//! ## Performance Optimizations
//!
//! - Single-pass analysis using running statistics to avoid memory allocation
//! - Correct vertex collection using deduplication to eliminate false duplicates
//! - Element-type-specific calculations for accurate metric computation
//! - Robust numerical methods avoiding floating-point conversion errors

use crate::mesh::{Mesh, Cell, ElementType};
use nalgebra::{Point3, RealField};
use num_traits::Float;
// Avoid ComplexField; use real-valued Float methods consistently
use serde::{Deserialize, Serialize};
use std::collections::HashSet;

/// Mesh quality metrics with comprehensive statistics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityMetrics<T: RealField> {
    /// Aspect ratio statistics
    pub aspect_ratio: QualityStatistics<T>,
    /// Skewness statistics
    pub skewness: QualityStatistics<T>,
    /// Orthogonality statistics  
    pub orthogonality: QualityStatistics<T>,
    /// Volume statistics
    pub volume: QualityStatistics<T>,
    /// Number of cells analyzed
    pub num_cells: usize,
    /// Number of cells that failed analysis
    pub failed_cells: usize,
}

/// Statistical summary computed using running algorithms for memory efficiency
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityStatistics<T: RealField> {
    /// Minimum value
    pub min: T,
    /// Maximum value
    pub max: T,
    /// Mean value
    pub mean: T,
    /// Standard deviation (computed using Welford's algorithm)
    pub std_dev: T,
    /// Number of samples
    pub count: usize,
}

/// Running statistics calculator using Welford's online algorithm
/// This allows single-pass computation without storing all values
#[derive(Debug, Clone)]
struct RunningStats<T: RealField> {
    count: usize,
    mean: T,
    m2: T,  // Sum of squares of differences from mean
    min: T,
    max: T,
}

impl<T: RealField + Float> RunningStats<T> {
    /// Create new running statistics tracker
    fn new() -> Self {
        Self {
            count: 0,
            mean: T::zero(),
            m2: T::zero(),
            min: T::infinity(),
            max: T::neg_infinity(),
        }
    }

    /// Update statistics with a new value using Welford's method
    fn update(&mut self, value: T) {
        self.count += 1;
        
        // Update min/max
        if value < self.min {
            self.min = value;
        }
        if value > self.max {
            self.max = value;
        }

        // Welford's algorithm for online mean and variance
        let delta = value - self.mean;
        self.mean = self.mean + delta / T::from(self.count).unwrap_or_else(|_| T::zero());
        let delta2 = value - self.mean;
        self.m2 = self.m2 + delta * delta2;
    }

    /// Finalize statistics computation
    fn finalize(self) -> QualityStatistics<T> {
        if self.count == 0 {
            return QualityStatistics {
                min: T::zero(),
                max: T::zero(),
                mean: T::zero(),
                std_dev: T::zero(),
                count: 0,
            };
        }

        let variance = if self.count > 1 {
            self.m2 / T::from(self.count - 1).unwrap_or_else(|_| T::zero())
        } else {
            T::zero()
        };

        QualityStatistics {
            min: self.min,
            max: self.max,
            mean: self.mean,
            std_dev: num_traits::Float::sqrt(variance),
            count: self.count,
        }
    }
}

/// Comprehensive mesh quality analyzer with element-type-specific calculations
#[derive(Debug)]
pub struct QualityAnalyzer<T: RealField> {
    /// Minimum acceptable aspect ratio
    pub min_aspect_ratio: T,
    /// Maximum acceptable skewness
    pub max_skewness: T,
    /// Minimum acceptable orthogonality
    pub min_orthogonality: T,
}

impl<T: RealField + Float> Default for QualityAnalyzer<T> {
    fn default() -> Self {
        Self {
            min_aspect_ratio: T::from(0.1).unwrap_or_else(|_| T::zero()),
            max_skewness: T::from(0.8).unwrap_or_else(|_| T::zero()),
            min_orthogonality: T::from(0.2).unwrap_or_else(|_| T::zero()),
        }
    }
}

impl<T: RealField + Float> QualityAnalyzer<T> {
    /// Create a new quality analyzer with custom thresholds
    pub fn new(min_aspect_ratio: T, max_skewness: T, min_orthogonality: T) -> Self {
        Self {
            min_aspect_ratio,
            max_skewness,
            min_orthogonality,
        }
    }

    /// Analyze mesh quality using single-pass algorithm for optimal performance
    /// 
    /// This method eliminates the multi-pass inefficiency by computing all statistics
    /// in a single iteration using running statistics algorithms.
    pub fn analyze(&self, mesh: &Mesh<T>) -> QualityMetrics<T> {
        // Initialize running statistics for each metric
        let mut aspect_ratio_stats = RunningStats::new();
        let mut skewness_stats = RunningStats::new();
        let mut orthogonality_stats = RunningStats::new();
        let mut volume_stats = RunningStats::new();
        
        let mut failed_cells = 0;

        // Single pass through all cells
        for cell in &mesh.cells {
            if let Some((ar, sk, orth, vol)) = self.analyze_cell(cell, mesh) {
                aspect_ratio_stats.update(ar);
                skewness_stats.update(sk);
                orthogonality_stats.update(orth);
                volume_stats.update(vol);
            } else {
                failed_cells += 1;
            }
        }

        QualityMetrics {
            aspect_ratio: aspect_ratio_stats.finalize(),
            skewness: skewness_stats.finalize(),
            orthogonality: orthogonality_stats.finalize(),
            volume: volume_stats.finalize(),
            num_cells: mesh.cells.len() - failed_cells,
            failed_cells,
        }
    }

    /// Analyze individual cell quality using correct vertex collection
    fn analyze_cell(&self, cell: &Cell, mesh: &Mesh<T>) -> Option<(T, T, T, T)> {
        // Use the correct unique_vertices method to eliminate duplicates
        let positions = cell.unique_vertices(mesh);
        
        if positions.len() < 3 {
            return None;
        }

        // Validate that the cell has the expected vertex count for its type
        if !cell.validate_vertex_count(mesh) {
            return None;
        }

        // Calculate metrics using element-type-specific implementations
        let aspect_ratio = self.calculate_aspect_ratio(cell, &positions);
        let skewness = self.calculate_skewness(cell, &positions);
        let orthogonality = self.calculate_orthogonality(cell, &positions);
        let volume = self.calculate_volume(cell, &positions);

        Some((aspect_ratio, skewness, orthogonality, volume))
    }

    /// Calculate aspect ratio using element-type-specific logic
    /// 
    /// Uses explicit element type dispatch instead of error-prone vertex count inference
    fn calculate_aspect_ratio(&self, cell: &Cell, vertices: &[&Point3<T>]) -> T {
        match cell.element_type {
            ElementType::Triangle => self.triangle_aspect_ratio(vertices),
            ElementType::Quadrilateral => self.quadrilateral_aspect_ratio(vertices),
            ElementType::Tetrahedron => self.tetrahedron_aspect_ratio(vertices),
            ElementType::Hexahedron => self.hexahedron_aspect_ratio(vertices),
            ElementType::Pyramid => self.pyramid_aspect_ratio(vertices),
            ElementType::Pentahedron => self.pentahedron_aspect_ratio(vertices),
            _ => T::one(), // Default for unsupported types
        }
    }

    /// Calculate skewness using element-type-specific logic with robust numerics
    fn calculate_skewness(&self, cell: &Cell, vertices: &[&Point3<T>]) -> T {
        match cell.element_type {
            ElementType::Triangle => self.triangle_skewness(vertices),
            ElementType::Quadrilateral => self.quadrilateral_skewness(vertices),
            ElementType::Tetrahedron => self.tetrahedron_skewness(vertices),
            ElementType::Hexahedron => self.hexahedron_skewness(vertices),
            ElementType::Pyramid => self.pyramid_skewness(vertices),
            ElementType::Pentahedron => self.pentahedron_skewness(vertices),
            _ => T::one(), // Worst quality for unsupported types
        }
    }

    /// Calculate orthogonality using element-type-specific logic
    fn calculate_orthogonality(&self, cell: &Cell, vertices: &[&Point3<T>]) -> T {
        match cell.element_type {
            ElementType::Triangle => T::one(), // Always orthogonal for 2D
            ElementType::Quadrilateral => self.quadrilateral_orthogonality(vertices),
            ElementType::Tetrahedron => self.tetrahedron_orthogonality(vertices),
            ElementType::Hexahedron => self.hexahedron_orthogonality(vertices),
            ElementType::Pyramid => self.pyramid_orthogonality(vertices),
            ElementType::Pentahedron => self.pentahedron_orthogonality(vertices),
            _ => T::zero(), // No orthogonality for unsupported types
        }
    }

    /// Calculate volume using element-type-specific formulas
    fn calculate_volume(&self, cell: &Cell, vertices: &[&Point3<T>]) -> T {
        match cell.element_type {
            ElementType::Triangle => self.triangle_area(vertices),
            ElementType::Quadrilateral => self.quadrilateral_area(vertices),
            ElementType::Tetrahedron => self.tetrahedron_volume(vertices),
            ElementType::Hexahedron => self.hexahedron_volume(vertices),
            ElementType::Pyramid => self.pyramid_volume(vertices),
            ElementType::Pentahedron => self.pentahedron_volume(vertices),
            _ => T::zero(),
        }
    }

    // Element-specific aspect ratio implementations

    fn triangle_aspect_ratio(&self, vertices: &[&Point3<T>]) -> T {
        if vertices.len() != 3 { return T::one(); }
        
        let edges = [
            (vertices[1] - vertices[0]).norm(),
            (vertices[2] - vertices[1]).norm(),
            (vertices[0] - vertices[2]).norm(),
        ];
        
        let min_edge = edges.iter().cloned().fold(T::infinity(), |a, b| Float::min(a, b));
        let max_edge = edges.iter().cloned().fold(T::neg_infinity(), |a, b| Float::max(a, b));
        
        if max_edge > T::zero() { min_edge / max_edge } else { T::one() }
    }

    fn quadrilateral_aspect_ratio(&self, vertices: &[&Point3<T>]) -> T {
        if vertices.len() != 4 { return T::one(); }
        
        let edges = [
            (vertices[1] - vertices[0]).norm(),
            (vertices[2] - vertices[1]).norm(),
            (vertices[3] - vertices[2]).norm(),
            (vertices[0] - vertices[3]).norm(),
        ];
        
        let min_edge = edges.iter().cloned().fold(T::infinity(), |a, b| Float::min(a, b));
        let max_edge = edges.iter().cloned().fold(T::neg_infinity(), |a, b| Float::max(a, b));
        
        if max_edge > T::zero() { min_edge / max_edge } else { T::one() }
    }

    fn tetrahedron_aspect_ratio(&self, vertices: &[&Point3<T>]) -> T {
        if vertices.len() != 4 { return T::one(); }
        
        let edges = [
            (vertices[1] - vertices[0]).norm(),
            (vertices[2] - vertices[0]).norm(),
            (vertices[3] - vertices[0]).norm(),
            (vertices[2] - vertices[1]).norm(),
            (vertices[3] - vertices[1]).norm(),
            (vertices[3] - vertices[2]).norm(),
        ];
        
        let min_edge = edges.iter().cloned().fold(T::infinity(), |a, b| Float::min(a, b));
        let max_edge = edges.iter().cloned().fold(T::neg_infinity(), |a, b| Float::max(a, b));
        
        if max_edge > T::zero() { min_edge / max_edge } else { T::one() }
    }

    fn hexahedron_aspect_ratio(&self, vertices: &[&Point3<T>]) -> T {
        if vertices.len() != 8 { return T::one(); }
        
        // Calculate edge lengths for a hexahedron (assuming standard ordering)
        let mut edges = Vec::with_capacity(12);
        
        // Bottom face edges (0-1-2-3)
        edges.extend_from_slice(&[
            (vertices[1] - vertices[0]).norm(),
            (vertices[2] - vertices[1]).norm(),
            (vertices[3] - vertices[2]).norm(),
            (vertices[0] - vertices[3]).norm(),
        ]);
        
        // Top face edges (4-5-6-7)
        edges.extend_from_slice(&[
            (vertices[5] - vertices[4]).norm(),
            (vertices[6] - vertices[5]).norm(),
            (vertices[7] - vertices[6]).norm(),
            (vertices[4] - vertices[7]).norm(),
        ]);
        
        // Vertical edges
        edges.extend_from_slice(&[
            (vertices[4] - vertices[0]).norm(),
            (vertices[5] - vertices[1]).norm(),
            (vertices[6] - vertices[2]).norm(),
            (vertices[7] - vertices[3]).norm(),
        ]);
        
        let min_edge = edges.iter().cloned().fold(T::infinity(), |a, b| Float::min(a, b));
        let max_edge = edges.iter().cloned().fold(T::neg_infinity(), |a, b| Float::max(a, b));
        
        if max_edge > T::zero() { min_edge / max_edge } else { T::one() }
    }

    fn pyramid_aspect_ratio(&self, vertices: &[&Point3<T>]) -> T {
        if vertices.len() != 5 { return T::one(); }
        
        // Base edges (assuming vertices 0-3 form the base)
        let mut edges = vec![
            (vertices[1] - vertices[0]).norm(),
            (vertices[2] - vertices[1]).norm(),
            (vertices[3] - vertices[2]).norm(),
            (vertices[0] - vertices[3]).norm(),
        ];
        
        // Apex edges (vertex 4 to base vertices)
        edges.extend_from_slice(&[
            (vertices[4] - vertices[0]).norm(),
            (vertices[4] - vertices[1]).norm(),
            (vertices[4] - vertices[2]).norm(),
            (vertices[4] - vertices[3]).norm(),
        ]);
        
        let min_edge = edges.iter().cloned().fold(T::infinity(), |a, b| Float::min(a, b));
        let max_edge = edges.iter().cloned().fold(T::neg_infinity(), |a, b| Float::max(a, b));
        
        if max_edge > T::zero() { min_edge / max_edge } else { T::one() }
    }

    fn pentahedron_aspect_ratio(&self, vertices: &[&Point3<T>]) -> T {
        if vertices.len() != 6 { return T::one(); }
        
        // Triangular prism edges
        let mut edges = Vec::with_capacity(9);
        
        // Bottom triangle edges (0-1-2)
        edges.extend_from_slice(&[
            (vertices[1] - vertices[0]).norm(),
            (vertices[2] - vertices[1]).norm(),
            (vertices[0] - vertices[2]).norm(),
        ]);
        
        // Top triangle edges (3-4-5)
        edges.extend_from_slice(&[
            (vertices[4] - vertices[3]).norm(),
            (vertices[5] - vertices[4]).norm(),
            (vertices[3] - vertices[5]).norm(),
        ]);
        
        // Vertical edges
        edges.extend_from_slice(&[
            (vertices[3] - vertices[0]).norm(),
            (vertices[4] - vertices[1]).norm(),
            (vertices[5] - vertices[2]).norm(),
        ]);
        
        let min_edge = edges.iter().cloned().fold(T::infinity(), |a, b| Float::min(a, b));
        let max_edge = edges.iter().cloned().fold(T::neg_infinity(), |a, b| Float::max(a, b));
        
        if max_edge > T::zero() { min_edge / max_edge } else { T::one() }
    }

    // Element-specific skewness implementations with robust numerics

    fn triangle_skewness(&self, vertices: &[&Point3<T>]) -> T {
        if vertices.len() != 3 { return T::one(); }
        
        // Calculate angles and check deviation from 60 degrees
        let edges = [
            vertices[1] - vertices[0],
            vertices[2] - vertices[1],
            vertices[0] - vertices[2],
        ];
        
        let mut max_deviation = T::zero();
        let ideal_angle = T::from(60.0_f64.to_radians()).expect("FIXME: Add proper error handling");
        
        for i in 0..3 {
            let e1 = &edges[i];
            let e2 = &edges[(i + 1) % 3];
            
            let cos_angle = -e1.dot(e2) / (e1.norm() * e2.norm());
            let cos_clamped = Float::min(Float::max(cos_angle, -T::one()), T::one());
            let angle = Float::acos(cos_clamped);
            
            let deviation = num_traits::Float::abs(angle - ideal_angle) / ideal_angle;
            max_deviation = Float::max(max_deviation, deviation);
        }
        
        Float::min(max_deviation, T::one())
    }

    fn quadrilateral_skewness(&self, vertices: &[&Point3<T>]) -> T {
        if vertices.len() != 4 { return T::one(); }
        
        // Check deviation of angles from 90 degrees
        let mut max_deviation = T::zero();
        let ideal_angle = T::from(90.0_f64.to_radians()).expect("FIXME: Add proper error handling");
        
        for i in 0..4 {
            let v1 = vertices[(i + 1) % 4] - vertices[i];
            let v2 = vertices[(i + 3) % 4] - vertices[i];
            
            let cos_angle = v1.dot(&v2) / (v1.norm() * v2.norm());
            let cos_clamped = Float::min(Float::max(cos_angle, -T::one()), T::one());
            let angle = Float::acos(cos_clamped);
            
            let deviation = num_traits::Float::abs(angle - ideal_angle) / ideal_angle;
            max_deviation = Float::max(max_deviation, deviation);
        }
        
        Float::min(max_deviation, T::one())
    }

    fn tetrahedron_skewness(&self, vertices: &[&Point3<T>]) -> T {
        if vertices.len() != 4 { return T::one(); }

        // Face-based skewness: max skewness across the four triangular faces
        let faces = [
            [vertices[0], vertices[1], vertices[2]],
            [vertices[0], vertices[2], vertices[3]],
            [vertices[0], vertices[3], vertices[1]],
            [vertices[1], vertices[3], vertices[2]],
        ];

        let mut max_face_skewness = T::zero();
        for face in &faces {
            let skew = self.triangle_skewness(face);
            max_face_skewness = Float::max(max_face_skewness, skew);
        }
        max_face_skewness
    }

    fn hexahedron_skewness(&self, vertices: &[&Point3<T>]) -> T {
        if vertices.len() != 8 { return T::one(); }
        
        // Check deviation of all corner angles from 90 degrees
        let mut max_deviation = T::zero();
        let ideal_angle = T::from(90.0_f64.to_radians()).expect("FIXME: Add proper error handling");
        
        // Check each corner (assuming standard hexahedron ordering)
        let corner_edges = [
            [(1, 0), (3, 0), (4, 0)], // Corner 0
            [(0, 1), (2, 1), (5, 1)], // Corner 1
            [(1, 2), (3, 2), (6, 2)], // Corner 2
            [(2, 3), (0, 3), (7, 3)], // Corner 3
            [(0, 4), (5, 4), (7, 4)], // Corner 4
            [(1, 5), (4, 5), (6, 5)], // Corner 5
            [(2, 6), (5, 6), (7, 6)], // Corner 6
            [(3, 7), (4, 7), (6, 7)], // Corner 7
        ];
        
        for corner in &corner_edges {
            for i in 0..3 {
                for j in (i+1)..3 {
                    let (a1, b1) = corner[i];
                    let (a2, b2) = corner[j];
                    
                    if a1 < vertices.len() && b1 < vertices.len() && 
                       a2 < vertices.len() && b2 < vertices.len() {
                        let e1 = vertices[a1] - vertices[b1];
                        let e2 = vertices[a2] - vertices[b2];
                        
                        let cos_angle = e1.dot(&e2) / (e1.norm() * e2.norm());
                        let cos_clamped = Float::min(Float::max(cos_angle, -T::one()), T::one());
                        let angle = Float::acos(cos_clamped);
                        
                        let deviation = num_traits::Float::abs(angle - ideal_angle) / ideal_angle;
                        max_deviation = Float::max(max_deviation, deviation);
                    }
                }
            }
        }
        
        Float::min(max_deviation, T::one())
    }

    fn pyramid_skewness(&self, vertices: &[&Point3<T>]) -> T {
        if vertices.len() != 5 { return T::one(); }
        
        // Check base quadrilateral skewness and apex angles
        let base_skewness = self.quadrilateral_skewness(&vertices[0..4]);
        
        // Check apex angles (should be close to some ideal)
        let mut apex_deviation = T::zero();
        for i in 0..4 {
            let e1 = vertices[i] - vertices[4];
            let e2 = vertices[(i + 1) % 4] - vertices[4];
            
            let cos_angle = e1.dot(&e2) / (e1.norm() * e2.norm());
            let cos_clamped = Float::min(Float::max(cos_angle, -T::one()), T::one());
            let angle = Float::acos(cos_clamped);
            
            // For a regular pyramid, apex angles should be equal
            let expected_angle = T::from(90.0_f64.to_radians()).expect("FIXME: Add proper error handling");
            let deviation = num_traits::Float::abs(angle - expected_angle) / expected_angle;
            apex_deviation = Float::max(apex_deviation, deviation);
        }
        
        Float::min(Float::max(base_skewness, apex_deviation), T::one())
    }

    fn pentahedron_skewness(&self, vertices: &[&Point3<T>]) -> T {
        if vertices.len() != 6 { return T::one(); }
        
        // Check triangular face skewness
        let bottom_skewness = self.triangle_skewness(&vertices[0..3]);
        let top_skewness = self.triangle_skewness(&vertices[3..6]);
        
        // Check rectangular face angles
        let mut rect_skewness = T::zero();
        let ideal_angle = T::from(90.0_f64.to_radians()).expect("FIXME: Add proper error handling");
        
        // Check the three rectangular faces
        for i in 0..3 {
            let rect_vertices = [
                vertices[i],
                vertices[(i + 1) % 3],
                vertices[((i + 1) % 3) + 3],
                vertices[i + 3],
            ];
            
            let face_skewness = self.quadrilateral_skewness(&rect_vertices);
            rect_skewness = Float::max(rect_skewness, face_skewness);
        }
        
        let a = Float::max(bottom_skewness, top_skewness);
        let b = Float::max(a, rect_skewness);
        Float::min(b, T::one())
    }

    // Element-specific orthogonality implementations

    fn quadrilateral_orthogonality(&self, vertices: &[&Point3<T>]) -> T {
        if vertices.len() != 4 { return T::zero(); }
        
        // Check how close adjacent edges are to perpendicular
        let mut min_orthogonality = T::one();
        
        for i in 0..4 {
            let e1 = vertices[(i + 1) % 4] - vertices[i];
            let e2 = vertices[(i + 3) % 4] - vertices[i];
            
            let cos_angle = e1.dot(&e2) / (e1.norm() * e2.norm());
            let orthogonality = T::one() - num_traits::Float::abs(cos_angle);
            min_orthogonality = Float::min(min_orthogonality, orthogonality);
        }
        
        Float::max(min_orthogonality, T::zero())
    }

    fn tetrahedron_orthogonality(&self, vertices: &[&Point3<T>]) -> T {
        if vertices.len() != 4 { return T::zero(); }
        
        // Calculate face normals and check their orthogonality
        let face_normals = [
            (vertices[1] - vertices[0]).cross(&(vertices[2] - vertices[0])),
            (vertices[1] - vertices[0]).cross(&(vertices[3] - vertices[0])),
            (vertices[2] - vertices[0]).cross(&(vertices[3] - vertices[0])),
            (vertices[2] - vertices[1]).cross(&(vertices[3] - vertices[1])),
        ];
        
        let mut min_orthogonality = T::one();
        
        for i in 0..face_normals.len() {
            for j in (i+1)..face_normals.len() {
                let norm_i = face_normals[i].norm();
                let norm_j = face_normals[j].norm();
                
                if norm_i > T::zero() && norm_j > T::zero() {
                    let cos_angle = face_normals[i].dot(&face_normals[j]) / (norm_i * norm_j);
                    let orthogonality = T::one() - num_traits::Float::abs(cos_angle);
                    min_orthogonality = Float::min(min_orthogonality, orthogonality);
                }
            }
        }
        
        Float::max(min_orthogonality, T::zero())
    }

    fn hexahedron_orthogonality(&self, vertices: &[&Point3<T>]) -> T {
        if vertices.len() != 8 { return T::zero(); }
        
        // Check orthogonality of face normals
        let face_normals = [
            // Bottom face (0-1-2-3)
            (vertices[1] - vertices[0]).cross(&(vertices[3] - vertices[0])),
            // Top face (4-5-6-7)
            (vertices[5] - vertices[4]).cross(&(vertices[7] - vertices[4])),
            // Front face (0-1-5-4)
            (vertices[1] - vertices[0]).cross(&(vertices[4] - vertices[0])),
            // Back face (2-3-7-6)
            (vertices[3] - vertices[2]).cross(&(vertices[6] - vertices[2])),
            // Left face (0-3-7-4)
            (vertices[3] - vertices[0]).cross(&(vertices[4] - vertices[0])),
            // Right face (1-2-6-5)
            (vertices[2] - vertices[1]).cross(&(vertices[5] - vertices[1])),
        ];
        
        let mut min_orthogonality = T::one();
        
        // Check orthogonality between adjacent faces
        let adjacent_pairs = [
            (0, 2), (0, 3), (0, 4), (0, 5), // Bottom with others
            (1, 2), (1, 3), (1, 4), (1, 5), // Top with others
            (2, 4), (2, 5), // Front with sides
            (3, 4), (3, 5), // Back with sides
        ];
        
        for &(i, j) in &adjacent_pairs {
            let norm_i = face_normals[i].norm();
            let norm_j = face_normals[j].norm();
            
            if norm_i > T::zero() && norm_j > T::zero() {
                let cos_angle = face_normals[i].dot(&face_normals[j]) / (norm_i * norm_j);
                let orthogonality = T::one() - num_traits::Float::abs(cos_angle);
                min_orthogonality = Float::min(min_orthogonality, orthogonality);
            }
        }
        
        Float::max(min_orthogonality, T::zero())
    }

    fn pyramid_orthogonality(&self, vertices: &[&Point3<T>]) -> T {
        if vertices.len() != 5 { return T::zero(); }
        
        // Check orthogonality between base and triangular faces
        let base_normal = (vertices[1] - vertices[0]).cross(&(vertices[3] - vertices[0]));
        
        let mut min_orthogonality = T::one();
        
        // Check each triangular face against the base
        for i in 0..4 {
            let tri_normal = (vertices[(i + 1) % 4] - vertices[i])
                .cross(&(vertices[4] - vertices[i]));
            
            let base_norm = base_normal.norm();
            let tri_norm = tri_normal.norm();
            
            if base_norm > T::zero() && tri_norm > T::zero() {
                let cos_angle = base_normal.dot(&tri_normal) / (base_norm * tri_norm);
                let orthogonality = T::one() - num_traits::Float::abs(cos_angle);
                min_orthogonality = Float::min(min_orthogonality, orthogonality);
            }
        }
        
        Float::max(min_orthogonality, T::zero())
    }

    fn pentahedron_orthogonality(&self, vertices: &[&Point3<T>]) -> T {
        if vertices.len() != 6 { return T::zero(); }
        
        // Check orthogonality between triangular faces and rectangular faces
        let bottom_normal = (vertices[1] - vertices[0]).cross(&(vertices[2] - vertices[0]));
        let top_normal = (vertices[4] - vertices[3]).cross(&(vertices[5] - vertices[3]));
        
        let mut min_orthogonality = T::one();
        
        // Check rectangular faces against triangular faces
        for i in 0..3 {
            let rect_normal = (vertices[(i + 1) % 3] - vertices[i])
                .cross(&(vertices[i + 3] - vertices[i]));
            
            // Check against bottom face
            let bottom_norm = bottom_normal.norm();
            let rect_norm = rect_normal.norm();
            
            if bottom_norm > T::zero() && rect_norm > T::zero() {
                let cos_angle = bottom_normal.dot(&rect_normal) / (bottom_norm * rect_norm);
                let orthogonality = T::one() - num_traits::Float::abs(cos_angle);
                min_orthogonality = Float::min(min_orthogonality, orthogonality);
            }
            
            // Check against top face
            let top_norm = top_normal.norm();
            if top_norm > T::zero() && rect_norm > T::zero() {
                let cos_angle = top_normal.dot(&rect_normal) / (top_norm * rect_norm);
                let orthogonality = T::one() - num_traits::Float::abs(cos_angle);
                min_orthogonality = Float::min(min_orthogonality, orthogonality);
            }
        }
        
        Float::max(min_orthogonality, T::zero())
    }

    // Element-specific volume/area calculations

    fn triangle_area(&self, vertices: &[&Point3<T>]) -> T {
        if vertices.len() != 3 { return T::zero(); }
        
        let v1 = vertices[1] - vertices[0];
        let v2 = vertices[2] - vertices[0];
        v1.cross(&v2).norm() / T::from(2.0).unwrap_or_else(|_| T::zero())
    }

    fn quadrilateral_area(&self, vertices: &[&Point3<T>]) -> T {
        if vertices.len() != 4 { return T::zero(); }
        
        // Split into two triangles and sum their areas
        let tri1_area = {
            let v1 = vertices[1] - vertices[0];
            let v2 = vertices[2] - vertices[0];
            v1.cross(&v2).norm() / T::from(2.0).unwrap_or_else(|_| T::zero())
        };
        
        let tri2_area = {
            let v1 = vertices[2] - vertices[0];
            let v2 = vertices[3] - vertices[0];
            v1.cross(&v2).norm() / T::from(2.0).unwrap_or_else(|_| T::zero())
        };
        
        tri1_area + tri2_area
    }

    fn tetrahedron_volume(&self, vertices: &[&Point3<T>]) -> T {
        if vertices.len() != 4 { return T::zero(); }
        
        let v1 = vertices[1] - vertices[0];
        let v2 = vertices[2] - vertices[0];
        let v3 = vertices[3] - vertices[0];
        
        num_traits::Float::abs(v1.cross(&v2).dot(&v3)) / T::from(6.0).unwrap_or_else(|_| T::zero())
    }

    fn hexahedron_volume(&self, vertices: &[&Point3<T>]) -> T {
        if vertices.len() != 8 { return T::zero(); }

        // Efficient 5-tetrahedra decomposition (assumes standard hexa ordering)
        let v = |i: usize| vertices[i];
        self.tetrahedron_volume(&[v(0), v(1), v(3), v(4)]) +
        self.tetrahedron_volume(&[v(1), v(2), v(3), v(6)]) +
        self.tetrahedron_volume(&[v(1), v(4), v(5), v(6)]) +
        self.tetrahedron_volume(&[v(3), v(4), v(6), v(7)]) +
        self.tetrahedron_volume(&[v(1), v(3), v(4), v(6)])
    }

    fn pyramid_volume(&self, vertices: &[&Point3<T>]) -> T {
        if vertices.len() != 5 { return T::zero(); }
        
        // Volume = (1/3) * base_area * height
        let base_area = self.quadrilateral_area(&vertices[0..4]);
        
        // Calculate height as distance from apex to base plane
        let base_normal = (vertices[1] - vertices[0]).cross(&(vertices[3] - vertices[0]));
        let base_normal_unit = base_normal / base_normal.norm();
        let height = num_traits::Float::abs((vertices[4] - vertices[0]).dot(&base_normal_unit));
        
        base_area * height / T::from(3.0).unwrap_or_else(|_| T::zero())
    }

    fn pentahedron_volume(&self, vertices: &[&Point3<T>]) -> T {
        if vertices.len() != 6 { return T::zero(); }
        let v = |i: usize| vertices[i];
        // Standard triangular prism decomposition: 3 tets
        let vol1 = self.tetrahedron_volume(&[v(0), v(1), v(2), v(5)]);
        let vol2 = self.tetrahedron_volume(&[v(0), v(1), v(3), v(5)]);
        let vol3 = self.tetrahedron_volume(&[v(1), v(3), v(4), v(5)]);
        vol1 + vol2 + vol3
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh::{Vertex, Face, Mesh, ElementType};
    use nalgebra::Point3;
    use approx::assert_relative_eq;

    /// Create a simple test mesh with a single tetrahedron
    fn create_tetrahedron_mesh() -> Mesh<f64> {
        let vertices = vec![
            Vertex { id: 0, position: Point3::new(0.0, 0.0, 0.0) },
            Vertex { id: 1, position: Point3::new(1.0, 0.0, 0.0) },
            Vertex { id: 2, position: Point3::new(0.0, 1.0, 0.0) },
            Vertex { id: 3, position: Point3::new(0.0, 0.0, 1.0) },
        ];

        let faces = vec![
            Face { id: 0, vertices: vec![0, 1, 2] }, // Bottom face
            Face { id: 1, vertices: vec![0, 1, 3] }, // Front face
            Face { id: 2, vertices: vec![1, 2, 3] }, // Right face
            Face { id: 3, vertices: vec![0, 2, 3] }, // Left face
        ];

        let cells = vec![
            Cell::new(0, vec![0, 1, 2, 3], ElementType::Tetrahedron),
        ];

        let mut mesh = Mesh {
            vertices,
            edges: vec![],
            faces,
            cells,
            topology: Default::default(),
        };
        mesh.update_topology();
        mesh
    }

    /// Create a test mesh with a hexahedron (cube)
    fn create_hexahedron_mesh() -> Mesh<f64> {
        let vertices = vec![
            // Bottom face
            Vertex { id: 0, position: Point3::new(0.0, 0.0, 0.0) },
            Vertex { id: 1, position: Point3::new(1.0, 0.0, 0.0) },
            Vertex { id: 2, position: Point3::new(1.0, 1.0, 0.0) },
            Vertex { id: 3, position: Point3::new(0.0, 1.0, 0.0) },
            // Top face
            Vertex { id: 4, position: Point3::new(0.0, 0.0, 1.0) },
            Vertex { id: 5, position: Point3::new(1.0, 0.0, 1.0) },
            Vertex { id: 6, position: Point3::new(1.0, 1.0, 1.0) },
            Vertex { id: 7, position: Point3::new(0.0, 1.0, 1.0) },
        ];

        let faces = vec![
            Face { id: 0, vertices: vec![0, 1, 2, 3] }, // Bottom
            Face { id: 1, vertices: vec![4, 5, 6, 7] }, // Top
            Face { id: 2, vertices: vec![0, 1, 5, 4] }, // Front
            Face { id: 3, vertices: vec![2, 3, 7, 6] }, // Back
            Face { id: 4, vertices: vec![0, 3, 7, 4] }, // Left
            Face { id: 5, vertices: vec![1, 2, 6, 5] }, // Right
        ];

        let cells = vec![
            Cell::new(0, vec![0, 1, 2, 3, 4, 5], ElementType::Hexahedron),
        ];

        let mut mesh = Mesh {
            vertices,
            edges: vec![],
            faces,
            cells,
            topology: Default::default(),
        };
        mesh.update_topology();
        mesh
    }

    #[test]
    fn test_unique_vertices_collection() {
        let mesh = create_tetrahedron_mesh();
        let cell = &mesh.cells[0];
        
        // Verify that unique_vertices returns exactly 4 vertices (no duplicates)
        let unique_vertices = cell.unique_vertices(&mesh);
        assert_eq!(unique_vertices.len(), 4, "Tetrahedron should have exactly 4 unique vertices");
        
        // Verify all expected positions are present
        let positions: Vec<_> = unique_vertices.iter().map(|v| *v).collect();
        assert!(positions.contains(&&Point3::new(0.0, 0.0, 0.0)));
        assert!(positions.contains(&&Point3::new(1.0, 0.0, 0.0)));
        assert!(positions.contains(&&Point3::new(0.0, 1.0, 0.0)));
        assert!(positions.contains(&&Point3::new(0.0, 0.0, 1.0)));
    }

    #[test]
    fn test_single_pass_analysis() {
        let mesh = create_tetrahedron_mesh();
        let analyzer = QualityAnalyzer::default();
        
        // Test that analysis completes successfully with single pass
        let metrics = analyzer.analyze(&mesh);
        
        assert_eq!(metrics.num_cells, 1);
        assert_eq!(metrics.failed_cells, 0);
        
        // Verify statistics are properly computed
        assert_eq!(metrics.aspect_ratio.count, 1);
        assert_eq!(metrics.skewness.count, 1);
        assert_eq!(metrics.orthogonality.count, 1);
        assert_eq!(metrics.volume.count, 1);
        
        // Verify volume is computed correctly for unit tetrahedron
        let expected_volume = 1.0 / 6.0; // Volume of tetrahedron with edge vectors (1,0,0), (0,1,0), (0,0,1)
        assert_relative_eq!(metrics.volume.mean, expected_volume, epsilon = 1e-10);
    }

    #[test]
    fn test_element_type_specific_calculations() {
        let tet_mesh = create_tetrahedron_mesh();
        let hex_mesh = create_hexahedron_mesh();
        let analyzer = QualityAnalyzer::default();
        
        // Test tetrahedron
        let tet_metrics = analyzer.analyze(&tet_mesh);
        assert_eq!(tet_metrics.num_cells, 1);
        
        // Test hexahedron
        let hex_metrics = analyzer.analyze(&hex_mesh);
        assert_eq!(hex_metrics.num_cells, 1);
        
        // Volume should be 1.0 for unit cube
        assert_relative_eq!(hex_metrics.volume.mean, 1.0, epsilon = 1e-6);
        
        // Aspect ratio should be 1.0 for perfect cube (all edges equal)
        assert_relative_eq!(hex_metrics.aspect_ratio.mean, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_robust_floating_point_handling() {
        let mut mesh = create_tetrahedron_mesh();
        
        // Create a nearly degenerate tetrahedron to test numerical robustness
        mesh.vertices[3].position = Point3::new(1e-15, 1e-15, 1e-15);
        
        let analyzer = QualityAnalyzer::default();
        let metrics = analyzer.analyze(&mesh);
        
        // Should not panic or produce NaN/infinite values
        assert!(metrics.aspect_ratio.mean.is_finite());
        assert!(metrics.skewness.mean.is_finite());
        assert!(metrics.orthogonality.mean.is_finite());
        assert!(metrics.volume.mean.is_finite());
    }

    #[test]
    fn test_explicit_element_types() {
        let tet_mesh = create_tetrahedron_mesh();
        let hex_mesh = create_hexahedron_mesh();
        
        // Verify element types are correctly set
        assert_eq!(tet_mesh.cells[0].element_type, ElementType::Tetrahedron);
        assert_eq!(hex_mesh.cells[0].element_type, ElementType::Hexahedron);
        
        // Verify validation works
        assert!(tet_mesh.cells[0].validate_vertex_count(&tet_mesh));
        assert!(hex_mesh.cells[0].validate_vertex_count(&hex_mesh));
        
        // Test expected vertex counts
        assert_eq!(ElementType::Tetrahedron.expected_vertex_count(), 4);
        assert_eq!(ElementType::Hexahedron.expected_vertex_count(), 8);
        assert_eq!(ElementType::Triangle.expected_vertex_count(), 3);
        assert_eq!(ElementType::Quadrilateral.expected_vertex_count(), 4);
    }

    #[test]
    fn test_running_statistics_accuracy() {
        // Test that running statistics (Welford's algorithm) produces correct results
        let mut stats = RunningStats::<f64>::new();
        
        let values = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        for value in &values {
            stats.update(*value);
        }
        
        let result = stats.finalize();
        
        // Mean should be 3.0
        assert_relative_eq!(result.mean, 3.0, epsilon = 1e-10);
        
        // Min should be 1.0, max should be 5.0
        assert_relative_eq!(result.min, 1.0, epsilon = 1e-10);
        assert_relative_eq!(result.max, 5.0, epsilon = 1e-10);
        
        // Standard deviation should be sqrt(2.5)
        let expected_std = (2.5_f64).sqrt();
        assert_relative_eq!(result.std_dev, expected_std, epsilon = 1e-10);
        
        assert_eq!(result.count, 5);
    }

    #[test]
    fn test_element_type_inference() {
        assert_eq!(ElementType::infer_from_vertex_count(3), ElementType::Triangle);
        assert_eq!(ElementType::infer_from_vertex_count(4), ElementType::Tetrahedron);
        assert_eq!(ElementType::infer_from_vertex_count(8), ElementType::Hexahedron);
        assert_eq!(ElementType::infer_from_vertex_count(100), ElementType::Polyhedron);
    }

    #[test]
    fn test_triangle_calculations() {
        // Create a simple triangle mesh
        let vertices = vec![
            Vertex { id: 0, position: Point3::new(0.0, 0.0, 0.0) },
            Vertex { id: 1, position: Point3::new(1.0, 0.0, 0.0) },
            Vertex { id: 2, position: Point3::new(0.0, 1.0, 0.0) },
        ];

        let faces = vec![
            Face { id: 0, vertices: vec![0, 1, 2] },
        ];

        let cells = vec![
            Cell::new(0, vec![0], ElementType::Triangle),
        ];

        let mesh = Mesh {
            vertices,
            edges: vec![],
            faces,
            cells,
            topology: Default::default(),
        };

        let analyzer = QualityAnalyzer::default();
        let metrics = analyzer.analyze(&mesh);
        
        // Area should be 0.5 for right triangle with legs of length 1
        assert_relative_eq!(metrics.volume.mean, 0.5, epsilon = 1e-10);
        
        // Perfect right triangle should have good aspect ratio
        assert!(metrics.aspect_ratio.mean > 0.5);
    }

    #[test]
    fn test_performance_with_multiple_cells() {
        // Create a mesh with multiple identical tetrahedra to test single-pass efficiency
        let base_mesh = create_tetrahedron_mesh();
        let mut vertices = Vec::new();
        let mut faces = Vec::new();
        let mut cells = Vec::new();
        
        // Create 100 copies of the tetrahedron
        for i in 0..100 {
            let offset = i * 4;
            
            // Add vertices
            for v in &base_mesh.vertices {
                vertices.push(Vertex {
                    id: v.id + offset,
                    position: v.position + nalgebra::Vector3::new(i as f64 * 2.0, 0.0, 0.0),
                });
            }
            
            // Add faces
            for f in &base_mesh.faces {
                faces.push(Face {
                    id: f.id + i * 4,
                    vertices: f.vertices.iter().map(|v| v + offset).collect(),
                });
            }
            
            // Add cell
            cells.push(Cell::new(
                i,
                (i * 4..(i + 1) * 4).collect(),
                ElementType::Tetrahedron,
            ));
        }
        
        let mesh = Mesh {
            vertices,
            edges: vec![],
            faces,
            cells,
            topology: Default::default(),
        };
        
        let analyzer = QualityAnalyzer::default();
        let metrics = analyzer.analyze(&mesh);
        
        // Should successfully analyze all 100 cells
        assert_eq!(metrics.num_cells, 100);
        assert_eq!(metrics.failed_cells, 0);
        
        // All cells are identical, so variance should be zero
        assert_relative_eq!(metrics.volume.std_dev, 0.0, epsilon = 1e-10);
        assert_relative_eq!(metrics.aspect_ratio.std_dev, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_mesh_validation() {
        let mut mesh = create_tetrahedron_mesh();
        
        // Should pass validation with correct element type
        assert!(mesh.validate().is_ok());
        
        // Break the mesh by changing element type to something invalid
        mesh.cells[0].element_type = ElementType::Hexahedron; // Wrong type for 4 vertices
        
        // Should fail validation
        assert!(mesh.validate().is_err());
    }
}
