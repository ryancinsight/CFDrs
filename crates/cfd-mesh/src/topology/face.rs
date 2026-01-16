//! Face representation for 2D elements in 3D meshes
//!
//! A face is a 2-simplex embedded in ℝ³, defined by an ordered list of vertex indices.
//! Faces form the boundaries of 3D cells and are fundamental to mesh topology.
//!
//! # Mathematical Background
//!
//! A face requires at least 3 non-collinear vertices to define a valid 2-simplex.
//! The vertex ordering follows counter-clockwise (CCW) winding convention when
//! viewed from the outside of the cell, defining the outward normal direction
//! via the right-hand rule: **n** = (**v₁** - **v₀**) × (**v₂** - **v₀**).
//!
//! # Topological vs Geometric Invariants
//!
//! This module validates **topological invariants** (index-based):
//! - Minimum 3 vertices (simplex condition)
//! - All vertex indices must be distinct (non-degenerate)
//!
//! **Geometric invariants** require vertex coordinates and are validated at `Mesh` level:
//! - Non-collinearity (area > 0)
//! - Planarity (for n > 3 vertices)
//! - Convexity (optional, application-dependent)
//! - Correct winding direction
//!
//! # Complexity Guarantees
//!
//! Let n = number of vertices in face:
//! - Construction: O(1) unchecked, O(n) validated
//! - `contains()`: O(n) linear search
//! - `edges()`: O(n) allocation + iteration
//! - `canonical_edges()`: O(n) iteration, zero allocation
//! - `validate()`: O(n) time, O(n) space for duplicate detection
//! - `shared_vertices()`: O(n + m) with HashSet optimization
//! - `shares_edge_with()`: O(n + m) with HashSet lookup
//!
//! # Preconditions
//!
//! Methods assume a topologically valid face (≥3 distinct vertices) unless
//! otherwise documented. Behavior is undefined for invalid faces.
//!
//! # References
//!
//! - Edelsbrunner, H. "Geometry and Topology for Mesh Generation" (2001), §2.1
//! - Shewchuk, J. "Triangle: Engineering a 2D Quality Mesh Generator" (1996)
//!
//! # Examples
//!
//! ```
//! use cfd_mesh::topology::Face;
//!
//! // Create a triangular face
//! let triangle = Face::triangle(0, 1, 2);
//! assert!(triangle.is_triangle());
//! assert_eq!(triangle.vertex_count(), 3);
//!
//! // Create a quadrilateral face
//! let quad = Face::quad(0, 1, 2, 3);
//! assert!(quad.is_quad());
//!
//! // Validate face topology
//! assert!(triangle.validate().is_ok());
//! ```

use crate::error::{FaceError, FaceResult};
use crate::topology::Edge;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashSet};

/// Minimum number of vertices required for a valid face (2-simplex)
pub const MIN_FACE_VERTICES: usize = 3;

/// Face defined by vertex indices
///
/// Represents a 2D polygonal face in a 3D mesh, storing vertex indices
/// in counter-clockwise order (when viewed from outside the cell).
///
/// # Fields
///
/// - `vertices`: Ordered vertex indices defining the face boundary
/// - `global_id`: Optional global identifier for distributed mesh partitioning
/// - `partition_id`: Optional MPI rank that owns this face
///
/// # Invariants
///
/// - `vertices.len() >= 3` (simplex condition)
/// - All elements in `vertices` are distinct (non-degenerate)
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Face {
    /// Ordered list of vertex indices (CCW winding convention)
    pub vertices: Vec<usize>,
    /// Global ID for distributed meshes (MPI partitioning)
    pub global_id: Option<usize>,
    /// Partition ID (MPI rank) that owns this face
    pub partition_id: Option<usize>,
}

impl Face {
    /// Create a triangular face from three vertex indices.
    ///
    /// # Arguments
    ///
    /// * `v0`, `v1`, `v2` - Vertex indices in CCW order
    ///
    /// # Examples
    ///
    /// ```
    /// use cfd_mesh::topology::Face;
    ///
    /// let face = Face::triangle(0, 1, 2);
    /// assert_eq!(face.vertex_count(), 3);
    /// assert!(face.is_triangle());
    /// ```
    #[must_use]
    pub fn triangle(v0: usize, v1: usize, v2: usize) -> Self {
        Self {
            vertices: vec![v0, v1, v2],
            global_id: None,
            partition_id: None,
        }
    }

    /// Create a quadrilateral face from four vertex indices.
    ///
    /// # Arguments
    ///
    /// * `v0`, `v1`, `v2`, `v3` - Vertex indices in CCW order
    ///
    /// # Examples
    ///
    /// ```
    /// use cfd_mesh::topology::Face;
    ///
    /// let face = Face::quad(0, 1, 2, 3);
    /// assert_eq!(face.vertex_count(), 4);
    /// assert!(face.is_quad());
    /// ```
    #[must_use]
    pub fn quad(v0: usize, v1: usize, v2: usize, v3: usize) -> Self {
        Self {
            vertices: vec![v0, v1, v2, v3],
            global_id: None,
            partition_id: None,
        }
    }

    /// Create a face from a vector of vertex indices.
    ///
    /// # Note
    ///
    /// This constructor does not validate the input. Use [`Face::try_new`]
    /// for validated construction that enforces invariants.
    ///
    /// # Arguments
    ///
    /// * `vertices` - Ordered vertex indices (CCW winding)
    ///
    /// # Examples
    ///
    /// ```
    /// use cfd_mesh::topology::Face;
    ///
    /// let pentagon = Face::from_vertices(vec![0, 1, 2, 3, 4]);
    /// assert_eq!(pentagon.vertex_count(), 5);
    /// ```
    #[must_use]
    pub fn from_vertices(vertices: Vec<usize>) -> Self {
        Self {
            vertices,
            global_id: None,
            partition_id: None,
        }
    }

    /// Create a validated face from vertex indices.
    ///
    /// # Errors
    ///
    /// Returns [`FaceError::InsufficientVertices`] if fewer than 3 vertices provided.
    /// Returns [`FaceError::DuplicateVertex`] if any vertex index appears more than once.
    ///
    /// # Examples
    ///
    /// ```
    /// use cfd_mesh::topology::Face;
    /// use cfd_mesh::error::FaceError;
    ///
    /// // Valid triangle
    /// let face = Face::try_new(vec![0, 1, 2]).unwrap();
    /// assert!(face.is_triangle());
    ///
    /// // Invalid: too few vertices
    /// let err = Face::try_new(vec![0, 1]).unwrap_err();
    /// assert!(matches!(err, FaceError::InsufficientVertices { .. }));
    ///
    /// // Invalid: duplicate vertex
    /// let err = Face::try_new(vec![0, 1, 1]).unwrap_err();
    /// assert!(matches!(err, FaceError::DuplicateVertex { .. }));
    /// ```
    pub fn try_new(vertices: Vec<usize>) -> FaceResult<Self> {
        let face = Self {
            vertices,
            global_id: None,
            partition_id: None,
        };
        face.validate()?;
        Ok(face)
    }

    /// Create a validated triangular face.
    ///
    /// # Errors
    ///
    /// Returns [`FaceError::DuplicateVertex`] if any vertices are equal.
    ///
    /// # Examples
    ///
    /// ```
    /// use cfd_mesh::topology::Face;
    /// use cfd_mesh::error::FaceError;
    ///
    /// let face = Face::try_triangle(0, 1, 2).unwrap();
    /// assert!(face.is_triangle());
    ///
    /// let err = Face::try_triangle(0, 0, 1).unwrap_err();
    /// assert!(matches!(err, FaceError::DuplicateVertex { .. }));
    /// ```
    pub fn try_triangle(v0: usize, v1: usize, v2: usize) -> FaceResult<Self> {
        Self::try_new(vec![v0, v1, v2])
    }

    /// Create a validated quadrilateral face.
    ///
    /// # Errors
    ///
    /// Returns [`FaceError::DuplicateVertex`] if any vertices are equal.
    ///
    /// # Examples
    ///
    /// ```
    /// use cfd_mesh::topology::Face;
    /// use cfd_mesh::error::FaceError;
    ///
    /// let face = Face::try_quad(0, 1, 2, 3).unwrap();
    /// assert!(face.is_quad());
    ///
    /// let err = Face::try_quad(0, 1, 2, 0).unwrap_err();
    /// assert!(matches!(err, FaceError::DuplicateVertex { .. }));
    /// ```
    pub fn try_quad(v0: usize, v1: usize, v2: usize, v3: usize) -> FaceResult<Self> {
        Self::try_new(vec![v0, v1, v2, v3])
    }

    /// Set distributed mesh properties using builder pattern.
    ///
    /// # Arguments
    ///
    /// * `global_id` - Global identifier across all MPI partitions
    /// * `partition_id` - MPI rank that owns this face
    ///
    /// # Examples
    ///
    /// ```
    /// use cfd_mesh::topology::Face;
    ///
    /// let face = Face::triangle(0, 1, 2)
    ///     .with_distributed_info(100, 2);
    /// assert_eq!(face.global_id, Some(100));
    /// assert_eq!(face.partition_id, Some(2));
    /// ```
    #[must_use]
    pub fn with_distributed_info(mut self, global_id: usize, partition_id: usize) -> Self {
        self.global_id = Some(global_id);
        self.partition_id = Some(partition_id);
        self
    }

    /// Validate face topology invariants.
    ///
    /// # Invariants Checked
    ///
    /// 1. Minimum vertex count (≥3 for 2-simplex)
    /// 2. No duplicate vertex indices
    ///
    /// # Errors
    ///
    /// Returns [`FaceError::InsufficientVertices`] if `vertices.len() < 3`.
    /// Returns [`FaceError::DuplicateVertex`] if any vertex appears multiple times.
    ///
    /// # Examples
    ///
    /// ```
    /// use cfd_mesh::topology::Face;
    ///
    /// let valid = Face::triangle(0, 1, 2);
    /// assert!(valid.validate().is_ok());
    ///
    /// let invalid = Face::from_vertices(vec![0, 1]);
    /// assert!(invalid.validate().is_err());
    /// ```
    pub fn validate(&self) -> FaceResult<()> {
        // Check minimum vertex count (simplex condition)
        if self.vertices.len() < MIN_FACE_VERTICES {
            return Err(FaceError::InsufficientVertices {
                count: self.vertices.len(),
                minimum: MIN_FACE_VERTICES,
            });
        }

        // Check for duplicate vertices
        if let Some((vertex, positions)) = self.find_duplicate_vertex() {
            return Err(FaceError::DuplicateVertex { vertex, positions });
        }

        Ok(())
    }

    /// Check if face has duplicate vertex indices.
    ///
    /// # Examples
    ///
    /// ```
    /// use cfd_mesh::topology::Face;
    ///
    /// let valid = Face::triangle(0, 1, 2);
    /// assert!(!valid.has_duplicate_vertices());
    ///
    /// let invalid = Face::from_vertices(vec![0, 1, 1, 2]);
    /// assert!(invalid.has_duplicate_vertices());
    /// ```
    #[must_use]
    pub fn has_duplicate_vertices(&self) -> bool {
        self.find_duplicate_vertex().is_some()
    }

    /// Find the first duplicate vertex and its positions.
    ///
    /// Returns `None` if all vertices are unique.
    /// Uses BTreeMap for deterministic ordering (smallest vertex index first).
    fn find_duplicate_vertex(&self) -> Option<(usize, Vec<usize>)> {
        let mut seen: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for (pos, &vertex) in self.vertices.iter().enumerate() {
            seen.entry(vertex).or_default().push(pos);
        }
        // BTreeMap iterates in sorted key order - deterministic
        for (vertex, positions) in seen {
            if positions.len() > 1 {
                return Some((vertex, positions));
            }
        }
        None
    }

    /// Number of vertices in the face.
    ///
    /// # Examples
    ///
    /// ```
    /// use cfd_mesh::topology::Face;
    ///
    /// assert_eq!(Face::triangle(0, 1, 2).vertex_count(), 3);
    /// assert_eq!(Face::quad(0, 1, 2, 3).vertex_count(), 4);
    /// ```
    #[must_use]
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }

    /// Check if face contains a specific vertex.
    ///
    /// # Examples
    ///
    /// ```
    /// use cfd_mesh::topology::Face;
    ///
    /// let face = Face::triangle(0, 1, 2);
    /// assert!(face.contains(1));
    /// assert!(!face.contains(5));
    /// ```
    #[must_use]
    pub fn contains(&self, vertex: usize) -> bool {
        self.vertices.contains(&vertex)
    }

    /// Get edges of the face as vertex index tuples.
    ///
    /// Returns edges in order around the face perimeter, where each edge
    /// connects consecutive vertices (with wraparound for the last edge).
    ///
    /// # Complexity
    ///
    /// O(n) allocation. Use [`edges_iter`] for zero-allocation iteration.
    ///
    /// # Examples
    ///
    /// ```
    /// use cfd_mesh::topology::Face;
    ///
    /// let face = Face::triangle(0, 1, 2);
    /// let edges = face.edges();
    /// assert_eq!(edges, vec![(0, 1), (1, 2), (2, 0)]);
    /// ```
    #[must_use]
    pub fn edges(&self) -> Vec<(usize, usize)> {
        self.edges_iter().collect()
    }

    /// Iterate over edges without allocation.
    ///
    /// Returns edges as (start, end) tuples in winding order around the face.
    /// Preserves edge direction (not canonically ordered).
    ///
    /// # Complexity
    ///
    /// O(n) iteration, zero allocation.
    ///
    /// # Examples
    ///
    /// ```
    /// use cfd_mesh::topology::Face;
    ///
    /// let face = Face::triangle(0, 1, 2);
    /// let edge_count = face.edges_iter().count();
    /// assert_eq!(edge_count, 3);
    ///
    /// // Find edge containing vertex 1
    /// let edges_with_v1: Vec<_> = face.edges_iter()
    ///     .filter(|(a, b)| *a == 1 || *b == 1)
    ///     .collect();
    /// assert_eq!(edges_with_v1.len(), 2);
    /// ```
    pub fn edges_iter(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        let n = self.vertices.len();
        (0..n).map(move |i| (self.vertices[i], self.vertices[(i + 1) % n]))
    }

    /// Get edges as canonically-ordered Edge structs.
    ///
    /// Unlike [`edges`], this returns [`Edge`] structs with canonical vertex
    /// ordering (smaller index first), enabling consistent hashing and comparison.
    ///
    /// # Examples
    ///
    /// ```
    /// use cfd_mesh::topology::{Face, Edge};
    ///
    /// let face = Face::triangle(2, 0, 1);
    /// let edges: Vec<Edge> = face.canonical_edges().collect();
    ///
    /// // Edges are canonically ordered (smaller vertex first)
    /// assert!(edges.contains(&Edge::new(0, 2)));
    /// assert!(edges.contains(&Edge::new(0, 1)));
    /// assert!(edges.contains(&Edge::new(1, 2)));
    /// ```
    pub fn canonical_edges(&self) -> impl Iterator<Item = Edge> + '_ {
        let n = self.vertices.len();
        (0..n).map(move |i| {
            let v1 = self.vertices[i];
            let v2 = self.vertices[(i + 1) % n];
            Edge::new(v1, v2)
        })
    }

    /// Find a shared edge with another face, if one exists.
    ///
    /// Returns the shared edge in canonical form if the faces share exactly
    /// two consecutive vertices (an edge).
    ///
    /// # Examples
    ///
    /// ```
    /// use cfd_mesh::topology::{Face, Edge};
    ///
    /// let face1 = Face::triangle(0, 1, 2);
    /// let face2 = Face::triangle(1, 2, 3);
    ///
    /// let shared = face1.shares_edge_with(&face2);
    /// assert_eq!(shared, Some(Edge::new(1, 2)));
    ///
    /// let face3 = Face::triangle(4, 5, 6);
    /// assert_eq!(face1.shares_edge_with(&face3), None);
    /// ```
    #[must_use]
    pub fn shares_edge_with(&self, other: &Face) -> Option<Edge> {
        let self_edges: std::collections::HashSet<_> = self.canonical_edges().collect();
        other.canonical_edges().find(|e| self_edges.contains(e))
    }

    /// Get vertices shared with another face.
    ///
    /// # Complexity
    ///
    /// O(n + m) where n and m are the vertex counts of self and other.
    ///
    /// # Examples
    ///
    /// ```
    /// use cfd_mesh::topology::Face;
    ///
    /// let face1 = Face::triangle(0, 1, 2);
    /// let face2 = Face::triangle(1, 2, 3);
    ///
    /// let shared = face1.shared_vertices(&face2);
    /// assert_eq!(shared.len(), 2);
    /// assert!(shared.contains(&1));
    /// assert!(shared.contains(&2));
    /// ```
    #[must_use]
    pub fn shared_vertices(&self, other: &Face) -> Vec<usize> {
        let other_set: HashSet<_> = other.vertices.iter().copied().collect();
        self.vertices
            .iter()
            .filter(|v| other_set.contains(v))
            .copied()
            .collect()
    }

    /// Get adjacent vertices for a given vertex in the face.
    ///
    /// Returns `(previous, next)` vertices in the face's winding order.
    /// Returns `None` if the vertex is not in the face.
    ///
    /// # Examples
    ///
    /// ```
    /// use cfd_mesh::topology::Face;
    ///
    /// let face = Face::triangle(0, 1, 2);
    ///
    /// assert_eq!(face.adjacent_vertices(1), Some((0, 2)));
    /// assert_eq!(face.adjacent_vertices(0), Some((2, 1)));
    /// assert_eq!(face.adjacent_vertices(5), None);
    /// ```
    #[must_use]
    pub fn adjacent_vertices(&self, vertex: usize) -> Option<(usize, usize)> {
        let pos = self.vertices.iter().position(|&v| v == vertex)?;
        let n = self.vertices.len();
        let prev = self.vertices[(pos + n - 1) % n];
        let next = self.vertices[(pos + 1) % n];
        Some((prev, next))
    }

    /// Create a new face with reversed winding order.
    ///
    /// Reversing the winding order flips the face normal direction.
    /// Used for correcting face orientation in mesh processing.
    ///
    /// # Mathematical Note
    ///
    /// If original normal **n** = (**v₁**-**v₀**) × (**v₂**-**v₀**),
    /// reversed face has normal **n'** = -**n**.
    ///
    /// # Examples
    ///
    /// ```
    /// use cfd_mesh::topology::Face;
    ///
    /// let face = Face::triangle(0, 1, 2);
    /// let reversed = face.winding_reversed();
    ///
    /// assert_eq!(reversed.vertices, vec![2, 1, 0]);
    /// ```
    #[must_use]
    pub fn winding_reversed(&self) -> Self {
        let mut reversed = self.vertices.clone();
        reversed.reverse();
        Self {
            vertices: reversed,
            global_id: self.global_id,
            partition_id: self.partition_id,
        }
    }

    /// Check if face is a triangle (3 vertices).
    ///
    /// # Examples
    ///
    /// ```
    /// use cfd_mesh::topology::Face;
    ///
    /// assert!(Face::triangle(0, 1, 2).is_triangle());
    /// assert!(!Face::quad(0, 1, 2, 3).is_triangle());
    /// ```
    #[must_use]
    pub fn is_triangle(&self) -> bool {
        self.vertices.len() == 3
    }

    /// Check if face is a quadrilateral (4 vertices).
    ///
    /// # Examples
    ///
    /// ```
    /// use cfd_mesh::topology::Face;
    ///
    /// assert!(Face::quad(0, 1, 2, 3).is_quad());
    /// assert!(!Face::triangle(0, 1, 2).is_quad());
    /// ```
    #[must_use]
    pub fn is_quad(&self) -> bool {
        self.vertices.len() == 4
    }

    /// Check if this is a valid face (quick validation).
    ///
    /// Equivalent to `self.validate().is_ok()` but more ergonomic
    /// for conditional checks.
    ///
    /// # Examples
    ///
    /// ```
    /// use cfd_mesh::topology::Face;
    ///
    /// assert!(Face::triangle(0, 1, 2).is_valid());
    /// assert!(!Face::from_vertices(vec![0, 1]).is_valid());
    /// ```
    #[must_use]
    pub fn is_valid(&self) -> bool {
        self.validate().is_ok()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ========================================================================
    // Constructor Tests
    // ========================================================================

    // TODO: Implement comprehensive face topology validation testing framework
    // DEPENDENCIES: Add rigorous validation for mesh topology consistency and geometric correctness
    // BLOCKED BY: Limited understanding of mesh topology validation requirements and edge cases
    // PRIORITY: Medium - Important for ensuring mesh generation robustness

    #[test]
    fn test_triangle_constructor() {
        let face = Face::triangle(0, 1, 2);
        assert_eq!(face.vertices, vec![0, 1, 2]);
        assert_eq!(face.global_id, None);
        assert_eq!(face.partition_id, None);
        assert!(face.is_triangle());
        assert!(!face.is_quad());
    }

    #[test]
    fn test_quad_constructor() {
        let face = Face::quad(0, 1, 2, 3);
        assert_eq!(face.vertices, vec![0, 1, 2, 3]);
        assert!(face.is_quad());
        assert!(!face.is_triangle());
    }

    #[test]
    fn test_from_vertices_pentagon() {
        let face = Face::from_vertices(vec![0, 1, 2, 3, 4]);
        assert_eq!(face.vertex_count(), 5);
        assert!(!face.is_triangle());
        assert!(!face.is_quad());
    }

    #[test]
    fn test_try_new_valid() {
        let face = Face::try_new(vec![0, 1, 2]).unwrap();
        assert!(face.is_triangle());
    }

    #[test]
    fn test_try_new_insufficient_vertices() {
        let err = Face::try_new(vec![0, 1]).unwrap_err();
        assert!(matches!(
            err,
            FaceError::InsufficientVertices {
                count: 2,
                minimum: 3
            }
        ));
    }

    #[test]
    fn test_try_new_empty() {
        let err = Face::try_new(vec![]).unwrap_err();
        assert!(matches!(
            err,
            FaceError::InsufficientVertices {
                count: 0,
                minimum: 3
            }
        ));
    }

    #[test]
    fn test_try_new_duplicate_vertex() {
        let err = Face::try_new(vec![0, 1, 1]).unwrap_err();
        match err {
            FaceError::DuplicateVertex { vertex, positions } => {
                assert_eq!(vertex, 1);
                assert_eq!(positions, vec![1, 2]);
            }
            _ => panic!("Expected DuplicateVertex error"),
        }
    }

    #[test]
    fn test_try_triangle_valid() {
        let face = Face::try_triangle(0, 1, 2).unwrap();
        assert_eq!(face.vertices, vec![0, 1, 2]);
    }

    #[test]
    fn test_try_triangle_duplicate() {
        let err = Face::try_triangle(0, 0, 1).unwrap_err();
        assert!(matches!(err, FaceError::DuplicateVertex { .. }));
    }

    #[test]
    fn test_try_quad_valid() {
        let face = Face::try_quad(0, 1, 2, 3).unwrap();
        assert_eq!(face.vertices, vec![0, 1, 2, 3]);
    }

    #[test]
    fn test_try_quad_duplicate() {
        let err = Face::try_quad(0, 1, 2, 0).unwrap_err();
        assert!(matches!(err, FaceError::DuplicateVertex { .. }));
    }

    // ========================================================================
    // Distributed Mesh Support Tests
    // ========================================================================

    #[test]
    fn test_distributed_info() {
        let face = Face::triangle(0, 1, 2).with_distributed_info(100, 2);
        assert_eq!(face.global_id, Some(100));
        assert_eq!(face.partition_id, Some(2));
    }

    // ========================================================================
    // Validation Tests
    // ========================================================================

    #[test]
    fn test_validate_valid_triangle() {
        let face = Face::triangle(0, 1, 2);
        assert!(face.validate().is_ok());
    }

    #[test]
    fn test_validate_valid_quad() {
        let face = Face::quad(0, 1, 2, 3);
        assert!(face.validate().is_ok());
    }

    #[test]
    fn test_validate_insufficient_vertices() {
        let face = Face::from_vertices(vec![0, 1]);
        let err = face.validate().unwrap_err();
        assert!(matches!(err, FaceError::InsufficientVertices { .. }));
    }

    #[test]
    fn test_validate_single_vertex() {
        let face = Face::from_vertices(vec![0]);
        assert!(!face.is_valid());
    }

    #[test]
    fn test_has_duplicate_vertices_true() {
        let face = Face::from_vertices(vec![0, 1, 1, 2]);
        assert!(face.has_duplicate_vertices());
    }

    #[test]
    fn test_has_duplicate_vertices_false() {
        let face = Face::triangle(0, 1, 2);
        assert!(!face.has_duplicate_vertices());
    }

    #[test]
    fn test_duplicate_detection_deterministic() {
        // Multiple duplicates: vertex 1 appears at [1,2], vertex 3 appears at [3,4]
        // BTreeMap ensures smallest vertex (1) is reported first
        let face = Face::from_vertices(vec![0, 1, 1, 3, 3, 5]);
        let err = face.validate().unwrap_err();
        match err {
            FaceError::DuplicateVertex { vertex, positions } => {
                assert_eq!(vertex, 1); // Smallest duplicate vertex
                assert_eq!(positions, vec![1, 2]);
            }
            _ => panic!("Expected DuplicateVertex error"),
        }
    }

    #[test]
    fn test_all_same_vertex() {
        // Edge case: all vertices are the same
        let face = Face::from_vertices(vec![0, 0, 0]);
        let err = face.validate().unwrap_err();
        match err {
            FaceError::DuplicateVertex { vertex, positions } => {
                assert_eq!(vertex, 0);
                assert_eq!(positions, vec![0, 1, 2]);
            }
            _ => panic!("Expected DuplicateVertex error"),
        }
    }

    #[test]
    fn test_is_valid_true() {
        let face = Face::triangle(0, 1, 2);
        assert!(face.is_valid());
    }

    #[test]
    fn test_is_valid_false() {
        let face = Face::from_vertices(vec![0, 1]);
        assert!(!face.is_valid());
    }

    // ========================================================================
    // Edge Tests
    // ========================================================================

    #[test]
    fn test_edges_triangle() {
        let face = Face::triangle(0, 1, 2);
        let edges = face.edges();
        assert_eq!(edges, vec![(0, 1), (1, 2), (2, 0)]);
    }

    #[test]
    fn test_edges_quad() {
        let face = Face::quad(0, 1, 2, 3);
        let edges = face.edges();
        assert_eq!(edges, vec![(0, 1), (1, 2), (2, 3), (3, 0)]);
    }

    #[test]
    fn test_edges_iter_zero_allocation() {
        let face = Face::triangle(0, 1, 2);
        // edges_iter returns iterator, no allocation until collect
        let count = face.edges_iter().count();
        assert_eq!(count, 3);
    }

    #[test]
    fn test_edges_iter_equivalence() {
        let face = Face::quad(0, 1, 2, 3);
        let edges_vec = face.edges();
        let edges_iter: Vec<_> = face.edges_iter().collect();
        assert_eq!(edges_vec, edges_iter);
    }

    #[test]
    fn test_edges_iter_filter() {
        let face = Face::triangle(0, 1, 2);
        // Filter edges containing vertex 1
        let edges_with_v1: Vec<_> = face
            .edges_iter()
            .filter(|(a, b)| *a == 1 || *b == 1)
            .collect();
        assert_eq!(edges_with_v1, vec![(0, 1), (1, 2)]);
    }

    #[test]
    fn test_canonical_edges_ordering() {
        // Vertices in descending order
        let face = Face::triangle(2, 1, 0);
        let edges: Vec<Edge> = face.canonical_edges().collect();

        // All edges should have smaller index first
        for edge in &edges {
            assert!(edge.start <= edge.end);
        }

        // Should contain edges (1,2), (0,1), (0,2)
        assert!(edges.contains(&Edge::new(1, 2)));
        assert!(edges.contains(&Edge::new(0, 1)));
        assert!(edges.contains(&Edge::new(0, 2)));
    }

    #[test]
    fn test_canonical_edges_count() {
        let face = Face::quad(0, 1, 2, 3);
        assert_eq!(face.canonical_edges().count(), 4);
    }

    // ========================================================================
    // Topology Query Tests
    // ========================================================================

    #[test]
    fn test_contains_vertex() {
        let face = Face::triangle(0, 1, 2);
        assert!(face.contains(0));
        assert!(face.contains(1));
        assert!(face.contains(2));
        assert!(!face.contains(3));
    }

    #[test]
    fn test_vertex_count() {
        assert_eq!(Face::triangle(0, 1, 2).vertex_count(), 3);
        assert_eq!(Face::quad(0, 1, 2, 3).vertex_count(), 4);
        assert_eq!(
            Face::from_vertices(vec![0, 1, 2, 3, 4, 5]).vertex_count(),
            6
        );
    }

    #[test]
    fn test_shares_edge_with_adjacent() {
        let face1 = Face::triangle(0, 1, 2);
        let face2 = Face::triangle(1, 2, 3);

        let shared = face1.shares_edge_with(&face2);
        assert_eq!(shared, Some(Edge::new(1, 2)));
    }

    #[test]
    fn test_shares_edge_with_none() {
        let face1 = Face::triangle(0, 1, 2);
        let face2 = Face::triangle(3, 4, 5);

        assert_eq!(face1.shares_edge_with(&face2), None);
    }

    #[test]
    fn test_shares_edge_with_single_vertex_shared() {
        let face1 = Face::triangle(0, 1, 2);
        let face2 = Face::triangle(2, 3, 4);

        // Only vertex 2 is shared, not an edge
        assert_eq!(face1.shares_edge_with(&face2), None);
    }

    #[test]
    fn test_shared_vertices() {
        let face1 = Face::triangle(0, 1, 2);
        let face2 = Face::triangle(1, 2, 3);

        let shared = face1.shared_vertices(&face2);
        assert_eq!(shared.len(), 2);
        assert!(shared.contains(&1));
        assert!(shared.contains(&2));
    }

    #[test]
    fn test_shared_vertices_none() {
        let face1 = Face::triangle(0, 1, 2);
        let face2 = Face::triangle(3, 4, 5);

        let shared = face1.shared_vertices(&face2);
        assert!(shared.is_empty());
    }

    #[test]
    fn test_adjacent_vertices() {
        let face = Face::triangle(0, 1, 2);

        assert_eq!(face.adjacent_vertices(0), Some((2, 1)));
        assert_eq!(face.adjacent_vertices(1), Some((0, 2)));
        assert_eq!(face.adjacent_vertices(2), Some((1, 0)));
        assert_eq!(face.adjacent_vertices(5), None);
    }

    #[test]
    fn test_adjacent_vertices_quad() {
        let face = Face::quad(0, 1, 2, 3);

        assert_eq!(face.adjacent_vertices(0), Some((3, 1)));
        assert_eq!(face.adjacent_vertices(2), Some((1, 3)));
    }

    // ========================================================================
    // Winding Order Tests
    // ========================================================================

    #[test]
    fn test_winding_reversed() {
        let face = Face::triangle(0, 1, 2);
        let reversed = face.winding_reversed();

        assert_eq!(reversed.vertices, vec![2, 1, 0]);
    }

    #[test]
    fn test_winding_reversed_preserves_metadata() {
        let face = Face::triangle(0, 1, 2).with_distributed_info(100, 2);
        let reversed = face.winding_reversed();

        assert_eq!(reversed.global_id, Some(100));
        assert_eq!(reversed.partition_id, Some(2));
    }

    #[test]
    fn test_winding_reversed_involutory() {
        let face = Face::triangle(0, 1, 2);
        let double_reversed = face.winding_reversed().winding_reversed();

        assert_eq!(face.vertices, double_reversed.vertices);
    }

    #[test]
    fn test_winding_reversed_quad() {
        let face = Face::quad(0, 1, 2, 3);
        let reversed = face.winding_reversed();

        assert_eq!(reversed.vertices, vec![3, 2, 1, 0]);
    }

    // ========================================================================
    // Equality and Hash Tests
    // ========================================================================

    #[test]
    fn test_face_equality() {
        let face1 = Face::triangle(0, 1, 2).with_distributed_info(100, 2);
        let face2 = Face::triangle(0, 1, 2).with_distributed_info(100, 2);
        let face3 = Face::triangle(0, 1, 3).with_distributed_info(100, 2);

        assert_eq!(face1, face2);
        assert_ne!(face1, face3);
    }

    #[test]
    fn test_face_hash_consistency() {
        use std::collections::hash_map::DefaultHasher;
        use std::hash::{Hash, Hasher};

        let face1 = Face::triangle(0, 1, 2);
        let face2 = Face::triangle(0, 1, 2);

        let mut hasher1 = DefaultHasher::new();
        let mut hasher2 = DefaultHasher::new();
        face1.hash(&mut hasher1);
        face2.hash(&mut hasher2);

        assert_eq!(hasher1.finish(), hasher2.finish());
    }

    // ========================================================================
    // Edge Case Tests
    // ========================================================================

    #[test]
    fn test_large_vertex_indices() {
        let face = Face::triangle(usize::MAX - 2, usize::MAX - 1, usize::MAX);
        assert!(face.is_valid());
        assert_eq!(face.vertex_count(), 3);
    }

    #[test]
    fn test_many_vertices() {
        let vertices: Vec<usize> = (0..100).collect();
        let face = Face::try_new(vertices).unwrap();
        assert_eq!(face.vertex_count(), 100);
        assert_eq!(face.canonical_edges().count(), 100);
    }
}
