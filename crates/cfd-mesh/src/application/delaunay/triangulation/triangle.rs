//! Triangle data structure for the Delaunay triangulation.
//!
//! Each triangle stores three vertex indices (CCW order) and three adjacency
//! links to neighbouring triangles.  The adjacency convention is:
//!
//! ```text
//!        v[0]
//!       / \
//!  e[2]/   \e[0]        adj[i] = triangle across edge opposite v[i]
//!     /     \
//! v[2]───────v[1]
//!       e[1]
//! ```
//!
//! - Edge `e[0]` is `v[1]→v[2]` (opposite `v[0]`), adjacent to `adj[0]`.
//! - Edge `e[1]` is `v[2]→v[0]` (opposite `v[1]`), adjacent to `adj[1]`.
//! - Edge `e[2]` is `v[0]→v[1]` (opposite `v[2]`), adjacent to `adj[2]`.

use std::fmt;

use crate::application::delaunay::pslg::vertex::PslgVertexId;

/// Strongly-typed triangle index.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct TriangleId(pub u32);

impl TriangleId {
    /// Create from raw `u32`.
    #[inline]
    pub fn new(raw: u32) -> Self {
        Self(raw)
    }

    /// Create from `usize`.
    #[inline]
    pub fn from_usize(n: usize) -> Self {
        Self(n as u32)
    }

    /// Raw index.
    #[inline]
    pub fn raw(self) -> u32 {
        self.0
    }

    /// As `usize`.
    #[inline]
    pub fn idx(self) -> usize {
        self.0 as usize
    }
}

impl From<usize> for TriangleId {
    #[inline]
    fn from(n: usize) -> Self {
        Self(n as u32)
    }
}

impl fmt::Display for TriangleId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "t{}", self.0)
    }
}

/// Sentinel value: no adjacent triangle (boundary).
pub const GHOST_TRIANGLE: TriangleId = TriangleId(u32::MAX);

/// A triangle in the triangulation.
///
/// # Invariants
///
/// - `vertices[0..3]` are in CCW order (verified via `orient_2d`).
/// - `adj[i]` is the triangle sharing edge opposite `vertices[i]`,
///   or [`GHOST_TRIANGLE`] if that edge is on the convex hull boundary.
/// - If `alive` is `false`, the triangle has been logically deleted
///   (tombstone) and will be skipped during iteration.
#[derive(Clone, Debug)]
pub struct Triangle {
    /// Three vertex indices in counter-clockwise order.
    pub vertices: [PslgVertexId; 3],
    /// Adjacent triangle across each edge. `adj[i]` is across the edge
    /// opposite `vertices[i]`.
    pub adj: [TriangleId; 3],
    /// Whether this triangle is alive (not tombstoned).
    pub alive: bool,
    /// Whether each edge is a constraint (PSLG segment).
    /// `constrained[i]` corresponds to edge opposite `vertices[i]`.
    pub constrained: [bool; 3],
}

impl Triangle {
    /// Create a new triangle with the given vertices and no adjacency.
    pub fn new(v0: PslgVertexId, v1: PslgVertexId, v2: PslgVertexId) -> Self {
        Self {
            vertices: [v0, v1, v2],
            adj: [GHOST_TRIANGLE; 3],
            alive: true,
            constrained: [false; 3],
        }
    }

    /// Return the local edge index (0, 1, or 2) of the edge opposite vertex `v`.
    ///
    /// Returns `None` if `v` is not a vertex of this triangle.
    #[inline]
    pub fn opposite_edge(&self, v: PslgVertexId) -> Option<usize> {
        self.vertices.iter().position(|&vi| vi == v)
    }

    /// Return the two vertices of edge `i` (the edge opposite `vertices[i]`).
    #[inline]
    pub fn edge_vertices(&self, edge: usize) -> (PslgVertexId, PslgVertexId) {
        let a = self.vertices[(edge + 1) % 3];
        let b = self.vertices[(edge + 2) % 3];
        (a, b)
    }

    /// Check if this triangle contains vertex `v`.
    #[inline]
    pub fn contains_vertex(&self, v: PslgVertexId) -> bool {
        self.vertices.contains(&v)
    }

    /// Return the local index (0, 1, 2) of vertex `v`, or `None`.
    #[inline]
    pub fn vertex_index(&self, v: PslgVertexId) -> Option<usize> {
        self.vertices.iter().position(|&vi| vi == v)
    }

    /// Return the edge index shared with triangle `other_tid`.
    #[inline]
    pub fn shared_edge(&self, other_tid: TriangleId) -> Option<usize> {
        self.adj.iter().position(|&a| a == other_tid)
    }
}
