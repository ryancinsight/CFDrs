//! Vertex type for the PSLG / Delaunay triangulation.
//!
//! Each vertex carries a 2-D position suitable for exact Shewchuk predicates.
//!
//! # Design
//!
//! `PslgVertexId` is a strongly-typed `u32` newtype following the same pattern
//! as [`crate::domain::core::index::VertexId`].  It indexes into the flat vertex array
//! stored in [`super::graph::Pslg`].

use std::fmt;

use crate::domain::core::scalar::Real;

/// A 2-D vertex position for the Delaunay triangulation.
///
/// Stored contiguously in the PSLG vertex pool.  Coordinates are `Real`
/// (default `f64`) to maintain full precision in all exact-predicate calls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PslgVertex {
    /// X-coordinate.
    pub x: Real,
    /// Y-coordinate.
    pub y: Real,
}

impl PslgVertex {
    /// Create a new vertex at `(x, y)`.
    #[inline]
    pub fn new(x: Real, y: Real) -> Self {
        Self { x, y }
    }

    /// Squared Euclidean distance to another vertex.
    #[inline]
    pub fn dist_sq(&self, other: &Self) -> Real {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        dx * dx + dy * dy
    }

    /// Euclidean distance to another vertex.
    #[inline]
    pub fn dist(&self, other: &Self) -> Real {
        self.dist_sq(other).sqrt()
    }

    /// Midpoint between `self` and `other`.
    #[inline]
    pub fn midpoint(&self, other: &Self) -> Self {
        Self {
            x: 0.5 * (self.x + other.x),
            y: 0.5 * (self.y + other.y),
        }
    }

    /// Convert to a `nalgebra::Point2<Real>` for predicate calls.
    #[inline]
    pub fn to_point2(&self) -> nalgebra::Point2<Real> {
        nalgebra::Point2::new(self.x, self.y)
    }
}

impl From<[Real; 2]> for PslgVertex {
    #[inline]
    fn from(arr: [Real; 2]) -> Self {
        Self {
            x: arr[0],
            y: arr[1],
        }
    }
}

impl From<(Real, Real)> for PslgVertex {
    #[inline]
    fn from((x, y): (Real, Real)) -> Self {
        Self { x, y }
    }
}

/// Strongly-typed vertex index into the PSLG vertex array.
///
/// Follows the same newtype-index pattern as [`crate::domain::core::index::VertexId`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct PslgVertexId(pub u32);

impl PslgVertexId {
    /// Create from a raw `u32`.
    #[inline]
    pub fn new(raw: u32) -> Self {
        Self(raw)
    }

    /// Create from `usize`.
    #[inline]
    pub fn from_usize(n: usize) -> Self {
        Self(n as u32)
    }

    /// Raw `u32` index.
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

impl From<usize> for PslgVertexId {
    #[inline]
    fn from(n: usize) -> Self {
        Self(n as u32)
    }
}

impl fmt::Display for PslgVertexId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "v{}", self.0)
    }
}

/// Sentinel value representing "no vertex" / the super-triangle ghost vertex.
pub const GHOST_VERTEX: PslgVertexId = PslgVertexId(u32::MAX);
