//! Element types and operations for mesh operations domain

use serde::{Deserialize, Serialize};

/// Element types - Single Source of Truth
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ElementType {
    // 1D Elements
    /// Line element (2 nodes)
    Line,
    /// Quadratic line element (3 nodes)
    Line3,

    // 2D Elements
    /// Triangle element (3 nodes)
    Triangle,
    /// Quadratic triangle element (6 nodes)
    Triangle6,
    /// Quadrilateral element (4 nodes)
    Quadrilateral,
    /// Quadratic quadrilateral element (9 nodes)
    Quadrilateral9,

    // 3D Elements
    /// Tetrahedral element (4 nodes)
    Tetrahedron,
    /// Quadratic tetrahedral element (10 nodes)
    Tetrahedron10,
    /// Hexahedral element (8 nodes)
    Hexahedron,
    /// Quadratic hexahedral element (20 nodes)
    Hexahedron20,
    /// Pyramid element (5 nodes)
    Pyramid,
    /// Prism/Wedge element (6 nodes)
    Prism,
}

impl ElementType {
    /// Get the number of nodes for this element type
    #[must_use] pub fn num_nodes(&self) -> usize {
        match self {
            Self::Line => 2,
            Self::Line3 => 3,
            Self::Triangle => 3,
            Self::Triangle6 => 6,
            Self::Quadrilateral => 4,
            Self::Quadrilateral9 => 9,
            Self::Tetrahedron => 4,
            Self::Tetrahedron10 => 10,
            Self::Hexahedron => 8,
            Self::Hexahedron20 => 20,
            Self::Pyramid => 5,
            Self::Prism => 6,
        }
    }

    /// Get the dimension of this element type
    #[must_use] pub fn dimension(&self) -> usize {
        match self {
            Self::Line | Self::Line3 => 1,
            Self::Triangle | Self::Triangle6 | Self::Quadrilateral | Self::Quadrilateral9 => 2,
            Self::Tetrahedron
            | Self::Tetrahedron10
            | Self::Hexahedron
            | Self::Hexahedron20
            | Self::Pyramid
            | Self::Prism => 3,
        }
    }

    /// Check if this is a linear element type
    #[must_use] pub fn is_linear(&self) -> bool {
        matches!(
            self,
            Self::Line
                | Self::Triangle
                | Self::Quadrilateral
                | Self::Tetrahedron
                | Self::Hexahedron
                | Self::Pyramid
                | Self::Prism
        )
    }

    /// Check if this is a quadratic element type
    #[must_use] pub fn is_quadratic(&self) -> bool {
        matches!(
            self,
            Self::Line3
                | Self::Triangle6
                | Self::Quadrilateral9
                | Self::Tetrahedron10
                | Self::Hexahedron20
        )
    }
}
