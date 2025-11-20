//! Edge connectivity between vertices

use serde::{Deserialize, Serialize};

/// Edge connecting two vertices
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Edge {
    /// Start vertex index
    pub start: usize,
    /// End vertex index
    pub end: usize,
}

impl Edge {
    /// Create a new edge with canonical vertex ordering
    #[must_use]
    pub fn new(v1: usize, v2: usize) -> Self {
        // Canonical ordering for consistent hashing
        if v1 < v2 {
            Self { start: v1, end: v2 }
        } else {
            Self { start: v2, end: v1 }
        }
    }

    /// Get the vertices of this edge
    #[must_use]
    pub fn vertices(&self) -> [usize; 2] {
        [self.start, self.end]
    }

    /// Check if edge contains a vertex
    #[must_use]
    pub fn contains(&self, vertex: usize) -> bool {
        self.start == vertex || self.end == vertex
    }

    /// Get the other vertex of an edge given one vertex
    #[must_use]
    pub fn other(&self, vertex: usize) -> Option<usize> {
        if self.start == vertex {
            Some(self.end)
        } else if self.end == vertex {
            Some(self.start)
        } else {
            None
        }
    }
}
