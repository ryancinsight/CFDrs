//! Coarsening strategies for AMG hierarchy construction.
//!
//! Provides multiple coarsening algorithms (Ruge-Stüben, aggregation, Falgout,
//! PMIS, HMIS) and quality analysis via algebraic distance measures.

mod algorithms;
mod quality;

pub use algorithms::*;
pub use quality::*;

use super::SparseMatrix;
use eunomia::RealField;

/// Result of coarsening operation
#[derive(Debug, Clone)]
pub struct CoarseningResult<T: RealField + Copy> {
    /// Indices of coarse points (C-points)
    pub coarse_points: Vec<usize>,
    /// Mapping from fine points to coarse points (None for F-points)
    pub fine_to_coarse_map: Vec<Option<usize>>,
    /// Strength of connection matrix
    pub strength_matrix: SparseMatrix<T>,
}

#[cfg(test)]
mod tests;
