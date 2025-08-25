//! Sparse matrix utilities and efficient assembly operations.
//!
//! This module provides zero-copy sparse matrix operations optimized for CFD applications
//! with support for parallel assembly and iterator-based construction.

mod builder;
mod patterns;
mod operations;
mod assembly;

pub use builder::{MatrixEntry, SparseMatrixBuilder};
pub use patterns::SparsePatterns;
pub use operations::SparseMatrixExt;
pub use assembly::ParallelAssembly;

// Re-export the core sparse matrix type
pub use nalgebra_sparse::CsrMatrix as SparseMatrix;

#[cfg(test)]
mod tests;