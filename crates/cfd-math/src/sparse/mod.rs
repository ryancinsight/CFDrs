//! Sparse matrix utilities and efficient assembly operations.
//!
//! This module provides zero-copy sparse matrix operations optimized for CFD applications
//! with support for parallel assembly and iterator-based construction.

mod assembly;
mod builder;
mod operations;
mod patterns;

pub use assembly::ParallelAssembly;
pub use builder::{MatrixEntry, SparseMatrixBuilder};
pub use operations::{spmv, spmv_parallel, spmv_f32_simd, SparseMatrixExt};
pub use patterns::SparsePatterns;

// Re-export the core sparse matrix type
pub use nalgebra_sparse::CsrMatrix as SparseMatrix;

#[cfg(test)]
mod tests;
