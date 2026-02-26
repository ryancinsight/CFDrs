//! Sparse matrix utilities and efficient assembly operations.
//!
//! This module provides zero-copy sparse matrix operations optimized for CFD applications
//! with support for parallel assembly and iterator-based construction.
//!
//! ## Theorem — CSR SpMV Complexity
//!
//! **Theorem (CSR SpMV Time & Space Complexity)**: For a sparse matrix A ∈ ℝ^{m×n}
//! stored in Compressed Sparse Row (CSR) format with nnz non-zeros:
//!
//! - **Time**: O(nnz) per matrix-vector product Ax — each non-zero contributes
//!   exactly one multiply-add to the output, with no wasted work.
//! - **Space**: O(n + nnz) — row pointer array of length n+1, column index and
//!   value arrays of length nnz.
//!
//! **Theorem (CSR vs Dense)**: For matrices with density d = nnz/(m·n) < 1/2,
//! CSR SpMV is faster than dense GEMV, which costs O(m·n) time and O(m·n) space.
//! CFD discretization matrices typically have d = O(1/n) (e.g., 5-7 non-zeros per row
//! for FD/FVM stencils), giving O(n) sparse vs O(n²) dense.
//!
//! ## Theorem — Sparse Cholesky Fill Bound
//!
//! **Theorem (Cholesky Fill-in, Rose 1972)**: The fill-in of sparse Cholesky
//! factorization is bounded by the elimination tree structure. For a 1D Laplacian
//! tridiagonal matrix, no fill-in occurs. For general unstructured meshes, bandwidth
//! reduction via Reverse Cuthill-McKee (RCM) ordering minimises fill-in.
//!
//! ## Theorem — Parallel Assembly Commutativity
//!
//! **Theorem (Atomic Accumulation)**: Parallel element-wise assembly of the global
//! stiffness matrix via atomic fetch-add converges to the same result as sequential
//! assembly, regardless of thread ordering, by commutativity and associativity of
//! floating-point addition (up to non-deterministic rounding order).
//!
//! ## Invariants
//! - CSR row pointer array is strictly monotonically non-decreasing.
//! - All column indices within a row are distinct and sorted in ascending order.
//! - Diagonal entries are explicitly stored for preconditioner construction.

mod assembly;
mod builder;
mod operations;
mod patterns;

pub use assembly::ParallelAssembly;
pub use builder::{MatrixEntry, SparseMatrixBuilder};
pub use operations::{sparse_sparse_mul, spmv, spmv_parallel, SparseMatrixExt};
pub use patterns::SparsePatterns;

// Re-export the core sparse matrix type
pub use nalgebra_sparse::CsrMatrix as SparseMatrix;

#[cfg(test)]
mod tests;
