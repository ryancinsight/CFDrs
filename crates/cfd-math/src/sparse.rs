//! Sparse matrix utilities and efficient assembly operations.
//!
//! This module provides zero-copy sparse matrix operations optimized for CFD applications
//! with support for parallel assembly and iterator-based construction.

use cfd_core::error::{Error, Result};
use nalgebra::{DVector, RealField};
use nalgebra_sparse::{CooMatrix, CsrMatrix};
use rayon::prelude::*;
use std::collections::HashMap;

/// Sparse matrix wrapper with additional CFD-specific operations
pub type SparseMatrix<T> = nalgebra_sparse::CsrMatrix<T>;

/// Entry for sparse matrix assembly
#[derive(Debug, Clone, Copy)]
pub struct MatrixEntry<T: RealField + Copy> {
    /// Row index
    pub row: usize,
    /// Column index
    pub col: usize,
    /// Value
    pub value: T,
}

impl<T: RealField + Copy> MatrixEntry<T> {
    /// Create a new matrix entry
    pub fn new(row: usize, col: usize, value: T) -> Self {
        Self { row, col, value }
    }
}

/// Sparse matrix builder with efficient assembly
pub struct SparseMatrixBuilder<T: RealField + Copy> {
    rows: usize,
    cols: usize,
    entries: Vec<MatrixEntry<T>>,
    allow_duplicates: bool,
}

impl<T: RealField + Copy> SparseMatrixBuilder<T> {
    /// Create a new builder
    pub fn new(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            entries: Vec::new(),
            allow_duplicates: false,
        }
    }

    /// Create with estimated capacity
    pub fn with_capacity(rows: usize, cols: usize, capacity: usize) -> Self {
        Self {
            rows,
            cols,
            entries: Vec::with_capacity(capacity),
            allow_duplicates: false,
        }
    }

    /// Allow duplicate entries (will be summed)
    pub fn allow_duplicates(mut self, allow: bool) -> Self {
        self.allow_duplicates = allow;
        self
    }

    /// Add a single entry
    pub fn add_entry(&mut self, row: usize, col: usize, value: T) -> Result<()> {
        if row >= self.rows || col >= self.cols {
            return Err(Error::InvalidConfiguration(
                format!("Entry ({}, {}) out of bounds for {}x{} matrix", row, col, self.rows, self.cols)
            ));
        }
        self.entries.push(MatrixEntry::new(row, col, value));
        Ok(())
    }

    /// Add multiple entries using iterator
    pub fn add_entries<I>(&mut self, entries: I) -> Result<()>
    where
        I: IntoIterator<Item = MatrixEntry<T>>,
    {
        for entry in entries {
            self.add_entry(entry.row, entry.col, entry.value)?;
        }
        Ok(())
    }

    /// Add entries from triplet format
    pub fn add_triplets<I>(&mut self, triplets: I) -> Result<()>
    where
        I: IntoIterator<Item = (usize, usize, T)>,
    {
        for (row, col, value) in triplets {
            self.add_entry(row, col, value)?;
        }
        Ok(())
    }

    /// Build the sparse matrix using COO format for efficiency
    pub fn build(self) -> Result<SparseMatrix<T>> {
        if self.entries.is_empty() {
            return Ok(CsrMatrix::zeros(self.rows, self.cols));
        }

        // Use COO matrix for efficient assembly
        let mut coo = CooMatrix::new(self.rows, self.cols);

        if self.allow_duplicates {
            // Add entries directly (duplicates will be summed)
            for entry in &self.entries {
                coo.push(entry.row, entry.col, entry.value);
            }
        } else {
            // Combine duplicate entries manually for better control
            let mut entry_map: HashMap<(usize, usize), T> = HashMap::new();

            for entry in &self.entries {
                let key = (entry.row, entry.col);
                entry_map.entry(key)
                    .and_modify(|v| *v += entry.value)
                    .or_insert(entry.value);
            }

            for ((row, col), value) in entry_map {
                coo.push(row, col, value);
            }
        }

        Ok(CsrMatrix::from(&coo))
    }

    /// Build matrix in parallel for large systems
    pub fn build_parallel(self) -> Result<SparseMatrix<T>>
    where
        T: Send + Sync,
    {
        if self.entries.is_empty() {
            return Ok(CsrMatrix::zeros(self.rows, self.cols));
        }

        // Zero-copy parallel aggregation using advanced iterator patterns
        let entry_map: HashMap<(usize, usize), T> = self.entries
            .par_iter()
            .map(|entry| ((entry.row, entry.col), entry.value))
            .fold(
                HashMap::new,
                |mut acc, (key, value)| {
                    *acc.entry(key).or_insert_with(T::zero) += value;
                    acc
                }
            )
            .reduce(
                HashMap::new,
                |mut acc1, acc2| {
                    for (key, value) in acc2 {
                        *acc1.entry(key).or_insert_with(T::zero) += value;
                    }
                    acc1
                }
            );

        let mut coo = CooMatrix::new(self.rows, self.cols);
        for ((row, col), value) in entry_map {
            coo.push(row, col, value);
        }

        Ok(CsrMatrix::from(&coo))
    }

    /// Get current number of entries
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Check if builder is empty
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

/// Extension trait for sparse matrix operations
pub trait SparseMatrixExt<T: RealField + Copy> {
    /// Matrix-vector multiplication with zero-copy optimization
    fn matvec_inplace(&self, x: &DVector<T>, y: &mut DVector<T>);

    /// Compute matrix norm
    fn norm_frobenius(&self) -> T;

    /// Get matrix statistics
    fn stats(&self) -> MatrixStats<T>;

    /// Extract diagonal as vector
    fn diagonal(&self) -> DVector<T>;

    /// Check if matrix is symmetric (within tolerance)
    fn is_symmetric(&self, tolerance: T) -> bool;
    
    /// Get condition number estimate (ratio of largest to smallest eigenvalue magnitude)
    fn condition_number_estimate(&self) -> T;
}

impl<T: RealField + Copy> SparseMatrixExt<T> for CsrMatrix<T> {
    fn matvec_inplace(&self, x: &DVector<T>, y: &mut DVector<T>) {
        assert_eq!(self.ncols(), x.len());
        assert_eq!(self.nrows(), y.len());

        y.fill(T::zero());

        // Zero-copy iterator-based approach with SIMD-friendly patterns
        self.row_iter()
            .enumerate()
            .for_each(|(row_idx, row)| {
                y[row_idx] = row.col_indices()
                    .iter()
                    .zip(row.values())
                    .map(|(&col_idx, value)| *value * x[col_idx])
                    .fold(T::zero(), |acc, val| acc + val);
            });
    }

    fn norm_frobenius(&self) -> T {
        // Zero-copy computation using iterator combinators
        self.values()
            .iter()
            .map(|v| *v * *v)
            .fold(T::zero(), |acc, v| acc + v)
            .sqrt()
    }

    fn stats(&self) -> MatrixStats<T> {
        let nnz = self.nnz();
        let density = T::from_usize(nnz).unwrap_or_else(|| T::zero()) /
                     T::from_usize(self.nrows() * self.ncols()).unwrap_or_else(|| T::one());

        let (min_val, max_val) = self.values()
            .iter()
            .fold((T::max_value().expect("CRITICAL: Add proper error handling"), T::min_value().expect("CRITICAL: Add proper error handling")),
                  |(min, max), val| {
                      (min.min(*val), max.max(*val))
                  });

        MatrixStats {
            rows: self.nrows(),
            cols: self.ncols(),
            nnz,
            density,
            min_value: min_val,
            max_value: max_val,
        }
    }

    fn diagonal(&self) -> DVector<T> {
        let n = self.nrows().min(self.ncols());
        let mut diag = DVector::zeros(n);

        // Extract diagonal using row iteration
        for (row_idx, row) in self.row_iter().enumerate().take(n) {
            for (col_idx, value) in row.col_indices().iter().zip(row.values()) {
                if *col_idx == row_idx {
                    diag[row_idx] = *value;
                    break;
                }
            }
        }

        diag
    }

    fn is_symmetric(&self, tolerance: T) -> bool {
        if self.nrows() != self.ncols() {
            return false;
        }

        // Tuned symmetry check: only check non-zero entries and their symmetric counterparts
        // This leverages the sparse structure to avoid O(n²) complexity
        use std::collections::HashMap;

        // Build a map of (row, col) -> value for efficient lookup
        let mut entries: HashMap<(usize, usize), T> = HashMap::new();

        for (row_idx, row) in self.row_iter().enumerate() {
            for (col_idx, value) in row.col_indices().iter().zip(row.values()) {
                entries.insert((row_idx, *col_idx), *value);
            }
        }

        // Check symmetry only for existing entries
        for ((i, j), value) in &entries {
            // Get the symmetric entry A[j,i]
            let symmetric_value = entries.get(&(*j, *i)).cloned().unwrap_or_else(T::zero);

            // Check if A[i,j] ≈ A[j,i]
            let diff = (*value - symmetric_value).abs();
            if diff > tolerance {
                return false;
            }
        }

        true
    }

    fn condition_number_estimate(&self) -> T {
        // Use Gershgorin circle theorem for eigenvalue bounds
        let (max_eigen, min_eigen) = (0..self.nrows())
            .map(|i| {
                let diag = self.get(i, i).unwrap_or_else(T::zero);
                let radius = self.row(i)
                    .filter(|(j, _)| *j != i)
                    .map(|(_, v)| v.abs())
                    .fold(T::zero(), |acc, v| acc + v);
                (diag.abs() + radius, (diag.abs() - radius).max(T::zero()))
            })
            .fold((T::zero(), T::from_f64(1e10).unwrap_or_else(|| T::one())),
                |(max_val, min_val), (upper, lower)| {
                    (max_val.max(upper), if lower > T::zero() { min_val.min(lower) } else { min_val })
                });
        
        if min_eigen > T::zero() {
            max_eigen / min_eigen
        } else {
            T::from_f64(1e10).unwrap_or_else(|| T::one()) // Return large number for singular matrix
        }
    }
}

/// Matrix statistics
#[derive(Debug, Clone)]
pub struct MatrixStats<T: RealField + Copy> {
    /// Number of rows
    pub rows: usize,
    /// Number of columns
    pub cols: usize,
    /// Number of non-zeros
    pub nnz: usize,
    /// Matrix density (nnz / (rows * cols))
    pub density: T,
    /// Minimum value
    pub min_value: T,
    /// Maximum value
    pub max_value: T,
}

/// Utility functions for common sparse matrix patterns
pub struct SparsePatterns;

impl SparsePatterns {
    /// Create a tridiagonal matrix
    pub fn tridiagonal<T: RealField + Copy>(
        n: usize,
        lower: T,
        diagonal: T,
        upper: T,
    ) -> Result<SparseMatrix<T>> {
        let mut builder = SparseMatrixBuilder::with_capacity(n, n, 3 * n - 2);

        for i in 0..n {
            // Diagonal element
            builder.add_entry(i, i, diagonal)?;

            // Lower diagonal
            if i > 0 {
                builder.add_entry(i, i - 1, lower)?;
            }

            // Upper diagonal
            if i < n - 1 {
                builder.add_entry(i, i + 1, upper)?;
            }
        }

        builder.build()
    }

    /// Create a 5-point stencil matrix for 2D finite differences
    pub fn five_point_stencil<T: RealField + Copy>(
        nx: usize,
        ny: usize,
        dx: T,
        dy: T,
    ) -> Result<SparseMatrix<T>> {
        let n = nx * ny;
        let mut builder = SparseMatrixBuilder::with_capacity(n, n, 5 * n);

        let dx2_inv = T::one() / (dx * dx);
        let dy2_inv = T::one() / (dy * dy);
        let center_coeff = -T::from_f64(2.0).unwrap_or_else(|| T::zero()) * (dx2_inv + dy2_inv);

        // Use iterator with cartesian product for 2D grid traversal
        (0..ny).flat_map(|j| (0..nx).map(move |i| (i, j)))
            .try_for_each(|(i, j)| {
                let idx = j * nx + i;

                // Center point
                builder.add_entry(idx, idx, center_coeff)?;

                // Neighbors - compute indices only when valid
                if i > 0 {
                    builder.add_entry(idx, idx - 1, dx2_inv)?;  // Left
                }
                if i < nx - 1 {
                    builder.add_entry(idx, idx + 1, dx2_inv)?;  // Right
                }
                if j > 0 {
                    builder.add_entry(idx, idx - nx, dy2_inv)?;  // Bottom
                }
                if j < ny - 1 {
                    builder.add_entry(idx, idx + nx, dy2_inv)?;  // Top
                }
                Ok::<(), cfd_core::error::Error>(())
            })?;

        builder.build()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use nalgebra::{DVector, DMatrix};

    // Helper function to convert sparse matrix to dense for testing
    fn sparse_to_dense<T: RealField + Copy>(sparse: &CsrMatrix<T>) -> DMatrix<T> {
        let mut dense = DMatrix::zeros(sparse.nrows(), sparse.ncols());

        for (row_idx, row) in sparse.row_iter().enumerate() {
            for (col_idx, value) in row.col_indices().iter().zip(row.values()) {
                dense[(row_idx, *col_idx)] = *value;
            }
        }

        dense
    }

    #[test]
    fn test_sparse_matrix_builder() -> Result<()> {
        let mut builder = SparseMatrixBuilder::<f64>::new(3, 3);
        
        builder.add_entry(0, 0, 1.0)?;
        builder.add_entry(0, 1, 2.0)?;
        builder.add_entry(1, 1, 3.0)?;
        builder.add_entry(2, 2, 4.0)?;
        
        let matrix = builder.build()?;
        
        assert_eq!(matrix.nrows(), 3);
        assert_eq!(matrix.ncols(), 3);
        assert_eq!(matrix.nnz(), 4);
        
        assert_relative_eq!(matrix.get(0, 0).unwrap_or(0.0), 1.0);
        assert_relative_eq!(matrix.get(0, 1).unwrap_or(0.0), 2.0);
        assert_relative_eq!(matrix.get(1, 1).unwrap_or(0.0), 3.0);
        assert_relative_eq!(matrix.get(2, 2).unwrap_or(0.0), 4.0);
        
        Ok(())
    }

    #[test]
    fn test_duplicate_entries_accumulate() -> Result<()> {
        let mut builder = SparseMatrixBuilder::<f64>::new(2, 2);
        
        builder.add_entry(0, 0, 1.0)?;
        builder.add_entry(0, 0, 2.0)?;  // Duplicate entry should accumulate
        builder.add_entry(1, 1, 3.0)?;
        
        let matrix = builder.build()?;
        
        assert_relative_eq!(matrix.get(0, 0).unwrap_or(0.0), 3.0); // 1.0 + 2.0
        assert_relative_eq!(matrix.get(1, 1).unwrap_or(0.0), 3.0);
        
        Ok(())
    }

    #[test]
    fn test_add_triplets() -> Result<()> {
        let mut builder = SparseMatrixBuilder::<f64>::new(3, 3);
        
        let triplets = vec![
            (0, 0, 1.0),
            (1, 1, 2.0),
            (2, 2, 3.0),
        ];
        
        builder.add_triplets(triplets)?;
        let matrix = builder.build()?;
        
        assert_eq!(matrix.nnz(), 3);
        assert_relative_eq!(matrix.get(0, 0).unwrap_or(0.0), 1.0);
        assert_relative_eq!(matrix.get(1, 1).unwrap_or(0.0), 2.0);
        assert_relative_eq!(matrix.get(2, 2).unwrap_or(0.0), 3.0);
        
        Ok(())
    }

    #[test]
    fn test_matrix_vector_multiply() -> Result<()> {
        let mut builder = SparseMatrixBuilder::<f64>::new(3, 3);
        
        builder.add_triplets(vec![
            (0, 0, 2.0),
            (0, 1, 1.0),
            (1, 1, 3.0),
            (2, 2, 4.0),
        ])?;
        
        let matrix = builder.build()?;
        let x = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        
        let y = &matrix * &x;
        
        assert_relative_eq!(y[0], 4.0);  // 2*1 + 1*2
        assert_relative_eq!(y[1], 6.0);  // 3*2
        assert_relative_eq!(y[2], 12.0); // 4*3
        
        Ok(())
    }

    #[test]
    fn test_transpose() -> Result<()> {
        let mut builder = SparseMatrixBuilder::<f64>::new(2, 3);
        
        builder.add_triplets(vec![
            (0, 0, 1.0),
            (0, 2, 2.0),
            (1, 1, 3.0),
        ])?;
        
        let matrix = builder.build()?;
        let transposed = matrix.transpose();
        
        assert_eq!(transposed.nrows(), 3);
        assert_eq!(transposed.ncols(), 2);
        
        assert_relative_eq!(transposed.get(0, 0).unwrap_or(0.0), 1.0);
        assert_relative_eq!(transposed.get(2, 0).unwrap_or(0.0), 2.0);
        assert_relative_eq!(transposed.get(1, 1).unwrap_or(0.0), 3.0);
        
        Ok(())
    }

    #[test]
    fn test_diagonal_extraction() -> Result<()> {
        let mut builder = SparseMatrixBuilder::<f64>::new(3, 3);
        
        builder.add_triplets(vec![
            (0, 0, 1.0),
            (1, 1, 2.0),
            (2, 2, 3.0),
            (0, 1, 4.0),  // Off-diagonal
        ])?;
        
        let matrix = builder.build()?;
        let diag = matrix.diagonal();
        
        assert_eq!(diag.len(), 3);
        assert_relative_eq!(diag[0], 1.0);
        assert_relative_eq!(diag[1], 2.0);
        assert_relative_eq!(diag[2], 3.0);
        
        Ok(())
    }

    #[test]
    fn test_row_iteration() -> Result<()> {
        let mut builder = SparseMatrixBuilder::<f64>::new(2, 3);
        
        builder.add_triplets(vec![
            (0, 0, 1.0),
            (0, 2, 2.0),
            (1, 1, 3.0),
        ])?;
        
        let matrix = builder.build()?;
        
        let row0: Vec<_> = matrix.row(0).collect();
        assert_eq!(row0.len(), 2);
        assert!(row0.contains(&(0, 1.0)));
        assert!(row0.contains(&(2, 2.0)));
        
        let row1: Vec<_> = matrix.row(1).collect();
        assert_eq!(row1.len(), 1);
        assert!(row1.contains(&(1, 3.0)));
        
        Ok(())
    }

    #[test]
    fn test_sparse_matrix_builder_with_duplicates() -> Result<()> {
        let mut builder = SparseMatrixBuilder::new(2, 2).allow_duplicates(true);

        builder.add_triplets(vec![
            (0, 0, 1.0),
            (0, 0, 2.0),  // Duplicate
            (1, 1, 3.0),
        ])?;

        let matrix = builder.build()?;
        
        // Duplicates should accumulate
        assert_relative_eq!(matrix.get(0, 0).unwrap_or(0.0), 3.0);
        
        Ok(())
    }

    #[test]
    fn test_sparse_matrix_builder_triplets() -> Result<()> {
        let mut builder = SparseMatrixBuilder::new(2, 2);
        let triplets = vec![(0, 0, 1.0), (0, 1, 2.0), (1, 1, 3.0)];

        builder.add_triplets(triplets)?;
        let matrix = builder.build()?;

        assert_eq!(matrix.nnz(), 3);
        let dense = sparse_to_dense(&matrix);
        assert_eq!(dense[(0, 0)], 1.0);
        assert_eq!(dense[(0, 1)], 2.0);
        assert_eq!(dense[(1, 1)], 3.0);
        Ok(())
    }

    #[test]
    fn test_sparse_matrix_builder_bounds_check() -> Result<()> {
        let mut builder = SparseMatrixBuilder::new(2, 2);

        // Should fail for out-of-bounds indices
        assert!(builder.add_entry(2, 0, 1.0).is_err());
        assert!(builder.add_entry(0, 2, 1.0).is_err());
        Ok(())
    }

    #[test]
    fn test_matrix_vector_multiplication() -> Result<()> {
        let mut builder = SparseMatrixBuilder::new(3, 3);

        // Create a simple matrix:
        // [1, 2, 0]
        // [0, 3, 1]
        // [0, 0, 4]
        builder.add_triplets(vec![
            (0, 0, 1.0), (0, 1, 2.0),
            (1, 1, 3.0), (1, 2, 1.0),
            (2, 2, 4.0)
        ])?;

        let matrix = builder.build()?;
        let x = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let mut y = DVector::zeros(3);

        matrix.matvec_inplace(&x, &mut y);

        // Expected: [1*1 + 2*2, 3*2 + 1*3, 4*3] = [5, 9, 12]
        assert_relative_eq!(y[0], 5.0, epsilon = 1e-10);
        assert_relative_eq!(y[1], 9.0, epsilon = 1e-10);
        assert_relative_eq!(y[2], 12.0, epsilon = 1e-10);
        Ok(())
    }

    #[test]
    fn test_matrix_stats() -> Result<()> {
        let mut builder = SparseMatrixBuilder::new(3, 3);
        builder.add_triplets(vec![
            (0, 0, 1.0), (0, 1, -2.0),
            (1, 1, 3.0), (2, 2, 0.5)
        ])?;

        let matrix = builder.build()?;
        let stats = matrix.stats();

        assert_eq!(stats.rows, 3);
        assert_eq!(stats.cols, 3);
        assert_eq!(stats.nnz, 4);
        assert_relative_eq!(stats.density, 4.0 / 9.0, epsilon = 1e-10);
        assert_relative_eq!(stats.min_value, -2.0, epsilon = 1e-10);
        assert_relative_eq!(stats.max_value, 3.0, epsilon = 1e-10);
        Ok(())
    }

    #[test]
    fn test_diagonal_extraction() -> Result<()> {
        let mut builder = SparseMatrixBuilder::new(3, 3);
        builder.add_triplets(vec![
            (0, 0, 1.0), (0, 1, 2.0),
            (1, 1, 3.0), (2, 2, 4.0)
        ])?;

        let matrix = builder.build()?;
        let diag = matrix.diagonal();

        assert_eq!(diag.len(), 3);
        assert_relative_eq!(diag[0], 1.0, epsilon = 1e-10);
        assert_relative_eq!(diag[1], 3.0, epsilon = 1e-10);
        assert_relative_eq!(diag[2], 4.0, epsilon = 1e-10);
        Ok(())
    }

    #[test]
    fn test_symmetry_check() -> Result<()> {
        // Symmetric matrix
        let mut builder = SparseMatrixBuilder::new(3, 3);
        builder.add_triplets(vec![
            (0, 0, 1.0), (0, 1, 2.0),
            (1, 0, 2.0), (1, 1, 3.0),
            (2, 2, 4.0)
        ])?;

        let matrix = builder.build()?;
        assert!(matrix.is_symmetric(1e-10));

        // Non-symmetric matrix
        let mut builder2 = SparseMatrixBuilder::new(2, 2);
        builder2.add_triplets(vec![
            (0, 0, 1.0), (0, 1, 2.0),
            (1, 0, 3.0), (1, 1, 4.0)  // (0,1) != (1,0)
        ])?;

        let matrix2 = builder2.build()?;
        assert!(!matrix2.is_symmetric(1e-10));
        Ok(())
    }

    #[test]
    fn test_tridiagonal_pattern() -> Result<()> {
        let matrix = SparsePatterns::tridiagonal(4, -1.0, 2.0, -1.0)?;
        
        assert_eq!(matrix.nrows(), 4);
        assert_eq!(matrix.ncols(), 4);
        
        // Check diagonal
        assert_relative_eq!(matrix.get(0, 0).unwrap_or(0.0), 2.0);
        assert_relative_eq!(matrix.get(1, 1).unwrap_or(0.0), 2.0);
        
        // Check off-diagonals
        assert_relative_eq!(matrix.get(0, 1).unwrap_or(0.0), -1.0);
        assert_relative_eq!(matrix.get(1, 0).unwrap_or(0.0), -1.0);
        
        Ok(())
    }

    #[test]
    fn test_five_point_stencil() -> Result<()> {
        let matrix = SparsePatterns::five_point_stencil(3, 3, 1.0, 1.0)?;
        
        // Check that it creates a 9x9 matrix (3x3 grid)
        assert_eq!(matrix.nrows(), 9);
        assert_eq!(matrix.ncols(), 9);
        
        // Center point should have 4 connections + diagonal
        let center = 4; // (1,1) in 3x3 grid
        let center_row: Vec<_> = matrix.row(center).collect();
        assert_eq!(center_row.len(), 5); // center + 4 neighbors
        
        Ok(())
    }

    #[test]
    fn test_frobenius_norm() -> Result<()> {
        let mut builder = SparseMatrixBuilder::new(2, 2);
        builder.add_triplets(vec![
            (0, 0, 3.0), (0, 1, 4.0),
            (1, 0, 0.0), (1, 1, 0.0)
        ])?;

        let matrix = builder.build()?;
        let norm = matrix.norm_frobenius();

        // ||A||_F = sqrt(3^2 + 4^2) = 5.0
        assert_relative_eq!(norm, 5.0, epsilon = 1e-10);
        Ok(())
    }
}