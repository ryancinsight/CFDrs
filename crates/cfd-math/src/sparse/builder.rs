//! Sparse matrix builder with efficient assembly

use cfd_core::error::{Error, Result};
use nalgebra::RealField;
use nalgebra_sparse::{CooMatrix, CsrMatrix};
use rayon::prelude::*;
use std::collections::HashMap;

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
            return Err(Error::InvalidConfiguration(format!(
                "Entry ({}, {}) out of bounds for {}x{} matrix",
                row, col, self.rows, self.cols
            )));
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
    pub fn build(self) -> Result<CsrMatrix<T>> {
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
                entry_map
                    .entry(key)
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
    pub fn build_parallel(self) -> Result<CsrMatrix<T>>
    where
        T: Send + Sync,
    {
        if self.entries.is_empty() {
            return Ok(CsrMatrix::zeros(self.rows, self.cols));
        }

        // Zero-copy parallel aggregation using advanced iterator patterns
        let entry_map: HashMap<(usize, usize), T> = self
            .entries
            .par_iter()
            .map(|entry| ((entry.row, entry.col), entry.value))
            .fold(HashMap::new, |mut acc, (key, value)| {
                acc.entry(key).and_modify(|v| *v += value).or_insert(value);
                acc
            })
            .reduce(HashMap::new, |mut acc, map| {
                for (key, value) in map {
                    acc.entry(key).and_modify(|v| *v += value).or_insert(value);
                }
                acc
            });

        // Build COO matrix from aggregated entries
        let mut coo = CooMatrix::new(self.rows, self.cols);
        for ((row, col), value) in entry_map {
            coo.push(row, col, value);
        }

        Ok(CsrMatrix::from(&coo))
    }

    /// Get the number of entries
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Check if builder is empty
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Clear all entries
    pub fn clear(&mut self) {
        self.entries.clear();
    }
}
