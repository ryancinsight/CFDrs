//! Sparse matrix builder with efficient assembly

use cfd_core::error::{Error, Result};
use nalgebra::{DVector, RealField};
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
    /// DOFs with strong Dirichlet enforcement: maps DOF index -> (diag_value, prescribed_value).
    /// During build, rows are replaced by `diag_value * I`, and column contributions
    /// from Dirichlet DOFs are eliminated into the RHS vector.
    dirichlet_dofs: HashMap<usize, (T, T)>,
}

impl<T: RealField + Copy> SparseMatrixBuilder<T> {
    /// Create a new builder
    pub fn new(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            entries: Vec::new(),
            allow_duplicates: false,
            dirichlet_dofs: HashMap::new(),
        }
    }

    /// Read-only access to the accumulated entries
    pub fn entries(&self) -> &[MatrixEntry<T>] {
        &self.entries
    }

    /// Create with estimated capacity
    pub fn with_capacity(rows: usize, cols: usize, capacity: usize) -> Self {
        Self {
            rows,
            cols,
            entries: Vec::with_capacity(capacity),
            allow_duplicates: false,
            dirichlet_dofs: HashMap::new(),
        }
    }

    /// Allow duplicate entries (will be summed)
    pub fn allow_duplicates(mut self, allow: bool) -> Self {
        self.allow_duplicates = allow;
        self
    }

    /// Mark a DOF for strong Dirichlet enforcement with column elimination.
    ///
    /// During `build()`, all entries in this row are discarded and the diagonal
    /// is set to `diag_value`. Column entries from Dirichlet DOFs in non-Dirichlet
    /// rows are eliminated into the RHS: `rhs[j] -= K[j,i] * prescribed_value`.
    ///
    /// The caller must also set `rhs[dof] = diag_value * prescribed_value`.
    pub fn set_dirichlet_row(&mut self, row: usize, diag_value: T, prescribed_value: T) {
        self.dirichlet_dofs.insert(row, (diag_value, prescribed_value));
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

    /// Build the sparse matrix with proper Dirichlet enforcement including
    /// column elimination. For each Dirichlet DOF i with prescribed value g_i:
    ///   - Row i is replaced by [0...0 diag_value 0...0]
    ///   - For each non-Dirichlet row j with entry K[j,i]:
    ///     rhs[j] -= K[j,i] * g_i  (column elimination into RHS)
    ///     K[j,i] is set to zero
    ///
    /// This is the standard FEM Dirichlet BC enforcement that maintains the
    /// correct coupling between constrained and free DOFs.
    pub fn build_with_rhs(self, rhs: &mut DVector<T>) -> Result<CsrMatrix<T>> {
        if self.entries.is_empty() && self.dirichlet_dofs.is_empty() {
            return Ok(CsrMatrix::zeros(self.rows, self.cols));
        }

        let mut coo = CooMatrix::new(self.rows, self.cols);

        // Combine duplicate entries and apply Dirichlet enforcement
        let mut entry_map: HashMap<(usize, usize), T> = HashMap::new();

        for entry in &self.entries {
            let row_is_dirichlet = self.dirichlet_dofs.contains_key(&entry.row);
            let col_is_dirichlet = self.dirichlet_dofs.get(&entry.col);

            if row_is_dirichlet {
                // Skip: row will be replaced by diagonal
                continue;
            }

            if let Some(&(_diag, prescribed_val)) = col_is_dirichlet {
                // Column elimination: move contribution to RHS
                // rhs[row] -= K[row, col] * prescribed_value
                if entry.row < rhs.len() {
                    rhs[entry.row] -= entry.value * prescribed_val;
                }
                // Don't add to matrix (column is zeroed)
                continue;
            }

            // Normal entry: accumulate
            let key = (entry.row, entry.col);
            entry_map
                .entry(key)
                .and_modify(|v| *v += entry.value)
                .or_insert(entry.value);
        }

        for ((row, col), value) in entry_map {
            coo.push(row, col, value);
        }

        // Insert clean diagonal entries for Dirichlet-constrained rows
        for (&row, &(diag_val, _)) in &self.dirichlet_dofs {
            coo.push(row, row, diag_val);
        }

        Ok(CsrMatrix::from(&coo))
    }

    /// Build without column elimination (legacy behavior for non-FEM use cases).
    /// Only zeros Dirichlet rows, does NOT eliminate columns.
    pub fn build(self) -> Result<CsrMatrix<T>> {
        if self.entries.is_empty() && self.dirichlet_dofs.is_empty() {
            return Ok(CsrMatrix::zeros(self.rows, self.cols));
        }

        let mut coo = CooMatrix::new(self.rows, self.cols);

        // Combine duplicate entries manually for better control
        let mut entry_map: HashMap<(usize, usize), T> = HashMap::new();

        for entry in &self.entries {
            if self.dirichlet_dofs.contains_key(&entry.row) {
                continue; // will be replaced by diagonal below
            }
            let key = (entry.row, entry.col);
            entry_map
                .entry(key)
                .and_modify(|v| *v += entry.value)
                .or_insert(entry.value);
        }

        for ((row, col), value) in entry_map {
            coo.push(row, col, value);
        }

        // Insert clean diagonal entries for Dirichlet-constrained rows
        for (&row, &(diag_val, _)) in &self.dirichlet_dofs {
            coo.push(row, row, diag_val);
        }

        Ok(CsrMatrix::from(&coo))
    }

    /// Build matrix in parallel for large systems
    pub fn build_parallel(self) -> Result<CsrMatrix<T>>
    where
        T: Send + Sync,
    {
        if self.entries.is_empty() && self.dirichlet_dofs.is_empty() {
            return Ok(CsrMatrix::zeros(self.rows, self.cols));
        }

        // Zero-copy parallel aggregation using advanced iterator patterns
        // Filter out Dirichlet-constrained rows
        let dirichlet = &self.dirichlet_dofs;
        let entry_map: HashMap<(usize, usize), T> = self
            .entries
            .par_iter()
            .filter(|entry| !dirichlet.contains_key(&entry.row))
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

        // Insert clean diagonal entries for Dirichlet-constrained rows
        for (&row, &(diag_val, _)) in &self.dirichlet_dofs {
            coo.push(row, row, diag_val);
        }

        Ok(CsrMatrix::from(&coo))
    }

    /// Get the number of entries
    #[must_use]
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Check if builder is empty
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Clear all entries
    pub fn clear(&mut self) {
        self.entries.clear();
    }
}
