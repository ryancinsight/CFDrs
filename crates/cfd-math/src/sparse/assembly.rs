//! Parallel assembly utilities for sparse matrices

use super::builder::{MatrixEntry, SparseMatrixBuilder};
use cfd_core::error::Result;
use nalgebra::RealField;
use nalgebra_sparse::CsrMatrix;
use rayon::prelude::*;
/// Parallel assembly utilities
pub struct ParallelAssembly;
impl ParallelAssembly {
    /// Assemble matrix from element contributions in parallel
    pub fn from_elements<T, F>(
        n: usize,
        elements: &[(Vec<usize>, Vec<Vec<T>>)],
        combine: F,
    ) -> Result<CsrMatrix<T>>
    where
        T: RealField + Copy + Send + Sync,
        F: Fn(&[Vec<T>]) -> Vec<T> + Sync,
    {
        // Parallel computation of element contributions
        let entries: Vec<MatrixEntry<T>> = elements
            .par_iter()
            .flat_map(|(indices, values)| {
                let combined = combine(values);
                let mut local_entries = Vec::new();
                for (i_local, &i_global) in indices.iter().enumerate() {
                    for (j_local, &j_global) in indices.iter().enumerate() {
                        let idx = i_local * indices.len() + j_local;
                        if idx < combined.len() {
                            local_entries.push(MatrixEntry::new(i_global, j_global, combined[idx]));
                        }
                    }
                }
                local_entries
            })
            .collect();
        let mut builder = SparseMatrixBuilder::new(n, n).allow_duplicates(true);
        builder.add_entries(entries)?;
        builder.build_parallel()
    }
    /// Assemble from domain decomposition
    pub fn from_domains<T>(
        domains: Vec<(Vec<usize>, CsrMatrix<T>)>,
        let mut builder = SparseMatrixBuilder::with_capacity(n, n, n * 10).allow_duplicates(true);
        // Add contributions from each domain
        for (indices, local_matrix) in domains {
            for i_local in 0..local_matrix.nrows() {
                let i_global = indices[i_local];
                let row = local_matrix.row(i_local);
                for (idx, &j_local) in row.col_indices().iter().enumerate() {
                    let j_global = indices[j_local];
                    let value = row.values()[idx];
                    builder.add_entry(i_global, j_global, value)?;
            }
        }
    /// Assemble block diagonal matrix
    pub fn block_diagonal<T>(blocks: Vec<CsrMatrix<T>>) -> Result<CsrMatrix<T>>
        T: RealField + Copy,
        if blocks.is_empty() {
            return Ok(CsrMatrix::zeros(0, 0));
        // Calculate total size
        let n: usize = blocks.iter().map(|b| b.nrows()).sum();
        let nnz: usize = blocks.iter().map(|b| b.nnz()).sum();
        let mut builder = SparseMatrixBuilder::with_capacity(n, n, nnz);
        let mut row_offset = 0;
        for block in blocks {
            for i in 0..block.nrows() {
                let row = block.row(i);
                for (idx, &j) in row.col_indices().iter().enumerate() {
                    builder.add_entry(row_offset + i, row_offset + j, row.values()[idx])?;
            row_offset += block.nrows();
        builder.build()
}
