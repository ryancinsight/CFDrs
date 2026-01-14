//! Overlapping Schwarz preconditioners for domain decomposition
//!
//! Schwarz methods decompose the domain into subdomains, solve local problems,
//! and combine the solutions. Overlapping subdomains improve convergence.

use super::IncompleteLU;
use crate::linear_solver::Preconditioner;
use cfd_core::error::{Error, Result};
use nalgebra::{DVector, RealField};
use nalgebra_sparse::CsrMatrix;
use num_traits::FromPrimitive;
use std::collections::HashMap;

/// Serial Overlapping Schwarz preconditioner
pub struct SerialSchwarzPreconditioner<T: RealField + Copy> {
    /// Local subdomain solvers (one per subdomain)
    local_solvers: Vec<IncompleteLU<T>>,
    /// Subdomain boundaries and mappings
    subdomain_map: Vec<Vec<usize>>,
    /// Overlap size between subdomains
    _overlap: usize,
}

impl<T: RealField + Copy + FromPrimitive> SerialSchwarzPreconditioner<T> {
    /// Create overlapping Schwarz preconditioner
    ///
    /// # Arguments
    ///
    /// * `matrix` - Global sparse matrix
    /// * `num_subdomains` - Number of subdomains to decompose domain into
    /// * `overlap` - Number of overlapping layers between subdomains
    ///
    /// # Returns
    ///
    /// Schwarz preconditioner with overlapping subdomains
    pub fn new(matrix: &CsrMatrix<T>, num_subdomains: usize, overlap: usize) -> Result<Self> {
        if num_subdomains < 2 {
            return Err(Error::InvalidConfiguration(
                "Need at least 2 subdomains for domain decomposition".to_string(),
            ));
        }

        let n = matrix.nrows();
        if n < num_subdomains {
            return Err(Error::InvalidConfiguration(format!(
                "Matrix size {n} too small for {num_subdomains} subdomains"
            )));
        }

        // Create subdomain partitioning (simple 1D decomposition for now)
        let subdomain_map = Self::create_subdomain_partitioning(n, num_subdomains, overlap);

        // Create local solvers for each subdomain
        let mut local_solvers = Vec::with_capacity(num_subdomains);

        for subdomain_indices in &subdomain_map {
            // Extract local subdomain matrix
            let local_matrix = Self::extract_subdomain_matrix(matrix, subdomain_indices)?;

            // Create ILU preconditioner for local subdomain
            let local_solver = IncompleteLU::new(&local_matrix)?;
            local_solvers.push(local_solver);
        }

        Ok(Self {
            local_solvers,
            subdomain_map,
            _overlap: overlap,
        })
    }

    /// Create subdomain partitioning with overlap
    fn create_subdomain_partitioning(
        n: usize,
        num_subdomains: usize,
        overlap: usize,
    ) -> Vec<Vec<usize>> {
        let mut subdomain_map = Vec::with_capacity(num_subdomains);

        // Simple 1D domain decomposition
        let base_size = n / num_subdomains;
        let remainder = n % num_subdomains;

        let mut start_idx = 0;

        for i in 0..num_subdomains {
            let mut end_idx = start_idx + base_size;
            if i < remainder {
                end_idx += 1;
            }

            // Add overlap to the left (except for first subdomain)
            let actual_start = if i > 0 {
                start_idx.saturating_sub(overlap)
            } else {
                start_idx
            };

            // Add overlap to the right (except for last subdomain)
            let actual_end = if i < num_subdomains - 1 {
                (end_idx + overlap).min(n)
            } else {
                end_idx.min(n)
            };

            let subdomain_indices: Vec<usize> = (actual_start..actual_end).collect();
            subdomain_map.push(subdomain_indices);

            start_idx = end_idx;
        }

        subdomain_map
    }

    /// Extract subdomain matrix from global matrix
    fn extract_subdomain_matrix(matrix: &CsrMatrix<T>, indices: &[usize]) -> Result<CsrMatrix<T>> {
        let subdomain_size = indices.len();

        // Create index mapping from global to local
        let mut global_to_local: HashMap<usize, usize> = HashMap::new();
        for (local_idx, &global_idx) in indices.iter().enumerate() {
            global_to_local.insert(global_idx, local_idx);
        }

        // Build local matrix
        let mut builder = crate::sparse::SparseMatrixBuilder::new(subdomain_size, subdomain_size);

        // Extract relevant entries from global matrix
        for (local_row, &global_row) in indices.iter().enumerate() {
            let row_start = matrix.row_offsets()[global_row];
            let row_end = matrix.row_offsets()[global_row + 1];

            for pos in row_start..row_end {
                let global_col = matrix.col_indices()[pos];
                let value = matrix.values()[pos];

                // Check if this column is in our subdomain
                if let Some(&local_col) = global_to_local.get(&global_col) {
                    builder.add_entry(local_row, local_col, value)?;
                }
            }
        }

        builder.build()
    }

    /// Apply Schwarz preconditioner (additive version)
    ///
    /// This implements the additive Schwarz method where all local solutions
    /// are computed independently and then summed.
    pub fn apply_additive(&self, r: &DVector<T>) -> Result<DVector<T>> {
        let mut result = DVector::zeros(r.len());

        // Apply each local subdomain solver
        for (subdomain_idx, local_solver) in self.local_solvers.iter().enumerate() {
            let subdomain_indices = &self.subdomain_map[subdomain_idx];

            // Extract local right-hand side
            let mut local_rhs = DVector::zeros(subdomain_indices.len());
            for (local_idx, &global_idx) in subdomain_indices.iter().enumerate() {
                local_rhs[local_idx] = r[global_idx];
            }

            // Solve local problem using preconditioner
            let mut local_solution = DVector::zeros(local_rhs.len());
            local_solver.apply_to(&local_rhs, &mut local_solution)?;

            // Add local solution to global result
            for (local_idx, &global_idx) in subdomain_indices.iter().enumerate() {
                result[global_idx] += local_solution[local_idx];
            }
        }

        Ok(result)
    }

    /// Apply Schwarz preconditioner (multiplicative version)
    ///
    /// This implements the multiplicative Schwarz method where subdomains
    /// are solved sequentially, updating the right-hand side.
    pub fn apply_multiplicative(&self, r: &DVector<T>) -> Result<DVector<T>> {
        let mut current_rhs = r.clone();
        let mut result = DVector::zeros(r.len());

        // Solve subdomains sequentially
        for (subdomain_idx, local_solver) in self.local_solvers.iter().enumerate() {
            let subdomain_indices = &self.subdomain_map[subdomain_idx];

            // Extract current local right-hand side
            let mut local_rhs = DVector::zeros(subdomain_indices.len());
            for (local_idx, &global_idx) in subdomain_indices.iter().enumerate() {
                local_rhs[local_idx] = current_rhs[global_idx];
            }

            // Solve local problem using preconditioner
            let mut local_solution = DVector::zeros(local_rhs.len());
            local_solver.apply_to(&local_rhs, &mut local_solution)?;

            // Update global solution and right-hand side for next subdomain
            for (local_idx, &global_idx) in subdomain_indices.iter().enumerate() {
                result[global_idx] += local_solution[local_idx];
                // Update RHS for overlapping regions (simple update)
                current_rhs[global_idx] -= local_solution[local_idx];
            }
        }

        Ok(result)
    }
}

impl<T: RealField + Copy> Preconditioner<T> for SerialSchwarzPreconditioner<T> {
    fn apply_to(&self, r: &DVector<T>, z: &mut DVector<T>) -> Result<()> {
        // Use additive Schwarz by default (more parallelizable)
        let result = self.apply_additive(r)?;
        z.copy_from(&result);
        Ok(())
    }
}
