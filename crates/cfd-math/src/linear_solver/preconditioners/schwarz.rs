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
use std::collections::{HashMap, HashSet, VecDeque};

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
        if matrix.ncols() != n {
            return Err(Error::InvalidConfiguration(
                "Schwarz preconditioner requires square matrix".to_string(),
            ));
        }

        if n < num_subdomains {
            return Err(Error::InvalidConfiguration(format!(
                "Matrix size {n} too small for {num_subdomains} subdomains"
            )));
        }

        // Support higher-dimensional and graph-based partitioning strategies.
        let subdomain_map = Self::create_subdomain_partitioning(matrix, num_subdomains, overlap);

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

    /// Compute BFS ordering of the matrix graph
    ///
    /// This uses a pseudo-peripheral node heuristic to improve the ordering quality
    /// (similar to Reverse Cuthill-McKee), reducing the "onion ring" effect in partitions.
    fn compute_bfs_ordering(matrix: &CsrMatrix<T>) -> Vec<usize> {
        let n = matrix.nrows();
        let mut ordering = Vec::with_capacity(n);
        let mut global_visited = vec![false; n];
        let mut probe_visited = vec![false; n];
        let mut queue = VecDeque::new();
        // Keep track of nodes in the current component to reset probe_visited efficiently
        let mut component_nodes = Vec::new();

        // Iterate through all nodes to handle disconnected components
        for i in 0..n {
            if global_visited[i] {
                continue;
            }

            // --- Step 1: Probe BFS to find a pseudo-peripheral node ---
            // Start from node i and find the "furthest" node in this component.
            let mut pseudo_peripheral_node = i;

            queue.push_back(i);
            probe_visited[i] = true;
            component_nodes.push(i);

            while let Some(node) = queue.pop_front() {
                pseudo_peripheral_node = node; // The last visited node is the furthest

                let row_start = matrix.row_offsets()[node];
                let row_end = matrix.row_offsets()[node + 1];

                for pos in row_start..row_end {
                    let neighbor = matrix.col_indices()[pos];
                    // Safety check for non-square matrices (though we check in new())
                    if neighbor < n && !probe_visited[neighbor] {
                        probe_visited[neighbor] = true;
                        queue.push_back(neighbor);
                        component_nodes.push(neighbor);
                    }
                }
            }

            // Reset probe_visited for next component usage (if any, though we use component_nodes)
            // Actually, we don't need to reset probe_visited immediately if we use it only for this component,
            // but we need to ensure we don't revisit nodes in global traversal.
            // Since components are disjoint, global_visited handles that.
            // But we need to reset probe_visited to reuse it?
            // Yes, because component_nodes only tracks THIS component.
            for &node in &component_nodes {
                probe_visited[node] = false;
            }
            component_nodes.clear();

            // --- Step 2: Main BFS from the pseudo-peripheral node ---
            // Generate the actual ordering starting from the better root.
            queue.push_back(pseudo_peripheral_node);
            global_visited[pseudo_peripheral_node] = true;

            while let Some(node) = queue.pop_front() {
                ordering.push(node);

                let row_start = matrix.row_offsets()[node];
                let row_end = matrix.row_offsets()[node + 1];

                for pos in row_start..row_end {
                    let neighbor = matrix.col_indices()[pos];
                    if neighbor < n && !global_visited[neighbor] {
                        global_visited[neighbor] = true;
                        queue.push_back(neighbor);
                    }
                }
            }
        }

        ordering
    }

    /// Expand subdomain by adding neighbors from the graph
    fn expand_subdomain(
        matrix: &CsrMatrix<T>,
        core_indices: &[usize],
        overlap: usize,
    ) -> Vec<usize> {
        if overlap == 0 {
            let mut result = core_indices.to_vec();
            result.sort_unstable();
            return result;
        }

        let mut included: HashSet<usize> = core_indices.iter().copied().collect();
        let mut frontier: VecDeque<usize> = core_indices.iter().copied().collect();

        for _ in 0..overlap {
            let mut next_frontier = VecDeque::new();

            while let Some(node) = frontier.pop_front() {
                let row_start = matrix.row_offsets()[node];
                let row_end = matrix.row_offsets()[node + 1];

                for pos in row_start..row_end {
                    let neighbor = matrix.col_indices()[pos];
                    if included.insert(neighbor) {
                        next_frontier.push_back(neighbor);
                    }
                }
            }

            frontier = next_frontier;
            if frontier.is_empty() {
                break;
            }
        }

        let mut result: Vec<usize> = included.into_iter().collect();
        result.sort_unstable();
        result
    }

    /// Create subdomain partitioning with overlap
    fn create_subdomain_partitioning(
        matrix: &CsrMatrix<T>,
        num_subdomains: usize,
        overlap: usize,
    ) -> Vec<Vec<usize>> {
        let n = matrix.nrows();
        let mut subdomain_map = Vec::with_capacity(num_subdomains);

        // Compute BFS ordering to linearize the graph
        // This effectively maps 2D/3D structures to 1D while preserving locality
        let ordering = Self::compute_bfs_ordering(matrix);

        let base_size = n / num_subdomains;
        let remainder = n % num_subdomains;

        let mut start_idx = 0;

        for i in 0..num_subdomains {
            let mut end_idx = start_idx + base_size;
            if i < remainder {
                end_idx += 1;
            }

            // Core subdomain from BFS ordering
            let core_indices = &ordering[start_idx..end_idx];

            // Expand subdomain with overlap
            let expanded_indices = Self::expand_subdomain(matrix, core_indices, overlap);
            subdomain_map.push(expanded_indices);

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sparse::SparsePatterns;

    #[test]
    fn test_2d_partitioning() -> Result<()> {
        let nx = 10;
        let ny = 10;
        let n = nx * ny;
        let matrix = SparsePatterns::five_point_stencil(nx, ny, 1.0, 1.0)?;

        let num_subdomains = 4;
        let overlap = 1;
        let preconditioner = SerialSchwarzPreconditioner::new(&matrix, num_subdomains, overlap)?;

        // Verify we have the correct number of subdomains
        assert_eq!(preconditioner.subdomain_map.len(), num_subdomains);

        // Verify coverage: union of all subdomains should be all nodes
        let mut covered = vec![false; n];
        for subdomain in &preconditioner.subdomain_map {
            for &idx in subdomain {
                covered[idx] = true;
            }
        }
        assert!(covered.iter().all(|&x| x), "All nodes should be covered by subdomains");

        // Verify overlap: sum of subdomain sizes should be > n
        let total_size: usize = preconditioner.subdomain_map.iter().map(|s| s.len()).sum();
        assert!(total_size > n, "Overlap should increase total size");

        Ok(())
    }
}
