//! Overlapping Schwarz preconditioners for domain decomposition
//!
//! Schwarz methods decompose the domain into subdomains, solve local problems,
//! and combine the solutions. Overlapping subdomains improve convergence.

use super::IncompleteLU;
use crate::linear_solver::Preconditioner;
use cfd_core::error::{Error, Result};
use eunomia::{NumericElement, RealField};
use leto::Array1;
use leto_ops::{CsrMatrix, Scalar as LetoScalar};
use std::collections::{HashMap, HashSet, VecDeque};

#[inline]
fn vector_len<T>(vector: &Array1<T>) -> usize {
    vector.shape()[0]
}

fn validate_vector_len<T>(name: &str, vector: &Array1<T>, expected: usize) -> Result<()> {
    let actual = vector_len(vector);
    if actual != expected {
        return Err(Error::InvalidConfiguration(format!(
            "{name} length mismatch: expected {expected}, got {actual}"
        )));
    }
    Ok(())
}

/// Serial Overlapping Schwarz preconditioner
pub struct SerialSchwarzPreconditioner<T: RealField + Copy + LetoScalar> {
    /// Local subdomain solvers (one per subdomain)
    local_solvers: Vec<IncompleteLU<T>>,
    /// Subdomain boundaries and mappings
    subdomain_map: Vec<Vec<usize>>,
    /// Global matrix dimension.
    system_size: usize,
    /// Overlap size between subdomains
    _overlap: usize,
}

impl<T: RealField + Copy + NumericElement + LetoScalar> SerialSchwarzPreconditioner<T> {
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
            system_size: n,
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

                let row = matrix.row(node);
                for &neighbor in row.col_indices() {
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

                let row = matrix.row(node);
                for &neighbor in row.col_indices() {
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
                let row = matrix.row(node);
                for &neighbor in row.col_indices() {
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

        let mut values = Vec::new();
        let mut col_indices = Vec::new();
        let mut row_ptr = Vec::with_capacity(subdomain_size + 1);
        row_ptr.push(0);

        // Extract relevant entries from global matrix
        for (local_row, &global_row) in indices.iter().enumerate() {
            let row = matrix.row(global_row);
            let mut local_entries = Vec::new();

            for (&global_col, &value) in row.col_indices().iter().zip(row.values()) {
                if let Some(&local_col) = global_to_local.get(&global_col) {
                    local_entries.push((local_col, value));
                }
            }

            local_entries.sort_by_key(|(local_col, _)| *local_col);
            for (local_col, value) in local_entries {
                col_indices.push(local_col);
                values.push(value);
            }
            debug_assert_eq!(row_ptr.len(), local_row + 1);
            row_ptr.push(values.len());
        }

        CsrMatrix::from_parts(values, col_indices, row_ptr, subdomain_size, subdomain_size).map_err(
            |error| {
                Error::InvalidConfiguration(format!(
                    "Schwarz subdomain Leto CSR construction failed: {error}"
                ))
            },
        )
    }

    /// Apply Schwarz preconditioner (additive version)
    ///
    /// This implements the additive Schwarz method where all local solutions
    /// are computed independently and then summed.
    pub fn apply_additive(&self, r: &Array1<T>) -> Result<Array1<T>> {
        validate_vector_len("Schwarz residual", r, self.system_size)?;
        let mut result = Array1::zeros([self.system_size]);

        // Apply each local subdomain solver
        for (subdomain_idx, local_solver) in self.local_solvers.iter().enumerate() {
            let subdomain_indices = &self.subdomain_map[subdomain_idx];

            // Extract local right-hand side
            let mut local_rhs = Array1::zeros([subdomain_indices.len()]);
            for (local_idx, &global_idx) in subdomain_indices.iter().enumerate() {
                local_rhs[local_idx] = r[global_idx];
            }

            // Solve local problem using preconditioner
            let mut local_solution = Array1::zeros([subdomain_indices.len()]);
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
    pub fn apply_multiplicative(&self, r: &Array1<T>) -> Result<Array1<T>> {
        validate_vector_len("Schwarz residual", r, self.system_size)?;
        let mut current_rhs = r.clone();
        let mut result = Array1::zeros([self.system_size]);

        // Solve subdomains sequentially
        for (subdomain_idx, local_solver) in self.local_solvers.iter().enumerate() {
            let subdomain_indices = &self.subdomain_map[subdomain_idx];

            // Extract current local right-hand side
            let mut local_rhs = Array1::zeros([subdomain_indices.len()]);
            for (local_idx, &global_idx) in subdomain_indices.iter().enumerate() {
                local_rhs[local_idx] = current_rhs[global_idx];
            }

            // Solve local problem using preconditioner
            let mut local_solution = Array1::zeros([subdomain_indices.len()]);
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

impl<T: RealField + Copy + NumericElement + LetoScalar> Preconditioner<T>
    for SerialSchwarzPreconditioner<T>
{
    fn apply_to(&self, r: &Array1<T>, z: &mut Array1<T>) -> Result<()> {
        validate_vector_len("Schwarz output", z, self.system_size)?;
        // Use additive Schwarz by default (more parallelizable)
        let result = self.apply_additive(r)?;
        for i in 0..self.system_size {
            z[i] = result[i];
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn identity_matrix(size: usize) -> Result<CsrMatrix<f64>> {
        let mut values = Vec::with_capacity(size);
        let mut col_indices = Vec::with_capacity(size);
        let mut row_ptr = Vec::with_capacity(size + 1);
        row_ptr.push(0);
        for idx in 0..size {
            col_indices.push(idx);
            values.push(1.0);
            row_ptr.push(values.len());
        }
        CsrMatrix::from_parts(values, col_indices, row_ptr, size, size).map_err(|error| {
            Error::InvalidConfiguration(format!("identity Leto CSR construction failed: {error}"))
        })
    }

    fn five_point_stencil_matrix(nx: usize, ny: usize) -> Result<CsrMatrix<f64>> {
        let n = nx * ny;
        let mut values = Vec::with_capacity(5 * n);
        let mut col_indices = Vec::with_capacity(5 * n);
        let mut row_ptr = Vec::with_capacity(n + 1);
        row_ptr.push(0);

        for j in 0..ny {
            for i in 0..nx {
                let idx = j * nx + i;
                if j > 0 {
                    col_indices.push(idx - nx);
                    values.push(1.0);
                }
                if i > 0 {
                    col_indices.push(idx - 1);
                    values.push(1.0);
                }
                col_indices.push(idx);
                values.push(-4.0);
                if i < nx - 1 {
                    col_indices.push(idx + 1);
                    values.push(1.0);
                }
                if j < ny - 1 {
                    col_indices.push(idx + nx);
                    values.push(1.0);
                }
                row_ptr.push(values.len());
            }
        }

        CsrMatrix::from_parts(values, col_indices, row_ptr, n, n).map_err(|error| {
            Error::InvalidConfiguration(format!(
                "five-point stencil Leto CSR construction failed: {error}"
            ))
        })
    }

    fn assert_array_eq(actual: &Array1<f64>, expected: &[f64]) {
        assert_eq!(vector_len(actual), expected.len());
        for (idx, expected_value) in expected.iter().enumerate() {
            assert_eq!(actual[idx], *expected_value, "entry {idx}");
        }
    }

    fn assert_invalid_configuration<T>(result: Result<T>, expected_message: &str) {
        match result {
            Err(Error::InvalidConfiguration(message)) => assert_eq!(message, expected_message),
            Err(error) => panic!("expected invalid configuration, got {error:?}"),
            Ok(_) => panic!("expected invalid configuration error"),
        }
    }

    #[test]
    fn test_2d_partitioning() -> Result<()> {
        let nx = 10;
        let ny = 10;
        let n = nx * ny;
        let matrix = five_point_stencil_matrix(nx, ny)?;

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
        assert!(
            covered.iter().all(|&x| x),
            "All nodes should be covered by subdomains"
        );

        // Verify overlap: sum of subdomain sizes should be > n
        let total_size: usize = preconditioner.subdomain_map.iter().map(|s| s.len()).sum();
        assert!(total_size > n, "Overlap should increase total size");

        Ok(())
    }

    #[test]
    fn schwarz_apply_methods_use_leto_vectors() -> Result<()> {
        let matrix = identity_matrix(4)?;
        let preconditioner = SerialSchwarzPreconditioner::new(&matrix, 2, 0)?;
        let residual =
            Array1::from_shape_vec([4], vec![2.0, -1.0, 0.5, 3.0]).expect("valid residual");

        let additive = preconditioner.apply_additive(&residual)?;
        assert_array_eq(&additive, &[2.0, -1.0, 0.5, 3.0]);

        let multiplicative = preconditioner.apply_multiplicative(&residual)?;
        assert_array_eq(&multiplicative, &[2.0, -1.0, 0.5, 3.0]);

        let mut output = Array1::zeros([4]);
        preconditioner.apply_to(&residual, &mut output)?;
        assert_array_eq(&output, &[2.0, -1.0, 0.5, 3.0]);

        Ok(())
    }

    #[test]
    fn schwarz_rejects_mismatched_leto_vector_lengths() -> Result<()> {
        let matrix = identity_matrix(4)?;
        let preconditioner = SerialSchwarzPreconditioner::new(&matrix, 2, 0)?;
        let residual =
            Array1::from_shape_vec([4], vec![2.0, -1.0, 0.5, 3.0]).expect("valid residual");
        let short_residual = Array1::zeros([3]);

        assert_invalid_configuration(
            preconditioner.apply_additive(&short_residual),
            "Schwarz residual length mismatch: expected 4, got 3",
        );
        assert_invalid_configuration(
            preconditioner.apply_multiplicative(&short_residual),
            "Schwarz residual length mismatch: expected 4, got 3",
        );

        let mut short_output = Array1::zeros([3]);
        assert_invalid_configuration(
            preconditioner.apply_to(&residual, &mut short_output),
            "Schwarz output length mismatch: expected 4, got 3",
        );

        let mut output = Array1::zeros([4]);
        assert_invalid_configuration(
            preconditioner.apply_to(&short_residual, &mut output),
            "Schwarz residual length mismatch: expected 4, got 3",
        );

        Ok(())
    }
}
