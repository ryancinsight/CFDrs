//! Global assembly for spectral element methods.
//!
//! This module provides functionality for assembling global matrices and vectors
//! from element-local contributions in spectral element methods.

use nalgebra::{DMatrix, DVector};

/// Represents a global sparse matrix in compressed sparse row (CSR) format
#[derive(Debug, Clone)]
pub struct SparseMatrixCSR {
    /// Number of rows
    nrows: usize,
    /// Number of columns
    ncols: usize,
    /// Row pointers (size nrows + 1)
    row_ptr: Vec<usize>,
    /// Column indices (size nnz)
    col_ind: Vec<usize>,
    /// Non-zero values (size nnz)
    values: Vec<f64>,
}

impl SparseMatrixCSR {
    /// Create a new empty sparse matrix
    pub fn new(nrows: usize, ncols: usize) -> Self {
        Self {
            nrows,
            ncols,
            row_ptr: vec![0; nrows + 1],
            col_ind: Vec::new(),
            values: Vec::new(),
        }
    }

    /// Get the number of non-zero elements
    pub fn nnz(&self) -> usize {
        self.values.len()
    }

    /// Get the number of rows
    pub fn nrows(&self) -> usize {
        self.nrows
    }

    /// Get the number of columns
    pub fn ncols(&self) -> usize {
        self.ncols
    }

    /// Get a reference to the row pointers
    pub fn row_ptr(&self) -> &[usize] {
        &self.row_ptr
    }

    /// Get a reference to the column indices
    pub fn col_ind(&self) -> &[usize] {
        &self.col_ind
    }

    /// Get a reference to the values
    pub fn values(&self) -> &[f64] {
        &self.values
    }

    /// Convert to a dense matrix (for debugging and small problems)
    pub fn to_dense(&self) -> DMatrix<f64> {
        let mut mat = DMatrix::zeros(self.nrows, self.ncols);

        for i in 0..self.nrows {
            for j in self.row_ptr[i]..self.row_ptr[i + 1] {
                mat[(i, self.col_ind[j])] = self.values[j];
            }
        }

        mat
    }

    /// Matrix-vector product: y = A * x
    pub fn multiply(&self, x: &[f64], y: &mut [f64]) {
        assert_eq!(x.len(), self.ncols);
        assert_eq!(y.len(), self.nrows);

        for i in 0..self.nrows {
            y[i] = 0.0;

            for j in self.row_ptr[i]..self.row_ptr[i + 1] {
                y[i] += self.values[j] * x[self.col_ind[j]];
            }
        }
    }

    /// Add a value to the matrix at (i,j)
    pub fn add_element(&mut self, i: usize, j: usize, value: f64) {
        assert!(i < self.nrows);
        assert!(j < self.ncols);

        // Find the position to insert
        let pos = match self.col_ind[self.row_ptr[i]..self.row_ptr[i + 1]].binary_search(&j) {
            Ok(idx) => self.row_ptr[i] + idx,
            Err(idx) => {
                // Insert new element
                let pos = self.row_ptr[i] + idx;
                self.col_ind.insert(pos, j);
                self.values.insert(pos, 0.0);

                // Update row pointers
                for k in i + 1..=self.nrows {
                    self.row_ptr[k] += 1;
                }

                pos
            }
        };

        self.values[pos] += value;
    }
}

/// Represents a global assembly for spectral element methods
pub struct GlobalAssembly {
    /// Global number of degrees of freedom
    num_dofs: usize,
    /// Element-to-DOF mapping
    element_dofs: Vec<Vec<usize>>,
    /// Global sparse matrix (stored in COO format during assembly)
    coo_rows: Vec<usize>,
    coo_cols: Vec<usize>,
    coo_vals: Vec<f64>,
    /// Global right-hand side vector
    rhs: Vec<f64>,
}

impl GlobalAssembly {
    /// Create a new global assembly
    pub fn new(num_dofs: usize, element_dofs: Vec<Vec<usize>>) -> Self {
        Self {
            num_dofs,
            element_dofs,
            coo_rows: Vec::new(),
            coo_cols: Vec::new(),
            coo_vals: Vec::new(),
            rhs: vec![0.0; num_dofs],
        }
    }

    /// Add a local element matrix to the global system
    pub fn add_element_matrix(
        &mut self,
        elem_idx: usize,
        local_mat: &DMatrix<f64>,
        local_rhs: Option<&DVector<f64>>,
    ) {
        let dofs = &self.element_dofs[elem_idx];
        let n = dofs.len();

        assert_eq!(local_mat.nrows(), n);
        assert_eq!(local_mat.ncols(), n);

        // Add matrix entries
        for i in 0..n {
            let gi = dofs[i];

            if let Some(rhs) = local_rhs {
                self.rhs[gi] += rhs[i];
            }

            for j in 0..n {
                let gj = dofs[j];
                let val = local_mat[(i, j)];

                if val != 0.0 {
                    self.coo_rows.push(gi);
                    self.coo_cols.push(gj);
                    self.coo_vals.push(val);
                }
            }
        }
    }

    /// Assemble the global sparse matrix in CSR format
    pub fn assemble_matrix(&self) -> SparseMatrixCSR {
        if self.coo_rows.is_empty() {
            return SparseMatrixCSR::new(self.num_dofs, self.num_dofs);
        }

        // 1. Collect and sort COO entries
        let mut entries: Vec<(usize, usize, f64)> = self
            .coo_rows
            .iter()
            .zip(self.coo_cols.iter())
            .zip(self.coo_vals.iter())
            .map(|((&r, &c), &v)| (r, c, v))
            .collect();

        // Sort by row, then by column
        entries.sort_by(|a, b| {
            if a.0 == b.0 {
                a.1.cmp(&b.1)
            } else {
                a.0.cmp(&b.0)
            }
        });

        // 2. Consolidate duplicates
        let mut consolidated = Vec::with_capacity(entries.len());
        if !entries.is_empty() {
            let (mut curr_r, mut curr_c, mut curr_v) = entries[0];
            for i in 1..entries.len() {
                let (r, c, v) = entries[i];
                if r == curr_r && c == curr_c {
                    curr_v += v;
                } else {
                    consolidated.push((curr_r, curr_c, curr_v));
                    curr_r = r;
                    curr_c = c;
                    curr_v = v;
                }
            }
            consolidated.push((curr_r, curr_c, curr_v));
        }

        // 3. Build CSR
        let mut mat = SparseMatrixCSR::new(self.num_dofs, self.num_dofs);
        mat.row_ptr = vec![0; self.num_dofs + 1];
        mat.col_ind = Vec::with_capacity(consolidated.len());
        mat.values = Vec::with_capacity(consolidated.len());

        let mut current_row = 0;
        for (r, c, v) in consolidated {
            while current_row < r {
                current_row += 1;
                mat.row_ptr[current_row] = mat.col_ind.len();
            }
            mat.col_ind.push(c);
            mat.values.push(v);
        }

        // Fill remaining row pointers
        while current_row < self.num_dofs {
            current_row += 1;
            mat.row_ptr[current_row] = mat.col_ind.len();
        }

        mat
    }

    /// Get the global right-hand side vector
    pub fn rhs(&self) -> &[f64] {
        &self.rhs
    }

    /// Apply boundary conditions to the global system
    pub fn apply_boundary_conditions<F>(
        &mut self,
        is_boundary: F,
        boundary_value: impl Fn(usize) -> f64,
    ) where
        F: Fn(usize) -> bool,
    {
        // First, identify boundary DOFs
        let mut boundary_dofs = Vec::new();

        for i in 0..self.num_dofs {
            if is_boundary(i) {
                boundary_dofs.push(i);
            }
        }

        // Apply boundary conditions to the matrix and RHS
        for &i in &boundary_dofs {
            // Zero out row i
            for k in 0..self.coo_vals.len() {
                if self.coo_rows[k] == i || self.coo_cols[k] == i {
                    if self.coo_rows[k] == self.coo_cols[k] {
                        // Diagonal entry
                        self.coo_vals[k] = 1.0;
                    } else {
                        // Off-diagonal entry
                        self.coo_vals[k] = 0.0;
                    }
                }
            }

            // Set RHS
            self.rhs[i] = boundary_value(i);
        }
    }
}

/// Represents a spectral element mesh for assembly
pub struct SpectralMesh {
    /// Number of elements
    num_elements: usize,
    /// Number of nodes per element
    nodes_per_element: usize,
    /// Element connectivity (element -> node indices)
    element_connectivity: Vec<Vec<usize>>,
    /// Node coordinates
    node_coords: Vec<f64>,
    /// Boundary node flags
    is_boundary: Vec<bool>,
}

impl SpectralMesh {
    /// Create a new 1D spectral element mesh
    pub fn new_1d(x_min: f64, x_max: f64, num_elements: usize, nodes_per_element: usize) -> Self {
        let num_nodes = (num_elements * (nodes_per_element - 1)) + 1;
        let dx = (x_max - x_min) / num_elements as f64;

        // Generate node coordinates
        let mut node_coords = Vec::with_capacity(num_nodes);
        let mut is_boundary = vec![false; num_nodes];

        // First node is a boundary
        node_coords.push(x_min);
        is_boundary[0] = true;

        // Interior nodes
        for e in 0..num_elements {
            let x0 = x_min + e as f64 * dx;

            // Interior nodes of the element
            for k in 1..nodes_per_element {
                let xi = -1.0 + 2.0 * (k as f64) / ((nodes_per_element - 1) as f64);
                let x = x0 + 0.5 * (xi + 1.0) * dx;
                node_coords.push(x);
            }
        }

        // Last node is a boundary
        is_boundary[num_nodes - 1] = true;

        // Set up element connectivity
        let mut element_connectivity = Vec::with_capacity(num_elements);

        for e in 0..num_elements {
            let mut conn = Vec::with_capacity(nodes_per_element);

            // First node of the element
            conn.push(e * (nodes_per_element - 1));

            // Interior nodes
            for k in 1..nodes_per_element {
                conn.push(e * (nodes_per_element - 1) + k);
            }

            element_connectivity.push(conn);
        }

        Self {
            num_elements,
            nodes_per_element,
            element_connectivity,
            node_coords,
            is_boundary,
        }
    }

    /// Get the number of elements
    pub fn num_elements(&self) -> usize {
        self.num_elements
    }

    /// Get the number of nodes per element
    pub fn nodes_per_element(&self) -> usize {
        self.nodes_per_element
    }

    /// Get the number of nodes
    pub fn num_nodes(&self) -> usize {
        self.node_coords.len()
    }

    /// Get the coordinates of a node
    pub fn node_coord(&self, node_idx: usize) -> f64 {
        self.node_coords[node_idx]
    }

    /// Get the element connectivity
    pub fn element_connectivity(&self, elem_idx: usize) -> &[usize] {
        &self.element_connectivity[elem_idx]
    }

    /// Check if a node is on the boundary
    pub fn is_boundary_node(&self, node_idx: usize) -> bool {
        self.is_boundary[node_idx]
    }

    /// Create a global assembly for this mesh
    pub fn create_assembly(&self) -> GlobalAssembly {
        GlobalAssembly::new(self.num_nodes(), self.element_connectivity.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_sparse_matrix() {
        let n = 5;
        let mut mat = SparseMatrixCSR::new(n, n);

        // Set up a tridiagonal matrix
        for i in 0..n {
            if i > 0 {
                mat.add_element(i, i - 1, -1.0);
            }

            mat.add_element(i, i, 2.0);

            if i < n - 1 {
                mat.add_element(i, i + 1, -1.0);
            }
        }

        // Convert to dense and check
        let dense = mat.to_dense();

        for i in 0..n {
            for j in 0..n {
                let expected = if i == j {
                    2.0
                } else if (i as i32 - j as i32).abs() == 1 {
                    -1.0
                } else {
                    0.0
                };

                assert_relative_eq!(dense[(i, j)], expected, epsilon = 1e-10);
            }
        }
    }

    #[test]
    fn test_global_assembly() {
        // Create a simple 1D mesh with 2 elements and 3 nodes per element
        let mesh = SpectralMesh::new_1d(0.0, 1.0, 2, 3);

        // Create global assembly
        let mut assembly = mesh.create_assembly();

        // Add element contributions (simple mass matrix for each element)
        for e in 0..mesh.num_elements() {
            let n = mesh.nodes_per_element();
            let mut local_mat = DMatrix::zeros(n, n);

            // Simple diagonal matrix for testing
            for i in 0..n {
                local_mat[(i, i)] = 1.0;
            }

            assembly.add_element_matrix(e, &local_mat, None);
        }

        // Assemble global matrix
        let mat = assembly.assemble_matrix();

        // Check that the matrix has the correct structure
        assert_eq!(mat.nrows(), 5); // 2 elements * 3 nodes - 1 shared node
        assert_eq!(mat.ncols(), 5);

        // The matrix should be identity except for the shared node (index 2)
        let dense = mat.to_dense();

        for i in 0..5 {
            for j in 0..5 {
                let expected = if i == j {
                    if i == 2 {
                        2.0 // Shared node
                    } else {
                        1.0
                    }
                } else {
                    0.0
                };

                assert_relative_eq!(dense[(i, j)], expected, epsilon = 1e-10);
            }
        }
    }

    #[test]
    fn test_boundary_conditions() {
        // Create a simple 1D mesh with 2 elements and 3 nodes per element
        let mesh = SpectralMesh::new_1d(0.0, 1.0, 2, 3);

        // Create global assembly
        let mut assembly = mesh.create_assembly();

        // Add element contributions (simple mass matrix for each element)
        for e in 0..mesh.num_elements() {
            let n = mesh.nodes_per_element();
            let mut local_mat = DMatrix::zeros(n, n);

            // Simple diagonal matrix for testing
            for i in 0..n {
                local_mat[(i, i)] = 1.0;
            }

            assembly.add_element_matrix(e, &local_mat, None);
        }

        // Apply boundary conditions (Dirichlet u=0 at x=0 and x=1)
        assembly.apply_boundary_conditions(
            |i| i == 0 || i == 4,               // Boundary nodes
            |i| if i == 0 { 0.0 } else { 1.0 }, // u(0)=0, u(1)=1
        );

        // Assemble global matrix
        let mat = assembly.assemble_matrix();
        let rhs = assembly.rhs();

        // Check boundary conditions
        assert_relative_eq!(rhs[0], 0.0, epsilon = 1e-10); // u(0) = 0
        assert_relative_eq!(rhs[4], 1.0, epsilon = 1e-10); // u(1) = 1

        // Check that the matrix has been modified correctly
        let dense = mat.to_dense();

        // First row should be [1, 0, 0, 0, 0]
        assert_relative_eq!(dense[(0, 0)], 1.0, epsilon = 1e-10);
        for j in 1..5 {
            assert_relative_eq!(dense[(0, j)], 0.0, epsilon = 1e-10);
        }

        // Last row should be [0, 0, 0, 0, 1]
        for j in 0..4 {
            assert_relative_eq!(dense[(4, j)], 0.0, epsilon = 1e-10);
        }
        assert_relative_eq!(dense[(4, 4)], 1.0, epsilon = 1e-10);
    }
}
