//! Matrix assembly for network flow equations

use crate::network::Network;
use cfd_core::error::{Error, Result};
use nalgebra::{DVector, RealField};
use nalgebra_sparse::{coo::CooMatrix, CsrMatrix};
use num_traits::FromPrimitive;

/// Matrix assembler for building the linear system from network equations
pub struct MatrixAssembler<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy> Default for MatrixAssembler<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy> MatrixAssembler<T> {
    /// Create a new matrix assembler
    #[must_use]
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField + Copy + FromPrimitive + Copy + Send + Sync + Copy> MatrixAssembler<T> {
    /// Assemble the linear system matrix and right-hand side vector
    ///
    /// This builds the system Ax = b where:
    /// - A is the conductance matrix
    /// - x is the pressure vector
    /// - b is the source/sink vector
    pub fn assemble(&self, network: &Network<T>) -> Result<(CsrMatrix<T>, DVector<T>)> {
        let n = network.node_count();
        let mut coo = CooMatrix::new(n, n);
        let mut rhs = DVector::zeros(n);

        // Invariants: units [A]=[Pa·s/m³], [b]=[Pa]
        // Enforce positivity of conductances and exact Dirichlet constraints
        // Gather boundary conditions for exact enforcement decisions
        let mut dirichlet_values: std::collections::HashMap<usize, T> = std::collections::HashMap::new();
        let mut neumann_sources: std::collections::HashMap<usize, T> = std::collections::HashMap::new();
        for (&node_idx, bc) in network.boundary_conditions() {
            let idx = node_idx.index();
            match bc {
                crate::network::BoundaryCondition::Dirichlet { value } => {
                    dirichlet_values.insert(idx, *value);
                }
                crate::network::BoundaryCondition::Neumann { gradient } => {
                    neumann_sources.insert(idx, *gradient);
                }
                _ => {
                    return Err(Error::Solver(format!(
                        "Unsupported boundary condition type encountered at node {idx}: {bc:?}"
                    )));
                }
            }
        }

        // Use a simple, sequential loop for better performance and clarity
        // Exact Dirichlet enforcement via row-replacement:
        // - Skip assembling row/column entries touching Dirichlet nodes
        // - Accumulate neighbor contributions into RHS for non-Dirichlet nodes
        for edge in network.edges_parallel() {
            let (i, j) = edge.nodes;
            let conductance = edge.conductance;

            // Basic validity: conductance must be positive
            if conductance <= T::zero() {
                return Err(Error::InvalidConfiguration(
                    "Non-positive conductance encountered in network assembly".to_string(),
                ));
            }

            let i_is_dir = dirichlet_values.contains_key(&i);
            let j_is_dir = dirichlet_values.contains_key(&j);

            match (i_is_dir, j_is_dir) {
                (false, false) => {
                    // Standard symmetric Laplacian contributions
                    coo.push(i, i, conductance);
                    coo.push(j, j, conductance);
                    coo.push(i, j, -conductance);
                    coo.push(j, i, -conductance);
                }
                (true, false) => {
                    // i fixed: remove column i and add contribution to RHS of j
                    coo.push(j, j, conductance);
                    let p_i = dirichlet_values[&i];
                    rhs[j] += conductance * p_i;
                }
                (false, true) => {
                    // j fixed: remove column j and add contribution to RHS of i
                    coo.push(i, i, conductance);
                    let p_j = dirichlet_values[&j];
                    rhs[i] += conductance * p_j;
                }
                (true, true) => {
                    // Both fixed: no equation to assemble between fixed nodes
                }
            }
        }

        // Apply Neumann sources to RHS
        for (idx, flow_rate) in neumann_sources.into_iter() {
            rhs[idx] += flow_rate;
        }

        // Inject identity rows for Dirichlet nodes with exact values
        for (idx, pressure) in dirichlet_values.into_iter() {
            coo.push(idx, idx, T::one());
            rhs[idx] = pressure;
        }

        let matrix = CsrMatrix::from(&coo);

        Ok((matrix, rhs))
    }
}
