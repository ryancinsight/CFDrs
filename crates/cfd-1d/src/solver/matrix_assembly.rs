//! Matrix assembly for network flow equations
//!
//! # Mathematical Foundation
//!
//! This module assembles the discrete Laplacian operator for the network flow equations:
//!
//! $$ \mathbf{A}\mathbf{x} = \mathbf{b} $$
//!
//! where $\mathbf{A}$ is the conductance matrix, $\mathbf{x}$ is the pressure vector, and $\mathbf{b}$
//! contains source terms and boundary condition contributions.
//!
//! ## Dirichlet Boundary Conditions
//!
//! Exact Dirichlet enforcement is achieved via **Row Replacement with Column Elimination**.
//! For a node $i$ with fixed pressure $P_i$:
//!
//! 1. **Row Replacement**: The equation for node $i$ becomes trivial: $1 \cdot x_i = P_i$.
//!    This is implemented by skipping edge contributions to row $i$ and injecting a diagonal 1 later.
//!
//! 2. **Column Elimination**: For every neighbor $j$ connected to $i$, the term $A_{ji} x_i$ is known ($A_{ji} P_i$).
//!    We move this to the RHS: $b_j \leftarrow b_j - A_{ji} P_i$.
//!    Since $A_{ji} = -G_{ij}$ (negative conductance), this becomes $b_j \leftarrow b_j + G_{ij} P_i$.
//!    The matrix entry $A_{ji}$ is then set to 0 (skipped).
//!
//! This preserves the symmetry of the active submatrix for the remaining unknowns, which is crucial for
//! solvers like Conjugate Gradient (CG).
//!
//! ## Positivity Invariants
//!
//! Conductances $G_{ij}$ must be strictly positive.
//! - $G_{ij} \le 0$ represents a physical impossibility (negative resistance) or a disconnected edge.
//! - Non-finite values (NaN/Inf) indicate numerical instability and are rejected immediately.

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

        // Invariants: units [A]=conductance [m³/(Pa·s)] or dimensionless, [b]=flow rate [m³/s] or pressure [Pa]
        // Interior rows: sum(G_ij * (P_i - P_j)) = Q_ext_i => G_ii*P_i + sum(G_ij*P_j) = Q_ext_i
        // Dirichlet rows: 1 * P_i = P_fixed_i
        let mut dirichlet_values: std::collections::HashMap<usize, T> =
            std::collections::HashMap::new();
        let mut neumann_sources: std::collections::HashMap<usize, T> =
            std::collections::HashMap::new();
        for (&node_idx, bc) in network.boundary_conditions() {
            if network.graph.node_weight(node_idx).is_none() {
                return Err(Error::InvalidConfiguration(format!(
                    "Boundary condition references missing node {}",
                    node_idx.index()
                )));
            }
            let idx: usize = node_idx.index();
            match bc {
                crate::network::BoundaryCondition::Dirichlet { value, .. } => {
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

        if dirichlet_values.is_empty() {
            return Err(Error::InvalidConfiguration(
                "At least one Dirichlet boundary condition is required".to_string(),
            ));
        }

        for node_idx in network.graph.node_indices() {
            if network
                .graph
                .neighbors_undirected(node_idx)
                .next()
                .is_none()
                && !dirichlet_values.contains_key(&node_idx.index())
            {
                return Err(Error::InvalidConfiguration(format!(
                    "Isolated node without Dirichlet condition: {}",
                    node_idx.index()
                )));
            }
        }

        // Use a simple, sequential loop for better performance and clarity
        // Exact Dirichlet enforcement via row-replacement:
        // - Skip assembling row/column entries touching Dirichlet nodes
        // - Accumulate neighbor contributions into RHS for non-Dirichlet nodes
        for edge in network.edges_parallel() {
            let (i, j) = edge.nodes;
            let conductance = edge.conductance;

            if i == j {
                return Err(Error::InvalidConfiguration(
                    "Self-loop edge detected in network assembly".to_string(),
                ));
            }

            if !conductance.is_finite() {
                return Err(Error::InvalidConfiguration(format!(
                    "Non-finite conductance encountered in network assembly: {conductance:?}"
                )));
            }

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
        for (idx, flow_rate) in neumann_sources {
            rhs[idx] += flow_rate;
        }

        // Inject identity rows for Dirichlet nodes with exact values
        for (idx, pressure) in dirichlet_values {
            coo.push(idx, idx, T::one());
            rhs[idx] = pressure;
        }

        let matrix = CsrMatrix::from(&coo);

        Ok((matrix, rhs))
    }
}
