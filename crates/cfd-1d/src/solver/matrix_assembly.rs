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
//! ## Theorem: Kirchhoff's Current Law (KCL)
//!
//! **Theorem**: At every internal (non-boundary) node $i$, the sum of all incoming
//! and outgoing mass flows must equal any external source $Q_{ext,i}$:
//!
//! $$ \sum_{j \in \mathcal{N}(i)} G_{ij} (P_i - P_j) = Q_{ext,i} $$
//!
//! where $P_i$ is the pressure at node $i$, $G_{ij}$ is the hydraulic conductance
//! of the edge connecting $i$ and $j$, and $\mathcal{N}(i)$ is the set of neighbors of $i$.
//!
//! ## Theorem: Network Laplacian
//!
//! **Theorem**: For a network with exclusively internal nodes, the assembled matrix $\mathbf{A}$
//! is the weighted Graph Laplacian $\mathbf{L} = \mathbf{D} - \mathbf{W}$, where $\mathbf{W}$
//! is the adjacency matrix of conductances $W_{ij} = G_{ij}$, and $\mathbf{D}$ is the diagonal
//! degree matrix $D_{ii} = \sum_{j} G_{ij}$.
//!
//! **Properties**:
//! 1. $\mathbf{A}$ is symmetric positive semi-definite (SPSD).
//! 2. The row sums of $\mathbf{A}$ are exactly zero: $\mathbf{A}\mathbf{1} = \mathbf{0}$.
//! 3. The null space is spanned by the constant vector $\mathbf{1}$, implying pressures
//!    are only unique up to an additive constant unless at least one Dirichlet boundary
//!    condition is specified.
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
use cfd_core::physics::fluid::FluidTrait;
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

/// Per-node Dirichlet and Neumann boundary condition values, indexed by node index.
struct BoundaryClassification<T> {
    dirichlet: Vec<Option<T>>,
    neumann: Vec<Option<T>>,
}

impl<T: RealField + Copy + FromPrimitive + Copy + Send + Sync + Copy> MatrixAssembler<T> {
    /// Classify boundary conditions into Dirichlet and Neumann arrays.
    ///
    /// Returns `(dirichlet_values, neumann_sources)` where each entry corresponds
    /// to a node index. Validates that at least one Dirichlet BC exists and that
    /// no isolated nodes lack a Dirichlet condition.
    fn classify_boundary_conditions<F: FluidTrait<T>>(
        network: &Network<T, F>,
        n: usize,
    ) -> Result<BoundaryClassification<T>> {
        let mut dirichlet_values: Vec<Option<T>> = vec![None; n];
        let mut neumann_sources: Vec<Option<T>> = vec![None; n];
        let mut has_dirichlet = false;

        for (&node_idx, bc) in network.boundary_conditions() {
            if network.graph.node_weight(node_idx).is_none() {
                return Err(Error::InvalidConfiguration(format!(
                    "Boundary condition references missing node {}",
                    node_idx.index()
                )));
            }
            let idx: usize = node_idx.index();
            if idx >= n {
                return Err(Error::InvalidConfiguration(format!(
                    "Boundary condition node index {idx} exceeds node count {n}. Deleted nodes are unsupported."
                )));
            }
            match bc {
                crate::network::BoundaryCondition::Dirichlet { value, .. } => {
                    dirichlet_values[idx] = Some(*value);
                    has_dirichlet = true;
                }
                crate::network::BoundaryCondition::Neumann { gradient } => {
                    neumann_sources[idx] = Some(*gradient);
                }
                _ => {
                    return Err(Error::Solver(format!(
                        "Unsupported boundary condition type encountered at node {idx}: {bc:?}"
                    )));
                }
            }
        }

        if !has_dirichlet {
            return Err(Error::InvalidConfiguration(
                "At least one Dirichlet boundary condition is required".to_string(),
            ));
        }

        for node_idx in network.graph.node_indices() {
            let idx = node_idx.index();
            if idx >= n {
                continue;
            }
            if network
                .graph
                .neighbors_undirected(node_idx)
                .next()
                .is_none()
                && dirichlet_values[idx].is_none()
            {
                return Err(Error::InvalidConfiguration(format!(
                    "Isolated node without Dirichlet condition: {idx}"
                )));
            }
        }

        Ok(BoundaryClassification {
            dirichlet: dirichlet_values,
            neumann: neumann_sources,
        })
    }

    /// Assemble the linear system matrix and right-hand side vector
    ///
    /// This builds the system Ax = b where:
    /// - A is the conductance matrix
    /// - x is the pressure vector
    /// - b is the source/sink vector
    pub fn assemble<F: FluidTrait<T>>(
        &self,
        network: &Network<T, F>,
    ) -> Result<(CsrMatrix<T>, DVector<T>)> {
        let n = network.node_count();
        let mut coo = CooMatrix::new(n, n);
        let mut rhs = DVector::zeros(n);

        // Invariants: units [A]=conductance [m³/(Pa·s)] or dimensionless, [b]=flow rate [m³/s] or pressure [Pa]
        // Interior rows: sum(G_ij * (P_i - P_j)) = Q_ext_i => G_ii*P_i + sum(G_ij*P_j) = Q_ext_i
        // Dirichlet rows: 1 * P_i = P_fixed_i
        let bc = Self::classify_boundary_conditions(network, n)?;
        let dirichlet_values = bc.dirichlet;
        let neumann_sources = bc.neumann;

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

            let i_is_dir = dirichlet_values[i].is_some();
            let j_is_dir = dirichlet_values[j].is_some();

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
                    let p_i = dirichlet_values[i].unwrap();
                    rhs[j] += conductance * p_i;
                }
                (false, true) => {
                    // j fixed: remove column j and add contribution to RHS of i
                    coo.push(i, i, conductance);
                    let p_j = dirichlet_values[j].unwrap();
                    rhs[i] += conductance * p_j;
                }
                (true, true) => {
                    // Both fixed: no equation to assemble between fixed nodes
                }
            }
        }

        // Apply Neumann sources to RHS
        for (idx, source) in neumann_sources.into_iter().enumerate() {
            if let Some(flow_rate) = source {
                rhs[idx] += flow_rate;
            }
        }

        // Inject identity rows for Dirichlet nodes with exact values
        for (idx, dir_val) in dirichlet_values.into_iter().enumerate() {
            if let Some(pressure) = dir_val {
                coo.push(idx, idx, T::one());
                rhs[idx] = pressure;
            }
        }

        let matrix = CsrMatrix::from(&coo);

        Ok((matrix, rhs))
    }
}
