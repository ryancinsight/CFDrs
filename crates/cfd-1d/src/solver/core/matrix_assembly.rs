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

use crate::domain::network::Network;
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

impl<T: RealField + Copy + FromPrimitive + Copy + Send + Sync + Copy> MatrixAssembler<T> {
    /// Classify boundary conditions into Dirichlet and Neumann arrays.
    ///
    /// Validates that at least one Dirichlet BC exists and that
    /// no isolated nodes lack a Dirichlet condition.
    pub fn classify_boundary_conditions_into<F: FluidTrait<T>>(
        network: &Network<T, F>,
        dirichlet_values: &mut [Option<T>],
        neumann_sources: &mut [Option<T>],
    ) -> Result<()> {
        let n = network.node_count();
        dirichlet_values.fill(None);
        neumann_sources.fill(None);

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
                crate::domain::network::BoundaryCondition::Dirichlet { value, .. } => {
                    dirichlet_values[idx] = Some(*value);
                    has_dirichlet = true;
                }
                crate::domain::network::BoundaryCondition::Neumann { gradient } => {
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

        Ok(())
    }

    /// Assemble the linear system matrix and right-hand side vector
    ///
    /// The boundary conditions must be pre-classified in the workspace
    /// before calling this function.
    pub fn assemble_into<F: FluidTrait<T>>(
        &self,
        network: &Network<T, F>,
        workspace: &mut crate::solver::core::workspace::SolverWorkspace<T>,
    ) -> Result<CsrMatrix<T>> {
        let n = network.node_count();
        let mut coo = CooMatrix::new(n, n);
        let rhs = &mut workspace.rhs;
        rhs.fill(T::zero());

        let dirichlet_values = &workspace.dirichlet_values;
        let neumann_sources = &workspace.neumann_sources;

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
        for (idx, source) in neumann_sources.iter().enumerate() {
            if let Some(flow_rate) = source {
                rhs[idx] += *flow_rate;
            }
        }

        // Inject identity rows for Dirichlet nodes with exact values
        for (idx, dir_val) in dirichlet_values.iter().enumerate() {
            if let Some(pressure) = dir_val {
                coo.push(idx, idx, T::one());
                rhs[idx] = *pressure;
            }
        }

        Ok(CsrMatrix::from(&coo))
    }

    /// Assemble the linear system matrix and right-hand side vector (allocates workspace)
    pub fn assemble<F: FluidTrait<T>>(
        &self,
        network: &Network<T, F>,
    ) -> Result<(CsrMatrix<T>, DVector<T>)> {
        let n = network.node_count();
        let mut workspace = crate::solver::core::workspace::SolverWorkspace::new(n, 1);
        Self::classify_boundary_conditions_into(
            network,
            &mut workspace.dirichlet_values,
            &mut workspace.neumann_sources,
        )?;
        let matrix = self.assemble_into(network, &mut workspace)?;
        Ok((matrix, workspace.rhs))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::network::{Edge, Network, Node};
    use cfd_core::physics::fluid::ConstantPropertyFluid;
    use cfd_schematics::domain::model::{EdgeKind, NodeKind};
    use petgraph::graph::Graph;

    /// Build a simple A--B--C line graph with conductances g1 (A-B) and g2 (B-C).
    fn abc_network(
        g1: f64,
        g2: f64,
        p_a: f64,
        p_c: f64,
        neumann_b: Option<f64>,
    ) -> Network<f64, ConstantPropertyFluid<f64>> {
        let mut graph = Graph::new();
        let a = graph.add_node(Node::<f64>::new("A".into(), NodeKind::Inlet));
        let b = graph.add_node(Node::<f64>::new("B".into(), NodeKind::Junction));
        let c = graph.add_node(Node::<f64>::new("C".into(), NodeKind::Outlet));

        graph.add_edge(
            a,
            b,
            Edge {
                id: "AB".into(),
                edge_type: EdgeKind::Pipe,
                flow_rate: 0.0,
                resistance: 1.0 / g1,
                quad_coeff: 0.0,
                area: 1.0,
            },
        );
        graph.add_edge(
            b,
            c,
            Edge {
                id: "BC".into(),
                edge_type: EdgeKind::Pipe,
                flow_rate: 0.0,
                resistance: 1.0 / g2,
                quad_coeff: 0.0,
                area: 1.0,
            },
        );

        let fluid = ConstantPropertyFluid::new(
            "test".into(), 1000.0, 1e-3, 4180.0, 0.6, 2.1e9,
        );
        let mut net = Network::new(graph, fluid);
        net.set_pressure(a, p_a);
        net.set_pressure(c, p_c);
        if let Some(q) = neumann_b {
            net.set_neumann_flow(b, q);
        }
        net
    }

    #[test]
    fn dirichlet_nodes_produce_identity_rows() {
        let net = abc_network(1.0, 2.0, 100.0, 0.0, None);
        let assembler = MatrixAssembler::<f64>::new();
        let (mat, rhs) = assembler.assemble(&net).unwrap();

        // Node A (index 0): identity row, RHS = 100
        assert!((mat.get_entry(0, 0).unwrap().into_value() - 1.0).abs() < 1e-12);
        assert!((rhs[0] - 100.0).abs() < 1e-12);

        // Node C (index 2): identity row, RHS = 0
        assert!((mat.get_entry(2, 2).unwrap().into_value() - 1.0).abs() < 1e-12);
        assert!((rhs[2] - 0.0).abs() < 1e-12);

        // Off-diagonals in Dirichlet rows should be zero
        assert!(mat.get_entry(0, 1).is_none() || mat.get_entry(0, 1).unwrap().into_value().abs() < 1e-12);
        assert!(mat.get_entry(2, 1).is_none() || mat.get_entry(2, 1).unwrap().into_value().abs() < 1e-12);
    }

    #[test]
    fn three_node_hand_computed_laplacian() {
        // A--B--C with g1=1, g2=2; A Dirichlet p=100, C Dirichlet p=0
        // Row B (interior): diagonal = g1+g2 = 3, off-diags eliminated by Dirichlet
        // RHS_B = g1*p_A + g2*p_C = 1*100 + 2*0 = 100
        let net = abc_network(1.0, 2.0, 100.0, 0.0, None);
        let assembler = MatrixAssembler::<f64>::new();
        let (mat, rhs) = assembler.assemble(&net).unwrap();

        // Row B (index 1): after Dirichlet column elimination, only diagonal remains
        // A(1,1) = g1 + g2 = 3.0
        let diag_b = mat.get_entry(1, 1).unwrap().into_value();
        assert!((diag_b - 3.0).abs() < 1e-12, "B diagonal should be 3.0, got {diag_b}");

        // Off-diagonal entries (1,0) and (1,2) should be absent (eliminated)
        let off_10 = mat.get_entry(1, 0).map(|e| e.into_value()).unwrap_or(0.0);
        let off_12 = mat.get_entry(1, 2).map(|e| e.into_value()).unwrap_or(0.0);
        assert!(off_10.abs() < 1e-12, "A(1,0) should be 0 after elimination, got {off_10}");
        assert!(off_12.abs() < 1e-12, "A(1,2) should be 0 after elimination, got {off_12}");

        // RHS[1] = g1*p_A + g2*p_C = 100
        assert!((rhs[1] - 100.0).abs() < 1e-12, "RHS[B] should be 100.0, got {}", rhs[1]);
    }

    #[test]
    fn neumann_source_appears_in_rhs() {
        let q_ext = 5.0;
        let net = abc_network(1.0, 2.0, 100.0, 0.0, Some(q_ext));
        let assembler = MatrixAssembler::<f64>::new();
        let (_mat, rhs) = assembler.assemble(&net).unwrap();

        // RHS[B] = g1*p_A + g2*p_C + q_ext = 100 + 0 + 5 = 105
        assert!((rhs[1] - 105.0).abs() < 1e-12, "Neumann source not in RHS: got {}", rhs[1]);
    }

    #[test]
    fn active_block_is_symmetric() {
        // Build a 4-node network: A-B, A-C, B-C, C-D with only A Dirichlet
        let mut graph = Graph::new();
        let a = graph.add_node(Node::<f64>::new("A".into(), NodeKind::Inlet));
        let b = graph.add_node(Node::<f64>::new("B".into(), NodeKind::Junction));
        let c = graph.add_node(Node::<f64>::new("C".into(), NodeKind::Junction));
        let d = graph.add_node(Node::<f64>::new("D".into(), NodeKind::Outlet));

        let mk_edge = |id: &str, g: f64| Edge {
            id: id.into(),
            edge_type: EdgeKind::Pipe,
            flow_rate: 0.0,
            resistance: 1.0 / g,
            quad_coeff: 0.0,
            area: 1.0,
        };
        graph.add_edge(a, b, mk_edge("AB", 1.0));
        graph.add_edge(a, c, mk_edge("AC", 2.0));
        graph.add_edge(b, c, mk_edge("BC", 3.0));
        graph.add_edge(c, d, mk_edge("CD", 4.0));

        let fluid = ConstantPropertyFluid::new(
            "test".into(), 1000.0, 1e-3, 4180.0, 0.6, 2.1e9,
        );
        let mut net = Network::new(graph, fluid);
        net.set_pressure(a, 100.0);
        net.set_pressure(d, 0.0);

        let assembler = MatrixAssembler::<f64>::new();
        let (mat, _rhs) = assembler.assemble(&net).unwrap();

        // Active nodes: B=1, C=2 (non-Dirichlet)
        let active = [1usize, 2];
        for &i in &active {
            for &j in &active {
                let aij = mat.get_entry(i, j).map(|e| e.into_value()).unwrap_or(0.0);
                let aji = mat.get_entry(j, i).map(|e| e.into_value()).unwrap_or(0.0);
                assert!(
                    (aij - aji).abs() < 1e-12,
                    "Symmetry violation: A[{i},{j}]={aij} != A[{j},{i}]={aji}"
                );
            }
        }
    }
}
