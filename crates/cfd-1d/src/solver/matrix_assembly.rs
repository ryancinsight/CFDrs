//! Matrix assembly for network flow equations

use crate::network::Network;
use cfd_core::{Error, Result};
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

        // Use a simple, sequential loop for better performance and clarity
        // The mutex-protected parallel version was actually slower due to lock contention
        for edge in network.edges_parallel() {
            let (i, j) = edge.nodes;
            let conductance = edge.conductance;

            // Add conductance terms to matrix
            coo.push(i, i, conductance);
            coo.push(j, j, conductance);
            coo.push(i, j, -conductance);
            coo.push(j, i, -conductance);
        }

        // Add boundary conditions
        for (node_idx, bc) in network.boundary_conditions() {
            let idx = node_idx.index();
            match bc {
                crate::network::BoundaryCondition::Dirichlet { value: pressure } => {
                    // Dirichlet boundary condition using the "big number" method
                    // This enforces P_i = pressure by making the diagonal term dominate
                    let large_number = T::from_f64(1.0e20)
                        .expect("Failed to represent large number for Dirichlet BC");

                    // Note: In COO format, we're adding to existing diagonal entries
                    // The large number will dominate the conductance terms
                    coo.push(idx, idx, large_number);
                    rhs[idx] = pressure * large_number;
                }
                crate::network::BoundaryCondition::Neumann {
                    gradient: flow_rate,
                } => {
                    // Neumann boundary condition
                    rhs[idx] = rhs[idx] + flow_rate;
                }
                _ => {
                    return Err(Error::Solver(format!(
                        "Unsupported boundary condition type encountered at node {}: {:?}",
                        idx, bc
                    )));
                }
            }
        }

        let matrix = CsrMatrix::from(&coo);

        Ok((matrix, rhs))
    }
}
