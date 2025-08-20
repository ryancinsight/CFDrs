//! Matrix assembly for network flow equations

use crate::network::Network;
use nalgebra::{RealField, DVector};
use nalgebra_sparse::{CsrMatrix, coo::CooMatrix};
use num_traits::FromPrimitive;
use cfd_core::Result;
use rayon::prelude::*;
use std::sync::Mutex;

/// Matrix assembler for building the linear system from network equations
pub struct MatrixAssembler<T: RealField> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField> MatrixAssembler<T> {
    /// Create a new matrix assembler
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField + FromPrimitive + Send + Sync> MatrixAssembler<T> {
    /// Assemble the linear system matrix and right-hand side vector
    /// 
    /// This builds the system Ax = b where:
    /// - A is the conductance matrix
    /// - x is the pressure vector
    /// - b is the source/sink vector
    pub fn assemble(&self, network: &Network<T>) -> Result<(CsrMatrix<T>, DVector<T>)> {
        let n = network.node_count();
        let coo_mutex = Mutex::new(CooMatrix::new(n, n));
        let mut rhs = DVector::zeros(n);

        // Parallel assembly of matrix entries
        network.edges_parallel().for_each(|edge| {
            let (i, j) = edge.nodes();
            let conductance = edge.conductance();
            
            let mut coo = coo_mutex.lock().unwrap();
            // Add conductance terms to matrix
            coo.push(i, i, conductance);
            coo.push(j, j, conductance);
            coo.push(i, j, -conductance);
            coo.push(j, i, -conductance);
        });

        // Add boundary conditions
        for (node_idx, bc) in network.boundary_conditions() {
            match bc {
                crate::network::BoundaryCondition::Pressure(p) => {
                    // Dirichlet boundary condition
                    let mut coo = coo_mutex.lock().unwrap();
                    // Set row to identity
                    coo.push(node_idx, node_idx, T::one());
                    rhs[node_idx] = p;
                }
                crate::network::BoundaryCondition::FlowRate(q) => {
                    // Neumann boundary condition
                    rhs[node_idx] += q;
                }
            }
        }

        let coo = coo_mutex.into_inner().unwrap();
        let matrix = CsrMatrix::from(&coo);
        
        Ok((matrix, rhs))
    }
}