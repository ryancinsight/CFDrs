//! Problem definition for 1D network flow analysis
//!
//! ## Theorem: Well-Posedness of the 1D Network Problem
//!
//! **Theorem**: A 1D network flow problem is strictly well-posed (admits a unique
//! pressure field $\mathbf{P} \in \mathbb{R}^{|\mathcal{V}|}$) if and only if:
//!
//! 1. **Graph Connectedness**: The network graph $\mathcal{G}$ is weakly connected.
//!    (Isolated sub-graphs require their own separate boundary conditions).
//! 2. **Dirichlet Uniqueness**: There exists at least one node $i \in \mathcal{V}$
//!    with a prescribed Dirichlet boundary condition (fixed pressure).
//! 3. **Positivity of Conductance**: All edges $e \in \mathcal{E}$ have a strictly
//!    positive, finite hydraulic conductance $G_e > 0$.
//!
//! If condition (2) is violated, the corresponding discrete Laplacian is exactly
//! singular with a nullspace rank of 1 (containing the constant vector $\mathbf{1}$).
//! In this physical state, only pressure *differentials* $\Delta P$ and corresponding
//! flows $\mathbf{Q}$ are uniquely defined; the absolute pressure floats.

use super::geometry::NetworkDomain;
use super::state::NetworkState;
use crate::domain::network::Network;
use cfd_core::abstractions::problem::Problem;
use cfd_core::error::Result;
use cfd_core::physics::boundary::BoundaryConditionSet;
use cfd_core::physics::fluid::{ConstantPropertyFluid, FluidTrait};
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Problem definition for 1D network flow analysis
///
/// This encapsulates the network state and configuration as a Problem that can be
/// solved using the core trait system, enabling polymorphism and plugin architecture.
#[derive(Debug, Clone)]
pub struct NetworkProblem<T: RealField + Copy, F: FluidTrait<T> = ConstantPropertyFluid<T>> {
    /// The network to solve
    pub network: Network<T, F>,
    /// Computational domain information
    pub domain: NetworkDomain<T>,
    /// Fluid properties
    pub fluid: F,
    /// Boundary conditions for the network
    boundary_conditions: BoundaryConditionSet<T>,
}

impl<T: RealField + Copy + FromPrimitive + Copy, F: FluidTrait<T> + Clone> NetworkProblem<T, F> {
    /// Create a new network problem
    pub fn new(network: Network<T, F>) -> Self {
        let node_count = network.node_count();
        let characteristic_length = network.characteristic_length();

        Self {
            domain: NetworkDomain::new(node_count, characteristic_length),
            fluid: network.fluid().clone(),
            boundary_conditions: BoundaryConditionSet::new(),
            network,
        }
    }
}

impl<T: RealField + Copy + FromPrimitive + Copy, F: FluidTrait<T> + Clone> Problem<T>
    for NetworkProblem<T, F>
{
    type Domain = NetworkDomain<T>;
    type State = NetworkState<T>;
    type Fluid = F;

    fn domain(&self) -> &NetworkDomain<T> {
        &self.domain
    }

    fn fluid(&self) -> &F {
        &self.fluid
    }

    fn boundary_conditions(&self) -> &BoundaryConditionSet<T> {
        &self.boundary_conditions
    }

    fn initial_state(&self) -> Result<NetworkState<T>> {
        Ok(NetworkState::from_network(&self.network))
    }
}
