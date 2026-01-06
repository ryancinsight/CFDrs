//! Problem definition for 1D network flow analysis

use super::domain::NetworkDomain;
use super::state::NetworkState;
use crate::network::Network;
use cfd_core::physics::boundary::BoundaryConditionSet;
use cfd_core::error::Result;
use cfd_core::physics::fluid::ConstantPropertyFluid;
use cfd_core::abstractions::problem::Problem;
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Problem definition for 1D network flow analysis
///
/// This encapsulates the network state and configuration as a Problem that can be
/// solved using the core trait system, enabling polymorphism and plugin architecture.
#[derive(Debug, Clone)]
pub struct NetworkProblem<T: RealField + Copy> {
    /// The network to solve
    pub network: Network<T>,
    /// Computational domain information
    pub domain: NetworkDomain<T>,
    /// Fluid properties
    fluid: ConstantPropertyFluid<T>,
    /// Boundary conditions for the network
    boundary_conditions: BoundaryConditionSet<T>,
}

impl<T: RealField + Copy + FromPrimitive + Copy> NetworkProblem<T> {
    /// Create a new network problem
    pub fn new(network: Network<T>) -> Self {
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

impl<T: RealField + Copy + FromPrimitive + Copy> Problem<T> for NetworkProblem<T> {
    type Domain = NetworkDomain<T>;
    type State = NetworkState<T>;

    fn domain(&self) -> &NetworkDomain<T> {
        &self.domain
    }

    fn fluid(&self) -> &ConstantPropertyFluid<T> {
        &self.fluid
    }

    fn boundary_conditions(&self) -> &BoundaryConditionSet<T> {
        &self.boundary_conditions
    }

    fn initial_state(&self) -> Result<NetworkState<T>> {
        Ok(NetworkState::from_network(&self.network))
    }
}
