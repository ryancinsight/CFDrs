//! Network state representation for 1D CFD

use crate::network::Network;
use cfd_core::physics::fluid::FluidTrait;
use nalgebra::{DVector, RealField};
/// State representation for a 1D network
#[derive(Debug, Clone)]
pub struct NetworkState<T: RealField + Copy> {
    /// Node pressures
    pub pressures: DVector<T>,
    /// Edge flow rates
    pub flow_rates: DVector<T>,
    /// Time for transient simulations
    pub time: T,
}

impl<T: RealField + Copy> NetworkState<T> {
    /// Create a new network state
    #[must_use]
    pub fn new(num_nodes: usize, num_edges: usize) -> Self {
        Self {
            pressures: DVector::zeros(num_nodes),
            flow_rates: DVector::zeros(num_edges),
            time: T::zero(),
        }
    }

    /// Create state from network
    pub fn from_network<F: FluidTrait<T>>(network: &Network<T, F>) -> Self {
        let num_nodes = network.node_count();
        let num_edges = network.edge_count();

        let mut pressures = DVector::zeros(num_nodes);
        for (node_idx, &pressure) in network.pressures() {
            if node_idx.index() < num_nodes {
                pressures[node_idx.index()] = pressure;
            }
        }

        let mut flow_rates = DVector::zeros(num_edges);
        for (edge_idx, &flow) in network.flow_rates() {
            if edge_idx.index() < num_edges {
                flow_rates[edge_idx.index()] = flow;
            }
        }

        Self {
            pressures,
            flow_rates,
            time: T::zero(),
        }
    }

    /// Get time
    pub fn time(&self) -> T {
        self.time
    }

    /// Set time
    pub fn set_time(&mut self, time: T) {
        self.time = time;
    }
}
