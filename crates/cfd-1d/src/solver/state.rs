//! Network state representation for 1D CFD

use crate::network::Network;
use nalgebra::{RealField, DVector};
use cfd_core::State;

/// State representation for a 1D network
#[derive(Debug, Clone)]
pub struct NetworkState<T: RealField> {
    /// Node pressures
    pub pressures: DVector<T>,
    /// Edge flow rates
    pub flow_rates: DVector<T>,
    /// Time for transient simulations
    pub time: T,
}

impl<T: RealField> NetworkState<T> {
    /// Create a new network state
    pub fn new(num_nodes: usize, num_edges: usize) -> Self {
        Self {
            pressures: DVector::zeros(num_nodes),
            flow_rates: DVector::zeros(num_edges),
            time: T::zero(),
        }
    }

    /// Create state from network
    pub fn from_network(network: &Network<T>) -> Self {
        Self {
            pressures: network.pressures().clone(),
            flow_rates: network.flow_rates().clone(),
            time: T::zero(),
        }
    }
}

impl<T: RealField> State<T> for NetworkState<T> {
    fn time(&self) -> T {
        self.time.clone()
    }

    fn set_time(&mut self, time: T) {
        self.time = time;
    }

    fn clone_state(&self) -> Box<dyn State<T>> {
        Box::new(self.clone())
    }
}