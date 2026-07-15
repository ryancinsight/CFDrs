//! Network state representation for 1D CFD
//!
//! ## Theorem: 1D Discrete State Space
//!
//! **Theorem**: The thermodynamic and hydrodynamic state of a 1D network at any
//! instant $t$ is fully described by the tuple $\Sigma(t) = (\mathbf{P}, \mathbf{Q})$, where:
//!
//! - $\mathbf{P} \in \mathbb{R}^{|\mathcal{V}|}$ is the vector of nodal pressures.
//! - $\mathbf{Q} \in \mathbb{R}^{|\mathcal{E}|}$ is the vector of volumetric flow rates
//!   through the network edges.
//!
//! **Corollary (Incompressibility)**: For an incompressible fluid without compliant
//! vessels, $\mathbf{P}$ alone strictly determines $\mathbf{Q}$ via local resistance
//! laws $Q_{ij} = G_{ij}(\mathbf{P}, \mathbf{Q}) \cdot (P_i - P_j)$.
//! Thus, the primary prognostic variable solved for is the nodal pressure field $\mathbf{P}$.

use crate::domain::network::Network;
use crate::scalar::Cfd1dScalar;
use cfd_core::physics::fluid::FluidTrait;
use leto::Array1;
/// State representation for a 1D network
#[derive(Debug, Clone)]
pub struct NetworkState<T: Cfd1dScalar + Copy> {
    /// Node pressures
    pub pressures: Array1<T>,
    /// Edge flow rates
    pub flow_rates: Array1<T>,
    /// Time for transient simulations
    pub time: T,
}

impl<T: Cfd1dScalar + Copy> NetworkState<T> {
    /// Create a new network state
    #[must_use]
    pub fn new(num_nodes: usize, num_edges: usize) -> Self {
        Self {
            pressures: Array1::from_elem([num_nodes], T::zero()),
            flow_rates: Array1::from_elem([num_edges], T::zero()),
            time: T::zero(),
        }
    }

    /// Create state from network
    pub fn from_network<F: FluidTrait<T>>(network: &Network<T, F>) -> Self {
        let num_nodes = network.node_count();
        let num_edges = network.edge_count();

        let mut pressures = vec![T::zero(); num_nodes];
        for (pressure, value) in pressures
            .iter_mut()
            .zip(network.pressures().iter().copied())
        {
            *pressure = value;
        }

        let mut flow_rates = vec![T::zero(); num_edges];
        for (flow_rate, value) in flow_rates
            .iter_mut()
            .zip(network.flow_rates().iter().copied())
        {
            *flow_rate = value;
        }

        Self {
            pressures: Array1::from_shape_vec([num_nodes], pressures)
                .expect("invariant: pressure vector length matches network node count"),
            flow_rates: Array1::from_shape_vec([num_edges], flow_rates)
                .expect("invariant: flow-rate vector length matches network edge count"),
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

#[cfg(test)]
mod tests {
    use super::NetworkState;

    #[test]
    fn new_allocates_leto_state_vectors_with_zero_values() {
        let state = NetworkState::<f64>::new(3, 2);

        assert_eq!(state.pressures.shape(), [3]);
        assert_eq!(state.flow_rates.shape(), [2]);
        assert_eq!(state.time(), 0.0);
        assert_eq!(
            (0..state.pressures.shape()[0])
                .map(|idx| state.pressures[idx])
                .collect::<Vec<_>>(),
            vec![0.0, 0.0, 0.0]
        );
        assert_eq!(
            (0..state.flow_rates.shape()[0])
                .map(|idx| state.flow_rates[idx])
                .collect::<Vec<_>>(),
            vec![0.0, 0.0]
        );
    }

    #[test]
    fn cloned_state_preserves_leto_vector_values() {
        let mut state = NetworkState::<f64>::new(2, 1);
        state.pressures[0] = 101_325.0;
        state.pressures[1] = 99_000.0;
        state.flow_rates[0] = 0.25;
        state.set_time(1.5);

        let cloned = state.clone();

        assert_eq!(cloned.pressures.shape(), [2]);
        assert_eq!(cloned.flow_rates.shape(), [1]);
        assert_eq!(cloned.pressures[0], 101_325.0);
        assert_eq!(cloned.pressures[1], 99_000.0);
        assert_eq!(cloned.flow_rates[0], 0.25);
        assert_eq!(cloned.time(), 1.5);
    }
}
