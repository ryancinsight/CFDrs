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
        for (i, &p) in network.pressures().iter().enumerate() {
            if i < num_nodes {
                pressures[i] = p;
            }
        }

        let mut flow_rates = DVector::zeros(num_edges);
        for (i, &q) in network.flow_rates().iter().enumerate() {
            if i < num_edges {
                flow_rates[i] = q;
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
