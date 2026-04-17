//! Network solver for multi-branch junction systems
//!
//! Solves complex branching networks (vascular-like structures) by iterating
//! two-way branch solutions through the network hierarchy.

use super::physics::{TwoWayBranchJunction, TwoWayBranchSolution};
use cfd_core::conversion::SafeFromF64;
use cfd_core::error::Error;
use cfd_core::physics::fluid::traits::Fluid as FluidTrait;
use cfd_core::physics::fluid::traits::NonNewtonianFluid;
use nalgebra::RealField;
use num_traits::{FromPrimitive, ToPrimitive};
use serde::{Deserialize, Serialize};

// ============================================================================
// Branching Network Configuration
// ============================================================================

/// Configuration for branching network solver
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BranchingNetworkConfig<T: RealField + Copy> {
    /// Inlet pressure [Pa]
    pub inlet_pressure: T,
    /// Inlet volumetric flow rate [m³/s]
    pub inlet_flow_rate: T,
    /// Outlet pressure [Pa] (typically 0 for gauge pressure)
    pub outlet_pressure: T,
    /// Maximum solver iterations (for implicit methods)
    pub max_iterations: usize,
    /// Convergence tolerance for flow conservation
    pub convergence_tolerance: T,
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + SafeFromF64> Default
    for BranchingNetworkConfig<T>
{
    fn default() -> Self {
        Self {
            inlet_pressure: T::from_f64_or_one(1000.0),
            inlet_flow_rate: T::from_f64_or_one(1e-6),
            outlet_pressure: T::zero(),
            max_iterations: 100,
            convergence_tolerance: T::from_f64_or_one(1e-6),
        }
    }
}

/// Explicit downstream continuation for cascaded branch-junction solves.
///
/// A `TwoWayBranchJunction` yields two physically distinct daughter states. Any
/// downstream junction must choose one of those states as its inlet; averaging
/// them destroys both mass conservation along the represented path and pressure
/// continuity. This enum makes that continuation explicit.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum DownstreamBranchRoute {
    /// Continue the cascade using daughter 1 outlet pressure and flow.
    Daughter1,
    /// Continue the cascade using daughter 2 outlet pressure and flow.
    Daughter2,
}

// ============================================================================
// Branching Network Solver
// ============================================================================

/// Solver for networks with multiple branch junctions
///
/// # Algorithm
///
/// For a single two-way branch junction:
/// 1. Assume flow split ratio based on daughter geometries
/// 2. Calculate Q_1, Q_2 from mass conservation
/// 3. Calculate pressure drops using apparent viscosity
/// 4. Return flow distribution and pressures
///
/// For networks with multiple branch junctions:
/// 1. Traverse network from inlet
/// 2. Solve each branch junction in series
/// 3. Pass downstream pressures upstream for convergence
/// 4. Iterate until pressure distribution converges
pub struct BranchingNetworkSolver<T: RealField + Copy> {
    config: BranchingNetworkConfig<T>,
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + SafeFromF64> BranchingNetworkSolver<T> {
    /// Create new network solver
    pub fn new(config: BranchingNetworkConfig<T>) -> Self {
        Self { config }
    }

    /// Solve a single two-way branch junction
    pub fn solve_two_way_branch<F: FluidTrait<T> + Copy>(
        &self,
        junction: &TwoWayBranchJunction<T>,
        fluid: F,
        inlet_pressure: T,
        inlet_flow_rate: T,
    ) -> Result<TwoWayBranchSolution<T>, Error> {
        junction.solve(
            fluid,
            inlet_flow_rate,
            inlet_pressure,
            T::from_f64_or_one(293.15),   // 20°C standard temperature [K]
            T::from_f64_or_one(101325.0), // 1 atm [Pa]
        )
    }

    /// Solve branching network by forward pass
    ///
    /// # Algorithm
    ///
    /// For a tree network (single path from inlet to each terminal):
    /// 1. Start at inlet with known P, Q
    /// 2. Pass through each branch junction
    /// 3. At each branch junction, distribute flow based on geometry
    /// 4. Calculate pressures in daughter branches
    /// 5. Use daughter pressures as inlet for next level
    ///
    /// # Time Complexity
    ///
    /// O(N) where N is number of branch junctions (single pass through tree)
    pub fn solve_network<F: FluidTrait<T> + Copy>(
        &self,
        branch_junctions: Vec<TwoWayBranchJunction<T>>,
        fluid: F,
    ) -> Result<Vec<TwoWayBranchSolution<T>>, Error> {
        if branch_junctions.len() > 1 {
            return Err(Error::InvalidConfiguration(
                "multi-junction branch cascades require explicit downstream routing; use solve_network_along_path"
                    .to_string(),
            ));
        }

        self.solve_network_along_path(&branch_junctions, &[], fluid)
    }

    /// Solve a cascaded sequence of branch junctions along an explicit daughter path.
    ///
    /// `downstream_routes[k]` selects which daughter state from junction `k`
    /// becomes the inlet state for junction `k + 1`.
    pub fn solve_network_along_path<F: FluidTrait<T> + Copy>(
        &self,
        branch_junctions: &[TwoWayBranchJunction<T>],
        downstream_routes: &[DownstreamBranchRoute],
        fluid: F,
    ) -> Result<Vec<TwoWayBranchSolution<T>>, Error> {
        let expected_routes = branch_junctions.len().saturating_sub(1);
        if downstream_routes.len() != expected_routes {
            return Err(Error::InvalidConfiguration(format!(
                "expected {expected_routes} downstream routes for {} branch junctions, got {}",
                branch_junctions.len(),
                downstream_routes.len()
            )));
        }

        let mut solutions = Vec::with_capacity(branch_junctions.len());
        let mut current_pressure = self.config.inlet_pressure;
        let mut current_flow = self.config.inlet_flow_rate;

        for (junction_index, branch_junction) in branch_junctions.iter().enumerate() {
            let solution = branch_junction.solve(
                fluid,
                current_flow,
                current_pressure,
                T::from_f64_or_one(293.15),
                T::from_f64_or_one(101325.0),
            )?;

            solutions.push(solution);

            if junction_index + 1 < branch_junctions.len() {
                let (next_pressure, next_flow) =
                    Self::downstream_state(&solution, downstream_routes[junction_index]);
                current_pressure = next_pressure;
                current_flow = next_flow;
            }
        }

        Ok(solutions)
    }

    /// Select the daughter outlet state that seeds the next junction in a routed cascade.
    fn downstream_state(
        solution: &TwoWayBranchSolution<T>,
        route: DownstreamBranchRoute,
    ) -> (T, T) {
        match route {
            DownstreamBranchRoute::Daughter1 => (solution.p_1, solution.q_1),
            DownstreamBranchRoute::Daughter2 => (solution.p_2, solution.q_2),
        }
    }

    /// Calculate total flow resistance of a two-way branch junction
    ///
    /// # Physics
    ///
    /// For branching networks, total resistance is NOT simply sum of resistances.
    /// Branch junctions create parallel paths:
    ///
    /// ```text
    /// 1/R_total = 1/R_1 + 1/R_2 + ... (parallel branches)
    /// R_total = R_parent + 1/(1/R_1 + 1/R_2) (serial + parallel)
    /// ```
    ///
    /// # Formula (for individual two-way branch junction)
    ///
    /// ```text
    /// R = ΔP / Q = (128 μ L) / (π D^4)
    /// ```
    pub fn two_way_branch_resistance<F: FluidTrait<T> + NonNewtonianFluid<T> + Copy>(
        &self,
        junction: &TwoWayBranchJunction<T>,
        fluid: F,
    ) -> (T, T, T) {
        // Use reference flow rate for viscosity evaluation
        let q_ref = self.config.inlet_flow_rate;

        let r_parent = TwoWayBranchJunction::pressure_drop(
            &fluid,
            q_ref,
            &junction.parent,
            T::from_f64_or_one(293.15),
            T::from_f64_or_one(101325.0),
        ) / (q_ref + T::from_f64_or_one(1e-15));
        let r_1 = TwoWayBranchJunction::pressure_drop(
            &fluid,
            q_ref / T::from_f64_or_one(2.0),
            &junction.daughter1,
            T::from_f64_or_one(293.15),
            T::from_f64_or_one(101325.0),
        ) / (q_ref / T::from_f64_or_one(2.0) + T::from_f64_or_one(1e-15));
        let r_2 = TwoWayBranchJunction::pressure_drop(
            &fluid,
            q_ref / T::from_f64_or_one(2.0),
            &junction.daughter2,
            T::from_f64_or_one(293.15),
            T::from_f64_or_one(101325.0),
        ) / (q_ref / T::from_f64_or_one(2.0) + T::from_f64_or_one(1e-15));

        (r_parent, r_1, r_2)
    }
}

// ============================================================================
// Network Analysis Results
// ============================================================================

/// Results from network analysis
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NetworkSolutionSummary<T: RealField + Copy> {
    /// Total number of branch junctions solved
    pub num_branch_junctions: usize,
    /// Total pressure drop across network [Pa]
    pub total_pressure_drop: T,
    /// Average pressure at branch junctions [Pa]
    pub average_junction_pressure: T,
    /// Flow conservation error (should be < 1e-10)
    pub mass_conservation_error: T,
    /// Number of iterations to convergence
    pub num_iterations: usize,
    /// Solver converged
    pub converged: bool,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::channel::Channel;
    use cfd_core::error::Error;

    #[test]
    fn test_branching_network_solver_creation() {
        let config = BranchingNetworkConfig::<f64>::default();
        let solver = BranchingNetworkSolver::new(config);

        assert_eq!(solver.config.max_iterations, 100);
    }

    #[test]
    fn test_two_way_branch_single_solve() {
        use crate::domain::channel::ChannelGeometry;

        let config = BranchingNetworkConfig {
            inlet_pressure: 1000.0,
            inlet_flow_rate: 1e-6,
            outlet_pressure: 0.0,
            max_iterations: 100,
            convergence_tolerance: 1e-6,
        };

        let solver = BranchingNetworkSolver::new(config);

        let parent_geom = ChannelGeometry::<f64>::circular(1.0e-2, 2.0e-3, 1e-6);
        let parent = Channel::new(parent_geom);

        let d1_geom = ChannelGeometry::<f64>::circular(1.0e-2, 1.58e-3, 1e-6);
        let d1 = Channel::new(d1_geom);

        let d2_geom = ChannelGeometry::<f64>::circular(1.0e-2, 1.58e-3, 1e-6);
        let d2 = Channel::new(d2_geom);

        let junction = TwoWayBranchJunction::new(parent, d1, d2, 0.5);
        let blood = cfd_core::physics::fluid::blood::CassonBlood::<f64>::normal_blood();

        let solution = solver
            .solve_two_way_branch(&junction, blood, 1000.0, 1e-6)
            .unwrap();

        assert!(solution.q_parent > 0.0);
        assert!(solution.p_parent >= 0.0);
    }

    #[test]
    fn test_solve_network_requires_explicit_routes_for_multiple_junctions() {
        use crate::domain::channel::ChannelGeometry;

        let solver = BranchingNetworkSolver::new(BranchingNetworkConfig::<f64>::default());

        let parent = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, 2.0e-3, 1e-6));
        let daughter1 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, 1.5e-3, 1e-6));
        let daughter2 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, 1.2e-3, 1e-6));
        let junction = TwoWayBranchJunction::new(parent, daughter1, daughter2, 0.35);
        let blood = cfd_core::physics::fluid::blood::CassonBlood::<f64>::normal_blood();

        let error = solver
            .solve_network(vec![junction.clone(), junction], blood)
            .expect_err("multi-junction path without routes must fail closed");

        match error {
            Error::InvalidConfiguration(message) => {
                assert!(message.contains("explicit downstream routing"));
            }
            other => panic!("unexpected error: {other:?}"),
        }
    }

    #[test]
    fn test_solve_network_along_path_propagates_selected_daughter_state() {
        use crate::domain::channel::ChannelGeometry;

        let solver = BranchingNetworkSolver::new(BranchingNetworkConfig::<f64>::default());

        let first_parent = Channel::new(ChannelGeometry::<f64>::circular(1.2e-2, 2.4e-3, 1e-6));
        let first_daughter1 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, 1.6e-3, 1e-6));
        let first_daughter2 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, 1.1e-3, 1e-6));
        let second_parent = Channel::new(ChannelGeometry::<f64>::circular(0.9e-2, 1.1e-3, 1e-6));
        let second_daughter1 = Channel::new(ChannelGeometry::<f64>::circular(0.8e-2, 0.8e-3, 1e-6));
        let second_daughter2 = Channel::new(ChannelGeometry::<f64>::circular(0.8e-2, 0.6e-3, 1e-6));

        let first_junction =
            TwoWayBranchJunction::new(first_parent, first_daughter1, first_daughter2, 0.3);
        let second_junction =
            TwoWayBranchJunction::new(second_parent, second_daughter1, second_daughter2, 0.4);
        let blood = cfd_core::physics::fluid::blood::CassonBlood::<f64>::normal_blood();

        let solutions = solver
            .solve_network_along_path(
                &[first_junction, second_junction],
                &[DownstreamBranchRoute::Daughter2],
                blood,
            )
            .expect("routed branch cascade");

        assert_eq!(solutions.len(), 2);
        assert_eq!(solutions[1].q_parent, solutions[0].q_2);
        assert_eq!(solutions[1].p_parent, solutions[0].p_2);
    }
}
