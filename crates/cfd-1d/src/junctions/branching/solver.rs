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
        junction.solve(fluid, inlet_flow_rate, inlet_pressure)
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
        let mut solutions = Vec::new();
        let mut current_pressure = self.config.inlet_pressure;
        let mut current_flow = self.config.inlet_flow_rate;

        for branch_junction in branch_junctions {
            let solution = branch_junction.solve(fluid, current_flow, current_pressure)?;

            // For next branch junction, use daughter pressures and flows
            // (In a real tree, this would route to specific daughter branches)
            current_pressure = (solution.p_1 + solution.p_2) / T::from_f64_or_one(2.0); // Average
            current_flow = (solution.q_1 + solution.q_2) / T::from_f64_or_one(2.0); // Average

            solutions.push(solution);
        }

        Ok(solutions)
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

        let r_parent = TwoWayBranchJunction::pressure_drop(&fluid, q_ref, &junction.parent)
            / (q_ref + T::from_f64_or_one(1e-15));
        let r_1 = TwoWayBranchJunction::pressure_drop(
            &fluid,
            q_ref / T::from_f64_or_one(2.0),
            &junction.daughter1,
        ) / (q_ref / T::from_f64_or_one(2.0) + T::from_f64_or_one(1e-15));
        let r_2 = TwoWayBranchJunction::pressure_drop(
            &fluid,
            q_ref / T::from_f64_or_one(2.0),
            &junction.daughter2,
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
    use crate::channel::Channel;

    #[test]
    fn test_branching_network_solver_creation() {
        let config = BranchingNetworkConfig::<f64>::default();
        let solver = BranchingNetworkSolver::new(config);

        assert_eq!(solver.config.max_iterations, 100);
    }

    #[test]
    fn test_two_way_branch_single_solve() {
        use crate::channel::ChannelGeometry;
        
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
}
