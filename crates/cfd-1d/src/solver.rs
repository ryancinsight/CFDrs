//! Network solver for 1D CFD analysis with robust linear system solving.
//!
//! This module provides a comprehensive solver for analyzing fluid flow in microfluidic
//! networks using sparse linear algebra and circuit analogies. The solver implements
//! the core CFD suite trait system for seamless integration with the broader framework.

use crate::network::{Network, BoundaryCondition};
use petgraph::visit::EdgeRef;
use cfd_core::{
    Result, Problem, Solver, Configurable, Validatable, 
    NetworkSolverConfig, Error as CoreError, Domain
};
use nalgebra::{RealField, DVector};
use nalgebra_sparse::{CsrMatrix, coo::CooMatrix};
use num_traits::cast::FromPrimitive;
use std::collections::HashMap;

// Re-export for convenience
pub use cfd_core::NetworkSolverConfig as SolverConfig;

/// 1D Network domain for the Problem trait
#[derive(Debug, Clone)]
pub struct NetworkDomain<T: RealField> {
    /// Number of nodes in the network
    pub node_count: usize,
    /// Network characteristic length scale
    pub characteristic_length: T,
}

impl<T: RealField + FromPrimitive> Domain<T> for NetworkDomain<T> {
    fn dimension(&self) -> usize {
        1 // 1D network
    }

    fn contains(&self, _point: &nalgebra::Point3<T>) -> bool {
        // For network domains, all points are conceptually "inside"
        true
    }

    fn bounding_box(&self) -> (nalgebra::Point3<T>, nalgebra::Point3<T>) {
        let zero = T::zero();
        let one = T::one();
        // Simple unit box for network domain
        (
            nalgebra::Point3::new(zero.clone(), zero.clone(), zero),
            nalgebra::Point3::new(one.clone(), one.clone(), one)
        )
    }

    fn volume(&self) -> T {
        // For 1D networks, "volume" is the characteristic length
        self.characteristic_length.clone()
    }
}

/// Problem definition for 1D network flow analysis
/// 
/// This encapsulates the network state and configuration as a Problem that can be
/// solved using the core trait system, enabling polymorphism and plugin architecture.
#[derive(Debug, Clone)]
pub struct NetworkProblem<T: RealField> {
    /// The network to solve
    pub network: Network<T>,
    /// Computational domain information
    pub domain: NetworkDomain<T>,
}

impl<T: RealField + FromPrimitive> Problem<T> for NetworkProblem<T> {
    type Domain = NetworkDomain<T>;
    type State = NetworkState<T>;

    fn domain(&self) -> &Self::Domain {
        &self.domain
    }

    fn fluid(&self) -> &cfd_core::Fluid<T> {
        self.network.fluid()
    }

    fn boundary_conditions(&self) -> &cfd_core::boundary::BoundaryConditionSet<T> {
        // For 1D networks, we use a simplified approach where boundary conditions
        // are handled directly in the network nodes rather than as a separate set
        // This is a placeholder implementation
        static EMPTY_BC_SET: std::sync::OnceLock<cfd_core::boundary::BoundaryConditionSet<f64>> = std::sync::OnceLock::new();
        
        // This is a workaround - in a full implementation, we'd properly convert
        // network boundary conditions to the core format or use a different approach
        unsafe {
            std::mem::transmute::<
                &cfd_core::boundary::BoundaryConditionSet<f64>,
                &cfd_core::boundary::BoundaryConditionSet<T>
            >(EMPTY_BC_SET.get_or_init(|| cfd_core::boundary::BoundaryConditionSet::new()))
        }
    }

    fn initial_state(&self) -> Result<Self::State> {
        Ok(NetworkState {
            pressures: DVector::zeros(self.network.node_count()),
            flow_rates: DVector::zeros(self.network.edge_count()),
        })
    }

    fn validate(&self) -> Result<()> {
        self.network.validate()?;
        Ok(())
    }
}

impl<T: RealField + FromPrimitive> NetworkProblem<T> {
    /// Create a new network problem from a network
    pub fn new(network: Network<T>) -> Self {
        let node_count = network.node_count();
        let characteristic_length = T::from_f64(1.0).unwrap(); // Default 1m
        
        Self {
            network,
            domain: NetworkDomain {
                node_count,
                characteristic_length,
            },
        }
    }
}

/// Network state containing solution variables
#[derive(Debug, Clone)]
pub struct NetworkState<T: RealField> {
    /// Node pressures
    pub pressures: DVector<T>,
    /// Edge flow rates
    pub flow_rates: DVector<T>,
}

/// Solution results from the robust linear solver
#[derive(Debug, Clone)]
pub struct SolutionResult<T: RealField> {
    /// Convergence status
    pub converged: bool,
    /// Number of linear solver iterations
    pub iterations: usize,
    /// Final residual norm
    pub residual: T,
    /// Solution time in seconds
    pub solve_time: f64,
    /// Condition number estimate (for numerical stability assessment)
    pub condition_number: Option<T>,
}

/// Robust network solver implementing core traits for integration with CFD suite
/// 
/// This solver uses sparse linear algebra to simultaneously solve for all node pressures,
/// providing much better convergence guarantees than simple iterative methods.
pub struct NetworkSolver<T: RealField> {
    config: NetworkSolverConfig<T>,
}

impl<T: RealField + FromPrimitive + num_traits::Float> NetworkSolver<T> {
    /// Create a new solver with default configuration
    pub fn new() -> Self {
        Self {
            config: NetworkSolverConfig::default(),
        }
    }

    /// Create a solver with custom configuration
    pub fn with_config(config: NetworkSolverConfig<T>) -> Self {
        Self { config }
    }

    /// Legacy interface for backward compatibility
    /// 
    /// This method provides the old interface while using the new robust solver internally.
    pub fn solve_steady_state(&self, network: &mut Network<T>) -> Result<SolutionResult<T>> {
        let problem = NetworkProblem::new(network.clone());
        let solution = self.solve_linear_system(&problem)?;
        
        // Update the network with the solution
        self.update_network_from_solution(network, &solution)?;
        
        Ok(SolutionResult {
            converged: true, // Linear solver either succeeds or fails
            iterations: 1,   // Direct solver
            residual: T::zero(), // Exact solution within machine precision
            solve_time: 0.0, // Would be measured in real implementation
            condition_number: None, // Could be computed for numerical stability
        })
    }

    /// Assemble the sparse linear system Ax = b for the network
    /// 
    /// This creates the conductance matrix A and source vector b representing
    /// the entire network simultaneously, enabling robust solution.
    fn assemble_system_matrix(
        &self,
        problem: &NetworkProblem<T>,
    ) -> Result<(CsrMatrix<T>, DVector<T>)> {
        let network = &problem.network;
        let num_nodes = network.graph.node_count();
        let mut coo = CooMatrix::new(num_nodes, num_nodes);
        let mut b = DVector::from_element(num_nodes, T::zero());

        // Create mapping from node indices to matrix rows/columns
        let node_indices: HashMap<_, _> = network.graph.node_indices()
            .enumerate()
            .map(|(i, idx)| (idx, i))
            .collect();

        // Build system matrix row by row
        for node_idx in network.graph.node_indices() {
            let node = &network.graph[node_idx];
            let i = node_indices[&node_idx];

            // Handle boundary conditions
            if let Some(bc) = &node.properties.boundary_condition {
                match bc {
                    BoundaryCondition::Dirichlet { value } => {
                        // Fixed pressure: A_ii = 1, b_i = pressure_value
                        coo.push(i, i, T::one());
                        b[i] = value.clone();
                        continue;
                    },
                    BoundaryCondition::PressureInlet { pressure } |
                    BoundaryCondition::PressureOutlet { pressure } => {
                        // Fixed pressure boundary conditions: A_ii = 1, b_i = pressure_value
                        coo.push(i, i, T::one());
                        b[i] = pressure.clone();
                        continue;
                    },
                    BoundaryCondition::Neumann { gradient } => {
                        // Flow rate boundary condition: add to source term
                        b[i] += gradient.clone();
                    },
                    BoundaryCondition::VolumeFlowInlet { flow_rate } => {
                        // Flow rate boundary condition: add to source term
                        b[i] += flow_rate.clone();
                    },
                    _ => {} // Handle other boundary condition types as needed
                }
            }

            // For internal nodes and Neumann BCs, build conductance matrix
            let mut diagonal_term = T::zero();

            // Process all edges connected to this node
            for edge_ref in network.graph.edges(node_idx) {
                let edge = edge_ref.weight();
                let neighbor_idx = if edge_ref.source() == node_idx {
                    edge_ref.target()
                } else {
                    edge_ref.source()
                };
                let j = node_indices[&neighbor_idx];

                // Conductance = 1/Resistance
                let resistance = edge.effective_resistance();
                if resistance <= T::zero() {
                    return Err(CoreError::InvalidConfiguration(
                        format!("Edge {} has non-positive resistance", edge.id)
                    ));
                }
                let conductance = T::one() / resistance;

                // Off-diagonal term: -G_ij (flow from i to j proportional to pressure difference)
                coo.push(i, j, -conductance.clone());

                // Accumulate diagonal term: sum of all conductances
                diagonal_term += conductance.clone();

                // Handle active components (pumps)
                if let Some(pressure_rise) = edge.pressure_rise() {
                    // Determine pump direction relative to current node
                    let sign = if edge_ref.target() == node_idx { 
                        T::one() 
                    } else { 
                        -T::one() 
                    };
                    // Add pump effect to source term
                    b[i] += sign * conductance * pressure_rise;
                }
            }

            // Set diagonal term
            coo.push(i, i, diagonal_term);
        }

        // Convert to efficient CSR format for solving
        let a = CsrMatrix::from(&coo);
        Ok((a, b))
    }

    /// Solve the linear system using robust sparse linear algebra
    fn solve_linear_system(&self, problem: &NetworkProblem<T>) -> Result<NetworkState<T>> {
        let (a, b) = self.assemble_system_matrix(problem)?;

        // Use nalgebra's dense solver for moderate-sized networks
        // In a full implementation, this would use iterative sparse solvers like GMRES
        // for large systems, but direct solvers work well for networks with < 1000 nodes
        
        // Convert sparse matrix to dense for nalgebra's robust LU solver
        let num_nodes = a.nrows();
        let mut a_dense = nalgebra::DMatrix::zeros(num_nodes, num_nodes);
        
        // Copy sparse matrix to dense format
        for (i, j, &value) in a.triplet_iter() {
            a_dense[(i, j)] = value;
        }
        
        let lu = a_dense.lu();
        
        let pressures = lu.solve(&b).ok_or_else(|| {
            CoreError::ConvergenceFailure(
                "Linear system is singular - check for disconnected components".to_string()
            )
        })?;

        // Compute flow rates from pressure solution
        let flow_rates = self.compute_flow_rates(problem, &pressures)?;

        Ok(NetworkState {
            pressures,
            flow_rates,
        })
    }

    /// Compute edge flow rates from node pressures using Ohm's law analogy
    fn compute_flow_rates(
        &self, 
        problem: &NetworkProblem<T>, 
        pressures: &DVector<T>
    ) -> Result<DVector<T>> {
        let network = &problem.network;
        let num_edges = network.graph.edge_count();
        let mut flow_rates = DVector::zeros(num_edges);

        // Create mapping from node indices to solution vector indices
        let node_indices: HashMap<_, _> = network.graph.node_indices()
            .enumerate()
            .map(|(i, idx)| (idx, i))
            .collect();

        for (edge_idx, edge_ref) in network.graph.edge_references().enumerate() {
            let edge = edge_ref.weight();
            let source_idx = node_indices[&edge_ref.source()];
            let target_idx = node_indices[&edge_ref.target()];

            // Pressure difference across edge
            let pressure_diff = pressures[source_idx].clone() - pressures[target_idx].clone();
            
            // Add pump pressure rise if present
            let effective_pressure_diff = if let Some(pressure_rise) = edge.pressure_rise() {
                pressure_diff + pressure_rise
            } else {
                pressure_diff
            };

            // Flow rate = pressure difference / resistance (Ohm's law analogy)
            let resistance = edge.effective_resistance();
            flow_rates[edge_idx] = effective_pressure_diff / resistance;
        }

        Ok(flow_rates)
    }

    /// Update network with solution from linear solver
    fn update_network_from_solution(
        &self,
        network: &mut Network<T>,
        solution: &NetworkState<T>
    ) -> Result<()> {
        // Update node pressures
        for (i, node_idx) in network.graph.node_indices().enumerate() {
            if let Some(node) = network.graph.node_weight_mut(node_idx) {
                node.pressure = Some(solution.pressures[i].clone());
            }
        }

        // Update edge flow rates
        for (i, edge_idx) in network.graph.edge_indices().enumerate() {
            if let Some(edge) = network.graph.edge_weight_mut(edge_idx) {
                edge.flow_rate = Some(solution.flow_rates[i].clone());
                // Compute pressure drop for reporting
                let resistance = edge.effective_resistance();
                edge.pressure_drop = Some(solution.flow_rates[i].clone() * resistance);
            }
        }

        Ok(())
    }
}

/// Implementation of core Solver trait for integration with CFD suite
impl<T: RealField + FromPrimitive + num_traits::Float> Solver<T> for NetworkSolver<T> {
    type Problem = NetworkProblem<T>;
    type Solution = NetworkState<T>;

    fn solve(&mut self, problem: &Self::Problem) -> Result<Self::Solution> {
        self.solve_linear_system(problem)
    }

    fn name(&self) -> &str {
        "1D Network Sparse Linear Solver"
    }
}

/// Implementation of Configurable trait for solver configuration management
impl<T: RealField> Configurable<T> for NetworkSolver<T> {
    type Config = NetworkSolverConfig<T>;

    fn config(&self) -> &Self::Config {
        &self.config
    }

    fn set_config(&mut self, config: Self::Config) {
        self.config = config;
    }
}

/// Implementation of Validatable trait for problem validation
impl<T: RealField + FromPrimitive + num_traits::Float> Validatable<T> for NetworkSolver<T> {
    type Problem = NetworkProblem<T>;

    fn validate_problem(&self, problem: &Self::Problem) -> Result<()> {
        let network = &problem.network;

        // Validate network structure
        problem.validate()?;

        // Check for at least one pressure boundary condition
        let has_pressure_bc = network.graph.node_weights().any(|node| {
            if let Some(bc) = &node.properties.boundary_condition {
                matches!(bc, 
                    BoundaryCondition::Dirichlet { .. } |
                    BoundaryCondition::PressureInlet { .. } |
                    BoundaryCondition::PressureOutlet { .. }
                )
            } else {
                false
            }
        });

        if !has_pressure_bc {
            return Err(CoreError::InvalidConfiguration(
                "Network must have at least one pressure boundary condition for well-posed problem".to_string()
            ));
        }

        // Check for positive resistances
        for edge in network.graph.edge_weights() {
            if edge.effective_resistance() <= T::zero() {
                return Err(CoreError::InvalidConfiguration(
                    format!("Edge {} has non-positive resistance", edge.id)
                ));
            }
        }

        // Check connectivity (ensure solvable system)
        if !network.is_connected() {
            return Err(CoreError::InvalidConfiguration(
                "Network has disconnected components - system matrix will be singular".to_string()
            ));
        }

        Ok(())
    }

    fn can_handle(&self, _problem: &Self::Problem) -> bool {
        // This solver can handle any valid 1D network problem
        true
    }
}

/// Default implementation
impl<T: RealField + FromPrimitive + num_traits::Float> Default for NetworkSolver<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::network::{NetworkBuilder, ChannelProperties};
    use cfd_core::SolverConfiguration;
    use approx::assert_relative_eq;

    #[test]
    fn test_solver_trait_implementation() {
        let mut solver = NetworkSolver::<f64>::new();
        
        // Test that the solver implements the required traits
        assert_eq!(solver.name(), "1D Network Sparse Linear Solver");
        assert!(solver.config().tolerance() > 0.0);
        
        // Create a simple network problem
        let network = NetworkBuilder::new()
            .add_inlet_pressure("inlet", 0.0, 0.0, 1000.0)
            .add_outlet_pressure("outlet", 1.0, 0.0, 0.0)
            .add_channel("ch1", "inlet", "outlet", ChannelProperties::new(100.0, 1.0, 1e-6))
            .build().unwrap();

        let problem = NetworkProblem::new(network);
        
        // Test validation
        assert!(solver.validate_problem(&problem).is_ok());
        assert!(solver.can_handle(&problem));
        
        // Test solving
        let solution = solver.solve(&problem).unwrap();
        assert_eq!(solution.pressures.len(), 2);
        assert_eq!(solution.flow_rates.len(), 1);
        
        // Check pressure boundary conditions are satisfied
        assert_relative_eq!(solution.pressures[0], 1000.0, epsilon = 1e-10);
        assert_relative_eq!(solution.pressures[1], 0.0, epsilon = 1e-10);
        
        // Check flow rate follows Ohm's law: Q = ΔP/R = 1000/100 = 10
        assert_relative_eq!(solution.flow_rates[0], 10.0, epsilon = 1e-10);
    }

    #[test]
    fn test_linear_system_assembly() {
        let solver = NetworkSolver::<f64>::new();
        
        // Create simple two-node network
        let network = NetworkBuilder::new()
            .add_inlet_pressure("inlet", 0.0, 0.0, 100.0)
            .add_outlet_pressure("outlet", 1.0, 0.0, 0.0)
            .add_channel("ch1", "inlet", "outlet", ChannelProperties::new(50.0, 1.0, 1e-6))
            .build().unwrap();

        let problem = NetworkProblem::new(network);
        let (a, b) = solver.assemble_system_matrix(&problem).unwrap();
        
        // Check matrix dimensions
        assert_eq!(a.nrows(), 2);
        assert_eq!(a.ncols(), 2);
        assert_eq!(b.len(), 2);
        
        // Check boundary conditions in source vector
        assert_relative_eq!(b[0], 100.0, epsilon = 1e-10); // Inlet pressure
        assert_relative_eq!(b[1], 0.0, epsilon = 1e-10);   // Outlet pressure
    }

    #[test]
    fn test_validation_failures() {
        let solver = NetworkSolver::<f64>::new();
        
        // Test network without pressure boundary conditions
        let network = NetworkBuilder::new()
            .add_junction("j1", 0.0, 0.0)
            .add_junction("j2", 1.0, 0.0)
            .add_channel("ch1", "j1", "j2", ChannelProperties::new(100.0, 1.0, 1e-6))
            .build_unchecked(); // Use unchecked to skip network validation

        let problem = NetworkProblem::new(network);
        assert!(solver.validate_problem(&problem).is_err());
    }

    #[test]
    fn test_backward_compatibility() {
        let solver = NetworkSolver::<f64>::new();
        
        let mut network = NetworkBuilder::new()
            .add_inlet_pressure("inlet", 0.0, 0.0, 1000.0)
            .add_outlet_pressure("outlet", 1.0, 0.0, 0.0)
            .add_channel("ch1", "inlet", "outlet", ChannelProperties::new(100.0, 1.0, 1e-6))
            .build().unwrap();

        // Test legacy interface still works
        let result = solver.solve_steady_state(&mut network).unwrap();
        assert!(result.converged);
        
        // Check that network was updated
        let inlet = network.get_node("inlet").unwrap();
        let outlet = network.get_node("outlet").unwrap();
        assert!(inlet.pressure.is_some());
        assert!(outlet.pressure.is_some());
        
        let channel = network.get_edge("ch1").unwrap();
        assert!(channel.flow_rate.is_some());
        assert!(channel.pressure_drop.is_some());
    }
}
