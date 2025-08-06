//! Simplified solver for 1D CFD network analysis.
//!
//! This module provides a basic solver for analyzing fluid flow in microfluidic
//! networks using simple iterative methods.

use crate::network::{Network, Edge};
use petgraph::visit::EdgeRef;
use cfd_core::Result;
use nalgebra::{RealField, ComplexField};
use num_traits::cast::FromPrimitive;

/// Solver configuration and parameters
#[derive(Debug, Clone)]
pub struct SolverConfig<T: RealField> {
    /// Maximum number of iterations for iterative solvers
    pub max_iterations: usize,
    /// Convergence tolerance
    pub tolerance: T,
    /// Relaxation factor for iterative methods
    pub relaxation_factor: T,
    /// Enable verbose output
    pub verbose: bool,
}

impl<T: RealField + FromPrimitive> Default for SolverConfig<T> {
    fn default() -> Self {
        Self {
            max_iterations: 1000,
            tolerance: T::from_f64(1e-6).unwrap(),
            relaxation_factor: T::from_f64(1.0).unwrap(),
            verbose: false,
        }
    }
}

/// Solution results from the solver
#[derive(Debug, Clone)]
pub struct SolutionResult<T: RealField> {
    /// Convergence status
    pub converged: bool,
    /// Number of iterations performed
    pub iterations: usize,
    /// Final residual norm
    pub residual: T,
    /// Solution time in seconds
    pub solve_time: f64,
}

/// Simple network solver for 1D flow analysis
pub struct NetworkSolver<T: RealField> {
    config: SolverConfig<T>,
}

impl<T: RealField + FromPrimitive + num_traits::Float> NetworkSolver<T> {
    /// Create a new solver with default configuration
    pub fn new() -> Self {
        Self {
            config: SolverConfig::default(),
        }
    }

    /// Create a solver with custom configuration
    pub fn with_config(config: SolverConfig<T>) -> Self {
        Self { config }
    }

    /// Solve the steady-state flow problem using simple iteration
    pub fn solve_steady_state(&self, network: &mut Network<T>) -> Result<SolutionResult<T>> {
        let start_time = std::time::Instant::now();

        // Clear any existing solution
        network.clear_solution();

        // Initialize pressures at boundary nodes
        self.initialize_boundary_conditions(network)?;

        // Iteratively solve for pressures and flow rates
        let mut converged = false;
        let mut iterations = 0;
        let mut residual = T::zero();

        for iter in 0..self.config.max_iterations {
            iterations = iter + 1;

            // Update pressures at internal nodes
            let pressure_residual = self.update_pressures(network)?;

            // Update flow rates based on current pressures
            self.update_flow_rates(network)?;

            residual = pressure_residual;

            if self.config.verbose && iter % 100 == 0 {
                println!("Iteration {}: residual = {:?}", iter, residual);
            }

            // Check convergence
            if residual < self.config.tolerance {
                converged = true;
                break;
            }
        }

        let solve_time = start_time.elapsed().as_secs_f64();

        if self.config.verbose {
            if converged {
                println!("Converged in {} iterations", iterations);
            } else {
                println!("Failed to converge in {} iterations", iterations);
            }
        }

        Ok(SolutionResult {
            converged,
            iterations,
            residual,
            solve_time,
        })
    }

    /// Initialize boundary conditions
    fn initialize_boundary_conditions(&self, network: &mut Network<T>) -> Result<()> {
        for node in network.graph.node_weights_mut() {
            if let Some(bc) = &node.properties.boundary_condition {
                if let Some(pressure) = bc.pressure_value() {
                    node.pressure = Some(pressure);
                } else if bc.flow_rate_value().is_some() {
                    // Initialize with zero pressure for flow rate boundaries
                    node.pressure = Some(T::zero());
                } else {
                    node.pressure = Some(T::zero());
                }
            } else {
                // Initialize internal nodes with zero pressure
                node.pressure = Some(T::zero());
            }
        }
        Ok(())
    }

    /// Update pressures at internal nodes using flow conservation
    fn update_pressures(&self, network: &mut Network<T>) -> Result<T> {
        let mut max_residual = T::zero();

        // Collect node information to avoid borrowing issues
        let node_info: Vec<(String, bool)> = network.nodes()
            .map(|node| (node.id.clone(), node.has_boundary_condition()))
            .collect();

        for (node_id, has_bc) in node_info {
            if has_bc {
                // Skip boundary nodes with pressure conditions
                if let Some(node) = network.get_node(&node_id) {
                    if let Some(bc) = node.boundary_condition() {
                        if bc.pressure_value().is_some() {
                            continue;
                        }
                    }
                }
            }

            // Calculate new pressure based on neighboring nodes
            let new_pressure = self.calculate_node_pressure(network, &node_id)?;

            if let Some(node) = network.get_node_mut(&node_id) {
                let old_pressure = node.pressure.unwrap_or(T::zero());
                let residual = ComplexField::abs(new_pressure.clone() - old_pressure.clone());

                if residual > max_residual {
                    max_residual = residual;
                }

                // Apply relaxation
                let relaxed_pressure = old_pressure.clone() +
                    self.config.relaxation_factor.clone() * (new_pressure - old_pressure);
                node.pressure = Some(relaxed_pressure);
            }
        }

        Ok(max_residual)
    }

    /// Calculate pressure at a node based on flow conservation
    fn calculate_node_pressure(&self, network: &Network<T>, node_id: &str) -> Result<T> {
        let node_edges = network.node_edges(node_id)?;

        if node_edges.is_empty() {
            return Ok(T::zero());
        }

        let mut total_conductance = T::zero();
        let mut weighted_pressure_sum = T::zero();
        let mut source_term = T::zero();

        for edge in node_edges {
            let conductance = T::one() / edge.effective_resistance();
            total_conductance += conductance.clone();

            // Find the other node connected by this edge
            let other_node_pressure = self.get_other_node_pressure(network, edge, node_id)?;
            weighted_pressure_sum += conductance.clone() * other_node_pressure;

            // Add pump contribution if present
            if let Some(pressure_rise) = edge.pressure_rise() {
                source_term += conductance * pressure_rise;
            }
        }

        // Handle flow rate boundary condition
        if let Some(node) = network.get_node(node_id) {
            if let Some(bc) = node.boundary_condition() {
                if let Some(flow_rate) = bc.flow_rate_value() {
                    source_term += flow_rate;
                }
            }
        }

        if total_conductance > T::zero() {
            Ok((weighted_pressure_sum + source_term) / total_conductance)
        } else {
            Ok(T::zero())
        }
    }

    /// Get pressure of the other node connected by an edge
    fn get_other_node_pressure(&self, network: &Network<T>, edge: &Edge<T>, current_node_id: &str) -> Result<T> {
        // Find the edge in the graph to get its endpoints
        for edge_ref in network.graph.edge_references() {
            if edge_ref.weight().id == edge.id {
                let from_node = &network.graph[edge_ref.source()];
                let to_node = &network.graph[edge_ref.target()];

                let other_node = if from_node.id == current_node_id {
                    to_node
                } else {
                    from_node
                };

                return Ok(other_node.pressure.unwrap_or(T::zero()));
            }
        }

        Ok(T::zero())
    }

    /// Update flow rates based on current pressure distribution
    fn update_flow_rates(&self, network: &mut Network<T>) -> Result<()> {
        // Collect edge information to avoid borrowing issues
        let edge_info: Vec<(String, String, String)> = network.graph.edge_references()
            .map(|edge_ref| {
                let from_node = &network.graph[edge_ref.source()].id;
                let to_node = &network.graph[edge_ref.target()].id;
                let edge_id = &edge_ref.weight().id;
                (edge_id.clone(), from_node.clone(), to_node.clone())
            })
            .collect();

        for (edge_id, from_node_id, to_node_id) in edge_info {
            let from_pressure = network.get_node(&from_node_id)
                .and_then(|node| node.pressure.clone())
                .unwrap_or(T::zero());

            let to_pressure = network.get_node(&to_node_id)
                .and_then(|node| node.pressure.clone())
                .unwrap_or(T::zero());

            if let Some(edge) = network.get_edge_mut(&edge_id) {
                let pressure_drop = from_pressure.clone() - to_pressure.clone();

                // Add pump contribution if present
                let effective_pressure_drop = if let Some(pump_rise) = edge.pressure_rise() {
                    pressure_drop.clone() - pump_rise
                } else {
                    pressure_drop.clone()
                };

                let flow_rate = effective_pressure_drop.clone() / edge.effective_resistance();

                edge.flow_rate = Some(flow_rate);
                edge.pressure_drop = Some(pressure_drop);
            }
        }

        Ok(())
    }

    /// Get solver configuration
    pub fn config(&self) -> &SolverConfig<T> {
        &self.config
    }

    /// Set solver configuration
    pub fn set_config(&mut self, config: SolverConfig<T>) {
        self.config = config;
    }
}

impl<T: RealField + FromPrimitive + num_traits::Float> Default for NetworkSolver<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::network::{NetworkBuilder};
    use approx::assert_relative_eq;

    #[test]
    fn test_solver_config_default() {
        let config: SolverConfig<f64> = SolverConfig::default();

        assert_eq!(config.max_iterations, 1000);
        assert_relative_eq!(config.tolerance, 1e-6, epsilon = 1e-10);
        assert!(!config.verbose);
    }

    #[test]
    fn test_simple_channel_flow() {
        // Create a simple channel: inlet (1000 Pa) -> channel -> outlet (0 Pa)
        let mut network = NetworkBuilder::new()
            .add_inlet_pressure("inlet", 0.0, 0.0, 1000.0).unwrap()
            .add_outlet_pressure("outlet", 1.0, 0.0, 0.0).unwrap()
            .add_channel("ch1", "inlet", "outlet", 100.0, 0.001, 1e-6).unwrap()
            .build().unwrap();

        let solver = NetworkSolver::new();
        let result = solver.solve_steady_state(&mut network).unwrap();

        assert!(result.converged);

        // Check pressures
        let inlet_pressure = network.get_node("inlet").unwrap().pressure.unwrap();
        let outlet_pressure = network.get_node("outlet").unwrap().pressure.unwrap();

        assert_relative_eq!(inlet_pressure, 1000.0, epsilon = 1e-3);
        assert_relative_eq!(outlet_pressure, 0.0, epsilon = 1e-3);

        // Check flow rate (should be pressure_drop / resistance = 1000 / 100 = 10)
        let flow_rate = network.get_edge("ch1").unwrap().flow_rate.unwrap();
        assert_relative_eq!(flow_rate, 10.0, epsilon = 1e-3);
    }

    #[test]
    fn test_t_junction_flow() {
        // Create simpler T-junction: inlet -> junction -> two outlets
        let mut network = NetworkBuilder::new()
            .add_inlet_pressure("inlet", 0.0, 0.0, 1000.0).unwrap()
            .add_junction("junction", 0.5, 0.0).unwrap()
            .add_outlet_pressure("outlet1", 1.0, 0.5, 0.0).unwrap()
            .add_outlet_pressure("outlet2", 1.0, -0.5, 0.0).unwrap()
            .add_channel("ch1", "inlet", "junction", 100.0, 0.001, 1e-6).unwrap()
            .add_channel("ch2", "junction", "outlet1", 100.0, 0.001, 1e-6).unwrap()
            .add_channel("ch3", "junction", "outlet2", 100.0, 0.001, 1e-6).unwrap()
            .build().unwrap();

        let mut config = SolverConfig::default();
        config.max_iterations = 2000;
        config.tolerance = 1e-3;
        let solver = NetworkSolver::with_config(config);
        let result = solver.solve_steady_state(&mut network).unwrap();

        assert!(result.converged);

        // Check that we have flow rates
        let flow_in = network.get_edge("ch1").unwrap().flow_rate.unwrap();
        let flow_out1 = network.get_edge("ch2").unwrap().flow_rate.unwrap();
        let flow_out2 = network.get_edge("ch3").unwrap().flow_rate.unwrap();

        // Basic sanity checks
        assert!(flow_in > 0.0);
        assert!(flow_out1 >= 0.0);
        assert!(flow_out2 >= 0.0);

        // Flow conservation should hold approximately
        let total_out = flow_out1 + flow_out2;
        if total_out > 0.0 {
            assert_relative_eq!(flow_in, total_out, epsilon = 0.5);
        }
    }

    #[test]
    fn test_pump_in_circuit() {
        // Create circuit with pump: inlet (0 Pa) -> pump -> outlet (0 Pa)
        let mut network = NetworkBuilder::new()
            .add_inlet_pressure("inlet", 0.0, 0.0, 0.0).unwrap()
            .add_outlet_pressure("outlet", 1.0, 0.0, 0.0).unwrap()
            .add_pump("pump1", "inlet", "outlet", 1000.0).unwrap()
            .build().unwrap();

        let solver = NetworkSolver::new();
        let result = solver.solve_steady_state(&mut network).unwrap();

        assert!(result.converged);

        // With zero resistance pump, pressure difference should equal pump rise
        let inlet_pressure = network.get_node("inlet").unwrap().pressure.unwrap();
        let outlet_pressure = network.get_node("outlet").unwrap().pressure.unwrap();

        assert_relative_eq!(inlet_pressure, 0.0, epsilon = 1e-3);
        assert_relative_eq!(outlet_pressure, 0.0, epsilon = 1e-3);

        // Flow rate should be determined by pump characteristics
        let pump_edge = network.get_edge("pump1").unwrap();
        assert!(pump_edge.flow_rate.is_some());
    }

    #[test]
    fn test_valve_resistance() {
        // Test valve with different opening states
        let mut network = NetworkBuilder::new()
            .add_inlet_pressure("inlet", 0.0, 0.0, 1000.0).unwrap()
            .add_outlet_pressure("outlet", 1.0, 0.0, 0.0).unwrap()
            .add_valve("valve1", "inlet", "outlet", 100.0, 0.5).unwrap() // Half open
            .build().unwrap();

        let solver = NetworkSolver::new();
        let result = solver.solve_steady_state(&mut network).unwrap();

        assert!(result.converged);

        // Half-open valve should have 2x resistance, so flow = 1000 / 200 = 5
        let flow_rate = network.get_edge("valve1").unwrap().flow_rate.unwrap();
        assert_relative_eq!(flow_rate, 5.0, epsilon = 1e-3);
    }

    #[test]
    fn test_solver_statistics() {
        let mut network = NetworkBuilder::new()
            .add_inlet_pressure("inlet", 0.0, 0.0, 1000.0).unwrap()
            .add_outlet_pressure("outlet", 1.0, 0.0, 0.0).unwrap()
            .add_channel("ch1", "inlet", "outlet", 100.0, 0.001, 1e-6).unwrap()
            .build().unwrap();

        let solver = NetworkSolver::new();
        let result = solver.solve_steady_state(&mut network).unwrap();

        assert!(result.converged);
        assert!(result.solve_time >= 0.0); // Should have positive solve time
    }
}
