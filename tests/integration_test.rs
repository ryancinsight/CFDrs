//! Integration tests for CFD Suite
//! Following SOLID principles

use cfd_suite::prelude::*;
use cfd_suite::core::Result;
use cfd_1d::solver::NetworkProblem;

#[test]
fn test_1d_network_integration() -> Result<()> {
    // Create a simple network
    let fluid = Fluid::<f64>::water()?;
    let mut network = Network::new(fluid);
    
    // Add nodes
    network.add_node(Node::new("inlet".to_string(), NodeType::Inlet));
    network.add_node(Node::new("outlet".to_string(), NodeType::Outlet));
    
    // Add edge
    let props = ChannelProperties::new(1.0, 0.01, 1e-6);
    network.add_edge("inlet", "outlet", props)?;
    
    // Set boundary conditions
    network.set_boundary_condition(
        "inlet",
        BoundaryCondition::PressureInlet { pressure: 101325.0 }
    )?;
    network.set_boundary_condition(
        "outlet", 
        BoundaryCondition::PressureOutlet { pressure: 101225.0 }
    )?;
    
    // Create and solve
    let solver = NetworkSolver::<f64>::new();
    let problem = NetworkProblem::new(network);
    let solution = solver.solve_network(&problem)?;
    
    // Verify solution
    assert!(solution.nodes().count() > 0);
    Ok(())
}

#[test]
fn test_2d_grid_creation() -> Result<()> {
    let grid = StructuredGrid2D::<f64>::new(10, 10, 0.0, 1.0, 0.0, 1.0)?;
    assert_eq!(grid.nx, 10);
    assert_eq!(grid.ny, 10);
    Ok(())
}

#[test]
fn test_fluid_properties() -> Result<()> {
    let water = Fluid::<f64>::water()?;
    assert!(water.density > 0.0);
    assert!(water.density < 2000.0); // Reasonable for water
    Ok(())
}

#[test]
fn test_solver_configuration() {
    let config = SolverConfig::<f64>::builder()
        .tolerance(1e-6)
        .max_iterations(1000)
        .build();
    assert_eq!(config.max_iterations(), 1000);
}

#[cfg(test)]
mod performance {
    use super::*;
    use std::time::Instant;
    
    #[test]
    fn test_solver_performance() -> Result<()> {
        let fluid = Fluid::<f64>::water()?;
        let mut network = Network::new(fluid);
        
        // Create larger network
        for i in 0..10 {
            network.add_node(Node::new(format!("node_{}", i), NodeType::Junction));
        }
        
        // Add edges
        for i in 0..9 {
            let props = ChannelProperties::new(1.0, 0.01, 1e-6);
            network.add_edge(&format!("node_{}", i), &format!("node_{}", i+1), props)?;
        }
        
        // Set boundary conditions
        network.set_boundary_condition(
            "node_0",
            BoundaryCondition::PressureInlet { pressure: 101325.0 }
        )?;
        network.set_boundary_condition(
            "node_9",
            BoundaryCondition::PressureOutlet { pressure: 101225.0 }
        )?;
        
        let solver = NetworkSolver::<f64>::new();
        let problem = NetworkProblem::new(network);
        
        let start = Instant::now();
        let _solution = solver.solve_network(&problem)?;
        let duration = start.elapsed();
        
        // Should complete in reasonable time
        assert!(duration.as_secs() < 1);
        Ok(())
    }
}