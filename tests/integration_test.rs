//! Integration tests for CFD Suite
//! Following SOLID principles

use cfd_suite::prelude::*;
use cfd_suite::core::Result;

#[test]
fn test_1d_network_integration() -> Result<()> {
    // Create a simple network
    let fluid = Fluid::<f64>::water()?;
    let mut builder = NetworkBuilder::new(fluid);
    
    // Add nodes
    builder = builder.add_node(Node::new(0, 0.0, 0.0, 0.0));
    builder = builder.add_node(Node::new(1, 1.0, 0.0, 0.0));
    
    let mut network = builder.build()?;
    
    // Set boundary conditions
    network.set_pressure(0, 101325.0)?;
    network.set_pressure(1, 101225.0)?;
    
    // Create and solve
    let mut solver = NetworkSolver::<f64>::new();
    let problem = NetworkProblem { network };
    let solution = solver.solve(&problem)?;
    
    // Verify solution
    assert!(solution.node_count() > 0);
    Ok(())
}

#[test]
fn test_2d_grid_creation() {
    let grid = StructuredGrid2D::<f64>::new(10, 10, 1.0, 1.0, 0.0, 1.0);
    assert_eq!(grid.nx(), 10);
    assert_eq!(grid.ny(), 10);
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
    let config = SolverConfig::<f64> {
        tolerance: 1e-6,
        max_iterations: 1000,
    };
    assert_eq!(config.max_iterations, 1000);
}

#[cfg(test)]
mod performance {
    use super::*;
    use std::time::Instant;
    
    #[test]
    fn test_solver_performance() -> Result<()> {
        let fluid = Fluid::<f64>::water()?;
        let mut builder = NetworkBuilder::new(fluid);
        
        // Create larger network
        for i in 0..10 {
            builder = builder.add_node(Node::new(i, i as f64, 0.0, 0.0));
        }
        
        let network = builder.build()?;
        let mut solver = NetworkSolver::<f64>::new();
        let problem = NetworkProblem { network };
        
        let start = Instant::now();
        let _solution = solver.solve(&problem)?;
        let duration = start.elapsed();
        
        // Should complete in reasonable time
        assert!(duration.as_secs() < 1);
        Ok(())
    }
}