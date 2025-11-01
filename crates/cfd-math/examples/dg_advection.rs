//! Example: Solving the 1D Advection Equation with Discontinuous Galerkin Methods
//!
//! This example demonstrates how to use the Discontinuous Galerkin (DG) methods
//! to solve the 1D linear advection equation:
//!
//! ∂u/∂t + c ∂u/∂x = 0
//!
//! with periodic boundary conditions and a smooth initial condition.

use cfd_math::high_order::dg::{
    DGOperator, DGOperatorParams, DGSolver, FluxType, TimeIntegratorFactory, TimeIntegration,
    TimeIntegrationParams,
};
use nalgebra::DVector;
use std::f64::consts::PI;
use std::fs::File;
use std::io::Write;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Problem parameters
    let order = 3;                  // Polynomial order
    let num_elements = 16;          // Number of elements
    let t_final = 1.0;              // Final time
    let cfl = 0.1;                  // CFL number
    
    // Set up the DG operator
    let params = DGOperatorParams::new()
        .with_volume_flux(FluxType::Upwind)
        .with_surface_flux(FluxType::Upwind);
    
    let dg_op = DGOperator::new(order, 1, Some(params))?;
    
    // Set up the time integrator
    let mut solver = DGSolver::new(
        dg_op,
        TimeIntegratorFactory::create(TimeIntegration::SSPRK3),
        TimeIntegrationParams::new(TimeIntegration::SSPRK3)
            .with_t_final(t_final)
            .with_cfl(cfl)
            .with_verbose(true),
    )?;
    
    // Initial condition: u(x,0) = sin(πx)
    let u0 = |x: f64| DVector::from_vec(vec![(PI * x).sin()]);
    
    // Initialize the solver
    solver.initialize(u0)?;
    
    // Set up output file
    let mut file = File::create("advection_solution.csv")?;
    writeln!(file, "x,u_initial,u_final,u_exact")?;
    
    // Save initial condition
    let x_values: Vec<f64> = (-100..=100)
        .map(|i| -1.0 + 2.0 * i as f64 / 200.0)
        .collect();
    
    for &x in &x_values {
        let u_init = (PI * x).sin();
        writeln!(file, "{},{},,", x, u_init)?;
    }
    
    // Solve the equation
    solver.solve(
        |t, u| {
            // Right-hand side function: du/dt = -du/dx
            Ok(-solver.operator.compute_derivative(u)?)
        },
        Some(|t, solver| {
            // Callback for monitoring progress
            println!("Time = {:.3}, dt = {:.3e}", t, solver.current_dt());
            Ok(())
        }),
    )?;
    
    // Save final solution and exact solution
    let mut file = std::fs::OpenOptions::new()
        .append(true)
        .open("advection_solution.csv")?;
    
    for &x in &x_values {
        let u_final = solver.evaluate(x)[0];
        let u_exact = (PI * (x - t_final)).sin();
        writeln!(file, "{},{},{},{}", x, "", u_final, u_exact)?;
    }
    
    // Compute and print error
    let (l2_error, l2_norm) = solver.compute_error(|x| {
        DVector::from_vec(vec![(PI * (x - t_final)).sin()])
    })?;
    
    println!("\nSimulation complete!");
    println!("Final time: {}", t_final);
    println!("L2 error: {:.3e}", l2_error);
    println!("Relative L2 error: {:.3e}", l2_error / l2_norm);
    println!("\nResults saved to advection_solution.csv");
    println!("To visualize the results, you can use the following Python code:");
    
    println!(
        r#"
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('advection_solution.csv', delimiter=',', skiprows=1)
x = data[:, 0]
u_initial = data[:, 1]
u_final = data[:, 2]
u_exact = data[:, 3]

plt.figure(figsize=(10, 6))
plt.plot(x, u_initial, 'b-', label='Initial condition')
plt.plot(x, u_final, 'r--', label='Numerical solution')
plt.plot(x, u_exact, 'k:', label='Exact solution')
plt.xlabel('x')
plt.ylabel('u')
plt.title('1D Advection Equation Solution')
plt.legend()
plt.grid(True)
plt.savefig('advection_solution.png')
plt.show()
"#
    );
    
    Ok(())
}
