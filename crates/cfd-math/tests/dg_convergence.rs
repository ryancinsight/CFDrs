//! Convergence tests for Discontinuous Galerkin methods

use cfd_math::high_order::dg::*;
use approx::assert_relative_eq;
use nalgebra::DVector;
use std::f64::consts::PI;

#[test]
fn test_convergence_linear_advection() -> Result<(), Box<dyn std::error::Error>> {
    let orders = [1, 2, 3, 4];
    let num_elements_list = [4, 8, 16, 32];
    let t_final = 1.0;
    let cfl = 0.1;
    
    let mut results = Vec::new();
    
    for &order in &orders {
        let mut prev_error = None;
        
        for &num_elements in &num_elements_list {
            let params = DGOperatorParams::new()
                .with_volume_flux(FluxType::Upwind)
                .with_surface_flux(FluxType::Upwind);
                
            let dg_op = DGOperator::new(order, 1, Some(params))?;
            
            let mut solver = DGSolver::new(
                dg_op,
                TimeIntegratorFactory::create(TimeIntegration::SSPRK3),
                TimeIntegrationParams::new(TimeIntegration::SSPRK3)
                    .with_t_final(t_final)
                    .with_cfl(cfl),
            )?;
            
            // Initial condition u(x,0) = sin(πx)
            solver.initialize(|x| DVector::from_vec(vec![(PI * x).sin()]))?;
            
            // Solve
            solver.solve(
                |t, u| {
                    // For linear advection: du/dt = -du/dx
                    Ok(-solver.operator.compute_derivative(u)?)
                },
                None::<fn(_, _) -> _>,
            )?;
            
            // Compute error
            let (error, _) = solver.compute_error(|x| {
                DVector::from_vec(vec![(PI * (x - t_final)).sin()])
            })?;
            
            if let Some(prev) = prev_error {
                let h = 2.0 / num_elements as f64;
                let rate = f64::log2(prev / error);
                results.push((order, h, error, rate));
                
                println!(
                    "Order {}: h = {:.3e}, error = {:.3e}, rate = {:.2}",
                    order, h, error, rate
                );
            }
            
            prev_error = Some(error);
        }
    }
    
    // Verify convergence rates
    for (order, h, error, rate) in &results {
        let expected_rate = *order as f64 + 1.0;
        assert!(
            *rate >= expected_rate * 0.8,  // Allow 20% tolerance
            "Convergence rate {} below expected {} for order {} (h = {})",
            rate, expected_rate, order, h
        );
    }
    
    Ok(())
}

#[test]
fn test_mass_conservation() -> Result<(), Box<dyn std::error::Error>> {
    let order = 3;
    let num_elements = 16;
    let t_final = 1.0;
    let cfl = 0.1;
    
    // Set up DG operator with periodic BCs
    let params = DGOperatorParams::new()
        .with_volume_flux(FluxType::Central)
        .with_surface_flux(FluxType::Central);
        
    let dg_op = DGOperator::new(order, 1, Some(params))?;
    
    let mut solver = DGSolver::new(
        dg_op,
        TimeIntegratorFactory::create(TimeIntegration::SSPRK3),
        TimeIntegrationParams::new(TimeIntegration::SSPRK3)
            .with_t_final(t_final)
            .with_cfl(cfl),
    )?;
    
    // Initial condition with non-zero mean
    let u0 = |x: f64| DVector::from_vec(vec![1.0 + 0.5 * (2.0 * PI * x).sin()]);
    
    solver.initialize(u0)?;
    
    // Get initial mass
    let mass_initial = solver.compute_mass()?;
    
    // Solve
    solver.solve(
        |t, u| {
            // For linear advection: du/dt = -du/dx
            Ok(-solver.operator.compute_derivative(u)?)
        },
        None::<fn(_, _) -> _>,
    )?;
    
    // Get final mass
    let mass_final = solver.compute_mass()?;
    
    // Mass should be conserved exactly for periodic BCs
    assert_abs_diff_eq!(
        mass_initial,
        mass_final,
        epsilon = 1e-10,
        abs_diff_all = 1e-10,
        relative = 1e-10
    );
    
    Ok(())
}
