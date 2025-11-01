//! Validation tests for Discontinuous Galerkin methods
//!
//! This module contains validation tests that compare numerical solutions
//! against analytical solutions or established benchmarks from the literature.

use crate::high_order::dg::*;
use approx::assert_relative_eq;
use nalgebra::{DMatrix, DVector};
use std::f64::consts::PI;

/// Validates the DG implementation against known analytical solutions
pub struct DGValidation;

impl DGValidation {
    /// Test case 1: Linear advection with smooth initial condition
    /// Solves u_t + u_x = 0 with u(x,0) = sin(πx) on [-1,1]
    /// Exact solution: u(x,t) = sin(π(x-t))
    pub fn linear_advection_smooth() -> Result<(), DGError> {
        let order = 4;
        let num_elements = 8;
        let t_final = 1.0;
        let cfl = 0.1;
        
        // Set up DG operator
        let params = DGOperatorParams::new()
            .with_volume_flux(FluxType::Upwind)
            .with_surface_flux(FluxType::Upwind);
            
        let mut dg_op = DGOperator::new(order, 1, Some(params))?;
        
        // Initial condition
        let u0 = |x: f64| DVector::from_vec(vec![(PI * x).sin()]);
        
        // Set up solver
        let mut solver = DGSolver::new(
            dg_op,
            TimeIntegratorFactory::create(TimeIntegration::SSPRK3),
            TimeIntegrationParams::new(TimeIntegration::SSPRK3)
                .with_t_final(t_final)
                .with_cfl(cfl),
        )?;
        
        solver.initialize(u0)?;
        
        // Solve
        solver.solve(
            |t, u| {
                // For linear advection: du/dt = -du/dx
                Ok(-solver.operator.compute_derivative(u)?)
            },
            None::<fn(_, _) -> _>,
        )?;
        
        // Compute error
        let (l2_error, _) = solver.compute_error(|x| {
            DVector::from_vec(vec![(PI * (x - t_final)).sin()])
        })?;
        
        // Check that the error is small
        assert!(l2_error < 1e-4, "L2 error {} too large", l2_error);
        
        Ok(())
    }
    
    /// Test case 2: Burgers' equation with smooth initial condition
    /// Solves u_t + (u²/2)_x = 0 with u(x,0) = 0.5 + sin(πx) on [-1,1]
    pub fn burgers_equation() -> Result<(), DGError> {
        let order = 3;
        let num_elements = 16;
        let t_final = 0.5;  // Before shock formation
        let cfl = 0.1;
        
        // Set up DG operator with Lax-Friedrichs flux
        let params = DGOperatorParams::new()
            .with_volume_flux(FluxType::LaxFriedrichs)
            .with_surface_flux(FluxType::LaxFriedrichs);
            
        let mut dg_op = DGOperator::new(order, 1, Some(params))?;
        
        // Initial condition
        let u0 = |x: f64| DVector::from_vec(vec![0.5 + (PI * x).sin()]);
        
        // Set up solver
        let mut solver = DGSolver::new(
            dg_op,
            TimeIntegratorFactory::create(TimeIntegration::SSPRK3),
            TimeIntegrationParams::new(TimeIntegration::SSPRK3)
                .with_t_final(t_final)
                .with_cfl(cfl),
        )?;
        
        solver.initialize(u0)?;
        
        // Solve
        solver.solve(
            |t, u| {
                // For Burgers' equation: du/dt = -u * du/dx
                let du_dx = solver.operator.compute_derivative(u)?;
                Ok(-u.component_mul(&du_dx))
            },
            None::<fn(_, _) -> _>,
        )?;
        
        // Check conservation properties
        let mass_initial = solver.initial_mass();
        let mass_final = solver.compute_mass()?;
        
        assert_relative_eq!(
            mass_initial, 
            mass_final, 
            epsilon = 1e-10,
            "Mass not conserved: initial = {}, final = {}", 
            mass_initial, 
            mass_final
        );
        
        Ok(())
    }
    
    /// Test case 3: Convergence test for linear advection
    /// Verifies that the method achieves the expected order of accuracy
    pub fn convergence_test() -> Result<Vec<(usize, f64)>, DGError> {
        let orders = [1, 2, 4];
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
                    
                let mut dg_op = DGOperator::new(order, 1, Some(params))?;
                
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
                    results.push((order, rate));
                    
                    println!(
                        "Order {}: h = {:.3e}, error = {:.3e}, rate = {:.2}",
                        order, h, error, rate
                    );
                }
                
                prev_error = Some(error);
            }
        }
        
        Ok(results)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    
    #[test]
    fn test_linear_advection_smooth() {
        DGValidation::linear_advection_smooth().unwrap();
    }
    
    #[test]
    fn test_burgers_equation() {
        DGValidation::burgers_equation().unwrap();
    }
    
    #[test]
    fn test_convergence() {
        let results = DGValidation::convergence_test().unwrap();
        
        // Check that the convergence rate is at least (order + 1/2) for each order
        for (order, rate) in results {
            let expected = order as f64 + 0.5;
            assert!(
                rate >= expected * 0.9,  // Allow 10% tolerance
                "Convergence rate {} below expected {} for order {}",
                rate, expected, order
            );
        }
    }
}
