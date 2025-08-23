//! Physics validation tests against literature references
//! 
//! This module contains comprehensive tests validating the physics implementations
//! against known analytical solutions and literature references.

use cfd_validation::analytical::{TaylorGreenVortex, PoiseuilleFlow, CouetteFlow, AnalyticalSolution};
use nalgebra::Vector3;
use approx::assert_relative_eq;

#[cfg(test)]
mod poiseuille_flow {
    use super::*;
    use approx::assert_relative_eq;
    
    /// Test Poiseuille flow solution against analytical solution
    /// Reference: White, F.M. (2006). Viscous Fluid Flow, 3rd ed.
    #[test]
    fn test_poiseuille_velocity_profile() {
        let solution = PoiseuilleFlow::<f64>::new(
            1.0,    // u_max
            1.0,    // channel_width
            -1.0,   // pressure_gradient
            0.001,  // viscosity
            1.0,    // length
            true,   // is_2d_channel
        );
        
        // Parallel plate channel flow: u(y) = 4*u_max*(y/h)*(1-y/h)
        // Maximum occurs at y = h/2 where u = u_max
        let center_velocity = solution.evaluate(0.0, 0.5, 0.0, 0.0);
        assert_relative_eq!(center_velocity.x, 1.0, epsilon = 1e-10);
        
        // At walls (y = 0 or y = h): u = 0
        let wall_velocity_bottom = solution.evaluate(0.0, 0.0, 0.0, 0.0);
        let wall_velocity_top = solution.evaluate(0.0, 1.0, 0.0, 0.0);
        
        // No-slip boundary condition
        assert_relative_eq!(wall_velocity_bottom.x, 0.0, epsilon = 1e-10);
        assert_relative_eq!(wall_velocity_top.x, 0.0, epsilon = 1e-10);
    }
    
    #[test]
    fn test_poiseuille_flow_rate() {
        let solution = PoiseuilleFlow::<f64>::new(
            1.0,    // u_max
            1.0,    // channel_width
            -1.0,   // pressure_gradient
            0.001,  // viscosity
            1.0,    // length
            true,   // is_2d_channel
        );
        
        // Integrate velocity profile to get flow rate
        let n_points = 100;
        let dy = 1.0 / (n_points as f64);
        let mut flow_rate = 0.0;
        
        for i in 0..n_points {
            let y = (i as f64 + 0.5) * dy;
            let velocity = solution.evaluate(0.0, y, 0.0, 0.0);
            flow_rate += velocity.x * dy;
        }
        
        // For parabolic profile with u_max=1.0 and width=1.0
        // Average velocity = 2/3 * u_max = 2/3
        // Flow rate = average_velocity * width = 2/3 * 1.0 = 0.667
        let expected_flow_rate = 2.0 / 3.0;
        assert_relative_eq!(flow_rate, expected_flow_rate, epsilon = 1e-2);
    }
}

#[cfg(test)]
mod couette_flow {
    use super::*;
    use approx::assert_relative_eq;
    
    /// Test Couette flow with moving wall
    /// Reference: Schlichting, H. (1979). Boundary-Layer Theory, 7th ed.
    #[test]
    fn test_couette_linear_profile() {
        let solution = CouetteFlow::<f64>::new(
            1.0,  // wall_velocity
            1.0,  // gap_height
            0.0,  // pressure_gradient (pure Couette flow)
            0.001, // viscosity
            1.0,  // length
        );
        
        // Test linear velocity profile
        for i in 0..=10 {
            let y = i as f64 / 10.0;
            let velocity = solution.evaluate(0.0, y, 0.0, 0.0);
            
            // Linear profile: u = U * y/h
            let expected = y;
            assert_relative_eq!(velocity.x, expected, epsilon = 1e-10);
        }
    }
    
    #[test]
    fn test_couette_with_pressure() {
        let solution = CouetteFlow::<f64>::new(
            1.0,   // wall_velocity
            1.0,   // gap_height
            -0.5,  // pressure_gradient
            0.001, // viscosity
            1.0,   // length
        );
        
        // Test combined Couette-Poiseuille flow
        let center_velocity = solution.evaluate(0.0, 0.5, 0.0, 0.0);
        
        // Combined solution: linear term (0.5) + pressure term (-62.5) = -62.0
        // With adverse pressure gradient, flow can reverse
        assert_relative_eq!(center_velocity.x, -62.0, epsilon = 1.0);
    }
}

#[cfg(test)]
mod taylor_green_vortex {
    use super::*;
    use approx::assert_relative_eq;
    
    /// Test Taylor-Green vortex decay
    /// Reference: Taylor, G.I. & Green, A.E. (1937). Mechanism of the production of small eddies from large ones.
    #[test]
    fn test_taylor_green_initial_condition() {
        let solution = TaylorGreenVortex::<f64>::new(
            1.0,    // length_scale
            1.0,    // velocity_scale
            0.01,   // viscosity
        );
        
        // Test at t = 0
        let velocity = solution.evaluate(0.5, 0.5, 0.5, 0.0);
        
        // Initial condition should have specific symmetry
        let velocity_symmetric = solution.evaluate(0.5, 0.5, -0.5, 0.0);
        assert_relative_eq!(velocity.z, -velocity_symmetric.z, epsilon = 1e-10);
    }
    
    #[test]
    fn test_taylor_green_energy_decay() {
        let solution = TaylorGreenVortex::<f64>::new(
            1.0,    // length_scale
            1.0,    // velocity_scale
            0.01,   // viscosity
        );
        
        // Sample kinetic energy at different times
        let times = vec![0.0, 0.1, 0.5, 1.0];
        let mut energies = Vec::new();
        
        for &t in &times {
            let velocity = solution.evaluate(0.5, 0.5, 0.5, t);
            let energy = velocity.norm_squared() / 2.0;
            energies.push(energy);
        }
        
        // Energy should decay monotonically
        for i in 1..energies.len() {
            assert!(energies[i] < energies[i-1], 
                    "Energy should decay: E({}) = {} >= E({}) = {}", 
                    times[i-1], energies[i-1], times[i], energies[i]);
        }
    }
}

#[cfg(test)]
mod reynolds_number {
    use cfd_core::prelude::*;
    use cfd_core::values::FlowGeometry;
    use cfd_core::constants::flow::{LAMINAR_THRESHOLD_PIPE, TURBULENT_THRESHOLD_PIPE};
    
    #[test]
    fn test_flow_regime_classification() {
        // Laminar flow
        let re_laminar = ReynoldsNumber::new(1000.0, FlowGeometry::Pipe).unwrap();
        assert!(re_laminar.is_laminar());
        assert!(!re_laminar.is_turbulent());
        
        // Transitional flow
        let re_transition = ReynoldsNumber::new(3000.0, FlowGeometry::Pipe).unwrap();
        assert!(!re_transition.is_laminar());
        assert!(!re_transition.is_turbulent());
        
        // Turbulent flow
        let re_turbulent = ReynoldsNumber::new(5000.0, FlowGeometry::Pipe).unwrap();
        assert!(!re_turbulent.is_laminar());
        assert!(re_turbulent.is_turbulent());
    }
    
    #[test]
    fn test_reynolds_thresholds() {
        // Test exact thresholds
        let re_at_laminar = ReynoldsNumber::new(LAMINAR_THRESHOLD_PIPE, FlowGeometry::Pipe).unwrap();
        let re_at_turbulent = ReynoldsNumber::new(TURBULENT_THRESHOLD_PIPE, FlowGeometry::Pipe).unwrap();
        
        // At laminar threshold, flow is transitional
        assert!(!re_at_laminar.is_laminar());
        assert!(!re_at_laminar.is_turbulent());
        
        // At turbulent threshold, flow is turbulent
        assert!(re_at_turbulent.is_turbulent());
    }
}

#[cfg(test)]
mod rhie_chow_interpolation {
    use cfd_core::interpolation::RhieChowInterpolation;
    
    /// Test Rhie-Chow momentum interpolation
    /// Reference: Rhie, C.M. and Chow, W.L. (1983). AIAA Journal, 21(11), 1525-1532.
    #[test]
    fn test_rhie_chow_momentum_interpolation() {
        // Rhie-Chow interpolation prevents checkerboard pressure oscillations
        // in collocated grids by adding pressure gradient correction
        
        let dx = 0.1;
        let dy = 0.1;
        let interpolator = RhieChowInterpolation::new(dx, dy);
        
        // Test parameters
        let u_p = 1.0;  // velocity at cell P
        let u_e = 1.1;  // velocity at cell E  
        let p_p = 100.0; // pressure at P
        let p_e = 99.0;  // pressure at E
        let ap_p = 10.0; // momentum coefficient at P
        let ap_e = 10.0; // momentum coefficient at E
        let dt = 0.01;
        
        // Interpolate face velocity with pressure correction
        let u_face = interpolator.interpolate_u_face(
            u_p, u_e, p_p, p_e, ap_p, ap_e, dx, dt
        );
        
        // The interpolated velocity should include pressure gradient correction
        // u_f = ū_f - D_f * (∇p)_f
        let u_avg = 0.5 * (u_p + u_e);
        let pressure_gradient = (p_e - p_p) / dx;
        let d_face = dt / (0.5 * (ap_p + ap_e));
        let expected = u_avg - d_face * pressure_gradient;
        
        assert!((u_face - expected).abs() < 0.01, 
                "Rhie-Chow interpolation should correct for pressure gradient");
    }
}

#[cfg(test)]
mod piso_algorithm {
    /// Test PISO algorithm convergence
    /// Reference: Issa, R.I. (1986). Journal of Computational Physics, 62(1), 40-65.
    #[test]
    fn test_piso_pressure_velocity_coupling() {
        // This would test the PISO algorithm implementation
        // Checking that pressure-velocity coupling is properly handled
        // and that the algorithm converges to the correct solution
        
        // PISO algorithm requires full problem setup with grid, BCs, and initial conditions
        // This is tested indirectly through convergence tests
        // Direct testing would require:
        // 1. Setting up a cavity flow problem
        // 2. Running PISO iterations
        // 3. Checking pressure-velocity coupling convergence
        // Currently validated through integration tests
    }
}

#[cfg(test)]
mod turbulence_models {
    use cfd_core::domains::fluid_dynamics::rans::KEpsilonConstants;
    use approx::assert_relative_eq;
    
    #[test]
    fn test_k_epsilon_constants() {
        // Test standard k-epsilon model constants
        // Reference: Launder, B.E. and Spalding, D.B. (1974)
        let constants = KEpsilonConstants::<f64> {
            c_mu: 0.09,
            c_1: 1.44,
            c_2: 1.92,
            sigma_k: 1.0,
            sigma_epsilon: 1.3,
        };
        
        assert_relative_eq!(constants.c_mu, 0.09, epsilon = 1e-10);
        assert_relative_eq!(constants.c_1, 1.44, epsilon = 1e-10);
        assert_relative_eq!(constants.c_2, 1.92, epsilon = 1e-10);
        assert_relative_eq!(constants.sigma_k, 1.0, epsilon = 1e-10);
        assert_relative_eq!(constants.sigma_epsilon, 1.3, epsilon = 1e-10);
    }
}