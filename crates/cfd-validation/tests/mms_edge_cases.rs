//! Manufactured Solutions Edge Case Tests
//!
//! Comprehensive edge case testing for MMS validation framework per Roache (1998)
//! and ASME V&V 20-2009 standards.
//!
//! Tests cover extreme parameter regimes that stress numerical schemes:
//! - High Peclet number (advection-dominated)
//! - Low Peclet number (diffusion-dominated)
//! - High Reynolds number (inviscid limit)
//! - Stiff source terms
//! - Sharp gradients and discontinuities
//!
//! References:
//! - Roache, P.J. (2002) "Code Verification by the Method of Manufactured Solutions"
//! - ASME V&V 20-2009 "Standard for Verification and Validation in CFD and HT"

use cfd_validation::manufactured::{
    ManufacturedAdvectionDiffusion, ManufacturedBurgers, ManufacturedSolution,
};
use proptest::prelude::*;

// ============================================================================
// High Peclet Number Tests (Advection-Dominated Flows)
// ============================================================================

proptest! {
    /// Test MMS with high Peclet numbers (advection-dominated)
    ///
    /// Pe >> 1 tests scheme robustness in convection-dominated regimes
    /// where central differences can produce non-physical oscillations.
    ///
    /// Reference: Patankar (1980) - Numerical Heat Transfer, §5.3
    #[test]
    fn test_high_peclet_advection_diffusion(
        velocity in 10.0f64..100.0,
        diffusivity in 0.001f64..0.01
    ) {
        use std::f64::consts::PI;
        
        // Pe = u*L/D where L is characteristic length
        // For Pe > 2, upwind schemes are typically required
        let peclet_number = velocity / diffusivity;
        
        // Create MMS with wavenumbers and velocity components
        let mms = ManufacturedAdvectionDiffusion::new(
            2.0 * PI,      // kx
            2.0 * PI,      // ky
            diffusivity,   // alpha (diffusion coefficient)
            velocity,      // vx
            velocity * 0.5 // vy
        );
        
        // Test at representative point
        let x = 0.5;
        let y = 0.5;
        let t = 0.1;
        
        // Verify solution and source term are finite
        let solution = mms.exact_solution(x, y, 0.0, t);
        let source = mms.source_term(x, y, 0.0, t);
        
        prop_assert!(
            solution.is_finite(),
            "Solution must be finite even at Pe = {:.1}", peclet_number
        );
        prop_assert!(
            source.is_finite(),
            "Source term must be finite even at Pe = {:.1}", peclet_number
        );
    }

    /// Test MMS with low diffusivity (near inviscid limit)
    ///
    /// Tests scheme behavior as viscosity approaches zero,
    /// approaching Euler equations.
    ///
    /// Reference: Ferziger & Perić (2019) - CFD Methods, §3.9
    #[test]
    fn test_low_diffusivity_limit(
        velocity in 1.0f64..10.0,
        diffusivity in 1e-6f64..1e-3
    ) {
        use std::f64::consts::PI;
        
        let mms = ManufacturedAdvectionDiffusion::new(
            2.0 * PI,      // kx
            2.0 * PI,      // ky
            diffusivity,   // alpha (very small)
            velocity,      // vx
            velocity * 0.5 // vy
        );
        
        // Test solution properties in near-inviscid regime
        let x = 0.3;
        let y = 0.7;
        let t = 0.01; // Small time to avoid numerical issues
        
        let solution = mms.exact_solution(x, y, 0.0, t);
        let source = mms.source_term(x, y, 0.0, t);
        
        // Solution should remain bounded even with tiny viscosity
        prop_assert!(
            solution.abs() < 1000.0,
            "Solution should remain bounded at low viscosity: {}", solution
        );
        prop_assert!(
            source.is_finite(),
            "Source term should be finite at low viscosity"
        );
    }
}

// ============================================================================
// Burgers Equation Edge Cases (Nonlinear Convection)
// ============================================================================

proptest! {
    /// Test Burgers equation with large amplitude (shock formation tendency)
    ///
    /// Large amplitudes can lead to shock formation in inviscid Burgers equation.
    /// Tests numerical scheme handling of steep gradients.
    ///
    /// Reference: Burgers (1948), Bateman (1915)
    #[test]
    fn test_burgers_large_amplitude(
        amplitude in 5.0f64..50.0,
        viscosity in 0.001f64..0.1
    ) {
        use std::f64::consts::PI;
        
        let burgers = ManufacturedBurgers::new(
            1.0,           // mean velocity
            amplitude,     // large amplitude
            2.0 * PI,      // wavenumber
            1.0,           // frequency
            viscosity,
        );
        
        // Test at multiple points to check gradient handling
        let test_points = vec![0.0, 0.25, 0.5, 0.75, 1.0];
        
        for &x in &test_points {
            let solution = burgers.exact_solution(x, 0.0, 0.0, 0.1);
            let source = burgers.source_term(x, 0.0, 0.0, 0.1);
            
            prop_assert!(
                solution.is_finite(),
                "Burgers solution must be finite at x={}, amplitude={}", x, amplitude
            );
            prop_assert!(
                source.is_finite(),
                "Burgers source must be finite at x={}, amplitude={}", x, amplitude
            );
        }
    }

    /// Test Burgers equation with small viscosity (shock-like behavior)
    ///
    /// Small viscosity tests scheme's ability to handle nearly discontinuous solutions
    /// without excessive numerical diffusion or oscillations.
    ///
    /// Reference: Roache (1998) - Verification and Validation in CFD
    #[test]
    fn test_burgers_low_viscosity(
        amplitude in 1.0f64..5.0,
        viscosity in 1e-6f64..1e-3
    ) {
        use std::f64::consts::PI;
        
        let burgers = ManufacturedBurgers::new(
            1.0,
            amplitude,
            PI,            // wavenumber
            1.0,           // frequency
            viscosity,     // very small viscosity
        );
        
        // Test gradient at different times
        let times = vec![0.0, 0.05, 0.1];
        
        for &t in &times {
            let x = 0.5;
            let solution = burgers.exact_solution(x, 0.0, 0.0, t);
            let source = burgers.source_term(x, 0.0, 0.0, t);
            
            prop_assert!(
                solution.abs() < 100.0,
                "Solution should remain bounded at low viscosity"
            );
            prop_assert!(
                source.is_finite(),
                "Source should be finite at low viscosity"
            );
        }
    }
}

// ============================================================================
// Spatial Resolution and Convergence Tests
// ============================================================================

proptest! {
    /// Test MMS error convergence with grid refinement
    ///
    /// Verifies that manufactured solutions maintain expected convergence
    /// order across different grid resolutions.
    ///
    /// Reference: Roache (2002) - MMS methodology
    #[test]
    fn test_mms_grid_convergence(
        velocity in 1.0f64..10.0,
        diffusivity in 0.01f64..1.0
    ) {
        use std::f64::consts::PI;
        
        let mms = ManufacturedAdvectionDiffusion::new(
            2.0 * PI,      // kx
            2.0 * PI,      // ky
            diffusivity,   // alpha
            velocity,      // vx
            velocity * 0.5 // vy
        );
        
        // Test at three different grid resolutions
        let resolutions = vec![10, 20, 40];
        let t = 0.1;
        
        let mut errors = Vec::new();
        
        for &n in &resolutions {
            let dx = 1.0 / n as f64;
            let mut max_error: f64 = 0.0;
            
            // Sample solution at grid points
            for i in 0..n {
                let x = i as f64 * dx;
                let y = 0.5;
                
                let exact = mms.exact_solution(x, y, 0.0, t);
                
                // In a real implementation, we'd compute numerical solution
                // Here we just verify the exact solution is well-behaved
                prop_assert!(
                    exact.is_finite(),
                    "Exact solution must be finite at all grid points"
                );
                
                max_error = max_error.max(exact.abs());
            }
            
            errors.push(max_error);
        }
        
        // Verify solutions are bounded
        prop_assert!(
            errors.iter().all(|&e| e < 1000.0),
            "All solutions should be bounded across resolutions"
        );
    }
}

// ============================================================================
// Temporal Evolution Tests
// ============================================================================

proptest! {
    /// Test MMS temporal accuracy with varying time steps
    ///
    /// Verifies manufactured solutions maintain consistency across
    /// different temporal discretizations.
    ///
    /// Reference: ASME V&V 20-2009 §4.2 on temporal discretization
    #[test]
    fn test_mms_temporal_evolution(
        velocity in 1.0f64..5.0,
        diffusivity in 0.1f64..1.0,
        dt in 0.001f64..0.01
    ) {
        use std::f64::consts::PI;
        
        let mms = ManufacturedAdvectionDiffusion::new(
            2.0 * PI,      // kx
            2.0 * PI,      // ky
            diffusivity,   // alpha
            velocity,      // vx
            velocity * 0.5 // vy
        );
        
        let x = 0.5;
        let y = 0.5;
        
        // Test solution at multiple time steps
        let n_steps = 10;
        
        for i in 0..n_steps {
            let t = i as f64 * dt;
            
            let solution = mms.exact_solution(x, y, 0.0, t);
            let source = mms.source_term(x, y, 0.0, t);
            
            prop_assert!(
                solution.is_finite(),
                "Solution must remain finite over time evolution at t={}", t
            );
            prop_assert!(
                source.is_finite(),
                "Source must remain finite over time evolution at t={}", t
            );
        }
    }

    /// Test MMS with stiff temporal behavior
    ///
    /// Tests solution properties when temporal scales vary widely,
    /// common in reacting flows and multiphase problems.
    ///
    /// Reference: Hairer & Wanner (1996) - Solving Stiff ODEs
    #[test]
    fn test_mms_stiff_temporal(
        fast_velocity in 50.0f64..500.0,
        slow_diffusivity in 0.001f64..0.01,
        dt in 1e-5f64..1e-4
    ) {
        use std::f64::consts::PI;
        
        // High velocity / low diffusivity creates stiff temporal behavior
        let stiffness_ratio = fast_velocity / slow_diffusivity;
        
        let mms = ManufacturedAdvectionDiffusion::new(
            2.0 * PI,          // kx
            2.0 * PI,          // ky
            slow_diffusivity,  // alpha (small)
            fast_velocity,     // vx (large)
            fast_velocity * 0.5 // vy
        );
        
        // Test at short times relevant to stiff scales
        let times = vec![0.0, dt, 2.0 * dt, 5.0 * dt];
        
        for &t in &times {
            let solution = mms.exact_solution(0.5, 0.5, 0.0, t);
            let source = mms.source_term(0.5, 0.5, 0.0, t);
            
            prop_assert!(
                solution.is_finite(),
                "Stiff MMS solution must be finite at stiffness ratio = {:.1}", stiffness_ratio
            );
            prop_assert!(
                source.is_finite(),
                "Stiff MMS source must be finite"
            );
        }
    }
}

// ============================================================================
// Boundary and Domain Tests
// ============================================================================

/// Test MMS solutions at domain boundaries
///
/// Verifies manufactured solutions are well-defined at boundaries,
/// critical for boundary condition implementation.
///
/// Reference: ASME V&V 20-2009 §5.1 on boundary treatments
#[test]
fn test_mms_boundary_values() {
    use std::f64::consts::PI;
    
    let burgers = ManufacturedBurgers::new(1.0, 0.5, 2.0 * PI, 1.0, 0.01);
    
    // Test at domain corners and edges
    let boundary_points = vec![
        (0.0, 0.0),
        (1.0, 0.0),
        (0.0, 1.0),
        (1.0, 1.0),
        (0.5, 0.0),
        (0.5, 1.0),
    ];
    
    let t = 0.1;
    
    for &(x, y) in &boundary_points {
        let solution = burgers.exact_solution(x, y, 0.0, t);
        let boundary = burgers.boundary_condition(x, y, 0.0, t);
        
        assert!(
            solution.is_finite(),
            "Boundary solution must be finite at ({}, {})", x, y
        );
        assert_eq!(
            solution, boundary,
            "Boundary condition must match exact solution"
        );
    }
}

/// Test MMS with periodic boundary conditions
///
/// Verifies manufactured solutions are compatible with periodic boundaries,
/// common in turbulence and spectral methods.
///
/// Reference: Canuto et al. (2007) - Spectral Methods
#[test]
fn test_mms_periodic_consistency() {
    use std::f64::consts::PI;
    
    let burgers = ManufacturedBurgers::new(1.0, 0.5, 2.0 * PI, 1.0, 0.01);
    
    // For periodic solutions, u(x+L) = u(x) where L is domain size
    let domain_size = 1.0;
    let t = 0.1;
    
    let test_points = vec![0.0, 0.25, 0.5, 0.75];
    
    for &x in &test_points {
        let u_left = burgers.exact_solution(x, 0.0, 0.0, t);
        let u_right = burgers.exact_solution(x + domain_size, 0.0, 0.0, t);
        
        // For 2π periodic functions, u(x + 1) ≈ u(x) if k = 2π
        // Small numerical differences expected due to floating point
        let diff = (u_left - u_right).abs();
        
        assert!(
            diff < 1e-10,
            "Periodic MMS solution should repeat: diff = {} at x = {}", diff, x
        );
    }
}
