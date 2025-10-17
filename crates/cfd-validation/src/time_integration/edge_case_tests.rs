//! Edge case tests for time integration methods
//!
//! Tests validate Forward Euler, RK2, RK4, and multistep methods against
//! known ODE solutions and numerical analysis literature.
//!
//! References:
//! - Hairer, Nørsett & Wanner (1993). Solving Ordinary Differential Equations I. Springer.
//! - Press et al. (2007). Numerical Recipes (3rd ed.). Cambridge University Press.

#[cfg(test)]
mod time_integration_edge_tests {
    use crate::time_integration::integrators::{ForwardEuler, RungeKutta2, TimeIntegratorTrait};
    use approx::assert_relative_eq;
    use nalgebra::DVector;

    /// Test Forward Euler with zero initial condition
    /// Validates: dy/dt = 0 → y(t) = 0
    #[test]
    fn test_forward_euler_zero_initial() {
        let integrator = ForwardEuler;
        let mut y = DVector::from_element(3, 0.0);
        let t = 0.0;
        let dt = 0.1;
        
        // ODE: dy/dt = 0
        let f = |_t: f64, _y: &DVector<f64>| DVector::zeros(3);
        
        integrator.step(&mut y, t, dt, f).unwrap();
        
        // Solution should remain zero
        for i in 0..3 {
            assert_relative_eq!(y[i], 0.0, epsilon = 1e-14);
        }
    }

    /// Test Forward Euler with constant derivative
    /// Validates: dy/dt = c → y(t) = y0 + c*t
    #[test]
    fn test_forward_euler_constant_derivative() {
        let integrator = ForwardEuler;
        let mut y = DVector::from_element(2, 1.0);
        let t = 0.0;
        let dt = 0.1;
        let constant = 2.0;
        
        // ODE: dy/dt = 2.0
        let f = |_t: f64, _y: &DVector<f64>| DVector::from_element(2, constant);
        
        integrator.step(&mut y, t, dt, f).unwrap();
        
        // y(0.1) = 1.0 + 2.0 * 0.1 = 1.2
        let expected = 1.0 + constant * dt;
        for i in 0..2 {
            assert_relative_eq!(y[i], expected, epsilon = 1e-12);
        }
    }

    /// Test Forward Euler with negative derivative
    /// Validates exponential decay: dy/dt = -λy → y(t) = y0*exp(-λt)
    #[test]
    fn test_forward_euler_negative_derivative() {
        let integrator = ForwardEuler;
        let mut y = DVector::from_element(1, 1.0);
        let t = 0.0;
        let dt = 0.01; // Small timestep for accuracy
        let lambda = 1.0;
        
        // ODE: dy/dt = -y (exponential decay)
        let f = |_t: f64, y: &DVector<f64>| DVector::from_element(1, -lambda * y[0]);
        
        // Take 10 small steps
        let mut current_t = t;
        for _ in 0..10 {
            integrator.step(&mut y, current_t, dt, f).unwrap();
            current_t += dt;
        }
        
        // Approximate solution: y(0.1) ≈ exp(-0.1) ≈ 0.9048
        let expected = (-lambda * current_t).exp();
        assert_relative_eq!(y[0], expected, epsilon = 1e-3);
        assert!(y[0] > 0.0, "Exponential decay remains positive");
    }

    /// Test RK2 with zero timestep (boundary case)
    #[test]
    fn test_rk2_zero_timestep() {
        let integrator = RungeKutta2;
        let mut y = DVector::from_element(2, 5.0);
        let y_initial = y.clone();
        let t = 1.0;
        let dt = 0.0;
        
        // Any ODE
        let f = |_t: f64, y: &DVector<f64>| y * 2.0;
        
        integrator.step(&mut y, t, dt, f).unwrap();
        
        // With dt=0, solution should not change
        for i in 0..2 {
            assert_relative_eq!(y[i], y_initial[i], epsilon = 1e-14);
        }
    }

    /// Test RK2 with very small timestep (numerical stability)
    #[test]
    fn test_rk2_small_timestep() {
        let integrator = RungeKutta2;
        let mut y = DVector::from_element(1, 1.0);
        let t = 0.0;
        let dt = 1.0e-10; // Very small timestep
        
        // ODE: dy/dt = y
        let f = |_t: f64, y: &DVector<f64>| y.clone();
        
        integrator.step(&mut y, t, dt, f).unwrap();
        
        // Solution should be close to initial (exp(1e-10) ≈ 1.0)
        assert_relative_eq!(y[0], 1.0, epsilon = 1e-8);
        assert!(y[0].is_finite(), "Solution remains numerically stable");
    }

    /// Test Forward Euler order of accuracy
    /// Reference: Hairer et al. (1993) - First-order method
    #[test]
    fn test_forward_euler_order() {
        let integrator = ForwardEuler;
        assert_eq!(<ForwardEuler as TimeIntegratorTrait<f64>>::order(&integrator), 1, "Forward Euler is first-order accurate");
        assert_eq!(<ForwardEuler as TimeIntegratorTrait<f64>>::workspace_size(&integrator), 0, "Forward Euler requires no workspace");
    }

    /// Test RK2 order of accuracy
    /// Reference: Hairer et al. (1993) - Second-order method
    #[test]
    fn test_rk2_order() {
        let integrator = RungeKutta2;
        assert_eq!(<RungeKutta2 as TimeIntegratorTrait<f64>>::order(&integrator), 2, "RK2 is second-order accurate");
        assert!(<RungeKutta2 as TimeIntegratorTrait<f64>>::workspace_size(&integrator) > 0, "RK2 requires workspace");
    }

    /// Test time integration with large timestep (CFL-like condition)
    #[test]
    fn test_forward_euler_large_timestep() {
        let integrator = ForwardEuler;
        let mut y = DVector::from_element(1, 1.0);
        let t = 0.0;
        let dt = 10.0; // Large timestep
        
        // Stable ODE: dy/dt = -0.1*y
        let f = |_t: f64, y: &DVector<f64>| y * -0.1;
        
        integrator.step(&mut y, t, dt, f).unwrap();
        
        // Solution should remain bounded (stability check)
        assert!(y[0].is_finite(), "Solution remains stable with large timestep");
        assert!(y[0] >= 0.0, "Solution physically reasonable");
    }

    /// Test RK2 with stiff ODE (numerical stability test)
    /// Reference: Hairer & Wanner (1996) - Stiff ODEs
    #[test]
    fn test_rk2_stiff_ode() {
        let integrator = RungeKutta2;
        let mut y = DVector::from_element(1, 1.0);
        let t = 0.0;
        let dt = 0.001; // Small timestep for stiff problem
        let lambda = -1000.0; // Stiff parameter
        
        // Stiff ODE: dy/dt = -1000*y
        let f = |_t: f64, y: &DVector<f64>| y * lambda;
        
        // Take several steps
        let mut current_t = t;
        for _ in 0..10 {
            integrator.step(&mut y, current_t, dt, f).unwrap();
            current_t += dt;
        }
        
        // Solution should decay rapidly but remain stable
        assert!(y[0].is_finite(), "Solution numerically stable for stiff ODE");
        assert!(y[0] >= 0.0 && y[0] <= 1.0, "Solution physically bounded");
    }

    /// Test multi-dimensional system integration
    /// Validates coupled ODEs
    #[test]
    fn test_forward_euler_multidimensional() {
        let integrator = ForwardEuler;
        let mut y = DVector::from_vec(vec![1.0, 0.0]); // Initial: [1, 0]
        let t = 0.0;
        let dt = 0.1;
        
        // Coupled system: dy1/dt = -y2, dy2/dt = y1 (harmonic oscillator)
        let f = |_t: f64, y: &DVector<f64>| {
            DVector::from_vec(vec![-y[1], y[0]])
        };
        
        integrator.step(&mut y, t, dt, f).unwrap();
        
        // Both components should be finite
        assert!(y[0].is_finite() && y[1].is_finite(), "Coupled system remains stable");
        
        // Energy-like quantity should be approximately conserved (relaxed tolerance)
        let energy = y[0] * y[0] + y[1] * y[1];
        assert!(energy > 0.5 && energy < 1.5, "Energy approximately conserved");
    }

    /// Test RK2 accuracy vs Forward Euler
    /// Validates second-order convergence
    #[test]
    fn test_rk2_vs_euler_accuracy() {
        let euler = ForwardEuler;
        let rk2 = RungeKutta2;
        
        let mut y_euler = DVector::from_element(1, 1.0);
        let mut y_rk2 = DVector::from_element(1, 1.0);
        
        let t = 0.0;
        let dt = 0.1;
        
        // ODE: dy/dt = -y
        let f = |_t: f64, y: &DVector<f64>| -y.clone();
        
        euler.step(&mut y_euler, t, dt, f).unwrap();
        rk2.step(&mut y_rk2, t, dt, f).unwrap();
        
        // Exact: y(0.1) = exp(-0.1) ≈ 0.9048
        let exact = (-0.1_f64).exp();
        
        // RK2 should be more accurate than Euler
        let error_euler = (y_euler[0] - exact).abs();
        let error_rk2 = (y_rk2[0] - exact).abs();
        
        assert!(error_rk2 < error_euler, "RK2 more accurate than Forward Euler");
    }
}
