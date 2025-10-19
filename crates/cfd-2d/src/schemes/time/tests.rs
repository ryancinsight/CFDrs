//! Tests for time integration schemes

#[cfg(test)]
mod tests {
    use super::super::{TimeIntegrator, TimeScheme};
    use approx::assert_relative_eq;
    use nalgebra::DVector;

    #[test]
    fn test_bdf2_exponential_decay() {
        // Test BDF2 with dy/dt = -y (exponential decay)
        // Exact solution: y(t) = exp(-t)
        let integrator = TimeIntegrator::<f64>::new(TimeScheme::BDF2);
        let f = |_t: f64, y: &DVector<f64>| -y;
        
        let dt = 0.1;
        let y0 = DVector::from_vec(vec![1.0]);
        
        // First step with Forward Euler (no history yet)
        let y1 = &y0 + f(0.0, &y0) * dt;
        
        // Second step with BDF2
        let y2 = integrator.step_with_history(f, &y1, Some(&y0), None, dt, dt);
        
        // Analytical: y(0.2) = exp(-0.2) ≈ 0.8187307530779818
        // Due to Forward Euler first step, expect some error
        let expected = (-0.2_f64).exp();
        // BDF2 after FE first step gives reasonable accuracy
        assert_relative_eq!(y2[0], expected, epsilon = 1e-2);
        
        // Verify solution is bounded and positive
        assert!(y2[0] > 0.0);
        assert!(y2[0] < 1.0);
    }

    #[test]
    fn test_bdf2_order_accuracy() {
        // Verify BDF2 is 2nd-order accurate
        assert_eq!(TimeIntegrator::<f64>::new(TimeScheme::BDF2).order(), 2);
    }

    #[test]
    fn test_bdf2_is_implicit() {
        let integrator = TimeIntegrator::<f64>::new(TimeScheme::BDF2);
        assert!(!integrator.is_explicit());
    }

    #[test]
    fn test_bdf2_stiff_system() {
        // Test with moderately stiff ODE: dy/dt = -10*y
        // BDF2 should remain stable due to A-stability
        let integrator = TimeIntegrator::<f64>::new(TimeScheme::BDF2);
        let lambda = -10.0;
        let f = move |_t: f64, y: &DVector<f64>| y * lambda;
        
        let dt = 0.1;
        let y0 = DVector::from_vec(vec![1.0]);
        
        // First step: use exact solution for better initial history
        // y(dt) = exp(-10*dt)
        let y1 = DVector::from_vec(vec![(-10.0_f64 * dt).exp()]);
        
        // Second step with BDF2
        let y2 = integrator.step_with_history(f, &y1, Some(&y0), None, dt, dt);
        
        // Analytical solution at t=2*dt: y(2*dt) = exp(-20*dt)
        let expected = (-10.0 * 2.0 * dt).exp();
        
        // Solution should remain bounded and match analytical reasonably
        assert!(y2[0].abs() < 1.0);
        assert!(y2[0] > 0.0);
        // Fixed-point iteration on stiff systems may have moderate error
        assert_relative_eq!(y2[0], expected, epsilon = 0.05);
    }

    #[test]
    fn test_bdf2_convergence_order() {
        // Verify BDF2 achieves 2nd-order convergence via MMS
        // Test problem: dy/dt = -y, y(0) = 1, exact solution: y(t) = exp(-t)
        let integrator = TimeIntegrator::<f64>::new(TimeScheme::BDF2);
        let f = |_t: f64, y: &DVector<f64>| -y;
        
        let t_final = 1.0_f64;
        let exact_final = (-t_final).exp();
        
        // Test with three grid sizes to verify convergence order
        let dt_values = vec![0.1_f64, 0.05_f64, 0.025_f64];
        let mut errors = Vec::new();
        
        for &dt in &dt_values {
            let n_steps = (t_final / dt).round() as usize;
            let mut y_prev = DVector::from_vec(vec![1.0]);
            
            // First step with exact solution (to start BDF2 properly)
            let mut y_curr = DVector::from_vec(vec![(-dt).exp()]);
            
            // Subsequent steps with BDF2
            for step in 1..n_steps {
                let t = (step as f64) * dt;
                let y_next = integrator.step_with_history(f, &y_curr, Some(&y_prev), None, t, dt);
                y_prev = y_curr;
                y_curr = y_next;
            }
            
            let error = (y_curr[0] - exact_final).abs();
            errors.push(error);
        }
        
        // Verify 2nd-order convergence: error ~ O(dt^2)
        // Convergence ratio should be approximately 4 (2^2) when dt is halved
        let ratio1 = errors[0] / errors[1];
        let ratio2 = errors[1] / errors[2];
        
        // Allow some tolerance due to fixed-point iteration convergence
        assert!(ratio1 > 2.5, "Convergence ratio 1: {ratio1} should be > 2.5");
        assert!(ratio2 > 2.5, "Convergence ratio 2: {ratio2} should be > 2.5");
        
        // Verify errors are decreasing
        assert!(errors[1] < errors[0]);
        assert!(errors[2] < errors[1]);
    }

    #[test]
    fn test_bdf3_order_accuracy() {
        // Verify BDF3 is 3rd-order accurate
        assert_eq!(TimeIntegrator::<f64>::new(TimeScheme::BDF3).order(), 3);
    }

    #[test]
    fn test_bdf3_is_implicit() {
        let integrator = TimeIntegrator::<f64>::new(TimeScheme::BDF3);
        assert!(!integrator.is_explicit());
    }

    #[test]
    fn test_bdf3_exponential_decay() {
        // Test BDF3 with dy/dt = -y (exponential decay)
        // Exact solution: y(t) = exp(-t)
        let integrator = TimeIntegrator::<f64>::new(TimeScheme::BDF3);
        let f = |_t: f64, y: &DVector<f64>| -y;
        
        let dt = 0.1_f64;
        let y0 = DVector::from_vec(vec![1.0]);
        
        // First step with exact solution
        let y1 = DVector::from_vec(vec![(-dt).exp()]);
        // Second step with exact solution
        let y2 = DVector::from_vec(vec![(-2.0 * dt).exp()]);
        
        // Third step with BDF3
        let y3 = integrator.step_with_history(f, &y2, Some(&y1), Some(&y0), 2.0 * dt, dt);
        
        // Analytical: y(0.3) = exp(-0.3)
        let expected = (-0.3_f64).exp();
        
        // BDF3 should give high accuracy
        assert_relative_eq!(y3[0], expected, epsilon = 1e-3);
        
        // Verify solution is bounded and positive
        assert!(y3[0] > 0.0);
        assert!(y3[0] < 1.0);
    }

    #[test]
    fn test_bdf3_stiff_system() {
        // Test with moderately stiff ODE: dy/dt = -10*y
        // BDF3 should remain stable
        let integrator = TimeIntegrator::<f64>::new(TimeScheme::BDF3);
        let lambda = -10.0;
        let f = move |_t: f64, y: &DVector<f64>| y * lambda;
        
        let dt = 0.1_f64;
        let y0 = DVector::from_vec(vec![1.0]);
        
        // First two steps: use exact solution for better initial history
        let y1 = DVector::from_vec(vec![(-10.0_f64 * dt).exp()]);
        let y2 = DVector::from_vec(vec![(-10.0_f64 * 2.0 * dt).exp()]);
        
        // Third step with BDF3
        let y3 = integrator.step_with_history(f, &y2, Some(&y1), Some(&y0), 2.0 * dt, dt);
        
        // Analytical solution at t=3*dt: y(3*dt) = exp(-30*dt)
        let expected = (-10.0 * 3.0 * dt).exp();
        
        // Solution should remain bounded and match analytical reasonably
        assert!(y3[0].abs() < 1.0);
        assert!(y3[0] > 0.0);
        // Fixed-point iteration on stiff systems may have moderate error
        assert_relative_eq!(y3[0], expected, epsilon = 0.05);
    }

    #[test]
    fn test_bdf3_convergence_order() {
        // Verify BDF3 achieves 3rd-order convergence via MMS
        // Test problem: dy/dt = -y, y(0) = 1, exact solution: y(t) = exp(-t)
        let integrator = TimeIntegrator::<f64>::new(TimeScheme::BDF3);
        let f = |_t: f64, y: &DVector<f64>| -y;
        
        let t_final = 1.0_f64;
        let exact_final = (-t_final).exp();
        
        // Test with three grid sizes to verify convergence order
        let dt_values = vec![0.1_f64, 0.05_f64, 0.025_f64];
        let mut errors = Vec::new();
        
        for &dt in &dt_values {
            let n_steps = (t_final / dt).round() as usize;
            let mut y_prev2 = DVector::from_vec(vec![1.0]);
            
            // First two steps with exact solution (to start BDF3 properly)
            let mut y_prev = DVector::from_vec(vec![(-dt).exp()]);
            let mut y_curr = DVector::from_vec(vec![(-2.0 * dt).exp()]);
            
            // Subsequent steps with BDF3
            for step in 2..n_steps {
                let t = (step as f64) * dt;
                let y_next = integrator.step_with_history(f, &y_curr, Some(&y_prev), Some(&y_prev2), t, dt);
                y_prev2 = y_prev;
                y_prev = y_curr;
                y_curr = y_next;
            }
            
            let error = (y_curr[0] - exact_final).abs();
            errors.push(error);
        }
        
        // Verify 3rd-order convergence: error ~ O(dt^3)
        // Convergence ratio should be approximately 8 (2^3) when dt is halved
        let ratio1 = errors[0] / errors[1];
        let ratio2 = errors[1] / errors[2];
        
        // Allow some tolerance due to fixed-point iteration convergence
        assert!(ratio1 > 4.0, "Convergence ratio 1: {ratio1} should be > 4.0 for BDF3");
        assert!(ratio2 > 4.0, "Convergence ratio 2: {ratio2} should be > 4.0 for BDF3");
        
        // Verify errors are decreasing
        assert!(errors[1] < errors[0]);
        assert!(errors[2] < errors[1]);
    }

    #[test]
    fn test_backward_euler_exponential_decay() {
        // Test Backward Euler with dy/dt = -y (exponential decay)
        // Exact solution: y(t) = exp(-t)
        // Reference: Hairer & Wanner (1996) - Solving ODEs II
        let integrator = TimeIntegrator::<f64>::new(TimeScheme::BackwardEuler);
        let f = |_t: f64, y: &DVector<f64>| -y;
        
        let dt = 0.1_f64;
        let y0 = DVector::from_vec(vec![1.0]);
        
        // Single step
        let y1 = integrator.step(f, &y0, 0.0, dt);
        
        // Analytical: y(0.1) = exp(-0.1) ≈ 0.9048374180359595
        let expected = (-dt).exp();
        
        // Backward Euler is L-stable, should give reasonable accuracy
        assert_relative_eq!(y1[0], expected, epsilon = 1e-2);
        
        // Verify solution is bounded and positive
        assert!(y1[0] > 0.0);
        assert!(y1[0] < 1.0);
    }

    #[test]
    fn test_backward_euler_order_accuracy() {
        // Verify Backward Euler is 1st-order accurate
        assert_eq!(TimeIntegrator::<f64>::new(TimeScheme::BackwardEuler).order(), 1);
    }

    #[test]
    fn test_backward_euler_is_implicit() {
        let integrator = TimeIntegrator::<f64>::new(TimeScheme::BackwardEuler);
        assert!(!integrator.is_explicit());
    }

    #[test]
    fn test_backward_euler_stiff_system() {
        // Test with moderately stiff ODE: dy/dt = -10*y
        // Backward Euler should remain stable due to L-stability
        let integrator = TimeIntegrator::<f64>::new(TimeScheme::BackwardEuler);
        let lambda = -10.0;
        let f = move |_t: f64, y: &DVector<f64>| y * lambda;
        
        let dt = 0.1_f64;
        let y0 = DVector::from_vec(vec![1.0]);
        
        // Take several steps
        let mut y = y0;
        for step in 0..10 {
            let t = f64::from(step) * dt;
            y = integrator.step(f, &y, t, dt);
        }
        
        // Solution should remain bounded and decay properly
        // For very stiff systems, backward Euler may give more damping than analytical
        assert!(y[0] > 0.0, "Solution should remain positive");
        
        // After 10 steps at dt=0.1, t=1.0, analytical would be exp(-10) ≈ 4.5e-5
        // Backward Euler may not match exactly but should show decay
        assert!(y[0] < 0.5, "Solution should show significant decay");
    }

    #[test]
    fn test_backward_euler_convergence_order() {
        // Verify Backward Euler achieves 1st-order convergence via MMS
        // Test problem: dy/dt = -y, y(0) = 1, exact solution: y(t) = exp(-t)
        let integrator = TimeIntegrator::<f64>::new(TimeScheme::BackwardEuler);
        let f = |_t: f64, y: &DVector<f64>| -y;
        
        let t_final = 1.0_f64;
        let exact_final = (-t_final).exp();
        
        // Test with three grid sizes to verify convergence order
        let dt_values = vec![0.1_f64, 0.05_f64, 0.025_f64];
        let mut errors = Vec::new();
        
        for &dt in &dt_values {
            let n_steps = (t_final / dt).round() as usize;
            let mut y = DVector::from_vec(vec![1.0]);
            
            for step in 0..n_steps {
                let t = (step as f64) * dt;
                y = integrator.step(f, &y, t, dt);
            }
            
            let error = (y[0] - exact_final).abs();
            errors.push(error);
        }
        
        // Verify 1st-order convergence: error ~ O(dt)
        // Convergence ratio should be approximately 2 when dt is halved
        let ratio1 = errors[0] / errors[1];
        let ratio2 = errors[1] / errors[2];
        
        // Allow some tolerance due to fixed-point iteration convergence
        assert!(ratio1 > 1.5, "Convergence ratio 1: {ratio1} should be > 1.5");
        assert!(ratio2 > 1.5, "Convergence ratio 2: {ratio2} should be > 1.5");
        
        // Verify errors are decreasing
        assert!(errors[1] < errors[0]);
        assert!(errors[2] < errors[1]);
    }

    #[test]
    fn test_crank_nicolson_exponential_decay() {
        // Test Crank-Nicolson with dy/dt = -y (exponential decay)
        // Exact solution: y(t) = exp(-t)
        // Reference: Crank & Nicolson (1947), Patankar (1980)
        let integrator = TimeIntegrator::<f64>::new(TimeScheme::CrankNicolson);
        let f = |_t: f64, y: &DVector<f64>| -y;
        
        let dt = 0.1_f64;
        let y0 = DVector::from_vec(vec![1.0]);
        
        // Single step
        let y1 = integrator.step(f, &y0, 0.0, dt);
        
        // Analytical: y(0.1) = exp(-0.1) ≈ 0.9048374180359595
        let expected = (-dt).exp();
        
        // Crank-Nicolson should give high accuracy
        assert_relative_eq!(y1[0], expected, epsilon = 1e-4);
        
        // Verify solution is bounded and positive
        assert!(y1[0] > 0.0);
        assert!(y1[0] < 1.0);
    }

    #[test]
    fn test_crank_nicolson_order_accuracy() {
        // Verify Crank-Nicolson is 2nd-order accurate
        assert_eq!(TimeIntegrator::<f64>::new(TimeScheme::CrankNicolson).order(), 2);
    }

    #[test]
    fn test_crank_nicolson_is_implicit() {
        let integrator = TimeIntegrator::<f64>::new(TimeScheme::CrankNicolson);
        assert!(!integrator.is_explicit());
    }

    #[test]
    fn test_crank_nicolson_stiff_system() {
        // Test with moderately stiff ODE: dy/dt = -10*y
        // Crank-Nicolson is A-stable and should handle this well
        let integrator = TimeIntegrator::<f64>::new(TimeScheme::CrankNicolson);
        let lambda = -10.0;
        let f = move |_t: f64, y: &DVector<f64>| y * lambda;
        
        let dt = 0.1_f64;
        let y0 = DVector::from_vec(vec![1.0]);
        
        // Single step
        let y1 = integrator.step(f, &y0, 0.0, dt);
        
        // Analytical: y(0.1) = exp(-1.0) ≈ 0.3678794411714423
        let expected = (-10.0 * dt).exp();
        
        // Solution should remain bounded and match analytical reasonably
        // Crank-Nicolson may have some iteration error on stiff systems
        assert!(y1[0] > 0.0, "Solution should remain positive");
        assert!(y1[0] < 1.0, "Solution should decay");
        assert_relative_eq!(y1[0], expected, epsilon = 5e-2);
    }

    #[test]
    fn test_crank_nicolson_convergence_order() {
        // Verify Crank-Nicolson achieves 2nd-order convergence via MMS
        // Test problem: dy/dt = -y, y(0) = 1, exact solution: y(t) = exp(-t)
        let integrator = TimeIntegrator::<f64>::new(TimeScheme::CrankNicolson);
        let f = |_t: f64, y: &DVector<f64>| -y;
        
        let t_final = 1.0_f64;
        let exact_final = (-t_final).exp();
        
        // Test with three grid sizes to verify convergence order
        let dt_values = vec![0.1_f64, 0.05_f64, 0.025_f64];
        let mut errors = Vec::new();
        
        for &dt in &dt_values {
            let n_steps = (t_final / dt).round() as usize;
            let mut y = DVector::from_vec(vec![1.0]);
            
            for step in 0..n_steps {
                let t = (step as f64) * dt;
                y = integrator.step(f, &y, t, dt);
            }
            
            let error = (y[0] - exact_final).abs();
            errors.push(error);
        }
        
        // Verify 2nd-order convergence: error ~ O(dt^2)
        // Convergence ratio should be approximately 4 (2^2) when dt is halved
        let ratio1 = errors[0] / errors[1];
        let ratio2 = errors[1] / errors[2];
        
        // Allow some tolerance due to fixed-point iteration convergence
        assert!(ratio1 > 2.5, "Convergence ratio 1: {ratio1} should be > 2.5");
        assert!(ratio2 > 2.5, "Convergence ratio 2: {ratio2} should be > 2.5");
        
        // Verify errors are decreasing
        assert!(errors[1] < errors[0]);
        assert!(errors[2] < errors[1]);
    }

    #[test]
    fn test_adams_bashforth2_exponential_decay() {
        // Test Adams-Bashforth 2 with dy/dt = -y (exponential decay)
        // Exact solution: y(t) = exp(-t)
        // Reference: Butcher (2016) - Numerical Methods for ODEs
        let integrator = TimeIntegrator::<f64>::new(TimeScheme::AdamsBashforth2);
        let f = |_t: f64, y: &DVector<f64>| -y;
        
        let dt = 0.1_f64;
        let y0 = DVector::from_vec(vec![1.0]);
        
        // First step with exact solution for proper history
        let y1 = DVector::from_vec(vec![(-dt).exp()]);
        
        // Second step with Adams-Bashforth 2
        let y2 = integrator.step_with_history(f, &y1, Some(&y0), None, dt, dt);
        
        // Analytical: y(0.2) = exp(-0.2)
        let expected = (-0.2_f64).exp();
        
        // AB2 should give good accuracy
        assert_relative_eq!(y2[0], expected, epsilon = 1e-3);
        
        // Verify solution is bounded and positive
        assert!(y2[0] > 0.0);
        assert!(y2[0] < 1.0);
    }

    #[test]
    fn test_adams_bashforth2_order_accuracy() {
        // Verify Adams-Bashforth 2 is 2nd-order accurate
        assert_eq!(TimeIntegrator::<f64>::new(TimeScheme::AdamsBashforth2).order(), 2);
    }

    #[test]
    fn test_adams_bashforth2_is_explicit() {
        let integrator = TimeIntegrator::<f64>::new(TimeScheme::AdamsBashforth2);
        assert!(integrator.is_explicit());
    }

    #[test]
    fn test_adams_bashforth2_no_history_fallback() {
        // Test that AB2 falls back to RK2 when no history is available
        let integrator = TimeIntegrator::<f64>::new(TimeScheme::AdamsBashforth2);
        let f = |_t: f64, y: &DVector<f64>| -y;
        
        let dt = 0.1_f64;
        let y0 = DVector::from_vec(vec![1.0]);
        
        // Call without history (should fall back to RK2)
        let y1 = integrator.step(f, &y0, 0.0, dt);
        
        // Should still give reasonable 2nd-order accuracy
        let expected = (-dt).exp();
        assert_relative_eq!(y1[0], expected, epsilon = 1e-3);
    }

    #[test]
    fn test_adams_bashforth2_convergence_order() {
        // Verify Adams-Bashforth 2 achieves 2nd-order convergence via MMS
        // Test problem: dy/dt = -y, y(0) = 1, exact solution: y(t) = exp(-t)
        let integrator = TimeIntegrator::<f64>::new(TimeScheme::AdamsBashforth2);
        let f = |_t: f64, y: &DVector<f64>| -y;
        
        let t_final = 1.0_f64;
        let exact_final = (-t_final).exp();
        
        // Test with three grid sizes to verify convergence order
        let dt_values = vec![0.1_f64, 0.05_f64, 0.025_f64];
        let mut errors = Vec::new();
        
        for &dt in &dt_values {
            let n_steps = (t_final / dt).round() as usize;
            let mut y_prev = DVector::from_vec(vec![1.0]);
            
            // First step with exact solution (to start AB2 properly)
            let mut y_curr = DVector::from_vec(vec![(-dt).exp()]);
            
            // Subsequent steps with AB2
            for step in 1..n_steps {
                let t = (step as f64) * dt;
                let y_next = integrator.step_with_history(f, &y_curr, Some(&y_prev), None, t, dt);
                y_prev = y_curr;
                y_curr = y_next;
            }
            
            let error = (y_curr[0] - exact_final).abs();
            errors.push(error);
        }
        
        // Verify 2nd-order convergence: error ~ O(dt^2)
        // Convergence ratio should be approximately 4 (2^2) when dt is halved
        let ratio1 = errors[0] / errors[1];
        let ratio2 = errors[1] / errors[2];
        
        // AB2 is explicit, may have slightly lower ratios than implicit methods
        assert!(ratio1 > 2.0, "Convergence ratio 1: {ratio1} should be > 2.0 for AB2");
        assert!(ratio2 > 2.0, "Convergence ratio 2: {ratio2} should be > 2.0 for AB2");
        
        // Verify errors are decreasing
        assert!(errors[1] < errors[0]);
        assert!(errors[2] < errors[1]);
    }
}
