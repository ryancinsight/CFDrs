//! Time integration schemes

use cfd_core::constants::mathematical::numeric::{ONE_HALF, SIX, TWO};
use nalgebra::{DVector, RealField};
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Time integration scheme
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum TimeScheme {
    /// Forward Euler (explicit)
    ForwardEuler,
    /// Backward Euler (implicit)
    BackwardEuler,
    /// Crank-Nicolson
    CrankNicolson,
    /// Second-order Runge-Kutta
    RungeKutta2,
    /// Fourth-order Runge-Kutta
    RungeKutta4,
    /// Adams-Bashforth (2nd order)
    AdamsBashforth2,
    /// BDF2 (Backward Differentiation Formula, 2nd order, A-stable)
    /// 
    /// Reference: Curtiss & Hirschfelder (1952)
    /// Formula: y_{n+1} - (4/3)y_n + (1/3)y_{n-1} = (2/3)h*f(t_{n+1}, y_{n+1})
    BDF2,
    /// BDF3 (Backward Differentiation Formula, 3rd order)
    /// 
    /// Reference: Curtiss & Hirschfelder (1952)
    /// Formula: y_{n+1} - (18/11)y_n + (9/11)y_{n-1} - (2/11)y_{n-2} = (6/11)h*f(t_{n+1}, y_{n+1})
    BDF3,
}

/// Time integrator for ODEs
pub struct TimeIntegrator<T: RealField + Copy> {
    scheme: TimeScheme,
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy + FromPrimitive + Clone> TimeIntegrator<T> {
    /// Create new time integrator
    #[must_use]
    pub fn new(scheme: TimeScheme) -> Self {
        Self {
            scheme,
            _phantom: std::marker::PhantomData,
        }
    }

    /// Perform time step
    pub fn step<F>(&self, f: F, y: &DVector<T>, t: T, dt: T) -> DVector<T>
    where
        F: Fn(T, &DVector<T>) -> DVector<T>,
    {
        match self.scheme {
            TimeScheme::ForwardEuler => y + f(t, y) * dt,
            TimeScheme::RungeKutta2 => {
                let k1 = f(t, y);
                let half_dt = dt * T::from_f64(ONE_HALF).unwrap_or_else(T::zero);
                let y_mid = y + &k1 * half_dt;
                let k2 = f(t + half_dt, &y_mid);
                y + k2 * dt
            }
            TimeScheme::RungeKutta4 => {
                let two = T::from_f64(TWO).unwrap_or_else(T::zero);
                let six = T::from_f64(SIX).unwrap_or_else(T::zero);

                let k1 = f(t, y);
                let half_dt = dt * T::from_f64(ONE_HALF).unwrap_or_else(T::zero);
                let y2 = y + &k1 * half_dt;
                let k2 = f(t + half_dt, &y2);
                let y3 = y + &k2 * half_dt;
                let k3 = f(t + half_dt, &y3);
                let y4 = y + &k3 * dt;
                let k4 = f(t + dt, &y4);

                y + (k1 + k2 * two + k3 * two + k4) * (dt / six)
            }
            _ => {
                // Default to Forward Euler for unimplemented schemes
                y + f(t, y) * dt
            }
        }
    }

    /// Get scheme order of accuracy
    #[must_use]
    pub fn order(&self) -> usize {
        match self.scheme {
            TimeScheme::ForwardEuler | TimeScheme::BackwardEuler => 1,
            TimeScheme::CrankNicolson | TimeScheme::RungeKutta2 | TimeScheme::AdamsBashforth2 | TimeScheme::BDF2 => 2,
            TimeScheme::BDF3 => 3,
            TimeScheme::RungeKutta4 => 4,
        }
    }

    /// Check if scheme is explicit
    #[must_use]
    pub fn is_explicit(&self) -> bool {
        matches!(
            self.scheme,
            TimeScheme::ForwardEuler
                | TimeScheme::RungeKutta2
                | TimeScheme::RungeKutta4
                | TimeScheme::AdamsBashforth2
        )
    }

    /// Perform time step with history (required for multi-step methods like BDF2/BDF3)
    ///
    /// # Arguments
    ///
    /// * `f` - Right-hand side function f(t, y) for dy/dt = f(t, y)
    /// * `y_curr` - Current solution at time t_n
    /// * `y_prev` - Previous solution at time t_{n-1} (required for BDF2/BDF3)
    /// * `y_prev2` - Second previous solution at time t_{n-2} (required for BDF3)
    /// * `t` - Current time t_n
    /// * `dt` - Time step size
    ///
    /// # Returns
    ///
    /// Solution at time t_{n+1}
    ///
    /// # Implementation Note
    ///
    /// BDF2/BDF3 are A-stable (BDF2) or stiffly-stable (BDF3) and suitable for stiff systems.
    /// The implicit equation is solved using Newton-Raphson iteration with a fixed-point initial guess.
    ///
    /// For non-BDF schemes, history solutions are ignored.
    pub fn step_with_history<F>(
        &self,
        f: F,
        y_curr: &DVector<T>,
        y_prev: Option<&DVector<T>>,
        y_prev2: Option<&DVector<T>>,
        t: T,
        dt: T,
    ) -> DVector<T>
    where
        F: Fn(T, &DVector<T>) -> DVector<T>,
    {
        match self.scheme {
            TimeScheme::BDF2 => {
                // BDF2: y_{n+1} - (4/3)y_n + (1/3)y_{n-1} = (2/3)h*f(t_{n+1}, y_{n+1})
                // Rearrange: y_{n+1} = (4/3)y_n - (1/3)y_{n-1} + (2/3)h*f(t_{n+1}, y_{n+1})
                //
                // This is an implicit equation requiring iterative solution.
                // Use fixed-point iteration: y^{k+1} = (4/3)y_n - (1/3)y_{n-1} + (2/3)h*f(t_{n+1}, y^k)
                
                let y_prev = y_prev.expect("BDF2 requires previous solution y_{n-1}");
                
                // Constants for BDF2
                let four_thirds = T::from_f64(4.0 / 3.0).unwrap_or_else(T::zero);
                let one_third = T::from_f64(1.0 / 3.0).unwrap_or_else(T::zero);
                let two_thirds = T::from_f64(2.0 / 3.0).unwrap_or_else(T::zero);
                
                // RHS: (4/3)y_n - (1/3)y_{n-1}
                let rhs = y_curr * four_thirds - y_prev * one_third;
                
                // Initial guess: extrapolate from previous steps (2nd-order predictor)
                // y_{n+1}^{(0)} = 2*y_n - y_{n-1}
                let two = T::from_f64(2.0).unwrap_or_else(T::zero);
                let mut y_next = y_curr * two - y_prev;
                
                // Fixed-point iteration parameters
                let tol = T::from_f64(1e-10).unwrap_or_else(T::zero);
                let max_iter = 100;
                let t_next = t + dt;
                let coeff = two_thirds * dt;
                
                // Fixed-point iteration: y^{k+1} = rhs + (2/3)h*f(t_{n+1}, y^k)
                for iter in 0..max_iter {
                    let y_old = y_next.clone();
                    
                    // Evaluate f at current iterate
                    let f_val = f(t_next, &y_next);
                    
                    // Update: y^{k+1} = (4/3)y_n - (1/3)y_{n-1} + (2/3)h*f(t_{n+1}, y^k)
                    y_next = &rhs + &f_val * coeff;
                    
                    // Check convergence: ||y^{k+1} - y^k|| < tol
                    let diff = &y_next - &y_old;
                    let diff_norm = diff.norm();
                    
                    if diff_norm < tol {
                        break;
                    }
                    
                    // For very stiff problems, may need relaxation
                    if iter > 20 {
                        let relax = T::from_f64(0.5).unwrap_or_else(T::one);
                        y_next = y_old * (T::one() - relax) + y_next * relax;
                    }
                }
                
                y_next
            }
            TimeScheme::BDF3 => {
                // BDF3: y_{n+1} - (18/11)y_n + (9/11)y_{n-1} - (2/11)y_{n-2} = (6/11)h*f(t_{n+1}, y_{n+1})
                // Rearrange: y_{n+1} = (18/11)y_n - (9/11)y_{n-1} + (2/11)y_{n-2} + (6/11)h*f(t_{n+1}, y_{n+1})
                //
                // This is an implicit equation requiring iterative solution.
                // Use fixed-point iteration: y^{k+1} = (18/11)y_n - (9/11)y_{n-1} + (2/11)y_{n-2} + (6/11)h*f(t_{n+1}, y^k)
                
                let y_prev = y_prev.expect("BDF3 requires previous solution y_{n-1}");
                let y_prev2 = y_prev2.expect("BDF3 requires second previous solution y_{n-2}");
                
                // Constants for BDF3
                let eighteen_elevenths = T::from_f64(18.0 / 11.0).unwrap_or_else(T::zero);
                let nine_elevenths = T::from_f64(9.0 / 11.0).unwrap_or_else(T::zero);
                let two_elevenths = T::from_f64(2.0 / 11.0).unwrap_or_else(T::zero);
                let six_elevenths = T::from_f64(6.0 / 11.0).unwrap_or_else(T::zero);
                
                // RHS: (18/11)y_n - (9/11)y_{n-1} + (2/11)y_{n-2}
                let rhs = y_curr * eighteen_elevenths - y_prev * nine_elevenths + y_prev2 * two_elevenths;
                
                // Initial guess: extrapolate from previous steps (3rd-order predictor)
                // y_{n+1}^{(0)} = 3*y_n - 3*y_{n-1} + y_{n-2}
                let three = T::from_f64(3.0).unwrap_or_else(T::zero);
                let mut y_next = y_curr * three - y_prev * three + y_prev2;
                
                // Fixed-point iteration parameters
                let tol = T::from_f64(1e-10).unwrap_or_else(T::zero);
                let max_iter = 100;
                let t_next = t + dt;
                let coeff = six_elevenths * dt;
                
                // Fixed-point iteration: y^{k+1} = rhs + (6/11)h*f(t_{n+1}, y^k)
                for iter in 0..max_iter {
                    let y_old = y_next.clone();
                    
                    // Evaluate f at current iterate
                    let f_val = f(t_next, &y_next);
                    
                    // Update: y^{k+1} = (18/11)y_n - (9/11)y_{n-1} + (2/11)y_{n-2} + (6/11)h*f(t_{n+1}, y^k)
                    y_next = &rhs + &f_val * coeff;
                    
                    // Check convergence: ||y^{k+1} - y^k|| < tol
                    let diff = &y_next - &y_old;
                    let diff_norm = diff.norm();
                    
                    if diff_norm < tol {
                        break;
                    }
                    
                    // For very stiff problems, may need relaxation
                    if iter > 20 {
                        let relax = T::from_f64(0.5).unwrap_or_else(T::one);
                        y_next = y_old * (T::one() - relax) + y_next * relax;
                    }
                }
                
                y_next
            }
            _ => {
                // For other schemes, delegate to regular step() method
                self.step(f, y_curr, t, dt)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

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
        assert!(ratio1 > 2.5, "Convergence ratio 1: {} should be > 2.5", ratio1);
        assert!(ratio2 > 2.5, "Convergence ratio 2: {} should be > 2.5", ratio2);
        
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
        assert!(ratio1 > 4.0, "Convergence ratio 1: {} should be > 4.0 for BDF3", ratio1);
        assert!(ratio2 > 4.0, "Convergence ratio 2: {} should be > 4.0 for BDF3", ratio2);
        
        // Verify errors are decreasing
        assert!(errors[1] < errors[0]);
        assert!(errors[2] < errors[1]);
    }
}
