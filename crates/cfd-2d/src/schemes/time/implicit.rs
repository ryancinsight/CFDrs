//! Implicit time integration methods

use cfd_core::constants::mathematical::numeric::{ONE_HALF, TWO};
use nalgebra::{DVector, RealField};
use num_traits::FromPrimitive;

/// Backward Euler step: y_{n+1} = y_n + dt*f(t_{n+1}, y_{n+1})
///
/// Implicit equation solved via fixed-point iteration.
///
/// Reference: Hairer & Wanner (1996) - Solving Ordinary Differential Equations II
pub fn backward_euler<T, F>(f: F, y: &DVector<T>, t: T, dt: T) -> DVector<T>
where
    T: RealField + Copy + FromPrimitive,
    F: Fn(T, &DVector<T>) -> DVector<T>,
{
    let t_next = t + dt;
    let tol = T::from_f64(1e-10).unwrap_or_else(T::zero);
    let max_iter = 100;
    
    // Initial guess: Forward Euler predictor
    let mut y_next = y + &f(t, y) * dt;
    
    // Fixed-point iteration: y^{k+1} = y_n + dt*f(t_{n+1}, y^k)
    for iter in 0..max_iter {
        let y_old = y_next.clone();
        
        // Evaluate f at current iterate
        let f_val = f(t_next, &y_next);
        
        // Update: y^{k+1} = y_n + dt*f(t_{n+1}, y^k)
        y_next = y + &f_val * dt;
        
        // Check convergence: ||y^{k+1} - y^k|| < tol
        let diff = &y_next - &y_old;
        let diff_norm = diff.norm();
        
        if diff_norm < tol {
            break;
        }
        
        // For very stiff problems, apply relaxation
        if iter > 20 {
            let relax = T::from_f64(0.5).unwrap_or_else(T::one);
            y_next = y_old * (T::one() - relax) + y_next * relax;
        }
    }
    
    y_next
}

/// Crank-Nicolson step (θ-method with θ=0.5):
/// y_{n+1} = y_n + (dt/2)*(f(t_n, y_n) + f(t_{n+1}, y_{n+1}))
///
/// Implicit equation solved via fixed-point iteration.
///
/// Reference: Crank & Nicolson (1947), Patankar (1980)
pub fn crank_nicolson<T, F>(f: F, y: &DVector<T>, t: T, dt: T) -> DVector<T>
where
    T: RealField + Copy + FromPrimitive,
    F: Fn(T, &DVector<T>) -> DVector<T>,
{
    let half = T::from_f64(ONE_HALF).unwrap_or_else(T::zero);
    let t_next = t + dt;
    let tol = T::from_f64(1e-10).unwrap_or_else(T::zero);
    let max_iter = 100;
    
    // Explicit part: dt/2 * f(t_n, y_n)
    let explicit_part = f(t, y) * (dt * half);
    
    // Initial guess: Forward Euler predictor
    let mut y_next = y + &explicit_part * T::from_f64(TWO).unwrap_or_else(T::one);
    
    // Fixed-point iteration: y^{k+1} = y_n + (dt/2)*(f(t_n, y_n) + f(t_{n+1}, y^k))
    for iter in 0..max_iter {
        let y_old = y_next.clone();
        
        // Implicit part: dt/2 * f(t_{n+1}, y^k)
        let implicit_part = f(t_next, &y_next) * (dt * half);
        
        // Update: y^{k+1} = y_n + (dt/2)*(f_n + f_{n+1})
        y_next = y + &explicit_part + &implicit_part;
        
        // Check convergence: ||y^{k+1} - y^k|| < tol
        let diff = &y_next - &y_old;
        let diff_norm = diff.norm();
        
        if diff_norm < tol {
            break;
        }
        
        // For very stiff problems, apply relaxation
        if iter > 20 {
            let relax = T::from_f64(0.5).unwrap_or_else(T::one);
            y_next = y_old * (T::one() - relax) + y_next * relax;
        }
    }
    
    y_next
}
