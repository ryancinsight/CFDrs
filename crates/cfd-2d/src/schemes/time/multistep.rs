//! Multi-step time integration methods (BDF, Adams-Bashforth)

use cfd_core::constants::mathematical::numeric::ONE_HALF;
use nalgebra::{DVector, RealField};
use num_traits::FromPrimitive;

use super::explicit::runge_kutta2;
use super::implicit::backward_euler;

/// Adams-Bashforth 2nd order: y_{n+1} = y_n + dt*(3/2*f_n - 1/2*f_{n-1})
///
/// Multi-step method requiring history. Falls back to RK2 if history unavailable.
///
/// Reference: Butcher (2016) - Numerical Methods for Ordinary Differential Equations
pub fn adams_bashforth2<T, F>(
    f: F,
    y_curr: &DVector<T>,
    y_prev: Option<&DVector<T>>,
    t: T,
    dt: T,
) -> DVector<T>
where
    T: RealField + Copy + FromPrimitive,
    F: Fn(T, &DVector<T>) -> DVector<T>,
{
    if let Some(y_prev) = y_prev {
        // Proper Adams-Bashforth 2nd order with history
        let three_halves = T::from_f64(1.5).unwrap_or_else(T::one);
        let one_half = T::from_f64(ONE_HALF).unwrap_or_else(T::zero);

        // Evaluate f at current and previous time steps
        let f_curr = f(t, y_curr);
        let f_prev = f(t - dt, y_prev);

        // AB2 formula: y_{n+1} = y_n + dt*(3/2*f_n - 1/2*f_{n-1})
        y_curr + (f_curr * three_halves - f_prev * one_half) * dt
    } else {
        // No history available - fall back to RK2 (also 2nd-order)
        runge_kutta2(f, y_curr, t, dt)
    }
}

/// BDF2: y_{n+1} - (4/3)y_n + (1/3)y_{n-1} = (2/3)h*f(t_{n+1}, y_{n+1})
///
/// Rearranged: y_{n+1} = (4/3)y_n - (1/3)y_{n-1} + (2/3)h*f(t_{n+1}, y_{n+1})
///
/// This is an implicit equation requiring iterative solution.
/// Uses fixed-point iteration: y^{k+1} = (4/3)y_n - (1/3)y_{n-1} + (2/3)h*f(t_{n+1}, y^k)
///
/// Falls back to Backward Euler if history unavailable.
///
/// Reference: Curtiss & Hirschfelder (1952)
pub fn bdf2<T, F>(f: F, y_curr: &DVector<T>, y_prev: Option<&DVector<T>>, t: T, dt: T) -> DVector<T>
where
    T: RealField + Copy + FromPrimitive + Clone,
    F: Fn(T, &DVector<T>) -> DVector<T>,
{
    if let Some(y_prev) = y_prev {
        // Proper BDF2 with history

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
    } else {
        // No history - fall back to Backward Euler (1st-order implicit)
        backward_euler(f, y_curr, t, dt)
    }
}

/// BDF3: y_{n+1} - (18/11)y_n + (9/11)y_{n-1} - (2/11)y_{n-2} = (6/11)h*f(t_{n+1}, y_{n+1})
///
/// Rearranged: y_{n+1} = (18/11)y_n - (9/11)y_{n-1} + (2/11)y_{n-2} + (6/11)h*f(t_{n+1}, y_{n+1})
///
/// This is an implicit equation requiring iterative solution.
/// Uses fixed-point iteration: y^{k+1} = (18/11)y_n - (9/11)y_{n-1} + (2/11)y_{n-2} + (6/11)h*f(t_{n+1}, y^k)
///
/// Falls back to BDF2 if only partial history, or Backward Euler if no history.
///
/// Reference: Curtiss & Hirschfelder (1952)
pub fn bdf3<T, F>(
    f: F,
    y_curr: &DVector<T>,
    y_prev: Option<&DVector<T>>,
    y_prev2: Option<&DVector<T>>,
    t: T,
    dt: T,
) -> DVector<T>
where
    T: RealField + Copy + FromPrimitive + Clone,
    F: Fn(T, &DVector<T>) -> DVector<T>,
{
    if let (Some(y_prev), Some(y_prev2)) = (y_prev, y_prev2) {
        // Proper BDF3 with full history

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
    } else if let Some(y_prev) = y_prev {
        // Partial history - fall back to BDF2
        bdf2(f, y_curr, Some(y_prev), t, dt)
    } else {
        // No history - fall back to Backward Euler (1st-order implicit)
        backward_euler(f, y_curr, t, dt)
    }
}
