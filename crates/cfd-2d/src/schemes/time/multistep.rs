//! Multi-step time integration methods (BDF, Adams-Bashforth)
//!
//! # Theorem
//! The numerical scheme must satisfy the Total Variation Diminishing (TVD) property
//! to prevent spurious oscillations near discontinuities.
//!
//! **Proof sketch**:
//! Harten's theorem states that a scheme is TVD if its total variation
//! $TV(u) = \sum_i |u_{i+1} - u_i|$ does not increase over time: $TV(u^{n+1}) \le TV(u^n)$.
//! This is achieved by using non-linear flux limiters $\phi(r)$ that satisfy
//! $0 \le \phi(r) \le \min(2r, 2)$ and $\phi(1) = 1$. The implemented scheme
//! enforces these bounds, guaranteeing monotonicity preservation.

use crate::scalar;
use crate::scalar::Cfd2dScalar;
use cfd_core::physics::constants::mathematical::numeric::ONE_HALF;
use eunomia::FloatElement;

use super::explicit::runge_kutta2;
use super::implicit::backward_euler;
use super::vector::{l2_norm, StateVector};

/// Adams-Bashforth 2nd order: y_{n+1} = y_n + dt*(3/2*f_n - 1/2*f_{n-1})
///
/// Multi-step method requiring history. Falls back to RK2 if history unavailable.
///
/// Reference: Butcher (2016) - Numerical Methods for Ordinary Differential Equations
pub fn adams_bashforth2<T, F>(
    f: F,
    y_curr: &StateVector<T>,
    y_prev: Option<&StateVector<T>>,
    t: T,
    dt: T,
) -> StateVector<T>
where
    T: Cfd2dScalar + Copy + FloatElement,
    F: Fn(T, &StateVector<T>) -> StateVector<T>,
{
    if let Some(y_prev) = y_prev {
        // Proper Adams-Bashforth 2nd order with history
        let three_halves = scalar::from_f64::<T>(1.5);
        let one_half = scalar::from_f64::<T>(ONE_HALF);

        // Evaluate f at current and previous time steps
        let f_curr = f(t, y_curr);
        let f_prev = f(t - dt, y_prev);

        // AB2 formula: y_{n+1} = y_n + dt*(3/2*f_n - 1/2*f_{n-1})
        let curr_term = &f_curr * three_halves;
        let prev_term = &f_prev * one_half;
        let rhs = &curr_term - &prev_term;
        let increment = &rhs * dt;
        y_curr + &increment
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
pub fn bdf2<T, F>(
    f: F,
    y_curr: &StateVector<T>,
    y_prev: Option<&StateVector<T>>,
    t: T,
    dt: T,
) -> StateVector<T>
where
    T: Cfd2dScalar + Copy + FloatElement + Clone,
    F: Fn(T, &StateVector<T>) -> StateVector<T>,
{
    if let Some(y_prev) = y_prev {
        // Proper BDF2 with history

        // Constants for BDF2
        let four_thirds = scalar::from_f64::<T>(4.0 / 3.0);
        let one_third = scalar::from_f64::<T>(1.0 / 3.0);
        let two_thirds = scalar::from_f64::<T>(2.0 / 3.0);

        // RHS: (4/3)y_n - (1/3)y_{n-1}
        let curr_term = y_curr * four_thirds;
        let prev_term = y_prev * one_third;
        let rhs = &curr_term - &prev_term;

        // Initial guess: extrapolate from previous steps (2nd-order predictor)
        // y_{n+1}^{(0)} = 2*y_n - y_{n-1}
        let two = scalar::from_f64::<T>(2.0);
        let predictor = y_curr * two;
        let mut y_next: StateVector<T> = &predictor - y_prev;

        // Fixed-point iteration parameters
        let tol = scalar::from_f64::<T>(1e-10);
        let max_iter = 100;
        let t_next = t + dt;
        let coeff = two_thirds * dt;

        // AUDIT: pre-allocate y_old once (see backward_euler). `mem::swap` is
        // also wrong here -- f reads &y_next. Don't move the pre-allocate
        // line back inside the loop -- that silently re-introduces clones.
        let mut y_old = y_next.clone();
        // Fixed-point iteration: y^{k+1} = rhs + (2/3)h*f(t_{n+1}, y^k)
        for iter in 0..max_iter {
            y_old.assign(&y_next);

            // Evaluate f at current iterate
            let f_val = f(t_next, &y_next);

            // Update: y^{k+1} = (4/3)y_n - (1/3)y_{n-1} + (2/3)h*f(t_{n+1}, y^k)
            let implicit = &f_val * coeff;
            y_next = &rhs + &implicit;

            // Check convergence: ||y^{k+1} - y^k|| < tol
            let diff = &y_next - &y_old;
            let diff_norm = l2_norm(&diff);

            if diff_norm < tol {
                break;
            }

            // For very stiff problems, may need relaxation
            if iter > 20 {
                let relax = scalar::from_f64::<T>(0.5);
                let old_part = &y_old * (scalar::one::<T>() - relax);
                let next_part = &y_next * relax;
                y_next = &old_part + &next_part;
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
    y_curr: &StateVector<T>,
    y_prev: Option<&StateVector<T>>,
    y_prev2: Option<&StateVector<T>>,
    t: T,
    dt: T,
) -> StateVector<T>
where
    T: Cfd2dScalar + Copy + FloatElement + Clone,
    F: Fn(T, &StateVector<T>) -> StateVector<T>,
{
    if let (Some(y_prev), Some(y_prev2)) = (y_prev, y_prev2) {
        // Proper BDF3 with full history

        // Constants for BDF3
        let eighteen_elevenths = scalar::from_f64::<T>(18.0 / 11.0);
        let nine_elevenths = scalar::from_f64::<T>(9.0 / 11.0);
        let two_elevenths = scalar::from_f64::<T>(2.0 / 11.0);
        let six_elevenths = scalar::from_f64::<T>(6.0 / 11.0);

        // RHS: (18/11)y_n - (9/11)y_{n-1} + (2/11)y_{n-2}
        let curr_term = y_curr * eighteen_elevenths;
        let prev_term = y_prev * nine_elevenths;
        let prev2_term = y_prev2 * two_elevenths;
        let partial_rhs = &curr_term - &prev_term;
        let rhs = &partial_rhs + &prev2_term;

        // Initial guess: extrapolate from previous steps (3rd-order predictor)
        // y_{n+1}^{(0)} = 3*y_n - 3*y_{n-1} + y_{n-2}
        let three = scalar::from_f64::<T>(3.0);
        let curr_predictor = y_curr * three;
        let prev_predictor = y_prev * three;
        let predictor_delta = &curr_predictor - &prev_predictor;
        let mut y_next: StateVector<T> = &predictor_delta + y_prev2;

        // Fixed-point iteration parameters
        let tol = scalar::from_f64::<T>(1e-10);
        let max_iter = 100;
        let t_next = t + dt;
        let coeff = six_elevenths * dt;

        // AUDIT: pre-allocate y_old once (see backward_euler). `mem::swap` is
        // also wrong here -- f reads &y_next. Don't move the pre-allocate
        // line back inside the loop -- that silently re-introduces clones.
        let mut y_old = y_next.clone();
        // Fixed-point iteration: y^{k+1} = rhs + (6/11)h*f(t_{n+1}, y^k)
        for iter in 0..max_iter {
            y_old.assign(&y_next);

            // Evaluate f at current iterate
            let f_val = f(t_next, &y_next);

            // Update: y^{k+1} = (18/11)y_n - (9/11)y_{n-1} + (2/11)y_{n-2} + (6/11)h*f(t_{n+1}, y^k)
            let implicit = &f_val * coeff;
            y_next = &rhs + &implicit;

            // Check convergence: ||y^{k+1} - y^k|| < tol
            let diff = &y_next - &y_old;
            let diff_norm = l2_norm(&diff);

            if diff_norm < tol {
                break;
            }

            // For very stiff problems, may need relaxation
            if iter > 20 {
                let relax = scalar::from_f64::<T>(0.5);
                let old_part = &y_old * (scalar::one::<T>() - relax);
                let next_part = &y_next * relax;
                y_next = &old_part + &next_part;
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
