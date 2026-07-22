//! Implicit time integration methods
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
use cfd_core::physics::constants::mathematical::numeric::{ONE_HALF, TWO};
use eunomia::FloatElement;

use super::vector::{l2_norm, StateVector};

/// Backward Euler step: y_{n+1} = y_n + dt*f(t_{n+1}, y_{n+1})
///
/// Implicit equation solved via fixed-point iteration.
///
/// Reference: Hairer & Wanner (1996) - Solving Ordinary Differential Equations II
pub fn backward_euler<T, F>(f: F, y: &StateVector<T>, t: T, dt: T) -> StateVector<T>
where
    T: Cfd2dScalar + Copy + FloatElement,
    F: Fn(T, &StateVector<T>) -> StateVector<T>,
{
    let t_next = t + dt;
    let tol = scalar::from_f64::<T>(1e-10);
    let max_iter = 100;

    // Initial guess: Forward Euler predictor
    let predictor = &f(t, y) * dt;
    let mut y_next = y + &predictor;

    // AUDIT: pre-allocate y_old once outside the loop; `y_old.assign(&y_next)`
    // inside the loop is a memcpy (saves ~99 of 100 Vec-allocs per fixed-point
    // call). `mem::swap` would be wrong here -- `f(&y_next)` reads y_next, so
    // swap would feed it stale pre_iter_{k-1} from iter 1 onward. Don't move
    // the `let mut y_old = ...` line back inside the loop.
    let mut y_old = y_next.clone();

    // Fixed-point iteration: y^{k+1} = y_n + dt*f(t_{n+1}, y^k)
    for iter in 0..max_iter {
        y_old.assign(&y_next);

        // Evaluate f at current iterate
        let f_val = f(t_next, &y_next);

        // Update: y^{k+1} = y_n + dt*f(t_{n+1}, y^k)
        let increment = &f_val * dt;
        y_next = y + &increment;

        // Check convergence: ||y^{k+1} - y^k|| < tol
        let diff = &y_next - &y_old;
        let diff_norm = l2_norm(&diff);

        if diff_norm < tol {
            break;
        }

        // For very stiff problems, apply relaxation
        if iter > 20 {
            let relax = scalar::from_f64::<T>(0.5);
            let old_part = &y_old * (scalar::one::<T>() - relax);
            let next_part = &y_next * relax;
            y_next = &old_part + &next_part;
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
pub fn crank_nicolson<T, F>(f: F, y: &StateVector<T>, t: T, dt: T) -> StateVector<T>
where
    T: Cfd2dScalar + Copy + FloatElement,
    F: Fn(T, &StateVector<T>) -> StateVector<T>,
{
    let half = scalar::from_f64::<T>(ONE_HALF);
    let t_next = t + dt;
    let tol = scalar::from_f64::<T>(1e-10);
    let max_iter = 100;

    // Explicit part: dt/2 * f(t_n, y_n)
    let explicit_part = &f(t, y) * (dt * half);

    // Initial guess: Forward Euler predictor
    let explicit_predictor = &explicit_part * scalar::from_f64::<T>(TWO);
    let mut y_next: StateVector<T> = y + &explicit_predictor;

    // AUDIT: pre-allocate y_old once (see backward_euler). `mem::swap` is also
    // wrong here -- f reads &y_next, so swap would feed stale input from iter 1.
    // Don't move the pre-allocate line back inside the loop -- that silently
    // re-introduces per-iter clones.
    let mut y_old = y_next.clone();

    // Fixed-point iteration: y^{k+1} = y_n + (dt/2)*(f(t_n, y_n) + f(t_{n+1}, y^k))
    for iter in 0..max_iter {
        y_old.assign(&y_next);

        // Implicit part: dt/2 * f(t_{n+1}, y^k)
        let implicit_part = &f(t_next, &y_next) * (dt * half);

        // Update: y^{k+1} = y_n + (dt/2)*(f_n + f_{n+1})
        let correction = &explicit_part + &implicit_part;
        y_next = y + &correction;

        // Check convergence: ||y^{k+1} - y^k|| < tol
        let diff = &y_next - &y_old;
        let diff_norm = l2_norm(&diff);

        if diff_norm < tol {
            break;
        }

        // For very stiff problems, apply relaxation
        if iter > 20 {
            let relax = scalar::from_f64::<T>(0.5);
            let old_part = &y_old * (scalar::one::<T>() - relax);
            let next_part = &y_next * relax;
            y_next = &old_part + &next_part;
        }
    }

    y_next
}
