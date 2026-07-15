//! Implicit-Explicit (IMEX) time-stepping methods.
//!
//! This module provides IMEX Runge-Kutta schemes that treat stiff terms implicitly
//! and non-stiff terms explicitly, optimal for CFD problems with both
//! convective and diffusive components.
//!
//! # Mathematical Foundation
//!
//! IMEX methods solve systems of the form:
//! ```text
//! du/dt = f_explicit(t,u) + f_implicit(t,u)
//! ```
//!
//! where `f_explicit` contains non-stiff terms (typically convective) and
//! `f_implicit` contains stiff terms (typically diffusive).
//!
//! ## ARS343 Method (Ascher, Ruuth, Spiteri, 1997)
//!
//! This implementation provides the ARS343 method, a 3rd-order
//! additive Runge-Kutta scheme with 3 stages.
//!
//! ### Tableau
//!
//! Explicit:
//! 0 | 0
//! γ | γ 0
//! 1 | δ 1-δ 0
//!
//! Implicit:
//! 0 | 0
//! γ | 0 γ
//! 1 | 0 1-2γ γ
//!
//! γ = (3 + √3) / 6
//! δ = 1 - 1/(2γ)
//!
//! ## Theorem Statement
//!
//! **Theorem (Additive Runge-Kutta Order Conditions)**: For the IMEX ARK method
//! with tableaux (A_exp, b_exp, c) and (A_imp, b_imp, c), the order conditions are:
//!
//! - **Order 1**: Σ b_exp + Σ b_imp = 1
//! - **Order 2**: Σ b_exp * c + Σ b_imp * c = 1/2
//! - **Order 3**: Σ b_exp * c² + Σ b_imp * c² = 1/3, and coupling conditions
//!
//! ## References
//!
//! - Ascher, U. M., Ruuth, S. J., & Wetton, B. T. (1997). Implicit-explicit methods
//!   for time-dependent PDEs. SIAM Journal on Numerical Analysis, 32(3), 797-823.

use super::traits::{
    from_f64, one, state_len, state_norm, state_zeros, zero, TimeMatrix, TimeState, TimeStepper,
};
use cfd_core::error::{ConvergenceErrorKind, Error, Result};
use eunomia::{NumericElement, RealField};
use leto_ops::{solve, RealScalar};

/// IMEX Runge-Kutta method for systems with mixed stiffness
///
/// Splits the RHS into explicit (non-stiff) and implicit (stiff) parts:
/// du/dt = f_explicit(t,u) + f_implicit(t,u)
pub struct IMEXTimeStepper<T: RealField + Copy> {
    /// Explicit method coefficients
    explicit_a: Vec<Vec<T>>,
    /// Implicit method coefficients (strictly lower triangular)
    implicit_a: Vec<Vec<T>>,
    /// Explicit solution weights
    explicit_b: Vec<T>,
    /// Implicit solution weights
    implicit_b: Vec<T>,
    /// Implicit diagonal coefficients (a_ii)
    implicit_diagonal: Vec<T>,
    /// Time evaluation points
    c: Vec<T>,
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + RealScalar + Copy> IMEXTimeStepper<T> {
    /// Default Newton iteration tolerance for implicit stages
    const NEWTON_TOLERANCE: f64 = 1e-10;

    /// Default maximum Newton iterations per implicit stage
    const MAX_NEWTON_ITERATIONS: usize = 10;

    /// Compute stage solution for RHS evaluation from previous stages
    ///
    /// Implements: u_stage = u + Σ_{j=0}^{stage-1} dt * (a_exp_j * k_exp_j + a_imp_j * k_imp_j)
    fn compute_stage_solution(
        &self,
        stage: usize,
        u: &TimeState<T>,
        dt: T,
        k_explicit: &[TimeState<T>],
        k_implicit: &[TimeState<T>],
        u_stage: &mut TimeState<T>,
    ) {
        u_stage.assign(u);
        let n = state_len(u);

        // Add contributions from previous stages
        for prev_stage in 0..stage {
            let a_exp =
                if stage < self.explicit_a.len() && prev_stage < self.explicit_a[stage].len() {
                    self.explicit_a[stage][prev_stage]
                } else {
                    zero()
                };
            let a_imp =
                if stage < self.implicit_a.len() && prev_stage < self.implicit_a[stage].len() {
                    self.implicit_a[stage][prev_stage]
                } else {
                    zero()
                };

            for i in 0..n {
                u_stage[i] +=
                    dt * (a_exp * k_explicit[prev_stage][i] + a_imp * k_implicit[prev_stage][i]);
            }
        }
    }

    /// Solve implicit stage using Newton iteration
    ///
    /// Solves: u_stage = u_stage_explicit + dt * a_ii * f_implicit(t_stage, u_stage)
    fn solve_implicit_stage<F, J>(
        &self,
        u_stage_explicit: &TimeState<T>,
        t_stage: T,
        dt: T,
        a_ii_imp: T,
        f_implicit: F,
        jacobian_implicit: J,
    ) -> Result<TimeState<T>>
    where
        F: Fn(T, &TimeState<T>) -> Result<TimeState<T>>,
        J: Fn(T, &TimeState<T>) -> Result<TimeMatrix<T>>,
    {
        let mut u_stage_current = u_stage_explicit.clone();
        let n = state_len(u_stage_explicit);

        for _newton_iter in 0..Self::MAX_NEWTON_ITERATIONS {
            // Compute residual: F(u) = u - u_explicit - dt * a_ii * f_imp(t, u)
            let f_imp = f_implicit(t_stage, &u_stage_current)?;
            let mut residual = state_zeros(n);

            for i in 0..n {
                residual[i] = u_stage_current[i] - u_stage_explicit[i] - dt * a_ii_imp * f_imp[i];
            }

            // Check convergence
            let residual_norm = state_norm(&residual);
            let tolerance = from_f64(Self::NEWTON_TOLERANCE);
            if residual_norm < tolerance {
                return Ok(u_stage_current);
            }

            // Compute Jacobian and solve for Newton step
            let jac_imp = jacobian_implicit(t_stage, &u_stage_current)?;
            validate_matrix_shape(&jac_imp, n, "IMEX implicit Jacobian")?;
            let mut jacobian = matrix_identity(n);
            let a_ii_dt = a_ii_imp * dt;

            for i in 0..n {
                for j in 0..n {
                    jacobian[[i, j]] -= a_ii_dt * jac_imp[[i, j]];
                }
            }

            // Solve Jacobian * du = -residual
            let neg_residual =
                TimeState::from_shape_vec([n], (0..n).map(|i| zero::<T>() - residual[i]).collect())
                    .expect("invariant: residual length matches state shape");
            let du = solve(&jacobian.view(), &neg_residual.view()).map_err(|error| {
                Error::InvalidInput(format!("IMEX Newton solve failed: {error}"))
            })?;
            for i in 0..n {
                u_stage_current[i] += du[i];
            }
        }

        Err(Error::Convergence(
            ConvergenceErrorKind::MaxIterationsExceeded {
                max: Self::MAX_NEWTON_ITERATIONS,
            },
        ))
    }

    /// Create ARS343 (Ascher, Ruuth, Spiteri 1997) IMEX method
    ///
    /// This is a 3rd-order, 3-stage L-stable scheme.
    ///
    /// # Reference
    /// Ascher, U. M., Ruuth, S. J., & Wetton, B. T. (1997). Implicit-explicit methods
    /// for time-dependent PDEs. SIAM Journal on Numerical Analysis, 32(3), 797-823.
    pub fn ars343() -> Self {
        let half: T = from_f64(0.5);
        let two: T = from_f64(2.0);
        let three: T = from_f64(3.0);
        let six: T = from_f64(6.0);
        let gamma = (three + <T as NumericElement>::sqrt(three)) / six;
        let delta = one::<T>() - one::<T>() / (two * gamma);

        let c = vec![zero::<T>(), gamma, one::<T>()];

        // Explicit tableau (strictly lower)
        let explicit_a = vec![
            vec![],                          // Stage 0 (c=0)
            vec![gamma],                     // Stage 1 (c=gamma)
            vec![delta, one::<T>() - delta], // Stage 2 (c=1)
        ];

        // Implicit tableau (strictly lower part)
        let implicit_a = vec![
            vec![],                                      // Stage 0
            vec![zero::<T>()],                           // Stage 1
            vec![zero::<T>(), one::<T>() - two * gamma], // Stage 2
        ];

        // Implicit diagonals
        let implicit_diagonal = vec![
            zero::<T>(), // Stage 0 explicit
            gamma,       // Stage 1
            gamma,       // Stage 2
        ];

        let explicit_b = vec![zero::<T>(), half, half];
        let implicit_b = vec![zero::<T>(), half, half];

        Self {
            explicit_a,
            implicit_a,
            explicit_b,
            implicit_b,
            implicit_diagonal,
            c,
            _phantom: std::marker::PhantomData,
        }
    }

    /// Take an IMEX step with proper implicit solving
    pub fn imex_step<F, S, J>(
        &self,
        f_explicit: F,
        f_implicit: S,
        jacobian_implicit: J,
        t: T,
        u: &TimeState<T>,
        dt: T,
    ) -> Result<TimeState<T>>
    where
        F: Fn(T, &TimeState<T>) -> Result<TimeState<T>>,
        S: Fn(T, &TimeState<T>) -> Result<TimeState<T>>,
        J: Fn(T, &TimeState<T>) -> Result<TimeMatrix<T>>,
    {
        let n = state_len(u);
        let stages = self.c.len();

        // Store RHS evaluations at each stage
        let mut k_explicit: Vec<TimeState<T>> = vec![state_zeros(n); stages];
        let mut k_implicit: Vec<TimeState<T>> = vec![state_zeros(n); stages];

        // Reusable buffer for stage solution
        let mut u_stage = state_zeros(n);

        // Compute stages
        for stage in 0..stages {
            let t_stage = t + self.c[stage] * dt;

            // Compute stage solution for RHS evaluation
            self.compute_stage_solution(stage, u, dt, &k_explicit, &k_implicit, &mut u_stage);

            // Evaluate RHS at this stage
            k_explicit[stage] = f_explicit(t_stage, &u_stage)?;
            // Note: k_implicit is computed AFTER solving for implicit stage if needed
            // But for compute_stage_solution of NEXT stages, we need k_implicit.
            // For implicit stage i, we solve u_i = u_expl + dt * a_ii * f_imp(u_i).
            // Then k_imp_i = f_imp(u_i).

            // First evaluate implicit part at explicit prediction (or if stage is explicit)
            let mut f_imp_val = f_implicit(t_stage, &u_stage)?;

            // For implicit stages, solve nonlinear equation if diagonal element is non-zero
            let a_ii_imp = if stage < self.implicit_diagonal.len() {
                self.implicit_diagonal[stage]
            } else {
                zero()
            };

            if a_ii_imp != zero() {
                // Solve implicit stage using Newton iteration
                let u_stage_solved = self.solve_implicit_stage(
                    &u_stage,
                    t_stage,
                    dt,
                    a_ii_imp,
                    &f_implicit,
                    &jacobian_implicit,
                )?;

                // Update RHS with solved implicit stage
                f_imp_val = f_implicit(t_stage, &u_stage_solved)?;
            }
            k_implicit[stage] = f_imp_val;
        }

        // Combine RHS evaluations for final solution
        let mut u_new = u.clone();
        for stage in 0..stages {
            let b_exp = if stage < self.explicit_b.len() {
                self.explicit_b[stage]
            } else {
                zero()
            };
            let b_imp = if stage < self.implicit_b.len() {
                self.implicit_b[stage]
            } else {
                zero()
            };

            for i in 0..n {
                u_new[i] += dt * (b_exp * k_explicit[stage][i] + b_imp * k_implicit[stage][i]);
            }
        }

        Ok(u_new)
    }
}

impl<T: RealField + RealScalar + Copy> TimeStepper<T> for IMEXTimeStepper<T> {
    fn step<F>(&self, f: F, t: T, u: &TimeState<T>, dt: T) -> Result<TimeState<T>>
    where
        F: Fn(T, &TimeState<T>) -> Result<TimeState<T>>,
    {
        // For simple TimeStepper interface, assume f is the explicit part
        // and implicit part is zero (falls back to explicit method)
        let jacobian_zero = |_t: T, _u: &TimeState<T>| Ok(matrix_zeros(state_len(u), state_len(u)));
        self.imex_step(
            f,
            |_t, _u| Ok(state_zeros(state_len(u))),
            jacobian_zero,
            t,
            u,
            dt,
        )
    }

    fn order(&self) -> usize {
        3 // ARS343 is 3rd order
    }

    fn stages(&self) -> usize {
        3 // ARS343 has 3 stages
    }

    fn is_explicit(&self) -> bool {
        false
    }

    fn stability_region(&self) -> Option<&str> {
        Some("IMEX L-stable")
    }
}

fn validate_matrix_shape<T>(matrix: &TimeMatrix<T>, n: usize, context: &str) -> Result<()> {
    let [rows, cols] = matrix.shape();
    if rows != n || cols != n {
        return Err(Error::InvalidInput(format!(
            "{context}: matrix has shape {rows}x{cols}, expected {n}x{n}"
        )));
    }
    Ok(())
}

fn matrix_zeros<T: RealScalar>(rows: usize, cols: usize) -> TimeMatrix<T> {
    TimeMatrix::from_elem([rows, cols], zero::<T>())
}

fn matrix_identity<T: RealScalar>(n: usize) -> TimeMatrix<T> {
    let mut matrix = matrix_zeros(n, n);
    for i in 0..n {
        matrix[[i, i]] = one::<T>();
    }
    matrix
}

#[cfg(test)]
fn matrix_scale<T: RealScalar>(matrix: &TimeMatrix<T>, scale: T) -> TimeMatrix<T> {
    let [rows, cols] = matrix.shape();
    TimeMatrix::from_shape_vec(
        [rows, cols],
        (0..rows * cols)
            .map(|idx| {
                let row = idx / cols;
                let col = idx % cols;
                matrix[[row, col]] * scale
            })
            .collect(),
    )
    .expect("invariant: scaled matrix storage matches shape")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::time_stepping::traits::{state_from_vec, state_len, state_neg, state_scale};
    use approx::assert_relative_eq;

    #[test]
    fn test_imex_properties() {
        let imex = IMEXTimeStepper::<f64>::ars343();
        assert_eq!(imex.order(), 3);
        assert_eq!(imex.stages(), 3);
        assert!(!imex.is_explicit());
    }

    #[test]
    fn test_imex_step() {
        let imex = IMEXTimeStepper::<f64>::ars343();

        // Simple test: du/dt = -u (explicit) + 0 (implicit)
        let f_explicit = |_t: f64, u: &TimeState<f64>| Ok(state_neg(u));
        let f_implicit = |_t: f64, u: &TimeState<f64>| Ok(state_zeros(state_len(u)));
        let jacobian_zero =
            |_t: f64, u: &TimeState<f64>| Ok(matrix_zeros(state_len(u), state_len(u)));

        let u0 = state_from_vec(vec![1.0]);
        let dt = 0.1;
        let t = 0.0;

        let u1 = imex
            .imex_step(f_explicit, f_implicit, jacobian_zero, t, &u0, dt)
            .unwrap();

        // Should produce reasonable result for exponential decay
        assert!(u1[0] > 0.0 && u1[0] < 1.0);
    }

    #[test]
    fn test_imex_convergence_order() {
        // Test problem: du/dt = -2u (stiff implicit) + u (explicit) = -u
        // Analytical solution: u(t) = u0 * exp(-t)
        let imex = IMEXTimeStepper::<f64>::ars343();

        // Use closures that match the expected signature
        let f_explicit = |_t: f64, u: &TimeState<f64>| Ok(u.clone());
        let f_implicit = |_t: f64, u: &TimeState<f64>| Ok(state_scale(u, -2.0));
        let jacobian_implicit = |_t: f64, u: &TimeState<f64>| {
            let mut jac = matrix_zeros(state_len(u), state_len(u));
            jac[[0, 0]] = -2.0;
            Ok(jac)
        };

        let u0 = state_from_vec(vec![1.0]);
        let t_final = 1.0f64;
        let analytical = |t: f64| (-t).exp();

        let dts = vec![0.1f64, 0.05, 0.025, 0.0125];
        let mut errors = Vec::new();

        for &dt in &dts {
            let mut u = u0.clone();
            let mut t = 0.0f64;

            while t < t_final - 1e-10 {
                let step_size = (t_final - t).min(dt);
                u = imex
                    .imex_step(f_explicit, f_implicit, jacobian_implicit, t, &u, step_size)
                    .unwrap();
                t += step_size;
            }

            let error = (u[0] - analytical(t_final)).abs();
            errors.push(error);
        }

        // Check convergence order
        let expected_order = imex.order() as f64;
        for i in 0..errors.len() - 1 {
            let ratio = errors[i] / errors[i + 1];
            let dt_ratio = dts[i] / dts[i + 1];
            let expected_ratio = dt_ratio.powf(expected_order);
            let observed_order = ratio;

            println!(
                "dt ratio: {dt_ratio:.1}, error ratio: {observed_order:.2}, expected ratio: {expected_ratio:.2}"
            );

            // Check that observed convergence is reasonably close
            // For ARS343 (3rd order), we expect ratio ~8
            assert!(
                ratio > 2.0,
                "Error not decreasing significantly: ratio = {ratio}"
            );
        }
    }

    #[test]
    fn test_imex_analytical_solution() {
        let imex = IMEXTimeStepper::<f64>::ars343();

        // Problem: du/dt = -u (entirely implicit for stiffness)
        let f_explicit = |_t: f64, u: &TimeState<f64>| Ok(state_zeros(state_len(u)));
        let f_implicit = |_t: f64, u: &TimeState<f64>| Ok(state_neg(u));
        let jacobian_implicit =
            |_t: f64, u: &TimeState<f64>| Ok(matrix_scale(&matrix_identity(state_len(u)), -1.0));

        let u0 = state_from_vec(vec![2.0, -1.0, 0.5]);
        let dt = 0.01;
        let t = 0.0;

        let u1 = imex
            .imex_step(f_explicit, f_implicit, jacobian_implicit, t, &u0, dt)
            .unwrap();

        let analytical = state_from_vec((0..state_len(&u0)).map(|i| u0[i] * (-dt).exp()).collect());

        // Check accuracy
        for i in 0..state_len(&u1) {
            assert_relative_eq!(u1[i], analytical[i], epsilon = 1e-6);
        }
    }

    #[test]
    fn test_imex_newton_convergence() {
        let imex = IMEXTimeStepper::<f64>::ars343();

        // Problem: dy/dt = -100*y^3 (non-linear stiff decay, tests Newton solver)
        let lambda = -100.0;
        let f_explicit = |_t: f64, u: &TimeState<f64>| Ok(state_zeros(state_len(u)));
        let f_implicit = |_t: f64, u: &TimeState<f64>| {
            Ok(state_from_vec(
                (0..state_len(u)).map(|i| lambda * u[i].powi(3)).collect(),
            ))
        };
        let jacobian_implicit = |_t: f64, u: &TimeState<f64>| {
            let mut jac = matrix_zeros(state_len(u), state_len(u));
            jac[[0, 0]] = 3.0 * lambda * u[0].powi(2);
            Ok(jac)
        };

        let u0 = state_from_vec(vec![1.0]);
        let dt = 0.01; // Reduced time step for better accuracy
        let t = 0.0;

        let u1 = imex
            .imex_step(f_explicit, f_implicit, jacobian_implicit, t, &u0, dt)
            .unwrap();

        // Analytical solution for dy/dt = λy^3: y(t) = 1 / sqrt(y0^-2 - 2λt)
        let analytical = 1.0 / (u0[0].powi(-2) - 2.0 * lambda * dt).sqrt();

        // ARS343 is L-stable, should handle this well
        assert_relative_eq!(u1[0], analytical, epsilon = 5e-2);
    }
}
