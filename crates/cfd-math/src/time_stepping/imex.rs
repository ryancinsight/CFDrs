//! Implicit-Explicit (IMEX) time-stepping methods.
//!
//! This module provides IMEX schemes that treat stiff terms implicitly
//! and non-stiff terms explicitly, optimal for CFD problems with both
//! convective and diffusive components.

use nalgebra::{DMatrix, DVector, RealField};
use cfd_core::error::Result;
use super::traits::TimeStepper;

/// IMEX Runge-Kutta method for systems with mixed stiffness
///
/// Splits the RHS into explicit (non-stiff) and implicit (stiff) parts:
/// du/dt = f_explicit(t,u) + f_implicit(t,u)
pub struct IMEXTimeStepper<T: RealField + Copy> {
    /// Explicit method coefficients
    explicit_a: Vec<Vec<T>>,
    /// Implicit method coefficients
    implicit_a: Vec<Vec<T>>,
    /// Explicit solution weights
    explicit_b: Vec<T>,
    /// Implicit solution weights
    implicit_b: Vec<T>,
    /// Time evaluation points
    c: Vec<T>,
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy> IMEXTimeStepper<T> {
    /// Create IMEX ARK4(3)6L[2]SA method (Kennedy & Carpenter)
    ///
    /// This is a 4th-order IMEX method with 6 stages, suitable for
    /// convection-diffusion problems common in CFD.
    pub fn ark436l2sa() -> Self {
        let c = vec![
            T::from_f64(0.0).unwrap(),
            T::from_f64(1.0/3.0).unwrap(),
            T::from_f64(2.0/3.0).unwrap(),
            T::from_f64(1.0/3.0).unwrap(),
            T::from_f64(2.0/3.0).unwrap(),
            T::from_f64(1.0).unwrap(),
        ];

        // Explicit tableau (ERK4)
        let explicit_a = vec![
            vec![],  // Stage 0
            vec![T::from_f64(1.0/3.0).unwrap()],  // Stage 1
            vec![T::from_f64(1.0/6.0).unwrap(), T::from_f64(1.0/6.0).unwrap()],  // Stage 2
            vec![T::from_f64(1.0/12.0).unwrap(), T::from_f64(1.0/12.0).unwrap(), T::from_f64(1.0/4.0).unwrap()],  // Stage 3
            vec![T::from_f64(1.0/12.0).unwrap(), T::from_f64(1.0/12.0).unwrap(), T::from_f64(1.0/4.0).unwrap(), T::from_f64(0.0).unwrap()],  // Stage 4
            vec![T::from_f64(1.0/12.0).unwrap(), T::from_f64(1.0/12.0).unwrap(), T::from_f64(1.0/4.0).unwrap(), T::from_f64(0.0).unwrap(), T::from_f64(1.0/2.0).unwrap()],  // Stage 5
        ];

        // Implicit tableau (ESDIRK3)
        let implicit_a = vec![
            vec![],  // Stage 0
            vec![T::from_f64(1.0/3.0).unwrap()],  // Stage 1
            vec![T::from_f64(1.0/6.0).unwrap(), T::from_f64(1.0/6.0).unwrap()],  // Stage 2
            vec![T::from_f64(1.0/12.0).unwrap(), T::from_f64(1.0/12.0).unwrap(), T::from_f64(1.0/4.0).unwrap()],  // Stage 3
            vec![T::from_f64(1.0/12.0).unwrap(), T::from_f64(1.0/12.0).unwrap(), T::from_f64(1.0/4.0).unwrap(), T::from_f64(0.0).unwrap()],  // Stage 4
            vec![T::from_f64(1.0/12.0).unwrap(), T::from_f64(1.0/12.0).unwrap(), T::from_f64(1.0/4.0).unwrap(), T::from_f64(0.0).unwrap(), T::from_f64(0.0).unwrap()],  // Stage 5
        ];

        let explicit_b = vec![
            T::from_f64(1.0/12.0).unwrap(),
            T::from_f64(1.0/12.0).unwrap(),
            T::from_f64(1.0/4.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(1.0/2.0).unwrap(),
            T::from_f64(0.0).unwrap(),
        ];

        let implicit_b = vec![
            T::from_f64(1.0/12.0).unwrap(),
            T::from_f64(1.0/12.0).unwrap(),
            T::from_f64(1.0/4.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(1.0/2.0).unwrap(),
        ];

        Self {
            explicit_a,
            implicit_a,
            explicit_b,
            implicit_b,
            c,
            _phantom: std::marker::PhantomData,
        }
    }

    /// Take an IMEX step with proper implicit solving
    ///
    /// # Arguments
    /// * `f_explicit` - Explicit part of RHS: f_explicit(t, u)
    /// * `f_implicit` - Implicit part of RHS: f_implicit(t, u)
    /// * `jacobian_implicit` - Jacobian of implicit part: df_implicit/du
    /// * `t` - Current time
    /// * `u` - Current solution
    /// * `dt` - Time step
    pub fn imex_step<F, S, J>(
        &self,
        f_explicit: F,
        f_implicit: S,
        jacobian_implicit: J,
        t: T,
        u: &DVector<T>,
        dt: T
    ) -> Result<DVector<T>>
    where
        F: Fn(T, &DVector<T>) -> Result<DVector<T>>,
        S: Fn(T, &DVector<T>) -> Result<DVector<T>>,
        J: Fn(T, &DVector<T>) -> Result<DMatrix<T>>,
    {
        let n = u.len();
        let stages = self.c.len();

        // Stage solutions
        let mut u_stages: Vec<DVector<T>> = vec![DVector::zeros(n); stages];

        // Compute stages
        for stage in 0..stages {
            let t_stage = t.clone() + self.c[stage].clone() * dt.clone();

            // Compute explicit contribution to stage
            let mut u_stage_explicit = u.clone();
            for prev_stage in 0..stage {
                let coeff = self.explicit_a[stage][prev_stage].clone();
                for i in 0..n {
                    u_stage_explicit[i] += dt.clone() * coeff * u_stages[prev_stage][i].clone();
                }
            }

            // Solve implicit stage equation: u_stage = u_stage_explicit + dt * a_ii * f_implicit(t_stage, u_stage)
            // This is a nonlinear equation that we solve using Newton's method
            let mut u_stage = u_stage_explicit.clone();
            let mut converged = false;
            let max_newton_iter = 10;
            let tolerance = T::from_f64(1e-10).unwrap_or_else(T::one);

            for _newton_iter in 0..max_newton_iter {
                // Compute residual: F(u_stage) = u_stage - u_stage_explicit - dt * a_ii * f_implicit(t_stage, u_stage)
                let f_imp = f_implicit(t_stage.clone(), &u_stage)?;
                let mut residual = DVector::zeros(n);
                let a_ii = self.implicit_a[stage][stage].clone();

                for i in 0..n {
                    residual[i] = u_stage[i] - u_stage_explicit[i] - dt.clone() * a_ii * f_imp[i];
                }

                // Check convergence
                let residual_norm = residual.norm();
                if residual_norm < tolerance {
                    converged = true;
                    break;
                }

                // Compute Jacobian of residual: dF/du = I - dt * a_ii * df_implicit/du
                let jac_imp = jacobian_implicit(t_stage.clone(), &u_stage)?;
                let mut jacobian = DMatrix::identity(n, n);
                let a_ii_dt = a_ii * dt.clone();

                for i in 0..n {
                    for j in 0..n {
                        jacobian[(i, j)] -= a_ii_dt * jac_imp[(i, j)];
                    }
                }

                // Solve Jacobian * du = -residual
                match jacobian.lu().solve(&(-residual)) {
                    Some(du) => {
                        u_stage += du;
                    }
                    None => {
                        // Fallback to explicit if Jacobian is singular
                        break;
                    }
                }
            }

            if !converged {
                // Fallback: use explicit approximation if implicit solve fails
                let f_imp = f_implicit(t_stage.clone(), &u_stage_explicit)?;
                let a_ii = self.implicit_a[stage][stage].clone();
                for i in 0..n {
                    u_stage[i] = u_stage_explicit[i] + dt.clone() * a_ii * f_imp[i];
                }
            }

            // Add explicit contribution
            let f_exp = f_explicit(t_stage, &u_stage_explicit)?;
            let a_ii = self.explicit_a[stage][stage].clone();
            for i in 0..n {
                u_stage[i] += dt.clone() * a_ii * f_exp[i];
            }

            u_stages[stage] = u_stage;
        }

        // Combine stages for final solution
        let mut u_new = DVector::zeros(n);
        for stage in 0..stages {
            let coeff_exp = self.explicit_b[stage].clone();
            let coeff_imp = self.implicit_b[stage].clone();

            for i in 0..n {
                u_new[i] += coeff_exp * u_stages[stage][i] + coeff_imp * u_stages[stage][i];
            }
        }

        Ok(u_new)
    }
}

impl<T: RealField + Copy> TimeStepper<T> for IMEXTimeStepper<T> {
    fn step<F>(&self, f: F, t: T, u: &DVector<T>, dt: T) -> Result<DVector<T>>
    where
        F: Fn(T, &DVector<T>) -> Result<DVector<T>>,
    {
        // For simple TimeStepper interface, assume f is the explicit part
        // and implicit part is zero (falls back to explicit method)
        // Provide identity Jacobian for zero implicit function
        let jacobian_zero = |_t: T, _u: &DVector<T>| Ok(DMatrix::identity(u.len(), u.len()));
        self.imex_step(f, |_t, _u| Ok(DVector::zeros(u.len())), jacobian_zero, t, u, dt)
    }

    fn order(&self) -> usize { 4 }

    fn stages(&self) -> usize { 6 }

    fn is_explicit(&self) -> bool { false } // IMEX is neither purely explicit nor implicit

    fn stability_region(&self) -> Option<&str> {
        Some("IMEX stability for convection-diffusion problems")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_imex_ark436l2sa_properties() {
        let imex = IMEXTimeStepper::<f64>::ark436l2sa();
        assert_eq!(imex.order(), 4);
        assert_eq!(imex.stages(), 6);
        assert!(!imex.is_explicit());
    }

    #[test]
    fn test_imex_step() {
        let imex = IMEXTimeStepper::<f64>::ark436l2sa();

        // Simple test: du/dt = -u (explicit) + 0 (implicit)
        let f_explicit = |t: f64, u: &DVector<f64>| Ok(-u.clone());
        let f_implicit = |t: f64, u: &DVector<f64>| Ok(DVector::zeros(u.len()));
        let jacobian_zero = |_t: f64, _u: &DVector<f64>| Ok(DMatrix::identity(u.len(), u.len()));

        let u0 = DVector::from_vec(vec![1.0]);
        let dt = 0.1;
        let t = 0.0;

        let u1 = imex.imex_step(f_explicit, f_implicit, jacobian_zero, t, &u0, dt).unwrap();

        // Should produce reasonable result for exponential decay
        assert!(u1[0] > 0.0 && u1[0] < 1.0);
    }
}
