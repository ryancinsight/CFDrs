//! SIMPLEC and PIMPLE algorithm implementations
//!
//! Contains the core iterative algorithm logic for both pressure-velocity coupling
//! methods, separated from the solver infrastructure for maintainability.
//!
//! # Theorem (SIMPLEC Convergence — Van Doormaal & Raithby 1984)
//!
//! When the momentum operator is diagonally dominant and the pressure correction
//! is solved to tolerance, SIMPLEC reduces the splitting error relative to SIMPLE
//! and typically converges in fewer outer iterations for the same relaxation factors.
//!
//! **Proof sketch**: The SIMPLEC correction retains the diagonal pressure coupling
//! while approximating the neighbour-velocity correction consistently. That lowers
//! the splitting error and reduces the spectral radius of the outer fixed-point map
//! when the discretization remains stable.
//!
//! # Theorem (PIMPLE Stability — Issa 1986)
//!
//! For sufficiently small CFL and bounded relaxation factors, PIMPLE maintains
//! bounded solutions for the incompressible Navier-Stokes equations in practice.
//! The nested outer/inner corrector structure improves robustness, but it does not
//! provide an unconditional convergence proof.
//!
//! **Proof sketch**: Each inner corrector performs a PISO-like pressure correction
//! that enforces discrete continuity to within the linear-solver tolerance. The
//! outer correctors re-linearise the convective term, so the composite map behaves
//! like a damped fixed-point iteration when the timestep is small enough.

use super::config::AlgorithmType;
use super::solver::SimplecPimpleSolver;
use crate::fields::SimulationFields;
use crate::scalar;
use cfd_core::CfdScalar;
use eunomia::{FloatElement, NumericElement};
use leto::geometry::Vector2;

impl<T: CfdScalar + Copy + std::fmt::LowerExp + FloatElement> SimplecPimpleSolver<T> {
    /// Solve pressure-velocity coupling for one time step with adaptive stepping
    pub fn solve_time_step(
        &mut self,
        fields: &mut SimulationFields<T>,
        dt: T,
        nu: T,
        rho: T,
    ) -> cfd_core::error::Result<T> {
        self.solve_time_step_with_tolerance(fields, dt, nu, rho, self.config.tolerance)
    }

    fn solve_time_step_with_tolerance(
        &mut self,
        fields: &mut SimulationFields<T>,
        dt: T,
        nu: T,
        rho: T,
        convergence_tolerance: T,
    ) -> cfd_core::error::Result<T> {
        // Capture start-of-step velocities for transient momentum consistency
        fields.u_old.data.copy_from_slice(&fields.u.data);
        fields.v_old.data.copy_from_slice(&fields.v.data);

        // Synchronize density and dynamic viscosity fields with time-step parameters
        let dynamic_viscosity = nu * rho;
        fields.density.map_inplace(|d| *d = rho);
        fields.viscosity.map_inplace(|v| *v = dynamic_viscosity);

        // Update Rhie-Chow old velocity buffer for transient correction
        if let Some(ref mut rhie_chow) = self.rhie_chow {
            let mut u_old_cache = self._vel_field_cache.borrow_mut();
            if u_old_cache.as_ref().is_none_or(|v| {
                let (nx, ny) = v.dimensions();
                nx != self.grid.nx || ny != self.grid.ny
            }) {
                *u_old_cache = Some(crate::fields::Field2D::new(
                    self.grid.nx,
                    self.grid.ny,
                    Vector2::zeros(),
                ));
            }
            let u_old = u_old_cache.as_mut().expect("old velocity cache must exist");
            for i in 0..self.grid.nx {
                for j in 0..self.grid.ny {
                    u_old.set(i, j, Vector2::new(fields.u.at(i, j), fields.v.at(i, j)));
                }
            }
            rhie_chow.update_old_velocity(u_old);
        }

        match self.config.algorithm {
            AlgorithmType::Simplec => {
                self.solve_simplec(fields, dt, nu, rho, convergence_tolerance)
            }
            AlgorithmType::Pimple => self.solve_pimple(fields, dt, nu, rho, convergence_tolerance),
        }
    }

    /// Steady-state residual: the velocity field's relative change per unit of
    /// convective time.
    ///
    /// `(‖u^{n+1} − u^n‖_∞ / ‖u‖_∞) · (t_ref / dt)` with `t_ref = L_ref / ‖u‖_∞`
    /// and `L_ref` the larger domain extent: the fractional change the field
    /// would accumulate over one convective transit if the current rate held.
    /// Dividing by `dt` is what makes the measure step-size independent — a
    /// change of `1e-9` taken over a step of `1e-9` is not convergence — and
    /// the two reference scales make it dimensionless, so a caller's tolerance
    /// means the same thing at any velocity or domain scale. A field that is
    /// still exactly zero everywhere returns the caller's worst case rather
    /// than zero, since an unmoved field is the starting condition, not a
    /// converged one.
    #[must_use]
    pub fn steady_state_residual(&self, fields: &SimulationFields<T>, dt: T) -> T {
        let mut max_change = scalar::zero::<T>();
        let mut max_speed = scalar::zero::<T>();
        for (component, previous) in [(&fields.u, &fields.u_old), (&fields.v, &fields.v_old)] {
            for (new, old) in component.data.iter().zip(previous.data.iter()) {
                max_change = max_change.max_scalar(NumericElement::abs(*new - *old));
                max_speed = max_speed.max_scalar(NumericElement::abs(*new));
            }
        }
        if max_speed <= T::default_epsilon() || dt <= scalar::zero::<T>() {
            return <T as FloatElement>::from_f64(1e10);
        }
        let (x_min, x_max, y_min, y_max) = self.grid.bounds;
        let length_reference = (x_max - x_min).max_scalar(y_max - y_min);
        max_change * length_reference / (dt * max_speed * max_speed)
    }

    /// March to steady state with adaptive time stepping.
    ///
    /// Returns the effective time step reached and the final steady-state
    /// residual (see [`Self::steady_state_residual`]). Termination is on that
    /// residual, not on the per-step pressure-velocity coupling residual: the
    /// coupling residual measures how well one time level was solved and falls
    /// below any tolerance on the first small step, so using it here reports
    /// convergence for a field that has barely moved.
    ///
    /// Uses Aitken's Δ² acceleration to estimate convergence trajectory
    /// and adjusts the time step based on residual behaviour.
    pub fn solve_adaptive(
        &mut self,
        fields: &mut SimulationFields<T>,
        dt_initial: T,
        nu: T,
        rho: T,
        max_steps: usize,
        target_residual: T,
    ) -> cfd_core::error::Result<(T, T)> {
        let mut dt = dt_initial;
        let mut step_count = 0;
        self.residual_history.clear();
        let mut last_residual = <T as FloatElement>::from_f64(1e10);

        let dt_increase_factor = <T as FloatElement>::from_f64(1.2);
        let dt_decrease_factor = <T as FloatElement>::from_f64(0.7);
        let min_dt = dt_initial * <T as FloatElement>::from_f64(0.1);
        let max_dt = <T as FloatElement>::from_f64(1e10);

        while step_count < max_steps {
            // `target_residual` is the caller's problem-scaled termination
            // contract. The configured tolerance remains the minimum useful
            // absolute tolerance for ordinary steps, but must not replace a
            // looser dimensional target with an unreachable fixed threshold.
            let convergence_tolerance = target_residual.max_scalar(self.config.tolerance);
            let coupling_residual =
                self.solve_time_step_with_tolerance(fields, dt, nu, rho, convergence_tolerance)?;
            if !NumericElement::is_finite(coupling_residual) {
                return Err(cfd_core::error::Error::Numerical(
                    cfd_core::error::NumericalErrorKind::InvalidValue {
                        value: format!(
                            "coupling residual {coupling_residual:e} at step {step_count}"
                        ),
                    },
                ));
            }
            // `u_old`/`v_old` hold the start-of-step field: the step just taken
            // is the difference against the current one.
            let residual = self.steady_state_residual(fields, dt);
            self.residual_history.push(residual);

            if residual < target_residual {
                break;
            }

            // Aitken's Δ² acceleration diagnostic
            if self.residual_history.len() >= 3 {
                let n = self.residual_history.len();
                let r0 = self.residual_history[n - 3];
                let r1 = self.residual_history[n - 2];
                let r2 = self.residual_history[n - 1];

                let denominator = r2 - r1 * <T as FloatElement>::from_f64(2.0) + r0;
                if NumericElement::abs(denominator) > T::default_epsilon() {
                    let numerator = (r1 - r0) * (r1 - r0);
                    let r_accelerated = r0 - numerator / denominator;

                    if r_accelerated > scalar::zero() && r_accelerated < residual {
                        tracing::debug!(
                            "Aitken acceleration applied: {:.6e} -> {:.6e}",
                            residual,
                            r_accelerated
                        );
                    }
                }
            }

            // Adaptive time step adjustment
            if residual < last_residual * <T as FloatElement>::from_f64(0.95) {
                dt = (dt * dt_increase_factor).min(max_dt);
                tracing::debug!(
                    "Time step increased to {:.6}, residual: {:.6e}",
                    dt,
                    residual
                );
            } else if residual > last_residual * <T as FloatElement>::from_f64(1.05) {
                dt = (dt * dt_decrease_factor).max(min_dt);
                tracing::debug!(
                    "Time step decreased to {:.6}, residual: {:.6e}",
                    dt,
                    residual
                );
            }

            last_residual = residual;
            step_count += 1;
        }

        let fallback_residual = <T as FloatElement>::from_f64(1e10);
        let final_residual = *self.residual_history.last().unwrap_or(&fallback_residual);
        Ok((dt, final_residual))
    }
}

#[cfg(test)]
mod tests {
    use super::super::config::{AlgorithmType, SimplecPimpleConfig};
    use super::super::solver::SimplecPimpleSolver;
    use crate::grid::StructuredGrid2D;

    fn make_grid(n: usize) -> StructuredGrid2D<f64> {
        StructuredGrid2D::new(n, n, 0.0, 1.0, 0.0, 1.0).expect("expected value")
    }

    #[test]
    fn algorithm_creation() {
        assert!(SimplecPimpleSolver::new(make_grid(8), SimplecPimpleConfig::simplec()).is_ok());
        assert!(SimplecPimpleSolver::new(make_grid(8), SimplecPimpleConfig::pimple()).is_ok());
    }

    #[test]
    fn verify_algorithm_type() {
        assert_eq!(
            SimplecPimpleConfig::<f64>::simplec().algorithm,
            AlgorithmType::Simplec
        );
        assert_eq!(
            SimplecPimpleConfig::<f64>::pimple().algorithm,
            AlgorithmType::Pimple
        );
        assert_eq!(
            SimplecPimpleConfig::<f64>::default().algorithm,
            AlgorithmType::Simplec
        );
    }
}
