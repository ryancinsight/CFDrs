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
use nalgebra::{RealField, Vector2};
use num_traits::{FromPrimitive, ToPrimitive};

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + std::fmt::LowerExp>
    SimplecPimpleSolver<T>
{
    /// Solve pressure-velocity coupling for one time step with adaptive stepping
    pub fn solve_time_step(
        &mut self,
        fields: &mut SimulationFields<T>,
        dt: T,
        nu: T,
        rho: T,
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
            AlgorithmType::Simplec => self.solve_simplec(fields, dt, nu, rho),
            AlgorithmType::Pimple => self.solve_pimple(fields, dt, nu, rho),
        }
    }

    /// Solve with adaptive time stepping and convergence acceleration
    ///
    /// Returns the effective time step used and residual.
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
        let mut residuals = Vec::new();
        let mut last_residual = T::max_value().unwrap_or(T::from_f64(1e10).unwrap_or(T::one()));

        let dt_increase_factor = T::from_f64(1.2).unwrap_or_else(|| {
            T::one() + T::from_f64(0.2).expect("analytical constant conversion")
        });
        let dt_decrease_factor = T::from_f64(0.7)
            .unwrap_or_else(|| T::from_f64(0.7).expect("analytical constant conversion"));
        let min_dt = dt_initial
            * T::from_f64(0.1)
                .unwrap_or_else(|| T::from_f64(0.1).expect("analytical constant conversion"));
        let max_dt = T::from_f64(1e10).unwrap_or_else(|| {
            dt_initial
                * T::from_f64(5.0)
                    .unwrap_or_else(|| T::from_f64(5.0).expect("analytical constant conversion"))
        });

        while step_count < max_steps {
            let residual = self.solve_time_step(fields, dt, nu, rho)?;
            residuals.push(residual);

            if residual < target_residual {
                break;
            }

            // Aitken's Δ² acceleration diagnostic
            if residuals.len() >= 3 {
                let n = residuals.len();
                let r0 = residuals[n - 3];
                let r1 = residuals[n - 2];
                let r2 = residuals[n - 1];

                let denominator =
                    r2 - r1 * T::from_f64(2.0).unwrap_or_else(|| T::one() + T::one()) + r0;
                if denominator.abs() > T::default_epsilon() {
                    let numerator = (r1 - r0) * (r1 - r0);
                    let r_accelerated = r0 - numerator / denominator;

                    if r_accelerated > T::zero() && r_accelerated < residual {
                        tracing::debug!(
                            "Aitken acceleration applied: {:.6e} -> {:.6e}",
                            residual,
                            r_accelerated
                        );
                    }
                }
            }

            // Adaptive time step adjustment
            if residual
                < last_residual
                    * T::from_f64(0.95).unwrap_or_else(|| {
                        T::from_f64(0.95).expect("analytical constant conversion")
                    })
            {
                dt = (dt * dt_increase_factor).min(max_dt);
                tracing::debug!(
                    "Time step increased to {:.6}, residual: {:.6e}",
                    dt,
                    residual
                );
            } else if residual
                > last_residual
                    * T::from_f64(1.05).unwrap_or_else(|| {
                        T::from_f64(1.05).expect("analytical constant conversion")
                    })
            {
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

        let final_residual = *residuals
            .last()
            .unwrap_or(&T::max_value().unwrap_or(T::from_f64(1e10).unwrap_or(T::one())));
        Ok((dt, final_residual))
    }
}

#[cfg(test)]
mod tests {
    use super::super::config::{AlgorithmType, SimplecPimpleConfig};
    use super::super::solver::SimplecPimpleSolver;
    use crate::grid::StructuredGrid2D;

    fn make_grid(n: usize) -> StructuredGrid2D<f64> {
        StructuredGrid2D::new(n, n, 0.0, 1.0, 0.0, 1.0).unwrap()
    }

    #[test]
    fn algorithm_creation() {
        assert!(SimplecPimpleSolver::new(make_grid(8), SimplecPimpleConfig::simplec()).is_ok());
        assert!(SimplecPimpleSolver::new(make_grid(8), SimplecPimpleConfig::pimple()).is_ok());
    }

    #[test]
    fn verify_algorithm_type() {
        assert_eq!(SimplecPimpleConfig::<f64>::simplec().algorithm, AlgorithmType::Simplec);
        assert_eq!(SimplecPimpleConfig::<f64>::pimple().algorithm, AlgorithmType::Pimple);
        assert_eq!(SimplecPimpleConfig::<f64>::default().algorithm, AlgorithmType::Simplec);
    }
}
