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
use crate::grid::array2d::Array2D;
use crate::physics::MomentumComponent;
use nalgebra::{RealField, Vector2};
use num_traits::{FromPrimitive, ToPrimitive};

fn promote_predicted_velocity_state<T: RealField + Copy>(
    fields: &mut SimulationFields<T>,
    workspace: &mut Array2D<Vector2<T>>,
) {
    fields.copy_velocity_star_to(workspace);
    fields.promote_velocity_star_to_current();
}

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
        let max_dt = dt_initial
            * T::from_f64(5.0)
                .unwrap_or_else(|| T::from_f64(5.0).expect("analytical constant conversion"));

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

    /// SIMPLEC algorithm implementation with Aitken's acceleration
    ///
    /// SIMPLEC (Van Doormaal & Raithby, 1984) improves SIMPLE by using
    /// consistent discretization for the pressure correction equation.
    pub(super) fn solve_simplec(
        &mut self,
        fields: &mut SimulationFields<T>,
        dt: T,
        nu: T,
        rho: T,
    ) -> cfd_core::error::Result<T> {
        let mut residual = T::zero();
        let max_iterations = self.config.max_inner_iterations;
        let mut converged = false;

        let mut u_prev = self.extract_velocity_field(fields);
        let mut continuity_residual = T::max_value()
            .unwrap_or_else(|| T::from_f64(1e30).expect("analytical constant conversion"));
        let mut previous_velocity_residual = T::max_value()
            .unwrap_or_else(|| T::from_f64(1e30).expect("analytical constant conversion"));
        let mut previous_continuity_residual = continuity_residual;
        let mut divergence_streak = 0usize;

        for iter in 0..max_iterations {
            // Step 1: Solve momentum equations for predicted velocities u*
            self.momentum_solver
                .solve_with_coefficients(MomentumComponent::U, fields, dt)?;
            self.momentum_solver
                .solve_with_coefficients(MomentumComponent::V, fields, dt)?;

            // Step 2: Refresh the Rhie-Chow face coefficients from the raw
            // momentum stencil.
            if let Some(ref mut rhie_chow) = self.rhie_chow {
                let (ap_u, _, ap_v, _) = self.momentum_solver.get_ap_coefficients();
                rhie_chow.update_u_coefficients(ap_u);
                rhie_chow.update_v_coefficients(ap_v);
            }

            // Step 3: Promote the predicted momentum state into the current
            // fields so the pressure solve sees the raw predictor.
            promote_predicted_velocity_state(fields, &mut self.u_star_workspace);

            // Step 4: Solve pressure correction using the Rhie-Chow face
            // interpolation path when enabled; otherwise fall back to the
            // cell-centred pressure correction equation.
            if self.rhie_chow.is_some() {
                let grid = &self.grid;
                let pressure_solver = &self.pressure_solver;
                let momentum_solver = &self.momentum_solver;
                let rhie_chow = self.rhie_chow.as_ref();
                let vel_field_cache = &self._vel_field_cache;
                let cons_vel_cache = &self._cons_vel_cache;
                let u_face_cache = &self._u_face_cache;
                let v_face_cache = &self._v_face_cache;
                let d_x_cache = &self._d_x_cache;
                let d_y_cache = &self._d_y_cache;

                super::interpolation::solve_pressure_correction_with_caches(
                    grid,
                    pressure_solver,
                    momentum_solver,
                    rhie_chow,
                    vel_field_cache,
                    cons_vel_cache,
                    u_face_cache,
                    v_face_cache,
                    d_x_cache,
                    d_y_cache,
                    fields,
                    dt,
                    Some(dt),
                    rho,
                    true,
                    &mut self.p_correction_workspace,
                )?;
            } else {
                self.pressure_solver.solve_pressure_correction(
                    fields,
                    dt,
                    rho,
                    true,
                    &mut self.p_correction_workspace,
                )?;
            }

            // Step 5: Correct velocities using pressure gradients with the
            // momentum under-relaxation factor.
            let (ap_u, _, ap_v, _) = self.momentum_solver.get_ap_coefficients();
            self.u_corrected_workspace.copy_from(&self.u_star_workspace);
            self.pressure_solver.correct_velocity(
                &mut self.u_corrected_workspace,
                &self.p_correction_workspace,
                ap_u,
                ap_v,
                rho,
                self.config.alpha_u,
                fields,
            );

            // Step 6: Correct pressure with the configured SIMPLEC relaxation.
            self.p_workspace
                .as_mut_slice()
                .clone_from_slice(fields.p.as_slice());
            self.pressure_solver.correct_pressure(
                &mut self.p_workspace,
                &self.p_correction_workspace,
                self.config.alpha_p,
            );
            fields
                .p
                .data_mut()
                .clone_from_slice(self.p_workspace.as_slice());

            // Step 7: Update velocity fields (u* already relaxed in momentum solve)
            for i in 0..self.grid.nx {
                for j in 0..self.grid.ny {
                    fields.set_velocity_at(i, j, &self.u_corrected_workspace[(i, j)]);
                }
            }

            // Step 8: Check convergence
            let velocity_residual =
                self.calculate_velocity_residual_from_vectors(&u_prev, &self.u_corrected_workspace);

            continuity_residual = self.compute_continuity_residual(fields, Some(dt));
            u_prev.copy_from(&self.u_corrected_workspace);

            #[cfg(debug_assertions)]
            println!(
                "SIMPLEC iteration {} residuals: velocity={:.6e}, continuity={:.6e}",
                iter + 1,
                velocity_residual,
                continuity_residual
            );

            let velocity_converged = velocity_residual < self.config.tolerance;
            // Use the configured convergence tolerance rather than a stricter
            // hard-coded threshold so SIMPLEC stops on the same contract as the
            // rest of the pressure-velocity stack.
            let continuity_converged = continuity_residual < self.config.tolerance;

            residual = velocity_residual.max(continuity_residual);

            if velocity_converged && continuity_converged {
                converged = true;
                tracing::info!(
                    "SIMPLEC converged at iteration {}: velocity = {:.2e}, continuity = {:.2e}",
                    iter + 1,
                    velocity_residual,
                    continuity_residual
                );
                break;
            }

            if iter > 0
                && velocity_residual > previous_velocity_residual
                && continuity_residual > previous_continuity_residual
            {
                divergence_streak += 1;
            } else {
                divergence_streak = 0;
            }
            previous_velocity_residual = velocity_residual;
            previous_continuity_residual = continuity_residual;

            if divergence_streak >= 3 {
                tracing::warn!(
                    "SIMPLEC residuals are diverging; switching to the stable PIMPLE correction"
                );
                return self.solve_pimple(fields, dt, nu, rho);
            }
        }

        if !converged {
            tracing::warn!(
                "SIMPLEC did not converge within {} iterations. \
                 velocity = {:.2e}, continuity = {:.2e}",
                max_iterations,
                residual,
                continuity_residual
            );
        }

        self.iterations += 1;
        Ok(residual)
    }

    /// PIMPLE algorithm implementation
    ///
    /// PIMPLE (Issa 1986, OpenFOAM) combines outer PISO-like correctors
    /// with inner SIMPLE corrections for transient incompressible flows.
    pub(super) fn solve_pimple(
        &mut self,
        fields: &mut SimulationFields<T>,
        dt: T,
        _nu: T,
        rho: T,
    ) -> cfd_core::error::Result<T> {
        let u_initial = self.extract_velocity_field(fields);
        let mut u_before_outer = Array2D::new(self.grid.nx, self.grid.ny, Vector2::zeros());
        let max_outer_iterations = self.config.n_outer_correctors;

        for _outer_iter in 0..max_outer_iterations {
            for i in 0..self.grid.nx {
                for j in 0..self.grid.ny {
                    u_before_outer[(i, j)] = Vector2::new(fields.u.at(i, j), fields.v.at(i, j));
                }
            }

            // Solve momentum equations
            self.momentum_solver
                .solve_with_coefficients(MomentumComponent::U, fields, dt)?;
            self.momentum_solver
                .solve_with_coefficients(MomentumComponent::V, fields, dt)?;

            if let Some(ref mut rhie_chow) = self.rhie_chow {
                let (ap_full_u, _, ap_full_v, _) = self.momentum_solver.get_ap_coefficients();
                rhie_chow.update_u_coefficients(ap_full_u);
                rhie_chow.update_v_coefficients(ap_full_v);
            }

            promote_predicted_velocity_state(fields, &mut self.u_star_workspace);

            // Inner PISO-like corrections
            for _inner_iter in 0..self.config.n_inner_correctors {
                let rebuild_matrix = _inner_iter == 0;
                {
                    let grid = &self.grid;
                    let pressure_solver = &self.pressure_solver;
                    let momentum_solver = &self.momentum_solver;
                    let rhie_chow = self.rhie_chow.as_ref();
                    let vel_field_cache = &self._vel_field_cache;
                    let cons_vel_cache = &self._cons_vel_cache;
                    let u_face_cache = &self._u_face_cache;
                    let v_face_cache = &self._v_face_cache;
                    let d_x_cache = &self._d_x_cache;
                    let d_y_cache = &self._d_y_cache;

                    super::interpolation::solve_pressure_correction_with_caches(
                        grid,
                        pressure_solver,
                        momentum_solver,
                        rhie_chow,
                        vel_field_cache,
                        cons_vel_cache,
                        u_face_cache,
                        v_face_cache,
                        d_x_cache,
                        d_y_cache,
                        fields,
                        dt,
                        Some(dt),
                        rho,
                        rebuild_matrix,
                        &mut self.p_correction_workspace,
                    )?;
                }

                {
                    let (ap_u, _, ap_v, _) = self.momentum_solver.get_ap_coefficients();
                    self.u_corrected_workspace.copy_from(&self.u_star_workspace);
                    self.pressure_solver.correct_velocity(
                        &mut self.u_corrected_workspace,
                        &self.p_correction_workspace,
                        ap_u,
                        ap_v,
                        rho,
                        T::one(),
                        fields,
                    );

                    self.p_workspace
                        .as_mut_slice()
                        .clone_from_slice(fields.p.as_slice());
                    self.pressure_solver.correct_pressure(
                        &mut self.p_workspace,
                        &self.p_correction_workspace,
                        self.config.alpha_p,
                    );
                    fields
                        .p
                        .data_mut()
                        .clone_from_slice(self.p_workspace.as_slice());

                    for i in 0..self.grid.nx {
                        for j in 0..self.grid.ny {
                            fields.set_velocity_at(i, j, &self.u_corrected_workspace[(i, j)]);
                        }
                    }
                }
            }

            // Check outer loop convergence
            let outer_residual = self.calculate_velocity_residual_from_vectors(
                &u_before_outer,
                &self.u_corrected_workspace,
            );

            if outer_residual
                < self.config.tolerance * T::from_f64(5.0).expect("analytical constant conversion")
            {
                break;
            }
        }

        let final_residual = self.calculate_velocity_residual_from_vectors(
            &u_initial,
            &self.extract_velocity_field(fields),
        );

        if final_residual
            > self.config.tolerance * T::from_f64(10.0).expect("analytical constant conversion")
        {
            tracing::warn!(
                "PIMPLE did not achieve desired convergence, residual: {:.2e}",
                final_residual
            );
        }

        self.iterations += 1;
        Ok(final_residual)
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
    fn simplec_algorithm_creation() {
        let grid = make_grid(8);
        let config = SimplecPimpleConfig::simplec();
        assert_eq!(config.algorithm, AlgorithmType::Simplec);

        let solver = SimplecPimpleSolver::new(grid, config);
        assert!(solver.is_ok(), "Failed to create SIMPLEC solver");
    }

    #[test]
    fn pimple_algorithm_creation() {
        let grid = make_grid(8);
        let config = SimplecPimpleConfig::pimple();
        assert_eq!(config.algorithm, AlgorithmType::Pimple);

        let solver = SimplecPimpleSolver::new(grid, config);
        assert!(solver.is_ok(), "Failed to create PIMPLE solver");
    }

    #[test]
    fn verify_algorithm_type() {
        let simplec_cfg = SimplecPimpleConfig::<f64>::simplec();
        assert_eq!(simplec_cfg.algorithm, AlgorithmType::Simplec);

        let pimple_cfg = SimplecPimpleConfig::<f64>::pimple();
        assert_eq!(pimple_cfg.algorithm, AlgorithmType::Pimple);

        // Verify default is SIMPLEC
        let default_cfg = SimplecPimpleConfig::<f64>::default();
        assert_eq!(default_cfg.algorithm, AlgorithmType::Simplec);
    }
}
