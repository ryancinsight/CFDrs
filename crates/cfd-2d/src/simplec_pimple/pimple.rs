use super::solver::SimplecPimpleSolver;
use crate::fields::SimulationFields;
use crate::grid::array2d::Array2D;
use crate::physics::MomentumComponent;
use nalgebra::{RealField, Vector2};
use num_traits::{FromPrimitive, ToPrimitive};

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + std::fmt::LowerExp>
    SimplecPimpleSolver<T>
{
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
        self.extrapolate_pressure_to_solids(fields);
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

            Self::promote_predicted_velocity_state(fields, &mut self.u_star_workspace);

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

                    let face_dt = if dt > T::from_f64(1.0).unwrap_or_else(T::one) {
                        None
                    } else {
                        Some(dt)
                    };
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
                        face_dt,
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

                    Self::apply_field_velocity_boundaries(
                        &self.grid,
                        self.momentum_solver.boundary_conditions(),
                        &mut self.u_corrected_workspace,
                    );

                    for i in 0..self.grid.nx {
                        for j in 0..self.grid.ny {
                            self.p_workspace[(i, j)] = fields.p.at(i, j);
                        }
                    }
                    self.pressure_solver.correct_pressure(
                        &mut self.p_workspace,
                        &self.p_correction_workspace,
                        self.config.alpha_p,
                    );
                    for i in 0..self.grid.nx {
                        for j in 0..self.grid.ny {
                            fields.p.set(i, j, self.p_workspace[(i, j)]);
                        }
                    }
                    self.extrapolate_pressure_to_solids(fields);

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
