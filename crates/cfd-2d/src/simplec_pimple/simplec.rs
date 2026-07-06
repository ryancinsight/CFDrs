use super::solver::SimplecPimpleSolver;
use crate::fields::SimulationFields;
use crate::grid::array2d::Array2D;
use crate::physics::MomentumComponent;
use crate::scalar;
use crate::scalar::Cfd2dScalar;
use eunomia::FloatElement;
use leto::geometry::Vector2;

impl<T: Cfd2dScalar + Copy + std::fmt::LowerExp + FloatElement> SimplecPimpleSolver<T> {
    /// SIMPLEC (Van Doormaal & Raithby, 1984) improves SIMPLE by using
    /// consistent discretization for the pressure correction equation.
    pub(super) fn solve_simplec(
        &mut self,
        fields: &mut SimulationFields<T>,
        dt: T,
        nu: T,
        rho: T,
    ) -> cfd_core::error::Result<T> {
        let mut residual = scalar::zero();
        let max_iterations = self.config.max_inner_iterations;
        let mut converged = false;

        let mut u_prev = self.extract_velocity_field(fields);
        let mut continuity_residual = scalar::from_f64(1e30);
        let mut previous_velocity_residual = scalar::from_f64(1e30);
        let mut previous_continuity_residual = continuity_residual;
        let mut divergence_streak = 0usize;

        self.extrapolate_pressure_to_solids(fields);

        for iter in 0..max_iterations {
            // Step 1: Solve momentum equations for predicted velocities u*
            self.momentum_solver
                .solve_with_coefficients(MomentumComponent::U, fields, dt)?;
            self.momentum_solver
                .solve_with_coefficients(MomentumComponent::V, fields, dt)?;

            // Step 2: Refresh the Rhie-Chow face coefficients from the raw
            // momentum stencil.
            if let Some(ref mut rhie_chow) = self.rhie_chow {
                let (_, ap_u_consistent, _, ap_v_consistent) =
                    self.momentum_solver.get_ap_coefficients();
                rhie_chow.update_u_coefficients(ap_u_consistent);
                rhie_chow.update_v_coefficients(ap_v_consistent);
            }

            // Step 3: Promote the predicted momentum state into the current
            // fields so the pressure solve sees the raw predictor.
            Self::promote_predicted_velocity_state(fields, &mut self.u_star_workspace);

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

                let face_dt = if dt > scalar::from_f64(1.0) {
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
                    true,
                    &mut self.p_correction_workspace,
                )?;
            } else {
                self.pressure_solver.solve_pressure_correction(
                    fields,
                    dt,
                    rho,
                    Some(self.momentum_solver.boundary_conditions()),
                    true,
                    &mut self.p_correction_workspace,
                )?;
            }

            // Step 5: Correct velocities using pressure gradients with the
            // momentum under-relaxation factor.
            let (_, ap_u_consistent, _, ap_v_consistent) =
                self.momentum_solver.get_ap_coefficients();
            self.u_corrected_workspace.copy_from(&self.u_star_workspace);
            self.pressure_solver.correct_velocity(
                &mut self.u_corrected_workspace,
                &self.p_correction_workspace,
                ap_u_consistent,
                ap_v_consistent,
                rho,
                self.config.alpha_u,
                fields,
            );

            Self::apply_field_velocity_boundaries(
                &self.grid,
                self.momentum_solver.boundary_conditions(),
                &mut self.u_corrected_workspace,
            );

            // Step 6: Correct pressure with the configured SIMPLEC relaxation.
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

    pub(super) fn apply_field_velocity_boundaries(
        grid: &crate::grid::StructuredGrid2D<T>,
        boundary_conditions: &std::collections::HashMap<
            String,
            cfd_core::physics::boundary::BoundaryCondition<T>,
        >,
        u_star: &mut Array2D<Vector2<T>>,
    ) {
        use cfd_core::physics::boundary::BoundaryCondition;
        let nx = grid.nx;
        let ny = grid.ny;

        // Apply boundary conditions to the boundaries of u_star
        // West boundary
        if let Some(bc) = boundary_conditions.get("west") {
            for j in 0..ny {
                match bc {
                    BoundaryCondition::Neumann { .. }
                    | BoundaryCondition::Outflow
                    | BoundaryCondition::PressureOutlet { .. }
                    | BoundaryCondition::PressureInlet { .. } => {
                        u_star[(0, j)] = u_star[(1, j)];
                    }
                    BoundaryCondition::Symmetry => {
                        u_star[(0, j)][0] = scalar::zero();
                        u_star[(0, j)][1] = u_star[(1, j)][1];
                    }
                    BoundaryCondition::Wall {
                        wall_type: cfd_core::physics::boundary::WallType::Slip,
                    } => {
                        u_star[(0, j)][0] = scalar::zero();
                        u_star[(0, j)][1] = u_star[(1, j)][1];
                    }
                    BoundaryCondition::Periodic { .. } => {
                        u_star[(0, j)] = u_star[(nx - 2, j)];
                    }
                    BoundaryCondition::VelocityInlet { velocity } => {
                        u_star[(0, j)] = Vector2::new(velocity[0], velocity[1]);
                    }
                    _ => {}
                }
            }
        }

        // East boundary
        if let Some(bc) = boundary_conditions.get("east") {
            for j in 0..ny {
                match bc {
                    BoundaryCondition::Neumann { .. }
                    | BoundaryCondition::Outflow
                    | BoundaryCondition::PressureOutlet { .. }
                    | BoundaryCondition::PressureInlet { .. } => {
                        u_star[(nx - 1, j)] = u_star[(nx - 2, j)];
                    }
                    BoundaryCondition::Symmetry => {
                        u_star[(nx - 1, j)][0] = scalar::zero();
                        u_star[(nx - 1, j)][1] = u_star[(nx - 2, j)][1];
                    }
                    BoundaryCondition::Wall {
                        wall_type: cfd_core::physics::boundary::WallType::Slip,
                    } => {
                        u_star[(nx - 1, j)][0] = scalar::zero();
                        u_star[(nx - 1, j)][1] = u_star[(nx - 2, j)][1];
                    }
                    BoundaryCondition::Periodic { .. } => {
                        u_star[(nx - 1, j)] = u_star[(1, j)];
                    }
                    BoundaryCondition::VelocityInlet { velocity } => {
                        u_star[(nx - 1, j)] = Vector2::new(velocity[0], velocity[1]);
                    }
                    _ => {}
                }
            }
        }

        // South boundary
        if let Some(bc) = boundary_conditions.get("south") {
            for i in 0..nx {
                match bc {
                    BoundaryCondition::Neumann { .. }
                    | BoundaryCondition::Outflow
                    | BoundaryCondition::PressureOutlet { .. }
                    | BoundaryCondition::PressureInlet { .. } => {
                        u_star[(i, 0)] = u_star[(i, 1)];
                    }
                    BoundaryCondition::Symmetry => {
                        u_star[(i, 0)][1] = scalar::zero();
                        u_star[(i, 0)][0] = u_star[(i, 1)][0];
                    }
                    BoundaryCondition::Wall {
                        wall_type: cfd_core::physics::boundary::WallType::Slip,
                    } => {
                        u_star[(i, 0)][1] = scalar::zero();
                        u_star[(i, 0)][0] = u_star[(i, 1)][0];
                    }
                    BoundaryCondition::Periodic { .. } => {
                        u_star[(i, 0)] = u_star[(i, ny - 2)];
                    }
                    BoundaryCondition::VelocityInlet { velocity } => {
                        u_star[(i, 0)] = Vector2::new(velocity[0], velocity[1]);
                    }
                    _ => {}
                }
            }
        }

        // North boundary
        if let Some(bc) = boundary_conditions.get("north") {
            for i in 0..nx {
                match bc {
                    BoundaryCondition::Neumann { .. }
                    | BoundaryCondition::Outflow
                    | BoundaryCondition::PressureOutlet { .. }
                    | BoundaryCondition::PressureInlet { .. } => {
                        u_star[(i, ny - 1)] = u_star[(i, ny - 2)];
                    }
                    BoundaryCondition::Symmetry => {
                        u_star[(i, ny - 1)][1] = scalar::zero();
                        u_star[(i, ny - 1)][0] = u_star[(i, ny - 2)][0];
                    }
                    BoundaryCondition::Wall {
                        wall_type: cfd_core::physics::boundary::WallType::Slip,
                    } => {
                        u_star[(i, ny - 1)][1] = scalar::zero();
                        u_star[(i, ny - 1)][0] = u_star[(i, ny - 2)][0];
                    }
                    BoundaryCondition::Periodic { .. } => {
                        u_star[(i, ny - 1)] = u_star[(i, 1)];
                    }
                    BoundaryCondition::VelocityInlet { velocity } => {
                        u_star[(i, ny - 1)] = Vector2::new(velocity[0], velocity[1]);
                    }
                    _ => {}
                }
            }
        }
    }
}
