use super::solver::{MomentumComponent, MomentumSolver};
use crate::fields::SimulationFields;
use crate::scalar;
use crate::scalar::Cfd2dScalar;
use cfd_math::linear_solver::preconditioners::IdentityPreconditioner;
use cfd_math::linear_solver::{DirectSparseSolver, IterativeLinearSolver};
use cfd_math::sparse::SparseMatrixBuilder;
use eunomia::{FloatElement, NumericElement};
use leto::Array1;
use leto_ops::norm_l2;

impl<T: Cfd2dScalar + Copy + FloatElement> MomentumSolver<T> {
    /// Solve momentum equation for specified component.
    ///
    /// Integrates with any configured turbulence model by updating
    /// estimated turbulence quantities based on flow physics.
    pub fn solve(
        &mut self,
        component: MomentumComponent,
        fields: &mut SimulationFields<T>,
        dt: T,
    ) -> cfd_core::error::Result<()> {
        self.solve_with_coefficients(component, fields, dt)?;
        Ok(())
    }

    /// Solve momentum equation and return coefficients for Rhie-Chow interpolation
    pub fn solve_with_coefficients(
        &mut self,
        component: MomentumComponent,
        fields: &mut SimulationFields<T>,
        dt: T,
    ) -> cfd_core::error::Result<()> {
        // Compute coefficients
        self.compute_coefficients_into(component, fields, dt)?;

        // Diagnostic: Check if coefficients are non-zero (helps debug false convergence)
        #[cfg(debug_assertions)]
        {
            let coeffs = match component {
                MomentumComponent::U => &self.coeffs_u,
                MomentumComponent::V => &self.coeffs_v,
            };
            let mut nonzero_count = 0;
            for j in 0..self.grid.ny {
                for i in 0..self.grid.nx {
                    if NumericElement::abs(coeffs.ap.at(i, j)) > T::default_epsilon() {
                        nonzero_count += 1;
                    }
                }
            }
            tracing::debug!(
                "Momentum coefficients: {}/{} non-zero entries",
                nonzero_count,
                self.grid.nx * self.grid.ny
            );
        }

        // Assemble linear system topology once and update values
        self.assemble_system_into(component, fields)?;

        let (matrix, rhs) = match component {
            MomentumComponent::U => (
                self.matrix_u.as_ref().unwrap(),
                self.rhs_u.as_ref().unwrap(),
            ),
            MomentumComponent::V => (
                self.matrix_v.as_ref().unwrap(),
                self.rhs_v.as_ref().unwrap(),
            ),
        };

        // Diagnostic: Check matrix and RHS statistics (helps debug false convergence)
        #[cfg(debug_assertions)]
        {
            let matrix_nnz = matrix.nnz();

            tracing::debug!(
                "Linear system: matrix {}x{}, {} nnz",
                matrix.nrows(),
                matrix.ncols(),
                matrix_nnz,
            );

            if matrix_nnz == 0 {
                tracing::error!("CRITICAL: Matrix has ZERO non-zero entries - system is empty!");
            }
        }

        // Solve linear system
        let mut solution = Array1::from_elem([matrix.nrows()], scalar::zero::<T>());
        let _rhs_norm =
            norm_l2(&rhs.view()).expect("invariant: momentum RHS Leto vector has a valid layout");
        let solve_result =
            self.linear_solver
                .solve(matrix, rhs, &mut solution, None::<&IdentityPreconditioner>);

        match solve_result {
            Ok(_) => {}
            Err(cfd_core::error::Error::Convergence(
                cfd_core::error::ConvergenceErrorKind::MaxIterationsExceeded { max },
            )) => {
                tracing::warn!(
                    component = ?component,
                    max_iterations = max,
                    "Momentum GMRES stalled; falling back to direct sparse solve"
                );
                let direct_solver = DirectSparseSolver::default();
                solution = direct_solver.solve(matrix, rhs).map_err(|error| {
                    cfd_core::error::Error::Solver(format!(
                        "Momentum direct sparse solve failed for {component:?}: {error}"
                    ))
                })?;
            }
            Err(cfd_core::error::Error::Convergence(
                cfd_core::error::ConvergenceErrorKind::Breakdown,
            )) => {
                tracing::warn!(
                    component = ?component,
                    "Momentum GMRES breakdown; falling back to direct sparse solve"
                );
                let direct_solver = DirectSparseSolver::default();
                solution = direct_solver.solve(matrix, rhs).map_err(|error| {
                    cfd_core::error::Error::Solver(format!(
                        "Momentum direct sparse solve failed for {component:?}: {error}"
                    ))
                })?;
            }
            Err(error) => return Err(error),
        }

        // Update velocity field
        self.update_velocity(component, fields, &solution);

        Ok(())
    }

    fn assemble_system_into(
        &mut self,
        component: MomentumComponent,
        fields: &SimulationFields<T>,
    ) -> cfd_core::error::Result<()> {
        let coeffs = match component {
            MomentumComponent::U => &self.coeffs_u,
            MomentumComponent::V => &self.coeffs_v,
        };
        let n = self.grid.nx * self.grid.ny;
        let nx = self.grid.nx;
        let ny = self.grid.ny;

        let should_rebuild = match component {
            MomentumComponent::U => self.matrix_u.is_none(),
            MomentumComponent::V => self.matrix_v.is_none(),
        };

        if match component {
            MomentumComponent::U => self.rhs_u.as_ref().is_none_or(|r| r.shape()[0] != n),
            MomentumComponent::V => self.rhs_v.as_ref().is_none_or(|r| r.shape()[0] != n),
        } {
            match component {
                MomentumComponent::U => {
                    self.rhs_u = Some(Array1::from_elem([n], scalar::zero::<T>()));
                }
                MomentumComponent::V => {
                    self.rhs_v = Some(Array1::from_elem([n], scalar::zero::<T>()));
                }
            }
        } else {
            match component {
                MomentumComponent::U => self.rhs_u.as_mut().unwrap().fill(scalar::zero::<T>()),
                MomentumComponent::V => self.rhs_v.as_mut().unwrap().fill(scalar::zero::<T>()),
            }
        }

        if should_rebuild {
            let mut builder = match component {
                MomentumComponent::U => self
                    .matrix_builder_u
                    .take()
                    .unwrap_or_else(|| SparseMatrixBuilder::new(n, n)),
                MomentumComponent::V => self
                    .matrix_builder_v
                    .take()
                    .unwrap_or_else(|| SparseMatrixBuilder::new(n, n)),
            };
            let mut rhs = Array1::from_elem([n], scalar::zero::<T>());
            let mask_data = fields.mask.data.as_slice();
            let ap_data = coeffs.ap.data.as_slice();
            let aw_data = coeffs.aw.data.as_slice();
            let ae_data = coeffs.ae.data.as_slice();
            let as_data = coeffs.as_.data.as_slice();
            let an_data = coeffs.an.data.as_slice();
            let source_data = coeffs.source.data.as_slice();

            for j in 0..ny {
                for i in 0..nx {
                    let idx = j * self.grid.nx + i;

                    // Check if this is a masked (solid) cell
                    if !mask_data[idx] {
                        // For solid cells: assemble identity equation φ = 0
                        builder.add_entry(idx, idx, T::one())?;
                        // RHS is already zero
                        continue;
                    }

                    // Check if this is a Dirichlet boundary node
                    if self.is_dirichlet_boundary(i, j, component) {
                        // For Dirichlet BC: assemble identity equation φ = bc_value
                        // This is handled in apply_momentum_boundaries, so skip coefficient assembly
                        builder.add_entry(idx, idx, T::one())?;
                        // RHS will be set by boundary condition handler
                        continue;
                    }

                    // Check if this is any other boundary node
                    if i == 0 || i == nx - 1 || j == 0 || j == ny - 1 {
                        continue;
                    }

                    // Central coefficient
                    builder.add_entry(idx, idx, ap_data[idx])?;

                    // Neighbor coefficients
                    if i > 0 {
                        builder.add_entry(idx, idx - 1, -aw_data[idx])?;
                    }
                    if i < self.grid.nx - 1 {
                        builder.add_entry(idx, idx + 1, -ae_data[idx])?;
                    }
                    if j > 0 {
                        builder.add_entry(idx, idx - self.grid.nx, -as_data[idx])?;
                    }
                    if j < self.grid.ny - 1 {
                        builder.add_entry(idx, idx + self.grid.nx, -an_data[idx])?;
                    }

                    // Source term including pressure gradient
                    rhs[idx] = source_data[idx];
                }
            }

            // Apply higher-order boundary conditions for improved near-wall accuracy
            // This provides better velocity gradients near walls for production accuracy
            super::boundary::apply_higher_order_wall_boundaries(
                &mut builder,
                &mut rhs,
                component,
                &self.boundary_conditions,
                &self.grid,
            )?;

            // Apply standard boundary conditions (sets RHS values for Dirichlet, modifies equations for Neumann)
            super::boundary::apply_momentum_boundaries(
                &mut builder,
                &mut rhs,
                component,
                &self.boundary_conditions,
                &self.grid,
            )?;

            match component {
                MomentumComponent::U => {
                    self.matrix_u = Some(builder.build()?);
                    self.matrix_builder_u = Some(SparseMatrixBuilder::new(n, n));
                    self.rhs_u = Some(rhs);
                }
                MomentumComponent::V => {
                    self.matrix_v = Some(builder.build()?);
                    self.matrix_builder_v = Some(SparseMatrixBuilder::new(n, n));
                    self.rhs_v = Some(rhs);
                }
            }
        } else {
            let mut matrix = match component {
                MomentumComponent::U => self.matrix_u.take().unwrap(),
                MomentumComponent::V => self.matrix_v.take().unwrap(),
            };
            matrix.values_mut().fill(T::zero());

            let mut rhs = Array1::from_elem([n], scalar::zero::<T>());

            let row_ptr = matrix.row_ptr().to_vec();
            let col_indices = matrix.col_indices().to_vec();
            let values = matrix.values_mut();
            let mut update_entry = |r: usize, c: usize, v: T| {
                let start = row_ptr[r];
                let end = row_ptr[r + 1];
                for idx in start..end {
                    if col_indices[idx] == c {
                        values[idx] += v;
                        break;
                    }
                }
            };

            let mask_data = fields.mask.data.as_slice();
            let ap_data = coeffs.ap.data.as_slice();
            let aw_data = coeffs.aw.data.as_slice();
            let ae_data = coeffs.ae.data.as_slice();
            let as_data = coeffs.as_.data.as_slice();
            let an_data = coeffs.an.data.as_slice();
            let source_data = coeffs.source.data.as_slice();

            for j in 0..ny {
                for i in 0..nx {
                    let idx = j * nx + i;

                    if !mask_data[idx] {
                        update_entry(idx, idx, T::one());
                        continue;
                    }

                    if self.is_dirichlet_boundary(i, j, component) {
                        update_entry(idx, idx, T::one());
                        continue;
                    }

                    // Check if this is any other boundary node
                    if i == 0 || i == nx - 1 || j == 0 || j == ny - 1 {
                        continue;
                    }

                    update_entry(idx, idx, ap_data[idx]);
                    if i > 0 {
                        update_entry(idx, idx - 1, -aw_data[idx]);
                    }
                    if i < nx - 1 {
                        update_entry(idx, idx + 1, -ae_data[idx]);
                    }
                    if j > 0 {
                        update_entry(idx, idx - nx, -as_data[idx]);
                    }
                    if j < ny - 1 {
                        update_entry(idx, idx + nx, -an_data[idx]);
                    }

                    rhs[idx] = source_data[idx];
                }
            }

            super::boundary::apply_higher_order_wall_boundaries(
                &mut matrix,
                &mut rhs,
                component,
                &self.boundary_conditions,
                &self.grid,
            )?;

            super::boundary::apply_momentum_boundaries(
                &mut matrix,
                &mut rhs,
                component,
                &self.boundary_conditions,
                &self.grid,
            )?;

            match component {
                MomentumComponent::U => {
                    self.matrix_u = Some(matrix);
                    self.rhs_u = Some(rhs);
                }
                MomentumComponent::V => {
                    self.matrix_v = Some(matrix);
                    self.rhs_v = Some(rhs);
                }
            }
        }

        Ok(())
    }

    fn update_velocity(
        &self,
        component: MomentumComponent,
        fields: &mut SimulationFields<T>,
        solution: &Array1<T>,
    ) {
        let alpha = self.velocity_relaxation;
        let one_minus_alpha = T::one() - alpha;

        for j in 0..self.grid.ny {
            for i in 0..self.grid.nx {
                let idx = j * self.grid.nx + i;
                let computed_value = solution[idx];

                match component {
                    MomentumComponent::U => {
                        if let Some(u) = fields.u.at_mut(i, j) {
                            let old_value = *u;
                            // Enforce zero velocity in masked (solid) cells
                            if fields.mask.at(i, j) {
                                // Under-relaxation: u_new = α * u_computed + (1-α) * u_old
                                let relaxed_value =
                                    alpha * computed_value + one_minus_alpha * old_value;
                                *u = relaxed_value;
                                if let Some(u_star) = fields.u_star.at_mut(i, j) {
                                    *u_star = relaxed_value;
                                }
                            } else {
                                if let Some(u_star) = fields.u_star.at_mut(i, j) {
                                    *u_star = T::zero();
                                }
                                *u = T::zero();
                            }
                        }
                    }
                    MomentumComponent::V => {
                        if let Some(v) = fields.v.at_mut(i, j) {
                            let old_value = *v;
                            // Enforce zero velocity in masked (solid) cells
                            if fields.mask.at(i, j) {
                                // Under-relaxation: v_new = α * v_computed + (1-α) * v_old
                                let relaxed_value =
                                    alpha * computed_value + one_minus_alpha * old_value;
                                *v = relaxed_value;
                                if let Some(v_star) = fields.v_star.at_mut(i, j) {
                                    *v_star = relaxed_value;
                                }
                            } else {
                                if let Some(v_star) = fields.v_star.at_mut(i, j) {
                                    *v_star = T::zero();
                                }
                                *v = T::zero();
                            }
                        }
                    }
                }
            }
        }
    }
}
