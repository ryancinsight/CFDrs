use crate::grid::array2d::Array2D;
use cfd_core::error::Result;
use cfd_core::physics::boundary::BoundaryCondition;
use nalgebra::RealField;
use std::collections::HashMap;

/// Energy equation solver for transporting thermal scalar fields.
pub struct EnergyEquationSolver<T: RealField + Copy> {
    /// Temperature field
    pub temperature: Array2D<T>,
    /// Thermal diffusivity field
    pub thermal_diffusivity: Array2D<T>,
    /// Heat source term
    pub heat_source: Array2D<T>,
    /// Grid dimensions
    pub nx: usize,
    /// Y dimension of the domain.
    pub ny: usize,
    /// Reusable work buffer for explicit time stepping (avoids per-call allocation)
    work_buffer: Array2D<T>,
    /// Whether to include viscous dissipation as an additional heat source.
    ///
    /// When enabled, the viscous dissipation function Phi/(rho*Cp) is added to the
    /// heat source before each time step. This requires velocity gradients and
    /// dynamic viscosity to be provided via [`Self::set_viscous_dissipation_params`].
    ///
    /// Default: `false` (no change to existing behavior).
    include_viscous_dissipation: bool,
    /// Dynamic viscosity field [Pa·s] for viscous dissipation computation.
    /// Only used when `include_viscous_dissipation` is true.
    viscous_dissipation_mu: Array2D<f64>,
    /// Density [kg/m³] for viscous dissipation normalization (Phi / (rho * Cp)).
    /// Only used when `include_viscous_dissipation` is true.
    viscous_dissipation_rho: f64,
    /// Specific heat capacity [J/(kg·K)] for viscous dissipation normalization.
    /// Only used when `include_viscous_dissipation` is true.
    viscous_dissipation_cp: f64,
}

impl<T: RealField + Copy> EnergyEquationSolver<T> {
    /// Create new energy equation solver
    pub fn new(nx: usize, ny: usize, initial_temperature: T, thermal_diffusivity: T) -> Self {
        Self {
            temperature: Array2D::new(nx, ny, initial_temperature),
            thermal_diffusivity: Array2D::new(nx, ny, thermal_diffusivity),
            heat_source: Array2D::new(nx, ny, T::zero()),
            nx,
            ny,
            work_buffer: Array2D::new(nx, ny, T::zero()),
            include_viscous_dissipation: false,
            viscous_dissipation_mu: Array2D::new(nx, ny, 0.0),
            viscous_dissipation_rho: 1000.0,
            viscous_dissipation_cp: 4186.0,
        }
    }

    /// Enable or disable viscous dissipation as an additional heat source.
    ///
    /// When enabled, the solver adds Phi/(rho*Cp) to the heat source at each
    /// interior point before the explicit time step, where Phi is the viscous
    /// dissipation function computed from velocity gradients.
    ///
    /// Default: `false`.
    pub fn set_include_viscous_dissipation(&mut self, enable: bool) {
        self.include_viscous_dissipation = enable;
    }

    /// Returns whether viscous dissipation is currently enabled.
    pub fn include_viscous_dissipation(&self) -> bool {
        self.include_viscous_dissipation
    }

    /// Set the fluid properties needed for viscous dissipation computation.
    ///
    /// # Arguments
    /// * `mu` - Uniform dynamic viscosity [Pa·s]
    /// * `rho` - Density [kg/m³]
    /// * `cp` - Specific heat capacity [J/(kg·K)]
    pub fn set_viscous_dissipation_params(&mut self, mu: f64, rho: f64, cp: f64) {
        self.viscous_dissipation_mu.fill(mu);
        self.viscous_dissipation_rho = rho;
        self.viscous_dissipation_cp = cp;
    }

    /// Solve energy equation using explicit time stepping
    pub fn solve_explicit(
        &mut self,
        u_velocity: &Array2D<T>,
        v_velocity: &Array2D<T>,
        dt: T,
        dx: T,
        dy: T,
        boundary_conditions: &HashMap<(usize, usize), BoundaryCondition<T>>,
    ) -> Result<()> {
        // Reuse pre-allocated work buffer instead of cloning temperature each call
        self.work_buffer.copy_from(&self.temperature);

        // Apply periodic boundary conditions first (before interior computation)
        // Left boundaries copy from right
        for j in 0..self.ny {
            if let Some(BoundaryCondition::Periodic { .. }) = boundary_conditions.get(&(0, j)) {
                self.work_buffer[(0, j)] = self.temperature[(self.nx - 1, j)];
            }
        }
        // Right boundaries copy from left
        for j in 0..self.ny {
            if let Some(BoundaryCondition::Periodic { .. }) =
                boundary_conditions.get(&(self.nx - 1, j))
            {
                self.work_buffer[(self.nx - 1, j)] = self.temperature[(0, j)];
            }
        }
        // Top boundaries copy from bottom
        for i in 0..self.nx {
            if let Some(BoundaryCondition::Periodic { .. }) = boundary_conditions.get(&(i, 0)) {
                self.work_buffer[(i, 0)] = self.temperature[(i, self.ny - 1)];
            }
        }
        // Bottom boundaries copy from top
        for i in 0..self.nx {
            if let Some(BoundaryCondition::Periodic { .. }) =
                boundary_conditions.get(&(i, self.ny - 1))
            {
                self.work_buffer[(i, self.ny - 1)] = self.temperature[(i, 0)];
            }
        }

        // Optionally add viscous dissipation to heat source.
        // This is controlled by a flag to avoid changing default behavior.
        // We compute Phi/(rho*Cp) from velocity gradients and add it to the
        // heat source field. The contribution is removed after the time step
        // to avoid accumulation across calls.
        //
        // Computation uses T-valued arithmetic throughout for generic
        // compatibility. The mu/rho/cp parameters are stored as f64 but
        // only consumed in the f64-specialized path (see below).
        let viscous_dissipation_contributions: Option<Array2D<T>> =
            if self.include_viscous_dissipation {
                let two = T::one() + T::one();
                let two_dx = two * dx;
                let two_dy = two * dy;
                let mut contribs = Array2D::new(self.nx, self.ny, T::zero());
                for i in 1..self.nx - 1 {
                    for j in 1..self.ny - 1 {
                        let du_dx = (u_velocity[(i + 1, j)] - u_velocity[(i - 1, j)]) / two_dx;
                        let du_dy = (u_velocity[(i, j + 1)] - u_velocity[(i, j - 1)]) / two_dy;
                        let dv_dx = (v_velocity[(i + 1, j)] - v_velocity[(i - 1, j)]) / two_dx;
                        let dv_dy = (v_velocity[(i, j + 1)] - v_velocity[(i, j - 1)]) / two_dy;

                        // Phi = 2*mu*[(du/dx)^2 + (dv/dy)^2] + mu*(du/dy + dv/dx)^2
                        // Compute in T-space using stored f64 mu converted via
                        // repeated addition (safe for any RealField).
                        // For simplicity, use the viscous_dissipation_2d function
                        // on f64 and store the contribution to add/remove.
                        let phi = two * (du_dx * du_dx + dv_dy * dv_dy)
                            + (du_dy + dv_dx) * (du_dy + dv_dx);
                        // phi is now the dissipation per unit viscosity (mu=1).
                        // We scale by mu/(rho*cp) but these are f64 stored params.
                        // Store the T-valued phi; the f64 scaling is applied below.
                        contribs[(i, j)] = phi;
                    }
                }
                Some(contribs)
            } else {
                None
            };

        // Apply viscous dissipation to heat_source (if enabled)
        if let Some(ref contribs) = viscous_dissipation_contributions {
            let rho_cp = self.viscous_dissipation_rho * self.viscous_dissipation_cp;
            if rho_cp > 0.0 {
                for i in 1..self.nx - 1 {
                    for j in 1..self.ny - 1 {
                        // Scale: mu / (rho * cp) * phi_normalized
                        // phi_normalized = 2*(du_dx^2 + dv_dy^2) + (du_dy+dv_dx)^2
                        // full contribution = mu * phi_normalized / (rho * cp)
                        let scale = self.viscous_dissipation_mu[(i, j)] / rho_cp;
                        // Convert scale factor from f64 to T via successive halvings.
                        // For f64 T this is exact; for f32 it rounds.
                        self.heat_source[(i, j)] += contribs[(i, j)] * T::from_subset(&scale);
                    }
                }
            }
        }

        // Interior points only - boundaries are handled separately
        for i in 1..self.nx - 1 {
            for j in 1..self.ny - 1 {
                // Skip any interior points that might have boundary conditions (unusual but possible)
                if boundary_conditions.contains_key(&(i, j)) {
                    continue;
                }

                let t = self.temperature[(i, j)];
                let alpha = self.thermal_diffusivity[(i, j)];

                // Face-interpolated temperatures (Central Differencing)
                let two = T::one() + T::one();
                let t_west = self.temperature[(i - 1, j)];
                let t_east = self.temperature[(i + 1, j)];
                let t_south = self.temperature[(i, j - 1)];
                let t_north = self.temperature[(i, j + 1)];

                // Face-interpolated velocities
                let u_west = (u_velocity[(i - 1, j)] + u_velocity[(i, j)]) / two;
                let u_east = (u_velocity[(i, j)] + u_velocity[(i + 1, j)]) / two;
                let v_south = (v_velocity[(i, j - 1)] + v_velocity[(i, j)]) / two;
                let v_north = (v_velocity[(i, j)] + v_velocity[(i, j + 1)]) / two;

                // Conservative convective fluxes across cell faces
                let f_conv_east = u_east * (t + t_east) / two;
                let f_conv_west = u_west * (t_west + t) / two;
                let f_conv_north = v_north * (t + t_north) / two;
                let f_conv_south = v_south * (t_south + t) / two;

                // Conservative diffusive fluxes across cell faces
                // (Using local alpha for the cell face gradient)
                let f_diff_east = alpha * (t_east - t) / dx;
                let f_diff_west = alpha * (t - t_west) / dx;
                let f_diff_north = alpha * (t_north - t) / dy;
                let f_diff_south = alpha * (t - t_south) / dy;

                // Flux balance (Finite Volume Method)
                let conv_term =
                    (f_conv_east - f_conv_west) / dx + (f_conv_north - f_conv_south) / dy;
                let diff_term =
                    (f_diff_east - f_diff_west) / dx + (f_diff_north - f_diff_south) / dy;

                // Explicit update
                self.work_buffer[(i, j)] =
                    t + dt * (-conv_term + diff_term + self.heat_source[(i, j)]);
            }
        }

        // Remove viscous dissipation contributions from heat_source to avoid
        // accumulation across successive calls to solve_explicit().
        if let Some(ref contribs) = viscous_dissipation_contributions {
            let rho_cp = self.viscous_dissipation_rho * self.viscous_dissipation_cp;
            if rho_cp > 0.0 {
                for i in 1..self.nx - 1 {
                    for j in 1..self.ny - 1 {
                        let scale = self.viscous_dissipation_mu[(i, j)] / rho_cp;
                        self.heat_source[(i, j)] -= contribs[(i, j)] * T::from_subset(&scale);
                    }
                }
            }
        }

        // Apply remaining boundary conditions (non-periodic)
        for (&(i, j), bc) in boundary_conditions {
            match bc {
                BoundaryCondition::Dirichlet { value, .. } => {
                    self.work_buffer[(i, j)] = *value;
                }
                BoundaryCondition::Neumann { gradient } => {
                    // Apply gradient boundary condition using interior points
                    if i == 0 {
                        self.work_buffer[(0, j)] = self.work_buffer[(1, j)] - *gradient * dx;
                    } else if i == self.nx - 1 {
                        self.work_buffer[(self.nx - 1, j)] =
                            self.work_buffer[(self.nx - 2, j)] + *gradient * dx;
                    }
                    if j == 0 {
                        self.work_buffer[(i, 0)] = self.work_buffer[(i, 1)] - *gradient * dy;
                    } else if j == self.ny - 1 {
                        self.work_buffer[(i, self.ny - 1)] =
                            self.work_buffer[(i, self.ny - 2)] + *gradient * dy;
                    }
                }
                BoundaryCondition::Symmetry => {
                    // Symmetry BC: zero normal gradient
                    if i == 0 {
                        self.work_buffer[(0, j)] = self.work_buffer[(1, j)];
                    } else if i == self.nx - 1 {
                        self.work_buffer[(self.nx - 1, j)] = self.work_buffer[(self.nx - 2, j)];
                    }
                    if j == 0 {
                        self.work_buffer[(i, 0)] = self.work_buffer[(i, 1)];
                    } else if j == self.ny - 1 {
                        self.work_buffer[(i, self.ny - 1)] = self.work_buffer[(i, self.ny - 2)];
                    }
                }
                _ => {}
            }
        }

        std::mem::swap(&mut self.temperature, &mut self.work_buffer);
        Ok(())
    }

    /// Calculate Nusselt number for heat transfer analysis
    pub fn nusselt_number(&self, wall_temp: T, bulk_temp: T, characteristic_length: T, dy: T) -> T {
        let dt_dy_wall = (self.temperature[(0, 1)] - self.temperature[(0, 0)]) / dy; // Use actual dy
        let h = self.thermal_diffusivity[(0, 0)] * dt_dy_wall / (wall_temp - bulk_temp);
        h * characteristic_length / self.thermal_diffusivity[(0, 0)]
    }
}
