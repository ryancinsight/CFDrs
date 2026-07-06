//! Cavitation Solver Integration with VOF Method
//!
//! This module provides comprehensive cavitation simulation capabilities
//! by integrating cavitation physics with the Volume of Fluid (VOF) method.
//!
//! # Theorem — Cavitation–VOF Mass Conservation
//!
//! The VOF transport equation with cavitation source term
//!
//! ```text
//! ∂α/∂t + ∇·(α u) = R_cav
//! ```
//!
//! conserves the total mass $\int_\Omega (\alpha \rho_l + (1-\alpha)\rho_v) \,d\Omega$
//! provided the cavitation source $R_{cav}$ satisfies
//!
//! ```text
//! ∫_Ω R_cav (ρ_l − ρ_v) dΩ = d/dt ∫_Ω ρ_v dΩ − d/dt ∫_Ω ρ_v dΩ
//! ```
//!
//! i.e., mass is transferred between phases without creation or destruction.
//!
//! ## Mathematical Foundation
//!
//! ### Cavitation-VOF Coupling
//! The cavitation model is coupled with VOF through the void fraction α:
//!
//! ```math
//! ∂α/∂t + ∇·(α u) = R_cav
//! ```
//!
//! where R_cav is the cavitation mass transfer rate from cavitation models.
//!
//! ### Multi-Phase Navier-Stokes with Cavitation
//! The momentum equation includes cavitation-induced momentum transfer:
//!
//! ```math
//! ∂(ρ u)/∂t + ∇·(ρ u ⊗ u) = -∇p + ∇·τ + F_σ + F_cav
//! ```
//!
//! where F_cav is the cavitation momentum source term.
//!
//! ## Implementation Features
//!
//! - **VOF-Cavitation Coupling**: Mass transfer integrated with interface tracking
//! - **Cavitation Damage Prediction**: Erosion rate calculation with material models
//! - **Bubble Dynamics**: Rayleigh-Plesset equation integration
//! - **Multi-Scale Modeling**: From microscopic bubbles to macroscopic damage
//! - **GPU Acceleration**: Optimized for large-scale cavitation simulations

/// Default grid spacing \[m] used when constructing the VOF solver grid.
const DEFAULT_GRID_SPACING: f64 = 0.01;

/// Polytropic index (ratio of specific heats, gamma = c_p / c_v) for
/// air at standard conditions.  Used in the Rayleigh-Plesset adiabatic
/// compression model for bubble-collapse frequency estimation.
#[allow(dead_code)]
const POLYTROPIC_INDEX_AIR: f64 = 1.4;

/// Cavitation inception field marker value.  When the local cavitation
/// number falls below the inception threshold *or* pressure drops below
/// vapor pressure, the inception field is set to this value (1.0 = active).
const INCEPTION_ACTIVE: f64 = 1.0;

/// Threshold applied to the inception field when counting cavitating
/// cells in statistics.  Cells with inception_field > this value are
/// counted as actively cavitating.
#[allow(dead_code)]
const INCEPTION_COUNT_THRESHOLD: f64 = 0.5;

use super::bubble_dynamics::BubbleDynamicsSolver;
use super::cavitation_types::{CavitationStatistics, CavitationVofConfig};
use crate::vof::solver::VofSolver;
use cfd_core::error::Result;
use cfd_core::physics::cavitation::rayleigh_plesset::RayleighPlesset;
use leto::{geometry::Vector3, Array2};

/// Dense row-major cavitation scalar field with logical shape `[nx, ny * nz]`.
pub type CavitationField = Array2<f64>;

#[inline]
fn field_len(field: &CavitationField) -> usize {
    field.size()
}

#[inline]
fn field_shape_matches(field: &CavitationField, rows: usize, columns: usize) -> bool {
    field.shape() == [rows, columns]
}

#[inline]
fn field_offset(column_count: usize, i: usize, col: usize) -> usize {
    i * column_count + col
}

#[inline]
fn alpha_offset(nx: usize, i: usize, col: usize) -> usize {
    col * nx + i
}

#[inline]
fn dense_field_slice(field: &CavitationField) -> &[f64] {
    field
        .as_slice()
        .expect("invariant: cavitation fields are dense row-major Leto Array2 values")
}

#[inline]
fn dense_field_slice_mut(field: &mut CavitationField) -> &mut [f64] {
    field
        .as_slice_mut()
        .expect("invariant: cavitation fields are dense row-major Leto Array2 values")
}

/// Cavitation-VOF solver
pub struct CavitationVofSolver {
    /// VOF solver for interface tracking
    vof_solver: VofSolver<f64>,
    /// Cavitation model
    cavitation_model: cfd_core::physics::cavitation::models::CavitationModel<f64>,
    /// Damage prediction model
    damage_model: Option<cfd_core::physics::cavitation::damage::CavitationDamage<f64>>,
    /// Bubble dynamics solver
    bubble_solver: Option<BubbleDynamicsSolver>,
    /// Configuration
    config: CavitationVofConfig,
    /// Cavitation inception field (tracks where cavitation starts)
    inception_field: CavitationField,
    /// Accumulated damage field
    damage_field: Option<CavitationField>,
    /// Bubble radius field for Rayleigh-Plesset
    bubble_radius_field: Option<CavitationField>,
    /// Passive scalar tracking for cavitation nuclei cascade
    nuclei_field: Option<CavitationField>,
    /// Reusable workspace for nuclei advection.
    nuclei_workspace: Option<CavitationField>,
    /// Reusable workspace for cavitation source terms.
    cavitation_source_workspace: CavitationField,
    /// Nuclei transport evaluator
    nuclei_transport: Option<cfd_core::physics::cavitation::nuclei_transport::NucleiTransport<f64>>,
}

impl CavitationVofSolver {
    /// Create new cavitation-VOF solver
    pub fn new(nx: usize, ny: usize, nz: usize, config: CavitationVofConfig) -> Result<Self> {
        if config.damage_model.is_some() && config.bubble_dynamics.is_none() {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "cavitation damage requires bubble dynamics to provide a physically grounded bubble radius".to_string(),
            ));
        }
        if let Some(nuclei) = &config.nuclei_transport {
            if nuclei.diffusion_coefficient < 0.0 {
                return Err(cfd_core::error::Error::InvalidConfiguration(
                    "nuclei diffusion coefficient must be nonnegative".to_string(),
                ));
            }
        }

        let vof_solver = VofSolver::create(
            config.vof_config.clone(),
            nx,
            ny,
            nz,
            DEFAULT_GRID_SPACING,
            DEFAULT_GRID_SPACING,
            DEFAULT_GRID_SPACING,
        );

        let bubble_solver = config.bubble_dynamics.as_ref().map(|bubble_config| {
            BubbleDynamicsSolver::new(
                bubble_config,
                nx,
                ny,
                nz,
                vof_solver.dx,
                vof_solver.dy,
                vof_solver.dz,
                config.liquid_density,
                config.liquid_blood_model.clone(),
                config.vapor_pressure,
            )
        });

        let inception_field = CavitationField::zeros([nx, ny * nz]);
        let damage_field = if config.damage_model.is_some() {
            Some(CavitationField::zeros([nx, ny * nz]))
        } else {
            None
        };

        let bubble_radius_field = if config.bubble_dynamics.is_some() {
            Some(CavitationField::zeros([nx, ny * nz]))
        } else {
            None
        };

        let nuclei_field = if config.nuclei_transport.is_some() {
            Some(CavitationField::zeros([nx, ny * nz]))
        } else {
            None
        };
        let nuclei_workspace = if config.nuclei_transport.is_some() {
            Some(CavitationField::zeros([nx, ny * nz]))
        } else {
            None
        };
        let cavitation_source_workspace = CavitationField::zeros([nx, ny * nz]);

        let nuclei_transport = config
            .nuclei_transport
            .map(cfd_core::physics::cavitation::nuclei_transport::NucleiTransport::new);

        Ok(Self {
            vof_solver,
            cavitation_model: config.cavitation_model.clone(),
            damage_model: config.damage_model.clone(),
            bubble_solver,
            config,
            inception_field,
            damage_field,
            bubble_radius_field,
            nuclei_field,
            nuclei_workspace,
            cavitation_source_workspace,
            nuclei_transport,
        })
    }

    /// Update cavitation-VOF solution for one time step
    pub fn step(
        &mut self,
        dt: f64,
        velocity_field: &[Vector3<f64>],
        pressure_field: &CavitationField,
        density_field: &CavitationField,
    ) -> Result<()> {
        self.validate_flow_field_dimensions(velocity_field, pressure_field, density_field)?;

        // 1. Update bubble dynamics if enabled
        self.update_bubble_dynamics(dt, velocity_field, pressure_field, density_field)?;

        // 2. Calculate cavitation mass transfer rates into the reusable workspace.
        let mut cavitation_source = CavitationField::zeros([0, 0]);
        std::mem::swap(
            &mut cavitation_source,
            &mut self.cavitation_source_workspace,
        );
        self.calculate_cavitation_source_into(
            dt,
            velocity_field,
            pressure_field,
            density_field,
            &mut cavitation_source,
        )?;

        // 2b. Advect and diffuse cavitation nuclei tracking the active cavitation
        self.update_nuclei_advection_diffusion(dt, velocity_field, &cavitation_source)?;

        // 3. Update VOF with cavitation source term
        let nx = self.vof_solver.nx;
        let column_count = self.vof_solver.ny * self.vof_solver.nz;
        let cavitation_source_data = dense_field_slice(&cavitation_source);
        for col in 0..column_count {
            for i in 0..nx {
                let alpha_idx = alpha_offset(nx, i, col);
                let field_idx = field_offset(column_count, i, col);
                self.vof_solver.alpha[alpha_idx] = (self.vof_solver.alpha[alpha_idx]
                    + cavitation_source_data[field_idx])
                    .clamp(0.0, 1.0);
            }
        }

        std::mem::swap(
            &mut cavitation_source,
            &mut self.cavitation_source_workspace,
        );

        // 4. Update VOF velocity and advance.
        self.vof_solver.velocity.clone_from_slice(velocity_field);
        self.vof_solver.advance(dt)?;

        // 5. Calculate cavitation damage if damage model is enabled
        self.update_damage(dt, pressure_field, density_field)?;

        Ok(())
    }

    /// # Theorem - Slice-Consistent Nuclei Transport Update
    ///
    /// If the nuclei transport state matches the solver grid, the explicit
    /// advection-diffusion-reaction update reads and writes only within the
    /// Leto row-major field slice allocated for `[nx, ny * nz]`. The final
    /// pointwise clamp enforces nonnegative nuclei values.
    ///
    /// **Proof sketch**: each physical cell maps to `field_offset(ny * nz, i,
    /// j + k * ny)` in the Leto field and `alpha_offset(nx, i, j + k * ny)` in
    /// the VOF flat vector. Neighbor coordinates are constructed only after
    /// checking the corresponding domain boundary. This keeps every access
    /// within the allocated domain, and the `max(0.0)` clamp enforces the
    /// nonnegative state constraint cell by cell.
    fn update_nuclei_advection_diffusion(
        &mut self,
        dt: f64,
        velocity_field: &[Vector3<f64>],
        cavitation_source: &CavitationField,
    ) -> Result<()> {
        let nx = self.vof_solver.nx;
        let ny = self.vof_solver.ny;
        let nz = self.vof_solver.nz;
        let dx = self.vof_solver.dx;
        let dy = self.vof_solver.dy;
        let dz = self.vof_solver.dz;
        let column_count = ny * nz;
        let grid_size = nx * column_count;

        if !field_shape_matches(cavitation_source, nx, column_count) {
            return Err(cfd_core::error::Error::DimensionMismatch {
                expected: grid_size,
                actual: field_len(cavitation_source),
            });
        }
        if velocity_field.len() != grid_size {
            return Err(cfd_core::error::Error::DimensionMismatch {
                expected: grid_size,
                actual: velocity_field.len(),
            });
        }

        let (Some(nuclei), Some(next_nuclei), Some(transport)) = (
            &mut self.nuclei_field,
            &mut self.nuclei_workspace,
            &self.nuclei_transport,
        ) else {
            return Ok(());
        };

        if !field_shape_matches(nuclei, nx, column_count)
            || !field_shape_matches(next_nuclei, nx, column_count)
        {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "nuclei transport workspace dimensions must match the solver grid".to_string(),
            ));
        }

        let diffusion_coefficient = transport.diffusion_coefficient();
        let inv_dx2 = if dx > 0.0 { 1.0 / (dx * dx) } else { 0.0 };
        let inv_dy2 = if dy > 0.0 { 1.0 / (dy * dy) } else { 0.0 };
        let inv_dz2 = if dz > 0.0 { 1.0 / (dz * dz) } else { 0.0 };
        let nuclei_data = dense_field_slice(nuclei);
        let next_data = dense_field_slice_mut(next_nuclei);
        let source_data = dense_field_slice(cavitation_source);

        next_data.copy_from_slice(nuclei_data);

        for k in 0..nz {
            let col_offset = k * ny;
            for j in 0..ny {
                let col = col_offset + j;
                for i in 0..nx {
                    let field_idx = field_offset(column_count, i, col);
                    let alpha_idx = alpha_offset(nx, i, col);
                    let u = velocity_field[alpha_idx];
                    let center = nuclei_data[field_idx];
                    let mut flux_x = 0.0;
                    let mut flux_y = 0.0;
                    let mut flux_z = 0.0;
                    let mut diffusion_x = 0.0;
                    let mut diffusion_y = 0.0;
                    let mut diffusion_z = 0.0;

                    // 1st-order upwind advection.
                    if u.x > 0.0 && i > 0 {
                        let neighbor = nuclei_data[field_offset(column_count, i - 1, col)];
                        flux_x = u.x * (center - neighbor) / dx;
                    } else if u.x < 0.0 && i + 1 < nx {
                        let neighbor = nuclei_data[field_offset(column_count, i + 1, col)];
                        flux_x = u.x * (neighbor - center) / dx;
                    }

                    if u.y > 0.0 {
                        if j > 0 {
                            let neighbor = nuclei_data[field_offset(column_count, i, col - 1)];
                            flux_y = u.y * (center - neighbor) / dy;
                        }
                    } else if u.y < 0.0 && j + 1 < ny {
                        let neighbor = nuclei_data[field_offset(column_count, i, col + 1)];
                        flux_y = u.y * (neighbor - center) / dy;
                    }

                    if u.z > 0.0 {
                        if k > 0 {
                            let neighbor = nuclei_data[field_offset(column_count, i, col - ny)];
                            flux_z = u.z * (center - neighbor) / dz;
                        }
                    } else if u.z < 0.0 && k + 1 < nz {
                        let neighbor = nuclei_data[field_offset(column_count, i, col + ny)];
                        flux_z = u.z * (neighbor - center) / dz;
                    }

                    if diffusion_coefficient > 0.0 {
                        if i > 0 {
                            diffusion_x += (nuclei_data[field_offset(column_count, i - 1, col)]
                                - center)
                                * inv_dx2;
                        }
                        if i + 1 < nx {
                            diffusion_x += (nuclei_data[field_offset(column_count, i + 1, col)]
                                - center)
                                * inv_dx2;
                        }

                        if j > 0 {
                            diffusion_y += (nuclei_data[field_offset(column_count, i, col - 1)]
                                - center)
                                * inv_dy2;
                        }
                        if j + 1 < ny {
                            diffusion_y += (nuclei_data[field_offset(column_count, i, col + 1)]
                                - center)
                                * inv_dy2;
                        }

                        if k > 0 {
                            diffusion_z += (nuclei_data[field_offset(column_count, i, col - ny)]
                                - center)
                                * inv_dz2;
                        }
                        if k + 1 < nz {
                            diffusion_z += (nuclei_data[field_offset(column_count, i, col + ny)]
                                - center)
                                * inv_dz2;
                        }

                        diffusion_x *= diffusion_coefficient;
                        diffusion_y *= diffusion_coefficient;
                        diffusion_z *= diffusion_coefficient;
                    }

                    let s_net = transport.calculate_net_reaction_rate(
                        center,
                        source_data[field_idx].max(0.0), // Only vapor generation acts as a nuclei source.
                    );

                    next_data[field_idx] += dt
                        * (-flux_x - flux_y - flux_z
                            + diffusion_x
                            + diffusion_y
                            + diffusion_z
                            + s_net);
                    next_data[field_idx] = next_data[field_idx].max(0.0);
                }
            }
        }

        std::mem::swap(nuclei, next_nuclei);
        Ok(())
    }

    fn update_bubble_dynamics(
        &mut self,
        dt: f64,
        velocity_field: &[Vector3<f64>],
        pressure_field: &CavitationField,
        density_field: &CavitationField,
    ) -> Result<()> {
        let nx = self.vof_solver.nx;
        let ny = self.vof_solver.ny;
        let nz = self.vof_solver.nz;

        if let (Some(bubble_solver), Some(radius_field)) =
            (&mut self.bubble_solver, &mut self.bubble_radius_field)
        {
            for i in 0..nx {
                for j in 0..ny {
                    for k in 0..nz {
                        let idx = self.vof_solver.index(i, j, k);
                        let col = j + k * ny;
                        let pressure = pressure_field[[i, col]];
                        let velocity = velocity_field[idx];
                        let density = density_field[[i, col]];

                        let radius = bubble_solver
                            .update_bubble(i, j, k, pressure, velocity, density, dt)?;
                        radius_field[[i, col]] = radius;
                    }
                }
            }
        }
        Ok(())
    }

    #[cfg(test)]
    fn calculate_cavitation_source(
        &mut self,
        dt: f64,
        velocity_field: &[Vector3<f64>],
        pressure_field: &CavitationField,
        density_field: &CavitationField,
    ) -> CavitationField {
        let nx = self.vof_solver.nx;
        let ny = self.vof_solver.ny;
        let nz = self.vof_solver.nz;

        let mut cavitation_source = CavitationField::zeros([nx, ny * nz]);
        self.calculate_cavitation_source_into(
            dt,
            velocity_field,
            pressure_field,
            density_field,
            &mut cavitation_source,
        )
        .expect("validated source calculation");
        cavitation_source
    }

    /// # Theorem - Feasible Cavitation Source Interval
    ///
    /// For each cell, the relaxed cavitation source remains within the
    /// physically feasible increment interval:
    /// `-void_fraction <= s <= max_void_fraction - void_fraction`.
    ///
    /// **Proof sketch**: the source is first computed from the physical
    /// mass-transfer law and relaxation time, then clamped by the remaining
    /// void-room on the vaporization branch and by the available void on the
    /// condensation branch. The clamp is pointwise, so the interval holds for
    /// every cell independently.
    fn calculate_cavitation_source_into(
        &mut self,
        dt: f64,
        velocity_field: &[Vector3<f64>],
        pressure_field: &CavitationField,
        density_field: &CavitationField,
        cavitation_source: &mut CavitationField,
    ) -> Result<()> {
        let nx = self.vof_solver.nx;
        let ny = self.vof_solver.ny;
        let nz = self.vof_solver.nz;
        let column_count = ny * nz;
        let grid_size = nx * column_count;

        if !field_shape_matches(pressure_field, nx, column_count) {
            return Err(cfd_core::error::Error::DimensionMismatch {
                expected: grid_size,
                actual: field_len(pressure_field),
            });
        }
        if !field_shape_matches(density_field, nx, column_count) {
            return Err(cfd_core::error::Error::DimensionMismatch {
                expected: grid_size,
                actual: field_len(density_field),
            });
        }
        if !field_shape_matches(cavitation_source, nx, column_count) {
            return Err(cfd_core::error::Error::DimensionMismatch {
                expected: grid_size,
                actual: field_len(cavitation_source),
            });
        }

        let pressure_data = dense_field_slice(pressure_field);
        let density_data = dense_field_slice(density_field);
        let source_data = dense_field_slice_mut(cavitation_source);
        let alpha = self.vof_solver.alpha.as_slice();
        let inception_threshold = self.config.inception_threshold;
        let max_void_fraction = self.config.max_void_fraction;
        let relaxation_time = self.config.relaxation_time;
        let vapor_density = self.config.vapor_density;
        let vapor_pressure_base = self.config.vapor_pressure;
        let nuclei_field = self.nuclei_field.as_ref();

        source_data.fill(0.0);

        for k in 0..nz {
            let col_offset = k * ny;
            for j in 0..ny {
                let col = col_offset + j;

                for i in 0..nx {
                    let field_idx = field_offset(column_count, i, col);
                    let alpha_idx = alpha_offset(nx, i, col);
                    let pressure = pressure_data[field_idx];
                    let void_fraction = alpha[alpha_idx].clamp(0.0, 1.0);
                    let density_liquid = density_data[field_idx];
                    let mut vapor_pressure = vapor_pressure_base;

                    // Cascade effect: modify vapor pressure assumption if pre-existing nuclei are advected.
                    if let Some(nuclei) = nuclei_field {
                        vapor_pressure = cfd_core::physics::cavitation::nuclei_transport::nuclei_adjusted_vapor_pressure(
                            vapor_pressure,
                            nuclei[[i, col]],
                        );
                    }

                    // Calculate cavitation number (relative to vapor pressure).
                    // Zero dynamic pressure maps to infinite cavitation number.
                    let ref_vel = velocity_field[alpha_idx].norm();
                    let cavitation_number =
                        Self::cavitation_number(pressure, vapor_pressure, density_liquid, ref_vel);

                    if pressure < vapor_pressure
                        || (cavitation_number.is_finite()
                            && cavitation_number < inception_threshold)
                    {
                        self.inception_field[[i, col]] = INCEPTION_ACTIVE;
                    }

                    let mass_transfer = self.cavitation_model.mass_transfer_rate(
                        pressure,
                        vapor_pressure,
                        void_fraction,
                        density_liquid,
                        vapor_density,
                    );

                    let alpha_source = mass_transfer / vapor_density;
                    let relaxed_rate = alpha_source * dt / (dt + relaxation_time);

                    source_data[field_idx] = relaxed_rate;

                    // Limit void fraction increment.
                    if source_data[field_idx] > 0.0 {
                        let room = (max_void_fraction - void_fraction).max(0.0);
                        if source_data[field_idx] > room {
                            source_data[field_idx] = room;
                        }
                    } else if source_data[field_idx] < 0.0 {
                        let available = void_fraction.max(0.0);
                        if source_data[field_idx].abs() > available {
                            source_data[field_idx] = -available;
                        }
                    }
                }
            }
        }

        Ok(())
    }

    fn validate_flow_field_dimensions(
        &self,
        velocity_field: &[Vector3<f64>],
        pressure_field: &CavitationField,
        density_field: &CavitationField,
    ) -> Result<()> {
        let nx = self.vof_solver.nx;
        let ny = self.vof_solver.ny;
        let nz = self.vof_solver.nz;

        if !field_shape_matches(pressure_field, nx, ny * nz) {
            return Err(cfd_core::error::Error::DimensionMismatch {
                expected: nx * ny * nz,
                actual: field_len(pressure_field),
            });
        }
        if !field_shape_matches(density_field, nx, ny * nz) {
            return Err(cfd_core::error::Error::DimensionMismatch {
                expected: nx * ny * nz,
                actual: field_len(density_field),
            });
        }
        if velocity_field.len() != nx * ny * nz {
            return Err(cfd_core::error::Error::DimensionMismatch {
                expected: nx * ny * nz,
                actual: velocity_field.len(),
            });
        }

        Ok(())
    }

    #[inline]
    fn bubble_population_weight(&self) -> f64 {
        self.config
            .bubble_dynamics
            .as_ref()
            .map_or(1.0, |bubble_cfg| {
                bubble_cfg.number_density
                    * self.vof_solver.dx
                    * self.vof_solver.dy
                    * self.vof_solver.dz
            })
    }

    #[inline]
    fn cavitation_number(
        pressure: f64,
        vapor_pressure: f64,
        density_liquid: f64,
        velocity_magnitude: f64,
    ) -> f64 {
        let denom = 0.5 * density_liquid * velocity_magnitude * velocity_magnitude;
        if denom > 0.0 {
            (pressure - vapor_pressure) / denom
        } else {
            f64::INFINITY
        }
    }

    fn update_damage(
        &mut self,
        dt: f64,
        pressure_field: &CavitationField,
        density_field: &CavitationField,
    ) -> Result<()> {
        let nx = self.vof_solver.nx;
        let ny = self.vof_solver.ny;
        let nz = self.vof_solver.nz;
        let column_count = ny * nz;
        let grid_size = nx * column_count;

        if !field_shape_matches(pressure_field, nx, column_count) {
            return Err(cfd_core::error::Error::DimensionMismatch {
                expected: grid_size,
                actual: field_len(pressure_field),
            });
        }
        if !field_shape_matches(density_field, nx, column_count) {
            return Err(cfd_core::error::Error::DimensionMismatch {
                expected: grid_size,
                actual: field_len(density_field),
            });
        }

        let (Some(damage_model), Some(damage_field), Some(bubble_solver)) = (
            &self.damage_model,
            &mut self.damage_field,
            &self.bubble_solver,
        ) else {
            return Ok(());
        };

        if !field_shape_matches(damage_field, nx, column_count) {
            return Err(cfd_core::error::Error::DimensionMismatch {
                expected: grid_size,
                actual: field_len(damage_field),
            });
        }
        let bubble_population = bubble_solver.population_weight();
        let pressure_data = dense_field_slice(pressure_field);
        let density_data = dense_field_slice(density_field);
        let damage_data = dense_field_slice_mut(damage_field);
        let alpha = self.vof_solver.alpha.as_slice();

        for k in 0..nz {
            let col_offset = k * ny;
            for j in 0..ny {
                let col = col_offset + j;

                for i in 0..nx {
                    let field_idx = field_offset(column_count, i, col);
                    let alpha_idx = alpha_offset(nx, i, col);
                    let void_fraction = alpha[alpha_idx];

                    if void_fraction > 0.0 {
                        let pressure = pressure_data[field_idx];
                        let impact_pressure = bubble_solver.collapse_pressure(
                            i,
                            j,
                            k,
                            density_data[field_idx],
                            self.config.sound_speed,
                        );
                        let impact_frequency = bubble_solver
                            .get_bubble_frequency(i, j, k, pressure)
                            * bubble_population;

                        let erosion_rate = damage_model.mdpr(impact_pressure, impact_frequency, dt);
                        damage_data[field_idx] += erosion_rate * void_fraction;
                    }
                }
            }
        }

        Ok(())
    }

    /// Get current volume fraction field
    pub fn volume_fraction(&self) -> CavitationField {
        let nx = self.vof_solver.nx;
        let ny = self.vof_solver.ny;
        let nz = self.vof_solver.nz;
        let column_count = ny * nz;
        let mut volume_fraction = CavitationField::zeros([nx, column_count]);
        let data = dense_field_slice_mut(&mut volume_fraction);
        for col in 0..column_count {
            for i in 0..nx {
                data[field_offset(column_count, i, col)] =
                    self.vof_solver.alpha[alpha_offset(nx, i, col)];
            }
        }
        volume_fraction
    }

    /// Set current volume fraction field
    pub fn set_volume_fraction(&mut self, volume_fraction: &CavitationField) -> Result<()> {
        let nx = self.vof_solver.nx;
        let ny = self.vof_solver.ny;
        let nz = self.vof_solver.nz;
        let grid_size = nx * ny * nz;

        if !field_shape_matches(volume_fraction, nx, ny * nz) {
            return Err(cfd_core::error::Error::DimensionMismatch {
                expected: grid_size,
                actual: field_len(volume_fraction),
            });
        }

        let column_count = ny * nz;
        let data = dense_field_slice(volume_fraction);
        for col in 0..column_count {
            for i in 0..nx {
                self.vof_solver.alpha[alpha_offset(nx, i, col)] =
                    data[field_offset(column_count, i, col)];
            }
        }
        Ok(())
    }

    /// Get cavitation inception field
    pub fn inception_field(&self) -> &CavitationField {
        &self.inception_field
    }

    /// Get accumulated damage field
    pub fn damage_field(&self) -> Option<&CavitationField> {
        self.damage_field.as_ref()
    }

    /// Get bubble radius field
    pub fn bubble_radius_field(&self) -> Option<&CavitationField> {
        self.bubble_radius_field.as_ref()
    }

    /// Calculate sonoluminescence energy field
    ///
    /// Estimates the radiated energy during bubble collapse using the
    /// Rayleigh-Plesset adiabatic compression model and Stefan-Boltzmann law.
    pub fn sonoluminescence_energy_field(
        &self,
        pressure_field: &CavitationField,
        ambient_temperature: f64,
        flash_duration: f64,
        emissivity: f64,
    ) -> Result<CavitationField> {
        let nx = self.vof_solver.nx;
        let ny = self.vof_solver.ny;
        let nz = self.vof_solver.nz;

        if !field_shape_matches(pressure_field, nx, ny * nz) {
            return Err(cfd_core::error::Error::DimensionMismatch {
                expected: nx * ny * nz,
                actual: field_len(pressure_field),
            });
        }
        let Some(radius_field) = &self.bubble_radius_field else {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "Bubble radius field not available; enable bubble dynamics".to_string(),
            ));
        };
        let Some(bubble_cfg) = &self.config.bubble_dynamics else {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "Bubble dynamics config not available".to_string(),
            ));
        };

        let mut energy = CavitationField::zeros([nx, ny * nz]);

        let rp = RayleighPlesset::<f64> {
            initial_radius: bubble_cfg.initial_radius,
            liquid_density: self.config.liquid_density,
            liquid_viscosity: self.config.liquid_blood_model.viscosity(0.0),
            surface_tension: bubble_cfg.surface_tension,
            vapor_pressure: self.config.vapor_pressure,
            polytropic_index: bubble_cfg.polytropic_exponent,
        };
        let bubble_population = self.bubble_population_weight();

        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let col = j + k * ny;
                    let r_collapse = radius_field[[i, col]];
                    if r_collapse <= 0.0 || r_collapse >= bubble_cfg.initial_radius {
                        continue;
                    }
                    let estimate = rp.estimate_sonoluminescence(
                        pressure_field[[i, col]],
                        ambient_temperature,
                        r_collapse,
                        emissivity,
                        flash_duration,
                    )?;
                    energy[[i, col]] = estimate.radiated_energy * bubble_population;
                }
            }
        }

        Ok(energy)
    }

    /// Calculate cavitation statistics
    pub fn cavitation_statistics(&self) -> CavitationStatistics {
        let total_cells = field_len(&self.inception_field);
        let cavitating_cells = self
            .inception_field
            .iter()
            .filter(|&&x| x > INCEPTION_COUNT_THRESHOLD)
            .count();
        let total_void_fraction: f64 = self.vof_solver.alpha.iter().sum();
        let max_void_fraction = self.vof_solver.alpha.iter().copied().fold(0.0, f64::max);

        let max_damage = self
            .damage_field
            .as_ref()
            .map_or(0.0, |field| field.iter().copied().fold(0.0, f64::max));

        CavitationStatistics {
            cavitation_fraction: cavitating_cells as f64 / total_cells as f64,
            total_void_fraction,
            max_void_fraction,
            max_damage,
            cavitating_cells,
            total_cells,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vof::config::VofConfig;
    use crate::vof::BubbleDynamicsConfig;
    use cfd_core::physics::cavitation::damage::CavitationDamage;
    use cfd_core::physics::cavitation::models::CavitationModel;
    use cfd_core::physics::fluid::BloodModel;

    /// Helper: build a minimal CavitationVofConfig with Kunz model.
    fn make_config() -> CavitationVofConfig {
        CavitationVofConfig {
            vof_config: VofConfig::default(),
            cavitation_model: CavitationModel::Kunz {
                vaporization_coeff: 100.0,
                condensation_coeff: 100.0,
            },
            damage_model: None,
            bubble_dynamics: None,
            inception_threshold: 1.0,
            max_void_fraction: 0.9,
            relaxation_time: 1e-4,
            vapor_pressure: 2300.0, // ~water at 20 C
            liquid_density: 1000.0,
            liquid_blood_model: BloodModel::Newtonian(1e-3),
            vapor_density: 0.017,
            sound_speed: 1500.0,
            nuclei_transport: None,
        }
    }

    #[test]
    fn cavitation_source_positive_when_pressure_below_vapor() {
        let config = make_config();
        let nx = 2;
        let ny = 2;
        let nz = 2;
        let grid = nx * ny * nz;
        let mut solver =
            CavitationVofSolver::new(nx, ny, nz, config.clone()).expect("solver construction");

        let dt = 1e-5;
        // Pressure well below vapor pressure -> vaporization expected.
        let pressure_field = CavitationField::from_elem([nx, ny * nz], 500.0); // 500 Pa << 2300 Pa
        let density_field = CavitationField::from_elem([nx, ny * nz], 1000.0);
        let velocity_field = vec![Vector3::new(1.0, 0.0, 0.0); grid];

        let source = solver.calculate_cavitation_source(
            dt,
            &velocity_field,
            &pressure_field,
            &density_field,
        );

        // At least one cell should have a positive (vaporization) source.
        let max_source = source.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        assert!(
            max_source > 0.0,
            "cavitation source should be positive when p < p_vapor, got max={}",
            max_source
        );
    }

    #[test]
    fn cavitation_number_is_infinite_without_dynamic_pressure() {
        let number = CavitationVofSolver::cavitation_number(1.0e5, 2.3e3, 1.0e3, 0.0);
        assert!(number.is_infinite() && number.is_sign_positive());
    }

    #[test]
    fn condensation_when_pressure_above_vapor() {
        let config = make_config();
        let nx = 2;
        let ny = 2;
        let nz = 2;
        let grid = nx * ny * nz;
        let mut solver =
            CavitationVofSolver::new(nx, ny, nz, config.clone()).expect("solver construction");

        // Pre-set some void fraction so condensation can occur.
        for a in solver.vof_solver.alpha.iter_mut() {
            *a = 0.5;
        }

        let dt = 1e-5;
        let pressure_field = CavitationField::from_elem([nx, ny * nz], 100_000.0); // well above vapor
        let density_field = CavitationField::from_elem([nx, ny * nz], 1000.0);
        let velocity_field = vec![Vector3::new(1.0, 0.0, 0.0); grid];

        let source = solver.calculate_cavitation_source(
            dt,
            &velocity_field,
            &pressure_field,
            &density_field,
        );

        // Condensation yields negative source (void fraction should decrease).
        let min_source = source.iter().copied().fold(f64::INFINITY, f64::min);
        assert!(
            min_source <= 0.0,
            "source should be non-positive (condensation) when p >> p_vapor, got min={}",
            min_source
        );
    }

    #[test]
    fn cavitation_solver_rejects_damage_without_bubble_dynamics() {
        let mut config = make_config();
        config.damage_model = Some(CavitationDamage {
            yield_strength: 1.0e6,
            ultimate_strength: 2.0e6,
            hardness: 3.0e6,
            fatigue_strength: 4.0e5,
            cycles: 1,
        });

        let err = CavitationVofSolver::new(2, 2, 2, config)
            .err()
            .expect("damage modeling must require explicit bubble dynamics");
        assert!(err.to_string().contains("bubble dynamics"));
    }

    #[test]
    fn nuclei_diffusion_spreads_without_flow_or_source() {
        let mut config = make_config();
        config.nuclei_transport = Some(
            cfd_core::physics::cavitation::nuclei_transport::NucleiTransportConfig {
                dissolution_time_s: f64::INFINITY,
                generation_rate_factor: 0.0,
                diffusion_coefficient: 1.0e-4,
            },
        );

        let mut solver = CavitationVofSolver::new(3, 3, 1, config).expect("solver construction");
        if let Some(nuclei) = solver.nuclei_field.as_mut() {
            nuclei.fill(0.0);
            nuclei[[1, 1]] = 1.0;
        }

        let dt = 1.0e-3;
        let pressure_field = CavitationField::from_elem([3, 3], 100_000.0);
        let density_field = CavitationField::from_elem([3, 3], 1000.0);
        let velocity_field = vec![Vector3::zeros(); 9];

        solver
            .step(dt, &velocity_field, &pressure_field, &density_field)
            .expect("step should remain stable under diffusion-only transport");

        let nuclei = solver.nuclei_field.as_ref().expect("nuclei field");
        let center = nuclei[[1, 1]];
        let neighbor_sum = nuclei[[0, 1]] + nuclei[[2, 1]] + nuclei[[1, 0]] + nuclei[[1, 2]];

        assert!(center < 1.0, "diffusion must decrease the central peak");
        assert!(
            neighbor_sum > 0.0,
            "diffusion must redistribute nuclei into neighboring cells"
        );
    }

    #[test]
    fn nuclei_generation_follows_positive_cavitation_source() {
        let mut config = make_config();
        config.nuclei_transport = Some(
            cfd_core::physics::cavitation::nuclei_transport::NucleiTransportConfig {
                dissolution_time_s: f64::INFINITY,
                generation_rate_factor: 1.0,
                diffusion_coefficient: 0.0,
            },
        );

        let mut solver = CavitationVofSolver::new(2, 1, 1, config).expect("solver construction");
        if let Some(nuclei) = solver.nuclei_field.as_mut() {
            nuclei.fill(0.0);
        }

        let dt = 1.0e-3;
        let pressure_field = CavitationField::from_shape_vec([2, 1], vec![500.0, 100_000.0])
            .expect("valid pressure field shape");
        let density_field = CavitationField::from_shape_vec([2, 1], vec![1000.0, 1000.0])
            .expect("valid density field shape");
        let velocity_field = vec![Vector3::zeros(); 2];

        solver
            .step(dt, &velocity_field, &pressure_field, &density_field)
            .expect("step should remain stable under source-driven transport");

        let nuclei = solver.nuclei_field.as_ref().expect("nuclei field");
        assert!(
            nuclei[[0, 0]] > 0.0,
            "active cavitation must generate nuclei in the vaporizing cell"
        );
        assert_eq!(
            nuclei[[1, 0]],
            0.0,
            "inactive cells must not accumulate nuclei from a nonpositive source"
        );
    }

    #[test]
    fn cavitation_source_workspace_is_reused_across_steps() {
        let config = make_config();
        let mut solver = CavitationVofSolver::new(2, 2, 1, config).expect("solver construction");
        let expected_len = field_len(&solver.cavitation_source_workspace);

        let pressure_field = CavitationField::from_elem([2, 2], 500.0);
        let density_field = CavitationField::from_elem([2, 2], 1000.0);
        let velocity_field = vec![Vector3::new(1.0, 0.0, 0.0); 4];

        solver
            .step(1.0e-5, &velocity_field, &pressure_field, &density_field)
            .expect("first cavitation step");
        assert_eq!(field_len(&solver.cavitation_source_workspace), expected_len);

        solver
            .step(1.0e-5, &velocity_field, &pressure_field, &density_field)
            .expect("second cavitation step");
        assert_eq!(field_len(&solver.cavitation_source_workspace), expected_len);
    }

    #[test]
    fn cavitation_source_clamps_to_feasible_bounds_at_extremes() {
        let config = make_config();
        let mut solver = CavitationVofSolver::new(2, 1, 1, config).expect("solver construction");

        solver.vof_solver.alpha[0] = solver.config.max_void_fraction;
        solver.vof_solver.alpha[1] = 0.0;

        let mut pressure_field = CavitationField::from_elem([2, 1], 100_000.0);
        pressure_field[[0, 0]] = 500.0;
        let density_field = CavitationField::from_elem([2, 1], 1000.0);
        let velocity_field = vec![Vector3::new(1.0, 0.0, 0.0); 2];

        let source = solver.calculate_cavitation_source(
            1.0e-5,
            &velocity_field,
            &pressure_field,
            &density_field,
        );

        assert_eq!(source[[0, 0]], 0.0);
        assert_eq!(source[[1, 0]], 0.0);
    }

    #[test]
    fn step_rejects_flow_field_dimension_mismatch_without_bubble_dynamics() {
        let config = make_config();
        let mut solver = CavitationVofSolver::new(2, 2, 1, config).expect("solver construction");

        let pressure_field = CavitationField::from_elem([2, 2], 500.0);
        let density_field = CavitationField::from_elem([2, 2], 1000.0);
        let velocity_field = vec![Vector3::new(1.0, 0.0, 0.0); 3];

        let err = solver
            .step(1.0e-5, &velocity_field, &pressure_field, &density_field)
            .expect_err("dimension mismatch should be rejected before cavitation update");

        assert!(matches!(
            err,
            cfd_core::error::Error::DimensionMismatch {
                expected: 4,
                actual: 3
            }
        ));
    }

    #[test]
    fn sonoluminescence_scales_with_number_density() {
        let mut low_density = make_config();
        low_density.bubble_dynamics = Some(BubbleDynamicsConfig {
            initial_radius: 2.0e-6,
            number_density: 1.0e6,
            polytropic_exponent: 1.4,
            surface_tension: 0.072,
        });

        let mut high_density = low_density.clone();
        high_density
            .bubble_dynamics
            .as_mut()
            .expect("bubble dynamics should be configured")
            .number_density = 2.0e6;

        let mut solver_low = CavitationVofSolver::new(1, 1, 1, low_density).unwrap();
        let mut solver_high = CavitationVofSolver::new(1, 1, 1, high_density).unwrap();

        if let Some(radius_field) = solver_low.bubble_radius_field.as_mut() {
            radius_field[[0, 0]] = 5.0e-7;
        }
        if let Some(radius_field) = solver_high.bubble_radius_field.as_mut() {
            radius_field[[0, 0]] = 5.0e-7;
        }

        let pressure_field = CavitationField::from_elem([1, 1], 1.0e5);
        let energy_low = solver_low
            .sonoluminescence_energy_field(&pressure_field, 293.15, 5.0e-11, 1.0)
            .unwrap();
        let energy_high = solver_high
            .sonoluminescence_energy_field(&pressure_field, 293.15, 5.0e-11, 1.0)
            .unwrap();

        let ratio = energy_high[[0, 0]] / energy_low[[0, 0]];
        assert!((ratio - 2.0).abs() < 1e-12 * 2.0);
    }

    #[test]
    fn volume_fraction_stays_clamped_after_step() {
        let config = make_config();
        let nx = 2;
        let ny = 2;
        let nz = 2;
        let grid = nx * ny * nz;
        let mut solver =
            CavitationVofSolver::new(nx, ny, nz, config.clone()).expect("solver construction");

        let dt = 1e-5;
        // Very low pressure to drive maximum vaporization.
        let pressure_field = CavitationField::from_elem([nx, ny * nz], 100.0);
        let density_field = CavitationField::from_elem([nx, ny * nz], 1000.0);
        let velocity_field = vec![Vector3::new(0.1, 0.0, 0.0); grid];

        // Run several steps.
        for _ in 0..20 {
            solver
                .step(dt, &velocity_field, &pressure_field, &density_field)
                .expect("step should succeed");
        }

        let vf = solver.volume_fraction();
        for val in vf.iter() {
            assert!(
                *val >= 0.0 && *val <= 1.0,
                "volume fraction {} out of [0,1] after stepping",
                val
            );
        }
    }

    #[test]
    fn damage_accumulates_below_previous_void_fraction_floor() {
        let mut config = make_config();
        config.damage_model = Some(CavitationDamage {
            yield_strength: 1.0e6,
            ultimate_strength: 2.0e6,
            hardness: 3.0e6,
            fatigue_strength: 4.0e5,
            cycles: 1,
        });
        config.bubble_dynamics = Some(BubbleDynamicsConfig {
            initial_radius: 2.0e-6,
            number_density: 1.0e12,
            polytropic_exponent: 1.4,
            surface_tension: 0.072,
        });

        let nx = 2;
        let ny = 1;
        let nz = 1;
        let mut solver = CavitationVofSolver::new(nx, ny, nz, config).expect("solver");

        solver.vof_solver.alpha[0] = 5.0e-3;
        solver.vof_solver.alpha[1] = 0.0;

        // Keep the pressure just above vapor pressure so the bubble frequency
        // remains positive while the low-void cell still survives the source update.
        let pressure_field = CavitationField::from_elem([nx, ny * nz], 2_300.1);
        let density_field = CavitationField::from_elem([nx, ny * nz], 1000.0);
        let velocity_field = vec![Vector3::new(0.1, 0.0, 0.0); nx * ny * nz];

        solver
            .step(1.0e-5, &velocity_field, &pressure_field, &density_field)
            .expect("step");

        let damage = solver.damage_field().expect("damage field");
        assert!(
            damage[[0, 0]] > 0.0,
            "damage must accumulate for sub-1% void fraction cells"
        );
        assert_eq!(
            damage[[1, 0]],
            0.0,
            "pure liquid cells must not accumulate cavitation damage"
        );
    }

    #[test]
    fn update_damage_rejects_pressure_dimension_mismatch() {
        let mut config = make_config();
        config.damage_model = Some(CavitationDamage {
            yield_strength: 1.0e6,
            ultimate_strength: 2.0e6,
            hardness: 3.0e6,
            fatigue_strength: 4.0e5,
            cycles: 1,
        });
        config.bubble_dynamics = Some(BubbleDynamicsConfig {
            initial_radius: 2.0e-6,
            number_density: 1.0e12,
            polytropic_exponent: 1.4,
            surface_tension: 0.072,
        });

        let mut solver = CavitationVofSolver::new(2, 1, 1, config).expect("solver");
        let pressure_field = CavitationField::from_elem([1, 1], 2_300.1);
        let density_field = CavitationField::from_elem([2, 1], 1000.0);

        let err = solver
            .update_damage(1.0e-5, &pressure_field, &density_field)
            .expect_err("pressure dimension mismatch should be rejected");

        assert!(matches!(
            err,
            cfd_core::error::Error::DimensionMismatch {
                expected: 2,
                actual: 1
            }
        ));
    }
}
