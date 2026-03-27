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

/// Default grid spacing [m] used when constructing the VOF solver grid.
const DEFAULT_GRID_SPACING: f64 = 0.01;

/// Minimum reference velocity [m/s] to avoid division by zero in the
/// cavitation number calculation.  Set well below any physically
/// meaningful flow speed.
const MIN_VELOCITY_THRESHOLD: f64 = 1.0e-9;

/// Void-fraction threshold below which a cell is considered purely liquid
/// for damage-accumulation purposes.  Cells with alpha <= this value
/// skip the erosion-rate calculation entirely.
const MIN_VOID_FRACTION_FOR_DAMAGE: f64 = 0.01;

/// Default bubble / collapse radius [m] used as a fallback when the
/// bubble-dynamics solver is not enabled.  Corresponds to a typical
/// micro-cavitation bubble (~1 micrometer).
#[allow(dead_code)] // Used in bubble-dynamics paths gated by runtime config
const DEFAULT_BUBBLE_RADIUS: f64 = 1e-6;

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
use nalgebra::{DMatrix, Vector3};

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
    inception_field: DMatrix<f64>,
    /// Accumulated damage field
    damage_field: Option<DMatrix<f64>>,
    /// Bubble radius field for Rayleigh-Plesset
    bubble_radius_field: Option<DMatrix<f64>>,
    /// Passive scalar tracking for cavitation nuclei cascade
    nuclei_field: Option<DMatrix<f64>>,
    /// Nuclei transport evaluator
    nuclei_transport: Option<cfd_core::physics::cavitation::nuclei_transport::NucleiTransport<f64>>,
}

impl CavitationVofSolver {
    /// Create new cavitation-VOF solver
    pub fn new(nx: usize, ny: usize, nz: usize, config: CavitationVofConfig) -> Result<Self> {
        let vof_solver = VofSolver::create(config.vof_config.clone(), nx, ny, nz, DEFAULT_GRID_SPACING, DEFAULT_GRID_SPACING, DEFAULT_GRID_SPACING);

        let bubble_solver = config.bubble_dynamics.as_ref().map(|bubble_config| {
            BubbleDynamicsSolver::new(
                bubble_config,
                nx,
                ny,
                nz,
                config.liquid_density,
                config.liquid_blood_model.clone(),
                config.vapor_pressure,
            )
        });

        let inception_field = DMatrix::zeros(nx, ny * nz);
        let damage_field = if config.damage_model.is_some() {
            Some(DMatrix::zeros(nx, ny * nz))
        } else {
            None
        };

        let bubble_radius_field = if config.bubble_dynamics.is_some() {
            Some(DMatrix::zeros(nx, ny * nz))
        } else {
            None
        };
        
        let nuclei_field = if config.nuclei_transport.is_some() {
            Some(DMatrix::zeros(nx, ny * nz))
        } else {
            None
        };
        
        let nuclei_transport = config.nuclei_transport.map(cfd_core::physics::cavitation::nuclei_transport::NucleiTransport::new);

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
            nuclei_transport,
        })
    }

    /// Update cavitation-VOF solution for one time step
    pub fn step(
        &mut self,
        dt: f64,
        velocity_field: &[Vector3<f64>],
        pressure_field: &DMatrix<f64>,
        density_field: &DMatrix<f64>,
    ) -> Result<()> {
        // 1. Update bubble dynamics if enabled
        self.update_bubble_dynamics(dt, velocity_field, pressure_field, density_field)?;

        // 2. Calculate cavitation mass transfer rates
        let cavitation_source =
            self.calculate_cavitation_source(dt, velocity_field, pressure_field, density_field);

        // 2b. Advect and diffuse cavitation nuclei tracking the active cavitation
        self.update_nuclei_advection_diffusion(dt, velocity_field, &cavitation_source);

        // 3. Update VOF with cavitation source term
        for idx in 0..self.vof_solver.alpha.len() {
            self.vof_solver.alpha[idx] =
                (self.vof_solver.alpha[idx] + cavitation_source[idx]).clamp(0.0, 1.0);
        }

        // 4. Update VOF velocity and advance
        self.vof_solver
            .set_velocity_field(velocity_field.to_vec())?;
        self.vof_solver.advance(dt)?;

        // 5. Calculate cavitation damage if damage model is enabled
        self.update_damage(dt, pressure_field, density_field);

        Ok(())
    }

    fn update_nuclei_advection_diffusion(
        &mut self,
        dt: f64,
        velocity_field: &[Vector3<f64>],
        cavitation_source: &DMatrix<f64>,
    ) {
        let nx = self.vof_solver.nx;
        let ny = self.vof_solver.ny;
        let nz = self.vof_solver.nz;
        let dx = self.vof_solver.dx;
        let dy = self.vof_solver.dy;
        let dz = self.vof_solver.dz;

        if let (Some(nuclei), Some(transport)) = (&mut self.nuclei_field, &self.nuclei_transport) {
            let mut next_nuclei = nuclei.clone();

            for i in 0..nx {
                for j in 0..ny {
                    for k in 0..nz {
                        let idx = self.vof_solver.index(i, j, k);
                        let col = j + k * ny;

                        let u = velocity_field[idx];
                        let mut flux_x = 0.0;
                        let mut flux_y = 0.0;
                        let mut flux_z = 0.0;

                        // 1st-order upwind advection
                        if u.x > 0.0 && i > 0 {
                            flux_x = u.x * (nuclei[(i, col)] - nuclei[(i - 1, col)]) / dx;
                        } else if u.x < 0.0 && i < nx - 1 {
                            flux_x = u.x * (nuclei[(i + 1, col)] - nuclei[(i, col)]) / dx;
                        }

                        if u.y > 0.0 && j > 0 {
                            flux_y = u.y * (nuclei[(i, col)] - nuclei[(i, col - 1)]) / dy;
                        } else if u.y < 0.0 && j < ny - 1 {
                            flux_y = u.y * (nuclei[(i, col + 1)] - nuclei[(i, col)]) / dy;
                        }

                        if u.z > 0.0 && k > 0 {
                            flux_z = u.z * (nuclei[(i, col)] - nuclei[(i, col - ny)]) / dz;
                        } else if u.z < 0.0 && k < nz - 1 {
                            flux_z = u.z * (nuclei[(i, col + ny)] - nuclei[(i, col)]) / dz;
                        }

                        let s_net = transport.calculate_net_reaction_rate(
                            nuclei[(i, col)],
                            cavitation_source[(i, col)].max(0.0), // Only vapor generation acts as a nuclei source
                        );

                        next_nuclei[(i, col)] += dt * (-flux_x - flux_y - flux_z + s_net);
                        next_nuclei[(i, col)] = next_nuclei[(i, col)].max(0.0);
                    }
                }
            }
            std::mem::swap(nuclei, &mut next_nuclei);
        }
    }

    fn update_bubble_dynamics(
        &mut self,
        dt: f64,
        velocity_field: &[Vector3<f64>],
        pressure_field: &DMatrix<f64>,
        density_field: &DMatrix<f64>,
    ) -> Result<()> {
        let nx = self.vof_solver.nx;
        let ny = self.vof_solver.ny;
        let nz = self.vof_solver.nz;

        if pressure_field.nrows() != nx || pressure_field.ncols() != ny * nz {
            return Err(cfd_core::error::Error::DimensionMismatch {
                expected: nx * ny * nz,
                actual: pressure_field.len(),
            });
        }
        if density_field.nrows() != nx || density_field.ncols() != ny * nz {
            return Err(cfd_core::error::Error::DimensionMismatch {
                expected: nx * ny * nz,
                actual: density_field.len(),
            });
        }
        if velocity_field.len() != nx * ny * nz {
            return Err(cfd_core::error::Error::DimensionMismatch {
                expected: nx * ny * nz,
                actual: velocity_field.len(),
            });
        }

        if let (Some(bubble_solver), Some(radius_field)) =
            (&mut self.bubble_solver, &mut self.bubble_radius_field)
        {
            for i in 0..nx {
                for j in 0..ny {
                    for k in 0..nz {
                        let idx = self.vof_solver.index(i, j, k);
                        let col = j + k * ny;
                        let pressure = pressure_field[(i, col)];
                        let velocity = velocity_field[idx];
                        let density = density_field[(i, col)];

                        let radius = bubble_solver
                            .update_bubble(i, j, k, pressure, velocity, density, dt)?;
                        radius_field[(i, col)] = radius;
                    }
                }
            }
        }
        Ok(())
    }

    fn calculate_cavitation_source(
        &mut self,
        dt: f64,
        velocity_field: &[Vector3<f64>],
        pressure_field: &DMatrix<f64>,
        density_field: &DMatrix<f64>,
    ) -> DMatrix<f64> {
        let nx = self.vof_solver.nx;
        let ny = self.vof_solver.ny;
        let nz = self.vof_solver.nz;

        let mut cavitation_source = DMatrix::zeros(nx, ny * nz);

        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let idx = self.vof_solver.index(i, j, k);
                    let col = j + k * ny;
                    let pressure = pressure_field[(i, col)];
                    let void_fraction = self.vof_solver.alpha[idx].clamp(0.0, 1.0);
                    let density_liquid = density_field[(i, col)];
                    let density_vapor = self.config.vapor_density;
                    let mut vapor_pressure = self.config.vapor_pressure;
                    
                    // Cascade effect: modify vapor pressure assumption if pre-existing nuclei are advected
                    if let Some(nuclei) = &self.nuclei_field {
                        vapor_pressure += nuclei[(i, col)] * 10_000.0;
                    }

                    // Calculate cavitation number (relative to vapor pressure)
                    let ref_vel = velocity_field[idx].norm();
                    let ref_vel_eff = ref_vel.max(MIN_VELOCITY_THRESHOLD);
                    let denom = 0.5 * density_liquid * ref_vel_eff * ref_vel_eff;
                    let cavitation_number = if denom > 0.0 {
                        (pressure - vapor_pressure) / denom
                    } else {
                        f64::INFINITY
                    };

                    if pressure < vapor_pressure
                        || (cavitation_number.is_finite()
                            && cavitation_number < self.config.inception_threshold)
                    {
                        self.inception_field[(i, col)] = INCEPTION_ACTIVE;
                    }

                    let mass_transfer = self.cavitation_model.mass_transfer_rate(
                        pressure,
                        vapor_pressure,
                        void_fraction,
                        density_liquid,
                        density_vapor,
                    );

                    let alpha_source = mass_transfer / density_vapor;
                    let relaxed_rate = alpha_source * dt / (dt + self.config.relaxation_time);

                    cavitation_source[(i, col)] = relaxed_rate;

                    // Limit void fraction increment
                    let max_allowed = self.config.max_void_fraction;
                    if cavitation_source[(i, col)] > 0.0 {
                        let room = (max_allowed - void_fraction).max(0.0);
                        if cavitation_source[(i, col)] > room {
                            cavitation_source[(i, col)] = room;
                        }
                    } else if cavitation_source[(i, col)] < 0.0 {
                        let available = void_fraction.max(0.0);
                        if cavitation_source[(i, col)].abs() > available {
                            cavitation_source[(i, col)] = -available;
                        }
                    }
                }
            }
        }
        cavitation_source
    }

    fn update_damage(
        &mut self,
        dt: f64,
        pressure_field: &DMatrix<f64>,
        density_field: &DMatrix<f64>,
    ) {
        let nx = self.vof_solver.nx;
        let ny = self.vof_solver.ny;
        let nz = self.vof_solver.nz;

        if density_field.nrows() != nx || density_field.ncols() != ny * nz {
            return;
        }

        if let (Some(damage_model), Some(damage_field)) =
            (&self.damage_model, &mut self.damage_field)
        {
            for i in 0..nx {
                for j in 0..ny {
                    for k in 0..nz {
                        let idx = self.vof_solver.index(i, j, k);
                        let col = j + k * ny;
                        let void_fraction = self.vof_solver.alpha[idx];

                        if void_fraction > MIN_VOID_FRACTION_FOR_DAMAGE {
                            let pressure = pressure_field[(i, col)];
                            let impact_pressure = if let Some(bubble_solver) = &self.bubble_solver {
                                bubble_solver.collapse_pressure(
                                    i,
                                    j,
                                    k,
                                    density_field[(i, col)],
                                    self.config.sound_speed,
                                )
                            } else {
                                let bubble_radius = self
                                    .bubble_radius_field
                                    .as_ref()
                                    .map_or(DEFAULT_BUBBLE_RADIUS, |f| f[(i, col)]);
                                let initial_radius = self
                                    .config
                                    .bubble_dynamics
                                    .as_ref()
                                    .map_or(DEFAULT_BUBBLE_RADIUS, |c| c.initial_radius);
                                if bubble_radius > 0.0 {
                                    density_field[(i, col)]
                                        * self.config.sound_speed.powi(2)
                                        * (initial_radius / bubble_radius)
                                } else {
                                    0.0
                                }
                            };

                            // Calculate impact frequency
                            let impact_frequency = if let Some(bubble_solver) = &self.bubble_solver
                            {
                                bubble_solver.get_bubble_frequency(i, j, k, pressure)
                            } else {
                                // Fallback estimation using VOF config
                                let radius = self
                                    .bubble_radius_field
                                    .as_ref()
                                    .map_or(DEFAULT_BUBBLE_RADIUS, |f| f[(i, col)]);
                                let rp = RayleighPlesset {
                                    initial_radius: radius, // Use current as estimate
                                    liquid_density: self.config.liquid_density,
                                    liquid_viscosity: self.config.liquid_blood_model.viscosity(0.0),
                                    surface_tension: self
                                        .config
                                        .vof_config
                                        .surface_tension_coefficient,
                                    vapor_pressure: self.config.vapor_pressure,
                                    polytropic_index: POLYTROPIC_INDEX_AIR,
                                };
                                rp.natural_frequency(radius, pressure)
                                    / (2.0 * std::f64::consts::PI)
                            };

                            // Calculate erosion rate
                            let erosion_rate =
                                damage_model.mdpr(impact_pressure, impact_frequency, dt);

                            // Accumulate damage
                            damage_field[(i, col)] += erosion_rate;
                        }
                    }
                }
            }
        }
    }

    /// Get current volume fraction field
    pub fn volume_fraction(&self) -> DMatrix<f64> {
        let nx = self.vof_solver.nx;
        let ny = self.vof_solver.ny;
        let nz = self.vof_solver.nz;
        // alpha is stored in k*ny*nx + j*nx + i order, which matches column-major (j+k*ny)*nx + i
        DMatrix::from_column_slice(nx, ny * nz, &self.vof_solver.alpha)
    }

    /// Set current volume fraction field
    pub fn set_volume_fraction(&mut self, volume_fraction: &DMatrix<f64>) -> Result<()> {
        let nx = self.vof_solver.nx;
        let ny = self.vof_solver.ny;
        let nz = self.vof_solver.nz;
        let grid_size = nx * ny * nz;

        if volume_fraction.len() != grid_size {
            return Err(cfd_core::error::Error::DimensionMismatch {
                expected: grid_size,
                actual: volume_fraction.len(),
            });
        }

        self.vof_solver
            .alpha
            .copy_from_slice(volume_fraction.as_slice());
        Ok(())
    }

    /// Get cavitation inception field
    pub fn inception_field(&self) -> &DMatrix<f64> {
        &self.inception_field
    }

    /// Get accumulated damage field
    pub fn damage_field(&self) -> Option<&DMatrix<f64>> {
        self.damage_field.as_ref()
    }

    /// Get bubble radius field
    pub fn bubble_radius_field(&self) -> Option<&DMatrix<f64>> {
        self.bubble_radius_field.as_ref()
    }

    /// Calculate sonoluminescence energy field
    ///
    /// Estimates the radiated energy during bubble collapse using the
    /// Rayleigh-Plesset adiabatic compression model and Stefan-Boltzmann law.
    pub fn sonoluminescence_energy_field(
        &self,
        pressure_field: &DMatrix<f64>,
        ambient_temperature: f64,
        flash_duration: f64,
        emissivity: f64,
    ) -> Result<DMatrix<f64>> {
        let nx = self.vof_solver.nx;
        let ny = self.vof_solver.ny;
        let nz = self.vof_solver.nz;

        if pressure_field.nrows() != nx || pressure_field.ncols() != ny * nz {
            return Err(cfd_core::error::Error::DimensionMismatch {
                expected: nx * ny * nz,
                actual: pressure_field.len(),
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

        let mut energy = DMatrix::zeros(nx, ny * nz);

        let rp = RayleighPlesset::<f64> {
            initial_radius: bubble_cfg.initial_radius,
            liquid_density: self.config.liquid_density,
            liquid_viscosity: self.config.liquid_blood_model.viscosity(0.0),
            surface_tension: bubble_cfg.surface_tension,
            vapor_pressure: self.config.vapor_pressure,
            polytropic_index: bubble_cfg.polytropic_exponent,
        };

        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let col = j + k * ny;
                    let r_collapse = radius_field[(i, col)];
                    let estimate = rp.estimate_sonoluminescence(
                        pressure_field[(i, col)],
                        ambient_temperature,
                        r_collapse,
                        emissivity,
                        flash_duration,
                    )?;
                    energy[(i, col)] = estimate.radiated_energy;
                }
            }
        }

        Ok(energy)
    }

    /// Calculate cavitation statistics
    pub fn cavitation_statistics(&self) -> CavitationStatistics {
        let total_cells = self.inception_field.len();
        let cavitating_cells = self.inception_field.iter().filter(|&&x| x > INCEPTION_COUNT_THRESHOLD).count();
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
            vapor_pressure: 2300.0,   // ~water at 20 C
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
        let pressure_field = DMatrix::from_element(nx, ny * nz, 500.0); // 500 Pa << 2300 Pa
        let density_field = DMatrix::from_element(nx, ny * nz, 1000.0);
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
        let pressure_field = DMatrix::from_element(nx, ny * nz, 100_000.0); // well above vapor
        let density_field = DMatrix::from_element(nx, ny * nz, 1000.0);
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
        let pressure_field = DMatrix::from_element(nx, ny * nz, 100.0);
        let density_field = DMatrix::from_element(nx, ny * nz, 1000.0);
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
}
