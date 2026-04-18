//! Bubble Dynamics Solver — Rayleigh-Plesset Integration
//!
//! This module implements the single-bubble Rayleigh-Plesset ordinary
//! differential equation on a per-cell basis for cavitation simulations.
//!
//! # Theorem — Rayleigh-Plesset Bubble Dynamics
//!
//! The radius $R(t)$ of a spherical cavitation bubble in an incompressible
//! liquid satisfies the Rayleigh-Plesset ODE:
//!
//! ```text
//! R R̈ + (3/2) Ṙ² = (1/ρ_l)[p_B(t) − p_∞(t) − 4μ Ṙ/R − 2σ/R]
//! ```
//!
//! where $p_B$ is the internal bubble pressure and $p_∞$ is the far-field
//! pressure.  The ODE is integrated with semi-implicit Euler to maintain
//! stability during violent collapse.
//!
//! **Proof sketch**: The ODE derives from the incompressible Euler equation
//! in spherical coordinates with the kinetic-energy integral evaluated from
//! $r = R$ to $r = \infty$.  Total mechanical energy is conserved up to
//! viscous dissipation and surface tension work.
//!
//! **Reference**: Plesset, M.S. (1949). J. Appl. Mech. 16:277–282.
//!
//! # Theorem — Non-Newtonian Apparent Viscosity at Bubble Wall
//!
//! The macroscopic shear rate $\dot{\gamma}$ at a spherical bubble wall expanding
//! or collapsing radially in a stationary fluid is analytically given by:
//!
//! ```text
//! \dot{\gamma}_{wall} = \sqrt{12} \left| \frac{\dot{R}}{R} \right|
//! ```
//!
//! **Proof sketch**: For purely radial flow $u_r = \dot{R} (R/r)^2$, the rate-of-strain
//! tensor components are non-zero only on the diagonal: $S_{rr} = -2\dot{R}R^2/r^3$ and
//! $S_{\theta\theta} = S_{\phi\phi} = \dot{R}R^2/r^3$. The second invariant
//! $S_{ij}S_{ij}$ evaluated at the wall $r=R$ is $3(\dot{R}/R)^2$. The shear rate
//! magnitude $\dot{\gamma} = \sqrt{2 S_{ij}S_{ij}} = \sqrt{12} |\dot{R}/R|$.
//!
//! This local shear rate dictates the apparent dynamic viscosity $\mu(\dot{\gamma})$
//! for non-Newtonian fluids like blood during the Rayleigh-Plesset integration.

use cfd_core::error::Result;
use cfd_core::physics::cavitation::rayleigh_plesset::RayleighPlesset;
use cfd_core::physics::fluid::BloodModel;
use nalgebra::Vector3;

/// Bubble dynamics configuration
#[derive(Debug, Clone)]
pub struct BubbleDynamicsConfig {
    /// Initial bubble radius (m)
    pub initial_radius: f64,
    /// Bubble number density (1/m³)
    ///
    /// This is converted to an expected bubble population per cell using the
    /// control-volume size so the per-cell damage and sonoluminescence outputs
    /// remain consistent with the modeled nuclei density.
    pub number_density: f64,
    /// Polytropic exponent for gas compression
    pub polytropic_exponent: f64,
    /// Surface tension coefficient (N/m)
    pub surface_tension: f64,
}

/// Bubble dynamics solver for Rayleigh-Plesset equation
pub struct BubbleDynamicsSolver {
    /// Grid dimensions.
    nx: usize,
    ny: usize,
    /// Bubble configurations per grid cell.
    configs: Vec<RayleighPlesset<f64>>,
    /// Current bubble radii per grid cell.
    radii: Vec<f64>,
    /// Current bubble wall velocities per grid cell.
    velocities: Vec<f64>,
    /// Expected bubble population represented by one coarse cell.
    bubble_population_weight: f64,
    /// Non-Newtonian blood model for apparent viscosity
    blood_model: BloodModel<f64>,
}

impl BubbleDynamicsSolver {
    /// Create new bubble dynamics solver
    pub fn new(
        config: &BubbleDynamicsConfig,
        nx: usize,
        ny: usize,
        nz: usize,
        dx: f64,
        dy: f64,
        dz: f64,
        liquid_density: f64,
        blood_model: BloodModel<f64>,
        vapor_pressure: f64,
    ) -> Self {
        let initial_viscosity = blood_model.viscosity(0.0);
        let rp = RayleighPlesset {
            initial_radius: config.initial_radius,
            liquid_density,
            liquid_viscosity: initial_viscosity,
            surface_tension: config.surface_tension,
            vapor_pressure,
            polytropic_index: config.polytropic_exponent,
        };
        let len = nx * ny * nz;
        let bubble_population_weight = config.number_density * dx * dy * dz;

        Self {
            nx,
            ny,
            configs: vec![rp; len],
            radii: vec![config.initial_radius; len],
            velocities: vec![0.0; len],
            bubble_population_weight,
            blood_model,
        }
    }

    #[inline]
    fn index(&self, i: usize, j: usize, k: usize) -> usize {
        k * self.ny * self.nx + j * self.nx + i
    }

    /// Expected bubble population represented by one coarse cell.
    #[must_use]
    pub(crate) fn population_weight(&self) -> f64 {
        self.bubble_population_weight
    }

    /// Update bubble dynamics for a specific cell
    pub fn update_bubble(
        &mut self,
        i: usize,
        j: usize,
        k: usize,
        pressure: f64,
        _velocity: Vector3<f64>,
        _density: f64,
        dt: f64,
    ) -> Result<f64> {
        let idx = self.index(i, j, k);
        let config = &mut self.configs[idx];
        let radius = self.radii[idx];
        let velocity = self.velocities[idx];

        // Calculate apparent viscosity at the bubble wall.
        // theorem: The shear rate at a spherical bubble wall expanding/collapsing radially
        // is given by γ̇_wall = √12 |Ṙ / R|.
        let shear_rate = 12.0_f64.sqrt() * (velocity.abs() / radius.max(1e-12));
        config.liquid_viscosity = self.blood_model.viscosity(shear_rate);

        // Calculate acceleration using Rayleigh-Plesset equation
        let acceleration = config.bubble_acceleration(radius, velocity, pressure);

        // Simple semi-implicit Euler integration
        let new_velocity = velocity + acceleration * dt;
        let new_radius = (radius + new_velocity * dt).max(1e-9);

        self.radii[idx] = new_radius;
        self.velocities[idx] = new_velocity;

        Ok(new_radius)
    }

    /// Get collapse pressure for damage calculation
    pub fn collapse_pressure(
        &self,
        i: usize,
        j: usize,
        k: usize,
        liquid_density: f64,
        sound_speed: f64,
    ) -> f64 {
        let idx = self.index(i, j, k);
        let config = &self.configs[idx];
        let radius = self.radii[idx];
        let initial_radius = config.initial_radius;

        if radius > 0.0 {
            // Impact pressure scales with (R_max / R_collapse)
            liquid_density * sound_speed * sound_speed * (initial_radius / radius)
        } else {
            0.0
        }
    }

    /// Get bubble natural frequency (Hz) for impact frequency estimation
    pub fn get_bubble_frequency(&self, i: usize, j: usize, k: usize, ambient_pressure: f64) -> f64 {
        let idx = self.index(i, j, k);
        let config = &self.configs[idx];
        let radius = self.radii[idx];

        let omega = config.natural_frequency(radius, ambient_pressure);
        omega / (2.0 * std::f64::consts::PI)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn collapse_pressure_and_frequency_use_initialized_bubble_state() {
        let config = BubbleDynamicsConfig {
            initial_radius: 2.0e-6,
            number_density: 1.0e12,
            polytropic_exponent: 1.4,
            surface_tension: 0.072,
        };
        let solver = BubbleDynamicsSolver::new(
            &config,
            1,
            1,
            1,
            1.0,
            1.0,
            1.0,
            1000.0,
            BloodModel::Newtonian(1.0e-3),
            2300.0,
        );

        let collapse_pressure = solver.collapse_pressure(0, 0, 0, 1000.0, 1500.0);
        let frequency = solver.get_bubble_frequency(0, 0, 0, 2500.0);

        assert!(collapse_pressure > 0.0);
        assert!(frequency > 0.0);
    }

    #[test]
    fn population_weight_reflects_cell_volume() {
        let config = BubbleDynamicsConfig {
            initial_radius: 2.0e-6,
            number_density: 2.5e12,
            polytropic_exponent: 1.4,
            surface_tension: 0.072,
        };
        let solver = BubbleDynamicsSolver::new(
            &config,
            2,
            2,
            2,
            0.01,
            0.02,
            0.03,
            1000.0,
            BloodModel::Newtonian(1.0e-3),
            2300.0,
        );

        let expected = config.number_density * 0.01 * 0.02 * 0.03;
        assert!((solver.population_weight() - expected).abs() < 1e-12 * expected.max(1.0));
    }
}
