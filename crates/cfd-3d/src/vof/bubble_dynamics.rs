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

use cfd_core::error::Result;
use cfd_core::physics::cavitation::rayleigh_plesset::RayleighPlesset;
use nalgebra::Vector3;
use std::collections::HashMap;

/// Bubble dynamics configuration
#[derive(Debug, Clone)]
pub struct BubbleDynamicsConfig {
    /// Initial bubble radius (m)
    pub initial_radius: f64,
    /// Bubble number density (1/m³)
    pub number_density: f64,
    /// Polytropic exponent for gas compression
    pub polytropic_exponent: f64,
    /// Surface tension coefficient (N/m)
    pub surface_tension: f64,
}

/// Bubble dynamics solver for Rayleigh-Plesset equation
pub struct BubbleDynamicsSolver {
    /// Bubble configurations per grid cell
    configs: HashMap<(usize, usize, usize), RayleighPlesset<f64>>,
    /// Current bubble radii per grid cell
    radii: HashMap<(usize, usize, usize), f64>,
    /// Current bubble wall velocities per grid cell
    velocities: HashMap<(usize, usize, usize), f64>,
}

impl BubbleDynamicsSolver {
    /// Create new bubble dynamics solver
    pub fn new(
        config: &BubbleDynamicsConfig,
        nx: usize,
        ny: usize,
        nz: usize,
        liquid_density: f64,
        liquid_viscosity: f64,
        vapor_pressure: f64,
    ) -> Self {
        let mut configs = HashMap::new();
        let mut radii = HashMap::new();
        let mut velocities = HashMap::new();

        // Initialize bubbles at each grid cell
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let rp = RayleighPlesset {
                        initial_radius: config.initial_radius,
                        liquid_density,
                        liquid_viscosity,
                        surface_tension: config.surface_tension,
                        vapor_pressure,
                        polytropic_index: config.polytropic_exponent,
                    };
                    configs.insert((i, j, k), rp);
                    radii.insert((i, j, k), config.initial_radius);
                    velocities.insert((i, j, k), 0.0);
                }
            }
        }

        Self {
            configs,
            radii,
            velocities,
        }
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
        let key = (i, j, k);
        let config = self
            .configs
            .get(&key)
            .ok_or_else(|| cfd_core::error::Error::Solver("Bubble config not found".to_string()))?;
        let radius = *self.radii.get(&key).unwrap_or(&config.initial_radius);
        let velocity = *self.velocities.get(&key).unwrap_or(&0.0);

        // Calculate acceleration using Rayleigh-Plesset equation
        let acceleration = config.bubble_acceleration(radius, velocity, pressure);

        // Simple semi-implicit Euler integration
        let new_velocity = velocity + acceleration * dt;
        let new_radius = (radius + new_velocity * dt).max(1e-9);

        self.radii.insert(key, new_radius);
        self.velocities.insert(key, new_velocity);

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
        let key = (i, j, k);
        let radius = *self.radii.get(&key).unwrap_or(&1e-6);
        let initial_radius = self.configs.get(&key).map_or(1e-6, |c| c.initial_radius);

        if radius > 0.0 {
            // Impact pressure scales with (R_max / R_collapse)
            liquid_density * sound_speed * sound_speed * (initial_radius / radius)
        } else {
            0.0
        }
    }

    /// Get bubble natural frequency (Hz) for impact frequency estimation
    pub fn get_bubble_frequency(&self, i: usize, j: usize, k: usize, ambient_pressure: f64) -> f64 {
        let key = (i, j, k);
        let config = self.configs.get(&key);
        let radius = self.radii.get(&key);

        if let (Some(config), Some(&r)) = (config, radius) {
            let omega = config.natural_frequency(r, ambient_pressure);
            omega / (2.0 * std::f64::consts::PI)
        } else {
            0.0
        }
    }
}
