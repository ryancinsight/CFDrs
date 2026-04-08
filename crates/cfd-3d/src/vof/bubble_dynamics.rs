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
        liquid_density: f64,
        blood_model: BloodModel<f64>,
        vapor_pressure: f64,
    ) -> Self {
        let mut configs = HashMap::new();
        let mut radii = HashMap::new();
        let mut velocities = HashMap::new();

        // Initialize bubbles at each grid cell
        let initial_viscosity = blood_model.viscosity(0.0);
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let rp = RayleighPlesset {
                        initial_radius: config.initial_radius,
                        liquid_density,
                        liquid_viscosity: initial_viscosity,
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
            blood_model,
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
            .get_mut(&key)
            .ok_or_else(|| cfd_core::error::Error::Solver("Bubble config not found".to_string()))?;
        let radius = *self
            .radii
            .get(&key)
            .expect("bubble radius is initialized for every grid cell");
        let velocity = *self
            .velocities
            .get(&key)
            .expect("bubble velocity is initialized for every grid cell");

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
        let config = self
            .configs
            .get(&key)
            .expect("bubble config is initialized for every grid cell");
        let radius = *self
            .radii
            .get(&key)
            .expect("bubble radius is initialized for every grid cell");
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
        let key = (i, j, k);
        let config = self
            .configs
            .get(&key)
            .expect("bubble config is initialized for every grid cell");
        let radius = *self
            .radii
            .get(&key)
            .expect("bubble radius is initialized for every grid cell");

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
            1000.0,
            BloodModel::Newtonian(1.0e-3),
            2300.0,
        );

        let collapse_pressure = solver.collapse_pressure(0, 0, 0, 1000.0, 1500.0);
        let frequency = solver.get_bubble_frequency(0, 0, 0, 2500.0);

        assert!(collapse_pressure > 0.0);
        assert!(frequency > 0.0);
    }
}
