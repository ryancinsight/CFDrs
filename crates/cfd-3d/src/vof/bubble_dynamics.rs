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

use aequitas::systems::si::quantities::{
    DynamicViscosity, Frequency, Length, MassDensity, NumberDensity, Pressure, SurfaceTension,
    Time, Velocity,
};
use cfd_core::error::{Error, Result};
use cfd_core::physics::cavitation::rayleigh_plesset::RayleighPlesset;
use cfd_core::physics::fluid::BloodModel;
use leto::geometry::Vector3;

/// Bubble dynamics configuration
#[derive(Debug, Clone)]
pub struct BubbleDynamicsConfig {
    /// Initial bubble radius (m)
    pub initial_radius: Length<f64>,
    /// Bubble number density (1/m³)
    ///
    /// This is converted to an expected bubble population per cell using the
    /// control-volume size so the per-cell damage and sonoluminescence outputs
    /// remain consistent with the modeled nuclei density.
    pub number_density: NumberDensity<f64>,
    /// Polytropic exponent for gas compression
    pub polytropic_exponent: f64,
    /// Surface tension coefficient (N/m)
    pub surface_tension: SurfaceTension<f64>,
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
        dx: Length<f64>,
        dy: Length<f64>,
        dz: Length<f64>,
        liquid_density: MassDensity<f64>,
        blood_model: BloodModel<f64>,
        vapor_pressure: Pressure<f64>,
    ) -> Self {
        let initial_viscosity = blood_model.viscosity(0.0);
        let rp = RayleighPlesset {
            initial_radius: config.initial_radius,
            liquid_density,
            liquid_viscosity: DynamicViscosity::from_base(initial_viscosity),
            surface_tension: config.surface_tension,
            vapor_pressure,
            polytropic_index: config.polytropic_exponent,
        };
        let len = nx * ny * nz;
        let bubble_population_weight =
            config.number_density.into_base() * dx.into_base() * dy.into_base() * dz.into_base();

        Self {
            nx,
            ny,
            configs: vec![rp; len],
            radii: vec![config.initial_radius.into_base(); len],
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
        pressure: Pressure<f64>,
        _velocity: Vector3<f64>,
        density: MassDensity<f64>,
        dt: Time<f64>,
    ) -> Result<Length<f64>> {
        let density = density.into_base();
        if !density.is_finite() || density <= 0.0 {
            return Err(Error::InvalidConfiguration(
                "local liquid density must be finite and positive".to_string(),
            ));
        }

        let idx = self.index(i, j, k);
        let base_viscosity = self.blood_model.viscosity(0.0);
        let config = &mut self.configs[idx];
        let radius = self.radii[idx];
        let velocity = self.velocities[idx];
        config.liquid_density = MassDensity::from_base(density);

        if radius <= 0.0 {
            config.liquid_viscosity = DynamicViscosity::from_base(base_viscosity);
            self.radii[idx] = 0.0;
            self.velocities[idx] = 0.0;
            return Ok(Length::from_base(0.0));
        }

        // Calculate apparent viscosity at the bubble wall.
        // theorem: The shear rate at a spherical bubble wall expanding/collapsing radially
        // is given by γ̇_wall = √12 |Ṙ / R|.
        let shear_rate = 12.0_f64.sqrt() * (velocity.abs() / radius);
        config.liquid_viscosity =
            DynamicViscosity::from_base(self.blood_model.viscosity(shear_rate));

        let (new_radius, new_velocity) =
            config.step_semi_implicit(radius, velocity, pressure.into_base(), dt.into_base())?;

        self.radii[idx] = new_radius;
        self.velocities[idx] = new_velocity;

        Ok(Length::from_base(new_radius))
    }

    /// Get collapse pressure for damage calculation
    pub fn collapse_pressure(
        &self,
        i: usize,
        j: usize,
        k: usize,
        liquid_density: MassDensity<f64>,
        sound_speed: Velocity<f64>,
    ) -> Pressure<f64> {
        let idx = self.index(i, j, k);
        let config = &self.configs[idx];
        let radius = self.radii[idx];
        let initial_radius = config.initial_radius.into_base();

        if radius > 0.0 {
            // Impact pressure scales with (R_max / R_collapse)
            Pressure::from_base(
                liquid_density.into_base()
                    * sound_speed.into_base()
                    * sound_speed.into_base()
                    * (initial_radius / radius),
            )
        } else {
            Pressure::from_base(0.0)
        }
    }

    /// Get bubble natural frequency (Hz) for impact frequency estimation
    pub fn get_bubble_frequency(
        &self,
        i: usize,
        j: usize,
        k: usize,
        ambient_pressure: Pressure<f64>,
    ) -> Frequency<f64> {
        let idx = self.index(i, j, k);
        let config = &self.configs[idx];
        let radius = self.radii[idx];

        let omega = config.natural_frequency(radius, ambient_pressure.into_base());
        Frequency::from_base(omega / (2.0 * std::f64::consts::PI))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use eunomia::assert_relative_eq;

    #[test]
    fn collapse_pressure_and_frequency_use_initialized_bubble_state() {
        let config = BubbleDynamicsConfig {
            initial_radius: Length::from_base(2.0e-6),
            number_density: NumberDensity::from_base(1.0e12),
            polytropic_exponent: 1.4,
            surface_tension: SurfaceTension::from_base(0.072),
        };
        let solver = BubbleDynamicsSolver::new(
            &config,
            1,
            1,
            1,
            Length::from_base(1.0),
            Length::from_base(1.0),
            Length::from_base(1.0),
            MassDensity::from_base(1000.0),
            BloodModel::Newtonian(1.0e-3),
            Pressure::from_base(2300.0),
        );

        let collapse_pressure = solver.collapse_pressure(
            0,
            0,
            0,
            MassDensity::from_base(1000.0),
            Velocity::from_base(1500.0),
        );
        let frequency = solver.get_bubble_frequency(0, 0, 0, Pressure::from_base(2500.0));

        assert!(collapse_pressure.into_base() > 0.0);
        assert!(frequency.into_base() > 0.0);
    }

    #[test]
    fn update_bubble_uses_local_liquid_density_and_rejects_nonphysical_inputs() {
        let config = BubbleDynamicsConfig {
            initial_radius: Length::from_base(2.0e-6),
            number_density: NumberDensity::from_base(1.0e12),
            polytropic_exponent: 1.4,
            surface_tension: SurfaceTension::from_base(0.072),
        };

        let make_solver = || {
            BubbleDynamicsSolver::new(
                &config,
                1,
                1,
                1,
                Length::from_base(1.0),
                Length::from_base(1.0),
                Length::from_base(1.0),
                MassDensity::from_base(1000.0),
                BloodModel::Newtonian(1.0e-3),
                Pressure::from_base(2300.0),
            )
        };

        let mut low_density = make_solver();
        let mut high_density = make_solver();

        low_density.radii[0] = 4.0e-6;
        high_density.radii[0] = 4.0e-6;

        let low = low_density
            .update_bubble(
                0,
                0,
                0,
                Pressure::from_base(1.0e5),
                Vector3::zeros(),
                MassDensity::from_base(800.0),
                Time::from_base(1.0e-7),
            )
            .expect("low-density update");
        let high = high_density
            .update_bubble(
                0,
                0,
                0,
                Pressure::from_base(1.0e5),
                Vector3::zeros(),
                MassDensity::from_base(1600.0),
                Time::from_base(1.0e-7),
            )
            .expect("high-density update");

        assert!(
            low < high,
            "lower liquid density must yield a stronger collapse response"
        );
        assert_relative_eq!(low_density.configs[0].liquid_density.into_base(), 800.0);
        assert_relative_eq!(high_density.configs[0].liquid_density.into_base(), 1600.0);

        let mut invalid_density = make_solver();
        let err = invalid_density
            .update_bubble(
                0,
                0,
                0,
                Pressure::from_base(1.0e5),
                Vector3::zeros(),
                MassDensity::from_base(0.0),
                Time::from_base(1.0e-7),
            )
            .unwrap_err();
        match err {
            Error::InvalidConfiguration(message) => {
                assert!(message.contains("density"));
            }
            other => panic!("expected invalid-configuration error, got {other:?}"),
        }
    }

    #[test]
    fn population_weight_reflects_cell_volume() {
        let config = BubbleDynamicsConfig {
            initial_radius: Length::from_base(2.0e-6),
            number_density: NumberDensity::from_base(2.5e12),
            polytropic_exponent: 1.4,
            surface_tension: SurfaceTension::from_base(0.072),
        };
        let solver = BubbleDynamicsSolver::new(
            &config,
            2,
            2,
            2,
            Length::from_base(0.01),
            Length::from_base(0.02),
            Length::from_base(0.03),
            MassDensity::from_base(1000.0),
            BloodModel::Newtonian(1.0e-3),
            Pressure::from_base(2300.0),
        );

        let expected = config.number_density.into_base() * 0.01 * 0.02 * 0.03;
        assert!((solver.population_weight() - expected).abs() < 1e-12 * expected.max(1.0));
    }

    #[test]
    fn collapsed_bubble_state_remains_absorbing_during_update() {
        let config = BubbleDynamicsConfig {
            initial_radius: Length::from_base(2.0e-6),
            number_density: NumberDensity::from_base(1.0e12),
            polytropic_exponent: 1.4,
            surface_tension: SurfaceTension::from_base(0.072),
        };
        let mut solver = BubbleDynamicsSolver::new(
            &config,
            1,
            1,
            1,
            Length::from_base(1.0),
            Length::from_base(1.0),
            Length::from_base(1.0),
            MassDensity::from_base(1000.0),
            BloodModel::Newtonian(1.0e-3),
            Pressure::from_base(2300.0),
        );

        solver.radii[0] = 0.0;
        solver.velocities[0] = 0.0;

        let radius = solver
            .update_bubble(
                0,
                0,
                0,
                Pressure::from_base(1.0e5),
                Vector3::zeros(),
                MassDensity::from_base(1000.0),
                Time::from_base(1.0e-7),
            )
            .unwrap();

        assert_eq!(radius.into_base(), 0.0);
        assert_eq!(solver.radii[0], 0.0);
        assert_eq!(solver.velocities[0], 0.0);
    }
}
