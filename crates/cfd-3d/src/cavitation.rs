//! Public cavitation closures for three-dimensional flow.
//!
//! The VOF solver owns dense fields and the full phase-transfer coupling. This
//! module owns the small closure seam used by callers that need an inception
//! threshold and a single-step collapse exposure without constructing a VOF
//! mesh. The dimensional contracts remain Aequitas quantities; the analytical
//! kernels in `cfd-core` extract base values at their formula boundary.

use aequitas::systems::si::quantities::{
    DynamicViscosity, Length, MassDensity, MassDensityRate, NumberDensity, Pressure,
    SurfaceTension, Time,
};
use cfd_core::error::{Error, Result};
use cfd_core::physics::cavitation::{
    models::CavitationModel, rayleigh_plesset::RayleighPlesset as CoreRayleighPlesset,
};

/// Cavitation-model trait.
///
/// The collapse rate is the non-negative mean liquid-mass-density transfer
/// rate caused by condensation at a pressure above the inception threshold.
/// Vaporization is represented by zero because this seam reports collapse
/// exposure only; the VOF phase-transfer model reports signed transfer.
pub trait Cavitation {
    /// Local cavity inception pressure.
    fn inception_pressure(&self) -> Pressure<f64>;

    /// Calculate the collapse exposure at the local pressure.
    fn collapse_rate(&self, pressure: Pressure<f64>) -> Result<MassDensityRate<f64>>;
}

/// Rayleigh-Plesset closure for a single bubble.
#[derive(Clone, Copy, Debug)]
pub struct RayleighPlesset {
    /// Initial bubble radius.
    pub initial_radius: Length<f64>,
    /// Liquid density.
    pub liquid_density: MassDensity<f64>,
    /// Liquid dynamic viscosity.
    pub liquid_viscosity: DynamicViscosity<f64>,
    /// Surface tension.
    pub surface_tension: SurfaceTension<f64>,
    /// Vapor pressure of the working fluid.
    pub vapor_pressure: Pressure<f64>,
    /// Polytropic index for the bubble gas.
    pub polytropic_index: f64,
}

impl Default for RayleighPlesset {
    fn default() -> Self {
        Self {
            initial_radius: Length::from_base(1.0e-6),
            liquid_density: MassDensity::from_base(cfd_core::physics::constants::physics::fluid::WATER_DENSITY),
            liquid_viscosity: DynamicViscosity::from_base(cfd_core::physics::constants::physics::fluid::WATER_VISCOSITY),
            surface_tension: SurfaceTension::from_base(cfd_core::physics::cavitation::constants::SURFACE_TENSION_WATER_VALUE),
            vapor_pressure: Pressure::from_base(cfd_core::physics::cavitation::constants::VAPOR_PRESSURE_WATER_20C_VALUE),
            polytropic_index: 1.4,
        }
    }
}

impl RayleighPlesset {
    fn core_model(&self) -> CoreRayleighPlesset<f64> {
        CoreRayleighPlesset {
            initial_radius: self.initial_radius,
            liquid_density: self.liquid_density,
            liquid_viscosity: self.liquid_viscosity,
            surface_tension: self.surface_tension,
            vapor_pressure: self.vapor_pressure,
            polytropic_index: self.polytropic_index,
        }
    }

    fn validate(&self) -> Result<()> {
        let fields = [
            self.initial_radius.into_base(),
            self.liquid_density.into_base(),
            self.liquid_viscosity.into_base(),
            self.surface_tension.into_base(),
            self.vapor_pressure.into_base(),
            self.polytropic_index,
        ];
        if !fields.iter().all(|value| value.is_finite()) {
            return Err(Error::InvalidConfiguration(
                "Rayleigh-Plesset closure values must be finite".to_string(),
            ));
        }
        if self.initial_radius.into_base() <= 0.0
            || self.liquid_density.into_base() <= 0.0
            || self.liquid_viscosity.into_base() < 0.0
            || self.surface_tension.into_base() < 0.0
            || self.vapor_pressure.into_base() < 0.0
            || self.polytropic_index <= 1.0
        {
            return Err(Error::InvalidConfiguration(
                "Rayleigh-Plesset closure values violate positivity constraints".to_string(),
            ));
        }
        Ok(())
    }
}

impl Cavitation for RayleighPlesset {
    fn inception_pressure(&self) -> Pressure<f64> {
        self.vapor_pressure
    }

    fn collapse_rate(&self, pressure: Pressure<f64>) -> Result<MassDensityRate<f64>> {
        self.validate()?;
        let pressure = pressure.into_base();
        if !pressure.is_finite() || pressure <= self.vapor_pressure.into_base() {
            return if pressure.is_finite() {
                Ok(MassDensityRate::from_base(0.0))
            } else {
                Err(Error::InvalidConfiguration(
                    "Collapse pressure must be finite".to_string(),
                ))
            };
        }

        let collapse_time = self.core_model().collapse_time(
            self.initial_radius.into_base(),
            pressure - self.vapor_pressure.into_base(),
        );
        let rate = self.liquid_density.into_base() / collapse_time;
        if !rate.is_finite() || rate < 0.0 {
            return Err(Error::InvalidConfiguration(
                "Rayleigh-Plesset collapse rate is non-finite".to_string(),
            ));
        }
        Ok(MassDensityRate::from_base(rate))
    }
}

/// Eulerian-Eulerian closure for cloud cavitation.
#[derive(Clone, Copy, Debug)]
pub struct EulerianEulerian {
    /// Threshold bubble number density.
    pub nucleation_density: NumberDensity<f64>,
    /// Vapor pressure of the working fluid.
    pub vapor_pressure: Pressure<f64>,
    /// Liquid density.
    pub liquid_density: MassDensity<f64>,
    /// Vapor density.
    pub vapor_density: MassDensity<f64>,
    /// Current vapor volume fraction.
    pub volume_fraction: f64,
    /// Representative bubble radius.
    pub bubble_radius: Length<f64>,
}

impl Default for EulerianEulerian {
    fn default() -> Self {
        Self {
            nucleation_density: NumberDensity::from_base(1.0e13),
            vapor_pressure: Pressure::from_base(cfd_core::physics::cavitation::constants::VAPOR_PRESSURE_WATER_20C_VALUE),
            liquid_density: MassDensity::from_base(cfd_core::physics::constants::physics::fluid::WATER_DENSITY),
            vapor_density: MassDensity::from_base(0.0173),
            volume_fraction: 0.5,
            bubble_radius: Length::from_base(1.0e-6),
        }
    }
}

impl EulerianEulerian {
    fn validate(&self) -> Result<()> {
        let fields = [
            self.nucleation_density.into_base(),
            self.vapor_pressure.into_base(),
            self.liquid_density.into_base(),
            self.vapor_density.into_base(),
            self.volume_fraction,
            self.bubble_radius.into_base(),
        ];
        if !fields.iter().all(|value| value.is_finite()) {
            return Err(Error::InvalidConfiguration(
                "Eulerian-Eulerian closure values must be finite".to_string(),
            ));
        }
        if self.nucleation_density.into_base() <= 0.0
            || self.vapor_pressure.into_base() < 0.0
            || self.liquid_density.into_base() <= 0.0
            || self.vapor_density.into_base() <= 0.0
            || !(0.0..=1.0).contains(&self.volume_fraction)
            || self.bubble_radius.into_base() <= 0.0
        {
            return Err(Error::InvalidConfiguration(
                "Eulerian-Eulerian closure values violate physical bounds".to_string(),
            ));
        }
        Ok(())
    }
}

impl Cavitation for EulerianEulerian {
    fn inception_pressure(&self) -> Pressure<f64> {
        self.vapor_pressure
    }

    fn collapse_rate(&self, pressure: Pressure<f64>) -> Result<MassDensityRate<f64>> {
        self.validate()?;
        let pressure_base = pressure.into_base();
        if !pressure_base.is_finite() {
            return Err(Error::InvalidConfiguration(
                "Collapse pressure must be finite".to_string(),
            ));
        }
        let model = CavitationModel::SchnerrSauer {
            bubble_density: self.nucleation_density,
            initial_radius: self.bubble_radius,
        };
        let signed_rate = model.mass_transfer_rate(
            pressure,
            self.vapor_pressure,
            self.volume_fraction,
            self.liquid_density,
            self.vapor_density,
        );
        let collapse_rate = (-signed_rate.into_base()).max(0.0);
        Ok(MassDensityRate::from_base(collapse_rate))
    }
}

/// Accumulate one timestep of collapse exposure.
pub fn damage_step(rate: MassDensityRate<f64>, dt: Time<f64>) -> Result<MassDensity<f64>> {
    let rate = rate.into_base();
    let dt = dt.into_base();
    if !rate.is_finite() || !dt.is_finite() || dt <= 0.0 {
        return Err(Error::InvalidConfiguration(
            "Collapse exposure rate must be finite and time step must be positive".to_string(),
        ));
    }
    let increment = rate * dt;
    if !increment.is_finite() {
        return Err(Error::InvalidConfiguration(
            "Collapse exposure increment is non-finite".to_string(),
        ));
    }
    Ok(MassDensity::from_base(increment))
}

#[cfg(test)]
mod tests {
    use super::{damage_step, Cavitation, EulerianEulerian, RayleighPlesset};
    use aequitas::systems::si::quantities::{MassDensityRate, Pressure, Time};

    #[test]
    fn rayleigh_plesset_collapse_is_pressure_sensitive() {
        let model = RayleighPlesset::default();
        let no_collapse = model
            .collapse_rate(Pressure::from_base(2_339.0))
            .expect("valid threshold pressure");
        let collapse = model
            .collapse_rate(Pressure::from_base(100_000.0))
            .expect("valid collapse pressure");

        assert_eq!(no_collapse.into_base(), 0.0);
        assert!(collapse.into_base() > 0.0);
    }

    #[test]
    fn eulerian_collapse_uses_signed_phase_transfer_model() {
        let model = EulerianEulerian::default();
        let collapse = model
            .collapse_rate(Pressure::from_base(100_000.0))
            .expect("valid collapse pressure");
        let vaporization = model
            .collapse_rate(Pressure::from_base(500.0))
            .expect("valid vaporization pressure");

        assert!(collapse.into_base() > 0.0);
        assert_eq!(vaporization.into_base(), 0.0);
    }

    #[test]
    fn damage_step_preserves_mass_density_rate_and_time_values() {
        let increment = damage_step(MassDensityRate::from_base(4.0), Time::from_base(0.25))
            .expect("positive timestep");

        assert_eq!(increment.into_base(), 1.0);
    }
}
