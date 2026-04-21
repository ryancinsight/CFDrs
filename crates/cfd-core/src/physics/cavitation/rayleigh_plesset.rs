//! Rayleigh-Plesset bubble dynamics model.
//!
//! # Theorem — Rayleigh-Plesset Equation (Rayleigh 1917, Plesset 1949)
//!
//! A spherical gas/vapour bubble of radius $R(t)$ in an infinite incompressible
//! liquid of density $\rho_L$, viscosity $\mu$, and surface tension $\sigma$
//! satisfies the second-order ODE:
//!
//! ```text
//! ρ_L (R R̈ + 3/2 Ṙ²) = p_B(R) − p_∞(t) − 2σ/R − 4μṘ/R
//! ```
//!
//! where $p_B(R) = p_{v} + p_{g0}(R_0/R)^{3\kappa}$ is the internal bubble
//! pressure (vapour + polytropic gas), $p_\infty(t)$ is the far-field liquid
//! pressure, and $\kappa$ is the polytropic index.
//!
//! **Proof sketch.** Assume spherical symmetry and incompressibility so that
//! the radial velocity is $u_r = \dot{R}R^2/r^2$. Substituting into the radial
//! Navier-Stokes equation and integrating from $r = R$ to $r \to \infty$ gives
//! the inertial terms on the left. The normal-stress balance at the bubble
//! surface provides the pressure and surface-tension terms. Viscous dissipation
//! contributes the $4\mu\dot{R}/R$ term (Plesset & Prosperetti 1977).
//!
//! A collapsed bubble is treated as an absorbing state: when `R = 0`, the
//! canonical integrator returns `(R, \dot{R}) = (0, 0)` and the bubble remains
//! collapsed until a separate nucleation process re-seeds it.
//!
//! ## References
//!
//! - Rayleigh, Lord (1917). "On the pressure developed in a liquid during the
//!   collapse of a spherical cavity." *Phil. Mag.* 34:94–98.
//! - Plesset, M. S. (1949). "The dynamics of cavitation bubbles." *J. Appl. Mech.*
//!   16:277–282.
//! - Plesset, M. S. & Prosperetti, A. (1977). "Bubble dynamics and cavitation."
//!   *Ann. Rev. Fluid Mech.* 9:145–185.

use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

use crate::error::{Error, NumericalErrorKind, Result};

/// Rayleigh-Plesset bubble dynamics model
///
/// Describes the growth and collapse of a spherical bubble in an infinite liquid
/// based on Rayleigh (1917) and Plesset (1949).
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct RayleighPlesset<T: RealField + Copy> {
    /// Initial bubble radius (m)
    pub initial_radius: T,
    /// Liquid density (kg/m³)
    pub liquid_density: T,
    /// Liquid viscosity (Pa·s)
    pub liquid_viscosity: T,
    /// Surface tension (N/m)
    pub surface_tension: T,
    /// Vapor pressure (Pa)
    pub vapor_pressure: T,
    /// Polytropic index for gas content
    pub polytropic_index: T,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[allow(missing_docs)]
pub struct SonoluminescenceEstimate<T: RealField + Copy> {
    #[allow(missing_docs)]
    pub peak_temperature: T,
    #[allow(missing_docs)]
    pub peak_pressure: T,
    #[allow(missing_docs)]
    pub radiated_energy: T,
}

impl<T: RealField + Copy + FromPrimitive> RayleighPlesset<T> {
    /// Calculate bubble wall acceleration (d²R/dt²)
    ///
    /// Implements the full Rayleigh-Plesset equation:
    ///
    /// ```text
    /// R̈ = [p_B(R) − p_∞ − 2σ/R − 4μṘ/R] / (ρ_L R) − (3/2) Ṙ²/R
    /// ```
    ///
    /// where the internal bubble pressure is:
    ///
    /// ```text
    /// p_B(R) = p_v + p_g0 · (R_0/R)^{3κ}
    /// ```
    ///
    /// and the equilibrium gas pressure `p_g0` is determined from the
    /// initial condition: at `R = R_0` the bubble is stationary (Ṙ = R̈ = 0),
    /// and the normal-stress balance gives `p_g0 = p_∞ − p_v + 2σ/R_0`.
    ///
    /// **References**: Plesset & Prosperetti (1977), Brennen (1995) §2.2.
    pub fn bubble_acceleration(&self, radius: T, velocity: T, ambient_pressure: T) -> T {
        if radius <= T::zero() {
            return T::zero();
        }

        let two = T::from_f64(2.0).unwrap_or_else(|| T::one() + T::one());
        let three = T::from_f64(3.0).unwrap_or_else(|| T::one() + T::one() + T::one());
        let four = two * two;
        let three_halves = T::from_f64(1.5).unwrap_or_else(|| T::one());

        // Equilibrium gas pressure from normal-stress balance at R = R_0:
        //   p_v + p_g0 = p_∞ + 2σ/R_0  →  p_g0 = p_∞ − p_v + 2σ/R_0
        // Using ambient_pressure as the reference far-field pressure p_∞.
        let p_g0 = (ambient_pressure - self.vapor_pressure
            + two * self.surface_tension / self.initial_radius)
            .max(T::zero());

        // Polytropic gas contribution: p_g = p_g0 · (R_0/R)^{3κ}
        let ratio = self.initial_radius / radius;
        let exponent = three * self.polytropic_index;
        let p_gas = p_g0 * ratio.powf(exponent);

        // Internal bubble pressure: p_B(R) = p_v + p_gas
        let p_bubble = self.vapor_pressure + p_gas;

        // Pressure difference: p_B − p_∞
        let pressure_diff = p_bubble - ambient_pressure;

        // Rayleigh-Plesset equation terms:
        //   R̈ = (p_B − p_∞)/(ρR) − (3/2)Ṙ²/R − 4μṘ/(ρR²) − 2σ/(ρR²)
        let term1 = pressure_diff / (self.liquid_density * radius);
        let term2 = -three_halves * velocity * velocity / radius;
        let term3 =
            -four * self.liquid_viscosity * velocity / (self.liquid_density * radius * radius);
        let term4 = -two * self.surface_tension / (self.liquid_density * radius * radius);

        term1 + term2 + term3 + term4
    }

    #[allow(missing_docs)]
    pub fn step_semi_implicit(
        &self,
        radius: T,
        velocity: T,
        ambient_pressure: T,
        dt: T,
    ) -> Result<(T, T)> {
        if radius < T::zero() {
            return Err(Error::InvalidConfiguration(
                "Bubble radius must be non-negative".to_string(),
            ));
        }
        if dt <= T::zero() {
            return Err(Error::InvalidConfiguration(
                "Time step must be positive".to_string(),
            ));
        }
        if !ambient_pressure.is_finite() {
            return Err(Error::InvalidConfiguration(
                "Ambient pressure must be finite".to_string(),
            ));
        }

        if radius == T::zero() {
            return Ok((T::zero(), T::zero()));
        }

        let accel = self.bubble_acceleration(radius, velocity, ambient_pressure);
        if !accel.is_finite() {
            return Err(Error::Numerical(NumericalErrorKind::InvalidValue {
                value: "Rayleigh-Plesset acceleration is non-finite".to_string(),
            }));
        }

        let new_velocity = velocity + accel * dt;
        let new_radius = radius + new_velocity * dt;
        if new_radius <= T::zero() {
            return Ok((T::zero(), T::zero()));
        }

        Ok((new_radius, new_velocity))
    }

    #[allow(missing_docs)]
    pub fn estimate_sonoluminescence(
        &self,
        ambient_pressure: T,
        ambient_temperature: T,
        collapse_radius: T,
        emissivity: T,
        flash_duration: T,
    ) -> Result<SonoluminescenceEstimate<T>> {
        if self.initial_radius <= T::zero() {
            return Err(Error::InvalidConfiguration(
                "Initial bubble radius must be positive".to_string(),
            ));
        }
        if collapse_radius <= T::zero() {
            return Err(Error::InvalidConfiguration(
                "Collapse radius must be positive".to_string(),
            ));
        }
        if collapse_radius > self.initial_radius {
            return Err(Error::InvalidConfiguration(
                "Collapse radius must not exceed initial radius".to_string(),
            ));
        }
        if ambient_temperature <= T::zero() || !ambient_temperature.is_finite() {
            return Err(Error::InvalidConfiguration(
                "Ambient temperature must be finite and positive".to_string(),
            ));
        }
        if ambient_pressure <= T::zero() || !ambient_pressure.is_finite() {
            return Err(Error::InvalidConfiguration(
                "Ambient pressure must be finite and positive".to_string(),
            ));
        }
        if emissivity < T::zero() || emissivity > T::one() || !emissivity.is_finite() {
            return Err(Error::InvalidConfiguration(
                "Emissivity must be finite and within [0, 1]".to_string(),
            ));
        }
        if flash_duration <= T::zero() || !flash_duration.is_finite() {
            return Err(Error::InvalidConfiguration(
                "Flash duration must be finite and positive".to_string(),
            ));
        }
        if self.polytropic_index <= T::one() || !self.polytropic_index.is_finite() {
            return Err(Error::InvalidConfiguration(
                "Polytropic index must be finite and > 1".to_string(),
            ));
        }

        let three = T::from_f64(3.0).unwrap_or_else(|| T::one() + T::one() + T::one());
        let four = T::from_f64(4.0).unwrap_or_else(|| T::one() + T::one() + T::one() + T::one());
        let pi = T::from_f64(std::f64::consts::PI).unwrap_or_else(T::zero);

        let ratio = self.initial_radius / collapse_radius;
        if !ratio.is_finite() || ratio < T::one() {
            return Err(Error::Numerical(NumericalErrorKind::InvalidValue {
                value: "Invalid radius ratio for sonoluminescence estimate".to_string(),
            }));
        }

        let gamma = self.polytropic_index;
        let exponent_t = three * (gamma - T::one());
        let exponent_p = three * gamma;

        let peak_temperature = ambient_temperature * ratio.powf(exponent_t);
        let peak_pressure = ambient_pressure * ratio.powf(exponent_p);

        if !peak_temperature.is_finite() || !peak_pressure.is_finite() {
            return Err(Error::Numerical(NumericalErrorKind::InvalidValue {
                value: "Non-finite peak temperature/pressure in sonoluminescence estimate"
                    .to_string(),
            }));
        }

        let sigma_sb = T::from_f64(5.670_374_419e-8).unwrap_or_else(T::zero);
        if sigma_sb <= T::zero() {
            return Err(Error::InvalidConfiguration(
                "Failed to construct Stefan-Boltzmann constant".to_string(),
            ));
        }

        let area = four * pi * collapse_radius * collapse_radius;
        let radiated_energy =
            emissivity * sigma_sb * peak_temperature.powi(4) * area * flash_duration;

        if !radiated_energy.is_finite() || radiated_energy < T::zero() {
            return Err(Error::Numerical(NumericalErrorKind::InvalidValue {
                value: "Non-finite radiated energy in sonoluminescence estimate".to_string(),
            }));
        }

        Ok(SonoluminescenceEstimate {
            peak_temperature,
            peak_pressure,
            radiated_energy,
        })
    }

    /// Calculate bubble growth rate for inviscid case
    /// Based on Rayleigh (1917) equation for bubble growth without viscosity
    pub fn growth_rate_inviscid(&self, _radius: T, ambient_pressure: T) -> T {
        let pressure_difference = self.vapor_pressure - ambient_pressure;

        if pressure_difference <= T::zero() {
            return T::zero(); // No growth if ambient pressure exceeds vapor pressure
        }

        // Rayleigh equation for inviscid growth: dR/dt = sqrt(2/3 * Δp/ρ)
        // This is the exact solution for constant pressure difference
        let two_thirds = T::from_f64(2.0 / 3.0).unwrap_or_else(|| T::one());
        (two_thirds * pressure_difference / self.liquid_density).sqrt()
    }

    /// Calculate collapse time from Rayleigh collapse formula
    pub fn collapse_time(&self, initial_radius: T, pressure_difference: T) -> T {
        // Rayleigh collapse time: t_c = 0.915 * R_0 * sqrt(ρ/Δp)
        let coefficient = T::from_f64(0.915).unwrap_or_else(|| T::one());

        if pressure_difference > T::zero() {
            coefficient * initial_radius * (self.liquid_density / pressure_difference).sqrt()
        } else {
            T::from_f64(1e10).unwrap_or_else(|| T::one()) // Very large time for no pressure difference
        }
    }

    /// Calculate maximum bubble radius during growth (Rayleigh-Plesset)
    pub fn maximum_radius(&self, pressure_ratio: T) -> T {
        // R_max/R_0 = (p_∞/p_v)^(1/3γ) for isothermal growth
        let exponent =
            T::one() / (T::from_f64(3.0).unwrap_or_else(|| T::one()) * self.polytropic_index);
        self.initial_radius * pressure_ratio.powf(exponent)
    }

    /// Calculate bubble natural frequency
    pub fn natural_frequency(&self, radius: T, ambient_pressure: T) -> T {
        if radius <= T::zero() {
            return T::zero();
        }

        // ω_0 = 1/R * sqrt(3γ(p_0 - p_v)/ρ)
        let three = T::from_f64(3.0).unwrap_or_else(|| T::one() + T::one() + T::one());
        let pressure_term =
            three * self.polytropic_index * (ambient_pressure - self.vapor_pressure);

        if pressure_term > T::zero() {
            (pressure_term / self.liquid_density).sqrt() / radius
        } else {
            T::zero()
        }
    }

    /// Calculate critical Blake radius for unstable growth
    pub fn blake_critical_radius(&self, ambient_pressure: T) -> T {
        // R_c = 0.85 * (2σ/(p_∞ - p_v))
        let coefficient =
            T::from_f64(super::constants::BLAKE_CRITICAL_COEFFICIENT).unwrap_or_else(|| T::one());
        let two = T::from_f64(2.0).unwrap_or_else(|| T::one() + T::one());

        let pressure_diff = ambient_pressure - self.vapor_pressure;
        if pressure_diff > T::zero() {
            coefficient * two * self.surface_tension / pressure_diff
        } else {
            T::from_f64(1e10).unwrap_or_else(|| T::one()) // Very large radius
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn step_semi_implicit_preserves_positive_radius() {
        let rp = RayleighPlesset::<f64> {
            initial_radius: 1e-6,
            liquid_density: 998.0,
            liquid_viscosity: 1.002e-3,
            surface_tension: 0.0728,
            vapor_pressure: 2339.0,
            polytropic_index: 1.4,
        };

        let (r1, v1) = rp.step_semi_implicit(1e-6, 0.0, 1.0e5, 1e-7).unwrap();
        assert!(r1 > 0.0);
        assert!(v1.is_finite());
    }

    #[test]
    fn step_semi_implicit_treats_collapsed_bubble_as_absorbing_state() {
        let rp = RayleighPlesset::<f64> {
            initial_radius: 1e-6,
            liquid_density: 998.0,
            liquid_viscosity: 1.002e-3,
            surface_tension: 0.0728,
            vapor_pressure: 2339.0,
            polytropic_index: 1.4,
        };

        let (r1, v1) = rp.step_semi_implicit(0.0, 0.0, 1.0e5, 1e-7).unwrap();
        assert_eq!(r1, 0.0);
        assert_eq!(v1, 0.0);
    }

    #[test]
    fn sonoluminescence_energy_increases_with_stronger_collapse() {
        let rp = RayleighPlesset::<f64> {
            initial_radius: 10e-6,
            liquid_density: 998.0,
            liquid_viscosity: 1.002e-3,
            surface_tension: 0.0728,
            vapor_pressure: 2339.0,
            polytropic_index: 1.4,
        };

        let ambient_pressure = 1.0e5;
        let ambient_temperature = 293.15;
        let flash_duration = 50e-12;
        let emissivity = 1.0;

        let e_weak = rp
            .estimate_sonoluminescence(
                ambient_pressure,
                ambient_temperature,
                5e-6,
                emissivity,
                flash_duration,
            )
            .unwrap()
            .radiated_energy;

        let e_strong = rp
            .estimate_sonoluminescence(
                ambient_pressure,
                ambient_temperature,
                1e-6,
                emissivity,
                flash_duration,
            )
            .unwrap()
            .radiated_energy;

        assert!(e_strong > e_weak);
    }
}
