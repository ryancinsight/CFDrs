//! Nuclei Transport Model for Cavitation Re-nucleation Cascade
//!
//! # Theorem — Nuclei Advection-Diffusion
//!
//! The transport of uncollapsed cavitation nuclei is governed by the scalar advection-diffusion equation:
//!
//! ```math
//! \frac{\partial n}{\partial t} + \nabla \cdot (n u) = \nabla \cdot (D \nabla n) + S_{gen} - S_{diss}
//! ```
//!
//! where $n$ is the nuclei volume fraction, $u$ is the fluid velocity vector, $D$ is the turbulent diffusion
//! coefficient (or effective kinematic viscosity), $S_{gen}$ is the cavitation generation source term, and
//! $S_{diss} = \frac{n}{\tau_{diss}}$ is the dissolution sink term based on residence time.
//!
//! **Proof sketch**: Conservation of mass for a passively transported dispersed scalar field. The nuclei represent
//! remnants of collapsed macroscopic cavitation bubbles that advect with the mean flow. By the Reynolds transport
//! theorem, the rate of change of nuclei within a material volume equals the net generation minus diffusion out of
//! the parcel. The exponential decay model for dissolution assumes nuclei dissolve linearly proportional to their
//! concentration at a rate scaled by the characteristic relaxation time $\tau_{diss}$.

use nalgebra::RealField;
use serde::{Deserialize, Serialize};

/// Configuration parameters for the Nuclei Transport model
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct NucleiTransportConfig<T: RealField + Copy> {
    /// Characteristic dissolution time [s] ($\tau_{diss}$)
    pub dissolution_time_s: T,
    /// Proportional generation factor from active cavitation events
    pub generation_rate_factor: T,
    /// Diffusion coefficient for the scalar field [m^2/s]
    pub diffusion_coefficient: T,
}

impl<T: RealField + Copy> Default for NucleiTransportConfig<T> {
    fn default() -> Self {
        Self {
            dissolution_time_s: T::from_f64(0.05).unwrap_or_else(T::one), // default 50ms
            generation_rate_factor: T::from_f64(1.0).unwrap_or_else(T::one),
            diffusion_coefficient: T::from_f64(1e-6).unwrap_or_else(T::one),
        }
    }
}

/// Evaluates the nuclei source and sink terms at a point in the flow.
pub struct NucleiTransport<T: RealField + Copy> {
    config: NucleiTransportConfig<T>,
}

impl<T: RealField + Copy> NucleiTransport<T> {
    /// Create a new nuclei transport evaluator
    #[must_use]
    pub const fn new(config: NucleiTransportConfig<T>) -> Self {
        Self { config }
    }

    /// Calculate the dissolution sink term $S_{diss} = n / \tau_{diss}$
    /// Returns the rate of change ($dn/dt$).
    #[must_use]
    pub fn calculate_dissolution_rate(&self, current_nuclei_fraction: T) -> T {
        if self.config.dissolution_time_s > T::zero() {
            -current_nuclei_fraction / self.config.dissolution_time_s
        } else {
            T::zero()
        }
    }

    /// Calculate the generation source term $S_{gen}$
    /// based on the local macroscopic cavitation source rate.
    #[must_use]
    pub fn calculate_generation_rate(&self, macroscopic_cavitation_source: T) -> T {
        if macroscopic_cavitation_source > T::zero() {
            macroscopic_cavitation_source * self.config.generation_rate_factor
        } else {
            T::zero()
        }
    }

    /// Total scalar reaction rate (generation - dissolution).
    /// Does not include advection or diffusion operators.
    #[must_use]
    pub fn calculate_net_reaction_rate(
        &self,
        current_nuclei_fraction: T,
        macroscopic_cavitation_source: T,
    ) -> T {
        self.calculate_generation_rate(macroscopic_cavitation_source)
            + self.calculate_dissolution_rate(current_nuclei_fraction)
    }

    /// For 1D / network models: analytical integration of the dissolution decay
    /// over a given transit time $\Delta t$.
    ///
    /// $n(t + \Delta t) = n(t) \exp(-\Delta t / \tau_{diss})$
    #[must_use]
    pub fn advect_1d_dissolution(&self, upstream_nuclei: T, transit_time: T) -> T {
        if self.config.dissolution_time_s > T::zero() {
            let exponent = -transit_time / self.config.dissolution_time_s;
            upstream_nuclei * exponent.exp()
        } else {
            upstream_nuclei
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_dissolution_rate() {
        let config = NucleiTransportConfig {
            dissolution_time_s: 0.1, // 100ms
            generation_rate_factor: 1.0,
            diffusion_coefficient: 1e-6,
        };
        let transport = NucleiTransport::new(config);

        let rate = transport.calculate_dissolution_rate(0.5);
        // -0.5 / 0.1 = -5.0
        assert_relative_eq!(rate, -5.0);
    }

    #[test]
    fn test_advect_1d_dissolution() {
        let config = NucleiTransportConfig {
            dissolution_time_s: 0.05, // 50ms
            generation_rate_factor: 1.0,
            diffusion_coefficient: 1e-6,
        };
        let transport = NucleiTransport::new(config);

        let dt = 0.05; // 1 time constant
        let n_out = transport.advect_1d_dissolution(1.0, dt);
        assert_relative_eq!(n_out, std::f64::consts::E.recip()); // exp(-1)
    }
}
