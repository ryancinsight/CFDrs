//! Cavitation models for mass transfer rate calculations.

use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Parameters for ZGB cavitation model calculation
#[derive(Debug, Clone, Copy)]
pub struct ZgbParams<T: RealField + Copy> {
    /// Pressure field value
    pub pressure: T,
    /// Vapor pressure
    pub vapor_pressure: T,
    /// Void fraction
    pub void_fraction: T,
    /// Liquid density
    pub density_liquid: T,
    /// Vapor density
    pub density_vapor: T,
    /// Nucleation site volume fraction
    pub nucleation_fraction: T,
    /// Bubble radius
    pub bubble_radius: T,
    /// Vaporization coefficient
    pub f_vap: T,
    /// Condensation coefficient
    pub f_cond: T,
}

/// Cavitation model types
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum CavitationModel<T: RealField + Copy> {
    /// Kunz model (2000)
    Kunz {
        /// Vaporization coefficient
        vaporization_coeff: T,
        /// Condensation coefficient
        condensation_coeff: T,
    },
    /// Schnerr-Sauer model (2001)
    SchnerrSauer {
        /// Bubble number density (#/m³)
        bubble_density: T,
        /// Initial bubble radius (m)
        initial_radius: T,
    },
    /// Zwart-Gerber-Belamri model (2004)
    ZGB {
        /// Nucleation site volume fraction
        nucleation_fraction: T,
        /// Bubble radius (m)
        bubble_radius: T,
        /// Vaporization coefficient
        f_vap: T,
        /// Condensation coefficient
        f_cond: T,
    },
}

impl<T: RealField + FromPrimitive + Copy> CavitationModel<T> {
    /// Calculate mass transfer rate (kg/m³/s)
    pub fn mass_transfer_rate(
        &self,
        pressure: T,
        vapor_pressure: T,
        void_fraction: T,
        density_liquid: T,
        density_vapor: T,
    ) -> T {
        match self {
            CavitationModel::Kunz {
                vaporization_coeff,
                condensation_coeff,
            } => Self::kunz_model(
                pressure,
                vapor_pressure,
                void_fraction,
                density_liquid,
                density_vapor,
                *vaporization_coeff,
                *condensation_coeff,
            ),

            CavitationModel::SchnerrSauer {
                bubble_density,
                initial_radius: _,
            } => Self::schnerr_sauer_model(
                pressure,
                vapor_pressure,
                void_fraction,
                density_liquid,
                density_vapor,
                *bubble_density,
            ),

            CavitationModel::ZGB {
                nucleation_fraction,
                bubble_radius,
                f_vap,
                f_cond,
            } => Self::zgb_model(ZgbParams {
                pressure,
                vapor_pressure,
                void_fraction,
                density_liquid,
                density_vapor,
                nucleation_fraction: *nucleation_fraction,
                bubble_radius: *bubble_radius,
                f_vap: *f_vap,
                f_cond: *f_cond,
            }),
        }
    }

    fn kunz_model(
        pressure: T,
        vapor_pressure: T,
        void_fraction: T,
        density_liquid: T,
        density_vapor: T,
        vaporization_coeff: T,
        condensation_coeff: T,
    ) -> T {
        let pressure_diff = pressure - vapor_pressure;
        let half = T::from_f64(0.5).unwrap_or_else(|| T::one() / (T::one() + T::one()));

        if pressure_diff < T::zero() {
            // Vaporization
            vaporization_coeff * density_vapor * (T::one() - void_fraction) * pressure_diff.abs()
                / (half * density_liquid)
        } else {
            // Condensation
            let rate = condensation_coeff * density_vapor * void_fraction * pressure_diff
                / (half * density_liquid);
            -rate
        }
    }

    fn schnerr_sauer_model(
        pressure: T,
        vapor_pressure: T,
        void_fraction: T,
        density_liquid: T,
        density_vapor: T,
        bubble_density: T,
    ) -> T {
        let three = T::from_f64(3.0).unwrap_or_else(|| T::one() + T::one() + T::one());
        let four_pi = T::from_f64(4.0 * std::f64::consts::PI).unwrap_or_else(|| T::one());
        let two_thirds = T::from_f64(2.0 / 3.0).unwrap_or_else(|| T::one());
        let one_third = T::from_f64(1.0 / 3.0).unwrap_or_else(|| T::one());

        let alpha = void_fraction.max(T::zero()).min(T::one());
        let n_b = bubble_density;

        // Calculate bubble radius from void fraction
        // R_B = [(3α)/(4πn_B(1-α))]^(1/3)
        let denominator = four_pi * n_b * (T::one() - alpha);
        let epsilon = T::from_f64(1e-12).unwrap_or_else(|| T::one());

        if denominator > epsilon && alpha > epsilon {
            let radius_cubed = three * alpha / denominator;
            let radius = radius_cubed.powf(one_third);

            if radius > epsilon {
                // Mass transfer rate
                let pressure_diff = pressure - vapor_pressure;
                let sign = if pressure_diff < T::zero() {
                    T::one()
                } else {
                    -T::one()
                };

                sign * three
                    * alpha
                    * (T::one() - alpha)
                    * density_vapor
                    * (two_thirds * pressure_diff.abs() / density_liquid).sqrt()
                    / radius
            } else {
                T::zero()
            }
        } else {
            T::zero()
        }
    }

    fn zgb_model(params: ZgbParams<T>) -> T {
        let pressure_diff = params.pressure - params.vapor_pressure;
        let three = T::from_f64(3.0).unwrap_or_else(|| T::one() + T::one() + T::one());
        let two_thirds = T::from_f64(2.0 / 3.0).unwrap_or_else(|| T::one());
        let epsilon = T::from_f64(1e-12).unwrap_or_else(|| T::one());

        if params.bubble_radius < epsilon {
            return T::zero();
        }

        if pressure_diff < T::zero() {
            // Vaporization
            params.f_vap
                * three
                * params.nucleation_fraction
                * (T::one() - params.void_fraction.max(T::zero()).min(T::one()))
                * params.density_vapor
                * (two_thirds * pressure_diff.abs() / params.density_liquid).sqrt()
                / params.bubble_radius
        } else {
            // Condensation
            let rate = params.f_cond
                * three
                * params.void_fraction.max(T::zero()).min(T::one())
                * params.density_vapor
                * (two_thirds * pressure_diff.max(T::zero()) / params.density_liquid).sqrt()
                / params.bubble_radius;
            -rate
        }
    }

    /// Get default parameters for each model type
    #[must_use]
    pub fn default_parameters(model_type: &str) -> Option<Self> {
        match model_type {
            "kunz" => Some(CavitationModel::Kunz {
                vaporization_coeff: T::from_f64(100.0)?,
                condensation_coeff: T::from_f64(100.0)?,
            }),
            "schnerr_sauer" => Some(CavitationModel::SchnerrSauer {
                bubble_density: T::from_f64(1e13)?, // #/m³
                initial_radius: T::from_f64(1e-6)?, // m
            }),
            "zgb" => Some(CavitationModel::ZGB {
                nucleation_fraction: T::from_f64(5e-4)?,
                bubble_radius: T::from_f64(1e-6)?, // m
                f_vap: T::from_f64(50.0)?,
                f_cond: T::from_f64(0.01)?,
            }),
            _ => None,
        }
    }
}
