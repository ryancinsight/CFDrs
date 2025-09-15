//! Cavitation models for mass transfer rate calculations.

use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

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
            } => self.kunz_model(
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
            } => self.schnerr_sauer_model(
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
            } => self.zgb_model(
                pressure,
                vapor_pressure,
                void_fraction,
                density_liquid,
                density_vapor,
                *nucleation_fraction,
                *bubble_radius,
                *f_vap,
                *f_cond,
            ),
        }
    }

    fn kunz_model(
        &self,
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
        &self,
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

        let alpha = void_fraction;
        let n_b = bubble_density;

        // Calculate bubble radius from void fraction
        // R_B = [(3α)/(4πn_B(1-α))]^(1/3)
        let denominator = four_pi * n_b * (T::one() - alpha);

        if denominator > T::from_f64(1e-10).unwrap_or_else(|| T::one()) {
            let radius_cubed = three * alpha / denominator;
            let radius = radius_cubed.powf(one_third);

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
    }

    fn zgb_model(
        &self,
        pressure: T,
        vapor_pressure: T,
        void_fraction: T,
        density_liquid: T,
        density_vapor: T,
        nucleation_fraction: T,
        bubble_radius: T,
        f_vap: T,
        f_cond: T,
    ) -> T {
        let pressure_diff = pressure - vapor_pressure;
        let three = T::from_f64(3.0).unwrap_or_else(|| T::one() + T::one() + T::one());
        let two_thirds = T::from_f64(2.0 / 3.0).unwrap_or_else(|| T::one());

        if pressure_diff < T::zero() {
            // Vaporization
            f_vap
                * three
                * nucleation_fraction
                * (T::one() - void_fraction)
                * density_vapor
                * (two_thirds * pressure_diff.abs() / density_liquid).sqrt()
                / bubble_radius
        } else {
            // Condensation
            let rate = f_cond
                * three
                * void_fraction
                * density_vapor
                * (two_thirds * pressure_diff / density_liquid).sqrt()
                / bubble_radius;
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
