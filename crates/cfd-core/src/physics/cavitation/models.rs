//! Cavitation models for mass transfer rate calculations.

use aequitas::systems::si::quantities::{
    Dimensionless, Length, MassDensity, MassDensityRate, NumberDensity, Pressure,
};
use eunomia::{FloatElement, NumericElement};
use serde::{Deserialize, Serialize};

/// Parameters for ZGB cavitation model calculation
#[derive(Debug, Clone, Copy)]
pub struct ZgbParams<T: FloatElement + Copy> {
    /// Pressure field value
    pub pressure: Pressure<T>,
    /// Vapor pressure
    pub vapor_pressure: Pressure<T>,
    /// Void fraction
    pub void_fraction: T,
    /// Liquid density
    pub density_liquid: MassDensity<T>,
    /// Vapor density
    pub density_vapor: MassDensity<T>,
    /// Nucleation site volume fraction
    pub nucleation_fraction: Dimensionless<T>,
    /// Bubble radius
    pub bubble_radius: Length<T>,
    /// Vaporization coefficient
    pub f_vap: T,
    /// Condensation coefficient
    pub f_cond: T,
}

/// Cavitation model types
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum CavitationModel<T: FloatElement + Copy> {
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
        bubble_density: NumberDensity<T>,
        /// Initial bubble radius (m)
        initial_radius: Length<T>,
    },
    /// Zwart-Gerber-Belamri model (2004)
    ZGB {
        /// Nucleation site volume fraction
        nucleation_fraction: Dimensionless<T>,
        /// Bubble radius (m)
        bubble_radius: Length<T>,
        /// Vaporization coefficient
        f_vap: T,
        /// Condensation coefficient
        f_cond: T,
    },
}

impl<T: FloatElement + Copy> CavitationModel<T> {
    /// Calculate mass transfer rate (kg/m³/s)
    pub fn mass_transfer_rate(
        &self,
        pressure: Pressure<T>,
        vapor_pressure: Pressure<T>,
        void_fraction: T,
        density_liquid: MassDensity<T>,
        density_vapor: MassDensity<T>,
    ) -> MassDensityRate<T> {
        let pressure = pressure.into_base();
        let vapor_pressure = vapor_pressure.into_base();
        let density_liquid = density_liquid.into_base();
        let density_vapor = density_vapor.into_base();
        let rate = match self {
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
                bubble_density.into_base(),
            ),

            CavitationModel::ZGB {
                nucleation_fraction,
                bubble_radius,
                f_vap,
                f_cond,
            } => Self::zgb_model(ZgbParams {
                pressure: Pressure::from_base(pressure),
                vapor_pressure: Pressure::from_base(vapor_pressure),
                void_fraction,
                density_liquid: MassDensity::from_base(density_liquid),
                density_vapor: MassDensity::from_base(density_vapor),
                nucleation_fraction: *nucleation_fraction,
                bubble_radius: *bubble_radius,
                f_vap: *f_vap,
                f_cond: *f_cond,
            }),
        };
        MassDensityRate::from_base(rate)
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
        let half = <T as FloatElement>::from_f64(0.5);

        if pressure_diff < <T as NumericElement>::ZERO {
            // Vaporization
            vaporization_coeff
                * density_vapor
                * (<T as NumericElement>::ONE - void_fraction)
                * <T as NumericElement>::abs(pressure_diff)
                / (half * density_liquid)
        } else {
            // Condensation
            let rate = condensation_coeff * density_vapor * void_fraction * pressure_diff
                / (half * density_liquid);
            <T as NumericElement>::ZERO - rate
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
        let three = <T as FloatElement>::from_f64(3.0);
        let four_pi = <T as FloatElement>::from_f64(4.0 * std::f64::consts::PI);
        let two_thirds = <T as FloatElement>::from_f64(2.0 / 3.0);
        let one_third = <T as FloatElement>::from_f64(1.0 / 3.0);

        let alpha = void_fraction
            .max_scalar(<T as NumericElement>::ZERO)
            .min_scalar(<T as NumericElement>::ONE);
        let n_b = bubble_density;

        // Calculate bubble radius from void fraction
        // R_B = [(3α)/(4πn_B(1-α))]^(1/3)
        let denominator = four_pi * n_b * (<T as NumericElement>::ONE - alpha);
        let epsilon = <T as FloatElement>::from_f64(1e-12);

        if denominator > epsilon && alpha > epsilon {
            let radius_cubed = three * alpha / denominator;
            let radius = <T as FloatElement>::powf(radius_cubed, one_third);

            if radius > epsilon {
                // Mass transfer rate
                let pressure_diff = pressure - vapor_pressure;
                let sign = if pressure_diff < <T as NumericElement>::ZERO {
                    <T as NumericElement>::ONE
                } else {
                    <T as NumericElement>::ZERO - <T as NumericElement>::ONE
                };

                sign * three
                    * alpha
                    * (<T as NumericElement>::ONE - alpha)
                    * density_vapor
                    * <T as NumericElement>::sqrt(
                        two_thirds * <T as NumericElement>::abs(pressure_diff) / density_liquid,
                    )
                    / radius
            } else {
                <T as NumericElement>::ZERO
            }
        } else {
            <T as NumericElement>::ZERO
        }
    }

    fn zgb_model(params: ZgbParams<T>) -> T {
        let pressure_diff = params.pressure.into_base() - params.vapor_pressure.into_base();
        let three = <T as FloatElement>::from_f64(3.0);
        let two_thirds = <T as FloatElement>::from_f64(2.0 / 3.0);
        let epsilon = <T as FloatElement>::from_f64(1e-12);

        if params.bubble_radius.into_base() < epsilon {
            return <T as NumericElement>::ZERO;
        }

        if pressure_diff < <T as NumericElement>::ZERO {
            // Vaporization
            params.f_vap
                * three
                * params.nucleation_fraction.into_base()
                * (<T as NumericElement>::ONE
                    - params
                        .void_fraction
                        .max_scalar(<T as NumericElement>::ZERO)
                        .min_scalar(<T as NumericElement>::ONE))
                * params.density_vapor.into_base()
                * <T as NumericElement>::sqrt(
                    two_thirds * <T as NumericElement>::abs(pressure_diff)
                        / params.density_liquid.into_base(),
                )
                / params.bubble_radius.into_base()
        } else {
            // Condensation
            let rate = params.f_cond
                * three
                * params
                    .void_fraction
                    .max_scalar(<T as NumericElement>::ZERO)
                    .min_scalar(<T as NumericElement>::ONE)
                * params.density_vapor.into_base()
                * <T as NumericElement>::sqrt(
                    two_thirds * pressure_diff.max_scalar(<T as NumericElement>::ZERO)
                        / params.density_liquid.into_base(),
                )
                / params.bubble_radius.into_base();
            <T as NumericElement>::ZERO - rate
        }
    }

    /// Get default parameters for each model type
    #[must_use]
    pub fn default_parameters(model_type: &str) -> Option<Self> {
        match model_type {
            "kunz" => Some(CavitationModel::Kunz {
                vaporization_coeff: <T as FloatElement>::from_f64(100.0),
                condensation_coeff: <T as FloatElement>::from_f64(100.0),
            }),
            "schnerr_sauer" => Some(CavitationModel::SchnerrSauer {
                bubble_density: NumberDensity::from_base(<T as FloatElement>::from_f64(1e13)),
                initial_radius: Length::from_base(<T as FloatElement>::from_f64(1e-6)),
            }),
            "zgb" => Some(CavitationModel::ZGB {
                nucleation_fraction: Dimensionless::from_base(<T as FloatElement>::from_f64(5e-4)),
                bubble_radius: Length::from_base(<T as FloatElement>::from_f64(1e-6)),
                f_vap: <T as FloatElement>::from_f64(50.0),
                f_cond: <T as FloatElement>::from_f64(0.01),
            }),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::CavitationModel;
    use aequitas::systems::si::quantities::{
        Dimensionless, Length, MassDensity, NumberDensity, Pressure,
    };

    #[test]
    fn kunz_model_vaporizes_below_vapor_pressure_and_condenses_above() {
        let model = CavitationModel::<f64>::Kunz {
            vaporization_coeff: 2.0,
            condensation_coeff: 3.0,
        };

        let vaporization = model
            .mass_transfer_rate(
                Pressure::from_base(90.0),
                Pressure::from_base(100.0),
                0.25,
                MassDensity::from_base(1_000.0),
                MassDensity::from_base(1.0),
            )
            .into_base();
        let condensation = model
            .mass_transfer_rate(
                Pressure::from_base(110.0),
                Pressure::from_base(100.0),
                0.25,
                MassDensity::from_base(1_000.0),
                MassDensity::from_base(1.0),
            )
            .into_base();

        assert!((vaporization - 0.03).abs() <= 1.0e-14);
        assert!((condensation + 0.015).abs() <= 1.0e-14);
    }

    #[test]
    fn schnerr_sauer_zero_bubble_density_has_zero_rate() {
        let model = CavitationModel::<f64>::SchnerrSauer {
            bubble_density: NumberDensity::from_base(0.0),
            initial_radius: Length::from_base(1.0e-6),
        };

        assert_eq!(
            model
                .mass_transfer_rate(
                    Pressure::from_base(90.0),
                    Pressure::from_base(100.0),
                    0.5,
                    MassDensity::from_base(1_000.0),
                    MassDensity::from_base(1.0),
                )
                .into_base(),
            0.0
        );
    }

    #[test]
    fn zgb_rejects_degenerate_bubble_radius() {
        let model = CavitationModel::<f64>::ZGB {
            nucleation_fraction: Dimensionless::from_base(5.0e-4),
            bubble_radius: Length::from_base(1.0e-14),
            f_vap: 50.0,
            f_cond: 0.01,
        };

        assert_eq!(
            model
                .mass_transfer_rate(
                    Pressure::from_base(90.0),
                    Pressure::from_base(100.0),
                    0.5,
                    MassDensity::from_base(1_000.0),
                    MassDensity::from_base(1.0),
                )
                .into_base(),
            0.0
        );
    }

    #[test]
    fn default_parameters_construct_known_models() {
        assert!(matches!(
            CavitationModel::<f64>::default_parameters("kunz"),
            Some(CavitationModel::Kunz { .. })
        ));
        assert!(matches!(
            CavitationModel::<f64>::default_parameters("schnerr_sauer"),
            Some(CavitationModel::SchnerrSauer { .. })
        ));
        assert!(matches!(
            CavitationModel::<f64>::default_parameters("zgb"),
            Some(CavitationModel::ZGB { .. })
        ));
        assert!(CavitationModel::<f64>::default_parameters("unknown").is_none());
    }
}
