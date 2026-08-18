//! Channel components for microfluidic networks
//!
//! # Hagen-Poiseuille Theorem (Circular Channel)
//!
//! For fully-developed laminar flow in a circular pipe (Re < 2300):
//! ```text
//! R = 128 μ L / (π D⁴)   [Pa·s/m³]
//! ```
//! This is the exact analytical solution to the Stokes equations with no-slip BC.
//! Derivation: integrate the parabolic velocity profile `u(r)` over the cross-section
//! to obtain Q, then ΔP = R·Q where R = 8μL/(πr⁴) = 128μL/(πD⁴).
//!
//! # Generalised Hydraulic Diameter
//!
//! For non-circular sections, the hydraulic diameter is:
//! ```text
//! D_h = 4 A / P_wet
//! ```
//! where A is the cross-sectional area and P_wet is the wetted perimeter.
//! For a rectangle (W×H): `D_h = 2WH/(W+H)`.
//! For a circle: `D_h = D` (exact).
//!
//! # Rectangular Shah-London Resistance
//!
//! For rectangular ducts the laminar resistance uses the Shah-London (1978) Poiseuille
//! number Po(α) = f·Re, which depends on the aspect ratio α = min(W,H)/max(W,H):
//! ```text
//! R = Po(α) · μ · L / (2 A · D_h²)   [Pa·s/m³]
//! ```
//! The `RectangularChannelModel` in `resistance::models` computes Po(α) accurately.

use super::{real_from_f64, Component};
use crate::physics::resistance::models::ResistanceModel;
use aequitas::systems::si::quantities::{Area, Length, Volume};
use cfd_core::conversion::SafeFromF64;
use cfd_core::error::Result;
use cfd_core::physics::fluid::ConstantPropertyFluid;
use cfd_core::CfdScalar;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Rectangular microchannel component
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RectangularChannel<T: CfdScalar + Copy> {
    /// Channel length \[m]
    pub length: Length<T>,
    /// Channel width \[m]
    pub width: Length<T>,
    /// Channel height \[m]
    pub height: Length<T>,
    /// Surface roughness \[m]
    pub roughness: Length<T>,
    /// Additional parameters
    pub parameters: HashMap<String, T>,
}

impl<T: CfdScalar + Copy + SafeFromF64> RectangularChannel<T> {
    /// Create a new rectangular channel
    pub fn new(
        length: Length<T>,
        width: Length<T>,
        height: Length<T>,
        roughness: Length<T>,
    ) -> Self {
        Self {
            length,
            width,
            height,
            roughness,
            parameters: HashMap::new(),
        }
    }

    /// Create a square channel
    pub fn square(length: Length<T>, side: Length<T>, roughness: Length<T>) -> Self {
        Self::new(length, side, side, roughness)
    }

    /// Get cross-sectional area
    pub fn area(&self) -> Area<T> {
        Area::from_base(self.width.into_base() * self.height.into_base())
    }

    /// Get hydraulic diameter
    pub fn hydraulic_diameter(&self) -> Length<T> {
        let two = T::ONE + T::ONE;
        Length::from_base(
            two * self.area().into_base() / (self.width.into_base() + self.height.into_base()),
        )
    }

    /// Get aspect ratio (width/height)
    pub fn aspect_ratio(&self) -> T {
        self.width.into_base() / self.height.into_base()
    }
}

impl<T: CfdScalar + Copy + SafeFromF64> Component<T> for RectangularChannel<T> {
    fn resistance(&self, fluid: &ConstantPropertyFluid<T>) -> T {
        // Use the validated RectangularChannelModel for consistency
        let model = crate::physics::resistance::models::RectangularChannelModel {
            width: self.width.into_base(),
            height: self.height.into_base(),
            length: self.length.into_base(),
        };

        // FlowConditions with zero flow rate to get purely laminar R
        let conditions = crate::physics::resistance::models::FlowConditions::new(T::ZERO);

        model
            .calculate_resistance(fluid, &conditions)
            .unwrap_or_else(|_| {
                // Fallback for safety, though model should be valid
                real_from_f64(1e12)
            })
    }

    fn coefficients(&self, fluid: &ConstantPropertyFluid<T>) -> (T, T) {
        let r = self.resistance(fluid);
        (r, T::ZERO)
    }

    fn component_type(&self) -> &'static str {
        "RectangularChannel"
    }

    fn parameters(&self) -> &HashMap<String, T> {
        &self.parameters
    }

    fn set_parameter(&mut self, key: &str, value: T) -> Result<()> {
        match key {
            "length" => self.length = Length::from_base(value),
            "width" => self.width = Length::from_base(value),
            "height" => self.height = Length::from_base(value),
            "roughness" => self.roughness = Length::from_base(value),
            _ => {
                self.parameters.insert(key.to_string(), value);
            }
        }
        Ok(())
    }

    fn volume(&self) -> Option<Volume<T>> {
        Some(Volume::from_base(
            self.length.into_base() * self.area().into_base(),
        ))
    }
}

/// Circular microchannel component
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CircularChannel<T: CfdScalar + Copy> {
    /// Channel length \[m]
    pub length: Length<T>,
    /// Channel diameter \[m]
    pub diameter: Length<T>,
    /// Surface roughness \[m]
    pub roughness: Length<T>,
    /// Additional parameters
    pub parameters: HashMap<String, T>,
}

impl<T: CfdScalar + Copy + SafeFromF64> CircularChannel<T> {
    /// Create a new circular channel
    pub fn new(length: Length<T>, diameter: Length<T>, roughness: Length<T>) -> Self {
        Self {
            length,
            diameter,
            roughness,
            parameters: HashMap::new(),
        }
    }

    /// Get cross-sectional area
    pub fn area(&self) -> Area<T> {
        let pi = T::pi();
        Area::from_base(
            pi * self.diameter.into_base() * self.diameter.into_base()
                / (T::ONE + T::ONE + T::ONE + T::ONE),
        )
    }

    /// Get hydraulic diameter (equals diameter for circular channels)
    pub fn hydraulic_diameter(&self) -> Length<T> {
        self.diameter
    }
}

impl<T: CfdScalar + Copy + SafeFromF64> Component<T> for CircularChannel<T> {
    fn resistance(&self, fluid: &ConstantPropertyFluid<T>) -> T {
        // Use Hagen-Poiseuille model for laminar circular flow
        let model = crate::physics::resistance::models::HagenPoiseuilleModel {
            diameter: self.diameter.into_base(),
            length: self.length.into_base(),
        };
        let conditions = crate::physics::resistance::models::FlowConditions::new(T::ZERO);
        model
            .calculate_resistance(fluid, &conditions)
            .unwrap_or_else(|_| real_from_f64(1e12))
    }

    fn coefficients(&self, fluid: &ConstantPropertyFluid<T>) -> (T, T) {
        // By default, components return laminar resistance.
        // For turbulent flow, one should use the ResistanceCalculator or a model
        // that supports flow-dependent coefficients.
        (self.resistance(fluid), T::ZERO)
    }

    fn component_type(&self) -> &'static str {
        "CircularChannel"
    }

    fn parameters(&self) -> &HashMap<String, T> {
        &self.parameters
    }

    fn set_parameter(&mut self, key: &str, value: T) -> Result<()> {
        match key {
            "length" => self.length = Length::from_base(value),
            "diameter" => self.diameter = Length::from_base(value),
            "roughness" => self.roughness = Length::from_base(value),
            _ => {
                self.parameters.insert(key.to_string(), value);
            }
        }
        Ok(())
    }

    fn volume(&self) -> Option<Volume<T>> {
        Some(Volume::from_base(
            self.length.into_base() * self.area().into_base(),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cfd_core::physics::fluid::database::water_20c;
    use eunomia::assert_relative_eq;

    #[test]
    fn test_rectangular_channel_resistance_matches_formula() {
        let fluid = water_20c::<f64>().expect("expected value");
        let width = 1e-3;
        let height = 1e-3;
        let length = 0.1;
        let chan = RectangularChannel::new(
            Length::from_base(length),
            Length::from_base(width),
            Length::from_base(height),
            Length::from_base(0.0),
        );
        let r = chan.resistance(&fluid);

        // For square channel, Bahrami exact fit yields slightly different Po (~56.4-56.9)
        // compared to the hardcoded Shah-London polynomial estimation of 56.91.
        let mu = fluid.viscosity.into_base();
        let area = width * height;
        let dh = 2.0 * width * height / (width + height);
        let r_expected = 56.91 * mu * length / (2.0 * area * dh * dh);
        assert_relative_eq!(r, r_expected, max_relative = 0.01);
    }

    #[test]
    fn test_rectangular_channel_hydraulic_diameter() {
        let chan = RectangularChannel::<f64>::new(
            Length::from_base(0.1),
            Length::from_base(2e-3),
            Length::from_base(1e-3),
            Length::from_base(0.0),
        );
        let expected = 2.0 * 2e-3 * 1e-3 / (2e-3 + 1e-3);
        assert_relative_eq!(
            chan.hydraulic_diameter().into_base(),
            expected,
            epsilon = 1e-15
        );
    }

    #[test]
    fn test_rectangular_channel_volume() {
        let chan = RectangularChannel::<f64>::new(
            Length::from_base(0.05),
            Length::from_base(1e-3),
            Length::from_base(5e-4),
            Length::from_base(0.0),
        );
        let expected = 0.05 * 1e-3 * 5e-4;
        assert_relative_eq!(
            chan.volume().expect("expected value").into_base(),
            expected,
            epsilon = 1e-20
        );
    }

    #[test]
    fn test_circular_channel_resistance_matches_hagen_poiseuille() {
        let fluid = water_20c::<f64>().expect("expected value");
        let d = 1e-3;
        let l = 0.1;
        let chan = CircularChannel::new(
            Length::from_base(l),
            Length::from_base(d),
            Length::from_base(0.0),
        );
        let r = chan.resistance(&fluid);

        // Hagen-Poiseuille: R = 128 mu L / (pi D^4)
        let expected = 128.0 * fluid.viscosity.into_base() * l / (std::f64::consts::PI * d.powi(4));
        assert_relative_eq!(r, expected, epsilon = 1e-6);
    }

    #[test]
    fn test_circular_channel_hydraulic_diameter_equals_diameter() {
        let chan = CircularChannel::<f64>::new(
            Length::from_base(0.1),
            Length::from_base(2e-3),
            Length::from_base(0.0),
        );
        assert_relative_eq!(chan.hydraulic_diameter().into_base(), 2e-3, epsilon = 1e-15);
    }

    #[test]
    fn test_circular_channel_volume() {
        let d = 1e-3_f64;
        let l = 0.1_f64;
        let chan = CircularChannel::<f64>::new(
            Length::from_base(l),
            Length::from_base(d),
            Length::from_base(0.0),
        );
        let expected = l * std::f64::consts::PI * d * d / 4.0;
        assert_relative_eq!(
            chan.volume().expect("expected value").into_base(),
            expected,
            epsilon = 1e-15
        );
    }
}
