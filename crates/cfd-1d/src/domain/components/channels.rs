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

use super::Component;
use crate::physics::resistance::models::ResistanceModel;
use cfd_core::error::Result;
use cfd_core::physics::fluid::ConstantPropertyFluid;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Rectangular microchannel component
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RectangularChannel<T: RealField + Copy> {
    /// Channel length [m]
    pub length: T,
    /// Channel width [m]
    pub width: T,
    /// Channel height [m]
    pub height: T,
    /// Surface roughness [m]
    pub roughness: T,
    /// Additional parameters
    pub parameters: HashMap<String, T>,
}

impl<T: RealField + Copy + FromPrimitive> RectangularChannel<T> {
    /// Create a new rectangular channel
    pub fn new(length: T, width: T, height: T, roughness: T) -> Self {
        Self {
            length,
            width,
            height,
            roughness,
            parameters: HashMap::new(),
        }
    }

    /// Create a square channel
    pub fn square(length: T, side: T, roughness: T) -> Self {
        Self::new(length, side, side, roughness)
    }

    /// Get cross-sectional area
    pub fn area(&self) -> T {
        self.width * self.height
    }

    /// Get hydraulic diameter
    pub fn hydraulic_diameter(&self) -> T {
        let two = T::one() + T::one();
        two * self.area() / (self.width + self.height)
    }

    /// Get aspect ratio (width/height)
    pub fn aspect_ratio(&self) -> T {
        self.width / self.height
    }
}

impl<T: RealField + Copy + FromPrimitive> Component<T> for RectangularChannel<T> {
    fn resistance(&self, fluid: &ConstantPropertyFluid<T>) -> T {
        // Use the validated RectangularChannelModel for consistency
        let model = crate::physics::resistance::models::RectangularChannelModel {
            width: self.width,
            height: self.height,
            length: self.length,
        };

        // FlowConditions with zero flow rate to get purely laminar R
        let conditions = crate::physics::resistance::models::FlowConditions::new(T::zero());

        model
            .calculate_resistance(fluid, &conditions)
            .unwrap_or_else(|_| {
                // Fallback for safety, though model should be valid
                T::from_f64(1e12).expect("Mathematical constant conversion compromised")
            })
    }

    fn coefficients(&self, fluid: &ConstantPropertyFluid<T>) -> (T, T) {
        let r = self.resistance(fluid);
        (r, T::zero())
    }

    fn component_type(&self) -> &'static str {
        "RectangularChannel"
    }

    fn parameters(&self) -> &HashMap<String, T> {
        &self.parameters
    }

    fn set_parameter(&mut self, key: &str, value: T) -> Result<()> {
        match key {
            "length" => self.length = value,
            "width" => self.width = value,
            "height" => self.height = value,
            "roughness" => self.roughness = value,
            _ => {
                self.parameters.insert(key.to_string(), value);
            }
        }
        Ok(())
    }

    fn volume(&self) -> Option<T> {
        Some(self.length * self.area())
    }
}

/// Circular microchannel component
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CircularChannel<T: RealField + Copy> {
    /// Channel length [m]
    pub length: T,
    /// Channel diameter [m]
    pub diameter: T,
    /// Surface roughness [m]
    pub roughness: T,
    /// Additional parameters
    pub parameters: HashMap<String, T>,
}

impl<T: RealField + Copy + FromPrimitive> CircularChannel<T> {
    /// Create a new circular channel
    pub fn new(length: T, diameter: T, roughness: T) -> Self {
        Self {
            length,
            diameter,
            roughness,
            parameters: HashMap::new(),
        }
    }

    /// Get cross-sectional area
    pub fn area(&self) -> T {
        let pi = T::pi();
        pi * self.diameter * self.diameter / (T::one() + T::one() + T::one() + T::one())
    }

    /// Get hydraulic diameter (equals diameter for circular channels)
    pub fn hydraulic_diameter(&self) -> T {
        self.diameter
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use cfd_core::physics::fluid::database::water_20c;

    #[test]
    fn test_rectangular_channel_resistance_matches_formula() {
        let fluid = water_20c::<f64>().unwrap();
        let width = 1e-3;
        let height = 1e-3;
        let length = 0.1;
        let chan = RectangularChannel::new(length, width, height, 0.0);
        let r = chan.resistance(&fluid);

        // For square channel, Shah-London: Po = 56.91, R = Po*mu*L/(2*A*Dh^2)
        let mu = fluid.viscosity;
        let area = width * height;
        let dh = 2.0 * width * height / (width + height);
        let r_expected = 56.91 * mu * length / (2.0 * area * dh * dh);
        assert_relative_eq!(r, r_expected, max_relative = 0.005);
    }

    #[test]
    fn test_rectangular_channel_hydraulic_diameter() {
        let chan = RectangularChannel::<f64>::new(0.1, 2e-3, 1e-3, 0.0);
        let expected = 2.0 * 2e-3 * 1e-3 / (2e-3 + 1e-3);
        assert_relative_eq!(chan.hydraulic_diameter(), expected, epsilon = 1e-15);
    }

    #[test]
    fn test_rectangular_channel_volume() {
        let chan = RectangularChannel::<f64>::new(0.05, 1e-3, 5e-4, 0.0);
        let expected = 0.05 * 1e-3 * 5e-4;
        assert_relative_eq!(chan.volume().unwrap(), expected, epsilon = 1e-20);
    }

    #[test]
    fn test_circular_channel_resistance_matches_hagen_poiseuille() {
        let fluid = water_20c::<f64>().unwrap();
        let d = 1e-3;
        let l = 0.1;
        let chan = CircularChannel::new(l, d, 0.0);
        let r = chan.resistance(&fluid);

        // Hagen-Poiseuille: R = 128 mu L / (pi D^4)
        let expected = 128.0 * fluid.viscosity * l / (std::f64::consts::PI * d.powi(4));
        assert_relative_eq!(r, expected, epsilon = 1e-6);
    }

    #[test]
    fn test_circular_channel_hydraulic_diameter_equals_diameter() {
        let chan = CircularChannel::<f64>::new(0.1, 2e-3, 0.0);
        assert_relative_eq!(chan.hydraulic_diameter(), 2e-3, epsilon = 1e-15);
    }

    #[test]
    fn test_circular_channel_volume() {
        let d = 1e-3_f64;
        let l = 0.1_f64;
        let chan = CircularChannel::<f64>::new(l, d, 0.0);
        let expected = l * std::f64::consts::PI * d * d / 4.0;
        assert_relative_eq!(chan.volume().unwrap(), expected, epsilon = 1e-15);
    }
}

impl<T: RealField + Copy + FromPrimitive> Component<T> for CircularChannel<T> {
    fn resistance(&self, fluid: &ConstantPropertyFluid<T>) -> T {
        // Use Hagen-Poiseuille model for laminar circular flow
        let model = crate::physics::resistance::models::HagenPoiseuilleModel {
            diameter: self.diameter,
            length: self.length,
        };
        let conditions = crate::physics::resistance::models::FlowConditions::new(T::zero());
        model
            .calculate_resistance(fluid, &conditions)
            .unwrap_or_else(|_| {
                T::from_f64(1e12).expect("Mathematical constant conversion compromised")
            })
    }

    fn coefficients(&self, fluid: &ConstantPropertyFluid<T>) -> (T, T) {
        // By default, components return laminar resistance.
        // For turbulent flow, one should use the ResistanceCalculator or a model
        // that supports flow-dependent coefficients.
        (self.resistance(fluid), T::zero())
    }

    fn component_type(&self) -> &'static str {
        "CircularChannel"
    }

    fn parameters(&self) -> &HashMap<String, T> {
        &self.parameters
    }

    fn set_parameter(&mut self, key: &str, value: T) -> Result<()> {
        match key {
            "length" => self.length = value,
            "diameter" => self.diameter = value,
            "roughness" => self.roughness = value,
            _ => {
                self.parameters.insert(key.to_string(), value);
            }
        }
        Ok(())
    }

    fn volume(&self) -> Option<T> {
        Some(self.length * self.area())
    }
}
