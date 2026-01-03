//! Channel components for microfluidic networks

use super::{constants, Component};
use crate::resistance::models::ResistanceModel;
use cfd_core::error::Result;
use cfd_core::fluid::Fluid;
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
        let two = T::from_f64(2.0).unwrap_or_else(T::zero);
        two * self.area() / (self.width + self.height)
    }

    /// Get aspect ratio (width/height)
    pub fn aspect_ratio(&self) -> T {
        self.width / self.height
    }

    /// Calculate friction factor for laminar flow
    #[allow(dead_code)]
    fn friction_factor_laminar(&self) -> T {
        let alpha = self.aspect_ratio();
        let one = T::one();
        let c1 = T::from_f64(constants::RECT_CHANNEL_C1).unwrap_or_else(T::zero);
        let c2 = T::from_f64(constants::RECT_CHANNEL_C2).unwrap_or_else(T::zero);

        if alpha >= one {
            // Wide channel approximation
            c1
        } else {
            // Tall channel (alpha < 1)
            let inv_alpha = one / alpha;
            c2 / inv_alpha
        }
    }
}

impl<T: RealField + Copy + FromPrimitive> Component<T> for RectangularChannel<T> {
    fn resistance(&self, fluid: &Fluid<T>) -> T {
        // Use the validated RectangularChannelModel for consistency
        let model = crate::resistance::models::RectangularChannelModel {
            width: self.width,
            height: self.height,
            length: self.length,
        };
        
        // FlowConditions with zero flow rate to get purely laminar R
        let conditions = crate::resistance::models::FlowConditions::new(T::zero());
        
        model.calculate_resistance(fluid, &conditions).unwrap_or_else(|_| {
            // Fallback for safety, though model should be valid
            T::from_f64(1e12).unwrap_or_else(T::zero)
        })
    }

    fn coefficients(&self, fluid: &Fluid<T>) -> (T, T) {
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
        let pi = T::from_f64(std::f64::consts::PI).unwrap_or_else(T::zero);
        pi * self.diameter * self.diameter / T::from_f64(4.0).unwrap_or_else(T::zero)
    }

    /// Get hydraulic diameter (equals diameter for circular channels)
    pub fn hydraulic_diameter(&self) -> T {
        self.diameter
    }
}

impl<T: RealField + Copy + FromPrimitive> Component<T> for CircularChannel<T> {
    fn resistance(&self, fluid: &Fluid<T>) -> T {
        // Use Hagen-Poiseuille model for laminar circular flow
        let model = crate::resistance::models::HagenPoiseuilleModel {
            diameter: self.diameter,
            length: self.length,
        };
        let conditions = crate::resistance::models::FlowConditions::new(T::zero());
        model.calculate_resistance(fluid, &conditions).unwrap_or_else(|_| {
            T::from_f64(1e12).unwrap_or_else(T::zero)
        })
    }

    fn coefficients(&self, fluid: &Fluid<T>) -> (T, T) {
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
