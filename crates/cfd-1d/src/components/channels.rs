//! Channel components for microfluidic networks

use super::{constants, Component};
use cfd_core::{Fluid, Result};
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};
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
impl<T: RealField + Copy + FromPrimitive + Float> RectangularChannel<T> {
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
    /// Get cross-sectional area
    pub fn area(&self) -> T {
        self.width * self.height
    /// Get hydraulic diameter
    pub fn hydraulic_diameter(&self) -> T {
        let two = T::from_f64(2.0).unwrap_or_else(T::zero);
        two * self.area() / (self.width + self.height)
    /// Get aspect ratio (width/height)
    pub fn aspect_ratio(&self) -> T {
        self.width / self.height
    /// Calculate friction factor for laminar flow
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
impl<T: RealField + Copy + FromPrimitive + Float> Component<T> for RectangularChannel<T> {
    fn resistance(&self, fluid: &Fluid<T>) -> T {
        let dh = self.hydraulic_diameter();
        let area = self.area();
        let f = self.friction_factor_laminar();
        // Use default temperature of 20°C (293.15 K)
        let temperature = T::from_f64(293.15).unwrap_or_else(T::zero);
        let kinematic_viscosity = fluid
            .dynamic_viscosity(temperature)
            .unwrap_or_else(|_| T::from_f64(0.001).unwrap_or_else(T::one))
            / fluid.density;
        let resistance = f * self.length * kinematic_viscosity / (area * dh * dh);
        // Ensure positive resistance
        if resistance > T::zero() {
            resistance
            T::from_f64(1e-12).unwrap_or_else(T::zero)
    fn component_type(&self) -> &str {
        "RectangularChannel"
    fn parameters(&self) -> &HashMap<String, T> {
        &self.parameters
    fn set_parameter(&mut self, key: &str, value: T) -> Result<()> {
        match key {
            "length" => self.length = value,
            "width" => self.width = value,
            "height" => self.height = value,
            "roughness" => self.roughness = value,
            _ => {
                self.parameters.insert(key.to_string(), value);
            }
        Ok(())
    fn volume(&self) -> Option<T> {
        Some(self.length * self.area())
/// Circular microchannel component
pub struct CircularChannel<T: RealField + Copy> {
    /// Channel diameter [m]
    pub diameter: T,
impl<T: RealField + Copy + FromPrimitive + Float> CircularChannel<T> {
    /// Create a new circular channel
    pub fn new(length: T, diameter: T, roughness: T) -> Self {
            diameter,
        let pi = T::from_f64(std::f64::consts::PI).unwrap_or_else(T::zero);
        pi * self.diameter * self.diameter / T::from_f64(4.0).unwrap_or_else(T::zero)
    /// Get hydraulic diameter (equals diameter for circular channels)
        self.diameter
impl<T: RealField + Copy + FromPrimitive + Float> Component<T> for CircularChannel<T> {
        // Hagen-Poiseuille equation for circular channels
        // R = 128 * μ * L / (π * D^4)
        let c128 = T::from_f64(128.0).unwrap_or_else(T::zero);
        let viscosity = fluid
            .unwrap_or_else(|_| T::from_f64(0.001).unwrap_or_else(T::one));
        let d4 = self.diameter * self.diameter * self.diameter * self.diameter;
        c128 * viscosity * self.length / (pi * d4)
        "CircularChannel"
            "diameter" => self.diameter = value,
