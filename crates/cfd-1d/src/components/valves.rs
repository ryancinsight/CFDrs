//! Valve components for microfluidic networks

use super::Component;
use cfd_core::{Fluid, Result};
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
/// Valve type enumeration
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum ValveType {
    /// Normally open valve
    NormallyOpen,
    /// Normally closed valve
    NormallyClosed,
    /// Check valve (one-way)
    Check,
    /// Proportional control valve
    Proportional,
}
/// Microvalve component
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Microvalve<T: RealField + Copy> {
    /// Flow coefficient [m³/s/Pa^0.5]
    pub cv: T,
    /// Opening fraction (0=closed, 1=open)
    pub opening: T,
    /// Additional parameters
    pub parameters: HashMap<String, T>,
impl<T: RealField + Copy + FromPrimitive + Float> Microvalve<T> {
    /// Create a new microvalve
    pub fn new(cv: T) -> Self {
        Self {
            cv,
            opening: T::one(),
            parameters: HashMap::new(),
        }
    }
impl<T: RealField + Copy + FromPrimitive + Float> Component<T> for Microvalve<T> {
    fn resistance(&self, _fluid: &Fluid<T>) -> T {
        if self.opening <= T::zero() {
            // Closed valve - infinite resistance
            T::from_f64(1e12).unwrap_or_else(T::one)
        } else {
            // Resistance inversely proportional to opening
            let base_resistance = T::one() / (self.cv * self.cv);
            base_resistance / (self.opening * self.opening)
    fn component_type(&self) -> &str {
        "Microvalve"
    fn parameters(&self) -> &HashMap<String, T> {
        &self.parameters
    fn set_parameter(&mut self, key: &str, value: T) -> Result<()> {
        match key {
            "cv" => self.cv = value,
            "opening" => {
                self.opening = if value < T::zero() {
                    T::zero()
                } else if value > T::one() {
                    T::one()
                } else {
                    value
                };
            }
            _ => {
                self.parameters.insert(key.to_string(), value);
        Ok(())
    fn is_active(&self) -> bool {
        true
