//! Traits and common types for resistance models.

use cfd_core::error::Result;
use cfd_core::fluid::ConstantPropertyFluid;
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;

/// Trait for hydraulic resistance models
pub trait ResistanceModel<T: RealField + Copy> {
    /// Calculate hydraulic resistance [Pa·s/m³]
    fn calculate_resistance(&self, fluid: &Fluid<T>, conditions: &FlowConditions<T>) -> Result<T>;

    /// Get model name
    fn model_name(&self) -> &str;

    /// Get applicable Reynolds number range
    fn reynolds_range(&self) -> (T, T);

    /// Check if model is applicable for given conditions
    fn is_applicable(&self, conditions: &FlowConditions<T>) -> bool {
        let (re_min, re_max) = self.reynolds_range();
        if let Some(re) = &conditions.reynolds_number {
            *re >= re_min && *re <= re_max
        } else {
            // Cannot determine applicability without Reynolds number
            false
        }
    }
}

/// Flow conditions for resistance calculations
#[derive(Debug, Clone)]
pub struct FlowConditions<T: RealField + Copy> {
    /// Reynolds number
    pub reynolds_number: Option<T>,
    /// Flow velocity [m/s]
    pub velocity: Option<T>,
    /// Flow rate [m³/s]
    pub flow_rate: Option<T>,
    /// Temperature [K]
    pub temperature: T,
    /// Pressure [Pa]
    pub pressure: T,
}

impl<T: RealField + Copy + FromPrimitive> FlowConditions<T> {
    /// Create new flow conditions with default temperature and pressure
    pub fn new(velocity: T) -> Self {
        use cfd_core::constants::physics::thermo::{P_ATM, T_STANDARD};

        Self {
            reynolds_number: None,
            velocity: Some(velocity),
            flow_rate: None,
            temperature: T::from_f64(T_STANDARD).unwrap_or_else(T::one),
            pressure: T::from_f64(P_ATM).unwrap_or_else(T::one),
        }
    }
}
