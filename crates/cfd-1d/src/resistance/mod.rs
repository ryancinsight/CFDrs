//! Hydraulic resistance models for 1D CFD networks.
//!
//! This module provides comprehensive resistance modeling for various
//! microfluidic components and flow conditions, including analytical
//! solutions and empirical correlations.

pub mod calculator;
pub mod conditions;
pub mod factory;
pub mod geometry;
pub mod models;
// Re-export main types
pub use calculator::ResistanceCalculator;
pub use conditions::FlowConditions;
pub use factory::{CombinationMethod, ResistanceModelFactory};
pub use geometry::ChannelGeometry;
pub use models::{
    DarcyWeisbachModel, EntranceEffectsModel, HagenPoiseuilleModel, RectangularChannelModel,
};
use cfd_core::error::Result;
use cfd_core::fluid::Fluid;
use nalgebra::RealField;
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
            true // Assume applicable if Re is unknown
        }
    }
}
