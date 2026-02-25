//! Parameter managers for different domains
//!
//! This module provides domain-specific parameter managers that implement
//! the `ParameterManager` trait. Each manager handles parameters for a specific
//! aspect of the microfluidic design system.

mod arc;
mod collision;
mod geometry;
mod serpentine;
mod symmetry;

pub use self::arc::ArcParameterManager;
pub use self::collision::CollisionParameterManager;
pub use self::geometry::GeometryParameterManager;
pub use self::serpentine::SerpentineParameterManager;
pub use self::symmetry::SymmetryParameterManager;

use crate::state_management::{
    errors::ParameterResult,
    parameters::ParameterMetadata,
    validation::ValidationRuleSet,
};
use std::fmt::Debug;

/// Core trait for all parameter managers
pub trait ParameterManager: Debug + Send + Sync {
    /// Get parameter value by name
    fn get_parameter(&self, name: &str) -> ParameterResult<Box<dyn std::any::Any>>;

    /// Set parameter value by name
    fn set_parameter(
        &mut self,
        name: &str,
        value: Box<dyn std::any::Any>,
        reason: &str,
    ) -> ParameterResult<()>;

    /// Get all parameter names managed by this manager
    fn parameter_names(&self) -> Vec<String>;

    /// Check if a parameter exists
    fn has_parameter(&self, name: &str) -> bool;

    /// Validate all parameters
    fn validate_all(&self) -> ParameterResult<()>;

    /// Get parameter metadata
    fn get_metadata(&self, name: &str) -> ParameterResult<&ParameterMetadata>;

    /// Get manager domain name
    fn domain_name(&self) -> &str;

    /// Reset all parameters to defaults
    fn reset_all(&mut self, reason: &str) -> ParameterResult<()>;

    /// Get validation rules for this manager
    fn validation_rules(&self) -> &ValidationRuleSet;
}
