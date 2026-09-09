//! Domain-specific error types for the microfluidic schematic design library.
//!
//! All error types are centralized in `cfd_core::error`. This module provides
//! re-exports and convenience constructors for schematic-domain errors.

// Re-export all error types from cfd-core
pub use cfd_core::error::{
    AdaptationErrorKind, ConfigurationErrorKind, ConstraintErrorKind, DependencyErrorKind, Error,
    GeometryErrorKind, ParameterErrorKind, RegistryErrorKind, Result, StrategyErrorKind,
    ValidationErrorKind, VisualizationErrorKind,
};

// Backward-compatible type aliases (old names → Kind types)
/// Geometry error type (alias for `GeometryErrorKind`)
pub type GeometryError = GeometryErrorKind;
/// Configuration error type (alias for `ConfigurationErrorKind`)
pub type ConfigurationError = ConfigurationErrorKind;
/// Visualization error type (alias for `VisualizationErrorKind`)
pub type VisualizationError = VisualizationErrorKind;
/// Strategy error type (alias for `StrategyErrorKind`)
pub type StrategyError = StrategyErrorKind;
/// Scheme error type (alias for `Error`)
pub type SchemeError = Error;

// Backward-compatible result type aliases using domain-specific error kinds
/// Result type for scheme operations
pub type SchemeResult<T> = std::result::Result<T, Error>;
/// Result type for geometry operations
pub type GeometryResult<T> = std::result::Result<T, GeometryErrorKind>;
/// Result type for configuration operations
pub type ConfigurationResult<T> = std::result::Result<T, ConfigurationErrorKind>;
/// Result type for visualization operations
pub type VisualizationResult<T> = std::result::Result<T, VisualizationErrorKind>;
/// Result type for strategy operations
pub type StrategyResult<T> = std::result::Result<T, StrategyErrorKind>;

// Convenience constructors are provided as inherent methods on each Kind type
// in cfd_core::error (e.g. GeometryErrorKind::invalid_box_dimensions(),
// ConfigurationErrorKind::invalid_arc_config(), etc.). Callers use them
// directly via the type aliases above.
