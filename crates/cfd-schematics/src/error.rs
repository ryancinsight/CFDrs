//! Domain-specific error types for the microfluidic schematic design library.
//!
//! All error types are centralized in `cfd_core::error`. This module provides
//! re-exports and convenience constructors for schematic-domain errors.

// Re-export all error types from cfd-core
pub use cfd_core::error::{
    AdaptationErrorKind, ConfigurationErrorKind, ConstraintErrorKind, DependencyErrorKind,
    Error, GeometryErrorKind, ParameterErrorKind, RegistryErrorKind, Result, StrategyErrorKind,
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

/// Convenience constructors for `GeometryErrorKind`
pub trait GeometryErrorExt {
    /// Create an invalid point error
    fn invalid_point(x: f64, y: f64) -> Self;
    /// Create an invalid box dimensions error
    fn invalid_box_dimensions(width: f64, height: f64) -> Self;
    /// Create an insufficient space error
    fn insufficient_space(required: f64, available: f64) -> Self;
}

impl GeometryErrorExt for Error {
    fn invalid_point(x: f64, y: f64) -> Self {
        Error::Geometry(GeometryErrorKind::InvalidPoint { x, y })
    }

    fn invalid_box_dimensions(width: f64, height: f64) -> Self {
        Error::Geometry(GeometryErrorKind::InvalidBoxDimensions { width, height })
    }

    fn insufficient_space(required: f64, available: f64) -> Self {
        Error::Geometry(GeometryErrorKind::InsufficientSpace { required, available })
    }
}

/// Convenience constructors for `ConfigurationErrorKind`
pub trait ConfigurationErrorExt {
    /// Create an invalid geometry config error
    fn invalid_geometry_config(field: &str, value: f64, constraint: &str) -> Self;
    /// Create an invalid serpentine config error
    fn invalid_serpentine_config(field: &str, value: f64, constraint: &str) -> Self;
    /// Create an invalid arc config error
    fn invalid_arc_config(field: &str, value: f64, constraint: &str) -> Self;
    /// Create an invalid frustum config error
    fn invalid_frustum_config(field: &str, value: f64, constraint: &str) -> Self;
    /// Create an invalid generation config error
    fn invalid_generation_config(field: &str, constraint: &str) -> Self;
}

impl ConfigurationErrorExt for Error {
    fn invalid_geometry_config(field: &str, value: f64, constraint: &str) -> Self {
        Error::Configuration(ConfigurationErrorKind::InvalidGeometryConfig {
            field: field.to_string(),
            value,
            constraint: constraint.to_string(),
        })
    }

    fn invalid_serpentine_config(field: &str, value: f64, constraint: &str) -> Self {
        Error::Configuration(ConfigurationErrorKind::InvalidSerpentineConfig {
            field: field.to_string(),
            value,
            constraint: constraint.to_string(),
        })
    }

    fn invalid_arc_config(field: &str, value: f64, constraint: &str) -> Self {
        Error::Configuration(ConfigurationErrorKind::InvalidArcConfig {
            field: field.to_string(),
            value,
            constraint: constraint.to_string(),
        })
    }

    fn invalid_frustum_config(field: &str, value: f64, constraint: &str) -> Self {
        Error::Configuration(ConfigurationErrorKind::InvalidFrustumConfig {
            field: field.to_string(),
            value,
            constraint: constraint.to_string(),
        })
    }

    fn invalid_generation_config(field: &str, constraint: &str) -> Self {
        Error::Configuration(ConfigurationErrorKind::InvalidGenerationConfig {
            field: field.to_string(),
            constraint: constraint.to_string(),
        })
    }
}

/// Convenience constructors for `VisualizationErrorKind`
pub trait VisualizationErrorExt {
    /// Create a rendering error
    fn rendering_error(message: &str) -> Self;
    /// Create an invalid output path error
    fn invalid_output_path(path: &str, reason: &str) -> Self;
    /// Create a file error
    fn file_error(message: &str) -> Self;
    /// Create an unsupported format error
    fn unsupported_format(format: &str, message: &str) -> Self;
}

impl VisualizationErrorExt for Error {
    fn rendering_error(message: &str) -> Self {
        Error::Visualization(VisualizationErrorKind::RenderingError {
            message: message.to_string(),
        })
    }

    fn invalid_output_path(path: &str, reason: &str) -> Self {
        Error::Visualization(VisualizationErrorKind::InvalidOutputPath {
            path: path.to_string(),
            reason: reason.to_string(),
        })
    }

    fn file_error(message: &str) -> Self {
        Error::Visualization(VisualizationErrorKind::FileError {
            message: message.to_string(),
        })
    }

    fn unsupported_format(format: &str, message: &str) -> Self {
        Error::Visualization(VisualizationErrorKind::UnsupportedFormat {
            format: format.to_string(),
            message: message.to_string(),
        })
    }
}

/// Convenience constructors for `StrategyErrorKind`
pub trait StrategyErrorExt {
    /// Create a strategy creation failed error
    fn strategy_creation_failed(channel_type: &str, reason: &str) -> Self;
    /// Create an execution failed error
    fn execution_failed(from: (f64, f64), to: (f64, f64), reason: &str) -> Self;
}

impl StrategyErrorExt for Error {
    fn strategy_creation_failed(channel_type: &str, reason: &str) -> Self {
        Error::Strategy(StrategyErrorKind::StrategyCreationFailed {
            channel_type: channel_type.to_string(),
            reason: reason.to_string(),
        })
    }

    fn execution_failed(from: (f64, f64), to: (f64, f64), reason: &str) -> Self {
        Error::Strategy(StrategyErrorKind::ExecutionFailed {
            from_x: from.0,
            from_y: from.1,
            to_x: to.0,
            to_y: to.1,
            reason: reason.to_string(),
        })
    }
}

/// Convenience constructors for `ParameterErrorKind`
pub trait ParameterErrorExt {
    /// Create a not-found error
    fn param_not_found(name: &str, domain: &str) -> Self;
    /// Create an invalid-value error
    fn param_invalid_value(name: &str, value: &dyn std::fmt::Debug, constraint: &str) -> Self;
    /// Create a type-mismatch error
    fn param_type_mismatch(name: &str, expected: &str, actual: &str) -> Self;
    /// Create a read-only error
    fn param_read_only(name: &str) -> Self;
}

impl ParameterErrorExt for Error {
    fn param_not_found(name: &str, domain: &str) -> Self {
        Error::Parameter(ParameterErrorKind::NotFound {
            name: name.to_string(),
            domain: domain.to_string(),
        })
    }

    fn param_invalid_value(name: &str, value: &dyn std::fmt::Debug, constraint: &str) -> Self {
        Error::Parameter(ParameterErrorKind::InvalidValue {
            name: name.to_string(),
            value: format!("{value:?}"),
            constraint: constraint.to_string(),
        })
    }

    fn param_type_mismatch(name: &str, expected: &str, actual: &str) -> Self {
        Error::Parameter(ParameterErrorKind::TypeMismatch {
            name: name.to_string(),
            expected: expected.to_string(),
            actual: actual.to_string(),
        })
    }

    fn param_read_only(name: &str) -> Self {
        Error::Parameter(ParameterErrorKind::ReadOnly {
            name: name.to_string(),
        })
    }
}
