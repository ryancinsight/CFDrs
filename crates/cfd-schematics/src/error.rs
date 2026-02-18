//! error.rs - Domain-Specific Error Types
//!
//! This module provides comprehensive error types for the microfluidic schematic
//! design library. Each error type is designed to provide clear, actionable
//! information to help users diagnose and fix issues.

use crate::geometry::Point2D;
use thiserror::Error;

/// Main error type for the scheme library
///
/// This enum encompasses all possible errors that can occur during
/// microfluidic schematic design operations.
#[derive(Error, Debug)]
pub enum SchemeError {
    /// Errors related to geometry generation and validation
    #[error("Geometry error: {0}")]
    Geometry(#[from] GeometryError),

    /// Errors related to configuration validation
    #[error("Configuration error: {0}")]
    Configuration(#[from] ConfigurationError),

    /// Errors related to visualization and rendering
    #[error("Visualization error: {0}")]
    Visualization(#[from] VisualizationError),

    /// Errors related to channel type strategy operations
    #[error("Strategy error: {0}")]
    Strategy(#[from] StrategyError),

    /// Legacy simulation errors (kept for backward compatibility)
    #[error("Simulation error: {0}")]
    Simulation(#[from] SimulationError),
}

/// Errors related to geometry generation and validation
#[derive(Error, Debug)]
pub enum GeometryError {
    /// Invalid box dimensions provided
    #[error(
        "Invalid box dimensions: width={width}, height={height}. Both dimensions must be positive."
    )]
    InvalidBoxDimensions {
        /// Width of the box
        width: f64,
        /// Height of the box
        height: f64,
    },

    /// Invalid point coordinates
    #[error("Invalid point coordinates: ({x}, {y}). Coordinates must be finite numbers.")]
    InvalidPoint { x: f64, y: f64 },

    /// Insufficient space for channel generation
    #[error("Insufficient space for channel generation. Required space: {required}, available: {available}")]
    InsufficientSpace { required: f64, available: f64 },

    /// Invalid split pattern
    #[error("Invalid split pattern: {reason}")]
    InvalidSplitPattern { reason: String },

    /// Node creation failed
    #[error("Failed to create node at position ({x}, {y}): {reason}")]
    NodeCreationFailed { x: f64, y: f64, reason: String },

    /// Channel creation failed
    #[error("Failed to create channel from node {from_id} to node {to_id}: {reason}")]
    ChannelCreationFailed {
        from_id: usize,
        to_id: usize,
        reason: String,
    },

    /// Invalid channel path
    #[error("Invalid channel path: {reason}. Path must contain at least 2 points.")]
    InvalidChannelPath { reason: String },

    /// Overlapping channels detected
    #[error("Overlapping channels detected between points ({x1}, {y1}) and ({x2}, {y2})")]
    OverlappingChannels { x1: f64, y1: f64, x2: f64, y2: f64 },
}

/// Errors related to configuration validation
#[derive(Error, Debug)]
pub enum ConfigurationError {
    /// Invalid geometry configuration
    #[error("Invalid geometry configuration: {field} = {value}. {constraint}")]
    InvalidGeometryConfig {
        field: String,
        value: f64,
        constraint: String,
    },

    /// Invalid serpentine configuration
    #[error("Invalid serpentine configuration: {field} = {value}. {constraint}")]
    InvalidSerpentineConfig {
        field: String,
        value: f64,
        constraint: String,
    },

    /// Invalid arc configuration
    #[error("Invalid arc configuration: {field} = {value}. {constraint}")]
    InvalidArcConfig {
        field: String,
        value: f64,
        constraint: String,
    },

    /// Invalid frustum configuration
    #[error("Invalid frustum configuration: {field} = {value}. {constraint}")]
    InvalidFrustumConfig {
        field: String,
        value: f64,
        constraint: String,
    },

    /// Invalid geometry generation configuration
    #[error("Invalid geometry generation configuration: {field}. {constraint}")]
    InvalidGenerationConfig { field: String, constraint: String },

    /// Conflicting configuration values
    #[error("Conflicting configuration values: {conflict}")]
    ConflictingValues { conflict: String },

    /// Missing required configuration
    #[error("Missing required configuration: {field}")]
    MissingConfiguration { field: String },
}

/// Errors related to visualization and rendering
#[derive(Error, Debug)]
pub enum VisualizationError {
    /// File I/O error during visualization
    #[error("File I/O error: {message}")]
    FileError { message: String },

    /// Invalid output path
    #[error("Invalid output path: '{path}'. {reason}")]
    InvalidOutputPath { path: String, reason: String },

    /// Rendering backend error
    #[error("Rendering backend error: {message}")]
    RenderingError { message: String },

    /// Invalid visualization parameters
    #[error("Invalid visualization parameters: {parameter} = {value}. {constraint}")]
    InvalidParameters {
        parameter: String,
        value: String,
        constraint: String,
    },

    /// Empty channel system
    #[error("Cannot visualize empty channel system")]
    EmptyChannelSystem,

    /// Unsupported output format
    #[error("Unsupported output format: {format}. {message}")]
    UnsupportedFormat { format: String, message: String },

    /// Coordinate transformation error
    #[error("Coordinate transformation error: {message}")]
    CoordinateTransformError { message: String },
}

/// Errors related to channel type strategy operations
#[derive(Error, Debug)]
pub enum StrategyError {
    /// Strategy creation failed
    #[error("Failed to create strategy for channel type: {channel_type}. Reason: {reason}")]
    StrategyCreationFailed {
        channel_type: String,
        reason: String,
    },

    /// Invalid strategy parameters
    #[error("Invalid strategy parameters: {parameter} = {value}. {constraint}")]
    InvalidParameters {
        parameter: String,
        value: String,
        constraint: String,
    },

    /// Strategy execution failed
    #[error("Strategy execution failed for channel from ({from_x}, {from_y}) to ({to_x}, {to_y}): {reason}")]
    ExecutionFailed {
        from_x: f64,
        from_y: f64,
        to_x: f64,
        to_y: f64,
        reason: String,
    },

    /// Unsupported channel type
    #[error("Unsupported channel type: {channel_type}")]
    UnsupportedChannelType { channel_type: String },
}

/// Legacy simulation errors (kept for backward compatibility)
#[derive(Error, Debug)]
pub enum SimulationError {
    #[error("The simulation's linear system could not be solved. This may be due to disconnected channels or other geometry issues.")]
    /// Linear system could not be solved due to geometry issues
    LinearSystemError,
}

/// Convenient result type for scheme operations
pub type SchemeResult<T> = Result<T, SchemeError>;

/// Convenient result type for geometry operations
pub type GeometryResult<T> = Result<T, GeometryError>;

/// Convenient result type for configuration operations
pub type ConfigurationResult<T> = Result<T, ConfigurationError>;

/// Convenient result type for visualization operations
pub type VisualizationResult<T> = Result<T, VisualizationError>;

/// Convenient result type for strategy operations
pub type StrategyResult<T> = Result<T, StrategyError>;

impl GeometryError {
    /// Create an invalid point error
    #[must_use]
    pub const fn invalid_point(point: Point2D) -> Self {
        Self::InvalidPoint {
            x: point.0,
            y: point.1,
        }
    }

    /// Create an invalid box dimensions error
    #[must_use]
    pub const fn invalid_box_dimensions(width: f64, height: f64) -> Self {
        Self::InvalidBoxDimensions { width, height }
    }

    /// Create an insufficient space error
    #[must_use]
    pub const fn insufficient_space(required: f64, available: f64) -> Self {
        Self::InsufficientSpace {
            required,
            available,
        }
    }
}

impl ConfigurationError {
    /// Create an invalid geometry config error
    #[must_use]
    pub fn invalid_geometry_config(field: &str, value: f64, constraint: &str) -> Self {
        Self::InvalidGeometryConfig {
            field: field.to_string(),
            value,
            constraint: constraint.to_string(),
        }
    }

    /// Create an invalid serpentine config error
    #[must_use]
    pub fn invalid_serpentine_config(field: &str, value: f64, constraint: &str) -> Self {
        Self::InvalidSerpentineConfig {
            field: field.to_string(),
            value,
            constraint: constraint.to_string(),
        }
    }

    /// Create an invalid arc config error
    #[must_use]
    pub fn invalid_arc_config(field: &str, value: f64, constraint: &str) -> Self {
        Self::InvalidArcConfig {
            field: field.to_string(),
            value,
            constraint: constraint.to_string(),
        }
    }

    /// Create an invalid frustum config error
    #[must_use]
    pub fn invalid_frustum_config(field: &str, value: f64, constraint: &str) -> Self {
        Self::InvalidFrustumConfig {
            field: field.to_string(),
            value,
            constraint: constraint.to_string(),
        }
    }

    /// Create an invalid geometry generation config error
    #[must_use]
    pub fn invalid_generation_config(field: &str, constraint: &str) -> Self {
        Self::InvalidGenerationConfig {
            field: field.to_string(),
            constraint: constraint.to_string(),
        }
    }
}

impl VisualizationError {
    /// Create a file error
    #[must_use]
    pub fn file_error(message: &str) -> Self {
        Self::FileError {
            message: message.to_string(),
        }
    }

    /// Create an invalid output path error
    #[must_use]
    pub fn invalid_output_path(path: &str, reason: &str) -> Self {
        Self::InvalidOutputPath {
            path: path.to_string(),
            reason: reason.to_string(),
        }
    }

    /// Create a rendering error
    #[must_use]
    pub fn rendering_error(message: &str) -> Self {
        Self::RenderingError {
            message: message.to_string(),
        }
    }

    /// Create an unsupported format error
    #[must_use]
    pub fn unsupported_format(format: &str, message: &str) -> Self {
        Self::UnsupportedFormat {
            format: format.to_string(),
            message: message.to_string(),
        }
    }
}

impl StrategyError {
    /// Create a strategy creation failed error
    #[must_use]
    pub fn strategy_creation_failed(channel_type: &str, reason: &str) -> Self {
        Self::StrategyCreationFailed {
            channel_type: channel_type.to_string(),
            reason: reason.to_string(),
        }
    }

    /// Create an execution failed error
    #[must_use]
    pub fn execution_failed(from: Point2D, to: Point2D, reason: &str) -> Self {
        Self::ExecutionFailed {
            from_x: from.0,
            from_y: from.1,
            to_x: to.0,
            to_y: to.1,
            reason: reason.to_string(),
        }
    }
}
