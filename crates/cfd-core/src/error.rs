//! Error types and result aliases for the CFD simulation suite.
//!
//! This module provides a comprehensive error handling system with:
//! - Structured error types for specific failure modes
//! - Automatic conversion from external error types
//! - Context extension trait for adding error context

use std::fmt;
use thiserror::Error;

/// Core error type for CFD operations
#[derive(Debug, Error)]
pub enum Error {
    /// Invalid input parameters
    #[error("Invalid input: {0}")]
    InvalidInput(String),

    /// Invalid configuration
    #[error("Invalid configuration: {0}")]
    InvalidConfiguration(String),

    /// Numerical computation error
    #[error("Numerical error: {0}")]
    Numerical(NumericalErrorKind),

    /// GPU compute error
    #[error("GPU compute error: {0}")]
    GpuCompute(String),

    /// Convergence failure
    #[error("Convergence failed: {0}")]
    Convergence(ConvergenceErrorKind),

    /// Plugin-related errors
    #[error("Plugin error: {0}")]
    Plugin(PluginErrorKind),

    /// Solver-specific errors
    #[error("Solver error: {0}")]
    Solver(String),

    /// Numeric conversion error
    #[error("Conversion error: {0}")]
    ConversionError(String),

    /// I/O errors
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// Serialization/deserialization errors
    #[error("Serialization error: {0}")]
    Serialization(String),

    /// Unsupported operation
    #[error("Unsupported operation: {0}")]
    UnsupportedOperation(String),

    /// Physical invariant violation
    #[error("Physical invariant violation: {0}")]
    PhysicsViolation(String),

    /// Dimension mismatch
    #[error("Dimension mismatch: expected {expected}, got {actual}")]
    DimensionMismatch {
        /// Expected dimension
        expected: usize,
        /// Actual dimension
        actual: usize,
    },

    /// Index out of bounds
    #[error("Index out of bounds: {index} >= {size}")]
    IndexOutOfBounds {
        /// Index that was accessed
        index: usize,
        /// Size of the collection
        size: usize,
    },

    /// Boundary condition error
    #[error("Boundary error: {0}")]
    Boundary(BoundaryErrorKind),

    /// Geometry generation and validation error
    #[error("Geometry error: {0}")]
    Geometry(GeometryErrorKind),

    /// Configuration validation error
    #[error("Configuration error: {0}")]
    Configuration(ConfigurationErrorKind),

    /// Visualization and rendering error
    #[error("Visualization error: {0}")]
    Visualization(VisualizationErrorKind),

    /// Strategy execution error
    #[error("Strategy error: {0}")]
    Strategy(StrategyErrorKind),

    /// Parameter error
    #[error("Parameter error: {0}")]
    Parameter(ParameterErrorKind),

    /// Validation error
    #[error("Validation error: {0}")]
    Validation(ValidationErrorKind),

    /// Registry error
    #[error("Registry error: {0}")]
    Registry(RegistryErrorKind),

    /// Constraint error
    #[error("Constraint error: {0}")]
    Constraint(ConstraintErrorKind),

    /// Dependency error
    #[error("Dependency error: {0}")]
    Dependency(DependencyErrorKind),

    /// Adaptation error
    #[error("Adaptation error: {0}")]
    Adaptation(AdaptationErrorKind),

    /// Resistance calculation error
    #[error("Resistance calculation error: {0}")]
    ResistanceCalculation(ResistanceCalculationErrorKind),

    /// Not implemented
    #[error("Not implemented: {0}")]
    NotImplemented(String),

    /// Generic error with context
    #[error("{context}: {source}")]
    WithContext {
        /// Context description
        context: String,
        /// Underlying error
        #[source]
        source: Box<Error>,
    },
}

// ── Plugin error variants ────────────────────────────────────────────────────

/// Plugin error variants
#[derive(Debug, Clone)]
pub enum PluginErrorKind {
    /// Plugin not found
    NotFound {
        /// Name of the plugin
        name: String,
    },
    /// Plugin already registered
    AlreadyRegistered {
        /// Name of the plugin
        name: String,
    },
    /// Plugin initialization failed
    InitializationFailed {
        /// Name of the plugin
        name: String,
        /// Reason for failure
        reason: String,
    },
    /// Plugin execution failed
    ExecutionFailed {
        /// Name of the plugin
        name: String,
        /// Reason for failure
        reason: String,
    },
    /// Dependency not satisfied
    DependencyNotSatisfied {
        /// Plugin name
        plugin: String,
        /// Missing dependency
        dependency: String,
    },
    /// Circular dependency detected
    CircularDependency {
        /// Chain of dependencies forming the cycle
        chain: Vec<String>,
    },
    /// Invalid plugin configuration
    InvalidConfiguration {
        /// Name of the plugin
        name: String,
        /// Reason for invalid configuration
        reason: String,
    },
    /// Dependency missing
    DependencyMissing {
        /// Plugin name
        plugin: String,
        /// Missing dependency
        dependency: String,
    },
    /// Plugin has dependents
    HasDependents {
        /// Plugin name
        plugin: String,
        /// List of dependent plugins
        dependents: Vec<String>,
    },
}

impl fmt::Display for PluginErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotFound { name } => write!(f, "Plugin '{name}' not found"),
            Self::AlreadyRegistered { name } => write!(f, "Plugin '{name}' already registered"),
            Self::InitializationFailed { name, reason } => {
                write!(f, "Plugin '{name}' initialization failed: {reason}")
            }
            Self::ExecutionFailed { name, reason } => {
                write!(f, "Plugin '{name}' execution failed: {reason}")
            }
            Self::DependencyNotSatisfied { plugin, dependency } => write!(
                f,
                "Plugin '{plugin}' dependency '{dependency}' not satisfied"
            ),
            Self::CircularDependency { chain } => {
                write!(f, "Circular dependency detected: {}", chain.join(" -> "))
            }
            Self::InvalidConfiguration { name, reason } => {
                write!(f, "Invalid configuration for plugin '{name}': {reason}")
            }
            Self::DependencyMissing { plugin, dependency } => {
                write!(f, "Plugin '{plugin}' missing dependency '{dependency}'")
            }
            Self::HasDependents { plugin, dependents } => {
                write!(
                    f,
                    "Plugin '{}' has dependents: {}",
                    plugin,
                    dependents.join(", ")
                )
            }
        }
    }
}

// ── Numerical error variants ─────────────────────────────────────────────────

/// Numerical error variants
#[derive(Debug, Clone)]
pub enum NumericalErrorKind {
    /// Division by zero
    DivisionByZero,
    /// Invalid value (NaN or Inf)
    InvalidValue {
        /// String representation of the invalid value
        value: String,
    },
    /// Underflow occurred
    Underflow {
        /// Value that caused underflow
        value: f64,
    },
    /// Overflow occurred
    Overflow {
        /// Value that caused overflow
        value: f64,
    },
    /// Matrix is singular
    SingularMatrix,
    /// Matrix is not positive definite
    NotPositiveDefinite,
    /// Invalid tolerance
    InvalidTolerance {
        /// The invalid tolerance value
        tolerance: f64,
    },
    /// Insufficient precision
    InsufficientPrecision {
        /// Precision achieved
        achieved: f64,
        /// Precision required
        required: f64,
    },
}

impl fmt::Display for NumericalErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DivisionByZero => write!(f, "Division by zero"),
            Self::InvalidValue { value } => write!(f, "Invalid numerical value: {value}"),
            Self::Underflow { value } => write!(f, "Numerical underflow: {value:.2e}"),
            Self::Overflow { value } => write!(f, "Numerical overflow: {value:.2e}"),
            Self::SingularMatrix => write!(f, "Matrix is singular"),
            Self::NotPositiveDefinite => write!(f, "Matrix is not positive definite"),
            Self::InvalidTolerance { tolerance } => {
                write!(f, "Invalid tolerance: {tolerance:.2e}")
            }
            Self::InsufficientPrecision { achieved, required } => write!(
                f,
                "Insufficient precision: achieved {achieved:.2e}, required {required:.2e}"
            ),
        }
    }
}

// ── Convergence error variants ───────────────────────────────────────────────

/// Convergence error variants
#[derive(Debug, Clone)]
pub enum ConvergenceErrorKind {
    /// Maximum iterations exceeded
    MaxIterationsExceeded {
        /// Maximum iteration limit
        max: usize,
    },
    /// Residual did not decrease
    StagnatedResidual {
        /// Current residual value
        residual: f64,
    },
    /// Solution diverged
    Diverged {
        /// Norm of the diverging solution
        norm: f64,
    },
    /// NaN or Inf detected
    InvalidValue,
    /// Algorithm breakdown (e.g., zero inner product in Krylov methods)
    Breakdown,
}

impl fmt::Display for ConvergenceErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MaxIterationsExceeded { max } => {
                write!(f, "Maximum iterations ({max}) exceeded")
            }
            Self::StagnatedResidual { residual } => {
                write!(f, "Residual stagnated at {residual:.2e}")
            }
            Self::Diverged { norm } => write!(f, "Solution diverged with norm {norm:.2e}"),
            Self::InvalidValue => write!(f, "Invalid value (NaN or Inf) detected"),
            Self::Breakdown => write!(f, "Algorithm breakdown (zero inner product)"),
        }
    }
}

// ── Boundary error variants ─────────────────────────────────────────────────

/// Boundary condition error variants
#[derive(Debug, Clone)]
pub enum BoundaryErrorKind {
    /// Insufficient interior values for requested stencil order
    InsufficientStencil {
        /// Required number of interior values
        required: usize,
        /// Stencil order
        order: usize,
        /// Actual number provided
        actual: usize,
    },

    /// Unsupported stencil order
    UnsupportedOrder(
        /// Stencil order number
        usize,
    ),

    /// Robin condition singularity
    RobinSingularity {
        /// Denominator value
        value: f64,
    },

    /// Invalid boundary region
    InvalidRegion(
        /// Description of the invalid region
        String,
    ),

    /// Dimension mismatch (string-based, for boundary metadata)
    DimensionMismatch {
        /// Expected dimensions
        expected: String,
        /// Actual dimensions
        actual: String,
    },
}

impl fmt::Display for BoundaryErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InsufficientStencil {
                required,
                order,
                actual,
            } => write!(
                f,
                "Insufficient interior values: need {required} for order {order} scheme, got {actual}"
            ),
            Self::UnsupportedOrder(order) => {
                write!(f, "Stencil order {order} is unsupported (supported: 1-4)")
            }
            Self::RobinSingularity { value } => {
                write!(f, "Robin condition denominator near zero: β + α*dx = {value}")
            }
            Self::InvalidRegion(msg) => write!(f, "Invalid boundary region: {msg}"),
            Self::DimensionMismatch { expected, actual } => {
                write!(f, "Dimension mismatch: expected {expected}, got {actual}")
            }
        }
    }
}

// ── Geometry error variants ─────────────────────────────────────────────────

/// Geometry generation and validation error variants
#[derive(Debug, Clone)]
pub enum GeometryErrorKind {
    /// Invalid box dimensions
    InvalidBoxDimensions {
        /// Width of the box
        width: f64,
        /// Height of the box
        height: f64,
    },
    /// Invalid point coordinates
    InvalidPoint {
        /// X coordinate
        x: f64,
        /// Y coordinate
        y: f64,
    },
    /// Insufficient space for generation
    InsufficientSpace {
        /// Required space
        required: f64,
        /// Available space
        available: f64,
    },
    /// Invalid split pattern
    InvalidSplitPattern {
        /// Reason for failure
        reason: String,
    },
    /// Node creation failed
    NodeCreationFailed {
        /// X position
        x: f64,
        /// Y position
        y: f64,
        /// Reason for failure
        reason: String,
    },
    /// Channel creation failed
    ChannelCreationFailed {
        /// Source node ID
        from_id: usize,
        /// Destination node ID
        to_id: usize,
        /// Reason for failure
        reason: String,
    },
    /// Invalid channel path
    InvalidChannelPath {
        /// Reason for failure
        reason: String,
    },
    /// Overlapping channels detected
    OverlappingChannels {
        /// First point X
        x1: f64,
        /// First point Y
        y1: f64,
        /// Second point X
        x2: f64,
        /// Second point Y
        y2: f64,
    },
}

impl GeometryErrorKind {
    /// Create an invalid box dimensions error
    #[must_use]
    pub const fn invalid_box_dimensions(width: f64, height: f64) -> Self {
        Self::InvalidBoxDimensions { width, height }
    }

    /// Create an invalid point error
    #[must_use]
    pub const fn invalid_point(x: f64, y: f64) -> Self {
        Self::InvalidPoint { x, y }
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

impl fmt::Display for GeometryErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidBoxDimensions { width, height } => write!(
                f,
                "Invalid box dimensions: width={width}, height={height}. Both must be positive."
            ),
            Self::InvalidPoint { x, y } => {
                write!(f, "Invalid point coordinates: ({x}, {y}). Must be finite.")
            }
            Self::InsufficientSpace {
                required,
                available,
            } => write!(
                f,
                "Insufficient space: required {required}, available {available}"
            ),
            Self::InvalidSplitPattern { reason } => write!(f, "Invalid split pattern: {reason}"),
            Self::NodeCreationFailed { x, y, reason } => {
                write!(f, "Failed to create node at ({x}, {y}): {reason}")
            }
            Self::ChannelCreationFailed {
                from_id,
                to_id,
                reason,
            } => write!(
                f,
                "Failed to create channel from {from_id} to {to_id}: {reason}"
            ),
            Self::InvalidChannelPath { reason } => write!(f, "Invalid channel path: {reason}"),
            Self::OverlappingChannels { x1, y1, x2, y2 } => write!(
                f,
                "Overlapping channels between ({x1},{y1}) and ({x2},{y2})"
            ),
        }
    }
}

// ── Configuration error variants ────────────────────────────────────────────

/// Configuration validation error variants
#[derive(Debug, Clone)]
pub enum ConfigurationErrorKind {
    /// Invalid geometry configuration
    InvalidGeometryConfig {
        /// Configuration field name
        field: String,
        /// Invalid value
        value: f64,
        /// Constraint description
        constraint: String,
    },
    /// Invalid serpentine configuration
    InvalidSerpentineConfig {
        /// Configuration field name
        field: String,
        /// Invalid value
        value: f64,
        /// Constraint description
        constraint: String,
    },
    /// Invalid arc configuration
    InvalidArcConfig {
        /// Configuration field name
        field: String,
        /// Invalid value
        value: f64,
        /// Constraint description
        constraint: String,
    },
    /// Invalid frustum configuration
    InvalidFrustumConfig {
        /// Configuration field name
        field: String,
        /// Invalid value
        value: f64,
        /// Constraint description
        constraint: String,
    },
    /// Invalid geometry generation configuration
    InvalidGenerationConfig {
        /// Configuration field name
        field: String,
        /// Constraint description
        constraint: String,
    },
    /// Conflicting configuration values
    ConflictingValues {
        /// Description of the conflict
        conflict: String,
    },
    /// Missing required configuration
    MissingConfiguration {
        /// Name of the missing field
        field: String,
    },
}

impl ConfigurationErrorKind {
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

    /// Create an invalid generation config error
    #[must_use]
    pub fn invalid_generation_config(field: &str, constraint: &str) -> Self {
        Self::InvalidGenerationConfig {
            field: field.to_string(),
            constraint: constraint.to_string(),
        }
    }
}

impl fmt::Display for ConfigurationErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidGeometryConfig {
                field,
                value,
                constraint,
            } => write!(
                f,
                "Invalid geometry config: {field} = {value}. {constraint}"
            ),
            Self::InvalidSerpentineConfig {
                field,
                value,
                constraint,
            } => write!(
                f,
                "Invalid serpentine config: {field} = {value}. {constraint}"
            ),
            Self::InvalidArcConfig {
                field,
                value,
                constraint,
            } => write!(f, "Invalid arc config: {field} = {value}. {constraint}"),
            Self::InvalidFrustumConfig {
                field,
                value,
                constraint,
            } => write!(f, "Invalid frustum config: {field} = {value}. {constraint}"),
            Self::InvalidGenerationConfig { field, constraint } => {
                write!(f, "Invalid generation config: {field}. {constraint}")
            }
            Self::ConflictingValues { conflict } => {
                write!(f, "Conflicting configuration values: {conflict}")
            }
            Self::MissingConfiguration { field } => {
                write!(f, "Missing required configuration: {field}")
            }
        }
    }
}

// ── Visualization error variants ────────────────────────────────────────────

/// Visualization and rendering error variants
#[derive(Debug, Clone)]
pub enum VisualizationErrorKind {
    /// File I/O error during visualization
    FileError {
        /// Error message
        message: String,
    },
    /// Invalid output path
    InvalidOutputPath {
        /// Path that was invalid
        path: String,
        /// Reason for failure
        reason: String,
    },
    /// Rendering backend error
    RenderingError {
        /// Error message
        message: String,
    },
    /// Invalid visualization parameters
    InvalidParameters {
        /// Parameter name
        parameter: String,
        /// Parameter value
        value: String,
        /// Constraint description
        constraint: String,
    },
    /// Empty channel system
    EmptyChannelSystem,
    /// Unsupported output format
    UnsupportedFormat {
        /// Format name
        format: String,
        /// Additional message
        message: String,
    },
    /// Coordinate transformation error
    CoordinateTransformError {
        /// Error message
        message: String,
    },
}

impl VisualizationErrorKind {
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

impl fmt::Display for VisualizationErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::FileError { message } => write!(f, "File I/O error: {message}"),
            Self::InvalidOutputPath { path, reason } => {
                write!(f, "Invalid output path: '{path}'. {reason}")
            }
            Self::RenderingError { message } => write!(f, "Rendering error: {message}"),
            Self::InvalidParameters {
                parameter,
                value,
                constraint,
            } => write!(f, "Invalid parameters: {parameter} = {value}. {constraint}"),
            Self::EmptyChannelSystem => write!(f, "Cannot visualize empty channel system"),
            Self::UnsupportedFormat { format, message } => {
                write!(f, "Unsupported format: {format}. {message}")
            }
            Self::CoordinateTransformError { message } => {
                write!(f, "Coordinate transform error: {message}")
            }
        }
    }
}

// ── Strategy error variants ─────────────────────────────────────────────────

/// Strategy execution error variants
#[derive(Debug, Clone)]
pub enum StrategyErrorKind {
    /// Strategy creation failed
    StrategyCreationFailed {
        /// Channel type
        channel_type: String,
        /// Reason for failure
        reason: String,
    },
    /// Invalid strategy parameters
    InvalidParameters {
        /// Parameter name
        parameter: String,
        /// Parameter value
        value: String,
        /// Constraint description
        constraint: String,
    },
    /// Strategy execution failed
    ExecutionFailed {
        /// Source X coordinate
        from_x: f64,
        /// Source Y coordinate
        from_y: f64,
        /// Destination X coordinate
        to_x: f64,
        /// Destination Y coordinate
        to_y: f64,
        /// Reason for failure
        reason: String,
    },
    /// Unsupported channel type
    UnsupportedChannelType {
        /// Channel type name
        channel_type: String,
    },
}

impl StrategyErrorKind {
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
    pub fn execution_failed(from_x: f64, from_y: f64, to_x: f64, to_y: f64, reason: &str) -> Self {
        Self::ExecutionFailed {
            from_x,
            from_y,
            to_x,
            to_y,
            reason: reason.to_string(),
        }
    }
}

impl fmt::Display for StrategyErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::StrategyCreationFailed {
                channel_type,
                reason,
            } => write!(f, "Failed to create strategy for {channel_type}: {reason}"),
            Self::InvalidParameters {
                parameter,
                value,
                constraint,
            } => write!(
                f,
                "Invalid strategy parameters: {parameter} = {value}. {constraint}"
            ),
            Self::ExecutionFailed {
                from_x,
                from_y,
                to_x,
                to_y,
                reason,
            } => write!(
                f,
                "Strategy execution failed from ({from_x},{from_y}) to ({to_x},{to_y}): {reason}"
            ),
            Self::UnsupportedChannelType { channel_type } => {
                write!(f, "Unsupported channel type: {channel_type}")
            }
        }
    }
}

// ── Parameter error variants ────────────────────────────────────────────────

/// Parameter error variants
#[derive(Debug, Clone)]
pub enum ParameterErrorKind {
    /// Parameter not found
    NotFound {
        /// Parameter name
        name: String,
        /// Domain name
        domain: String,
    },
    /// Invalid parameter value
    InvalidValue {
        /// Parameter name
        name: String,
        /// Invalid value
        value: String,
        /// Constraint description
        constraint: String,
    },
    /// Parameter type mismatch
    TypeMismatch {
        /// Parameter name
        name: String,
        /// Expected type
        expected: String,
        /// Actual type
        actual: String,
    },
    /// Parameter is read-only
    ReadOnly {
        /// Parameter name
        name: String,
    },
    /// Parameter dependency not satisfied
    DependencyNotSatisfied {
        /// Parameter name
        name: String,
        /// Missing dependency
        dependency: String,
    },
    /// Circular dependency detected
    CircularDependency {
        /// Parameter name
        name: String,
    },
    /// Adaptation failed
    AdaptationFailed {
        /// Parameter name
        name: String,
        /// Reason for failure
        reason: String,
    },
}

impl ParameterErrorKind {
    /// Create a not-found error
    #[must_use]
    pub fn not_found(name: &str, domain: &str) -> Self {
        Self::NotFound {
            name: name.to_string(),
            domain: domain.to_string(),
        }
    }

    /// Create an invalid-value error
    pub fn invalid_value(name: &str, value: &dyn fmt::Debug, constraint: &str) -> Self {
        Self::InvalidValue {
            name: name.to_string(),
            value: format!("{value:?}"),
            constraint: constraint.to_string(),
        }
    }

    /// Create a type-mismatch error
    #[must_use]
    pub fn type_mismatch(name: &str, expected: &str, actual: &str) -> Self {
        Self::TypeMismatch {
            name: name.to_string(),
            expected: expected.to_string(),
            actual: actual.to_string(),
        }
    }

    /// Create a read-only error
    #[must_use]
    pub fn read_only(name: &str) -> Self {
        Self::ReadOnly {
            name: name.to_string(),
        }
    }

    /// Create a dependency-not-satisfied error
    #[must_use]
    pub fn dependency_not_satisfied(name: &str, dependency: &str) -> Self {
        Self::DependencyNotSatisfied {
            name: name.to_string(),
            dependency: dependency.to_string(),
        }
    }

    /// Create a circular-dependency error
    #[must_use]
    pub fn circular_dependency(name: &str) -> Self {
        Self::CircularDependency {
            name: name.to_string(),
        }
    }

    /// Create an adaptation-failed error
    #[must_use]
    pub fn adaptation_failed(name: &str, reason: &str) -> Self {
        Self::AdaptationFailed {
            name: name.to_string(),
            reason: reason.to_string(),
        }
    }
}

impl fmt::Display for ParameterErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotFound { name, domain } => {
                write!(f, "Parameter '{name}' not found in domain '{domain}'")
            }
            Self::InvalidValue {
                name,
                value,
                constraint,
            } => write!(f, "Invalid value for '{name}': {value}. {constraint}"),
            Self::TypeMismatch {
                name,
                expected,
                actual,
            } => write!(
                f,
                "Type mismatch for '{name}': expected {expected}, got {actual}"
            ),
            Self::ReadOnly { name } => write!(f, "Parameter '{name}' is read-only"),
            Self::DependencyNotSatisfied { name, dependency } => write!(
                f,
                "Parameter '{name}' dependency '{dependency}' not satisfied"
            ),
            Self::CircularDependency { name } => {
                write!(f, "Circular dependency involving parameter '{name}'")
            }
            Self::AdaptationFailed { name, reason } => {
                write!(f, "Failed to adapt parameter '{name}': {reason}")
            }
        }
    }
}

// ── Validation error variants ───────────────────────────────────────────────

/// Validation error variants
#[derive(Debug, Clone)]
pub enum ValidationErrorKind {
    /// Validation rule failed
    RuleFailed {
        /// Field that failed validation
        field: String,
        /// Error message
        message: String,
    },
    /// Multiple validation failures
    Multiple {
        /// List of failures
        failures: Vec<ValidationErrorKind>,
    },
    /// Custom validation error
    Custom {
        /// Error message
        message: String,
    },
    /// Constraint violation
    ConstraintViolation {
        /// Constraint description
        constraint: String,
    },
}

impl ValidationErrorKind {
    /// Create a rule-failed error
    #[must_use]
    pub fn rule_failed(field: &str, message: &str) -> Self {
        Self::RuleFailed {
            field: field.to_string(),
            message: message.to_string(),
        }
    }

    /// Create a multiple-failures error
    #[must_use]
    pub const fn multiple(failures: Vec<Self>) -> Self {
        Self::Multiple { failures }
    }

    /// Create a custom error
    #[must_use]
    pub fn custom(message: &str) -> Self {
        Self::Custom {
            message: message.to_string(),
        }
    }

    /// Create a constraint-violation error
    #[must_use]
    pub fn constraint_violation(constraint: &str) -> Self {
        Self::ConstraintViolation {
            constraint: constraint.to_string(),
        }
    }
}

impl fmt::Display for ValidationErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::RuleFailed { field, message } => {
                write!(f, "Validation failed for '{field}': {message}")
            }
            Self::Multiple { failures } => write!(
                f,
                "Multiple validation failures ({} errors)",
                failures.len()
            ),
            Self::Custom { message } => write!(f, "Validation error: {message}"),
            Self::ConstraintViolation { constraint } => {
                write!(f, "Constraint violation: {constraint}")
            }
        }
    }
}

// ── Registry error variants ─────────────────────────────────────────────────

/// Registry error variants
#[derive(Debug, Clone)]
pub enum RegistryErrorKind {
    /// Manager not found
    ManagerNotFound {
        /// Domain name
        domain: String,
    },
    /// Manager already registered
    ManagerAlreadyRegistered {
        /// Domain name
        domain: String,
    },
    /// Registry is locked
    RegistryLocked,
    /// Initialization failed
    InitializationFailed {
        /// Reason for failure
        reason: String,
    },
    /// Update conflict
    UpdateConflict {
        /// Parameter name
        parameter: String,
        /// Reason for conflict
        reason: String,
    },
}

impl fmt::Display for RegistryErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ManagerNotFound { domain } => {
                write!(f, "Parameter manager for domain '{domain}' not found")
            }
            Self::ManagerAlreadyRegistered { domain } => write!(
                f,
                "Parameter manager for domain '{domain}' already registered"
            ),
            Self::RegistryLocked => write!(f, "Registry is locked"),
            Self::InitializationFailed { reason } => {
                write!(f, "Registry initialization failed: {reason}")
            }
            Self::UpdateConflict { parameter, reason } => {
                write!(f, "Update conflict for '{parameter}': {reason}")
            }
        }
    }
}

// ── Constraint error variants ───────────────────────────────────────────────

/// Constraint error variants
#[derive(Debug, Clone)]
pub enum ConstraintErrorKind {
    /// Range constraint violation
    RangeViolation {
        /// Violating value
        value: String,
        /// Minimum allowed
        min: String,
        /// Maximum allowed
        max: String,
    },
    /// Set constraint violation
    SetViolation {
        /// Violating value
        value: String,
        /// Allowed values
        allowed: Vec<String>,
    },
    /// Custom constraint violation
    CustomViolation {
        /// Error message
        message: String,
    },
    /// Constraint composition error
    CompositionError {
        /// Reason for failure
        reason: String,
    },
}

impl ConstraintErrorKind {
    /// Create a range-violation error
    pub fn range_violation(
        value: &dyn fmt::Debug,
        min: &dyn fmt::Debug,
        max: &dyn fmt::Debug,
    ) -> Self {
        Self::RangeViolation {
            value: format!("{value:?}"),
            min: format!("{min:?}"),
            max: format!("{max:?}"),
        }
    }

    /// Create a set-violation error
    pub fn set_violation(value: &dyn fmt::Debug, allowed: &[&dyn fmt::Debug]) -> Self {
        Self::SetViolation {
            value: format!("{value:?}"),
            allowed: allowed.iter().map(|v| format!("{v:?}")).collect(),
        }
    }

    /// Create a custom-violation error
    #[must_use]
    pub fn custom_violation(message: &str) -> Self {
        Self::CustomViolation {
            message: message.to_string(),
        }
    }
}

impl fmt::Display for ConstraintErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::RangeViolation { value, min, max } => {
                write!(f, "Value {value} is outside range [{min}, {max}]")
            }
            Self::SetViolation { value, allowed } => {
                write!(f, "Value {value} is not in allowed set: {allowed:?}")
            }
            Self::CustomViolation { message } => write!(f, "Constraint violation: {message}"),
            Self::CompositionError { reason } => {
                write!(f, "Constraint composition error: {reason}")
            }
        }
    }
}

// ── Dependency error variants ───────────────────────────────────────────────

/// Dependency error variants
#[derive(Debug, Clone)]
pub enum DependencyErrorKind {
    /// Dependency cycle detected
    CycleDetected {
        /// Cycle path
        cycle: Vec<String>,
    },
    /// Missing dependency
    MissingDependency {
        /// Parameter that has the dependency
        parameter: String,
        /// Missing dependency name
        dependency: String,
    },
    /// Dependency resolution failed
    ResolutionFailed {
        /// Parameter name
        parameter: String,
        /// Reason for failure
        reason: String,
    },
    /// Invalid dependency graph
    InvalidGraph {
        /// Reason for failure
        reason: String,
    },
}

impl fmt::Display for DependencyErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::CycleDetected { cycle } => write!(f, "Dependency cycle: {}", cycle.join(" -> ")),
            Self::MissingDependency {
                parameter,
                dependency,
            } => write!(f, "Missing dependency '{dependency}' for '{parameter}'"),
            Self::ResolutionFailed { parameter, reason } => write!(
                f,
                "Dependency resolution failed for '{parameter}': {reason}"
            ),
            Self::InvalidGraph { reason } => write!(f, "Invalid dependency graph: {reason}"),
        }
    }
}

// ── Adaptation error variants ───────────────────────────────────────────────

/// Adaptation error variants
#[derive(Debug, Clone)]
pub enum AdaptationErrorKind {
    /// Invalid context provided
    InvalidContext {
        /// Reason for failure
        reason: String,
    },
    /// Adaptation calculation failed
    CalculationFailed {
        /// Parameter name
        parameter: String,
        /// Reason for failure
        reason: String,
    },
    /// Result value is invalid
    InvalidResult {
        /// Invalid value
        value: String,
        /// Constraint description
        constraint: String,
    },
    /// Dependency not available
    DependencyMissing {
        /// Missing dependency name
        dependency: String,
    },
}

impl fmt::Display for AdaptationErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidContext { reason } => write!(f, "Invalid adaptation context: {reason}"),
            Self::CalculationFailed { parameter, reason } => {
                write!(f, "Adaptation failed for '{parameter}': {reason}")
            }
            Self::InvalidResult { value, constraint } => {
                write!(f, "Adaptation result '{value}' violates: {constraint}")
            }
            Self::DependencyMissing { dependency } => {
                write!(f, "Required dependency '{dependency}' not available")
            }
        }
    }
}

// ── Result type and helpers ─────────────────────────────────────────────────

/// Result type alias for CFD operations
pub type Result<T> = std::result::Result<T, Error>;

/// Extension trait for adding context to errors
pub trait ErrorContext<T> {
    /// Add context to an error
    ///
    /// # Errors
    /// Returns the original error with added context
    fn context(self, msg: impl Into<String>) -> Result<T>;

    /// Add context with a closure (lazy evaluation)
    ///
    /// # Errors
    /// Returns the original error with context added via the closure
    fn with_context<F>(self, f: F) -> Result<T>
    where
        F: FnOnce() -> String;
}

impl<T> ErrorContext<T> for Result<T> {
    fn context(self, msg: impl Into<String>) -> Result<T> {
        self.map_err(|e| Error::WithContext {
            context: msg.into(),
            source: Box::new(e),
        })
    }

    fn with_context<F>(self, f: F) -> Result<T>
    where
        F: FnOnce() -> String,
    {
        self.map_err(|e| Error::WithContext {
            context: f(),
            source: Box::new(e),
        })
    }
}

/// Helper function to convert Option to Result
///
/// # Errors
/// Returns `Error::InvalidInput` if the option is None
pub fn require<T>(opt: Option<T>, msg: impl Into<String>) -> Result<T> {
    opt.ok_or_else(|| Error::InvalidInput(msg.into()))
}

impl From<&str> for Error {
    fn from(msg: &str) -> Self {
        Error::InvalidInput(msg.to_string())
    }
}

// ── std::error::Error impls for all Kind types ─────────────────────────────

impl std::error::Error for GeometryErrorKind {}
impl std::error::Error for ConfigurationErrorKind {}
impl std::error::Error for VisualizationErrorKind {}
impl std::error::Error for StrategyErrorKind {}
impl std::error::Error for ParameterErrorKind {}
impl std::error::Error for ValidationErrorKind {}
impl std::error::Error for RegistryErrorKind {}
impl std::error::Error for ConstraintErrorKind {}
impl std::error::Error for DependencyErrorKind {}
impl std::error::Error for AdaptationErrorKind {}
impl std::error::Error for BoundaryErrorKind {}
impl std::error::Error for NumericalErrorKind {}
impl std::error::Error for ConvergenceErrorKind {}
impl std::error::Error for PluginErrorKind {}

// ── Resistance calculation error variants ───────────────────────────────────

/// Resistance calculation error variants
#[derive(Debug, Clone)]
pub enum ResistanceCalculationErrorKind {
    /// Hydraulic diameter is missing for a component
    MissingHydraulicDiameter,
    /// Resistance model calculation failed
    ModelError(
        /// Error message
        String,
    ),
    /// Invalid flow conditions
    InvalidFlowConditions(
        /// Error message
        String,
    ),
    /// Numerical computation error
    NumericalError(
        /// Error message
        String,
    ),
}

impl fmt::Display for ResistanceCalculationErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingHydraulicDiameter => {
                write!(f, "Hydraulic diameter is missing for component")
            }
            Self::ModelError(msg) => write!(f, "Resistance model calculation failed: {msg}"),
            Self::InvalidFlowConditions(msg) => write!(f, "Invalid flow conditions: {msg}"),
            Self::NumericalError(msg) => write!(f, "Numerical computation error: {msg}"),
        }
    }
}

impl std::error::Error for ResistanceCalculationErrorKind {}

impl From<String> for Error {
    fn from(msg: String) -> Self {
        Error::InvalidInput(msg)
    }
}

// ── From impls for all Kind types → Error ───────────────────────────────────

impl From<GeometryErrorKind> for Error {
    fn from(kind: GeometryErrorKind) -> Self {
        Error::Geometry(kind)
    }
}
impl From<ConfigurationErrorKind> for Error {
    fn from(kind: ConfigurationErrorKind) -> Self {
        Error::Configuration(kind)
    }
}
impl From<VisualizationErrorKind> for Error {
    fn from(kind: VisualizationErrorKind) -> Self {
        Error::Visualization(kind)
    }
}
impl From<StrategyErrorKind> for Error {
    fn from(kind: StrategyErrorKind) -> Self {
        Error::Strategy(kind)
    }
}
impl From<ParameterErrorKind> for Error {
    fn from(kind: ParameterErrorKind) -> Self {
        Error::Parameter(kind)
    }
}
impl From<ValidationErrorKind> for Error {
    fn from(kind: ValidationErrorKind) -> Self {
        Error::Validation(kind)
    }
}
impl From<RegistryErrorKind> for Error {
    fn from(kind: RegistryErrorKind) -> Self {
        Error::Registry(kind)
    }
}
impl From<ConstraintErrorKind> for Error {
    fn from(kind: ConstraintErrorKind) -> Self {
        Error::Constraint(kind)
    }
}
impl From<DependencyErrorKind> for Error {
    fn from(kind: DependencyErrorKind) -> Self {
        Error::Dependency(kind)
    }
}
impl From<AdaptationErrorKind> for Error {
    fn from(kind: AdaptationErrorKind) -> Self {
        Error::Adaptation(kind)
    }
}
impl From<ResistanceCalculationErrorKind> for Error {
    fn from(kind: ResistanceCalculationErrorKind) -> Self {
        Error::ResistanceCalculation(kind)
    }
}
impl From<BoundaryErrorKind> for Error {
    fn from(kind: BoundaryErrorKind) -> Self {
        Error::Boundary(kind)
    }
}

// Cross-kind conversions used by state_management
impl From<ConstraintErrorKind> for ParameterErrorKind {
    fn from(err: ConstraintErrorKind) -> Self {
        Self::InvalidValue {
            name: String::new(),
            value: String::new(),
            constraint: err.to_string(),
        }
    }
}
impl From<ValidationErrorKind> for ParameterErrorKind {
    fn from(err: ValidationErrorKind) -> Self {
        Self::InvalidValue {
            name: String::new(),
            value: String::new(),
            constraint: err.to_string(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_context() {
        let result: Result<()> = Err(Error::InvalidInput("test".into()));
        let with_context = result.context("Additional context");
        assert!(with_context.is_err());
        let error_msg = format!("{}", with_context.unwrap_err());
        assert!(error_msg.contains("Additional context"));
    }

    #[test]
    fn test_require() {
        let some_value = Some(42);
        assert_eq!(require(some_value, "missing").unwrap(), 42);

        let none_value: Option<i32> = None;
        assert!(require(none_value, "missing").is_err());
    }
}
