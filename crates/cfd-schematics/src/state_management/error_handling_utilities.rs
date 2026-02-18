//! Error handling utilities for the state management system
//!
//! This module provides comprehensive error handling utilities that improve
//! error context, validation, and recovery mechanisms throughout the state
//! management system, following SOLID, CUPID, and GRASP design principles.

use crate::{
    state_management::{
        adaptive::{AdaptationError, ChannelGenerationContext},
        errors::{ConstraintError, ParameterError, StateManagementError, ValidationError},
        parameters::ConfigurableParameter,
    },
    // Removed unused error imports
};
use std::fmt::Debug;

/// Enhanced error context for parameter operations
#[derive(Debug, Clone)]
pub struct ParameterErrorContext {
    /// Parameter name
    pub parameter_name: String,

    /// Parameter domain/category
    pub domain: String,

    /// Operation being performed
    pub operation: String,

    /// Current parameter value (if available)
    pub current_value: Option<String>,

    /// Expected value or constraint
    pub expected: Option<String>,

    /// Additional context information
    pub context: Vec<(String, String)>,
}

impl ParameterErrorContext {
    /// Create a new parameter error context
    #[must_use]
    pub fn new(parameter_name: &str, domain: &str, operation: &str) -> Self {
        Self {
            parameter_name: parameter_name.to_string(),
            domain: domain.to_string(),
            operation: operation.to_string(),
            current_value: None,
            expected: None,
            context: Vec::new(),
        }
    }

    /// Add current value information
    #[must_use]
    pub fn with_current_value(mut self, value: &str) -> Self {
        self.current_value = Some(value.to_string());
        self
    }

    /// Add expected value information
    #[must_use]
    pub fn with_expected(mut self, expected: &str) -> Self {
        self.expected = Some(expected.to_string());
        self
    }

    /// Add context information
    #[must_use]
    pub fn with_context(mut self, key: &str, value: &str) -> Self {
        self.context.push((key.to_string(), value.to_string()));
        self
    }

    /// Generate a detailed error message
    #[must_use]
    pub fn to_detailed_message(&self) -> String {
        let mut message = format!(
            "Parameter '{}' in domain '{}' failed during operation '{}'",
            self.parameter_name, self.domain, self.operation
        );

        if let Some(current) = &self.current_value {
            message.push_str(&format!(" (current value: {current})"));
        }

        if let Some(expected) = &self.expected {
            message.push_str(&format!(" (expected: {expected})"));
        }

        if !self.context.is_empty() {
            message.push_str(" [Context: ");
            let context_parts: Vec<String> = self
                .context
                .iter()
                .map(|(k, v)| format!("{k}={v}"))
                .collect();
            message.push_str(&context_parts.join(", "));
            message.push(']');
        }

        message
    }
}

/// Enhanced validation result with detailed error information
#[derive(Debug)]
pub struct ValidationResult {
    /// Whether validation passed
    pub is_valid: bool,

    /// Validation errors (if any)
    pub errors: Vec<ValidationError>,

    /// Validation warnings (non-fatal issues)
    pub warnings: Vec<String>,

    /// Validation context
    pub context: ParameterErrorContext,
}

impl ValidationResult {
    /// Create a successful validation result
    #[must_use]
    pub const fn success(context: ParameterErrorContext) -> Self {
        Self {
            is_valid: true,
            errors: Vec::new(),
            warnings: Vec::new(),
            context,
        }
    }

    /// Create a failed validation result
    #[must_use]
    pub const fn failure(context: ParameterErrorContext, errors: Vec<ValidationError>) -> Self {
        Self {
            is_valid: false,
            errors,
            warnings: Vec::new(),
            context,
        }
    }

    /// Add a warning to the validation result
    #[must_use]
    pub fn with_warning(mut self, warning: String) -> Self {
        self.warnings.push(warning);
        self
    }

    /// Convert to a Result type
    pub fn into_result(self) -> Result<Vec<String>, StateManagementError> {
        if self.is_valid {
            Ok(self.warnings)
        } else {
            // Create a comprehensive error message
            let error_messages: Vec<String> = self
                .errors
                .iter()
                .map(std::string::ToString::to_string)
                .collect();

            let detailed_message = format!(
                "{}\nValidation errors: {}",
                self.context.to_detailed_message(),
                error_messages.join("; ")
            );

            Err(StateManagementError::Validation(
                ValidationError::rule_failed(&self.context.parameter_name, &detailed_message),
            ))
        }
    }
}

/// Enhanced parameter validation utilities
pub struct ParameterValidator;

impl ParameterValidator {
    /// Validate a parameter with enhanced error context
    pub fn validate_with_context<T>(
        parameter: &ConfigurableParameter<T>,
        context: ParameterErrorContext,
    ) -> ValidationResult
    where
        T: Clone + Debug + PartialEq + PartialOrd + 'static,
    {
        let mut errors = Vec::new();
        let mut warnings = Vec::new();

        // Validate current value
        if let Err(constraint_error) = parameter.validate() {
            errors.push(ValidationError::constraint_violation(
                &constraint_error.to_string(),
            ));
        }

        // Check for potential issues
        if parameter.is_modified() && parameter.get_change_history().len() > 50 {
            warnings.push(format!(
                "Parameter '{}' has been modified {} times, consider reviewing change frequency",
                context.parameter_name,
                parameter.get_change_history().len()
            ));
        }

        let mut result = if errors.is_empty() {
            ValidationResult::success(context)
        } else {
            ValidationResult::failure(context, errors)
        };

        for warning in warnings {
            result = result.with_warning(warning);
        }

        result
    }

    /// Validate adaptation context
    pub fn validate_adaptation_context(
        context: &ChannelGenerationContext,
        _parameter_name: &str,
    ) -> Result<(), AdaptationError> {
        // Validate basic context properties
        if context.total_branches == 0 {
            return Err(AdaptationError::InvalidContext {
                reason: "Total branches must be greater than zero".to_string(),
            });
        }

        if context.box_dims.0 <= 0.0 || context.box_dims.1 <= 0.0 {
            return Err(AdaptationError::InvalidContext {
                reason: "Box dimensions must be positive".to_string(),
            });
        }

        if context.geometry_config.channel_width <= 0.0 {
            return Err(AdaptationError::InvalidContext {
                reason: "Channel width must be positive".to_string(),
            });
        }

        // Validate neighbor information if present
        if let Some(neighbors) = &context.neighbor_info {
            for (i, &neighbor_y) in neighbors.iter().enumerate() {
                if neighbor_y < 0.0 || neighbor_y > context.box_dims.1 {
                    return Err(AdaptationError::InvalidContext {
                        reason: format!(
                            "Neighbor {} position ({}) is outside box bounds (0, {})",
                            i, neighbor_y, context.box_dims.1
                        ),
                    });
                }
            }
        }

        // Validate endpoints
        let (start, end) = context.channel_endpoints;
        if start != (0.0, 0.0) || end != (0.0, 0.0) {
            if start.0 < 0.0
                || start.0 > context.box_dims.0
                || start.1 < 0.0
                || start.1 > context.box_dims.1
            {
                return Err(AdaptationError::InvalidContext {
                    reason: format!(
                        "Start point ({}, {}) is outside box bounds",
                        start.0, start.1
                    ),
                });
            }

            if end.0 < 0.0
                || end.0 > context.box_dims.0
                || end.1 < 0.0
                || end.1 > context.box_dims.1
            {
                return Err(AdaptationError::InvalidContext {
                    reason: format!("End point ({}, {}) is outside box bounds", end.0, end.1),
                });
            }
        }

        Ok(())
    }
}

/// Error recovery utilities
pub struct ErrorRecovery;

impl ErrorRecovery {
    /// Attempt to recover from a parameter error by using default values
    pub fn recover_parameter_error<T>(
        error: &ParameterError,
        default_value: T,
        parameter_name: &str,
    ) -> Result<T, ParameterError>
    where
        T: Clone + Debug,
    {
        match error {
            ParameterError::NotFound { .. } => {
                #[cfg(debug_assertions)]
                eprintln!("Warning: Parameter '{parameter_name}' not found, using default value");
                Ok(default_value)
            }
            ParameterError::InvalidValue { .. } => {
                #[cfg(debug_assertions)]
                eprintln!("Warning: Parameter '{parameter_name}' has invalid value, using default");
                Ok(default_value)
            }
            ParameterError::ReadOnly { .. } => {
                // Cannot recover from read-only errors
                Err(error.clone())
            }
            _ => {
                // For other errors, try using default
                #[cfg(debug_assertions)]
                eprintln!(
                    "Warning: Parameter '{parameter_name}' error, attempting recovery with default"
                );
                Ok(default_value)
            }
        }
    }

    /// Attempt to recover from adaptation errors
    pub fn recover_adaptation_error<T>(
        error: &AdaptationError,
        base_value: T,
        parameter_name: &str,
    ) -> Result<T, AdaptationError>
    where
        T: Clone + Debug,
    {
        match error {
            AdaptationError::InvalidContext { .. } | AdaptationError::DependencyMissing { .. } => {
                #[cfg(debug_assertions)]
                eprintln!("Warning: Adaptation failed for '{parameter_name}', using base value");
                Ok(base_value)
            }
            AdaptationError::CalculationFailed { .. } | AdaptationError::InvalidResult { .. } => {
                // These are more serious errors that shouldn't be silently recovered
                Err(error.clone())
            }
        }
    }
}

/// Comprehensive error reporting utilities
pub struct ErrorReporter;

impl ErrorReporter {
    /// Generate a comprehensive error report for debugging
    #[must_use]
    pub fn generate_error_report(
        error: &StateManagementError,
        context: Option<&ParameterErrorContext>,
    ) -> String {
        let mut report = String::new();

        report.push_str("=== State Management Error Report ===\n");
        report.push_str(&format!("Error Type: {error}\n"));

        if let Some(ctx) = context {
            report.push_str(&format!("Parameter: {}\n", ctx.parameter_name));
            report.push_str(&format!("Domain: {}\n", ctx.domain));
            report.push_str(&format!("Operation: {}\n", ctx.operation));

            if let Some(current) = &ctx.current_value {
                report.push_str(&format!("Current Value: {current}\n"));
            }

            if let Some(expected) = &ctx.expected {
                report.push_str(&format!("Expected: {expected}\n"));
            }

            if !ctx.context.is_empty() {
                report.push_str("Additional Context:\n");
                for (key, value) in &ctx.context {
                    report.push_str(&format!("  {key}: {value}\n"));
                }
            }
        }

        report.push_str("=====================================\n");

        report
    }

    /// Log error with appropriate level based on severity
    pub fn log_error(error: &StateManagementError, context: Option<&ParameterErrorContext>) {
        let report = Self::generate_error_report(error, context);

        match error {
            StateManagementError::Parameter(ParameterError::ReadOnly { .. })
            | StateManagementError::Constraint(ConstraintError::RangeViolation { .. }) => {
                // These are user errors, log as warnings
                #[cfg(debug_assertions)]
                eprintln!("Warning: {report}");
            }
            _ => {
                // These are system errors, log as errors
                #[cfg(debug_assertions)]
                eprintln!("Error: {report}");
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::GeometryConfig;

    #[test]
    fn test_parameter_error_context_creation() {
        let context = ParameterErrorContext::new("amplitude", "serpentine", "validation")
            .with_current_value("5.0")
            .with_expected("positive number")
            .with_context("channel_length", "100.0");

        let message = context.to_detailed_message();
        assert!(message.contains("amplitude"));
        assert!(message.contains("serpentine"));
        assert!(message.contains("validation"));
        assert!(message.contains("5.0"));
        assert!(message.contains("positive number"));
        assert!(message.contains("channel_length=100.0"));
    }

    #[test]
    fn test_validation_result_success() {
        let context = ParameterErrorContext::new("test", "domain", "operation");
        let result = ValidationResult::success(context);

        assert!(result.is_valid);
        assert!(result.errors.is_empty());
        assert!(result.into_result().is_ok());
    }

    #[test]
    fn test_adaptation_context_validation() {
        let context = ChannelGenerationContext::new(
            GeometryConfig::default(),
            (100.0, 50.0),
            4,
            Some(&[10.0, 20.0, 30.0, 40.0]),
        );

        let result = ParameterValidator::validate_adaptation_context(&context, "test_param");
        assert!(result.is_ok());
    }

    #[test]
    fn test_error_recovery() {
        let error = ParameterError::not_found("missing_param", "test_domain");
        let recovered = ErrorRecovery::recover_parameter_error(&error, 42.0, "missing_param");

        assert!(recovered.is_ok());
        assert_eq!(recovered.unwrap(), 42.0);
    }
}
