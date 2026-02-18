//! Parameter constraint system
//!
//! This module provides a flexible constraint system for validating parameter values.
//! It supports range constraints, set constraints, custom validation functions,
//! and composite constraints that can be combined for complex validation logic.

use crate::state_management::errors::{ConstraintError, ParameterResult};
use std::fmt::Debug;

/// Trait for types that can be used in constraints
pub trait ConstraintValue: Clone + Debug + PartialEq + PartialOrd {}

impl<T> ConstraintValue for T where T: Clone + Debug + PartialEq + PartialOrd {}

/// Parameter constraint system for validation
#[derive(Debug, Clone)]
pub enum ParameterConstraints<T>
where
    T: ConstraintValue,
{
    /// Value must be within specified range (inclusive)
    Range { min: T, max: T },

    /// Value must be one of the specified values
    Set(Vec<T>),

    /// Custom validation function with parameters
    Custom {
        name: String,
        validator: fn(&T) -> Result<(), String>,
    },

    /// Length range constraint for strings
    LengthRange { min: usize, max: usize },

    /// Contains pattern constraint for strings
    Contains { pattern: String },

    /// Composite constraint (all must pass)
    All(Vec<ParameterConstraints<T>>),

    /// Alternative constraint (at least one must pass)
    Any(Vec<ParameterConstraints<T>>),

    /// No constraints (always valid)
    None,
}

impl<T> ParameterConstraints<T>
where
    T: ConstraintValue + 'static,
{
    /// Create a range constraint
    pub const fn range(min: T, max: T) -> Self {
        Self::Range { min, max }
    }

    /// Create a set constraint
    #[must_use]
    pub const fn set(values: Vec<T>) -> Self {
        Self::Set(values)
    }

    /// Create a custom constraint
    pub fn custom(name: &str, validator: fn(&T) -> Result<(), String>) -> Self {
        Self::Custom {
            name: name.to_string(),
            validator,
        }
    }

    /// Create a composite constraint where all must pass
    #[must_use]
    pub const fn all(constraints: Vec<Self>) -> Self {
        Self::All(constraints)
    }

    /// Create a composite constraint where any can pass
    #[must_use]
    pub const fn any(constraints: Vec<Self>) -> Self {
        Self::Any(constraints)
    }

    /// Create no constraints
    #[must_use]
    pub const fn none() -> Self {
        Self::None
    }

    /// Validate a value against this constraint
    pub fn validate(&self, value: &T) -> ParameterResult<()> {
        match self {
            Self::Range { min, max } => {
                if value < min || value > max {
                    Err(ConstraintError::range_violation(value, min, max).into())
                } else {
                    Ok(())
                }
            }

            Self::Set(allowed_values) => {
                if allowed_values.contains(value) {
                    Ok(())
                } else {
                    let allowed_debug: Vec<&dyn Debug> =
                        allowed_values.iter().map(|v| v as &dyn Debug).collect();
                    Err(ConstraintError::set_violation(value, &allowed_debug).into())
                }
            }

            Self::Custom { name, validator } => validator(value)
                .map_err(|msg| ConstraintError::custom_violation(&format!("{name}: {msg}")).into()),

            Self::LengthRange { min, max } => {
                // This only works for String types
                if let Some(string_val) = (value as &dyn std::any::Any).downcast_ref::<String>() {
                    let len = string_val.len();
                    if len >= *min && len <= *max {
                        Ok(())
                    } else {
                        Err(ConstraintError::custom_violation(&format!(
                            "length {len} is not between {min} and {max}"
                        ))
                        .into())
                    }
                } else {
                    Ok(()) // Skip validation for non-string types
                }
            }

            Self::Contains { pattern } => {
                // This only works for String types
                if let Some(string_val) = (value as &dyn std::any::Any).downcast_ref::<String>() {
                    if string_val.contains(pattern) {
                        Ok(())
                    } else {
                        Err(ConstraintError::custom_violation(&format!(
                            "string does not contain '{pattern}'"
                        ))
                        .into())
                    }
                } else {
                    Ok(()) // Skip validation for non-string types
                }
            }

            Self::All(constraints) => {
                for constraint in constraints {
                    constraint.validate(value)?;
                }
                Ok(())
            }

            Self::Any(constraints) => {
                if constraints.is_empty() {
                    return Ok(());
                }

                let mut errors = Vec::new();
                for constraint in constraints {
                    match constraint.validate(value) {
                        Ok(()) => return Ok(()),
                        Err(e) => errors.push(e),
                    }
                }

                // If we get here, all constraints failed
                Err(ConstraintError::custom_violation(&format!(
                    "All alternative constraints failed: {errors:?}"
                ))
                .into())
            }

            Self::None => Ok(()),
        }
    }

    /// Check if this constraint is satisfied by a value
    pub fn is_satisfied_by(&self, value: &T) -> bool {
        self.validate(value).is_ok()
    }

    /// Get a human-readable description of this constraint
    pub fn description(&self) -> String {
        match self {
            Self::Range { min, max } => format!("must be between {min:?} and {max:?}"),
            Self::Set(values) => format!("must be one of {values:?}"),
            Self::Custom { name, .. } => format!("must satisfy {name}"),
            Self::LengthRange { min, max } => format!("length must be between {min} and {max}"),
            Self::Contains { pattern } => format!("must contain '{pattern}'"),
            Self::All(constraints) => {
                let descriptions: Vec<String> = constraints.iter().map(Self::description).collect();
                format!("must satisfy all: [{}]", descriptions.join(", "))
            }
            Self::Any(constraints) => {
                let descriptions: Vec<String> = constraints.iter().map(Self::description).collect();
                format!("must satisfy any: [{}]", descriptions.join(" OR "))
            }
            Self::None => "no constraints".to_string(),
        }
    }
}

/// Common constraint builders for numeric types
impl ParameterConstraints<f64> {
    /// Create a positive constraint (> 0.0)
    #[must_use]
    pub fn positive() -> Self {
        Self::custom("positive", |v| {
            if *v > 0.0 {
                Ok(())
            } else {
                Err("must be positive".to_string())
            }
        })
    }

    /// Create a non-negative constraint (>= 0.0)
    #[must_use]
    pub fn non_negative() -> Self {
        Self::custom("non-negative", |v| {
            if *v >= 0.0 {
                Ok(())
            } else {
                Err("must be non-negative".to_string())
            }
        })
    }

    /// Create a non-zero constraint (!= 0.0)
    #[must_use]
    pub fn non_zero() -> Self {
        Self::custom("non-zero", |v| {
            if v.abs() > f64::EPSILON {
                Ok(())
            } else {
                Err("must be non-zero".to_string())
            }
        })
    }

    /// Create a normalized constraint (0.0 <= value <= 1.0)
    #[must_use]
    pub const fn normalized() -> Self {
        Self::range(0.0, 1.0)
    }

    /// Create a percentage constraint (0.0 <= value <= 100.0)
    #[must_use]
    pub const fn percentage() -> Self {
        Self::range(0.0, 100.0)
    }
}

/// Common constraint builders for integer types
impl ParameterConstraints<usize> {
    /// Create a positive constraint (> 0)
    #[must_use]
    pub fn positive() -> Self {
        Self::custom("positive", |v| {
            if *v > 0 {
                Ok(())
            } else {
                Err("must be positive".to_string())
            }
        })
    }
}

/// Common constraint builders for u32 types
impl ParameterConstraints<u32> {
    /// Create a positive constraint (> 0)
    #[must_use]
    pub fn positive() -> Self {
        Self::custom("positive", |v| {
            if *v > 0 {
                Ok(())
            } else {
                Err("must be positive".to_string())
            }
        })
    }

    /// Create a non-zero constraint (!= 0)
    #[must_use]
    pub fn non_zero() -> Self {
        Self::positive() // For usize, positive is the same as non-zero
    }

    /// Create a power of two constraint
    #[must_use]
    pub fn power_of_two() -> Self {
        Self::custom("power-of-two", |v| {
            if *v > 0 && (*v).is_power_of_two() {
                Ok(())
            } else {
                Err("must be a power of two".to_string())
            }
        })
    }
}

/// Common constraint builders for string types
impl ParameterConstraints<String> {
    /// Create a non-empty constraint
    #[must_use]
    pub fn non_empty() -> Self {
        Self::custom("non-empty", |v| {
            if v.is_empty() {
                Err("must not be empty".to_string())
            } else {
                Ok(())
            }
        })
    }

    /// Create a length constraint
    #[must_use]
    pub const fn length_range(min: usize, max: usize) -> Self {
        Self::LengthRange { min, max }
    }

    /// Create a pattern constraint (simple contains check)
    #[must_use]
    pub fn contains(pattern: &str) -> Self {
        Self::Contains {
            pattern: pattern.to_string(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_range_constraint() {
        let constraint = ParameterConstraints::range(0.0, 10.0);

        assert!(constraint.validate(&5.0).is_ok());
        assert!(constraint.validate(&0.0).is_ok());
        assert!(constraint.validate(&10.0).is_ok());
        assert!(constraint.validate(&-1.0).is_err());
        assert!(constraint.validate(&11.0).is_err());
    }

    #[test]
    fn test_set_constraint() {
        let constraint = ParameterConstraints::set(vec![1, 2, 3]);

        assert!(constraint.validate(&1).is_ok());
        assert!(constraint.validate(&2).is_ok());
        assert!(constraint.validate(&3).is_ok());
        assert!(constraint.validate(&4).is_err());
    }

    #[test]
    fn test_positive_constraint() {
        let constraint = ParameterConstraints::<f64>::positive();

        assert!(constraint.validate(&1.0).is_ok());
        assert!(constraint.validate(&0.1).is_ok());
        assert!(constraint.validate(&0.0).is_err());
        assert!(constraint.validate(&-1.0).is_err());
    }

    #[test]
    fn test_composite_constraint() {
        let constraint = ParameterConstraints::all(vec![
            ParameterConstraints::<f64>::positive(),
            ParameterConstraints::range(0.0, 100.0),
        ]);

        assert!(constraint.validate(&50.0).is_ok());
        assert!(constraint.validate(&0.0).is_err()); // Fails positive
        assert!(constraint.validate(&150.0).is_err()); // Fails range
    }

    #[test]
    fn test_constraint_description() {
        let constraint = ParameterConstraints::range(0.0, 10.0);
        let description = constraint.description();
        assert!(description.contains("between"));
        assert!(description.contains("0"));
        assert!(description.contains("10"));
    }
}
