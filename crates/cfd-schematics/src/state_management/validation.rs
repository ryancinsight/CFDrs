//! Parameter validation system
//!
//! This module provides a comprehensive validation framework for parameter
//! relationships, dependencies, and cross-parameter constraints.

use crate::state_management::errors::{ParameterResult, ValidationError};
use std::collections::{HashMap, HashSet};
use std::fmt::Debug;

/// Type alias for complex validation function signatures
pub type ValidationFunction = fn(&HashMap<String, Box<dyn std::any::Any>>) -> Result<(), String>;

/// Validation rule for parameter relationships
pub trait ValidationRule: Debug + Send + Sync {
    /// Validate parameters and return any errors
    fn validate(&self, parameters: &HashMap<String, Box<dyn std::any::Any>>)
        -> ParameterResult<()>;

    /// Get the name of this validation rule
    fn name(&self) -> &str;

    /// Get the parameters this rule depends on
    fn dependencies(&self) -> Vec<String>;

    /// Get a description of what this rule validates
    fn description(&self) -> String;
}

/// Set of validation rules for a parameter domain
#[derive(Debug)]
pub struct ValidationRuleSet {
    /// Individual validation rules
    rules: Vec<Box<dyn ValidationRule>>,

    /// Whether validation is enabled
    enabled: bool,

    /// Cache of validation results
    cache: HashMap<String, bool>,
}

impl ValidationRuleSet {
    /// Create a new validation rule set
    #[must_use]
    pub fn new() -> Self {
        Self {
            rules: Vec::new(),
            enabled: true,
            cache: HashMap::new(),
        }
    }

    /// Add a validation rule
    #[must_use]
    pub fn add_rule<R: ValidationRule + 'static>(mut self, rule: R) -> Self {
        self.rules.push(Box::new(rule));
        self.cache.clear(); // Clear cache when rules change
        self
    }

    /// Enable or disable validation
    pub fn set_enabled(&mut self, enabled: bool) {
        self.enabled = enabled;
        if !enabled {
            self.cache.clear();
        }
    }

    /// Check if validation is enabled
    #[must_use]
    pub const fn is_enabled(&self) -> bool {
        self.enabled
    }

    /// Validate all rules against parameters
    ///
    /// # Errors
    ///
    /// Returns an error if any validation rule fails. The error contains details
    /// about the first validation failure encountered.
    ///
    /// # Panics
    ///
    /// This method will panic if the error collection is empty when errors are expected.
    /// This should never happen under normal usage.
    pub fn validate_all(
        &self,
        parameters: &HashMap<String, Box<dyn std::any::Any>>,
    ) -> ParameterResult<()> {
        if !self.enabled {
            return Ok(());
        }

        let mut errors = Vec::new();

        for rule in &self.rules {
            if let Err(e) = rule.validate(parameters) {
                errors.push(e);
            }
        }

        if errors.is_empty() {
            Ok(())
        } else if errors.len() == 1 {
            Err(errors.into_iter().next().unwrap())
        } else {
            Err(ValidationError::multiple(
                errors
                    .into_iter()
                    .map(|e| ValidationError::custom(&e.to_string()))
                    .collect(),
            )
            .into())
        }
    }

    /// Get all rule names
    #[must_use]
    pub fn rule_names(&self) -> Vec<&str> {
        self.rules.iter().map(|r| r.name()).collect()
    }

    /// Get all dependencies across all rules
    #[must_use]
    pub fn all_dependencies(&self) -> HashSet<String> {
        self.rules.iter().flat_map(|r| r.dependencies()).collect()
    }

    /// Clear validation cache
    pub fn clear_cache(&mut self) {
        self.cache.clear();
    }
}

impl Default for ValidationRuleSet {
    fn default() -> Self {
        Self::new()
    }
}

/// Range validation rule for numeric parameters
#[derive(Debug)]
pub struct RangeValidationRule {
    name: String,
    parameter_name: String,
    min_value: f64,
    max_value: f64,
}

impl RangeValidationRule {
    /// Create a new range validation rule
    #[must_use]
    pub fn new(name: &str, parameter_name: &str, min_value: f64, max_value: f64) -> Self {
        Self {
            name: name.to_string(),
            parameter_name: parameter_name.to_string(),
            min_value,
            max_value,
        }
    }
}

impl ValidationRule for RangeValidationRule {
    fn validate(
        &self,
        parameters: &HashMap<String, Box<dyn std::any::Any>>,
    ) -> ParameterResult<()> {
        if let Some(param) = parameters.get(&self.parameter_name) {
            if let Some(value) = param.downcast_ref::<f64>() {
                if *value < self.min_value || *value > self.max_value {
                    return Err(ValidationError::rule_failed(
                        &self.parameter_name,
                        &format!(
                            "Value {} is outside range [{}, {}]",
                            value, self.min_value, self.max_value
                        ),
                    )
                    .into());
                }
            }
        }
        Ok(())
    }

    fn name(&self) -> &str {
        &self.name
    }

    fn dependencies(&self) -> Vec<String> {
        vec![self.parameter_name.clone()]
    }

    fn description(&self) -> String {
        format!(
            "Validates that {} is between {} and {}",
            self.parameter_name, self.min_value, self.max_value
        )
    }
}

/// Relationship validation rule for parameter dependencies
#[derive(Debug)]
pub struct RelationshipValidationRule {
    name: String,
    primary_param: String,
    dependent_param: String,
    validator: fn(f64, f64) -> bool,
    error_message: String,
}

impl RelationshipValidationRule {
    /// Create a new relationship validation rule
    pub fn new(
        name: &str,
        primary_param: &str,
        dependent_param: &str,
        validator: fn(f64, f64) -> bool,
        error_message: &str,
    ) -> Self {
        Self {
            name: name.to_string(),
            primary_param: primary_param.to_string(),
            dependent_param: dependent_param.to_string(),
            validator,
            error_message: error_message.to_string(),
        }
    }

    /// Create a rule where dependent must be less than primary
    #[must_use]
    pub fn less_than(name: &str, primary: &str, dependent: &str) -> Self {
        Self::new(
            name,
            primary,
            dependent,
            |p, d| d < p,
            &format!("{dependent} must be less than {primary}"),
        )
    }

    /// Create a rule where dependent must be greater than primary
    #[must_use]
    pub fn greater_than(name: &str, primary: &str, dependent: &str) -> Self {
        Self::new(
            name,
            primary,
            dependent,
            |p, d| d > p,
            &format!("{dependent} must be greater than {primary}"),
        )
    }

    /// Create a rule where dependent must be a multiple of primary
    #[must_use]
    pub fn multiple_of(name: &str, primary: &str, dependent: &str) -> Self {
        Self::new(
            name,
            primary,
            dependent,
            |p, d| p > 0.0 && (d % p).abs() < 1e-10,
            &format!("{dependent} must be a multiple of {primary}"),
        )
    }
}

impl ValidationRule for RelationshipValidationRule {
    fn validate(
        &self,
        parameters: &HashMap<String, Box<dyn std::any::Any>>,
    ) -> ParameterResult<()> {
        let primary_value = parameters
            .get(&self.primary_param)
            .and_then(|p| p.downcast_ref::<f64>());
        let dependent_value = parameters
            .get(&self.dependent_param)
            .and_then(|p| p.downcast_ref::<f64>());

        if let (Some(&primary), Some(&dependent)) = (primary_value, dependent_value) {
            if !(self.validator)(primary, dependent) {
                return Err(ValidationError::rule_failed(
                    &format!("{}+{}", self.primary_param, self.dependent_param),
                    &self.error_message,
                )
                .into());
            }
        }

        Ok(())
    }

    fn name(&self) -> &str {
        &self.name
    }

    fn dependencies(&self) -> Vec<String> {
        vec![self.primary_param.clone(), self.dependent_param.clone()]
    }

    fn description(&self) -> String {
        format!(
            "Validates relationship between {} and {}: {}",
            self.primary_param, self.dependent_param, self.error_message
        )
    }
}

/// Custom validation rule with user-defined logic
#[derive(Debug)]
pub struct CustomValidationRule {
    name: String,
    dependencies: Vec<String>,
    validator: ValidationFunction,
    description: String,
}

impl CustomValidationRule {
    /// Create a new custom validation rule
    pub fn new(
        name: &str,
        dependencies: Vec<String>,
        validator: ValidationFunction,
        description: &str,
    ) -> Self {
        Self {
            name: name.to_string(),
            dependencies,
            validator,
            description: description.to_string(),
        }
    }
}

impl ValidationRule for CustomValidationRule {
    fn validate(
        &self,
        parameters: &HashMap<String, Box<dyn std::any::Any>>,
    ) -> ParameterResult<()> {
        (self.validator)(parameters)
            .map_err(|msg| ValidationError::custom(&format!("{}: {}", self.name, msg)).into())
    }

    fn name(&self) -> &str {
        &self.name
    }

    fn dependencies(&self) -> Vec<String> {
        self.dependencies.clone()
    }

    fn description(&self) -> String {
        self.description.clone()
    }
}

/// Common validation rules for microfluidic parameters
pub struct CommonValidationRules;

impl CommonValidationRules {
    /// Create validation rules for serpentine parameters
    #[must_use]
    pub fn serpentine_rules() -> ValidationRuleSet {
        ValidationRuleSet::new()
            .add_rule(RangeValidationRule::new(
                "amplitude_range",
                "amplitude",
                0.1,
                100.0,
            ))
            .add_rule(RangeValidationRule::new(
                "wavelength_range",
                "wavelength_factor",
                0.1,
                10.0,
            ))
            .add_rule(RelationshipValidationRule::less_than(
                "amplitude_wavelength_ratio",
                "wavelength_factor",
                "amplitude",
            ))
    }

    /// Create validation rules for arc parameters
    #[must_use]
    pub fn arc_rules() -> ValidationRuleSet {
        ValidationRuleSet::new()
            .add_rule(RangeValidationRule::new(
                "curvature_range",
                "curvature_factor",
                0.0,
                2.0,
            ))
            .add_rule(RangeValidationRule::new(
                "smoothness_range",
                "smoothness",
                5.0,
                100.0,
            ))
    }

    /// Create validation rules for geometry parameters
    #[must_use]
    pub fn geometry_rules() -> ValidationRuleSet {
        ValidationRuleSet::new()
            .add_rule(RangeValidationRule::new(
                "wall_clearance_range",
                "wall_clearance",
                0.01,
                10.0,
            ))
            .add_rule(RangeValidationRule::new(
                "channel_width_range",
                "channel_width",
                0.01,
                50.0,
            ))
            .add_rule(RelationshipValidationRule::less_than(
                "clearance_width_ratio",
                "channel_width",
                "wall_clearance",
            ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::any::Any;

    fn create_test_parameters() -> HashMap<String, Box<dyn Any>> {
        let mut params = HashMap::new();
        params.insert("amplitude".to_string(), Box::new(1.5) as Box<dyn Any>);
        params.insert(
            "wavelength_factor".to_string(),
            Box::new(2.0) as Box<dyn Any>,
        );
        params.insert("channel_width".to_string(), Box::new(1.0) as Box<dyn Any>);
        params.insert("wall_clearance".to_string(), Box::new(0.5) as Box<dyn Any>);
        params
    }

    #[test]
    fn test_range_validation_rule() {
        let rule = RangeValidationRule::new("test_range", "amplitude", 1.0, 10.0);
        let params = create_test_parameters();

        assert!(rule.validate(&params).is_ok());

        let mut bad_params = HashMap::new();
        bad_params.insert("amplitude".to_string(), Box::new(15.0) as Box<dyn Any>);
        bad_params.insert(
            "wavelength_factor".to_string(),
            Box::new(2.0) as Box<dyn Any>,
        );
        bad_params.insert("channel_width".to_string(), Box::new(1.0) as Box<dyn Any>);
        bad_params.insert("wall_clearance".to_string(), Box::new(0.5) as Box<dyn Any>);
        assert!(rule.validate(&bad_params).is_err());
    }

    #[test]
    fn test_relationship_validation_rule() {
        let rule = RelationshipValidationRule::less_than(
            "test_relationship",
            "wavelength_factor",
            "amplitude",
        );
        let params = create_test_parameters();

        // amplitude (1.5) < wavelength_factor (2.0), so this should pass
        assert!(rule.validate(&params).is_ok());

        let mut bad_params = HashMap::new();
        bad_params.insert("amplitude".to_string(), Box::new(3.0) as Box<dyn Any>);
        bad_params.insert(
            "wavelength_factor".to_string(),
            Box::new(2.0) as Box<dyn Any>,
        );
        bad_params.insert("channel_width".to_string(), Box::new(1.0) as Box<dyn Any>);
        bad_params.insert("wall_clearance".to_string(), Box::new(0.5) as Box<dyn Any>);
        // amplitude (3.0) > wavelength_factor (2.0), so this should fail
        assert!(rule.validate(&bad_params).is_err());
    }

    #[test]
    fn test_validation_rule_set() {
        let rule_set = CommonValidationRules::serpentine_rules();
        let params = create_test_parameters();

        assert!(rule_set.validate_all(&params).is_ok());
        assert!(rule_set.is_enabled());
        assert!(!rule_set.rule_names().is_empty());
    }
}
