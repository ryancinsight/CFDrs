//! Symmetry parameter manager.
//!
//! Manages parameters that control bilateral mirror symmetry enforcement
//! for millifluidic channel layouts. Parameters mirror the
//! [`BilateralSymmetryConfig`](crate::state_management::bilateral_symmetry::BilateralSymmetryConfig)
//! struct used at runtime.

use crate::state_management::{
    constraints::ParameterConstraints,
    errors::{ParameterError, ParameterResult},
    parameters::{ConfigurableParameter, ParameterMetadata},
    validation::ValidationRuleSet,
};

use super::ParameterManager;

/// Parameter manager for bilateral mirror symmetry enforcement
///
/// Manages the following parameters:
///
/// | Name | Type | Default | Range | Units |
/// |------|------|---------|-------|-------|
/// | `enable_vertical_symmetry` | bool | true | — | — |
/// | `enable_horizontal_symmetry` | bool | true | — | — |
/// | `symmetry_tolerance` | f64 | 1e-6 | 1e-12–1.0 | mm |
/// | `enable_adaptive_symmetry` | bool | true | — | — |
/// | `enforcement_strength` | f64 | 1.0 | 0.0–1.0 | — |
#[derive(Debug)]
pub struct SymmetryParameterManager {
    /// Enable perfect vertical centerline symmetry (splits mirror merges)
    enable_vertical_symmetry: ConfigurableParameter<bool>,

    /// Enable perfect horizontal centerline symmetry (upper mirrors lower)
    enable_horizontal_symmetry: ConfigurableParameter<bool>,

    /// Positional tolerance for symmetry validation
    symmetry_tolerance: ConfigurableParameter<f64>,

    /// Enable context-aware adaptive symmetry adjustments
    enable_adaptive_symmetry: ConfigurableParameter<bool>,

    /// Symmetry enforcement strength (0.0 = weak, 1.0 = strict)
    enforcement_strength: ConfigurableParameter<f64>,

    /// Validation rules
    validation_rules: ValidationRuleSet,
}

impl SymmetryParameterManager {
    /// Create a new symmetry parameter manager with defaults from `BilateralSymmetryConfig`
    #[must_use]
    pub fn new() -> Self {
        let enable_vertical_symmetry = ConfigurableParameter::new(
            true,
            ParameterConstraints::none(),
            ParameterMetadata::new(
                "enable_vertical_symmetry",
                "Enable perfect vertical centerline mirror symmetry (splits mirror merges)",
                "symmetry_parameters",
            ),
        );

        let enable_horizontal_symmetry = ConfigurableParameter::new(
            true,
            ParameterConstraints::none(),
            ParameterMetadata::new(
                "enable_horizontal_symmetry",
                "Enable perfect horizontal centerline mirror symmetry (upper mirrors lower)",
                "symmetry_parameters",
            ),
        );

        let symmetry_tolerance = ConfigurableParameter::new(
            1e-6,
            ParameterConstraints::all(vec![
                ParameterConstraints::<f64>::positive(),
                ParameterConstraints::range(1e-12, 1.0),
            ]),
            ParameterMetadata::new(
                "symmetry_tolerance",
                "Positional tolerance for symmetry validation checks",
                "symmetry_parameters",
            )
            .with_units("mm"),
        );

        let enable_adaptive_symmetry = ConfigurableParameter::new(
            true,
            ParameterConstraints::none(),
            ParameterMetadata::new(
                "enable_adaptive_symmetry",
                "Enable context-aware adaptive symmetry adjustments",
                "symmetry_parameters",
            ),
        );

        let enforcement_strength = ConfigurableParameter::new(
            1.0,
            ParameterConstraints::all(vec![
                ParameterConstraints::non_negative(),
                ParameterConstraints::normalized(),
            ]),
            ParameterMetadata::new(
                "enforcement_strength",
                "Symmetry enforcement strength (0.0 = weak, 1.0 = strict)",
                "symmetry_parameters",
            ),
        );

        let validation_rules = ValidationRuleSet::new();

        Self {
            enable_vertical_symmetry,
            enable_horizontal_symmetry,
            symmetry_tolerance,
            enable_adaptive_symmetry,
            enforcement_strength,
            validation_rules,
        }
    }

    /// Check if vertical symmetry is enabled
    #[must_use]
    pub const fn is_vertical_symmetry_enabled(&self) -> bool {
        *self.enable_vertical_symmetry.get_raw_value()
    }

    /// Check if horizontal symmetry is enabled
    #[must_use]
    pub const fn is_horizontal_symmetry_enabled(&self) -> bool {
        *self.enable_horizontal_symmetry.get_raw_value()
    }

    /// Get symmetry tolerance
    #[must_use]
    pub const fn get_symmetry_tolerance(&self) -> f64 {
        *self.symmetry_tolerance.get_raw_value()
    }

    /// Get enforcement strength
    #[must_use]
    pub const fn get_enforcement_strength(&self) -> f64 {
        *self.enforcement_strength.get_raw_value()
    }
}

impl Default for SymmetryParameterManager {
    fn default() -> Self {
        Self::new()
    }
}

const SYMMETRY_PARAMS: &[&str] = &[
    "enable_vertical_symmetry",
    "enable_horizontal_symmetry",
    "symmetry_tolerance",
    "enable_adaptive_symmetry",
    "enforcement_strength",
];

impl ParameterManager for SymmetryParameterManager {
    fn get_parameter(&self, name: &str) -> ParameterResult<Box<dyn std::any::Any>> {
        match name {
            "enable_vertical_symmetry" => {
                Ok(Box::new(*self.enable_vertical_symmetry.get_raw_value()))
            }
            "enable_horizontal_symmetry" => {
                Ok(Box::new(*self.enable_horizontal_symmetry.get_raw_value()))
            }
            "symmetry_tolerance" => Ok(Box::new(*self.symmetry_tolerance.get_raw_value())),
            "enable_adaptive_symmetry" => {
                Ok(Box::new(*self.enable_adaptive_symmetry.get_raw_value()))
            }
            "enforcement_strength" => Ok(Box::new(*self.enforcement_strength.get_raw_value())),
            _ => Err(ParameterError::not_found(name, "symmetry")),
        }
    }

    fn set_parameter(
        &mut self,
        name: &str,
        value: Box<dyn std::any::Any>,
        reason: &str,
    ) -> ParameterResult<()> {
        match name {
            "enable_vertical_symmetry" => {
                if let Some(val) = value.downcast_ref::<bool>() {
                    self.enable_vertical_symmetry.set_value(*val, reason)
                } else {
                    Err(ParameterError::type_mismatch(name, "bool", "unknown"))
                }
            }
            "enable_horizontal_symmetry" => {
                if let Some(val) = value.downcast_ref::<bool>() {
                    self.enable_horizontal_symmetry.set_value(*val, reason)
                } else {
                    Err(ParameterError::type_mismatch(name, "bool", "unknown"))
                }
            }
            "symmetry_tolerance" => {
                if let Some(val) = value.downcast_ref::<f64>() {
                    self.symmetry_tolerance.set_value(*val, reason)
                } else {
                    Err(ParameterError::type_mismatch(name, "f64", "unknown"))
                }
            }
            "enable_adaptive_symmetry" => {
                if let Some(val) = value.downcast_ref::<bool>() {
                    self.enable_adaptive_symmetry.set_value(*val, reason)
                } else {
                    Err(ParameterError::type_mismatch(name, "bool", "unknown"))
                }
            }
            "enforcement_strength" => {
                if let Some(val) = value.downcast_ref::<f64>() {
                    self.enforcement_strength.set_value(*val, reason)
                } else {
                    Err(ParameterError::type_mismatch(name, "f64", "unknown"))
                }
            }
            _ => Err(ParameterError::not_found(name, "symmetry")),
        }
    }

    fn parameter_names(&self) -> Vec<String> {
        SYMMETRY_PARAMS.iter().map(|s| (*s).to_string()).collect()
    }

    fn has_parameter(&self, name: &str) -> bool {
        SYMMETRY_PARAMS.contains(&name)
    }

    fn validate_all(&self) -> ParameterResult<()> {
        self.enable_vertical_symmetry.validate()?;
        self.enable_horizontal_symmetry.validate()?;
        self.symmetry_tolerance.validate()?;
        self.enable_adaptive_symmetry.validate()?;
        self.enforcement_strength.validate()?;
        Ok(())
    }

    fn get_metadata(&self, name: &str) -> ParameterResult<&ParameterMetadata> {
        match name {
            "enable_vertical_symmetry" => Ok(self.enable_vertical_symmetry.metadata()),
            "enable_horizontal_symmetry" => Ok(self.enable_horizontal_symmetry.metadata()),
            "symmetry_tolerance" => Ok(self.symmetry_tolerance.metadata()),
            "enable_adaptive_symmetry" => Ok(self.enable_adaptive_symmetry.metadata()),
            "enforcement_strength" => Ok(self.enforcement_strength.metadata()),
            _ => Err(ParameterError::not_found(name, "symmetry")),
        }
    }

    fn domain_name(&self) -> &'static str {
        "symmetry"
    }

    fn reset_all(&mut self, reason: &str) -> ParameterResult<()> {
        self.enable_vertical_symmetry.reset(reason)?;
        self.enable_horizontal_symmetry.reset(reason)?;
        self.symmetry_tolerance.reset(reason)?;
        self.enable_adaptive_symmetry.reset(reason)?;
        self.enforcement_strength.reset(reason)?;
        Ok(())
    }

    fn validation_rules(&self) -> &ValidationRuleSet {
        &self.validation_rules
    }
}
