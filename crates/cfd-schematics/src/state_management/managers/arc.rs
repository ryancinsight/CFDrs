//! Arc channel parameter manager.

use crate::state_management::{
    constraints::ParameterConstraints,
    errors::{ParameterError, ParameterResult},
    parameters::{ConfigurableParameter, ParameterMetadata},
    validation::ValidationRuleSet,
};

use super::ParameterManager;

/// Parameter manager for arc channel parameters
#[derive(Debug)]
pub struct ArcParameterManager {
    /// Base curvature factor
    curvature_factor: ConfigurableParameter<f64>,

    /// Number of smoothness points
    smoothness: ConfigurableParameter<usize>,

    /// Curvature direction (-1.0 to 1.0)
    curvature_direction: ConfigurableParameter<f64>,

    /// Minimum separation distance for collision prevention
    min_separation_distance: ConfigurableParameter<f64>,

    /// Maximum curvature reduction factor
    max_curvature_reduction: ConfigurableParameter<f64>,

    /// Enable collision prevention
    enable_collision_prevention: ConfigurableParameter<bool>,

    /// Enable adaptive curvature
    enable_adaptive_curvature: ConfigurableParameter<bool>,

    /// Validation rules
    validation_rules: ValidationRuleSet,
}

impl ArcParameterManager {
    /// Create a new arc parameter manager with default values
    #[must_use]
    pub fn new() -> Self {
        let curvature_factor = ConfigurableParameter::new(
            0.5,
            ParameterConstraints::all(vec![
                ParameterConstraints::non_negative(),
                ParameterConstraints::range(0.0, 2.0),
            ]),
            ParameterMetadata::new(
                "curvature_factor",
                "Base curvature factor for arc generation",
                "arc_parameters",
            )
            .affects_others(),
        );

        let smoothness = ConfigurableParameter::new(
            20usize,
            ParameterConstraints::all(vec![
                ParameterConstraints::<usize>::positive(),
                ParameterConstraints::range(5, 100),
            ]),
            ParameterMetadata::new(
                "smoothness",
                "Number of points for arc smoothness",
                "arc_parameters",
            ),
        );

        let curvature_direction = ConfigurableParameter::new(
            0.0,
            ParameterConstraints::range(-1.0, 1.0),
            ParameterMetadata::new(
                "curvature_direction",
                "Direction of arc curvature (-1.0 to 1.0, 0.0 for auto)",
                "arc_parameters",
            ),
        );

        let min_separation_distance = ConfigurableParameter::new(
            2.0,
            ParameterConstraints::<f64>::positive(),
            ParameterMetadata::new(
                "min_separation_distance",
                "Minimum separation distance for collision prevention",
                "collision_parameters",
            )
            .with_units("mm"),
        );

        let max_curvature_reduction = ConfigurableParameter::new(
            0.8,
            ParameterConstraints::all(vec![
                ParameterConstraints::<f64>::positive(),
                ParameterConstraints::normalized(),
            ]),
            ParameterMetadata::new(
                "max_curvature_reduction",
                "Maximum allowed curvature reduction factor",
                "collision_parameters",
            ),
        );

        let enable_collision_prevention = ConfigurableParameter::new(
            true,
            ParameterConstraints::none(),
            ParameterMetadata::new(
                "enable_collision_prevention",
                "Enable collision prevention for arc channels",
                "collision_parameters",
            ),
        );

        let enable_adaptive_curvature = ConfigurableParameter::new(
            true,
            ParameterConstraints::none(),
            ParameterMetadata::new(
                "enable_adaptive_curvature",
                "Enable adaptive curvature based on context",
                "adaptive_parameters",
            ),
        );

        let validation_rules = ValidationRuleSet::new();

        Self {
            curvature_factor,
            smoothness,
            curvature_direction,
            min_separation_distance,
            max_curvature_reduction,
            enable_collision_prevention,
            enable_adaptive_curvature,
            validation_rules,
        }
    }

    /// Get curvature factor
    #[must_use]
    pub const fn get_curvature_factor(&self) -> f64 {
        *self.curvature_factor.get_raw_value()
    }

    /// Get smoothness points
    #[must_use]
    pub const fn get_smoothness(&self) -> usize {
        *self.smoothness.get_raw_value()
    }

    /// Check if collision prevention is enabled
    #[must_use]
    pub const fn is_collision_prevention_enabled(&self) -> bool {
        *self.enable_collision_prevention.get_raw_value()
    }

    /// Check if adaptive curvature is enabled
    #[must_use]
    pub const fn is_adaptive_curvature_enabled(&self) -> bool {
        *self.enable_adaptive_curvature.get_raw_value()
    }
}

impl Default for ArcParameterManager {
    fn default() -> Self {
        Self::new()
    }
}

impl ParameterManager for ArcParameterManager {
    fn get_parameter(&self, name: &str) -> ParameterResult<Box<dyn std::any::Any>> {
        match name {
            "curvature_factor" => Ok(Box::new(*self.curvature_factor.get_raw_value())),
            "smoothness" => Ok(Box::new(*self.smoothness.get_raw_value())),
            "curvature_direction" => Ok(Box::new(*self.curvature_direction.get_raw_value())),
            "min_separation_distance" => {
                Ok(Box::new(*self.min_separation_distance.get_raw_value()))
            }
            "max_curvature_reduction" => {
                Ok(Box::new(*self.max_curvature_reduction.get_raw_value()))
            }
            "enable_collision_prevention" => {
                Ok(Box::new(*self.enable_collision_prevention.get_raw_value()))
            }
            "enable_adaptive_curvature" => {
                Ok(Box::new(*self.enable_adaptive_curvature.get_raw_value()))
            }
            _ => Err(ParameterError::not_found(name, "arc")),
        }
    }

    fn set_parameter(
        &mut self,
        name: &str,
        value: Box<dyn std::any::Any>,
        reason: &str,
    ) -> ParameterResult<()> {
        match name {
            "curvature_factor" => {
                if let Some(val) = value.downcast_ref::<f64>() {
                    self.curvature_factor.set_value(*val, reason)
                } else {
                    Err(ParameterError::type_mismatch(name, "f64", "unknown"))
                }
            }
            "smoothness" => {
                if let Some(val) = value.downcast_ref::<usize>() {
                    self.smoothness.set_value(*val, reason)
                } else {
                    Err(ParameterError::type_mismatch(name, "usize", "unknown"))
                }
            }
            "enable_collision_prevention" => {
                if let Some(val) = value.downcast_ref::<bool>() {
                    self.enable_collision_prevention.set_value(*val, reason)
                } else {
                    Err(ParameterError::type_mismatch(name, "bool", "unknown"))
                }
            }
            _ => Err(ParameterError::not_found(name, "arc")),
        }
    }

    fn parameter_names(&self) -> Vec<String> {
        vec![
            "curvature_factor".to_string(),
            "smoothness".to_string(),
            "curvature_direction".to_string(),
            "min_separation_distance".to_string(),
            "max_curvature_reduction".to_string(),
            "enable_collision_prevention".to_string(),
            "enable_adaptive_curvature".to_string(),
        ]
    }

    fn has_parameter(&self, name: &str) -> bool {
        self.parameter_names().contains(&name.to_string())
    }

    fn validate_all(&self) -> ParameterResult<()> {
        self.curvature_factor.validate()?;
        self.smoothness.validate()?;
        self.curvature_direction.validate()?;
        self.min_separation_distance.validate()?;
        self.max_curvature_reduction.validate()?;
        self.enable_collision_prevention.validate()?;
        self.enable_adaptive_curvature.validate()?;
        Ok(())
    }

    fn get_metadata(&self, name: &str) -> ParameterResult<&ParameterMetadata> {
        match name {
            "curvature_factor" => Ok(self.curvature_factor.metadata()),
            "smoothness" => Ok(self.smoothness.metadata()),
            "curvature_direction" => Ok(self.curvature_direction.metadata()),
            "min_separation_distance" => Ok(self.min_separation_distance.metadata()),
            "max_curvature_reduction" => Ok(self.max_curvature_reduction.metadata()),
            "enable_collision_prevention" => Ok(self.enable_collision_prevention.metadata()),
            "enable_adaptive_curvature" => Ok(self.enable_adaptive_curvature.metadata()),
            _ => Err(ParameterError::not_found(name, "arc")),
        }
    }

    fn domain_name(&self) -> &'static str {
        "arc"
    }

    fn reset_all(&mut self, reason: &str) -> ParameterResult<()> {
        self.curvature_factor.reset(reason)?;
        self.smoothness.reset(reason)?;
        self.curvature_direction.reset(reason)?;
        self.min_separation_distance.reset(reason)?;
        self.max_curvature_reduction.reset(reason)?;
        self.enable_collision_prevention.reset(reason)?;
        self.enable_adaptive_curvature.reset(reason)?;
        Ok(())
    }

    fn validation_rules(&self) -> &ValidationRuleSet {
        &self.validation_rules
    }
}
