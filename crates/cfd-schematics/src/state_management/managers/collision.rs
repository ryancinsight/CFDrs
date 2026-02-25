//! Collision detection parameter manager.
//!
//! Manages parameters that control channel-to-channel and channel-to-wall
//! collision detection and avoidance behaviour. Parameters mirror the
//! [`CollisionParameters`](crate::geometry::collision_detection::CollisionParameters)
//! struct used at runtime.

use crate::state_management::{
    constraints::ParameterConstraints,
    errors::{ParameterError, ParameterResult},
    parameters::{ConfigurableParameter, ParameterMetadata},
    validation::ValidationRuleSet,
};

use super::ParameterManager;

/// Parameter manager for collision detection and avoidance
///
/// Manages the following parameters:
///
/// | Name | Type | Default | Range | Units |
/// |------|------|---------|-------|-------|
/// | `min_channel_distance` | f64 | 0.5 | 0.01–50.0 | mm |
/// | `min_wall_distance` | f64 | 0.5 | 0.01–50.0 | mm |
/// | `safety_margin_factor` | f64 | 1.1 | 1.0–5.0 | — |
/// | `max_reduction_factor` | f64 | 0.8 | 0.01–1.0 | — |
/// | `detection_sensitivity` | f64 | 0.1 | 0.001–1.0 | — |
/// | `enable_neighbor_detection` | bool | true | — | — |
/// | `enable_wall_detection` | bool | true | — | — |
#[derive(Debug)]
pub struct CollisionParameterManager {
    /// Minimum distance between channels
    min_channel_distance: ConfigurableParameter<f64>,

    /// Minimum distance from walls
    min_wall_distance: ConfigurableParameter<f64>,

    /// Safety margin multiplier applied to clearance checks
    safety_margin_factor: ConfigurableParameter<f64>,

    /// Maximum curvature/amplitude reduction allowed to avoid collisions
    max_reduction_factor: ConfigurableParameter<f64>,

    /// Sensitivity of the collision detection sweep
    detection_sensitivity: ConfigurableParameter<f64>,

    /// Enable neighbour-based collision detection
    enable_neighbor_detection: ConfigurableParameter<bool>,

    /// Enable wall-based collision detection
    enable_wall_detection: ConfigurableParameter<bool>,

    /// Validation rules
    validation_rules: ValidationRuleSet,
}

impl CollisionParameterManager {
    /// Create a new collision parameter manager with defaults from the constants registry
    #[must_use]
    pub fn new() -> Self {
        let min_channel_distance = ConfigurableParameter::new(
            0.5,
            ParameterConstraints::all(vec![
                ParameterConstraints::<f64>::positive(),
                ParameterConstraints::range(0.01, 50.0),
            ]),
            ParameterMetadata::new(
                "min_channel_distance",
                "Minimum centre-to-centre distance between channels",
                "collision_parameters",
            )
            .with_units("mm")
            .affects_others(),
        );

        let min_wall_distance = ConfigurableParameter::new(
            0.5,
            ParameterConstraints::all(vec![
                ParameterConstraints::<f64>::positive(),
                ParameterConstraints::range(0.01, 50.0),
            ]),
            ParameterMetadata::new(
                "min_wall_distance",
                "Minimum distance from channels to wall boundaries",
                "collision_parameters",
            )
            .with_units("mm")
            .affects_others(),
        );

        let safety_margin_factor = ConfigurableParameter::new(
            1.1,
            ParameterConstraints::range(1.0, 5.0),
            ParameterMetadata::new(
                "safety_margin_factor",
                "Multiplicative safety margin applied to clearance checks",
                "collision_parameters",
            ),
        );

        let max_reduction_factor = ConfigurableParameter::new(
            0.8,
            ParameterConstraints::all(vec![
                ParameterConstraints::<f64>::positive(),
                ParameterConstraints::range(0.01, 1.0),
            ]),
            ParameterMetadata::new(
                "max_reduction_factor",
                "Maximum curvature/amplitude reduction factor for collision avoidance",
                "collision_parameters",
            ),
        );

        let detection_sensitivity = ConfigurableParameter::new(
            0.1,
            ParameterConstraints::all(vec![
                ParameterConstraints::<f64>::positive(),
                ParameterConstraints::range(0.001, 1.0),
            ]),
            ParameterMetadata::new(
                "detection_sensitivity",
                "Collision detection sweep sensitivity",
                "collision_parameters",
            ),
        );

        let enable_neighbor_detection = ConfigurableParameter::new(
            true,
            ParameterConstraints::none(),
            ParameterMetadata::new(
                "enable_neighbor_detection",
                "Enable channel-to-channel collision detection",
                "collision_parameters",
            ),
        );

        let enable_wall_detection = ConfigurableParameter::new(
            true,
            ParameterConstraints::none(),
            ParameterMetadata::new(
                "enable_wall_detection",
                "Enable channel-to-wall collision detection",
                "collision_parameters",
            ),
        );

        let validation_rules = ValidationRuleSet::new();

        Self {
            min_channel_distance,
            min_wall_distance,
            safety_margin_factor,
            max_reduction_factor,
            detection_sensitivity,
            enable_neighbor_detection,
            enable_wall_detection,
            validation_rules,
        }
    }

    /// Get the minimum channel distance
    #[must_use]
    pub const fn get_min_channel_distance(&self) -> f64 {
        *self.min_channel_distance.get_raw_value()
    }

    /// Get the minimum wall distance
    #[must_use]
    pub const fn get_min_wall_distance(&self) -> f64 {
        *self.min_wall_distance.get_raw_value()
    }

    /// Get the safety margin factor
    #[must_use]
    pub const fn get_safety_margin_factor(&self) -> f64 {
        *self.safety_margin_factor.get_raw_value()
    }

    /// Check if neighbour detection is enabled
    #[must_use]
    pub const fn is_neighbor_detection_enabled(&self) -> bool {
        *self.enable_neighbor_detection.get_raw_value()
    }

    /// Check if wall detection is enabled
    #[must_use]
    pub const fn is_wall_detection_enabled(&self) -> bool {
        *self.enable_wall_detection.get_raw_value()
    }
}

impl Default for CollisionParameterManager {
    fn default() -> Self {
        Self::new()
    }
}

const COLLISION_PARAMS: &[&str] = &[
    "min_channel_distance",
    "min_wall_distance",
    "safety_margin_factor",
    "max_reduction_factor",
    "detection_sensitivity",
    "enable_neighbor_detection",
    "enable_wall_detection",
];

impl ParameterManager for CollisionParameterManager {
    fn get_parameter(&self, name: &str) -> ParameterResult<Box<dyn std::any::Any>> {
        match name {
            "min_channel_distance" => Ok(Box::new(*self.min_channel_distance.get_raw_value())),
            "min_wall_distance" => Ok(Box::new(*self.min_wall_distance.get_raw_value())),
            "safety_margin_factor" => Ok(Box::new(*self.safety_margin_factor.get_raw_value())),
            "max_reduction_factor" => Ok(Box::new(*self.max_reduction_factor.get_raw_value())),
            "detection_sensitivity" => Ok(Box::new(*self.detection_sensitivity.get_raw_value())),
            "enable_neighbor_detection" => {
                Ok(Box::new(*self.enable_neighbor_detection.get_raw_value()))
            }
            "enable_wall_detection" => Ok(Box::new(*self.enable_wall_detection.get_raw_value())),
            _ => Err(ParameterError::not_found(name, "collision")),
        }
    }

    fn set_parameter(
        &mut self,
        name: &str,
        value: Box<dyn std::any::Any>,
        reason: &str,
    ) -> ParameterResult<()> {
        match name {
            "min_channel_distance" => {
                if let Some(val) = value.downcast_ref::<f64>() {
                    self.min_channel_distance.set_value(*val, reason)
                } else {
                    Err(ParameterError::type_mismatch(name, "f64", "unknown"))
                }
            }
            "min_wall_distance" => {
                if let Some(val) = value.downcast_ref::<f64>() {
                    self.min_wall_distance.set_value(*val, reason)
                } else {
                    Err(ParameterError::type_mismatch(name, "f64", "unknown"))
                }
            }
            "safety_margin_factor" => {
                if let Some(val) = value.downcast_ref::<f64>() {
                    self.safety_margin_factor.set_value(*val, reason)
                } else {
                    Err(ParameterError::type_mismatch(name, "f64", "unknown"))
                }
            }
            "max_reduction_factor" => {
                if let Some(val) = value.downcast_ref::<f64>() {
                    self.max_reduction_factor.set_value(*val, reason)
                } else {
                    Err(ParameterError::type_mismatch(name, "f64", "unknown"))
                }
            }
            "detection_sensitivity" => {
                if let Some(val) = value.downcast_ref::<f64>() {
                    self.detection_sensitivity.set_value(*val, reason)
                } else {
                    Err(ParameterError::type_mismatch(name, "f64", "unknown"))
                }
            }
            "enable_neighbor_detection" => {
                if let Some(val) = value.downcast_ref::<bool>() {
                    self.enable_neighbor_detection.set_value(*val, reason)
                } else {
                    Err(ParameterError::type_mismatch(name, "bool", "unknown"))
                }
            }
            "enable_wall_detection" => {
                if let Some(val) = value.downcast_ref::<bool>() {
                    self.enable_wall_detection.set_value(*val, reason)
                } else {
                    Err(ParameterError::type_mismatch(name, "bool", "unknown"))
                }
            }
            _ => Err(ParameterError::not_found(name, "collision")),
        }
    }

    fn parameter_names(&self) -> Vec<String> {
        COLLISION_PARAMS.iter().map(|s| (*s).to_string()).collect()
    }

    fn has_parameter(&self, name: &str) -> bool {
        COLLISION_PARAMS.contains(&name)
    }

    fn validate_all(&self) -> ParameterResult<()> {
        self.min_channel_distance.validate()?;
        self.min_wall_distance.validate()?;
        self.safety_margin_factor.validate()?;
        self.max_reduction_factor.validate()?;
        self.detection_sensitivity.validate()?;
        self.enable_neighbor_detection.validate()?;
        self.enable_wall_detection.validate()?;
        Ok(())
    }

    fn get_metadata(&self, name: &str) -> ParameterResult<&ParameterMetadata> {
        match name {
            "min_channel_distance" => Ok(self.min_channel_distance.metadata()),
            "min_wall_distance" => Ok(self.min_wall_distance.metadata()),
            "safety_margin_factor" => Ok(self.safety_margin_factor.metadata()),
            "max_reduction_factor" => Ok(self.max_reduction_factor.metadata()),
            "detection_sensitivity" => Ok(self.detection_sensitivity.metadata()),
            "enable_neighbor_detection" => Ok(self.enable_neighbor_detection.metadata()),
            "enable_wall_detection" => Ok(self.enable_wall_detection.metadata()),
            _ => Err(ParameterError::not_found(name, "collision")),
        }
    }

    fn domain_name(&self) -> &'static str {
        "collision"
    }

    fn reset_all(&mut self, reason: &str) -> ParameterResult<()> {
        self.min_channel_distance.reset(reason)?;
        self.min_wall_distance.reset(reason)?;
        self.safety_margin_factor.reset(reason)?;
        self.max_reduction_factor.reset(reason)?;
        self.detection_sensitivity.reset(reason)?;
        self.enable_neighbor_detection.reset(reason)?;
        self.enable_wall_detection.reset(reason)?;
        Ok(())
    }

    fn validation_rules(&self) -> &ValidationRuleSet {
        &self.validation_rules
    }
}
