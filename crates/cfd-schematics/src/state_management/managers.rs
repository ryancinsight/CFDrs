//! Parameter managers for different domains
//!
//! This module provides domain-specific parameter managers that implement
//! the `ParameterManager` trait. Each manager handles parameters for a specific
//! aspect of the microfluidic design system.

use crate::config::{OptimizationProfile, WaveShape};
use crate::state_management::{
    adaptive::{
        BranchCountDensityAdapter, ChannelGenerationContext, DistanceBasedAmplitudeAdapter,
        LengthBasedWavelengthAdapter,
    },
    constraints::ParameterConstraints,
    errors::{ParameterError, ParameterResult},
    parameters::{ConfigurableParameter, ParameterMetadata},
    validation::ValidationRuleSet,
};
use std::collections::HashMap;
use std::fmt::Debug;

/// Core trait for all parameter managers
pub trait ParameterManager: Debug + Send + Sync {
    /// Get parameter value by name
    fn get_parameter(&self, name: &str) -> ParameterResult<Box<dyn std::any::Any>>;

    /// Set parameter value by name
    fn set_parameter(
        &mut self,
        name: &str,
        value: Box<dyn std::any::Any>,
        reason: &str,
    ) -> ParameterResult<()>;

    /// Get all parameter names managed by this manager
    fn parameter_names(&self) -> Vec<String>;

    /// Check if a parameter exists
    fn has_parameter(&self, name: &str) -> bool;

    /// Validate all parameters
    fn validate_all(&self) -> ParameterResult<()>;

    /// Get parameter metadata
    fn get_metadata(&self, name: &str) -> ParameterResult<&ParameterMetadata>;

    /// Get manager domain name
    fn domain_name(&self) -> &str;

    /// Reset all parameters to defaults
    fn reset_all(&mut self, reason: &str) -> ParameterResult<()>;

    /// Get validation rules for this manager
    fn validation_rules(&self) -> &ValidationRuleSet;
}

/// Parameter manager for serpentine channel parameters
#[derive(Debug)]
pub struct SerpentineParameterManager {
    /// Wave amplitude parameter
    amplitude: ConfigurableParameter<f64>,

    /// Wavelength factor parameter
    wavelength_factor: ConfigurableParameter<f64>,

    /// Wave frequency multiplier
    frequency_multiplier: ConfigurableParameter<f64>,

    /// Phase offset for wave generation
    phase_offset: ConfigurableParameter<f64>,

    /// Gaussian width factor for envelope
    gaussian_width_factor: ConfigurableParameter<f64>,

    /// Wave density factor
    wave_density_factor: ConfigurableParameter<f64>,

    /// Fill factor for amplitude scaling
    fill_factor: ConfigurableParameter<f64>,

    /// Wave shape type
    wave_shape: ConfigurableParameter<WaveShape>,

    /// Optimization profile
    optimization_profile: ConfigurableParameter<OptimizationProfile>,

    /// Target fill ratio for optimization
    target_fill_ratio: ConfigurableParameter<f64>,

    /// Validation rules
    validation_rules: ValidationRuleSet,
}

impl SerpentineParameterManager {
    /// Create a new serpentine parameter manager with default values
    #[must_use]
    pub fn new() -> Self {
        // Create amplitude parameter with adaptive behavior
        let amplitude = ConfigurableParameter::new(
            5.0,
            ParameterConstraints::all(vec![
                ParameterConstraints::<f64>::positive(),
                ParameterConstraints::range(0.1, 100.0),
            ]),
            ParameterMetadata::new(
                "amplitude",
                "Base amplitude for serpentine wave generation",
                "wave_parameters",
            )
            .with_units("mm")
            .affects_others(),
        )
        .with_adaptive_behavior(DistanceBasedAmplitudeAdapter::default());

        // Create wavelength factor parameter with adaptive behavior
        let wavelength_factor = ConfigurableParameter::new(
            2.0,
            ParameterConstraints::all(vec![
                ParameterConstraints::<f64>::positive(),
                ParameterConstraints::range(0.1, 10.0),
            ]),
            ParameterMetadata::new(
                "wavelength_factor",
                "Factor for calculating wavelength relative to channel width",
                "wave_parameters",
            )
            .affects_others(),
        )
        .with_adaptive_behavior(LengthBasedWavelengthAdapter::default());

        // Create frequency multiplier parameter
        let frequency_multiplier = ConfigurableParameter::new(
            1.0,
            ParameterConstraints::all(vec![
                ParameterConstraints::<f64>::positive(),
                ParameterConstraints::range(0.1, 5.0),
            ]),
            ParameterMetadata::new(
                "frequency_multiplier",
                "Multiplier for wave frequency calculation",
                "wave_parameters",
            ),
        );

        // Create phase offset parameter
        let phase_offset = ConfigurableParameter::new(
            0.0,
            ParameterConstraints::range(-std::f64::consts::PI, std::f64::consts::PI),
            ParameterMetadata::new(
                "phase_offset",
                "Phase offset for wave generation in radians",
                "wave_parameters",
            )
            .with_units("radians"),
        );

        // Create gaussian width factor parameter
        let gaussian_width_factor = ConfigurableParameter::new(
            0.3,
            ParameterConstraints::all(vec![
                ParameterConstraints::<f64>::positive(),
                ParameterConstraints::range(0.1, 1.0),
            ]),
            ParameterMetadata::new(
                "gaussian_width_factor",
                "Factor for Gaussian envelope width calculation",
                "envelope_parameters",
            ),
        );

        // Create wave density factor parameter with adaptive behavior
        let wave_density_factor = ConfigurableParameter::new(
            2.0,
            ParameterConstraints::all(vec![
                ParameterConstraints::<f64>::positive(),
                ParameterConstraints::range(0.5, 10.0),
            ]),
            ParameterMetadata::new(
                "wave_density_factor",
                "Factor controlling wave density along channel length",
                "wave_parameters",
            )
            .affects_others(),
        )
        .with_adaptive_behavior(BranchCountDensityAdapter::default());

        // Create fill factor parameter
        let fill_factor = ConfigurableParameter::new(
            0.8,
            ParameterConstraints::all(vec![
                ParameterConstraints::<f64>::positive(),
                ParameterConstraints::normalized(),
            ]),
            ParameterMetadata::new(
                "fill_factor",
                "Factor for scaling amplitude relative to available space",
                "scaling_parameters",
            ),
        );

        // Create wave shape parameter
        let wave_shape = ConfigurableParameter::new(
            WaveShape::Sine,
            ParameterConstraints::set(vec![WaveShape::Sine, WaveShape::Square]),
            ParameterMetadata::new(
                "wave_shape",
                "Shape of the wave function (sine or square)",
                "wave_parameters",
            ),
        );

        // Create optimization profile parameter
        let optimization_profile = ConfigurableParameter::new(
            OptimizationProfile::Balanced,
            ParameterConstraints::set(vec![
                OptimizationProfile::Fast,
                OptimizationProfile::Balanced,
                OptimizationProfile::Thorough,
            ]),
            ParameterMetadata::new(
                "optimization_profile",
                "Profile for optimization algorithm behavior",
                "optimization_parameters",
            ),
        );

        // Create target fill ratio parameter
        let target_fill_ratio = ConfigurableParameter::new(
            0.9,
            ParameterConstraints::all(vec![
                ParameterConstraints::<f64>::positive(),
                ParameterConstraints::normalized(),
            ]),
            ParameterMetadata::new(
                "target_fill_ratio",
                "Target fill ratio for optimization algorithms",
                "optimization_parameters",
            ),
        );

        // Create validation rules
        let validation_rules = ValidationRuleSet::new();

        Self {
            amplitude,
            wavelength_factor,
            frequency_multiplier,
            phase_offset,
            gaussian_width_factor,
            wave_density_factor,
            fill_factor,
            wave_shape,
            optimization_profile,
            target_fill_ratio,
            validation_rules,
        }
    }

    /// Get amplitude value with context adaptation
    #[must_use]
    pub fn get_amplitude(&self, context: Option<&ChannelGenerationContext>) -> f64 {
        self.amplitude.get_value(context)
    }

    /// Get wavelength factor value with context adaptation
    #[must_use]
    pub fn get_wavelength_factor(&self, context: Option<&ChannelGenerationContext>) -> f64 {
        self.wavelength_factor.get_value(context)
    }

    /// Get wave density factor value with context adaptation
    #[must_use]
    pub fn get_wave_density_factor(&self, context: Option<&ChannelGenerationContext>) -> f64 {
        self.wave_density_factor.get_value(context)
    }

    /// Set amplitude with validation
    pub fn set_amplitude(&mut self, value: f64, reason: &str) -> ParameterResult<()> {
        self.amplitude.set_value(value, reason)
    }

    /// Set wavelength factor with validation
    pub fn set_wavelength_factor(&mut self, value: f64, reason: &str) -> ParameterResult<()> {
        self.wavelength_factor.set_value(value, reason)
    }

    /// Get all wave parameters as a map for easy access
    #[must_use]
    pub fn get_wave_parameters(
        &self,
        context: Option<&ChannelGenerationContext>,
    ) -> HashMap<String, f64> {
        let mut params = HashMap::new();
        params.insert("amplitude".to_string(), self.get_amplitude(context));
        params.insert(
            "wavelength_factor".to_string(),
            self.get_wavelength_factor(context),
        );
        params.insert(
            "frequency_multiplier".to_string(),
            self.frequency_multiplier.get_value(context),
        );
        params.insert(
            "phase_offset".to_string(),
            self.phase_offset.get_value(context),
        );
        params.insert(
            "gaussian_width_factor".to_string(),
            self.gaussian_width_factor.get_value(context),
        );
        params.insert(
            "wave_density_factor".to_string(),
            self.get_wave_density_factor(context),
        );
        params.insert(
            "fill_factor".to_string(),
            self.fill_factor.get_value(context),
        );
        params.insert(
            "target_fill_ratio".to_string(),
            self.target_fill_ratio.get_value(context),
        );
        params
    }
}

impl Default for SerpentineParameterManager {
    fn default() -> Self {
        Self::new()
    }
}

impl ParameterManager for SerpentineParameterManager {
    fn get_parameter(&self, name: &str) -> ParameterResult<Box<dyn std::any::Any>> {
        match name {
            "amplitude" => Ok(Box::new(*self.amplitude.get_raw_value())),
            "wavelength_factor" => Ok(Box::new(*self.wavelength_factor.get_raw_value())),
            "frequency_multiplier" => Ok(Box::new(*self.frequency_multiplier.get_raw_value())),
            "phase_offset" => Ok(Box::new(*self.phase_offset.get_raw_value())),
            "gaussian_width_factor" => Ok(Box::new(*self.gaussian_width_factor.get_raw_value())),
            "wave_density_factor" => Ok(Box::new(*self.wave_density_factor.get_raw_value())),
            "fill_factor" => Ok(Box::new(*self.fill_factor.get_raw_value())),
            "wave_shape" => Ok(Box::new(*self.wave_shape.get_raw_value())),
            "optimization_profile" => Ok(Box::new(*self.optimization_profile.get_raw_value())),
            "target_fill_ratio" => Ok(Box::new(*self.target_fill_ratio.get_raw_value())),
            _ => Err(ParameterError::not_found(name, "serpentine")),
        }
    }

    fn set_parameter(
        &mut self,
        name: &str,
        value: Box<dyn std::any::Any>,
        reason: &str,
    ) -> ParameterResult<()> {
        match name {
            "amplitude" => {
                if let Some(val) = value.downcast_ref::<f64>() {
                    self.amplitude.set_value(*val, reason)
                } else {
                    Err(ParameterError::type_mismatch(name, "f64", "unknown"))
                }
            }
            "wavelength_factor" => {
                if let Some(val) = value.downcast_ref::<f64>() {
                    self.wavelength_factor.set_value(*val, reason)
                } else {
                    Err(ParameterError::type_mismatch(name, "f64", "unknown"))
                }
            }
            "fill_factor" => {
                if let Some(val) = value.downcast_ref::<f64>() {
                    self.fill_factor.set_value(*val, reason)
                } else {
                    Err(ParameterError::type_mismatch(name, "f64", "unknown"))
                }
            }
            _ => Err(ParameterError::not_found(name, "serpentine")),
        }
    }

    fn parameter_names(&self) -> Vec<String> {
        vec![
            "amplitude".to_string(),
            "wavelength_factor".to_string(),
            "frequency_multiplier".to_string(),
            "phase_offset".to_string(),
            "gaussian_width_factor".to_string(),
            "wave_density_factor".to_string(),
            "fill_factor".to_string(),
            "wave_shape".to_string(),
            "optimization_profile".to_string(),
            "target_fill_ratio".to_string(),
        ]
    }

    fn has_parameter(&self, name: &str) -> bool {
        self.parameter_names().contains(&name.to_string())
    }

    fn validate_all(&self) -> ParameterResult<()> {
        self.amplitude.validate()?;
        self.wavelength_factor.validate()?;
        self.frequency_multiplier.validate()?;
        self.phase_offset.validate()?;
        self.gaussian_width_factor.validate()?;
        self.wave_density_factor.validate()?;
        self.fill_factor.validate()?;
        self.wave_shape.validate()?;
        self.optimization_profile.validate()?;
        self.target_fill_ratio.validate()?;
        Ok(())
    }

    fn get_metadata(&self, name: &str) -> ParameterResult<&ParameterMetadata> {
        match name {
            "amplitude" => Ok(self.amplitude.metadata()),
            "wavelength_factor" => Ok(self.wavelength_factor.metadata()),
            "frequency_multiplier" => Ok(self.frequency_multiplier.metadata()),
            "phase_offset" => Ok(self.phase_offset.metadata()),
            "gaussian_width_factor" => Ok(self.gaussian_width_factor.metadata()),
            "wave_density_factor" => Ok(self.wave_density_factor.metadata()),
            "fill_factor" => Ok(self.fill_factor.metadata()),
            "wave_shape" => Ok(self.wave_shape.metadata()),
            "optimization_profile" => Ok(self.optimization_profile.metadata()),
            "target_fill_ratio" => Ok(self.target_fill_ratio.metadata()),
            _ => Err(ParameterError::not_found(name, "serpentine")),
        }
    }

    fn domain_name(&self) -> &'static str {
        "serpentine"
    }

    fn reset_all(&mut self, reason: &str) -> ParameterResult<()> {
        self.amplitude.reset(reason)?;
        self.wavelength_factor.reset(reason)?;
        self.frequency_multiplier.reset(reason)?;
        self.phase_offset.reset(reason)?;
        self.gaussian_width_factor.reset(reason)?;
        self.wave_density_factor.reset(reason)?;
        self.fill_factor.reset(reason)?;
        self.wave_shape.reset(reason)?;
        self.optimization_profile.reset(reason)?;
        self.target_fill_ratio.reset(reason)?;
        Ok(())
    }

    fn validation_rules(&self) -> &ValidationRuleSet {
        &self.validation_rules
    }
}

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

/// Parameter manager for geometry parameters
#[derive(Debug)]
pub struct GeometryParameterManager {
    /// Wall clearance parameter
    wall_clearance: ConfigurableParameter<f64>,

    /// Channel width parameter
    channel_width: ConfigurableParameter<f64>,

    /// Channel height parameter
    channel_height: ConfigurableParameter<f64>,

    /// Validation rules
    validation_rules: ValidationRuleSet,
}

impl GeometryParameterManager {
    /// Create a new geometry parameter manager
    #[must_use]
    pub fn new() -> Self {
        let wall_clearance = ConfigurableParameter::new(
            0.5,
            ParameterConstraints::all(vec![
                ParameterConstraints::<f64>::positive(),
                ParameterConstraints::range(0.01, 10.0),
            ]),
            ParameterMetadata::new(
                "wall_clearance",
                "Minimum distance between channels and walls",
                "geometry_parameters",
            )
            .with_units("mm")
            .affects_others(),
        );

        let channel_width = ConfigurableParameter::new(
            1.0,
            ParameterConstraints::all(vec![
                ParameterConstraints::<f64>::positive(),
                ParameterConstraints::range(0.01, 50.0),
            ]),
            ParameterMetadata::new(
                "channel_width",
                "Width of microfluidic channels",
                "geometry_parameters",
            )
            .with_units("mm")
            .affects_others(),
        );

        let channel_height = ConfigurableParameter::new(
            0.5,
            ParameterConstraints::all(vec![
                ParameterConstraints::<f64>::positive(),
                ParameterConstraints::range(0.01, 25.0),
            ]),
            ParameterMetadata::new(
                "channel_height",
                "Height of microfluidic channels",
                "geometry_parameters",
            )
            .with_units("mm"),
        );

        let validation_rules = ValidationRuleSet::new();

        Self {
            wall_clearance,
            channel_width,
            channel_height,
            validation_rules,
        }
    }
}

impl Default for GeometryParameterManager {
    fn default() -> Self {
        Self::new()
    }
}

impl ParameterManager for GeometryParameterManager {
    fn get_parameter(&self, name: &str) -> ParameterResult<Box<dyn std::any::Any>> {
        match name {
            "wall_clearance" => Ok(Box::new(*self.wall_clearance.get_raw_value())),
            "channel_width" => Ok(Box::new(*self.channel_width.get_raw_value())),
            "channel_height" => Ok(Box::new(*self.channel_height.get_raw_value())),
            _ => Err(ParameterError::not_found(name, "geometry")),
        }
    }

    fn set_parameter(
        &mut self,
        name: &str,
        value: Box<dyn std::any::Any>,
        reason: &str,
    ) -> ParameterResult<()> {
        match name {
            "wall_clearance" => {
                if let Some(val) = value.downcast_ref::<f64>() {
                    self.wall_clearance.set_value(*val, reason)
                } else {
                    Err(ParameterError::type_mismatch(name, "f64", "unknown"))
                }
            }
            "channel_width" => {
                if let Some(val) = value.downcast_ref::<f64>() {
                    self.channel_width.set_value(*val, reason)
                } else {
                    Err(ParameterError::type_mismatch(name, "f64", "unknown"))
                }
            }
            "channel_height" => {
                if let Some(val) = value.downcast_ref::<f64>() {
                    self.channel_height.set_value(*val, reason)
                } else {
                    Err(ParameterError::type_mismatch(name, "f64", "unknown"))
                }
            }
            _ => Err(ParameterError::not_found(name, "geometry")),
        }
    }

    fn parameter_names(&self) -> Vec<String> {
        vec![
            "wall_clearance".to_string(),
            "channel_width".to_string(),
            "channel_height".to_string(),
        ]
    }

    fn has_parameter(&self, name: &str) -> bool {
        self.parameter_names().contains(&name.to_string())
    }

    fn validate_all(&self) -> ParameterResult<()> {
        self.wall_clearance.validate()?;
        self.channel_width.validate()?;
        self.channel_height.validate()?;
        Ok(())
    }

    fn get_metadata(&self, name: &str) -> ParameterResult<&ParameterMetadata> {
        match name {
            "wall_clearance" => Ok(self.wall_clearance.metadata()),
            "channel_width" => Ok(self.channel_width.metadata()),
            "channel_height" => Ok(self.channel_height.metadata()),
            _ => Err(ParameterError::not_found(name, "geometry")),
        }
    }

    fn domain_name(&self) -> &'static str {
        "geometry"
    }

    fn reset_all(&mut self, reason: &str) -> ParameterResult<()> {
        self.wall_clearance.reset(reason)?;
        self.channel_width.reset(reason)?;
        self.channel_height.reset(reason)?;
        Ok(())
    }

    fn validation_rules(&self) -> &ValidationRuleSet {
        &self.validation_rules
    }
}

/// Stub parameter managers for collision and symmetry (to be implemented later)
#[derive(Debug)]
pub struct CollisionParameterManager {
    validation_rules: ValidationRuleSet,
}

impl CollisionParameterManager {
    #[must_use]
    pub fn new() -> Self {
        Self {
            validation_rules: ValidationRuleSet::new(),
        }
    }
}

impl Default for CollisionParameterManager {
    fn default() -> Self {
        Self::new()
    }
}

impl ParameterManager for CollisionParameterManager {
    fn get_parameter(&self, _name: &str) -> ParameterResult<Box<dyn std::any::Any>> {
        Err(ParameterError::not_found("any", "collision"))
    }

    fn set_parameter(
        &mut self,
        name: &str,
        _value: Box<dyn std::any::Any>,
        _reason: &str,
    ) -> ParameterResult<()> {
        Err(ParameterError::not_found(name, "collision"))
    }

    fn parameter_names(&self) -> Vec<String> {
        vec![]
    }

    fn has_parameter(&self, _name: &str) -> bool {
        false
    }

    fn validate_all(&self) -> ParameterResult<()> {
        Ok(())
    }

    fn get_metadata(&self, name: &str) -> ParameterResult<&ParameterMetadata> {
        Err(ParameterError::not_found(name, "collision"))
    }

    fn domain_name(&self) -> &'static str {
        "collision"
    }

    fn reset_all(&mut self, _reason: &str) -> ParameterResult<()> {
        Ok(())
    }

    fn validation_rules(&self) -> &ValidationRuleSet {
        &self.validation_rules
    }
}

#[derive(Debug)]
pub struct SymmetryParameterManager {
    validation_rules: ValidationRuleSet,
}

impl SymmetryParameterManager {
    #[must_use]
    pub fn new() -> Self {
        Self {
            validation_rules: ValidationRuleSet::new(),
        }
    }
}

impl Default for SymmetryParameterManager {
    fn default() -> Self {
        Self::new()
    }
}

impl ParameterManager for SymmetryParameterManager {
    fn get_parameter(&self, _name: &str) -> ParameterResult<Box<dyn std::any::Any>> {
        Err(ParameterError::not_found("any", "symmetry"))
    }

    fn set_parameter(
        &mut self,
        name: &str,
        _value: Box<dyn std::any::Any>,
        _reason: &str,
    ) -> ParameterResult<()> {
        Err(ParameterError::not_found(name, "symmetry"))
    }

    fn parameter_names(&self) -> Vec<String> {
        vec![]
    }

    fn has_parameter(&self, _name: &str) -> bool {
        false
    }

    fn validate_all(&self) -> ParameterResult<()> {
        Ok(())
    }

    fn get_metadata(&self, name: &str) -> ParameterResult<&ParameterMetadata> {
        Err(ParameterError::not_found(name, "symmetry"))
    }

    fn domain_name(&self) -> &'static str {
        "symmetry"
    }

    fn reset_all(&mut self, _reason: &str) -> ParameterResult<()> {
        Ok(())
    }

    fn validation_rules(&self) -> &ValidationRuleSet {
        &self.validation_rules
    }
}
