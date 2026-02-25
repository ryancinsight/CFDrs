//! Serpentine channel parameter manager.

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

use super::ParameterManager;

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

        let wave_shape = ConfigurableParameter::new(
            WaveShape::Sine,
            ParameterConstraints::set(vec![WaveShape::Sine, WaveShape::Square]),
            ParameterMetadata::new(
                "wave_shape",
                "Shape of the wave function (sine or square)",
                "wave_parameters",
            ),
        );

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
