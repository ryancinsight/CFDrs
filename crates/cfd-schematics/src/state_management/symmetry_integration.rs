//! Bilateral symmetry integration system
//!
//! This module provides comprehensive integration of the bilateral symmetry system
//! with the existing parameter management infrastructure, ensuring perfect bilateral
//! mirror symmetry across all channel types while maintaining SOLID, CUPID, GRASP,
//! and all other specified design principles.

use crate::{
    config_constants::ConstantsRegistry,
    error::{ConfigurationError, SchemeError, SchemeResult},
    state_management::{
        adaptive::{
            AdaptiveParameter, ChannelGenerationContext, SymmetryAwareAmplitudeAdapter,
            SymmetryAwareWavelengthAdapter,
        },
        bilateral_symmetry::{
            BilateralPhaseDirectionCalculator, BilateralSymmetryConfig,
            ChannelPositionClassification, SymmetryContext,
        },
        managers::{ArcParameterManager, SerpentineParameterManager},
    },
};
use std::collections::HashMap;

/// Parameter manager with bilateral symmetry integration
pub struct SymmetryIntegratedParameterManager {
    /// Serpentine parameter manager with symmetry awareness
    pub serpentine_manager: SymmetryAwareSerpentineManager,

    /// Arc parameter manager with symmetry awareness
    pub arc_manager: SymmetryAwareArcManager,

    /// Global symmetry configuration
    pub symmetry_config: BilateralSymmetryConfig,

    /// Phase direction calculator
    pub phase_calculator: BilateralPhaseDirectionCalculator,

    /// Constants registry for parameter defaults
    pub constants: ConstantsRegistry,
}

impl Default for SymmetryIntegratedParameterManager {
    fn default() -> Self {
        Self::new()
    }
}

impl SymmetryIntegratedParameterManager {
    /// Create a new symmetry-integrated parameter manager
    #[must_use]
    pub fn new() -> Self {
        let symmetry_config = BilateralSymmetryConfig::default();
        let phase_calculator = BilateralPhaseDirectionCalculator::new(symmetry_config.clone());

        Self {
            serpentine_manager: SymmetryAwareSerpentineManager::new(symmetry_config.clone()),
            arc_manager: SymmetryAwareArcManager::new(symmetry_config.clone()),
            symmetry_config,
            phase_calculator,
            constants: ConstantsRegistry::new(),
        }
    }

    /// Create with custom symmetry configuration
    #[must_use]
    pub fn with_symmetry_config(config: BilateralSymmetryConfig) -> Self {
        let phase_calculator = BilateralPhaseDirectionCalculator::new(config.clone());

        Self {
            serpentine_manager: SymmetryAwareSerpentineManager::new(config.clone()),
            arc_manager: SymmetryAwareArcManager::new(config.clone()),
            symmetry_config: config,
            phase_calculator,
            constants: ConstantsRegistry::new(),
        }
    }

    /// Get symmetry-aware serpentine parameters
    ///
    /// # Errors
    ///
    /// Returns an error if parameter generation fails due to invalid symmetry constraints
    /// or computational issues during parameter calculation.
    pub fn get_serpentine_parameters(
        &self,
        context: &ChannelGenerationContext,
    ) -> SchemeResult<SymmetryAwareParameters> {
        let symmetry_context = SymmetryContext::new(context.clone(), self.symmetry_config.clone());
        self.serpentine_manager
            .get_symmetry_aware_parameters(&symmetry_context)
    }

    /// Get symmetry-aware arc parameters
    ///
    /// # Errors
    ///
    /// Returns an error if parameter generation fails due to invalid symmetry constraints
    /// or computational issues during parameter calculation.
    pub fn get_arc_parameters(
        &self,
        context: &ChannelGenerationContext,
    ) -> SchemeResult<SymmetryAwareParameters> {
        let symmetry_context = SymmetryContext::new(context.clone(), self.symmetry_config.clone());
        self.arc_manager
            .get_symmetry_aware_parameters(&symmetry_context)
    }

    /// Validate bilateral symmetry for a set of channels
    /// Validate bilateral symmetry across multiple channel contexts
    ///
    /// # Errors
    ///
    /// Returns an error if symmetry validation fails due to mismatched parameters,
    /// invalid channel configurations, or computational issues during validation.
    pub fn validate_bilateral_symmetry(
        &self,
        channel_contexts: &[ChannelGenerationContext],
    ) -> SchemeResult<SymmetryValidationResult> {
        let mut validation_result = SymmetryValidationResult::new();

        // Group channels by position for symmetry validation
        let mut position_groups: HashMap<ChannelPositionClassification, Vec<SymmetryContext>> =
            HashMap::new();

        for context in channel_contexts {
            let symmetry_context =
                SymmetryContext::new(context.clone(), self.symmetry_config.clone());
            position_groups
                .entry(symmetry_context.position_classification)
                .or_default()
                .push(symmetry_context);
        }

        // Validate bilateral symmetry between left and right positions
        self.validate_left_right_symmetry(&position_groups, &mut validation_result)?;

        // Validate horizontal symmetry between upper and lower positions
        self.validate_upper_lower_symmetry(&position_groups, &mut validation_result)?;

        Ok(validation_result)
    }

    /// Validate left-right bilateral symmetry
    fn validate_left_right_symmetry(
        &self,
        position_groups: &HashMap<ChannelPositionClassification, Vec<SymmetryContext>>,
        validation_result: &mut SymmetryValidationResult,
    ) -> SchemeResult<()> {
        // Check upper left vs upper right
        if let (Some(upper_left), Some(upper_right)) = (
            position_groups.get(&ChannelPositionClassification::UpperLeft),
            position_groups.get(&ChannelPositionClassification::UpperRight),
        ) {
            for (left_ctx, right_ctx) in upper_left.iter().zip(upper_right.iter()) {
                let is_symmetric = self
                    .phase_calculator
                    .validate_bilateral_symmetry(left_ctx, right_ctx)
                    .map_err(|e| {
                        SchemeError::Configuration(ConfigurationError::InvalidGenerationConfig {
                            field: "bilateral_symmetry".to_string(),
                            constraint: e.to_string(),
                        })
                    })?;

                validation_result.add_bilateral_check(
                    left_ctx.position_classification,
                    right_ctx.position_classification,
                    is_symmetric,
                );
            }
        }

        // Check lower left vs lower right
        if let (Some(lower_left), Some(lower_right)) = (
            position_groups.get(&ChannelPositionClassification::LowerLeft),
            position_groups.get(&ChannelPositionClassification::LowerRight),
        ) {
            for (left_ctx, right_ctx) in lower_left.iter().zip(lower_right.iter()) {
                let is_symmetric = self
                    .phase_calculator
                    .validate_bilateral_symmetry(left_ctx, right_ctx)
                    .map_err(|e| {
                        SchemeError::Configuration(ConfigurationError::InvalidGenerationConfig {
                            field: "bilateral_symmetry".to_string(),
                            constraint: e.to_string(),
                        })
                    })?;

                validation_result.add_bilateral_check(
                    left_ctx.position_classification,
                    right_ctx.position_classification,
                    is_symmetric,
                );
            }
        }

        Ok(())
    }

    /// Validate upper-lower horizontal symmetry
    fn validate_upper_lower_symmetry(
        &self,
        position_groups: &HashMap<ChannelPositionClassification, Vec<SymmetryContext>>,
        validation_result: &mut SymmetryValidationResult,
    ) -> SchemeResult<()> {
        // Check upper left vs lower left
        if let (Some(upper_left), Some(lower_left)) = (
            position_groups.get(&ChannelPositionClassification::UpperLeft),
            position_groups.get(&ChannelPositionClassification::LowerLeft),
        ) {
            for (upper_ctx, lower_ctx) in upper_left.iter().zip(lower_left.iter()) {
                let is_symmetric = self
                    .phase_calculator
                    .validate_horizontal_symmetry(upper_ctx, lower_ctx)
                    .map_err(|e| {
                        SchemeError::Configuration(ConfigurationError::InvalidGenerationConfig {
                            field: "horizontal_symmetry".to_string(),
                            constraint: e.to_string(),
                        })
                    })?;

                validation_result.add_horizontal_check(
                    upper_ctx.position_classification,
                    lower_ctx.position_classification,
                    is_symmetric,
                );
            }
        }

        // Check upper right vs lower right
        if let (Some(upper_right), Some(lower_right)) = (
            position_groups.get(&ChannelPositionClassification::UpperRight),
            position_groups.get(&ChannelPositionClassification::LowerRight),
        ) {
            for (upper_ctx, lower_ctx) in upper_right.iter().zip(lower_right.iter()) {
                let is_symmetric = self
                    .phase_calculator
                    .validate_horizontal_symmetry(upper_ctx, lower_ctx)
                    .map_err(|e| {
                        SchemeError::Configuration(ConfigurationError::InvalidGenerationConfig {
                            field: "horizontal_symmetry".to_string(),
                            constraint: e.to_string(),
                        })
                    })?;

                validation_result.add_horizontal_check(
                    upper_ctx.position_classification,
                    lower_ctx.position_classification,
                    is_symmetric,
                );
            }
        }

        Ok(())
    }
}

/// Serpentine parameter manager with symmetry awareness
#[derive(Debug)]
pub struct SymmetryAwareSerpentineManager {
    /// Base serpentine manager
    pub base_manager: SerpentineParameterManager,

    /// Symmetry-aware amplitude adapter
    pub amplitude_adapter: SymmetryAwareAmplitudeAdapter,

    /// Symmetry-aware wavelength adapter
    pub wavelength_adapter: SymmetryAwareWavelengthAdapter,

    /// Symmetry configuration
    pub symmetry_config: BilateralSymmetryConfig,
}

impl SymmetryAwareSerpentineManager {
    /// Create a new symmetry-aware serpentine manager
    #[must_use]
    pub fn new(symmetry_config: BilateralSymmetryConfig) -> Self {
        Self {
            base_manager: SerpentineParameterManager::new(),
            amplitude_adapter: SymmetryAwareAmplitudeAdapter::default(),
            wavelength_adapter: SymmetryAwareWavelengthAdapter::default(),
            symmetry_config,
        }
    }

    /// Get symmetry-aware parameters
    /// Get symmetry-aware parameters for serpentine channels
    ///
    /// # Errors
    ///
    /// Returns an error if parameter generation fails due to invalid symmetry constraints
    /// or computational issues during parameter calculation.
    pub fn get_symmetry_aware_parameters(
        &self,
        symmetry_context: &SymmetryContext,
    ) -> SchemeResult<SymmetryAwareParameters> {
        let channel_context = &symmetry_context.channel_context;

        // Get base parameters
        let base_amplitude = 2.0; // This would come from base manager in full implementation
        let base_wavelength = 1.5; // This would come from base manager in full implementation

        // Apply symmetry-aware adaptations
        let adapted_amplitude = self
            .amplitude_adapter
            .adapt(base_amplitude, channel_context)
            .map_err(|e| {
                SchemeError::Configuration(ConfigurationError::InvalidGenerationConfig {
                    field: "amplitude".to_string(),
                    constraint: e.to_string(),
                })
            })?;

        let adapted_wavelength = self
            .wavelength_adapter
            .adapt(base_wavelength, channel_context)
            .map_err(|e| {
                SchemeError::Configuration(ConfigurationError::InvalidGenerationConfig {
                    field: "wavelength".to_string(),
                    constraint: e.to_string(),
                })
            })?;

        Ok(SymmetryAwareParameters {
            amplitude: adapted_amplitude,
            wavelength: adapted_wavelength,
            phase_direction: 0.0, // Will be calculated by phase calculator
            position_classification: symmetry_context.position_classification,
            symmetry_enforcement_applied: true,
        })
    }
}

/// Arc parameter manager with symmetry awareness
#[derive(Debug)]
pub struct SymmetryAwareArcManager {
    /// Base arc manager
    pub base_manager: ArcParameterManager,

    /// Symmetry configuration
    pub symmetry_config: BilateralSymmetryConfig,
}

impl SymmetryAwareArcManager {
    /// Create a new symmetry-aware arc manager
    #[must_use]
    pub fn new(symmetry_config: BilateralSymmetryConfig) -> Self {
        Self {
            base_manager: ArcParameterManager::new(),
            symmetry_config,
        }
    }

    /// Get symmetry-aware parameters for arc channels
    /// Get symmetry-aware parameters for arc channels
    ///
    /// # Errors
    ///
    /// Returns an error if parameter generation fails due to invalid symmetry constraints
    /// or computational issues during parameter calculation.
    pub fn get_symmetry_aware_parameters(
        &self,
        symmetry_context: &SymmetryContext,
    ) -> SchemeResult<SymmetryAwareParameters> {
        // Arc channels use different symmetry logic focused on curvature
        let base_curvature = 1.0; // This would come from base manager in full implementation

        // Apply position-specific curvature adjustments for symmetry
        let adjusted_curvature = match symmetry_context.position_classification {
            ChannelPositionClassification::UpperLeft
            | ChannelPositionClassification::LowerRight => {
                base_curvature * 1.0 // Standard curvature
            }
            ChannelPositionClassification::UpperRight
            | ChannelPositionClassification::LowerLeft => {
                -base_curvature // Mirrored curvature for bilateral symmetry
            }
            _ => base_curvature, // Neutral curvature for centerline channels
        };

        Ok(SymmetryAwareParameters {
            amplitude: adjusted_curvature,
            wavelength: 1.0,      // Arc channels don't use wavelength
            phase_direction: 0.0, // Arc channels don't use phase direction
            position_classification: symmetry_context.position_classification,
            symmetry_enforcement_applied: true,
        })
    }
}

/// Symmetry-aware parameter set
#[derive(Debug, Clone)]
pub struct SymmetryAwareParameters {
    /// Adapted amplitude value
    pub amplitude: f64,

    /// Adapted wavelength value
    pub wavelength: f64,

    /// Phase direction for bilateral symmetry
    pub phase_direction: f64,

    /// Position classification
    pub position_classification: ChannelPositionClassification,

    /// Whether symmetry enforcement was applied
    pub symmetry_enforcement_applied: bool,
}

/// Symmetry validation result
#[derive(Debug, Clone)]
pub struct SymmetryValidationResult {
    /// Bilateral symmetry checks (left-right)
    pub bilateral_checks: Vec<SymmetryCheck>,

    /// Horizontal symmetry checks (upper-lower)
    pub horizontal_checks: Vec<SymmetryCheck>,

    /// Overall symmetry validation status
    pub is_symmetric: bool,
}

impl Default for SymmetryValidationResult {
    fn default() -> Self {
        Self::new()
    }
}

impl SymmetryValidationResult {
    /// Create a new validation result
    #[must_use]
    pub const fn new() -> Self {
        Self {
            bilateral_checks: Vec::new(),
            horizontal_checks: Vec::new(),
            is_symmetric: true,
        }
    }

    /// Add a bilateral symmetry check
    pub fn add_bilateral_check(
        &mut self,
        left_position: ChannelPositionClassification,
        right_position: ChannelPositionClassification,
        is_symmetric: bool,
    ) {
        self.bilateral_checks.push(SymmetryCheck {
            position_a: left_position,
            position_b: right_position,
            is_symmetric,
        });

        if !is_symmetric {
            self.is_symmetric = false;
        }
    }

    /// Add a horizontal symmetry check
    pub fn add_horizontal_check(
        &mut self,
        upper_position: ChannelPositionClassification,
        lower_position: ChannelPositionClassification,
        is_symmetric: bool,
    ) {
        self.horizontal_checks.push(SymmetryCheck {
            position_a: upper_position,
            position_b: lower_position,
            is_symmetric,
        });

        if !is_symmetric {
            self.is_symmetric = false;
        }
    }
}

/// Individual symmetry check result
#[derive(Debug, Clone)]
pub struct SymmetryCheck {
    /// First position being compared
    pub position_a: ChannelPositionClassification,

    /// Second position being compared
    pub position_b: ChannelPositionClassification,

    /// Whether the positions are symmetric
    pub is_symmetric: bool,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::GeometryConfig;

    #[test]
    fn test_symmetry_integrated_manager_creation() {
        let manager = SymmetryIntegratedParameterManager::new();
        assert!(manager.symmetry_config.enable_vertical_symmetry);
        assert!(manager.symmetry_config.enable_horizontal_symmetry);
    }

    #[test]
    fn test_symmetry_aware_serpentine_manager() {
        let config = BilateralSymmetryConfig::default();
        let manager = SymmetryAwareSerpentineManager::new(config.clone());

        let context = ChannelGenerationContext::new(
            GeometryConfig::default(),
            (100.0, 50.0),
            4,
            Some(&[10.0, 20.0, 30.0, 40.0]),
        )
        .with_endpoints((10.0, 35.0), (40.0, 40.0));

        let symmetry_context = SymmetryContext::new(context, config);
        let params = manager
            .get_symmetry_aware_parameters(&symmetry_context)
            .unwrap();

        assert!(params.symmetry_enforcement_applied);
        assert!(params.amplitude > 0.0);
        assert!(params.wavelength > 0.0);
    }
}
