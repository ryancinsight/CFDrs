//! Symmetry-aware adaptive parameter implementations
//!
//! Extends base adapters with bilateral symmetry enforcement for
//! amplitude and wavelength parameters across mirrored channel positions.

use crate::state_management::bilateral_symmetry::{
    BilateralSymmetryConfig, ChannelPositionClassification, SymmetryContext,
};

use super::adapters::{DistanceBasedAmplitudeAdapter, LengthBasedWavelengthAdapter};
use super::{AdaptationError, AdaptiveParameter, ChannelGenerationContext};

/// Enhanced symmetry-aware amplitude adapter
#[derive(Debug, Clone)]
pub struct SymmetryAwareAmplitudeAdapter {
    /// Base distance-based amplitude adapter
    pub base_adapter: DistanceBasedAmplitudeAdapter,

    /// Symmetry configuration
    pub symmetry_config: BilateralSymmetryConfig,

    /// Symmetry enforcement factor
    pub symmetry_enforcement_factor: f64,

    /// Enable cross-quadrant amplitude matching
    pub enable_cross_quadrant_matching: bool,
}

impl Default for SymmetryAwareAmplitudeAdapter {
    fn default() -> Self {
        Self {
            base_adapter: DistanceBasedAmplitudeAdapter::default(),
            symmetry_config: BilateralSymmetryConfig::default(),
            symmetry_enforcement_factor: 0.8,
            enable_cross_quadrant_matching: true,
        }
    }
}

impl AdaptiveParameter<f64, ChannelGenerationContext> for SymmetryAwareAmplitudeAdapter {
    fn adapt(
        &self,
        base_value: f64,
        context: &ChannelGenerationContext,
    ) -> Result<f64, AdaptationError> {
        // First apply base distance-based adaptation
        let distance_adapted = self.base_adapter.adapt(base_value, context)?;

        if !self.symmetry_config.enable_adaptive_symmetry {
            return Ok(distance_adapted);
        }

        // Create symmetry context for enhanced symmetry calculations
        let symmetry_context = SymmetryContext::new(context.clone(), self.symmetry_config.clone());

        // Apply symmetry-aware adjustments
        let symmetry_adjusted =
            self.apply_symmetry_adjustments(distance_adapted, &symmetry_context)?;

        Ok(symmetry_adjusted)
    }

    fn is_adaptive(&self) -> bool {
        true
    }

    fn adaptation_description(&self) -> String {
        format!(
            "symmetry-aware amplitude adaptation (enforcement: {}, cross-quadrant: {})",
            self.symmetry_enforcement_factor, self.enable_cross_quadrant_matching
        )
    }

    fn validate_context(&self, context: &ChannelGenerationContext) -> Result<(), AdaptationError> {
        // Validate base context
        self.base_adapter.validate_context(context)?;

        // Validate symmetry-specific requirements
        if self.symmetry_enforcement_factor < 0.0 || self.symmetry_enforcement_factor > 1.0 {
            return Err(AdaptationError::InvalidContext {
                reason: "Symmetry enforcement factor must be between 0.0 and 1.0".to_string(),
            });
        }

        Ok(())
    }
}

impl SymmetryAwareAmplitudeAdapter {
    /// Apply symmetry-aware amplitude adjustments
    fn apply_symmetry_adjustments(
        &self,
        base_amplitude: f64,
        symmetry_context: &SymmetryContext,
    ) -> Result<f64, AdaptationError> {
        let mut adjusted_amplitude = base_amplitude;

        // Apply position-specific symmetry adjustments
        match symmetry_context.position_classification {
            ChannelPositionClassification::UpperLeft
            | ChannelPositionClassification::UpperRight => {
                // Upper channels: ensure consistent amplitude for bilateral symmetry
                if self.enable_cross_quadrant_matching {
                    adjusted_amplitude *= self.symmetry_enforcement_factor.mul_add(0.1, 1.0);
                }
            }
            ChannelPositionClassification::LowerLeft
            | ChannelPositionClassification::LowerRight => {
                // Lower channels: mirror upper channel amplitude adjustments
                if self.enable_cross_quadrant_matching {
                    adjusted_amplitude *= self.symmetry_enforcement_factor.mul_add(0.1, 1.0);
                }
            }
            ChannelPositionClassification::OnHorizontalCenter => {
                // Channels on horizontal centerline: use neutral amplitude
                adjusted_amplitude *= self.symmetry_enforcement_factor.mul_add(-0.05, 1.0);
            }
            _ => {
                // Other positions: minimal adjustment
                adjusted_amplitude *= self.symmetry_enforcement_factor.mul_add(0.02, 1.0);
            }
        }

        // Ensure amplitude remains positive and finite
        if adjusted_amplitude <= 0.0
            || adjusted_amplitude.is_nan()
            || adjusted_amplitude.is_infinite()
        {
            return Err(AdaptationError::InvalidResult {
                value: adjusted_amplitude.to_string(),
                constraint: "Symmetry-adjusted amplitude must be positive and finite".to_string(),
            });
        }

        Ok(adjusted_amplitude)
    }
}

/// Enhanced symmetry-aware wavelength adapter
#[derive(Debug, Clone)]
pub struct SymmetryAwareWavelengthAdapter {
    /// Base length-based wavelength adapter
    pub base_adapter: LengthBasedWavelengthAdapter,

    /// Symmetry configuration
    pub symmetry_config: BilateralSymmetryConfig,

    /// Wavelength synchronization factor for symmetry
    pub wavelength_sync_factor: f64,

    /// Enable wavelength matching across mirror positions
    pub enable_wavelength_matching: bool,
}

impl Default for SymmetryAwareWavelengthAdapter {
    fn default() -> Self {
        Self {
            base_adapter: LengthBasedWavelengthAdapter::default(),
            symmetry_config: BilateralSymmetryConfig::default(),
            wavelength_sync_factor: 0.9,
            enable_wavelength_matching: true,
        }
    }
}

impl AdaptiveParameter<f64, ChannelGenerationContext> for SymmetryAwareWavelengthAdapter {
    fn adapt(
        &self,
        base_value: f64,
        context: &ChannelGenerationContext,
    ) -> Result<f64, AdaptationError> {
        // First apply base length-based adaptation
        let length_adapted = self.base_adapter.adapt(base_value, context)?;

        if !self.symmetry_config.enable_adaptive_symmetry || !self.enable_wavelength_matching {
            return Ok(length_adapted);
        }

        // Create symmetry context for enhanced symmetry calculations
        let symmetry_context = SymmetryContext::new(context.clone(), self.symmetry_config.clone());

        // Apply wavelength synchronization for perfect symmetry
        let synchronized_wavelength =
            self.apply_wavelength_synchronization(length_adapted, &symmetry_context)?;

        Ok(synchronized_wavelength)
    }

    fn is_adaptive(&self) -> bool {
        true
    }

    fn adaptation_description(&self) -> String {
        format!(
            "symmetry-aware wavelength adaptation (sync: {}, matching: {})",
            self.wavelength_sync_factor, self.enable_wavelength_matching
        )
    }

    fn validate_context(&self, context: &ChannelGenerationContext) -> Result<(), AdaptationError> {
        // Validate base context
        self.base_adapter.validate_context(context)?;

        // Validate symmetry-specific requirements
        if self.wavelength_sync_factor < 0.0 || self.wavelength_sync_factor > 1.0 {
            return Err(AdaptationError::InvalidContext {
                reason: "Wavelength sync factor must be between 0.0 and 1.0".to_string(),
            });
        }

        Ok(())
    }
}

impl SymmetryAwareWavelengthAdapter {
    /// Apply wavelength synchronization for perfect bilateral symmetry
    fn apply_wavelength_synchronization(
        &self,
        base_wavelength: f64,
        symmetry_context: &SymmetryContext,
    ) -> Result<f64, AdaptationError> {
        let mut synchronized_wavelength = base_wavelength;

        // Apply position-specific wavelength synchronization
        match symmetry_context.position_classification {
            ChannelPositionClassification::UpperLeft | ChannelPositionClassification::LowerLeft => {
                // Left side channels: use base wavelength as reference
                synchronized_wavelength *= self.wavelength_sync_factor;
            }
            ChannelPositionClassification::UpperRight
            | ChannelPositionClassification::LowerRight => {
                // Right side channels: synchronize with left side for perfect bilateral symmetry
                synchronized_wavelength *= self.wavelength_sync_factor;
            }
            ChannelPositionClassification::OnVerticalCenter => {
                // Channels on vertical centerline: use neutral wavelength
                synchronized_wavelength *= 1.0;
            }
            _ => {
                // Other positions: minimal synchronization
                synchronized_wavelength *= self.wavelength_sync_factor.mul_add(0.05, 0.95);
            }
        }

        // Ensure wavelength remains positive and finite
        if synchronized_wavelength <= 0.0
            || synchronized_wavelength.is_nan()
            || synchronized_wavelength.is_infinite()
        {
            return Err(AdaptationError::InvalidResult {
                value: synchronized_wavelength.to_string(),
                constraint: "Synchronized wavelength must be positive and finite".to_string(),
            });
        }

        Ok(synchronized_wavelength)
    }
}
