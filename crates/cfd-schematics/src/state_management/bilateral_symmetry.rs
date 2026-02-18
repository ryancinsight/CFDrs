//! Enhanced bilateral mirror symmetry system
//!
//! This module provides comprehensive bilateral mirror symmetry capabilities that ensure
//! perfect symmetry across both vertical and horizontal centerlines. It integrates with
//! the adaptive parameter system to provide context-aware symmetry adjustments while
//! following SOLID, CUPID, GRASP, and all other specified design principles.

#![allow(unused_imports)] // Allow unused imports for comprehensive framework

use crate::{
    geometry::Point2D, // May be used in future symmetry operations
    state_management::adaptive::{AdaptationError, AdaptiveParameter, ChannelGenerationContext},
};
use std::collections::HashMap;

/// Enhanced bilateral mirror symmetry configuration
#[derive(Debug, Clone)]
pub struct BilateralSymmetryConfig {
    /// Enable perfect vertical centerline symmetry (splits mirror merges)
    pub enable_vertical_symmetry: bool,

    /// Enable perfect horizontal centerline symmetry (upper mirrors lower)
    pub enable_horizontal_symmetry: bool,

    /// Symmetry tolerance for validation (in units)
    pub symmetry_tolerance: f64,

    /// Enable adaptive symmetry adjustments
    pub enable_adaptive_symmetry: bool,

    /// Symmetry enforcement strength (0.0 = weak, 1.0 = strict)
    pub enforcement_strength: f64,
}

impl Default for BilateralSymmetryConfig {
    fn default() -> Self {
        Self {
            enable_vertical_symmetry: true,
            enable_horizontal_symmetry: true,
            symmetry_tolerance: 1e-6,
            enable_adaptive_symmetry: true,
            enforcement_strength: 1.0,
        }
    }
}

/// Bilateral symmetry context for enhanced symmetry calculations
#[derive(Debug, Clone)]
pub struct SymmetryContext {
    /// Channel generation context
    pub channel_context: ChannelGenerationContext,

    /// Vertical centerline position
    pub vertical_centerline: f64,

    /// Horizontal centerline position
    pub horizontal_centerline: f64,

    /// Channel position relative to centerlines
    pub position_classification: ChannelPositionClassification,

    /// Symmetry configuration
    pub config: BilateralSymmetryConfig,
}

/// Classification of channel position relative to centerlines
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ChannelPositionClassification {
    /// Upper left quadrant (split, above horizontal center)
    UpperLeft,

    /// Upper right quadrant (merge, above horizontal center)
    UpperRight,

    /// Lower left quadrant (split, below horizontal center)
    LowerLeft,

    /// Lower right quadrant (merge, below horizontal center)
    LowerRight,

    /// On vertical centerline
    OnVerticalCenter,

    /// On horizontal centerline
    OnHorizontalCenter,

    /// At intersection of both centerlines
    AtIntersection,
}

impl SymmetryContext {
    /// Create a new symmetry context from channel generation context
    #[must_use]
    pub fn new(channel_context: ChannelGenerationContext, config: BilateralSymmetryConfig) -> Self {
        let vertical_centerline = channel_context.box_dims.0 / 2.0;
        let horizontal_centerline = channel_context.box_dims.1 / 2.0;

        let position_classification = Self::classify_channel_position(
            &channel_context,
            vertical_centerline,
            horizontal_centerline,
            config.symmetry_tolerance,
        );

        Self {
            channel_context,
            vertical_centerline,
            horizontal_centerline,
            position_classification,
            config,
        }
    }

    /// Classify channel position relative to centerlines
    fn classify_channel_position(
        context: &ChannelGenerationContext,
        vertical_center: f64,
        horizontal_center: f64,
        tolerance: f64,
    ) -> ChannelPositionClassification {
        let (start, end) = context.channel_endpoints;
        let center_x = f64::midpoint(start.0, end.0);
        let center_y = f64::midpoint(start.1, end.1);

        let on_vertical = (center_x - vertical_center).abs() < tolerance;
        let on_horizontal = (center_y - horizontal_center).abs() < tolerance;

        if on_vertical && on_horizontal {
            ChannelPositionClassification::AtIntersection
        } else if on_vertical {
            ChannelPositionClassification::OnVerticalCenter
        } else if on_horizontal {
            ChannelPositionClassification::OnHorizontalCenter
        } else {
            let is_left = center_x < vertical_center;
            let is_upper = center_y > horizontal_center;

            match (is_left, is_upper) {
                (true, true) => ChannelPositionClassification::UpperLeft,
                (false, true) => ChannelPositionClassification::UpperRight,
                (true, false) => ChannelPositionClassification::LowerLeft,
                (false, false) => ChannelPositionClassification::LowerRight,
            }
        }
    }

    /// Get the mirror position classification for perfect bilateral symmetry
    #[must_use]
    pub const fn get_mirror_position(&self) -> ChannelPositionClassification {
        match self.position_classification {
            ChannelPositionClassification::UpperLeft => ChannelPositionClassification::UpperRight,
            ChannelPositionClassification::UpperRight => ChannelPositionClassification::UpperLeft,
            ChannelPositionClassification::LowerLeft => ChannelPositionClassification::LowerRight,
            ChannelPositionClassification::LowerRight => ChannelPositionClassification::LowerLeft,
            // Centerline positions mirror to themselves
            other => other,
        }
    }

    /// Get the horizontal mirror position classification
    #[must_use]
    pub const fn get_horizontal_mirror_position(&self) -> ChannelPositionClassification {
        match self.position_classification {
            ChannelPositionClassification::UpperLeft => ChannelPositionClassification::LowerLeft,
            ChannelPositionClassification::UpperRight => ChannelPositionClassification::LowerRight,
            ChannelPositionClassification::LowerLeft => ChannelPositionClassification::UpperLeft,
            ChannelPositionClassification::LowerRight => ChannelPositionClassification::UpperRight,
            // Centerline positions mirror to themselves or their counterparts
            ChannelPositionClassification::OnVerticalCenter => {
                ChannelPositionClassification::OnVerticalCenter
            }
            other => other,
        }
    }
}

/// Enhanced phase direction calculator with perfect bilateral symmetry
#[derive(Debug, Clone)]
pub struct BilateralPhaseDirectionCalculator {
    /// Symmetry configuration
    pub config: BilateralSymmetryConfig,

    /// Phase direction mapping for each position classification
    pub phase_direction_map: HashMap<ChannelPositionClassification, f64>,
}

impl Default for BilateralPhaseDirectionCalculator {
    fn default() -> Self {
        let mut phase_direction_map = HashMap::new();

        // Perfect bilateral symmetry mapping:
        // Upper channels have positive phase, lower channels have negative phase
        // This ensures perfect mirror symmetry across horizontal centerline
        phase_direction_map.insert(ChannelPositionClassification::UpperLeft, 1.0);
        phase_direction_map.insert(ChannelPositionClassification::UpperRight, 1.0);
        phase_direction_map.insert(ChannelPositionClassification::LowerLeft, -1.0);
        phase_direction_map.insert(ChannelPositionClassification::LowerRight, -1.0);

        // Centerline channels use neutral phase
        phase_direction_map.insert(ChannelPositionClassification::OnVerticalCenter, 0.0);
        phase_direction_map.insert(ChannelPositionClassification::OnHorizontalCenter, 0.0);
        phase_direction_map.insert(ChannelPositionClassification::AtIntersection, 0.0);

        Self {
            config: BilateralSymmetryConfig::default(),
            phase_direction_map,
        }
    }
}

impl BilateralPhaseDirectionCalculator {
    /// Create a new phase direction calculator with custom configuration
    #[must_use]
    pub fn new(config: BilateralSymmetryConfig) -> Self {
        Self { config, ..Default::default() }
    }

    /// Calculate phase direction for perfect bilateral symmetry
    pub fn calculate_phase_direction(
        &self,
        context: &SymmetryContext,
    ) -> Result<f64, AdaptationError> {
        if !self.config.enable_vertical_symmetry && !self.config.enable_horizontal_symmetry {
            return Ok(0.0); // No symmetry enforcement
        }

        let base_phase = self
            .phase_direction_map
            .get(&context.position_classification)
            .copied()
            .unwrap_or(0.0);

        // Apply enforcement strength
        let enforced_phase = base_phase * self.config.enforcement_strength;

        // Validate phase direction
        if enforced_phase.is_nan() || enforced_phase.is_infinite() {
            return Err(AdaptationError::CalculationFailed {
                parameter: "phase_direction".to_string(),
                reason: "Phase direction calculation resulted in invalid value".to_string(),
            });
        }

        Ok(enforced_phase)
    }

    /// Validate bilateral symmetry for a pair of channels
    pub fn validate_bilateral_symmetry(
        &self,
        left_context: &SymmetryContext,
        right_context: &SymmetryContext,
    ) -> Result<bool, AdaptationError> {
        if !self.config.enable_vertical_symmetry {
            return Ok(true); // Symmetry validation disabled
        }

        let left_phase = self.calculate_phase_direction(left_context)?;
        let right_phase = self.calculate_phase_direction(right_context)?;

        // For perfect bilateral symmetry, corresponding channels should have the same phase
        let phase_difference = (left_phase - right_phase).abs();
        let is_symmetric = phase_difference < self.config.symmetry_tolerance;

        Ok(is_symmetric)
    }

    /// Validate horizontal symmetry for a pair of channels
    pub fn validate_horizontal_symmetry(
        &self,
        upper_context: &SymmetryContext,
        lower_context: &SymmetryContext,
    ) -> Result<bool, AdaptationError> {
        if !self.config.enable_horizontal_symmetry {
            return Ok(true); // Symmetry validation disabled
        }

        let upper_phase = self.calculate_phase_direction(upper_context)?;
        let lower_phase = self.calculate_phase_direction(lower_context)?;

        // For perfect horizontal symmetry, upper and lower channels should have opposite phases
        let expected_phase_difference = 2.0; // |1.0 - (-1.0)| = 2.0
        let actual_phase_difference = (upper_phase - lower_phase).abs();
        let difference_error = (actual_phase_difference - expected_phase_difference).abs();

        let is_symmetric = difference_error < self.config.symmetry_tolerance;

        Ok(is_symmetric)
    }
}

impl AdaptiveParameter<f64, SymmetryContext> for BilateralPhaseDirectionCalculator {
    fn adapt(&self, _base_value: f64, context: &SymmetryContext) -> Result<f64, AdaptationError> {
        self.calculate_phase_direction(context)
    }

    fn is_adaptive(&self) -> bool {
        self.config.enable_adaptive_symmetry
    }

    fn adaptation_description(&self) -> String {
        format!(
            "bilateral symmetry phase direction (vertical: {}, horizontal: {}, strength: {})",
            self.config.enable_vertical_symmetry,
            self.config.enable_horizontal_symmetry,
            self.config.enforcement_strength
        )
    }

    fn validate_context(&self, context: &SymmetryContext) -> Result<(), AdaptationError> {
        // Validate symmetry context
        if context.vertical_centerline <= 0.0 || context.horizontal_centerline <= 0.0 {
            return Err(AdaptationError::InvalidContext {
                reason: "Centerline positions must be positive".to_string(),
            });
        }

        if context.config.symmetry_tolerance <= 0.0 {
            return Err(AdaptationError::InvalidContext {
                reason: "Symmetry tolerance must be positive".to_string(),
            });
        }

        if context.config.enforcement_strength < 0.0 || context.config.enforcement_strength > 1.0 {
            return Err(AdaptationError::InvalidContext {
                reason: "Enforcement strength must be between 0.0 and 1.0".to_string(),
            });
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::GeometryConfig;

    fn create_test_context(start: Point2D, end: Point2D) -> SymmetryContext {
        let channel_context = ChannelGenerationContext::new(
            GeometryConfig::default(),
            (100.0, 50.0),
            4,
            Some(&[10.0, 20.0, 30.0, 40.0]),
        )
        .with_endpoints(start, end);

        SymmetryContext::new(channel_context, BilateralSymmetryConfig::default())
    }

    #[test]
    fn test_channel_position_classification() {
        // Test upper left quadrant
        let context = create_test_context((10.0, 35.0), (40.0, 40.0));
        assert_eq!(
            context.position_classification,
            ChannelPositionClassification::UpperLeft
        );

        // Test lower right quadrant
        let context = create_test_context((60.0, 10.0), (90.0, 20.0));
        assert_eq!(
            context.position_classification,
            ChannelPositionClassification::LowerRight
        );
    }

    #[test]
    fn test_bilateral_phase_direction_calculation() {
        let calculator = BilateralPhaseDirectionCalculator::default();

        // Test upper channel (should have positive phase)
        let upper_context = create_test_context((10.0, 35.0), (40.0, 40.0));
        let upper_phase = calculator
            .calculate_phase_direction(&upper_context)
            .unwrap();
        assert_eq!(upper_phase, 1.0);

        // Test lower channel (should have negative phase)
        let lower_context = create_test_context((10.0, 10.0), (40.0, 15.0));
        let lower_phase = calculator
            .calculate_phase_direction(&lower_context)
            .unwrap();
        assert_eq!(lower_phase, -1.0);
    }

    #[test]
    fn test_bilateral_symmetry_validation() {
        let calculator = BilateralPhaseDirectionCalculator::default();

        // Test symmetric pair (both upper)
        let left_context = create_test_context((10.0, 35.0), (40.0, 40.0));
        let right_context = create_test_context((60.0, 35.0), (90.0, 40.0));

        let is_symmetric = calculator
            .validate_bilateral_symmetry(&left_context, &right_context)
            .unwrap();
        assert!(is_symmetric);
    }

    #[test]
    fn test_horizontal_symmetry_validation() {
        let calculator = BilateralPhaseDirectionCalculator::default();

        // Test horizontally symmetric pair
        let upper_context = create_test_context((10.0, 35.0), (40.0, 40.0));
        let lower_context = create_test_context((10.0, 10.0), (40.0, 15.0));

        let is_symmetric = calculator
            .validate_horizontal_symmetry(&upper_context, &lower_context)
            .unwrap();
        assert!(is_symmetric);
    }
}
