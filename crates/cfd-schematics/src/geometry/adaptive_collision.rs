//! Adaptive collision detection and avoidance system
//!
//! This module provides advanced collision detection capabilities that adapt
//! their behavior based on channel generation context, following SOLID, CUPID,
//! and GRASP design principles for maximum flexibility and maintainability.

use crate::{
    config::ConstantsRegistry,
    geometry::collision_detection::CollisionDetectionResult,
    state_management::{
        adaptive::{
            AdaptationError, AdaptiveParameter, AdaptiveParameterCompat, ChannelGenerationContext,
        },
        constraints::ParameterConstraints,
        parameters::ConfigurableParameter,
    },
};
use std::fmt::Debug;

/// Adaptive collision detection strategy that adjusts behavior based on context
#[derive(Debug)]
pub struct AdaptiveCollisionDetector {
    /// Base collision parameters
    base_parameters: CollisionParameterSet,

    /// Adaptive behaviors for different parameters
    distance_adapter: DistanceBasedCollisionAdapter,
    sensitivity_adapter: ContextSensitivityAdapter,
    avoidance_adapter: AdaptiveAvoidanceStrategy,
}

/// Set of collision parameters with adaptive capabilities
#[derive(Debug)]
pub struct CollisionParameterSet {
    /// Minimum channel distance parameter
    pub min_channel_distance: ConfigurableParameter<f64>,

    /// Minimum wall distance parameter
    pub min_wall_distance: ConfigurableParameter<f64>,

    /// Safety margin factor parameter
    pub safety_margin_factor: ConfigurableParameter<f64>,

    /// Detection sensitivity parameter
    pub detection_sensitivity: ConfigurableParameter<f64>,

    /// Maximum reduction factor parameter
    pub max_reduction_factor: ConfigurableParameter<f64>,
}

impl CollisionParameterSet {
    /// Create a new parameter set with default values from constants registry
    #[must_use]
    pub fn from_constants_registry(_constants: &ConstantsRegistry) -> Self {
        Self {
            min_channel_distance: ConfigurableParameter::new(
                2.0, // This would come from constants in full implementation
                ParameterConstraints::range(0.5, 10.0),
                crate::state_management::parameters::ParameterMetadata::new(
                    "min_channel_distance",
                    "Minimum distance between channels for collision detection",
                    "collision_detection",
                )
                .with_units("mm"),
            ),

            min_wall_distance: ConfigurableParameter::new(
                1.0,
                ParameterConstraints::range(0.1, 5.0),
                crate::state_management::parameters::ParameterMetadata::new(
                    "min_wall_distance",
                    "Minimum distance from walls for collision detection",
                    "collision_detection",
                )
                .with_units("mm"),
            ),

            safety_margin_factor: ConfigurableParameter::new(
                1.2,
                ParameterConstraints::range(1.0, 3.0),
                crate::state_management::parameters::ParameterMetadata::new(
                    "safety_margin_factor",
                    "Safety margin multiplier for collision distances",
                    "collision_detection",
                )
                .with_units("factor"),
            ),

            detection_sensitivity: ConfigurableParameter::new(
                1.0,
                ParameterConstraints::range(0.1, 5.0),
                crate::state_management::parameters::ParameterMetadata::new(
                    "detection_sensitivity",
                    "Sensitivity factor for collision detection",
                    "collision_detection",
                )
                .with_units("factor"),
            ),

            max_reduction_factor: ConfigurableParameter::new(
                0.8,
                ParameterConstraints::range(0.1, 0.95),
                crate::state_management::parameters::ParameterMetadata::new(
                    "max_reduction_factor",
                    "Maximum reduction factor for collision avoidance",
                    "collision_detection",
                )
                .with_units("factor"),
            ),
        }
    }

    /// Get parameter values with optional adaptive context
    #[must_use]
    pub fn get_values(
        &self,
        context: Option<&ChannelGenerationContext>,
    ) -> CollisionParameterValues {
        CollisionParameterValues {
            min_channel_distance: self.min_channel_distance.get_value(context),
            min_wall_distance: self.min_wall_distance.get_value(context),
            safety_margin_factor: self.safety_margin_factor.get_value(context),
            detection_sensitivity: self.detection_sensitivity.get_value(context),
            max_reduction_factor: self.max_reduction_factor.get_value(context),
        }
    }
}

/// Concrete parameter values for collision detection
#[derive(Debug, Clone)]
pub struct CollisionParameterValues {
    pub min_channel_distance: f64,
    pub min_wall_distance: f64,
    pub safety_margin_factor: f64,
    pub detection_sensitivity: f64,
    pub max_reduction_factor: f64,
}

/// Adaptive collision distance calculator based on neighbor proximity
#[derive(Debug)]
pub struct DistanceBasedCollisionAdapter {
    /// Scaling factor for neighbor-based adjustments
    pub neighbor_scale_factor: f64,

    /// Minimum distance threshold
    pub min_distance_threshold: f64,

    /// Maximum distance adjustment factor
    pub max_adjustment_factor: f64,
}

impl Default for DistanceBasedCollisionAdapter {
    fn default() -> Self {
        let constants = crate::config::ConstantsRegistry::new();
        Self {
            neighbor_scale_factor: constants.get_neighbor_scale_factor(),
            min_distance_threshold: constants.get_min_distance_threshold(),
            max_adjustment_factor: constants.get_max_adjustment_factor(),
        }
    }
}

impl AdaptiveParameter<f64, ChannelGenerationContext> for DistanceBasedCollisionAdapter {
    fn adapt(
        &self,
        base_value: f64,
        context: &ChannelGenerationContext,
    ) -> Result<f64, AdaptationError> {
        context
            .min_neighbor_distance()
            .map_or(Ok(base_value), |min_neighbor_dist| {
                // Adjust collision distance based on neighbor proximity
                let proximity_factor = (min_neighbor_dist / self.min_distance_threshold)
                    .min(self.max_adjustment_factor);
                let adjusted_value = base_value * proximity_factor * self.neighbor_scale_factor;

                // Ensure minimum threshold
                Ok(adjusted_value.max(self.min_distance_threshold))
            })
    }

    fn is_adaptive(&self) -> bool {
        true
    }

    fn adaptation_description(&self) -> String {
        format!(
            "distance-based collision adaptation (scale: {}, threshold: {})",
            self.neighbor_scale_factor, self.min_distance_threshold
        )
    }
}

/// Context-sensitive collision detection sensitivity adapter
#[derive(Debug)]
pub struct ContextSensitivityAdapter {
    /// Branch count scaling factor
    pub branch_scale_factor: f64,

    /// Channel length scaling factor
    pub length_scale_factor: f64,

    /// Maximum sensitivity multiplier
    pub max_sensitivity: f64,
}

impl Default for ContextSensitivityAdapter {
    fn default() -> Self {
        Self {
            branch_scale_factor: 0.1,
            length_scale_factor: 0.01,
            max_sensitivity: 3.0,
        }
    }
}

impl AdaptiveParameter<f64, ChannelGenerationContext> for ContextSensitivityAdapter {
    fn adapt(
        &self,
        base_value: f64,
        context: &ChannelGenerationContext,
    ) -> Result<f64, AdaptationError> {
        let mut sensitivity = base_value;

        // Adjust based on branch count
        let branch_factor = (context.total_branches as f64).mul_add(self.branch_scale_factor, 1.0);
        sensitivity *= branch_factor;

        // Adjust based on channel length
        let channel_length = context.channel_length();
        let length_factor = channel_length.mul_add(self.length_scale_factor, 1.0);
        sensitivity *= length_factor;

        // Cap at maximum sensitivity
        Ok(sensitivity.min(self.max_sensitivity))
    }

    fn is_adaptive(&self) -> bool {
        true
    }

    fn adaptation_description(&self) -> String {
        format!(
            "context-sensitive detection (branch: {}, length: {})",
            self.branch_scale_factor, self.length_scale_factor
        )
    }
}

/// Adaptive collision avoidance strategy
#[derive(Debug)]
pub struct AdaptiveAvoidanceStrategy {
    /// Base reduction factor
    pub base_reduction: f64,

    /// Severity scaling factor
    pub severity_scale: f64,

    /// Context adjustment factor
    pub context_adjustment: f64,
}

impl Default for AdaptiveAvoidanceStrategy {
    fn default() -> Self {
        Self {
            base_reduction: 0.5,
            severity_scale: 0.3,
            context_adjustment: 0.2,
        }
    }
}

impl AdaptiveAvoidanceStrategy {
    /// Calculate adaptive reduction factor based on collision result and context
    #[must_use]
    pub fn calculate_reduction_factor(
        &self,
        collision_result: &CollisionDetectionResult,
        context: &ChannelGenerationContext,
        max_reduction: f64,
    ) -> f64 {
        // Base reduction from collision severity
        let severity_reduction = collision_result.severity_score * self.severity_scale;

        // Context-based adjustments
        let context_factor = self.calculate_context_factor(context);

        // Combine factors
        let total_reduction = self.base_reduction + severity_reduction + context_factor;

        // Cap at maximum allowed reduction
        total_reduction.min(max_reduction)
    }

    /// Calculate context-based adjustment factor
    fn calculate_context_factor(&self, context: &ChannelGenerationContext) -> f64 {
        let mut factor = 0.0;

        // Adjust based on branch density
        let branch_density = context.total_branches as f64 / 10.0; // Normalize to typical range
        factor += branch_density * self.context_adjustment;

        // Adjust based on neighbor proximity
        if let Some(min_neighbor_dist) = context.min_neighbor_distance() {
            if min_neighbor_dist < 5.0 {
                factor += (5.0 - min_neighbor_dist) * 0.1;
            }
        }

        // Adjust based on channel length
        let channel_length = context.channel_length();
        if channel_length > 50.0 {
            factor += 0.1; // Longer channels can handle more reduction
        }

        factor
    }
}

impl AdaptiveCollisionDetector {
    /// Create a new adaptive collision detector
    #[must_use]
    pub fn new(constants: &ConstantsRegistry) -> Self {
        Self {
            base_parameters: CollisionParameterSet::from_constants_registry(constants),
            distance_adapter: DistanceBasedCollisionAdapter::default(),
            sensitivity_adapter: ContextSensitivityAdapter::default(),
            avoidance_adapter: AdaptiveAvoidanceStrategy::default(),
        }
    }

    /// Get adaptive collision parameters for a given context
    #[must_use]
    pub fn get_adaptive_parameters(
        &self,
        context: &ChannelGenerationContext,
    ) -> CollisionParameterValues {
        let mut values = self.base_parameters.get_values(Some(context));

        // Apply adaptive adjustments with error handling
        values.min_channel_distance = self
            .distance_adapter
            .adapt_or_base(values.min_channel_distance, context);
        values.min_wall_distance = self
            .distance_adapter
            .adapt_or_base(values.min_wall_distance, context);
        values.detection_sensitivity = self
            .sensitivity_adapter
            .adapt_or_base(values.detection_sensitivity, context);

        values
    }

    /// Calculate adaptive avoidance strategy
    #[must_use]
    pub fn calculate_adaptive_avoidance(
        &self,
        collision_result: &CollisionDetectionResult,
        context: &ChannelGenerationContext,
    ) -> f64 {
        let max_reduction = self
            .base_parameters
            .max_reduction_factor
            .get_value(Some(context));
        self.avoidance_adapter
            .calculate_reduction_factor(collision_result, context, max_reduction)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::GeometryConfig;

    fn create_test_context() -> ChannelGenerationContext {
        ChannelGenerationContext::new(
            GeometryConfig::default(),
            (100.0, 50.0),
            4,
            Some(&[10.0, 20.0, 30.0, 40.0]),
        )
        .with_endpoints((0.0, 25.0), (100.0, 25.0))
    }

    #[test]
    fn test_adaptive_collision_detector_creation() {
        let constants = ConstantsRegistry::new();
        let detector = AdaptiveCollisionDetector::new(&constants);

        // Test that detector is created successfully
        assert!(detector.distance_adapter.is_adaptive());
        assert!(detector.sensitivity_adapter.is_adaptive());
    }

    #[test]
    fn test_distance_based_adaptation() {
        let adapter = DistanceBasedCollisionAdapter::default();
        let context = create_test_context();

        let base_distance = 2.0;
        let adapted_distance = adapter.adapt(base_distance, &context).unwrap();

        // Should adapt based on neighbor proximity
        assert!(adapted_distance >= adapter.min_distance_threshold);
    }

    #[test]
    fn test_context_sensitivity_adaptation() {
        let adapter = ContextSensitivityAdapter::default();
        let context = create_test_context();

        let base_sensitivity = 1.0;
        let adapted_sensitivity = adapter.adapt(base_sensitivity, &context).unwrap();

        // Should increase sensitivity based on context
        assert!(adapted_sensitivity >= base_sensitivity);
        assert!(adapted_sensitivity <= adapter.max_sensitivity);
    }
}
