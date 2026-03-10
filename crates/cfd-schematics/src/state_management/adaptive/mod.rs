//! Adaptive parameter behavior system
//!
//! This module provides traits and implementations for parameters that can
//! adapt their values based on context, such as channel generation context,
//! neighbor proximity, or geometric constraints.

mod adapters;
mod symmetry;

pub use adapters::{
    BranchCountDensityAdapter, DistanceBasedAmplitudeAdapter, LengthBasedWavelengthAdapter,
};
pub use symmetry::{SymmetryAwareAmplitudeAdapter, SymmetryAwareWavelengthAdapter};

use crate::config::GeometryConfig;
use crate::geometry::Point2D;
use std::fmt::Debug;

/// Context information for channel generation
#[derive(Debug, Clone)]
pub struct ChannelGenerationContext {
    /// Geometry configuration
    pub geometry_config: GeometryConfig,

    /// Bounding box dimensions (width, height)
    pub box_dims: (f64, f64),

    /// Total number of branches in the system
    pub total_branches: usize,

    /// Information about neighboring channels (y-coordinates)
    pub neighbor_info: Option<Vec<f64>>,

    /// Channel start and end points
    pub channel_endpoints: (Point2D, Point2D),

    /// Current channel index in the system
    pub channel_index: usize,

    /// Additional context data
    pub custom_data: std::collections::HashMap<String, f64>,
}

impl ChannelGenerationContext {
    /// Create a new channel generation context
    #[must_use]
    pub fn new(
        geometry_config: GeometryConfig,
        box_dims: (f64, f64),
        total_branches: usize,
        neighbor_info: Option<&[f64]>,
    ) -> Self {
        Self {
            geometry_config,
            box_dims,
            total_branches,
            neighbor_info: neighbor_info.map(<[f64]>::to_vec),
            channel_endpoints: ((0.0, 0.0), (0.0, 0.0)),
            channel_index: 0,
            custom_data: std::collections::HashMap::new(),
        }
    }

    /// Set channel endpoints
    #[must_use]
    pub const fn with_endpoints(mut self, start: Point2D, end: Point2D) -> Self {
        self.channel_endpoints = (start, end);
        self
    }

    /// Set channel index
    #[must_use]
    pub const fn with_index(mut self, index: usize) -> Self {
        self.channel_index = index;
        self
    }

    /// Add custom data
    #[must_use]
    pub fn with_custom_data(mut self, key: &str, value: f64) -> Self {
        self.custom_data.insert(key.to_string(), value);
        self
    }

    /// Get channel length
    #[must_use]
    pub fn channel_length(&self) -> f64 {
        let (start, end) = self.channel_endpoints;
        let dx = end.0 - start.0;
        let dy = end.1 - start.1;
        dx.hypot(dy)
    }

    /// Get channel center point
    #[must_use]
    pub const fn channel_center(&self) -> Point2D {
        let (start, end) = self.channel_endpoints;
        (f64::midpoint(start.0, end.0), f64::midpoint(start.1, end.1))
    }

    /// Check if channel is mostly horizontal
    #[must_use]
    pub fn is_horizontal(&self) -> bool {
        let (start, end) = self.channel_endpoints;
        let dx = (end.0 - start.0).abs();
        let dy = (end.1 - start.1).abs();
        dx > dy * 2.0 // Horizontal if dx is significantly larger than dy
    }

    /// Get minimum distance to neighbors
    ///
    /// Returns None if no neighbors exist or if distance calculation fails
    #[must_use]
    pub fn min_neighbor_distance(&self) -> Option<f64> {
        let center = self.channel_center();
        self.neighbor_info.as_ref().and_then(|neighbors| {
            neighbors
                .iter()
                .map(|&neighbor_y| (neighbor_y - center.1).abs())
                .filter(|&dist| dist > 1e-6) // Exclude self
                .min_by(|a, b| {
                    // Handle NaN values gracefully instead of panicking
                    if let Some(ordering) = a.partial_cmp(b) {
                        ordering
                    } else {
                        // Log warning for NaN values in debug builds
                        #[cfg(debug_assertions)]
                        ::tracing::info!(
                            "Warning: NaN encountered in neighbor distance calculation"
                        );
                        std::cmp::Ordering::Equal
                    }
                })
        })
    }

    /// Get distance to walls
    #[must_use]
    pub fn wall_distances(&self) -> (f64, f64, f64, f64) {
        let center = self.channel_center();
        let (width, height) = self.box_dims;
        let half_channel = self.geometry_config.channel_width / 2.0;

        (
            center.0 - half_channel,          // left
            width - center.0 - half_channel,  // right
            center.1 - half_channel,          // bottom
            height - center.1 - half_channel, // top
        )
    }

    /// Get minimum distance to any wall
    #[must_use]
    pub fn min_wall_distance(&self) -> f64 {
        let (left, right, bottom, top) = self.wall_distances();
        left.min(right).min(bottom).min(top)
    }
}

/// Error type for adaptation failures
#[derive(Debug, Clone)]
pub enum AdaptationError {
    /// Invalid context provided
    InvalidContext { reason: String },

    /// Adaptation calculation failed
    CalculationFailed { parameter: String, reason: String },

    /// Result value is invalid
    InvalidResult { value: String, constraint: String },

    /// Dependency not available
    DependencyMissing { dependency: String },
}

impl std::fmt::Display for AdaptationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidContext { reason } => {
                write!(f, "Invalid adaptation context: {reason}")
            }
            Self::CalculationFailed { parameter, reason } => {
                write!(
                    f,
                    "Adaptation calculation failed for parameter '{parameter}': {reason}"
                )
            }
            Self::InvalidResult { value, constraint } => {
                write!(
                    f,
                    "Adaptation result '{value}' violates constraint: {constraint}"
                )
            }
            Self::DependencyMissing { dependency } => {
                write!(
                    f,
                    "Required dependency '{dependency}' not available for adaptation"
                )
            }
        }
    }
}

impl std::error::Error for AdaptationError {}

/// Convenience extension trait for fallible adaptive parameters.
pub trait AdaptiveParameterCompat<T, Context>: AdaptiveParameter<T, Context> {
    /// Adapt with fallback to default value on error
    fn adapt_or_default(&self, base_value: T, context: &Context, default: T) -> T
    where
        T: Clone,
    {
        if let Ok(adapted) = self.adapt(base_value, context) {
            adapted
        } else {
            #[cfg(debug_assertions)]
            ::tracing::info!("Warning: Adaptation failed, using default value");
            default
        }
    }

    /// Adapt with fallback to base value on error
    fn adapt_or_base(&self, base_value: T, context: &Context) -> T
    where
        T: Clone,
    {
        if let Ok(adapted) = self.adapt(base_value.clone(), context) {
            adapted
        } else {
            #[cfg(debug_assertions)]
            ::tracing::info!("Warning: Adaptation failed, using base value");
            base_value
        }
    }
}

// Blanket implementation for all AdaptiveParameter types
impl<T, Context, A> AdaptiveParameterCompat<T, Context> for A where A: AdaptiveParameter<T, Context> {}

/// Trait for parameters that can adapt based on context with proper error handling
pub trait AdaptiveParameter<T, Context>: Debug + Send + Sync {
    /// Calculate adaptive value based on context
    ///
    /// # Errors
    ///
    /// Returns an error if adaptation fails due to invalid context, calculation errors,
    /// constraint violations, or missing dependencies.
    fn adapt(&self, base_value: T, context: &Context) -> Result<T, AdaptationError>;

    /// Check if adaptation is enabled
    fn is_adaptive(&self) -> bool;

    /// Get adaptation description
    fn adaptation_description(&self) -> String;

    /// Validate that the context is suitable for adaptation
    ///
    /// # Errors
    ///
    /// Returns an error if the context is invalid or unsuitable for adaptation.
    fn validate_context(&self, context: &Context) -> Result<(), AdaptationError> {
        // Default implementation - can be overridden
        let _ = context; // Suppress unused parameter warning
        Ok(())
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
    fn test_distance_based_amplitude_adapter() {
        let adapter = DistanceBasedAmplitudeAdapter::default();
        let context = create_test_context();

        let base_amplitude = 10.0;
        let adapted = adapter
            .adapt(base_amplitude, &context)
            .expect("structural invariant");

        // Should be scaled down due to constraints
        assert!(adapted < base_amplitude);
        assert!(adapted > 0.0);
    }

    #[test]
    fn test_branch_count_density_adapter() {
        let adapter = BranchCountDensityAdapter::default();
        let context = create_test_context();

        let base_density = 2.0;
        let adapted = adapter
            .adapt(base_density, &context)
            .expect("structural invariant");

        // Should be scaled based on branch count
        assert!(adapted > 0.0);
        assert!(AdaptiveParameter::is_adaptive(&adapter));
    }

    #[test]
    fn test_length_based_wavelength_adapter() {
        let adapter = LengthBasedWavelengthAdapter::default();
        let context = create_test_context();

        let base_wavelength = 1.0;
        let adapted = adapter
            .adapt(base_wavelength, &context)
            .expect("structural invariant");

        // Should adjust based on channel length
        assert!(adapted > 0.0);
        assert!(AdaptiveParameter::is_adaptive(&adapter));
    }
}
