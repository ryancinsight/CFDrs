//! Adaptive parameter behavior system
//!
//! This module provides traits and implementations for parameters that can
//! adapt their values based on context, such as channel generation context,
//! neighbor proximity, or geometric constraints.

use crate::{
    config::GeometryConfig,
    geometry::Point2D,
    state_management::bilateral_symmetry::{
        BilateralSymmetryConfig, ChannelPositionClassification, SymmetryContext,
    },
    // Removed unused error imports
};
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
                        eprintln!("Warning: NaN encountered in neighbor distance calculation");
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

/// Extension trait for backward compatibility with existing code
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
            eprintln!("Warning: Adaptation failed, using default value");
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
            eprintln!("Warning: Adaptation failed, using base value");
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

/// Distance-based adaptive behavior for amplitude parameters
#[derive(Debug, Clone)]
pub struct DistanceBasedAmplitudeAdapter {
    /// Enable neighbor-based scaling
    pub neighbor_scaling: bool,

    /// Enable wall-based scaling
    pub wall_scaling: bool,

    /// Scaling factor for neighbor proximity (0.0 to 1.0)
    pub neighbor_scale_factor: f64,

    /// Scaling factor for wall proximity (0.0 to 1.0)
    pub wall_scale_factor: f64,

    /// Minimum allowed amplitude ratio
    pub min_amplitude_ratio: f64,
}

impl Default for DistanceBasedAmplitudeAdapter {
    fn default() -> Self {
        Self {
            neighbor_scaling: true,
            wall_scaling: true,
            neighbor_scale_factor: 0.8,
            wall_scale_factor: 0.8,
            min_amplitude_ratio: 0.1,
        }
    }
}

impl AdaptiveParameter<f64, ChannelGenerationContext> for DistanceBasedAmplitudeAdapter {
    fn adapt(
        &self,
        base_value: f64,
        context: &ChannelGenerationContext,
    ) -> Result<f64, AdaptationError> {
        // Validate input parameters
        if base_value <= 0.0 || base_value.is_nan() || base_value.is_infinite() {
            return Err(AdaptationError::InvalidResult {
                value: base_value.to_string(),
                constraint: "Base amplitude must be positive and finite".to_string(),
            });
        }

        let mut scale_factor: f64 = 1.0;

        // Apply neighbor-based scaling
        if self.neighbor_scaling {
            if let Some(min_neighbor_dist) = context.min_neighbor_distance() {
                if min_neighbor_dist <= 0.0 {
                    return Err(AdaptationError::CalculationFailed {
                        parameter: "amplitude".to_string(),
                        reason: "Invalid neighbor distance".to_string(),
                    });
                }

                let neighbor_constraint = (min_neighbor_dist / 2.0) * self.neighbor_scale_factor;
                let neighbor_ratio = neighbor_constraint / base_value;
                scale_factor = scale_factor.min(neighbor_ratio);
            }
        }

        // Apply wall-based scaling
        if self.wall_scaling {
            let min_wall_dist = context.min_wall_distance();
            if min_wall_dist <= 0.0 {
                return Err(AdaptationError::CalculationFailed {
                    parameter: "amplitude".to_string(),
                    reason: "Invalid wall distance".to_string(),
                });
            }

            let wall_constraint = min_wall_dist * self.wall_scale_factor;
            let wall_ratio = wall_constraint / base_value;
            scale_factor = scale_factor.min(wall_ratio);
        }

        // Apply minimum ratio constraint
        scale_factor = scale_factor.max(self.min_amplitude_ratio);

        let result = base_value * scale_factor;

        // Validate result
        if result.is_nan() || result.is_infinite() || result <= 0.0 {
            return Err(AdaptationError::InvalidResult {
                value: result.to_string(),
                constraint: "Adapted amplitude must be positive and finite".to_string(),
            });
        }

        Ok(result)
    }

    fn is_adaptive(&self) -> bool {
        self.neighbor_scaling || self.wall_scaling
    }

    fn adaptation_description(&self) -> String {
        let mut parts = Vec::new();
        if self.neighbor_scaling {
            parts.push(format!("neighbor-aware ({}x)", self.neighbor_scale_factor));
        }
        if self.wall_scaling {
            parts.push(format!("wall-aware ({}x)", self.wall_scale_factor));
        }
        if parts.is_empty() {
            "no adaptation".to_string()
        } else {
            parts.join(", ")
        }
    }
}

/// Branch-count based adaptive behavior for density parameters
#[derive(Debug, Clone)]
pub struct BranchCountDensityAdapter {
    /// Base scaling exponent
    pub scaling_exponent: f64,

    /// Maximum scaling factor
    pub max_scale_factor: f64,

    /// Minimum scaling factor
    pub min_scale_factor: f64,
}

impl Default for BranchCountDensityAdapter {
    fn default() -> Self {
        Self {
            scaling_exponent: 0.75,
            max_scale_factor: 2.0,
            min_scale_factor: 0.5,
        }
    }
}

impl AdaptiveParameter<f64, ChannelGenerationContext> for BranchCountDensityAdapter {
    fn adapt(
        &self,
        base_value: f64,
        context: &ChannelGenerationContext,
    ) -> Result<f64, AdaptationError> {
        // Validate input parameters
        if base_value <= 0.0 || base_value.is_nan() || base_value.is_infinite() {
            return Err(AdaptationError::InvalidResult {
                value: base_value.to_string(),
                constraint: "Base value must be positive and finite".to_string(),
            });
        }

        if context.total_branches == 0 {
            return Err(AdaptationError::InvalidContext {
                reason: "Total branches must be greater than zero".to_string(),
            });
        }

        let branch_factor = (context.total_branches as f64)
            .powf(self.scaling_exponent)
            .max(1.0);

        if branch_factor.is_nan() || branch_factor.is_infinite() {
            return Err(AdaptationError::CalculationFailed {
                parameter: "branch_density".to_string(),
                reason: "Branch factor calculation resulted in NaN or infinite value".to_string(),
            });
        }

        let scale_factor = (1.0 / branch_factor)
            .max(self.min_scale_factor)
            .min(self.max_scale_factor);

        let result = base_value * scale_factor;

        // Validate result
        if result.is_nan() || result.is_infinite() || result <= 0.0 {
            return Err(AdaptationError::InvalidResult {
                value: result.to_string(),
                constraint: "Adapted value must be positive and finite".to_string(),
            });
        }

        Ok(result)
    }

    fn is_adaptive(&self) -> bool {
        true
    }

    fn adaptation_description(&self) -> String {
        format!(
            "branch-count scaling (exp: {}, range: {}-{})",
            self.scaling_exponent, self.min_scale_factor, self.max_scale_factor
        )
    }
}

/// Length-based adaptive behavior for wavelength parameters
#[derive(Debug, Clone)]
pub struct LengthBasedWavelengthAdapter {
    /// Target number of wavelengths per channel
    pub target_wavelengths: f64,

    /// Minimum wavelength factor
    pub min_wavelength_factor: f64,

    /// Maximum wavelength factor
    pub max_wavelength_factor: f64,
}

impl Default for LengthBasedWavelengthAdapter {
    fn default() -> Self {
        Self {
            target_wavelengths: 3.0,
            min_wavelength_factor: 0.5,
            max_wavelength_factor: 5.0,
        }
    }
}

impl AdaptiveParameter<f64, ChannelGenerationContext> for LengthBasedWavelengthAdapter {
    fn adapt(
        &self,
        base_value: f64,
        context: &ChannelGenerationContext,
    ) -> Result<f64, AdaptationError> {
        // Validate input parameters
        if base_value <= 0.0 || base_value.is_nan() || base_value.is_infinite() {
            return Err(AdaptationError::InvalidResult {
                value: base_value.to_string(),
                constraint: "Base wavelength factor must be positive and finite".to_string(),
            });
        }

        let channel_length = context.channel_length();
        if channel_length <= 0.0 {
            return Err(AdaptationError::InvalidContext {
                reason: "Channel length must be positive".to_string(),
            });
        }

        let channel_width = context.geometry_config.channel_width;
        if channel_width <= 0.0 {
            return Err(AdaptationError::InvalidContext {
                reason: "Channel width must be positive".to_string(),
            });
        }

        let base_wavelength = base_value * channel_width;

        if base_wavelength > 0.0 {
            let current_wavelengths = channel_length / base_wavelength;

            if current_wavelengths <= 0.0
                || current_wavelengths.is_nan()
                || current_wavelengths.is_infinite()
            {
                return Err(AdaptationError::CalculationFailed {
                    parameter: "wavelength".to_string(),
                    reason: "Invalid wavelength count calculation".to_string(),
                });
            }

            let adjustment_factor = self.target_wavelengths / current_wavelengths;

            if adjustment_factor.is_nan() || adjustment_factor.is_infinite() {
                return Err(AdaptationError::CalculationFailed {
                    parameter: "wavelength".to_string(),
                    reason: "Invalid adjustment factor calculation".to_string(),
                });
            }

            let adjusted_factor = (base_value * adjustment_factor)
                .max(self.min_wavelength_factor)
                .min(self.max_wavelength_factor);

            // Validate result
            if adjusted_factor.is_nan() || adjusted_factor.is_infinite() || adjusted_factor <= 0.0 {
                return Err(AdaptationError::InvalidResult {
                    value: adjusted_factor.to_string(),
                    constraint: "Adjusted wavelength factor must be positive and finite"
                        .to_string(),
                });
            }

            Ok(adjusted_factor)
        } else {
            Ok(base_value)
        }
    }

    fn is_adaptive(&self) -> bool {
        true
    }

    fn adaptation_description(&self) -> String {
        format!(
            "length-based wavelength (target: {} waves, range: {}-{})",
            self.target_wavelengths, self.min_wavelength_factor, self.max_wavelength_factor
        )
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
        let adapted = adapter.adapt(base_amplitude, &context).unwrap();

        // Should be scaled down due to constraints
        assert!(adapted < base_amplitude);
        assert!(adapted > 0.0);
    }

    #[test]
    fn test_branch_count_density_adapter() {
        let adapter = BranchCountDensityAdapter::default();
        let context = create_test_context();

        let base_density = 2.0;
        let adapted = adapter.adapt(base_density, &context).unwrap();

        // Should be scaled based on branch count
        assert!(adapted > 0.0);
        assert!(AdaptiveParameter::is_adaptive(&adapter));
    }

    #[test]
    fn test_length_based_wavelength_adapter() {
        let adapter = LengthBasedWavelengthAdapter::default();
        let context = create_test_context();

        let base_wavelength = 1.0;
        let adapted = adapter.adapt(base_wavelength, &context).unwrap();

        // Should adjust based on channel length
        assert!(adapted > 0.0);
        assert!(AdaptiveParameter::is_adaptive(&adapter));
    }
}

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
