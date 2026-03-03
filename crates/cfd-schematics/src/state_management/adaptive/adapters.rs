//! Base adaptive parameter implementations
//!
//! Distance-based amplitude adaptation, branch-count density scaling,
//! and length-based wavelength adaptation.

use super::{AdaptationError, AdaptiveParameter, ChannelGenerationContext};

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
