//! Envelope calculation strategies for serpentine channel amplitude modulation.

use crate::config::ConstantsRegistry;

/// Trait for calculating envelope functions for serpentine channels
///
/// This trait abstracts envelope calculation logic to eliminate code duplication
/// and provide a clean interface for different envelope types.
pub trait EnvelopeCalculator {
    /// Calculate envelope value at parameter t (0.0 to 1.0)
    fn calculate_envelope(&self, t: f64, context: &EnvelopeContext) -> f64;
}

/// Context for envelope calculations
#[derive(Debug, Clone)]
pub struct EnvelopeContext {
    /// Channel length
    pub channel_length: f64,
    /// Channel direction vector (dx, dy)
    pub direction: (f64, f64),
    /// Distance to nearest node
    pub node_distance: f64,
    /// Adaptive configuration
    pub adaptive_config: crate::config::AdaptiveSerpentineConfig,
    /// Gaussian width factor
    pub gaussian_width_factor: f64,
}

/// Smooth endpoint envelope calculator
///
/// Provides smooth transitions at channel endpoints using smoothstep function.
pub struct SmoothEndpointEnvelopeCalculator;

impl EnvelopeCalculator for SmoothEndpointEnvelopeCalculator {
    fn calculate_envelope(&self, t: f64, _context: &EnvelopeContext) -> f64 {
        // Smooth endpoint envelope using smoothstep function
        // This ensures zero amplitude and zero derivative at endpoints
        let smoothstep = |x: f64| x * x * 2.0f64.mul_add(-x, 3.0);

        let constants = ConstantsRegistry::new();
        let start_threshold = constants.get_smooth_endpoint_start_threshold();
        let end_threshold = constants.get_smooth_endpoint_end_threshold();

        // Create smooth transitions at both ends
        let start_transition = if t < start_threshold {
            smoothstep(t / start_threshold)
        } else {
            1.0
        };

        let end_transition = if t > end_threshold {
            smoothstep((1.0 - t) / (1.0 - end_threshold))
        } else {
            1.0
        };

        start_transition * end_transition
    }
}

/// Gaussian envelope calculator with adaptive behavior
///
/// Provides Gaussian-shaped envelope with adaptive parameters based on
/// channel characteristics and proximity to nodes/walls.
pub struct AdaptiveGaussianEnvelopeCalculator;

impl EnvelopeCalculator for AdaptiveGaussianEnvelopeCalculator {
    fn calculate_envelope(&self, t: f64, context: &EnvelopeContext) -> f64 {
        let dx = context.direction.0;
        let dy = context.direction.1;
        let channel_length = context.channel_length;
        let node_distance = context.node_distance;

        // Determine if this is primarily a horizontal channel
        let is_horizontal = dx.abs() > dy.abs();
        let horizontal_ratio = dx.abs() / node_distance;

        // For horizontal channels (middle sections), we want less aggressive tapering
        let middle_section_factor = if is_horizontal
            && horizontal_ratio > context.adaptive_config.horizontal_ratio_threshold
        {
            (1.0 - context.adaptive_config.middle_section_amplitude_factor).mul_add(
                horizontal_ratio,
                context.adaptive_config.middle_section_amplitude_factor,
            )
        } else {
            1.0
        };

        // Distance-based normalization
        let distance_normalization = if context.adaptive_config.enable_distance_based_scaling {
            (node_distance / context.adaptive_config.node_distance_normalization)
                .clamp(0.1, 1.0)
        } else {
            1.0
        };

        // Calculate effective sigma based on distance and section type
        let base_sigma = channel_length / context.gaussian_width_factor;
        let effective_sigma = base_sigma * distance_normalization * middle_section_factor;

        // Center the envelope
        let center = 0.5;

        // Create smooth dome-shaped envelope instead of sharp Gaussian peaks
        // This prevents self-intersection when channels curve back on themselves
        let dome_envelope = if (t - center).abs() < 0.45 {
            // Use raised cosine for the main dome (much smoother than Gaussian)
            let normalized_t = (t - center) / 0.45; // Scale to [-1, 1] range
            let cosine_factor = 0.5 * (1.0 + (std::f64::consts::PI * normalized_t).cos());

            // Apply effective sigma scaling to the dome
            let sigma_factor = (effective_sigma / channel_length).clamp(0.1, 0.3);
            let dome_width = 0.45 * sigma_factor / 0.2; // Scale dome width based on sigma

            if (t - center).abs() < dome_width {
                let dome_t = (t - center) / dome_width;
                0.5 * (1.0 + (std::f64::consts::PI * dome_t).cos())
            } else {
                // Smooth transition to zero
                let transition_factor = ((t - center).abs() - dome_width) / (0.45 - dome_width);
                let smoothstep = (transition_factor * transition_factor)
                    .mul_add(-2.0f64.mul_add(-transition_factor, 3.0), 1.0);
                cosine_factor * smoothstep * 0.1
            }
        } else {
            // Smooth transition to zero at edges using smoothstep
            let edge_distance = ((t - center).abs() - 0.45) / 0.05; // 0.05 is transition zone
            if (0.0..1.0).contains(&edge_distance) {
                let smoothstep = (edge_distance * edge_distance)
                    .mul_add(-2.0f64.mul_add(-edge_distance, 3.0), 1.0);
                smoothstep * 0.05 // Very small amplitude at edges
            } else {
                0.0
            }
        };

        // For middle sections, enhance the dome but keep it smooth
        if is_horizontal && horizontal_ratio > context.adaptive_config.horizontal_ratio_threshold {
            let plateau_width = context.adaptive_config.plateau_width_factor.min(0.3); // Limit plateau width
            let plateau_start = 0.5 - plateau_width / 2.0;
            let plateau_end = 0.5 + plateau_width / 2.0;

            if t >= plateau_start && t <= plateau_end {
                // In the plateau region, enhance the dome but keep it smooth
                let plateau_factor = 1.0 - ((t - 0.5).abs() / (plateau_width / 2.0));
                let enhanced_amplitude = (1.0 - context.adaptive_config.plateau_amplitude_factor)
                    .mul_add(
                        plateau_factor,
                        context.adaptive_config.plateau_amplitude_factor,
                    );
                dome_envelope.max(enhanced_amplitude * dome_envelope)
            } else {
                dome_envelope
            }
        } else {
            dome_envelope
        }
    }
}
