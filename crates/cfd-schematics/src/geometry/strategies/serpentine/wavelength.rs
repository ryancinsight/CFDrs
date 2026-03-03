//! Wavelength calculation and validation for serpentine channels.

use super::super::ChannelGenerationContext;
use super::SerpentineChannelStrategy;
use crate::config::{constants, ConstantsRegistry};

impl SerpentineChannelStrategy {
    /// Calculate wavelength-aware scaling factor with adaptive thresholds
    pub(super) fn calculate_wavelength_scaling_factor(
        &self,
        wavelength: f64,
        channel_width: f64,
    ) -> f64 {
        let min_separation = channel_width;
        let wavelength_ratio = wavelength / min_separation;

        if wavelength_ratio >= 3.0 {
            1.0
        } else if wavelength_ratio >= 2.0 {
            0.95
        } else if wavelength_ratio >= 1.5 {
            0.85
        } else if wavelength_ratio >= 1.2 {
            0.75
        } else {
            0.65
        }
    }

    /// Calculate density enhancement factor for better space utilization
    pub(super) fn calculate_density_enhancement_factor(
        &self,
        context: &ChannelGenerationContext,
    ) -> f64 {
        let box_area = context.box_dims.0 * context.box_dims.1;
        let branch_density = context.total_branches as f64 / box_area;

        if branch_density < 0.01 {
            1.08
        } else if branch_density < 0.02 {
            1.05
        } else if branch_density < 0.05 {
            1.02
        } else {
            1.0
        }
    }

    /// Validate and adjust wavelength for manufacturing constraints
    pub(super) fn validate_wavelength_for_diameter(
        &self,
        wavelength: f64,
        channel_width: f64,
    ) -> f64 {
        let min_separation = channel_width;
        let min_wavelength = min_separation * 3.0;
        wavelength.max(min_wavelength)
    }

    /// Calculate minimum wavelength from diameter-aware curvature bounds.
    fn calculate_curvature_safe_wavelength(&self, amplitude: f64, channel_diameter: f64) -> f64 {
        let safe_amplitude = amplitude.max(channel_diameter * 0.25);
        let min_radius_multiplier =
            if matches!(self.config.wave_shape, crate::config::WaveShape::Square) {
                1.5
            } else {
                1.0
            };
        let min_bend_radius = channel_diameter * min_radius_multiplier;

        // For y=A*sin(2πx/λ), curvature max at peaks is κ_max = A*(2π/λ)^2.
        // Enforcing R_min = 1/κ_max yields λ >= 2π*sqrt(A*R_min).
        2.0 * std::f64::consts::PI * (safe_amplitude * min_bend_radius).sqrt()
    }

    /// Adjust effective wavelength for smoother non-square waves with arc-like turns.
    pub(super) fn calculate_effective_wavelength(
        &self,
        base_wavelength: f64,
        amplitude: f64,
        channel_diameter: f64,
    ) -> f64 {
        const MAX_SINE_SLOPE: f64 = 50.0;
        const ARC_SPAN_MULTIPLIER: f64 = 0.15;

        let slope_limited_wavelength = if amplitude > 1e-6 {
            2.0 * std::f64::consts::PI * amplitude / MAX_SINE_SLOPE
        } else {
            base_wavelength
        };
        let arc_span_wavelength = amplitude * ARC_SPAN_MULTIPLIER;
        let curvature_safe = self.calculate_curvature_safe_wavelength(amplitude, channel_diameter);

        base_wavelength
            .max(slope_limited_wavelength)
            .max(arc_span_wavelength)
            .max(curvature_safe)
    }

    /// Compute adaptive square-wave sharpness from wavelength/diameter ratio.
    pub(super) fn calculate_effective_square_sharpness(
        &self,
        channel_diameter: f64,
        effective_wavelength: f64,
    ) -> f64 {
        let constants = ConstantsRegistry::new();
        let base_sharpness = constants.get_square_wave_sharpness();
        let ratio = effective_wavelength / channel_diameter.max(1e-6);
        let max_diameter_safe_sharpness = base_sharpness.min(3.0);
        let capped_sharpness = (ratio * 0.2).clamp(1.2, max_diameter_safe_sharpness);
        if matches!(self.config.wave_shape, crate::config::WaveShape::Square) {
            capped_sharpness
        } else {
            base_sharpness
        }
    }

    /// Compute adaptive path sample count from requested wave count.
    pub(super) fn calculate_required_wave_points(
        &self,
        half_periods: f64,
        base_points: usize,
    ) -> usize {
        let min_points_per_half_period =
            if matches!(self.config.wave_shape, crate::config::WaveShape::Square) {
                6.0
            } else {
                12.0
            };
        let required_points = (half_periods * min_points_per_half_period).ceil() as usize + 1;

        base_points.max(required_points).clamp(
            constants::MIN_SERPENTINE_POINTS,
            constants::MAX_SERPENTINE_POINTS,
        )
    }
}
