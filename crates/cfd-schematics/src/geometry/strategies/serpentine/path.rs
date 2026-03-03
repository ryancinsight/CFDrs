//! Serpentine path generation and geometric safety.

use super::super::envelope::{
    AdaptiveGaussianEnvelopeCalculator, EnvelopeCalculator, EnvelopeContext,
};
use super::super::ChannelGenerationContext;
use super::SerpentineChannelStrategy;
use crate::config::{ConstantsRegistry, SerpentineConfig};
use crate::geometry::optimization::optimize_serpentine_parameters;
use crate::geometry::Point2D;

impl SerpentineChannelStrategy {
    /// Estimate local bend radius from three consecutive points.
    fn estimate_local_bend_radius(a: Point2D, b: Point2D, c: Point2D) -> Option<f64> {
        let ab = ((b.0 - a.0).powi(2) + (b.1 - a.1).powi(2)).sqrt();
        let bc = ((c.0 - b.0).powi(2) + (c.1 - b.1).powi(2)).sqrt();
        let ac = ((c.0 - a.0).powi(2) + (c.1 - a.1).powi(2)).sqrt();
        if ab <= 1e-9 || bc <= 1e-9 || ac <= 1e-9 {
            return None;
        }

        let twice_area = ((b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0)).abs();
        if twice_area <= 1e-12 {
            return Some(f64::INFINITY);
        }

        Some((ab * bc * ac) / (2.0 * twice_area))
    }

    /// Smooth path points when local radius falls below a target threshold.
    fn enforce_min_bend_radius(&self, path: &mut [Point2D], min_radius: f64, iterations: usize) {
        if path.len() < 3 || min_radius <= 0.0 {
            return;
        }

        for _ in 0..iterations {
            let mut changed = false;
            for i in 1..path.len() - 1 {
                let a = path[i - 1];
                let b = path[i];
                let c = path[i + 1];

                if let Some(radius) = Self::estimate_local_bend_radius(a, b, c) {
                    if radius.is_finite() && radius < min_radius {
                        let midpoint = ((a.0 + c.0) * 0.5, (a.1 + c.1) * 0.5);
                        let blend = (1.0 - radius / min_radius).clamp(0.15, 0.65);
                        path[i] = (
                            b.0 + (midpoint.0 - b.0) * blend,
                            b.1 + (midpoint.1 - b.1) * blend,
                        );
                        changed = true;
                    }
                }
            }

            if !changed {
                break;
            }
        }
    }

    /// Generate a serpentine path between two points using zero-copy techniques
    pub(super) fn generate_serpentine_path(
        &self,
        p1: Point2D,
        p2: Point2D,
        context: &ChannelGenerationContext,
    ) -> Vec<Point2D> {
        let initial_wavelength =
            self.config.wavelength_factor * context.geometry_config.channel_width;
        let wavelength = self.validate_wavelength_for_diameter(
            initial_wavelength,
            context.geometry_config.channel_width,
        );
        let amplitude = self.calculate_adaptive_amplitude(p1, p2, context, wavelength);

        if amplitude <= 0.0 {
            return self.generate_straight_line_path(
                p1,
                p2,
                context.geometry_config.generation.serpentine_points,
            );
        }

        let dx = p2.0 - p1.0;
        let dy = p2.1 - p1.1;
        let channel_length = dx.hypot(dy);
        let _angle = dy.atan2(dx);

        let constants = ConstantsRegistry::new();
        let _branch_factor = (context.total_branches as f64)
            .powf(constants.get_branch_factor_exponent())
            .max(1.0);

        let initial_wavelength =
            self.config.wavelength_factor * context.geometry_config.channel_width;
        let base_wavelength = self.validate_wavelength_for_diameter(
            initial_wavelength,
            context.geometry_config.channel_width,
        );

        let initial_amplitude = self.calculate_adaptive_amplitude(p1, p2, context, base_wavelength);
        let start_guard = self.calculate_junction_guard_fraction(p1.1, channel_length, context);
        let end_guard = self.calculate_junction_guard_fraction(p2.1, channel_length, context);
        let channel_diameter = context.geometry_config.channel_width;
        let effective_wavelength = self.calculate_effective_wavelength(
            base_wavelength,
            initial_amplitude,
            channel_diameter,
        );
        let square_sharpness =
            self.calculate_effective_square_sharpness(channel_diameter, effective_wavelength);

        let length_based_periods =
            (channel_length / effective_wavelength) * self.config.wave_density_factor;
        let base_periods = length_based_periods.max(1.0);
        let requested_half_periods = (base_periods * 2.0).max(1.0);

        let base_points = context.geometry_config.generation.serpentine_points;
        let n_points = self.calculate_required_wave_points(requested_half_periods, base_points);
        let min_points_per_half_period =
            if matches!(self.config.wave_shape, crate::config::WaveShape::Square) {
                6.0
            } else {
                12.0
            };
        let max_resolvable_half_periods =
            ((n_points - 1) as f64 / min_points_per_half_period).max(1.0);
        let bounded_half_periods = requested_half_periods.min(max_resolvable_half_periods);
        let mut half_periods = bounded_half_periods.round().max(1.0);
        if !matches!(self.config.wave_shape, crate::config::WaveShape::Square)
            && half_periods < 2.0
            && max_resolvable_half_periods >= 2.0
        {
            half_periods = 2.0;
        }
        if !matches!(self.config.wave_shape, crate::config::WaveShape::Square) {
            let central_window = (1.0 - start_guard - end_guard).max(0.0);
            if central_window < 0.35 && max_resolvable_half_periods >= 4.0 {
                half_periods = half_periods.max(4.0).min(max_resolvable_half_periods);
            }
        }

        let phase_direction = self.calculate_wave_phase_direction(p1, p2, context);

        let mut path = Vec::with_capacity(n_points);
        let perp_x = -dy / channel_length;
        let perp_y = dx / channel_length;
        let smoothstep = |v: f64| {
            let x = v.clamp(0.0, 1.0);
            x * x * 2.0f64.mul_add(-x, 3.0)
        };

        for i in 0..n_points {
            let t = i as f64 / (n_points - 1) as f64;

            let base_x = t.mul_add(dx, p1.0);
            let base_y = t.mul_add(dy, p1.1);

            let node_distance = dx.hypot(dy);
            let envelope_context = EnvelopeContext {
                channel_length,
                direction: (dx, dy),
                node_distance,
                adaptive_config: self.config.adaptive_config,
                gaussian_width_factor: self.config.gaussian_width_factor,
            };
            let improved_envelope_calc = AdaptiveGaussianEnvelopeCalculator;
            let envelope = improved_envelope_calc.calculate_envelope(t, &envelope_context);

            let start_guard_scale = if t < start_guard {
                let x = smoothstep(t / start_guard);
                x * x * x
            } else {
                1.0
            };
            let end_guard_scale = if t > (1.0 - end_guard) {
                let x = smoothstep((1.0 - t) / end_guard);
                x * x * x
            } else {
                1.0
            };
            let guarded_envelope = envelope * start_guard_scale * end_guard_scale;

            let wave_phase = std::f64::consts::PI * half_periods * t;

            let phase_offset = if phase_direction < 0.0 {
                std::f64::consts::PI
            } else {
                0.0
            };

            let local_cap = self.calculate_local_amplitude_cap(base_y, perp_y.abs(), context);
            let mut capped_amplitude = (initial_amplitude * guarded_envelope).min(local_cap);

            if t < start_guard || t > (1.0 - end_guard) {
                let distance_to_nearest_endpoint = t.min(1.0 - t) * channel_length;
                let endpoint_distance_cap = distance_to_nearest_endpoint * 0.25;
                capped_amplitude = capped_amplitude.min(endpoint_distance_cap);
            }
            let wave_amplitude = capped_amplitude
                * self.calculate_wave_amplitude(wave_phase, phase_offset, square_sharpness);

            let x = wave_amplitude.mul_add(perp_x, base_x);
            let y = wave_amplitude.mul_add(perp_y, base_y);

            if i == 0 {
                path.push(p1);
            } else if i == n_points - 1 {
                path.push(p2);
            } else {
                path.push((x, y));
            }
        }

        let min_radius = if matches!(self.config.wave_shape, crate::config::WaveShape::Square) {
            channel_diameter
        } else {
            channel_diameter * 0.5
        };
        self.enforce_min_bend_radius(&mut path, min_radius, 4);

        path
    }

    /// Generate an optimized serpentine path between two points
    pub(super) fn generate_optimized_serpentine_path(
        &self,
        p1: Point2D,
        p2: Point2D,
        context: &ChannelGenerationContext,
    ) -> Vec<Point2D> {
        let initial_wavelength =
            self.config.wavelength_factor * context.geometry_config.channel_width;
        let wavelength = self.validate_wavelength_for_diameter(
            initial_wavelength,
            context.geometry_config.channel_width,
        );
        let amplitude = self.calculate_adaptive_amplitude(p1, p2, context, wavelength);

        if amplitude <= 0.0 {
            return self.generate_straight_line_path(
                p1,
                p2,
                context.geometry_config.generation.serpentine_points,
            );
        }

        let optimization_result = optimize_serpentine_parameters(
            p1,
            p2,
            context.geometry_config,
            &self.config,
            context.box_dims,
            context.neighbor_info,
        );

        let optimized_config = SerpentineConfig {
            wavelength_factor: optimization_result.params.wavelength_factor,
            wave_density_factor: optimization_result.params.wave_density_factor,
            fill_factor: optimization_result.params.fill_factor,
            gaussian_width_factor: self.config.gaussian_width_factor,
            wave_phase_direction: self.config.wave_phase_direction,
            wave_shape: self.config.wave_shape,
            optimization_enabled: false,
            target_fill_ratio: self.config.target_fill_ratio,
            optimization_profile: self.config.optimization_profile,
            adaptive_config: self.config.adaptive_config,
        };

        let temp_strategy = Self::new(optimized_config);
        temp_strategy.generate_serpentine_path(p1, p2, context)
    }

    /// Generate serpentine path for optimization purposes (public interface)
    #[must_use]
    pub fn generate_serpentine_path_for_optimization(
        &self,
        p1: Point2D,
        p2: Point2D,
        geometry_config: &crate::config::GeometryConfig,
        box_dims: (f64, f64),
        total_branches: usize,
        neighbor_info: Option<&[f64]>,
    ) -> Vec<Point2D> {
        let context =
            ChannelGenerationContext::new(geometry_config, box_dims, total_branches, neighbor_info);
        self.generate_serpentine_path(p1, p2, &context)
    }

    /// Generate a straight line path when serpentine amplitude is too small
    fn generate_straight_line_path(
        &self,
        p1: Point2D,
        p2: Point2D,
        n_points: usize,
    ) -> Vec<Point2D> {
        let dx = p2.0 - p1.0;
        let dy = p2.1 - p1.1;

        (0..n_points)
            .map(|i| {
                let t = i as f64 / (n_points - 1) as f64;
                (t.mul_add(dx, p1.0), t.mul_add(dy, p1.1))
            })
            .collect()
    }
}
