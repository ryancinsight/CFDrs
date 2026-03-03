//! Amplitude calculation and space metrics for serpentine channels.

use super::super::ChannelGenerationContext;
use super::SerpentineChannelStrategy;
use super::SpaceMetrics;
use crate::config::ConstantsRegistry;
use crate::geometry::Point2D;

impl SerpentineChannelStrategy {
    /// Calculate maximum safe amplitude using advanced adaptive algorithms
    pub(super) fn calculate_adaptive_amplitude(
        &self,
        p1: Point2D,
        p2: Point2D,
        context: &ChannelGenerationContext,
        wavelength: f64,
    ) -> f64 {
        let constants = ConstantsRegistry::new();
        let channel_width = context.geometry_config.channel_width;

        let space_metrics = self.analyze_space_metrics(p1, p2, context);
        if space_metrics.available_space <= constants.get_geometric_tolerance() {
            return 0.0;
        }

        let base_amplitude = space_metrics.available_space;
        let wavelength_factor = self.calculate_wavelength_scaling_factor(wavelength, channel_width);
        let density_factor = self.calculate_density_enhancement_factor(context);
        let mut amplitude =
            base_amplitude * wavelength_factor * density_factor * self.config.fill_factor;
        if !matches!(self.config.wave_shape, crate::config::WaveShape::Square) {
            amplitude *= 1.08;
        }

        amplitude = amplitude.min(base_amplitude * 0.99);

        let min_wave_amplitude = channel_width * 0.15;
        if amplitude < min_wave_amplitude {
            amplitude = min_wave_amplitude.min(base_amplitude * 0.95);
        }

        if amplitude <= constants.get_geometric_tolerance() {
            0.0
        } else {
            amplitude
        }
    }

    /// Space metrics analysis for amplitude calculation with dynamic space utilization
    fn analyze_space_metrics(
        &self,
        p1: Point2D,
        p2: Point2D,
        context: &ChannelGenerationContext,
    ) -> SpaceMetrics {
        let channel_center_y = f64::midpoint(p1.1, p2.1);
        let box_height = context.box_dims.1;
        let wall_clearance = context.geometry_config.wall_clearance;
        let channel_width = context.geometry_config.channel_width;

        let available_space = if self.config.adaptive_config.enable_neighbor_avoidance {
            if let Some(neighbor_info) = context.neighbor_info {
                if neighbor_info.is_empty() {
                    let safety_margin = wall_clearance + channel_width * 0.5;
                    let space_above = box_height - channel_center_y - safety_margin;
                    let space_below = channel_center_y - safety_margin;
                    space_above.min(space_below)
                } else {
                    self.calculate_dynamic_available_space(
                        channel_center_y,
                        box_height,
                        wall_clearance,
                        channel_width,
                        neighbor_info,
                    )
                }
            } else {
                let safety_margin = wall_clearance + channel_width * 0.5;
                let space_above = box_height - channel_center_y - safety_margin;
                let space_below = channel_center_y - safety_margin;
                space_above.min(space_below)
            }
        } else {
            let safety_margin = wall_clearance + channel_width * 0.5;
            let space_above = box_height - channel_center_y - safety_margin;
            let space_below = channel_center_y - safety_margin;
            space_above.min(space_below)
        };

        SpaceMetrics {
            available_space: available_space.max(0.0),
        }
    }

    /// Calculate dynamic available space with proper geometric constraints
    fn calculate_dynamic_available_space(
        &self,
        channel_center_y: f64,
        box_height: f64,
        wall_clearance: f64,
        channel_width: f64,
        neighbor_info: &[f64],
    ) -> f64 {
        let mut closest_above = box_height;
        let mut closest_below: f64 = 0.0;

        for &neighbor_y in neighbor_info {
            if (neighbor_y - channel_center_y).abs() > 0.1 {
                if neighbor_y > channel_center_y {
                    closest_above = closest_above.min(neighbor_y);
                } else {
                    closest_below = closest_below.max(neighbor_y);
                }
            }
        }

        let distance_to_neighbor_above = closest_above - channel_center_y;
        let distance_to_neighbor_below = channel_center_y - closest_below;
        let min_neighbor_distance = distance_to_neighbor_above.min(distance_to_neighbor_below);

        let wall_margin = wall_clearance + channel_width * 0.5;
        let distance_to_top_wall = box_height - channel_center_y - wall_margin;
        let distance_to_bottom_wall = channel_center_y - wall_margin;
        let min_wall_distance = distance_to_top_wall.min(distance_to_bottom_wall);

        let neighbor_constraint = if min_neighbor_distance < box_height {
            let branch_safety = channel_width;
            (min_neighbor_distance - branch_safety) * 0.5
        } else {
            f64::INFINITY
        };

        let wall_constraint = min_wall_distance;

        neighbor_constraint.min(wall_constraint).max(0.0)
    }

    /// Find the minimum centerline distance from a given y-position to neighbor centerlines.
    pub(super) fn min_neighbor_distance_at_y(
        &self,
        y: f64,
        neighbors: Option<&[f64]>,
    ) -> Option<f64> {
        neighbors.and_then(|neighbor_values| {
            neighbor_values
                .iter()
                .map(|neighbor_y| (neighbor_y - y).abs())
                .filter(|distance| *distance > 0.1)
                .min_by(f64::total_cmp)
        })
    }

    /// Compute endpoint guard length (fraction of channel length) based on local neighbor spacing.
    pub(super) fn calculate_junction_guard_fraction(
        &self,
        endpoint_y: f64,
        channel_length: f64,
        context: &ChannelGenerationContext,
    ) -> f64 {
        const MIN_GUARD_FRACTION: f64 = 0.10;
        const MAX_GUARD_FRACTION: f64 = 0.35;

        let mut guard_fraction = MIN_GUARD_FRACTION;

        let channel_diameter = context.geometry_config.channel_width;
        let physical_guard_length =
            context.geometry_config.wall_clearance + channel_diameter * 20.0;
        let physical_guard_fraction = if channel_length > 1e-6 {
            (physical_guard_length / channel_length).clamp(MIN_GUARD_FRACTION, 0.45)
        } else {
            0.45
        };
        guard_fraction = guard_fraction.max(physical_guard_fraction);

        if !self.config.adaptive_config.enable_neighbor_avoidance {
            return guard_fraction.min(0.45);
        }

        if let Some(min_neighbor_distance) =
            self.min_neighbor_distance_at_y(endpoint_y, context.neighbor_info)
        {
            let required_distance = channel_diameter * 1.2;
            if min_neighbor_distance <= required_distance {
                return MAX_GUARD_FRACTION.max(guard_fraction).min(0.45);
            }

            let tightness = (required_distance / min_neighbor_distance).clamp(0.0, 1.0);
            let neighbor_guard =
                (MAX_GUARD_FRACTION - MIN_GUARD_FRACTION).mul_add(tightness, MIN_GUARD_FRACTION);
            guard_fraction = guard_fraction.max(neighbor_guard);
        }

        guard_fraction.min(0.45)
    }

    /// Calculate a local amplitude cap at a given base point using wall and neighbor clearances.
    pub(super) fn calculate_local_amplitude_cap(
        &self,
        base_y: f64,
        perp_y_abs: f64,
        context: &ChannelGenerationContext,
    ) -> f64 {
        let normal_y_factor = perp_y_abs.max(0.1);
        let channel_diameter = context.geometry_config.channel_width;
        let mut amplitude_cap = f64::INFINITY;

        if self.config.adaptive_config.enable_wall_proximity_scaling {
            let wall_margin = context.geometry_config.wall_clearance + channel_diameter * 0.5;
            let wall_available =
                (base_y - wall_margin).min(context.box_dims.1 - base_y - wall_margin);
            amplitude_cap = amplitude_cap.min((wall_available / normal_y_factor).max(0.0));
        }

        if self.config.adaptive_config.enable_neighbor_avoidance {
            if let Some(min_neighbor_distance) =
                self.min_neighbor_distance_at_y(base_y, context.neighbor_info)
            {
                let neighbor_available =
                    ((min_neighbor_distance - channel_diameter) * 0.5).max(0.0);
                amplitude_cap = amplitude_cap.min(neighbor_available / normal_y_factor);
            }
        }

        amplitude_cap
    }
}
