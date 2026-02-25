//! Serpentine channel generation strategy.

use crate::config::{constants, ConstantsRegistry, GeometryConfig, SerpentineConfig};
use crate::geometry::optimization::optimize_serpentine_parameters;
use crate::geometry::{ChannelType, Point2D};
use crate::state_management::bilateral_symmetry::{
    BilateralPhaseDirectionCalculator, BilateralSymmetryConfig, SymmetryContext,
};

use super::envelope::{AdaptiveGaussianEnvelopeCalculator, EnvelopeCalculator, EnvelopeContext};
use super::{ChannelGenerationContext, ChannelTypeStrategy};

/// Space metrics for amplitude calculation
#[derive(Debug, Clone)]
struct SpaceMetrics {
    /// Available space for amplitude expansion
    available_space: f64,
}

/// Strategy for creating serpentine channels
#[derive(Debug, Clone)]
pub struct SerpentineChannelStrategy {
    config: SerpentineConfig,
}

impl SerpentineChannelStrategy {
    /// Create a new serpentine channel strategy with the given configuration
    ///
    /// # Arguments
    ///
    /// * `config` - Configuration parameters for serpentine channel generation
    ///
    /// # Examples
    ///
    /// ```rust
    /// use cfd_schematics::geometry::strategies::SerpentineChannelStrategy;
    /// use cfd_schematics::config::SerpentineConfig;
    ///
    /// let strategy = SerpentineChannelStrategy::new(SerpentineConfig::default());
    /// ```
    #[must_use]
    pub const fn new(config: SerpentineConfig) -> Self {
        Self { config }
    }
}

impl ChannelTypeStrategy for SerpentineChannelStrategy {
    fn create_channel(
        &self,
        from: Point2D,
        to: Point2D,
        geometry_config: &GeometryConfig,
        box_dims: (f64, f64),
        total_branches: usize,
        neighbor_info: Option<&[f64]>,
    ) -> ChannelType {
        let context =
            ChannelGenerationContext::new(geometry_config, box_dims, total_branches, neighbor_info);

        let path = if self.config.optimization_enabled {
            self.generate_optimized_serpentine_path(from, to, &context)
        } else {
            self.generate_serpentine_path(from, to, &context)
        };
        ChannelType::Serpentine { path }
    }
}

impl SerpentineChannelStrategy {
    /// Calculate wave amplitude based on wave shape and phase.
    fn calculate_wave_amplitude(
        &self,
        wave_phase: f64,
        phase_offset: f64,
        square_sharpness: f64,
    ) -> f64 {
        use crate::config::WaveShape;

        match self.config.wave_shape {
            WaveShape::Sine => {
                // Smooth sine wave
                (wave_phase + phase_offset).sin()
            }
            WaveShape::Square => {
                // Square wave with smooth transitions
                let sine_value = (wave_phase + phase_offset).sin();
                // Use tanh to create smooth square wave transitions.
                (square_sharpness * sine_value).tanh()
            }
        }
    }

    /// Calculate maximum safe amplitude using advanced adaptive algorithms
    fn calculate_adaptive_amplitude(
        &self,
        p1: Point2D,
        p2: Point2D,
        context: &ChannelGenerationContext,
        wavelength: f64,
    ) -> f64 {
        let constants = ConstantsRegistry::new();
        let channel_width = context.geometry_config.channel_width;

        // Dynamic space analysis with wall/neighbor constraints.
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
            // Non-square waves should better utilize the available lane height while
            // local caps still enforce wall/neighbor safety.
            amplitude *= 1.08;
        }

        // Clamp strictly to local available space so we do not violate branch/wall spacing.
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

        // Calculate available space considering neighbors and walls.
        // Neighbor-based reduction can be disabled for more aggressive wave amplitudes.
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
                // Single channel - use full box space minus wall clearances
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
        // Find the closest neighbors above and below
        let mut closest_above = box_height;
        let mut closest_below: f64 = 0.0;

        for &neighbor_y in neighbor_info {
            if (neighbor_y - channel_center_y).abs() > 0.1 {
                // Exclude self (within 0.1mm tolerance)
                if neighbor_y > channel_center_y {
                    closest_above = closest_above.min(neighbor_y);
                } else {
                    closest_below = closest_below.max(neighbor_y);
                }
            }
        }

        // Calculate distance to nearest neighbor
        let distance_to_neighbor_above = closest_above - channel_center_y;
        let distance_to_neighbor_below = channel_center_y - closest_below;
        let min_neighbor_distance = distance_to_neighbor_above.min(distance_to_neighbor_below);

        // Calculate distance to walls
        let wall_margin = wall_clearance + channel_width * 0.5;
        let distance_to_top_wall = box_height - channel_center_y - wall_margin;
        let distance_to_bottom_wall = channel_center_y - wall_margin;
        let min_wall_distance = distance_to_top_wall.min(distance_to_bottom_wall);

        // Available space is constrained by both neighbors and walls
        let neighbor_constraint = if min_neighbor_distance < box_height {
            let branch_safety = channel_width;
            (min_neighbor_distance - branch_safety) * 0.5
        } else {
            f64::INFINITY // No neighbors
        };

        let wall_constraint = min_wall_distance;

        // Use the most restrictive constraint
        neighbor_constraint.min(wall_constraint).max(0.0)
    }

    /// Find the minimum centerline distance from a given y-position to neighbor centerlines.
    fn min_neighbor_distance_at_y(&self, y: f64, neighbors: Option<&[f64]>) -> Option<f64> {
        neighbors.and_then(|neighbor_values| {
            neighbor_values
                .iter()
                .map(|neighbor_y| (neighbor_y - y).abs())
                .filter(|distance| *distance > 0.1)
                .min_by(f64::total_cmp)
        })
    }

    /// Compute endpoint guard length (fraction of channel length) based on local neighbor spacing.
    fn calculate_junction_guard_fraction(
        &self,
        endpoint_y: f64,
        channel_length: f64,
        context: &ChannelGenerationContext,
    ) -> f64 {
        const MIN_GUARD_FRACTION: f64 = 0.10;
        const MAX_GUARD_FRACTION: f64 = 0.35;

        let mut guard_fraction = MIN_GUARD_FRACTION;

        // Always enforce a physical guard zone around nodes for short segments
        // to prevent early wave bulging into split/merge connector geometry.
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
    fn calculate_local_amplitude_cap(
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

    /// Calculate wavelength-aware scaling factor with adaptive thresholds
    fn calculate_wavelength_scaling_factor(&self, wavelength: f64, channel_width: f64) -> f64 {
        let min_separation = channel_width;
        let wavelength_ratio = wavelength / min_separation;

        // Adaptive scaling based on wavelength ratio
        if wavelength_ratio >= 3.0 {
            1.0 // Full utilization for large wavelengths
        } else if wavelength_ratio >= 2.0 {
            0.95 // Near-full utilization
        } else if wavelength_ratio >= 1.5 {
            0.85 // Good utilization
        } else if wavelength_ratio >= 1.2 {
            0.75 // Moderate utilization
        } else {
            0.65 // Conservative but still aggressive
        }
    }

    /// Calculate density enhancement factor for better space utilization
    fn calculate_density_enhancement_factor(&self, context: &ChannelGenerationContext) -> f64 {
        let box_area = context.box_dims.0 * context.box_dims.1;
        let branch_density = context.total_branches as f64 / box_area;

        // More conservative enhancement for lower density layouts
        if branch_density < 0.01 {
            1.08 // 8% boost for sparse layouts
        } else if branch_density < 0.02 {
            1.05 // 5% boost for moderate layouts
        } else if branch_density < 0.05 {
            1.02 // 2% boost for dense layouts
        } else {
            1.0 // No boost for very dense layouts
        }
    }

    /// Validate and adjust wavelength for manufacturing constraints (enhanced)
    fn validate_wavelength_for_diameter(&self, wavelength: f64, channel_width: f64) -> f64 {
        // Conservative minimum wavelength calculation for proper serpentine spacing
        let min_separation = channel_width;
        let min_wavelength = min_separation * 3.0; // Increased for better serpentine channel spacing

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
    fn calculate_effective_wavelength(
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
    fn calculate_effective_square_sharpness(
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
    fn calculate_required_wave_points(&self, half_periods: f64, base_points: usize) -> usize {
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
                        // Move point toward neighbor midpoint to open local curvature radius.
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
    fn generate_serpentine_path(
        &self,
        p1: Point2D,
        p2: Point2D,
        context: &ChannelGenerationContext,
    ) -> Vec<Point2D> {
        // Check if amplitude is below threshold - if so, return straight line
        let initial_wavelength =
            self.config.wavelength_factor * context.geometry_config.channel_width;
        let wavelength = self.validate_wavelength_for_diameter(
            initial_wavelength,
            context.geometry_config.channel_width,
        );
        let amplitude = self.calculate_adaptive_amplitude(p1, p2, context, wavelength);

        if amplitude <= 0.0 {
            // Return straight line when amplitude is too small for meaningful serpentines
            return self.generate_straight_line_path(
                p1,
                p2,
                context.geometry_config.generation.serpentine_points,
            );
        }

        let dx = p2.0 - p1.0;
        let dy = p2.1 - p1.1;
        let channel_length = dx.hypot(dy); // More efficient than sqrt(dx*dx + dy*dy)
        let _angle = dy.atan2(dx);

        let constants = ConstantsRegistry::new();
        let _branch_factor = (context.total_branches as f64)
            .powf(constants.get_branch_factor_exponent())
            .max(1.0);

        // Calculate number of periods to ensure complete wave cycles
        let initial_wavelength =
            self.config.wavelength_factor * context.geometry_config.channel_width;
        let base_wavelength = self.validate_wavelength_for_diameter(
            initial_wavelength,
            context.geometry_config.channel_width,
        );

        // Calculate amplitude with advanced adaptive algorithms
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

        // For smooth endpoint transitions, use half-periods to ensure zero amplitude at endpoints
        // Scale the number of periods with channel length and ensure minimum complete cycles
        let length_based_periods =
            (channel_length / effective_wavelength) * self.config.wave_density_factor;
        let base_periods = length_based_periods.max(1.0); // Minimum 1 complete cycle
        let requested_half_periods = (base_periods * 2.0).max(1.0);

        // Increase point density when many waves are requested so non-square waves keep smooth crests.
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
            // Ensure at least one positive and one negative lobe for sine-like waves.
            half_periods = 2.0;
        }
        if !matches!(self.config.wave_shape, crate::config::WaveShape::Square) {
            let central_window = (1.0 - start_guard - end_guard).max(0.0);
            if central_window < 0.35 && max_resolvable_half_periods >= 4.0 {
                // In strongly guarded segments, keep at least two full sign changes in
                // the middle region so both positive and negative lobes remain visible.
                half_periods = half_periods.max(4.0).min(max_resolvable_half_periods);
            }
        }

        // Calculate wave phase direction for perfect mirror symmetry
        let phase_direction = self.calculate_wave_phase_direction(p1, p2, context);

        // Pre-allocate path with exact capacity
        let mut path = Vec::with_capacity(n_points);
        let perp_x = -dy / channel_length;
        let perp_y = dx / channel_length;
        let smoothstep = |v: f64| {
            let x = v.clamp(0.0, 1.0);
            x * x * 2.0f64.mul_add(-x, 3.0)
        };

        for i in 0..n_points {
            let t = i as f64 / (n_points - 1) as f64;

            // Base position along the line
            let base_x = t.mul_add(dx, p1.0);
            let base_y = t.mul_add(dy, p1.1);

            // Use the improved envelope that respects adaptive configuration
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

            // Keep amplitude suppressed near split/merge nodes until branches have enough separation.
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

            // Serpentine wave with half-periods to ensure zero amplitude at endpoints
            let wave_phase = std::f64::consts::PI * half_periods * t;

            // Apply phase direction correctly for bilateral mirror symmetry
            // phase_direction determines the initial phase offset, not frequency scaling
            let phase_offset = if phase_direction < 0.0 {
                std::f64::consts::PI // Negative phase: start with inverted sine wave (π phase)
            } else {
                0.0 // Neutral/positive phase: start with sine wave (0 phase)
            };

            let local_cap = self.calculate_local_amplitude_cap(base_y, perp_y.abs(), context);
            let mut capped_amplitude = (initial_amplitude * guarded_envelope).min(local_cap);

            // Additional endpoint-distance cap inside guarded node zones to keep
            // channels from bulging into split/merge connectors immediately.
            if t < start_guard || t > (1.0 - end_guard) {
                let distance_to_nearest_endpoint = t.min(1.0 - t) * channel_length;
                let endpoint_distance_cap = distance_to_nearest_endpoint * 0.25;
                capped_amplitude = capped_amplitude.min(endpoint_distance_cap);
            }
            let wave_amplitude = capped_amplitude
                * self.calculate_wave_amplitude(wave_phase, phase_offset, square_sharpness);

            let x = wave_amplitude.mul_add(perp_x, base_x);
            let y = wave_amplitude.mul_add(perp_y, base_y);

            // Ensure exact endpoint matching for first and last points to maintain precision
            // The smooth envelope should make wave_amplitude ≈ 0 at endpoints, but we ensure exactness
            if i == 0 {
                path.push(p1);
            } else if i == n_points - 1 {
                path.push(p2);
            } else {
                path.push((x, y));
            }
        }

        // Final geometric safety pass: enforce diameter-scaled minimum bend radius
        // to prevent sharp local points and self-overlap tendencies.
        let min_radius = if matches!(self.config.wave_shape, crate::config::WaveShape::Square) {
            // Square waves need a stronger radius floor to avoid knife-edge turns
            // at larger diameters.
            channel_diameter
        } else {
            channel_diameter * 0.5
        };
        self.enforce_min_bend_radius(&mut path, min_radius, 4);

        path
    }

    /// Generate an optimized serpentine path between two points
    fn generate_optimized_serpentine_path(
        &self,
        p1: Point2D,
        p2: Point2D,
        context: &ChannelGenerationContext,
    ) -> Vec<Point2D> {
        // Check if amplitude is below threshold - if so, return straight line
        let initial_wavelength =
            self.config.wavelength_factor * context.geometry_config.channel_width;
        let wavelength = self.validate_wavelength_for_diameter(
            initial_wavelength,
            context.geometry_config.channel_width,
        );
        let amplitude = self.calculate_adaptive_amplitude(p1, p2, context, wavelength);

        if amplitude <= 0.0 {
            // Return straight line when amplitude is too small for meaningful serpentines
            return self.generate_straight_line_path(
                p1,
                p2,
                context.geometry_config.generation.serpentine_points,
            );
        }

        // Run optimization to find best parameters
        let optimization_result = optimize_serpentine_parameters(
            p1,
            p2,
            context.geometry_config,
            &self.config,
            context.box_dims,
            context.neighbor_info,
        );

        // Create optimized configuration without full clone
        let optimized_config = SerpentineConfig {
            wavelength_factor: optimization_result.params.wavelength_factor,
            wave_density_factor: optimization_result.params.wave_density_factor,
            fill_factor: optimization_result.params.fill_factor,
            gaussian_width_factor: self.config.gaussian_width_factor,
            wave_phase_direction: self.config.wave_phase_direction,
            wave_shape: self.config.wave_shape,
            optimization_enabled: false, // Disable nested optimization
            target_fill_ratio: self.config.target_fill_ratio,
            optimization_profile: self.config.optimization_profile,
            adaptive_config: self.config.adaptive_config,
        };

        // Generate path with optimized parameters using temporary strategy
        let temp_strategy = Self::new(optimized_config);
        temp_strategy.generate_serpentine_path(p1, p2, context)
    }

    /// Generate serpentine path for optimization purposes (public interface)
    #[must_use]
    pub fn generate_serpentine_path_for_optimization(
        &self,
        p1: Point2D,
        p2: Point2D,
        geometry_config: &GeometryConfig,
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

    /// Calculate wave phase direction for perfect bilateral mirror symmetry using enhanced symmetry system
    fn calculate_wave_phase_direction(
        &self,
        p1: Point2D,
        p2: Point2D,
        context: &ChannelGenerationContext,
    ) -> f64 {
        // If wave phase direction is explicitly set, use it
        if self.config.wave_phase_direction.abs() > 1e-6 {
            return self.config.wave_phase_direction;
        }

        // Use enhanced bilateral symmetry system for perfect symmetry
        let symmetry_config = BilateralSymmetryConfig {
            enable_vertical_symmetry: true,
            enable_horizontal_symmetry: true,
            symmetry_tolerance: 1e-6,
            enable_adaptive_symmetry: true,
            enforcement_strength: 1.0,
        };

        // Create a temporary channel generation context for symmetry calculation
        let temp_context = crate::state_management::adaptive::ChannelGenerationContext::new(
            *context.geometry_config,
            context.box_dims,
            context.total_branches,
            context.neighbor_info,
        )
        .with_endpoints(p1, p2);

        // Create symmetry context
        let symmetry_context = SymmetryContext::new(temp_context, symmetry_config);

        // Use bilateral phase direction calculator for perfect symmetry
        let phase_calculator = BilateralPhaseDirectionCalculator::default();

        match phase_calculator.calculate_phase_direction(&symmetry_context) {
            Ok(phase_direction) => phase_direction,
            Err(_) => self.calculate_wave_phase_direction_legacy(p1, p2, context.box_dims),
        }
    }

    /// Legacy fallback for phase direction that preserves centerline mirror behavior.
    fn calculate_wave_phase_direction_legacy(
        &self,
        p1: Point2D,
        p2: Point2D,
        box_dims: (f64, f64),
    ) -> f64 {
        let center_y = box_dims.1 / 2.0;
        let channel_center_y = f64::midpoint(p1.1, p2.1);
        let tolerance = box_dims.1 * 0.01;

        if (channel_center_y - center_y).abs() <= tolerance {
            0.0
        } else if channel_center_y > center_y {
            1.0
        } else {
            -1.0
        }
    }
}
