//! Arc channel generation strategy.

use crate::config::{ArcConfig, ConstantsRegistry, GeometryConfig};
use crate::geometry::{ChannelType, Point2D};

use super::ChannelTypeStrategy;

/// Strategy for creating arc channels
#[derive(Debug, Clone)]
pub struct ArcChannelStrategy {
    config: ArcConfig,
}

impl ArcChannelStrategy {
    /// Create a new arc channel strategy with the given configuration
    ///
    /// # Arguments
    ///
    /// * `config` - Configuration parameters for arc channel generation
    ///
    /// # Examples
    ///
    /// ```rust
    /// use cfd_schematics::geometry::strategies::ArcChannelStrategy;
    /// use cfd_schematics::config::ArcConfig;
    ///
    /// let strategy = ArcChannelStrategy::new(ArcConfig::default());
    /// ```
    #[must_use]
    pub const fn new(config: ArcConfig) -> Self {
        Self { config }
    }
}

impl ChannelTypeStrategy for ArcChannelStrategy {
    fn create_channel(
        &self,
        from: Point2D,
        to: Point2D,
        geometry_config: &GeometryConfig,
        box_dims: (f64, f64),
        total_branches: usize,
        neighbor_info: Option<&[f64]>,
    ) -> ChannelType {
        let path = self.generate_arc_path_with_enhanced_symmetry(
            from,
            to,
            geometry_config,
            box_dims,
            total_branches,
            neighbor_info,
        );
        ChannelType::Arc { path }
    }
}

impl ArcChannelStrategy {
    /// Generate arc path with enhanced bilateral mirror symmetry
    fn generate_arc_path_with_enhanced_symmetry(
        &self,
        p1: Point2D,
        p2: Point2D,
        geometry_config: &GeometryConfig,
        box_dims: (f64, f64),
        total_branches: usize,
        neighbor_info: Option<&[f64]>,
    ) -> Vec<Point2D> {
        // Check if this is a center channel in a trifurcation that needs figure-8 pattern
        if self.is_center_trifurcation_channel(p1, p2, box_dims, total_branches) {
            return self.generate_figure_eight_path(
                p1,
                p2,
                geometry_config,
                box_dims,
                neighbor_info,
            );
        }

        if !self.config.enable_collision_prevention {
            return self.generate_arc_path_with_bilateral_symmetry(
                p1,
                p2,
                geometry_config,
                box_dims,
                total_branches,
                neighbor_info,
            );
        }

        // Calculate adaptive curvature factor based on proximity to neighbors
        let adaptive_curvature = self.calculate_adaptive_curvature(
            p1,
            p2,
            geometry_config,
            box_dims,
            total_branches,
            neighbor_info,
        );

        // Create temporary config with adaptive curvature and enhanced symmetry
        let adaptive_config = ArcConfig {
            curvature_factor: adaptive_curvature,
            ..self.config
        };

        // Generate path with adaptive curvature and bilateral symmetry
        let temp_strategy = Self::new(adaptive_config);
        temp_strategy.generate_arc_path_with_bilateral_symmetry(
            p1,
            p2,
            geometry_config,
            box_dims,
            total_branches,
            neighbor_info,
        )
    }

    /// Check if this is a center channel in a trifurcation that needs figure-8 pattern
    fn is_center_trifurcation_channel(
        &self,
        p1: Point2D,
        p2: Point2D,
        box_dims: (f64, f64),
        total_branches: usize,
    ) -> bool {
        // Only apply to trifurcations (3 or more branches)
        if total_branches < 3 {
            return false;
        }

        let (_, height) = box_dims;
        let center_y = height / 2.0;
        let tolerance = height * 0.1; // 10% tolerance for center detection

        // Check if both points are near the vertical center
        let p1_near_center = (p1.1 - center_y).abs() < tolerance;
        let p2_near_center = (p2.1 - center_y).abs() < tolerance;

        p1_near_center && p2_near_center
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

    /// Compute endpoint guard length (fraction of channel length) for figure-8 center channels.
    fn calculate_figure_eight_junction_guard_fraction(
        &self,
        endpoint_y: f64,
        channel_length: f64,
        geometry_config: &GeometryConfig,
        neighbor_info: Option<&[f64]>,
    ) -> f64 {
        const MIN_GUARD_FRACTION: f64 = 0.15;
        const MAX_GUARD_FRACTION: f64 = 0.40;

        let channel_diameter = geometry_config.channel_width;
        let physical_guard_length = geometry_config.wall_clearance + channel_diameter * 20.0;
        let mut guard_fraction = if channel_length > 1e-6 {
            (physical_guard_length / channel_length).clamp(MIN_GUARD_FRACTION, 0.45)
        } else {
            0.45
        };

        if let Some(min_neighbor_distance) =
            self.min_neighbor_distance_at_y(endpoint_y, neighbor_info)
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

    /// Compute a local figure-8 amplitude cap from wall and neighbor spacing at base y.
    fn calculate_local_figure_eight_amplitude_cap(
        &self,
        base_y: f64,
        normal_y_factor: f64,
        geometry_config: &GeometryConfig,
        box_dims: (f64, f64),
        neighbor_info: Option<&[f64]>,
    ) -> f64 {
        let clear_up = self.calculate_max_offset_toward_direction(
            base_y,
            1.0,
            geometry_config,
            box_dims,
            neighbor_info,
        );
        let clear_down = self.calculate_max_offset_toward_direction(
            base_y,
            -1.0,
            geometry_config,
            box_dims,
            neighbor_info,
        );
        let vertical_cap = clear_up.min(clear_down).max(0.0);
        vertical_cap / normal_y_factor.max(0.1)
    }

    /// Generate figure-8 (weave) path for center channels in trifurcations
    fn generate_figure_eight_path(
        &self,
        p1: Point2D,
        p2: Point2D,
        geometry_config: &GeometryConfig,
        box_dims: (f64, f64),
        neighbor_info: Option<&[f64]>,
    ) -> Vec<Point2D> {
        let constants = ConstantsRegistry::new();
        let num_points = self.config.smoothness + 2;
        let mut path = Vec::with_capacity(num_points);

        let dx = p2.0 - p1.0;
        let dy = p2.1 - p1.1;
        let distance = dx.hypot(dy);

        if distance < constants.get_geometric_tolerance() {
            return vec![p1, p2];
        }

        let (_, height) = box_dims;
        let base_amplitude = height * 0.12 * self.config.curvature_factor;
        let clear_up = self.calculate_max_offset_toward_direction(
            f64::midpoint(p1.1, p2.1),
            1.0,
            geometry_config,
            box_dims,
            neighbor_info,
        );
        let clear_down = self.calculate_max_offset_toward_direction(
            f64::midpoint(p1.1, p2.1),
            -1.0,
            geometry_config,
            box_dims,
            neighbor_info,
        );
        let amplitude = base_amplitude.min(clear_up.min(clear_down) * 0.7);
        if amplitude <= constants.get_geometric_tolerance() {
            return vec![p1, p2];
        }

        let perp_x = -dy / distance;
        let perp_y = dx / distance;
        let normal_y_factor = perp_y.abs();
        let start_guard = self.calculate_figure_eight_junction_guard_fraction(
            p1.1,
            distance,
            geometry_config,
            neighbor_info,
        );
        let end_guard = self.calculate_figure_eight_junction_guard_fraction(
            p2.1,
            distance,
            geometry_config,
            neighbor_info,
        );
        let smoothstep = |v: f64| {
            let x = v.clamp(0.0, 1.0);
            x * x * 2.0f64.mul_add(-x, 3.0)
        };

        // Generate smooth figure-8 pattern with two crossing arcs
        for i in 0..num_points {
            let t = i as f64 / (num_points - 1) as f64;
            let base_x = t.mul_add(dx, p1.0);
            let base_y = t.mul_add(dy, p1.1);

            // Create smooth figure-8 with two arcs that cross in the middle
            // First half curves one way, second half curves the opposite way
            let wave_shape = if t < 0.5 {
                // First arc: smooth curve upward then back to center
                let local_t = t * 2.0; // Scale to 0-1 for first half
                (std::f64::consts::PI * local_t).sin()
            } else {
                // Second arc: smooth curve downward then back to center
                let local_t = (t - 0.5) * 2.0; // Scale to 0-1 for second half
                -(std::f64::consts::PI * local_t).sin()
            };

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

            let local_cap = self.calculate_local_figure_eight_amplitude_cap(
                base_y,
                normal_y_factor,
                geometry_config,
                box_dims,
                neighbor_info,
            );
            let mut capped_amplitude =
                (amplitude * start_guard_scale * end_guard_scale).min(local_cap);

            // Endpoint-distance cap keeps center branch from immediately
            // bulging into split/merge connectors near junction nodes.
            let distance_to_nearest_endpoint = t.min(1.0 - t) * distance;
            let endpoint_distance_cap = distance_to_nearest_endpoint * 0.16;
            capped_amplitude = capped_amplitude.min(endpoint_distance_cap);

            let wave_offset = wave_shape * capped_amplitude;
            let x = wave_offset.mul_add(perp_x, base_x);
            let y = wave_offset.mul_add(perp_y, base_y);

            if i == 0 {
                path.push(p1);
            } else if i == num_points - 1 {
                path.push(p2);
            } else {
                path.push((x, y));
            }
        }

        path
    }

    /// Generate arc path with enhanced bilateral mirror symmetry
    fn generate_arc_path_with_bilateral_symmetry(
        &self,
        p1: Point2D,
        p2: Point2D,
        geometry_config: &GeometryConfig,
        box_dims: (f64, f64),
        total_branches: usize,
        neighbor_info: Option<&[f64]>,
    ) -> Vec<Point2D> {
        let constants = ConstantsRegistry::new();
        let num_points = self.config.smoothness + 2;

        let dx = p2.0 - p1.0;
        let dy = p2.1 - p1.1;
        let distance = dx.hypot(dy);

        // For very short channels or zero curvature, return straight line
        if distance < constants.get_geometric_tolerance()
            || self.config.curvature_factor < constants.get_geometric_tolerance()
        {
            return vec![p1, p2];
        }

        // Calculate enhanced arc direction with bilateral symmetry
        let arc_direction =
            self.calculate_bilateral_symmetric_arc_direction(p1, p2, box_dims, total_branches);

        // Generate symmetric arc path
        self.generate_symmetric_arc_with_direction(
            p1,
            p2,
            geometry_config,
            box_dims,
            neighbor_info,
            arc_direction,
            num_points,
        )
    }

    /// Calculate bilateral symmetric arc direction for enhanced symmetry
    fn calculate_bilateral_symmetric_arc_direction(
        &self,
        p1: Point2D,
        p2: Point2D,
        box_dims: (f64, f64),
        total_branches: usize,
    ) -> f64 {
        // If curvature direction is explicitly set, use it
        if self.config.curvature_direction.abs() > 1e-6 {
            return self.config.curvature_direction;
        }

        let (_length, height) = box_dims;
        let center_y = height / 2.0;

        // Calculate channel position relative to centers
        let _channel_center_x = f64::midpoint(p1.0, p2.0);
        let channel_center_y = f64::midpoint(p1.1, p2.1);
        // Determine if channel is peripheral or internal
        let is_peripheral = self.is_peripheral_channel(p1, p2, box_dims, total_branches);

        if is_peripheral {
            // Peripheral channels curve toward walls
            if channel_center_y > center_y {
                1.0 // Upper peripheral channels curve upward (toward top wall)
            } else {
                -1.0 // Lower peripheral channels curve downward (toward bottom wall)
            }
        } else {
            // Internal channels curve toward center
            if channel_center_y > center_y {
                -1.0 // Upper internal channels curve downward (toward center)
            } else {
                1.0 // Lower internal channels curve upward (toward center)
            }
        }
    }

    /// Check if channel is peripheral (outer) vs internal
    fn is_peripheral_channel(
        &self,
        p1: Point2D,
        p2: Point2D,
        box_dims: (f64, f64),
        total_branches: usize,
    ) -> bool {
        let (_, height) = box_dims;
        let center_y = height / 2.0;
        let channel_center_y = f64::midpoint(p1.1, p2.1);

        // For bifurcations (2 branches), both are peripheral
        if total_branches <= 2 {
            return true;
        }

        // For trifurcations and higher, determine based on distance from center
        let distance_from_center = (channel_center_y - center_y).abs();
        let threshold = (height / (total_branches as f64 + 1.0)).max(height * 0.1);

        distance_from_center > threshold
    }

    /// Generate symmetric arc with specific direction
    #[allow(clippy::too_many_arguments)]
    fn generate_symmetric_arc_with_direction(
        &self,
        p1: Point2D,
        p2: Point2D,
        geometry_config: &GeometryConfig,
        box_dims: (f64, f64),
        neighbor_info: Option<&[f64]>,
        arc_direction: f64,
        num_points: usize,
    ) -> Vec<Point2D> {
        let mut path = Vec::with_capacity(num_points);

        let dx = p2.0 - p1.0;
        let dy = p2.1 - p1.1;
        let distance = dx.hypot(dy);
        if distance <= ConstantsRegistry::new().get_geometric_tolerance() {
            return vec![p1, p2];
        }

        // Calculate perpendicular direction for arc curvature
        let perp_x = -dy / distance;
        let perp_y = dx / distance;

        // Apply directional multiplier
        let directed_perp_x = perp_x * arc_direction;
        let directed_perp_y = perp_y * arc_direction;

        // Arc height based on curvature factor and constrained by wall/neighbor spacing.
        let base_arc_height = distance * self.config.curvature_factor * 0.5;
        let channel_center_y = f64::midpoint(p1.1, p2.1);
        let max_offset = if arc_direction.abs() > 1e-6 {
            self.calculate_max_offset_toward_direction(
                channel_center_y,
                arc_direction,
                geometry_config,
                box_dims,
                neighbor_info,
            )
        } else {
            let up = self.calculate_max_offset_toward_direction(
                channel_center_y,
                1.0,
                geometry_config,
                box_dims,
                neighbor_info,
            );
            let down = self.calculate_max_offset_toward_direction(
                channel_center_y,
                -1.0,
                geometry_config,
                box_dims,
                neighbor_info,
            );
            up.min(down)
        };
        let arc_height = base_arc_height.min(max_offset * 0.9);
        if arc_height <= ConstantsRegistry::new().get_geometric_tolerance() {
            return vec![p1, p2];
        }

        // Generate smooth arc using quadratic Bezier curve
        for i in 0..num_points {
            let t = i as f64 / (num_points - 1) as f64;

            // Quadratic Bezier: B(t) = (1-t)²P₀ + 2(1-t)tP₁ + t²P₂
            let control_x = directed_perp_x.mul_add(arc_height, f64::midpoint(p1.0, p2.0));
            let control_y = directed_perp_y.mul_add(arc_height, f64::midpoint(p1.1, p2.1));

            let one_minus_t = 1.0 - t;
            let one_minus_t_sq = one_minus_t * one_minus_t;
            let t_sq = t * t;
            let two_t_one_minus_t = 2.0 * t * one_minus_t;

            let x = t_sq.mul_add(
                p2.0,
                one_minus_t_sq.mul_add(p1.0, two_t_one_minus_t * control_x),
            );
            let y = t_sq.mul_add(
                p2.1,
                one_minus_t_sq.mul_add(p1.1, two_t_one_minus_t * control_y),
            );

            path.push((x, y));
        }

        path
    }

    /// Compute the maximum safe lateral offset in a direction while respecting
    /// wall clearance and neighbor spacing.
    fn calculate_max_offset_toward_direction(
        &self,
        channel_center_y: f64,
        direction: f64,
        geometry_config: &GeometryConfig,
        box_dims: (f64, f64),
        neighbor_info: Option<&[f64]>,
    ) -> f64 {
        let wall_margin = geometry_config.wall_clearance + geometry_config.channel_width * 0.5;
        let (_, box_height) = box_dims;
        let wall_limit = if direction > 0.0 {
            box_height - channel_center_y - wall_margin
        } else {
            channel_center_y - wall_margin
        };

        let neighbor_limit = match neighbor_info {
            Some(neighbors) => neighbors
                .iter()
                .filter_map(|&neighbor_y| {
                    let delta = neighbor_y - channel_center_y;
                    if delta.abs() <= 0.1 {
                        return None;
                    }

                    if direction > 0.0 && delta > 0.0 {
                        Some(delta - geometry_config.channel_width)
                    } else if direction < 0.0 && delta < 0.0 {
                        Some(-delta - geometry_config.channel_width)
                    } else {
                        None
                    }
                })
                .fold(f64::INFINITY, f64::min),
            None => f64::INFINITY,
        };

        wall_limit.min(neighbor_limit).max(0.0)
    }

    /// Calculate adaptive curvature factor based on neighbor proximity
    fn calculate_adaptive_curvature(
        &self,
        p1: Point2D,
        p2: Point2D,
        geometry_config: &GeometryConfig,
        box_dims: (f64, f64),
        total_branches: usize,
        neighbor_info: Option<&[f64]>,
    ) -> f64 {
        let constants = ConstantsRegistry::new();
        if !self.config.enable_adaptive_curvature {
            return self.config.curvature_factor;
        }

        let dx = p2.0 - p1.0;
        let dy = p2.1 - p1.1;
        let channel_length = dx.hypot(dy);

        // Base curvature factor
        let mut adaptive_factor = self.config.curvature_factor;

        // Calculate proximity-based reduction
        let proximity_reduction = self.calculate_proximity_reduction(
            p1,
            p2,
            geometry_config.channel_width,
            box_dims,
            total_branches,
            neighbor_info,
        );

        // Apply proximity reduction with limits
        adaptive_factor *= (1.0 - proximity_reduction).max(self.config.max_curvature_reduction);

        // Additional safety check for very short channels
        if channel_length
            < geometry_config.channel_width * constants.get_short_channel_width_multiplier()
        {
            adaptive_factor *= constants.get_max_curvature_reduction_factor(); // Reduce curvature for very short channels
        }

        // Ensure we don't go below minimum curvature
        adaptive_factor.max(constants.get_min_curvature_factor())
    }

    /// Calculate proximity-based curvature reduction factor
    fn calculate_proximity_reduction(
        &self,
        p1: Point2D,
        p2: Point2D,
        channel_diameter: f64,
        box_dims: (f64, f64),
        total_branches: usize,
        neighbor_info: Option<&[f64]>,
    ) -> f64 {
        // If we don't have neighbor information, use branch density estimation
        let channel_center_y = f64::midpoint(p1.1, p2.1);
        let neighbor_distances: Vec<f64> = match neighbor_info {
            Some(neighbors) => neighbors
                .iter()
                .map(|neighbor_y| (neighbor_y - channel_center_y).abs())
                .filter(|distance| *distance > 0.1)
                .collect(),
            None => return self.estimate_density_based_reduction(p1, p2, box_dims, total_branches),
        };
        let mut max_reduction: f64 = 0.0;

        // Calculate channel midpoint for proximity calculations
        let _mid_x = f64::midpoint(p1.0, p2.0);
        let _mid_y = f64::midpoint(p1.1, p2.1);

        // Check proximity to each neighbor
        for neighbor_distance in neighbor_distances {
            if neighbor_distance < channel_diameter {
                // Calculate reduction factor based on how close the neighbor is
                let proximity_ratio = neighbor_distance / channel_diameter;
                let reduction = (1.0 - proximity_ratio).max(0.0);
                max_reduction = max_reduction.max(reduction);
            }
        }

        // Apply maximum reduction limit
        max_reduction.min(1.0 - self.config.max_curvature_reduction)
    }

    /// Estimate curvature reduction based on branch density
    fn estimate_density_based_reduction(
        &self,
        p1: Point2D,
        p2: Point2D,
        box_dims: (f64, f64),
        total_branches: usize,
    ) -> f64 {
        // Calculate effective area per branch
        let box_area = box_dims.0 * box_dims.1;
        let area_per_branch = box_area / total_branches as f64;

        // Calculate channel length
        let dx = p2.0 - p1.0;
        let dy = p2.1 - p1.1;
        let channel_length = dx.hypot(dy);

        // Estimate potential arc area
        let potential_arc_area = channel_length * channel_length * self.config.curvature_factor;

        // If potential arc area is large relative to available space, reduce curvature
        if potential_arc_area > area_per_branch * 0.5 {
            let density_ratio = potential_arc_area / (area_per_branch * 0.5);
            let reduction = (density_ratio - 1.0).clamp(0.0, 0.8);
            return reduction;
        }

        0.0 // No reduction needed
    }
}
