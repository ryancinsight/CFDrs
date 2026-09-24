use super::{ArcChannelStrategy, ConstantsRegistry, GeometryConfig, Point2D};

impl ArcChannelStrategy {
    /// Generate symmetric arc with specific direction
    #[allow(clippy::too_many_arguments)]
    pub(super) fn generate_symmetric_arc_with_direction(
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
    pub(super) fn calculate_max_offset_toward_direction(
        &self,
        channel_center_y: f64,
        direction: f64,
        geometry_config: &GeometryConfig,
        box_dims: (f64, f64),
        neighbor_info: Option<&[f64]>,
    ) -> f64 {
        let wall_margin =
            geometry_config.wall_clearance_mm() + geometry_config.channel_width_mm() * 0.5;
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
                        Some(delta - geometry_config.channel_width_mm())
                    } else if direction < 0.0 && delta < 0.0 {
                        Some(-delta - geometry_config.channel_width_mm())
                    } else {
                        None
                    }
                })
                .fold(f64::INFINITY, f64::min),
            None => f64::INFINITY,
        };

        wall_limit.min(neighbor_limit).max(0.0)
    }

}
