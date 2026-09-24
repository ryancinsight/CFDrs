use super::{ArcChannelStrategy, ConstantsRegistry, GeometryConfig, Point2D};

impl ArcChannelStrategy {
    /// Calculate adaptive curvature factor based on neighbor proximity
    pub(super) fn calculate_adaptive_curvature(
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
            geometry_config.channel_width_mm(),
            box_dims,
            total_branches,
            neighbor_info,
        );

        // Apply proximity reduction with limits
        adaptive_factor *= (1.0 - proximity_reduction).max(self.config.max_curvature_reduction);

        // Additional safety check for very short channels
        if channel_length
            < geometry_config.channel_width_mm() * constants.get_short_channel_width_multiplier()
        {
            adaptive_factor *= constants.get_max_curvature_reduction_factor();
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
