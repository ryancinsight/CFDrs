use super::{ArcChannelStrategy, ArcConfig, ConstantsRegistry, GeometryConfig, Point2D};

impl ArcChannelStrategy {
    /// Generate arc path with enhanced bilateral mirror symmetry
    pub(super) fn generate_arc_path_with_enhanced_symmetry(
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

}
