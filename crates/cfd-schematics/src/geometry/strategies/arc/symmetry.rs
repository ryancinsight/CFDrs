use super::{ArcChannelStrategy, Point2D};

impl ArcChannelStrategy {
    /// Calculate bilateral symmetric arc direction for enhanced symmetry
    pub(super) fn calculate_bilateral_symmetric_arc_direction(
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

}
