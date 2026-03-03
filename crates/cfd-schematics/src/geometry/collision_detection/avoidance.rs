//! Collision avoidance strategies for adaptive path reduction.

use super::{
    CollisionAvoidanceResult, CollisionContext, CollisionDetectionResult, CollisionDetectionSystem,
    CollisionParameters, CollisionSeverity,
};
use crate::{
    config::ConstantsRegistry, error::SchemeResult, geometry::Point2D,
    state_management::adaptive::ChannelGenerationContext,
};

impl CollisionDetectionSystem {
    /// Apply collision avoidance to a path
    ///
    /// # Errors
    ///
    /// Returns an error if collision detection fails or if path modification
    /// encounters computational issues during the avoidance process.
    pub fn apply_collision_avoidance(
        &mut self,
        path: &mut Vec<Point2D>,
        context: &CollisionContext,
    ) -> SchemeResult<CollisionAvoidanceResult> {
        let detection_result = self.detect_collisions(path, context)?;

        if !detection_result.has_collisions {
            return Ok(CollisionAvoidanceResult {
                applied: false,
                reduction_factor: 1.0,
                original_severity: 0.0,
                final_severity: 0.0,
            });
        }

        // Get adaptive parameters for avoidance strategy
        let params = Self::get_collision_parameters(Some(&context.channel_context));
        let original_severity = detection_result.severity_score;

        // Apply adaptive avoidance strategies
        let reduction_factor = Self::calculate_adaptive_reduction_factor(
            &detection_result,
            &params,
            &context.channel_context,
        );
        Self::apply_adaptive_path_reduction(path, context, reduction_factor);

        // Verify improvement with adaptive parameters
        let final_detection = self.detect_collisions(path, context)?;
        let final_severity = final_detection.severity_score;

        Ok(CollisionAvoidanceResult {
            applied: true,
            reduction_factor,
            original_severity,
            final_severity,
        })
    }

    /// Calculate adaptive reduction factor based on collision severity and context
    fn calculate_adaptive_reduction_factor(
        detection_result: &CollisionDetectionResult,
        params: &CollisionParameters,
        context: &ChannelGenerationContext,
    ) -> f64 {
        let neighbor_max = detection_result
            .neighbor_collisions
            .iter()
            .map(CollisionSeverity::get_severity)
            .fold(0.0, f64::max);

        let wall_max = detection_result
            .wall_collisions
            .iter()
            .map(CollisionSeverity::get_severity)
            .fold(0.0, f64::max);

        let max_severity = neighbor_max.max(wall_max);

        // Base reduction factor
        let mut reduction_factor = max_severity * params.max_reduction_factor;

        // Apply adaptive adjustments based on context
        if params.adaptive_enabled {
            // Adjust based on channel length - longer channels can handle more reduction
            let channel_length = context.channel_length();
            if channel_length > 50.0 {
                reduction_factor *= 1.1;
            } else if channel_length < 20.0 {
                reduction_factor *= 0.9;
            }

            // Adjust based on neighbor density
            if let Some(neighbor_info) = &context.neighbor_info {
                let neighbor_density = neighbor_info.len() as f64 / context.total_branches as f64;
                if neighbor_density > 0.8 {
                    // High density - be more aggressive with reduction
                    reduction_factor *= 1.2;
                }
            }

            // Adjust based on branch count
            let constants = ConstantsRegistry::new();
            let branch_factor =
                (context.total_branches as f64).powf(constants.get_branch_factor_exponent());
            if branch_factor > 2.0 {
                reduction_factor *= (branch_factor - 2.0).mul_add(0.1, 1.0);
            }
        }

        reduction_factor.min(params.max_reduction_factor)
    }

    /// Apply adaptive path reduction to avoid collisions
    fn apply_adaptive_path_reduction(
        path: &mut [Point2D],
        context: &CollisionContext,
        reduction_factor: f64,
    ) {
        if reduction_factor <= 0.0 {
            return;
        }

        let start = context.current_channel.start;
        let end = context.current_channel.end;
        let channel_context = &context.channel_context;

        // Apply different reduction strategies based on context
        if channel_context.total_branches > 8 {
            // High branch count - use more sophisticated reduction
            Self::apply_sophisticated_reduction(
                path,
                start,
                end,
                reduction_factor,
                channel_context,
            );
        } else {
            // Standard reduction for simpler cases
            Self::apply_standard_reduction(path, start, end, reduction_factor);
        }
    }

    /// Apply sophisticated reduction for complex channel systems
    fn apply_sophisticated_reduction(
        path: &mut [Point2D],
        start: Point2D,
        end: Point2D,
        reduction_factor: f64,
        context: &ChannelGenerationContext,
    ) {
        let path_len = path.len();

        // Use adaptive reduction that varies along the path
        for (i, point) in path.iter_mut().enumerate() {
            let t = i as f64 / (path_len - 1) as f64;

            // Calculate adaptive reduction factor based on position
            let position_factor = if (0.2..=0.8).contains(&t) {
                // Standard reduction in the middle
                reduction_factor
            } else {
                // Reduce more aggressively at endpoints
                reduction_factor * 1.2
            };

            // Consider neighbor proximity for local adjustments
            let local_reduction = if let Some(neighbor_info) = &context.neighbor_info {
                let current_y = point.1;
                let min_neighbor_dist = neighbor_info
                    .iter()
                    .map(|&ny| (current_y - ny).abs())
                    .fold(f64::INFINITY, f64::min);

                if min_neighbor_dist < 5.0 {
                    position_factor * 1.3 // More aggressive reduction near neighbors
                } else {
                    position_factor
                }
            } else {
                position_factor
            };

            let straight_x = t.mul_add(end.0 - start.0, start.0);
            let straight_y = t.mul_add(end.1 - start.1, start.1);

            // Apply adaptive interpolation
            let final_reduction = local_reduction.min(0.95); // Cap at 95% reduction
            point.0 = point
                .0
                .mul_add(1.0 - final_reduction, straight_x * final_reduction);
            point.1 = point
                .1
                .mul_add(1.0 - final_reduction, straight_y * final_reduction);
        }
    }

    /// Apply standard path reduction (legacy method)
    fn apply_standard_reduction(
        path: &mut [Point2D],
        start: Point2D,
        end: Point2D,
        reduction_factor: f64,
    ) {
        let path_len = path.len();
        for (i, point) in path.iter_mut().enumerate() {
            let t = i as f64 / (path_len - 1) as f64;
            let straight_x = t.mul_add(end.0 - start.0, start.0);
            let straight_y = t.mul_add(end.1 - start.1, start.1);

            // Interpolate between current path and straight line
            point.0 = point
                .0
                .mul_add(1.0 - reduction_factor, straight_x * reduction_factor);
            point.1 = point
                .1
                .mul_add(1.0 - reduction_factor, straight_y * reduction_factor);
        }
    }
}
