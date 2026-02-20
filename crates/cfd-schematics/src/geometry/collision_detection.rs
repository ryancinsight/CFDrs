//! Enhanced collision detection system with centralized parameter management
//!
//! This module provides a modular collision detection and avoidance system
//! that integrates with the centralized state management system and supports
//! both channel-to-channel and channel-to-wall boundary constraints.

use crate::{
    config::ConstantsRegistry,
    error::{ConfigurationError, SchemeError, SchemeResult},
    geometry::Point2D,
    state_management::{adaptive::ChannelGenerationContext, ParameterRegistry},
};

/// Enhanced collision detection context that integrates with adaptive parameter system
#[derive(Debug, Clone)]
pub struct CollisionContext {
    /// Channel generation context for adaptive behavior
    pub channel_context: ChannelGenerationContext,

    /// Information about neighboring channels
    pub neighbor_info: Vec<NeighborInfo>,

    /// Wall boundaries
    pub wall_boundaries: WallBoundaries,

    /// Current channel being processed
    pub current_channel: ChannelInfo,
}

impl CollisionContext {
    /// Create a new collision context from channel generation context
    #[must_use]
    pub const fn from_channel_context(
        channel_context: ChannelGenerationContext,
        neighbor_info: Vec<NeighborInfo>,
        wall_boundaries: WallBoundaries,
        current_channel: ChannelInfo,
    ) -> Self {
        Self {
            channel_context,
            neighbor_info,
            wall_boundaries,
            current_channel,
        }
    }

    /// Get the underlying channel generation context
    #[must_use]
    pub const fn channel_context(&self) -> &ChannelGenerationContext {
        &self.channel_context
    }
}

/// Information about a neighboring channel
#[derive(Debug, Clone)]
pub struct NeighborInfo {
    /// Y-coordinate of the neighbor
    pub y_position: f64,

    /// Channel width of the neighbor
    pub width: f64,

    /// Distance to this neighbor
    pub distance: f64,

    /// Whether this neighbor is active/relevant
    pub is_active: bool,
}

/// Wall boundary information
#[derive(Debug, Clone)]
pub struct WallBoundaries {
    /// Left wall x-coordinate
    pub left: f64,

    /// Right wall x-coordinate
    pub right: f64,

    /// Bottom wall y-coordinate
    pub bottom: f64,

    /// Top wall y-coordinate
    pub top: f64,
}

/// Information about the current channel
#[derive(Debug, Clone)]
pub struct ChannelInfo {
    /// Start point of the channel
    pub start: Point2D,

    /// End point of the channel
    pub end: Point2D,

    /// Channel width
    pub width: f64,

    /// Channel index in the system
    pub index: usize,
}

/// Collision detection and avoidance system
pub struct CollisionDetectionSystem {
    /// Parameter registry for collision parameters
    #[allow(dead_code)] // Part of comprehensive collision detection framework
    registry: ParameterRegistry,

    /// Cached collision parameters
    #[allow(dead_code)] // Part of comprehensive collision detection framework
    cached_params: Option<CollisionParameters>,
}

/// Enhanced collision detection parameters with adaptive behavior
#[derive(Debug, Clone)]
pub struct CollisionParameters {
    /// Minimum distance between channels (adaptive)
    pub min_channel_distance: f64,

    /// Minimum distance from walls (adaptive)
    pub min_wall_distance: f64,

    /// Safety margin factor (adaptive)
    pub safety_margin_factor: f64,

    /// Enable neighbor-based collision detection
    pub enable_neighbor_detection: bool,

    /// Enable wall-based collision detection
    pub enable_wall_detection: bool,

    /// Maximum reduction factor for collision avoidance (adaptive)
    pub max_reduction_factor: f64,

    /// Collision detection sensitivity (adaptive)
    pub detection_sensitivity: f64,

    /// Adaptive behavior enabled
    pub adaptive_enabled: bool,
}

impl CollisionParameters {
    /// Create parameters from constants registry with optional adaptive context
    #[must_use]
    pub fn from_constants_registry(
        constants: &ConstantsRegistry,
        context: Option<&ChannelGenerationContext>,
    ) -> Self {
        // Base parameters from constants registry
        let mut params = Self {
            min_channel_distance: constants.get_min_channel_distance(),
            min_wall_distance: constants.get_min_wall_distance(),
            safety_margin_factor: constants.get_safety_margin_factor(),
            enable_neighbor_detection: true,
            enable_wall_detection: true,
            max_reduction_factor: constants.get_max_reduction_factor(),
            detection_sensitivity: constants.get_detection_sensitivity(),
            adaptive_enabled: context.is_some(),
        };

        // Apply adaptive adjustments if context is provided
        if let Some(ctx) = context {
            params.apply_adaptive_adjustments(ctx, constants);
        }

        params
    }

    /// Apply adaptive adjustments based on context
    fn apply_adaptive_adjustments(
        &mut self,
        context: &ChannelGenerationContext,
        constants: &ConstantsRegistry,
    ) {
        // Adjust minimum distances based on neighbor proximity
        if let Some(min_neighbor_dist) = context.min_neighbor_distance() {
            // Reduce minimum distances when neighbors are close
            let proximity_factor = (min_neighbor_dist / constants.get_proximity_divisor()).clamp(
                constants.get_min_proximity_factor(),
                constants.get_max_proximity_factor(),
            );
            self.min_channel_distance *= proximity_factor;
            self.min_wall_distance *= proximity_factor;
        }

        // Adjust sensitivity based on branch count
        let branch_factor = constants.get_branch_factor_exponent();
        let branch_adjustment = (context.total_branches as f64).powf(branch_factor)
            / constants.get_branch_adjustment_divisor();
        self.detection_sensitivity *=
            (1.0 + branch_adjustment).min(constants.get_max_sensitivity_multiplier());

        // Adjust reduction factor based on channel length
        let channel_length = context.channel_length();
        if channel_length > constants.get_long_channel_threshold() {
            // Longer channels can tolerate more reduction
            self.max_reduction_factor = (self.max_reduction_factor
                * constants.get_long_channel_reduction_multiplier())
            .min(constants.get_max_reduction_limit());
        }
    }
}

impl CollisionDetectionSystem {
    /// Create a new collision detection system
    ///
    /// # Errors
    ///
    /// Returns an error if the parameter registry cannot be initialized with default values.
    pub fn new() -> SchemeResult<Self> {
        let registry = ParameterRegistry::with_defaults().map_err(|e| {
            SchemeError::Configuration(ConfigurationError::MissingConfiguration {
                field: e.to_string(),
            })
        })?;

        Ok(Self {
            registry,
            cached_params: None,
        })
    }

    /// Create with existing parameter registry
    #[must_use]
    pub const fn with_registry(registry: ParameterRegistry) -> Self {
        Self {
            registry,
            cached_params: None,
        }
    }

    /// Get collision parameters with adaptive behavior
    fn get_collision_parameters(context: Option<&ChannelGenerationContext>) -> CollisionParameters {
        let constants = ConstantsRegistry::new();

        // Always create fresh parameters to ensure adaptive behavior is applied
        CollisionParameters::from_constants_registry(&constants, context)
    }

    /// Detect collisions for a given path with adaptive parameter behavior
    ///
    /// # Errors
    ///
    /// Returns an error if collision detection fails due to invalid parameters or
    /// computational issues during the detection process.
    pub fn detect_collisions(
        &mut self,
        path: &[Point2D],
        context: &CollisionContext,
    ) -> SchemeResult<CollisionDetectionResult> {
        // Get adaptive parameters based on channel context
        let params = Self::get_collision_parameters(Some(&context.channel_context));

        let mut result = CollisionDetectionResult {
            has_collisions: false,
            neighbor_collisions: Vec::new(),
            wall_collisions: Vec::new(),
            severity_score: 0.0,
        };

        // Check neighbor collisions with adaptive parameters
        if params.enable_neighbor_detection {
            Self::detect_neighbor_collisions(path, context, &params, &mut result);
        }

        // Check wall collisions with adaptive parameters
        if params.enable_wall_detection {
            Self::detect_wall_collisions(path, context, &params, &mut result);
        }

        // Calculate overall severity with adaptive sensitivity
        result.severity_score = Self::calculate_severity_score(&result, &params);
        result.has_collisions = result.severity_score > 0.0;

        Ok(result)
    }

    /// Detect collisions with neighboring channels
    fn detect_neighbor_collisions(
        path: &[Point2D],
        context: &CollisionContext,
        params: &CollisionParameters,
        result: &mut CollisionDetectionResult,
    ) {
        let min_distance = params.min_channel_distance * params.safety_margin_factor;

        for neighbor in &context.neighbor_info {
            if !neighbor.is_active {
                continue;
            }

            // Check each point in the path against this neighbor
            for (i, &point) in path.iter().enumerate() {
                let distance_to_neighbor = (point.1 - neighbor.y_position).abs();
                let required_distance =
                    min_distance + f64::midpoint(context.current_channel.width, neighbor.width);

                if distance_to_neighbor < required_distance {
                    result.neighbor_collisions.push(NeighborCollision {
                        point_index: i,
                        neighbor_y: neighbor.y_position,
                        actual_distance: distance_to_neighbor,
                        required_distance,
                        severity: (required_distance - distance_to_neighbor) / required_distance,
                    });
                }
            }
        }

        // Function completed successfully
    }

    /// Detect collisions with walls
    fn detect_wall_collisions(
        path: &[Point2D],
        context: &CollisionContext,
        params: &CollisionParameters,
        result: &mut CollisionDetectionResult,
    ) {
        let min_distance = params.min_wall_distance * params.safety_margin_factor;
        let half_width = context.current_channel.width / 2.0;

        for (i, &point) in path.iter().enumerate() {
            // Check distance to each wall
            let distances = [
                point.0 - half_width - context.wall_boundaries.left, // Left wall
                context.wall_boundaries.right - point.0 - half_width, // Right wall
                point.1 - half_width - context.wall_boundaries.bottom, // Bottom wall
                context.wall_boundaries.top - point.1 - half_width,  // Top wall
            ];

            let wall_names = ["left", "right", "bottom", "top"];

            for (wall_idx, &distance) in distances.iter().enumerate() {
                if distance < min_distance {
                    result.wall_collisions.push(WallCollision {
                        point_index: i,
                        wall_name: wall_names[wall_idx].to_string(),
                        actual_distance: distance,
                        required_distance: min_distance,
                        severity: (min_distance - distance) / min_distance,
                    });
                }
            }
        }

        // Function completed successfully
    }

    /// Calculate overall collision severity score
    fn calculate_severity_score(
        result: &CollisionDetectionResult,
        _params: &CollisionParameters,
    ) -> f64 {
        let neighbor_severity: f64 = result.neighbor_collisions.iter().map(|c| c.severity).sum();

        let wall_severity: f64 = result.wall_collisions.iter().map(|c| c.severity).sum();

        neighbor_severity + wall_severity
    }

    /// Apply adaptive collision avoidance to a path
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

    /// Legacy method for backward compatibility
    #[allow(dead_code)] // Part of comprehensive collision detection framework
    fn calculate_reduction_factor(
        detection_result: &CollisionDetectionResult,
        params: &CollisionParameters,
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

        // Scale reduction factor based on severity
        let base_reduction = max_severity * params.max_reduction_factor;
        base_reduction.min(params.max_reduction_factor)
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

        // Function completed successfully
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

        // Function completed successfully
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

        // Function completed successfully
    }

    /// Legacy method for backward compatibility
    #[allow(dead_code)] // Part of comprehensive collision detection framework
    fn apply_path_reduction(
        path: &mut [Point2D],
        context: &CollisionContext,
        reduction_factor: f64,
    ) {
        Self::apply_standard_reduction(
            path,
            context.current_channel.start,
            context.current_channel.end,
            reduction_factor,
        );
    }
}

/// Result of collision detection
#[derive(Debug)]
pub struct CollisionDetectionResult {
    /// Whether any collisions were detected
    pub has_collisions: bool,

    /// Collisions with neighboring channels
    pub neighbor_collisions: Vec<NeighborCollision>,

    /// Collisions with walls
    pub wall_collisions: Vec<WallCollision>,

    /// Overall severity score
    pub severity_score: f64,
}

/// Collision with a neighboring channel
#[derive(Debug)]
pub struct NeighborCollision {
    /// Index of the colliding point in the path
    pub point_index: usize,

    /// Y-coordinate of the neighbor
    pub neighbor_y: f64,

    /// Actual distance to neighbor
    pub actual_distance: f64,

    /// Required minimum distance
    pub required_distance: f64,

    /// Collision severity (0.0 to 1.0)
    pub severity: f64,
}

/// Collision with a wall
#[derive(Debug)]
pub struct WallCollision {
    /// Index of the colliding point in the path
    pub point_index: usize,

    /// Name of the wall (left, right, top, bottom)
    pub wall_name: String,

    /// Actual distance to wall
    pub actual_distance: f64,

    /// Required minimum distance
    pub required_distance: f64,

    /// Collision severity (0.0 to 1.0)
    pub severity: f64,
}

/// Result of collision avoidance application
#[derive(Debug)]
pub struct CollisionAvoidanceResult {
    /// Whether avoidance was applied
    pub applied: bool,

    /// Reduction factor used
    pub reduction_factor: f64,

    /// Original collision severity
    pub original_severity: f64,

    /// Final collision severity after avoidance
    pub final_severity: f64,
}

/// Trait for objects that have collision severity
trait CollisionSeverity {
    fn get_severity(&self) -> f64;
}

impl CollisionSeverity for NeighborCollision {
    fn get_severity(&self) -> f64 {
        self.severity
    }
}

impl CollisionSeverity for WallCollision {
    fn get_severity(&self) -> f64 {
        self.severity
    }
}

impl Default for CollisionDetectionSystem {
    fn default() -> Self {
        Self::new().expect("Failed to create default CollisionDetectionSystem")
    }
}
