//! Enhanced collision detection system with centralized parameter management
//!
//! This module provides a modular collision detection and avoidance system
//! that integrates with the centralized state management system and supports
//! both channel-to-channel and channel-to-wall boundary constraints.

mod avoidance;

use crate::{
    config::ConstantsRegistry,
    error::SchemeResult,
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
pub struct CollisionDetectionSystem {}

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
        Ok(Self {})
    }

    /// Create with existing parameter registry
    #[must_use]
    pub fn with_registry(_registry: ParameterRegistry) -> Self {
        Self {}
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
