//! Channel type generation strategies
//!
//! This module implements the Strategy pattern for channel type generation,
//! providing a clean separation of concerns and enabling easy extension
//! with new channel types while adhering to SOLID principles.

use crate::config::{
    constants, ArcConfig, ChannelTypeConfig, FrustumConfig, GeometryConfig, SerpentineConfig,
};
use crate::config_constants::ConstantsRegistry;
use crate::geometry::optimization::optimize_serpentine_parameters;
use crate::geometry::{ChannelType, Point2D};
use crate::state_management::bilateral_symmetry::{
    BilateralPhaseDirectionCalculator, BilateralSymmetryConfig, SymmetryContext,
};

/// Context object for channel generation to reduce parameter coupling
///
/// This struct groups related parameters together, following the
/// Parameter Object pattern to improve method signatures and reduce coupling.
#[derive(Debug, Clone)]
pub struct ChannelGenerationContext<'a> {
    /// Geometry configuration parameters
    pub geometry_config: &'a GeometryConfig,
    /// Bounding box dimensions (width, height)
    pub box_dims: (f64, f64),
    /// Total number of branches in the system
    pub total_branches: usize,
    /// Information about neighboring channels for collision avoidance
    pub neighbor_info: Option<&'a [f64]>,
}

/// Space metrics for amplitude calculation
#[derive(Debug, Clone)]
struct SpaceMetrics {
    /// Available space for amplitude expansion
    available_space: f64,
}

impl<'a> ChannelGenerationContext<'a> {
    /// Create a new channel generation context
    #[must_use]
    pub const fn new(
        geometry_config: &'a GeometryConfig,
        box_dims: (f64, f64),
        total_branches: usize,
        neighbor_info: Option<&'a [f64]>,
    ) -> Self {
        Self {
            geometry_config,
            box_dims,
            total_branches,
            neighbor_info,
        }
    }
}

/// Trait for calculating envelope functions for serpentine channels
///
/// This trait abstracts envelope calculation logic to eliminate code duplication
/// and provide a clean interface for different envelope types.
pub trait EnvelopeCalculator {
    /// Calculate envelope value at parameter t (0.0 to 1.0)
    fn calculate_envelope(&self, t: f64, context: &EnvelopeContext) -> f64;
}

/// Context for envelope calculations
#[derive(Debug, Clone)]
pub struct EnvelopeContext {
    /// Channel length
    pub channel_length: f64,
    /// Channel direction vector (dx, dy)
    pub direction: (f64, f64),
    /// Distance to nearest node
    pub node_distance: f64,
    /// Adaptive configuration
    pub adaptive_config: crate::config::AdaptiveSerpentineConfig,
    /// Gaussian width factor
    pub gaussian_width_factor: f64,
}

/// Smooth endpoint envelope calculator
///
/// Provides smooth transitions at channel endpoints using smoothstep function.
pub struct SmoothEndpointEnvelopeCalculator;

impl EnvelopeCalculator for SmoothEndpointEnvelopeCalculator {
    fn calculate_envelope(&self, t: f64, _context: &EnvelopeContext) -> f64 {
        // Smooth endpoint envelope using smoothstep function
        // This ensures zero amplitude and zero derivative at endpoints
        let smoothstep = |x: f64| x * x * 2.0f64.mul_add(-x, 3.0);

        let constants = ConstantsRegistry::new();
        let start_threshold = constants.get_smooth_endpoint_start_threshold();
        let end_threshold = constants.get_smooth_endpoint_end_threshold();

        // Create smooth transitions at both ends
        let start_transition = if t < start_threshold {
            smoothstep(t / start_threshold)
        } else {
            1.0
        };

        let end_transition = if t > end_threshold {
            smoothstep((1.0 - t) / (1.0 - end_threshold))
        } else {
            1.0
        };

        start_transition * end_transition
    }
}

/// Gaussian envelope calculator with adaptive behavior
///
/// Provides Gaussian-shaped envelope with adaptive parameters based on
/// channel characteristics and proximity to nodes/walls.
pub struct AdaptiveGaussianEnvelopeCalculator;

impl EnvelopeCalculator for AdaptiveGaussianEnvelopeCalculator {
    fn calculate_envelope(&self, t: f64, context: &EnvelopeContext) -> f64 {
        let dx = context.direction.0;
        let dy = context.direction.1;
        let channel_length = context.channel_length;
        let node_distance = context.node_distance;

        // Determine if this is primarily a horizontal channel
        let is_horizontal = dx.abs() > dy.abs();
        let horizontal_ratio = dx.abs() / node_distance;

        // For horizontal channels (middle sections), we want less aggressive tapering
        let middle_section_factor = if is_horizontal
            && horizontal_ratio > context.adaptive_config.horizontal_ratio_threshold
        {
            (1.0 - context.adaptive_config.middle_section_amplitude_factor).mul_add(
                horizontal_ratio,
                context.adaptive_config.middle_section_amplitude_factor,
            )
        } else {
            1.0
        };

        // Distance-based normalization
        let distance_normalization = if context.adaptive_config.enable_distance_based_scaling {
            (node_distance / context.adaptive_config.node_distance_normalization)
                .clamp(0.1, 1.0)
        } else {
            1.0
        };

        // Calculate effective sigma based on distance and section type
        let base_sigma = channel_length / context.gaussian_width_factor;
        let effective_sigma = base_sigma * distance_normalization * middle_section_factor;

        // Center the envelope
        let center = 0.5;

        // Create smooth dome-shaped envelope instead of sharp Gaussian peaks
        // This prevents self-intersection when channels curve back on themselves
        let dome_envelope = if (t - center).abs() < 0.45 {
            // Use raised cosine for the main dome (much smoother than Gaussian)
            let normalized_t = (t - center) / 0.45; // Scale to [-1, 1] range
            let cosine_factor = 0.5 * (1.0 + (std::f64::consts::PI * normalized_t).cos());

            // Apply effective sigma scaling to the dome
            let sigma_factor = (effective_sigma / channel_length).clamp(0.1, 0.3);
            let dome_width = 0.45 * sigma_factor / 0.2; // Scale dome width based on sigma

            if (t - center).abs() < dome_width {
                let dome_t = (t - center) / dome_width;
                0.5 * (1.0 + (std::f64::consts::PI * dome_t).cos())
            } else {
                // Smooth transition to zero
                let transition_factor = ((t - center).abs() - dome_width) / (0.45 - dome_width);
                let smoothstep = (transition_factor * transition_factor)
                    .mul_add(-2.0f64.mul_add(-transition_factor, 3.0), 1.0);
                cosine_factor * smoothstep * 0.1
            }
        } else {
            // Smooth transition to zero at edges using smoothstep
            let edge_distance = ((t - center).abs() - 0.45) / 0.05; // 0.05 is transition zone
            if (0.0..1.0).contains(&edge_distance) {
                let smoothstep = (edge_distance * edge_distance)
                    .mul_add(-2.0f64.mul_add(-edge_distance, 3.0), 1.0);
                smoothstep * 0.05 // Very small amplitude at edges
            } else {
                0.0
            }
        };

        // For middle sections, enhance the dome but keep it smooth
        if is_horizontal && horizontal_ratio > context.adaptive_config.horizontal_ratio_threshold {
            let plateau_width = context.adaptive_config.plateau_width_factor.min(0.3); // Limit plateau width
            let plateau_start = 0.5 - plateau_width / 2.0;
            let plateau_end = 0.5 + plateau_width / 2.0;

            if t >= plateau_start && t <= plateau_end {
                // In the plateau region, enhance the dome but keep it smooth
                let plateau_factor = 1.0 - ((t - 0.5).abs() / (plateau_width / 2.0));
                let enhanced_amplitude = (1.0 - context.adaptive_config.plateau_amplitude_factor)
                    .mul_add(
                        plateau_factor,
                        context.adaptive_config.plateau_amplitude_factor,
                    );
                dome_envelope.max(enhanced_amplitude * dome_envelope)
            } else {
                dome_envelope
            }
        } else {
            dome_envelope
        }
    }
}

/// Strategy trait for generating different types of channels
///
/// This trait follows the Strategy pattern to enable different channel
/// generation algorithms while maintaining a consistent interface.
/// Each strategy is responsible for creating a specific type of channel
/// based on the provided points and configuration.
pub trait ChannelTypeStrategy {
    /// Create a channel type between two points
    ///
    /// # Arguments
    /// * `from` - Starting point of the channel
    /// * `to` - Ending point of the channel
    /// * `geometry_config` - General geometry configuration
    /// * `box_dims` - Dimensions of the containing box
    /// * `total_branches` - Total number of branches in the system (for scaling)
    /// * `neighbor_info` - Optional information about neighboring channels
    ///
    /// # Returns
    /// A `ChannelType` representing the generated channel
    fn create_channel(
        &self,
        from: Point2D,
        to: Point2D,
        geometry_config: &GeometryConfig,
        box_dims: (f64, f64),
        total_branches: usize,
        neighbor_info: Option<&[f64]>,
    ) -> ChannelType;
}

/// Strategy for creating straight channels
#[derive(Debug, Clone)]
pub struct StraightChannelStrategy;

impl ChannelTypeStrategy for StraightChannelStrategy {
    fn create_channel(
        &self,
        _from: Point2D,
        _to: Point2D,
        _geometry_config: &GeometryConfig,
        _box_dims: (f64, f64),
        _total_branches: usize,
        _neighbor_info: Option<&[f64]>,
    ) -> ChannelType {
        ChannelType::Straight
    }
}

/// Strategy for creating smooth straight channels with transition zones
#[derive(Debug, Clone)]
pub struct SmoothStraightChannelStrategy {
    /// Configuration for transition zones
    pub transition_config: SmoothTransitionConfig,
}

/// Configuration for smooth transition zones in straight channels
///
/// This configuration controls the smooth transition zones at the endpoints
/// of straight channels to eliminate sharp corners when connecting to other
/// channel types.
///
/// # Examples
///
/// ```rust
/// use scheme::geometry::strategies::SmoothTransitionConfig;
///
/// // Create with default values
/// let config = SmoothTransitionConfig::default();
///
/// // Create with custom values for subtle transitions
/// let subtle = SmoothTransitionConfig {
///     transition_length_factor: 0.1,
///     transition_amplitude_factor: 0.2,
///     transition_smoothness: 15,
///     wave_multiplier: 1.5,
/// };
/// ```
#[derive(Debug, Clone, Copy)]
pub struct SmoothTransitionConfig {
    /// Length of transition zone as fraction of channel length (0.0 to 0.5)
    pub transition_length_factor: f64,
    /// Maximum amplitude of transition waves relative to channel width (0.0 to 1.0)
    pub transition_amplitude_factor: f64,
    /// Number of points to use for transition smoothing (5 to 100)
    pub transition_smoothness: usize,
    /// Wave multiplier for transition waves (0.5 to 10.0, where 2.0 = one complete wave)
    pub wave_multiplier: f64,
}

impl Default for SmoothTransitionConfig {
    fn default() -> Self {
        let constants = ConstantsRegistry::new();
        Self {
            transition_length_factor: constants.get_default_transition_length_factor(),
            transition_amplitude_factor: constants.get_default_transition_amplitude_factor(),
            transition_smoothness: constants.get_default_transition_smoothness(),
            wave_multiplier: constants.get_default_wave_multiplier(),
        }
    }
}

impl SmoothTransitionConfig {
    /// Create a new smooth transition configuration with validation
    pub fn new(
        transition_length_factor: f64,
        transition_amplitude_factor: f64,
        transition_smoothness: usize,
        wave_multiplier: f64,
    ) -> Result<Self, String> {
        let config = Self {
            transition_length_factor,
            transition_amplitude_factor,
            transition_smoothness,
            wave_multiplier,
        };
        config.validate()?;
        Ok(config)
    }

    /// Validate the configuration parameters
    pub fn validate(&self) -> Result<(), String> {
        if self.transition_length_factor < 0.0 || self.transition_length_factor > 0.5 {
            return Err("transition_length_factor must be between 0.0 and 0.5".to_string());
        }

        if self.transition_amplitude_factor < 0.0 || self.transition_amplitude_factor > 1.0 {
            return Err("transition_amplitude_factor must be between 0.0 and 1.0".to_string());
        }

        if self.transition_smoothness < 5 || self.transition_smoothness > 100 {
            return Err("transition_smoothness must be between 5 and 100".to_string());
        }

        if self.wave_multiplier < 0.5 || self.wave_multiplier > 10.0 {
            return Err("wave_multiplier must be between 0.5 and 10.0".to_string());
        }

        Ok(())
    }

    /// Create a subtle transition configuration
    #[must_use]
    pub const fn subtle() -> Self {
        Self {
            transition_length_factor: 0.1,
            transition_amplitude_factor: 0.2,
            transition_smoothness: 15,
            wave_multiplier: 1.5,
        }
    }

    /// Create a pronounced transition configuration
    #[must_use]
    pub const fn pronounced() -> Self {
        Self {
            transition_length_factor: 0.25,
            transition_amplitude_factor: 0.5,
            transition_smoothness: 30,
            wave_multiplier: 3.0,
        }
    }

    /// Create a high-quality transition configuration for detailed work
    #[must_use]
    pub const fn high_quality() -> Self {
        Self {
            transition_length_factor: 0.2,
            transition_amplitude_factor: 0.4,
            transition_smoothness: 50,
            wave_multiplier: 2.0,
        }
    }

    /// Create a fast transition configuration for quick generation
    #[must_use]
    pub const fn fast() -> Self {
        Self {
            transition_length_factor: 0.15,
            transition_amplitude_factor: 0.3,
            transition_smoothness: 10,
            wave_multiplier: 2.0,
        }
    }
}

impl SmoothStraightChannelStrategy {
    #[must_use]
    pub const fn new(config: SmoothTransitionConfig) -> Self {
        Self {
            transition_config: config,
        }
    }
}

impl ChannelTypeStrategy for SmoothStraightChannelStrategy {
    fn create_channel(
        &self,
        from: Point2D,
        to: Point2D,
        geometry_config: &GeometryConfig,
        _box_dims: (f64, f64),
        _total_branches: usize,
        _neighbor_info: Option<&[f64]>,
    ) -> ChannelType {
        let path = self.generate_smooth_straight_path(from, to, geometry_config);
        ChannelType::SmoothStraight { path }
    }
}

impl SmoothStraightChannelStrategy {
    /// Generate a smooth straight path with transition zones at endpoints
    fn generate_smooth_straight_path(
        &self,
        p1: Point2D,
        p2: Point2D,
        geometry_config: &GeometryConfig,
    ) -> Vec<Point2D> {
        let dx = p2.0 - p1.0;
        let dy = p2.1 - p1.1;
        let channel_length = dx.hypot(dy);

        let constants = ConstantsRegistry::new();
        // For very short channels, just return straight line
        if channel_length
            < geometry_config.channel_width * constants.get_short_channel_width_multiplier()
        {
            return vec![p1, p2];
        }

        let transition_length = channel_length * self.transition_config.transition_length_factor;
        let max_amplitude =
            geometry_config.channel_width * self.transition_config.transition_amplitude_factor;

        // Calculate total points: transition + middle + transition
        let transition_points = self.transition_config.transition_smoothness;
        let middle_points = geometry_config.generation.smooth_straight_middle_points;
        let total_points = transition_points * 2 + middle_points;

        let mut path = Vec::with_capacity(total_points);

        // Perpendicular direction for wave displacement
        let perp_x = -dy / channel_length;
        let perp_y = dx / channel_length;

        for i in 0..total_points {
            let t = i as f64 / (total_points - 1) as f64;

            // Base position along the line
            let base_x = t.mul_add(dx, p1.0);
            let base_y = t.mul_add(dy, p1.1);

            // Calculate smooth transition amplitude
            let amplitude = self.calculate_transition_amplitude(
                t,
                transition_length / channel_length,
                max_amplitude,
            );

            // Apply small wave for smooth transition
            let wave_phase = std::f64::consts::PI * self.transition_config.wave_multiplier * t;
            let wave_amplitude = amplitude * wave_phase.sin();

            let x = wave_amplitude.mul_add(perp_x, base_x);
            let y = wave_amplitude.mul_add(perp_y, base_y);

            // Ensure exact endpoints
            if i == 0 {
                path.push(p1);
            } else if i == total_points - 1 {
                path.push(p2);
            } else {
                path.push((x, y));
            }
        }

        path
    }

    /// Calculate transition amplitude that smoothly goes to zero at endpoints
    fn calculate_transition_amplitude(
        &self,
        t: f64,
        transition_factor: f64,
        max_amplitude: f64,
    ) -> f64 {
        // Create smooth transitions at both ends
        let start_transition = if t < transition_factor {
            // Smooth ramp up from 0 to 1 using smoothstep
            let local_t = t / transition_factor;
            local_t * local_t * 2.0f64.mul_add(-local_t, 3.0)
        } else {
            1.0
        };

        let end_transition = if t > (1.0 - transition_factor) {
            // Smooth ramp down from 1 to 0 using smoothstep
            let local_t = (1.0 - t) / transition_factor;
            local_t * local_t * 2.0f64.mul_add(-local_t, 3.0)
        } else {
            1.0
        };

        max_amplitude * start_transition * end_transition
    }
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
    /// use scheme::geometry::strategies::SerpentineChannelStrategy;
    /// use scheme::config::SerpentineConfig;
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
                .min_by(|a, b| a.total_cmp(b))
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

        // Gaussian envelope parameters
        // Note: sigma and center are now calculated in the improved envelope function

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
    /// use scheme::geometry::strategies::ArcChannelStrategy;
    /// use scheme::config::ArcConfig;
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
                .min_by(|a, b| a.total_cmp(b))
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

        let (length, height) = box_dims;
        let _center_x = length / 2.0;
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

/// Factory for creating channel type strategies based on configuration
///
/// This factory implements the Factory pattern to create appropriate
/// channel type strategies based on the provided configuration.
/// It encapsulates the logic for determining which strategy to use
/// and handles complex configurations like Adaptive and `MixedByPosition`.
pub struct ChannelTypeFactory;

impl ChannelTypeFactory {
    /// Create a strategy based on the channel type configuration
    ///
    /// # Arguments
    /// * `config` - The channel type configuration
    /// * `from` - Starting point of the channel (for context-aware strategies)
    /// * `to` - Ending point of the channel (for context-aware strategies)
    /// * `box_dims` - Dimensions of the containing box (for context-aware strategies)
    ///
    /// # Returns
    /// A boxed trait object implementing `ChannelTypeStrategy`
    #[must_use]
    pub fn create_strategy(
        config: &ChannelTypeConfig,
        from: Point2D,
        to: Point2D,
        box_dims: (f64, f64),
    ) -> Box<dyn ChannelTypeStrategy> {
        match config {
            ChannelTypeConfig::AllStraight => Box::new(StraightChannelStrategy),

            ChannelTypeConfig::AllSmoothStraight(smooth_config) => {
                Box::new(SmoothStraightChannelStrategy::new(*smooth_config))
            }

            ChannelTypeConfig::AllSerpentine(serpentine_config) => {
                Box::new(SerpentineChannelStrategy::new(*serpentine_config))
            }

            ChannelTypeConfig::AllArcs(arc_config) => {
                Box::new(ArcChannelStrategy::new(*arc_config))
            }

            ChannelTypeConfig::AllFrustum(frustum_config) => {
                Box::new(FrustumChannelStrategy::new(*frustum_config))
            }

            ChannelTypeConfig::MixedByPosition {
                middle_zone_fraction,
                serpentine_config,
                arc_config,
            } => {
                let (length, _) = box_dims;
                let mid_x = length / 2.0;
                let channel_mid_x = f64::midpoint(from.0, to.0);
                let tolerance = length * middle_zone_fraction / 2.0;

                if (channel_mid_x - mid_x).abs() < tolerance {
                    Box::new(SerpentineChannelStrategy::new(*serpentine_config))
                } else if Self::is_angled_channel(from, to) {
                    Box::new(ArcChannelStrategy::new(*arc_config))
                } else {
                    Box::new(StraightChannelStrategy)
                }
            }

            ChannelTypeConfig::Adaptive {
                serpentine_config,
                arc_config,
                frustum_config,
            } => Self::create_adaptive_strategy(
                from,
                to,
                box_dims,
                *serpentine_config,
                *arc_config,
                *frustum_config,
            ),

            ChannelTypeConfig::SmoothSerpentineWithTransitions {
                serpentine_config,
                smooth_straight_config,
            } => Self::create_smooth_serpentine_strategy(
                from,
                to,
                box_dims,
                *serpentine_config,
                *smooth_straight_config,
            ),

            ChannelTypeConfig::Custom(func) => {
                // For custom functions, we create a wrapper strategy
                let channel_type = func(from, to, box_dims);
                Box::new(CustomChannelStrategy::new(channel_type))
            }
        }
    }

    /// Create an adaptive strategy based on channel characteristics
    fn create_adaptive_strategy(
        from: Point2D,
        to: Point2D,
        box_dims: (f64, f64),
        serpentine_config: SerpentineConfig,
        arc_config: ArcConfig,
        frustum_config: FrustumConfig,
    ) -> Box<dyn ChannelTypeStrategy> {
        let constants = ConstantsRegistry::new();
        let dx = to.0 - from.0;
        let dy = to.1 - from.1;
        let length = dx.hypot(dy);

        // Smart logic: use frustum for medium-length channels that could benefit from flow control,
        // serpentine for longer horizontal channels, arcs for angled channels, straight for short channels
        if length > box_dims.0 * constants.get_long_horizontal_threshold()
            && dy.abs() < dx.abs() * constants.get_horizontal_angle_threshold()
        {
            // Long horizontal channel - use serpentine for mixing
            Box::new(SerpentineChannelStrategy::new(serpentine_config))
        } else if length > box_dims.0 * constants.get_frustum_min_length_threshold()
            && length < box_dims.0 * constants.get_frustum_max_length_threshold()
            && dy.abs() < dx.abs() * constants.get_frustum_angle_threshold()
        {
            // Medium-length horizontal channel - use frustum for flow control
            Box::new(FrustumChannelStrategy::new(frustum_config))
        } else if Self::is_angled_channel(from, to)
            && length > box_dims.0 * constants.get_min_arc_length_threshold()
        {
            // Angled channel of reasonable length - use arc
            Box::new(ArcChannelStrategy::new(arc_config))
        } else {
            // Default to straight for short channels
            Box::new(StraightChannelStrategy)
        }
    }

    /// Create a smooth serpentine strategy with smooth straight junction connectors
    fn create_smooth_serpentine_strategy(
        from: Point2D,
        to: Point2D,
        box_dims: (f64, f64),
        serpentine_config: SerpentineConfig,
        smooth_straight_config: SmoothTransitionConfig,
    ) -> Box<dyn ChannelTypeStrategy> {
        let constants = ConstantsRegistry::new();
        let dx = to.0 - from.0;
        let dy = to.1 - from.1;
        let length = dx.hypot(dy);

        // Use serpentine for longer horizontal channels (branches)
        // Use smooth straight for shorter channels and junction connectors
        if length > box_dims.0 * constants.get_long_horizontal_threshold()
            && dy.abs() < dx.abs() * constants.get_horizontal_angle_threshold()
        {
            // Long horizontal channel - use serpentine
            Box::new(SerpentineChannelStrategy::new(serpentine_config))
        } else {
            // Junction connectors and short channels - use smooth straight
            Box::new(SmoothStraightChannelStrategy::new(smooth_straight_config))
        }
    }

    /// Check if a channel is significantly angled
    fn is_angled_channel(from: Point2D, to: Point2D) -> bool {
        let constants = ConstantsRegistry::new();
        let dx = to.0 - from.0;
        let dy = to.1 - from.1;

        if dx.abs() < constants.get_geometric_tolerance() {
            return dy.abs() > constants.get_geometric_tolerance(); // Vertical channel
        }

        let slope = dy / dx;
        slope.abs() > constants::strategy_thresholds::ANGLED_CHANNEL_SLOPE_THRESHOLD
    }
}

/// Strategy wrapper for custom channel type functions
///
/// This strategy wraps custom channel type functions to fit into
/// the strategy pattern while maintaining backward compatibility.
#[derive(Debug, Clone)]
pub struct CustomChannelStrategy {
    channel_type: ChannelType,
}

impl CustomChannelStrategy {
    /// Create a new direct channel strategy with the given channel type
    ///
    /// This strategy directly uses the provided channel type without any
    /// configuration-based selection logic.
    ///
    /// # Arguments
    ///
    /// * `channel_type` - The specific channel type to use for all channels
    ///
    /// # Examples
    ///
    /// ```rust
    /// use scheme::geometry::strategies::CustomChannelStrategy;
    /// use scheme::geometry::ChannelType;
    ///
    /// let strategy = CustomChannelStrategy::new(ChannelType::Straight);
    /// ```
    #[must_use]
    pub const fn new(channel_type: ChannelType) -> Self {
        Self { channel_type }
    }
}

impl ChannelTypeStrategy for CustomChannelStrategy {
    fn create_channel(
        &self,
        _from: Point2D,
        _to: Point2D,
        _geometry_config: &GeometryConfig,
        _box_dims: (f64, f64),
        _total_branches: usize,
        _neighbor_info: Option<&[f64]>,
    ) -> ChannelType {
        self.channel_type.clone()
    }
}

/// Strategy for generating frustum (tapered) channels with venturi throat functionality
///
/// This strategy creates channels with variable width along their length, featuring:
/// - Configurable inlet, throat, and outlet widths
/// - Multiple taper profiles (linear, exponential, smooth)
/// - Precise throat positioning
/// - Smooth width transitions
///
/// The frustum channel is ideal for applications requiring flow acceleration/deceleration,
/// such as venturi throats, flow focusing, or particle sorting.
///
/// # Design Principles
/// - **Single Responsibility**: Focused solely on frustum channel generation
/// - **Open/Closed**: Extensible for new taper profiles without modification
/// - **Dependency Inversion**: Depends on abstractions (`ChannelTypeStrategy` trait)
/// - **DRY**: Reuses common path generation patterns
/// - **KISS**: Simple, clear implementation
///
/// # Examples
///
/// ```rust
/// use scheme::geometry::strategies::FrustumChannelStrategy;
/// use scheme::config::FrustumConfig;
///
/// let config = FrustumConfig::default();
/// let strategy = FrustumChannelStrategy::new(config);
/// ```
#[derive(Debug, Clone)]
pub struct FrustumChannelStrategy {
    /// Configuration parameters for the frustum channel
    config: FrustumConfig,
}

impl FrustumChannelStrategy {
    /// Create a new frustum channel strategy with the specified configuration
    ///
    /// # Arguments
    /// * `config` - Configuration parameters for the frustum channel
    ///
    /// # Returns
    /// * `Self` - A new frustum channel strategy instance
    ///
    /// # Examples
    ///
    /// ```rust
    /// use scheme::geometry::strategies::FrustumChannelStrategy;
    /// use scheme::config::FrustumConfig;
    ///
    /// let config = FrustumConfig::default();
    /// let strategy = FrustumChannelStrategy::new(config);
    /// ```
    #[must_use]
    pub const fn new(config: FrustumConfig) -> Self {
        Self { config }
    }

    /// Generate the centerline path for the frustum channel using iterator combinators
    ///
    /// Creates a straight line path from start to end point with the specified
    /// number of points for smooth width transitions.
    ///
    /// # Arguments
    /// * `from` - Starting point of the channel
    /// * `to` - Ending point of the channel
    ///
    /// # Returns
    /// * `Vec<Point2D>` - The centerline path points
    fn generate_centerline_path(&self, from: Point2D, to: Point2D) -> Vec<Point2D> {
        let dx = to.0 - from.0;
        let dy = to.1 - from.1;
        let inv_smoothness = 1.0 / (self.config.smoothness - 1) as f64;

        // Use iterator combinators for zero-copy generation
        (0..self.config.smoothness)
            .map(|i| {
                let t = i as f64 * inv_smoothness;
                (dx.mul_add(t, from.0), dy.mul_add(t, from.1))
            })
            .collect()
    }

    /// Generate width values corresponding to each point in the path
    ///
    /// Calculates the width at each point along the channel using the
    /// configured taper profile.
    ///
    /// # Arguments
    /// * `path_length` - Number of points in the path
    ///
    /// # Returns
    /// * `Vec<f64>` - Width values for each path point
    fn generate_width_profile(&self, path_length: usize) -> Vec<f64> {
        let mut widths = Vec::with_capacity(path_length);

        for i in 0..path_length {
            let t = i as f64 / (path_length - 1) as f64;
            let width = self.config.width_at_position(t);
            widths.push(width);
        }

        widths
    }
}

impl ChannelTypeStrategy for FrustumChannelStrategy {
    fn create_channel(
        &self,
        from: Point2D,
        to: Point2D,
        _geometry_config: &GeometryConfig,
        _box_dims: (f64, f64),
        _total_branches: usize,
        _neighbor_info: Option<&[f64]>,
    ) -> ChannelType {
        // Generate the centerline path
        let path = self.generate_centerline_path(from, to);

        // Generate the width profile
        let widths = self.generate_width_profile(path.len());

        ChannelType::Frustum {
            path,
            widths,
            inlet_width: self.config.inlet_width,
            throat_width: self.config.throat_width,
            outlet_width: self.config.outlet_width,
            taper_profile: self.config.taper_profile,
            throat_position: self.config.throat_position,
            has_venturi_throat: true,
        }
    }
}
