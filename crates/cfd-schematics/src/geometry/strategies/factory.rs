//! Factory for creating channel type strategies based on configuration.

use crate::config::{
    constants, ArcConfig, ChannelTypeConfig, ConstantsRegistry, FrustumConfig, SerpentineConfig,
};
use crate::geometry::Point2D;

use super::straight::SmoothTransitionConfig;
use super::{
    ArcChannelStrategy, ChannelTypeStrategy, FrustumChannelStrategy, SerpentineChannelStrategy,
    SmoothStraightChannelStrategy, StraightChannelStrategy,
};

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
                Box::new(super::CustomChannelStrategy::new(channel_type))
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
