use crate::config::channel::arc::ArcConfig;
use crate::config::channel::frustum::FrustumConfig;
use crate::config::channel::serpentine::SerpentineConfig;
use crate::config::constants::primitives as constants;
use crate::geometry::strategies::SmoothTransitionConfig;
use crate::geometry::ChannelType;

/// Configuration for selecting channel types in microfluidic schematics
///
/// This enum provides different strategies for determining what type of channel
/// (straight, serpentine, or arc) to use for each connection in the system.
///
/// # Examples
///
/// ```rust
/// use cfd_schematics::config::{ChannelTypeConfig, SerpentineConfig, ArcConfig, FrustumConfig};
///
/// // All channels will be straight lines
/// let straight_config = ChannelTypeConfig::AllStraight;
///
/// // All channels will be serpentine with default parameters
/// let serpentine_config = ChannelTypeConfig::AllSerpentine(SerpentineConfig::default());
///
/// // All channels will be frustum (tapered) with default parameters
/// let frustum_config = ChannelTypeConfig::AllFrustum(FrustumConfig::default());
///
/// // Adaptive selection based on channel characteristics
/// let adaptive_config = ChannelTypeConfig::Adaptive {
///     serpentine_config: SerpentineConfig::default(),
///     arc_config: ArcConfig::default(),
///     frustum_config: FrustumConfig::default(),
/// };
/// ```
#[derive(Debug, Clone, Copy)]
pub enum ChannelTypeConfig {
    /// All channels will be straight lines
    AllStraight,
    /// All channels will be smooth straight lines with transition zones
    AllSmoothStraight(SmoothTransitionConfig),
    /// All channels will be serpentine with the specified configuration
    AllSerpentine(SerpentineConfig),
    /// All channels will be arcs with the specified configuration
    AllArcs(ArcConfig),
    /// All channels will be frustum (tapered) with the specified configuration
    AllFrustum(FrustumConfig),
    /// Channels are selected based on their position in the layout
    MixedByPosition {
        /// Fraction of the box width that defines the middle zone for serpentine channels (0.0 to 1.0)
        middle_zone_fraction: f64,
        /// Configuration for serpentine channels in the middle zone
        serpentine_config: SerpentineConfig,
        /// Configuration for arc channels outside the middle zone
        arc_config: ArcConfig,
    },
    /// Adaptive channel type selection based on channel characteristics
    Adaptive {
        /// Configuration for serpentine channels when selected by adaptive algorithm
        serpentine_config: SerpentineConfig,
        /// Configuration for arc channels when selected by adaptive algorithm
        arc_config: ArcConfig,
        /// Configuration for frustum channels when selected by adaptive algorithm
        frustum_config: FrustumConfig,
    },
    /// Smooth serpentine channels with smooth straight junction connectors
    SmoothSerpentineWithTransitions {
        /// Configuration for serpentine channels in branches
        serpentine_config: SerpentineConfig,
        /// Configuration for smooth straight channels in junction connectors
        smooth_straight_config: SmoothTransitionConfig,
    },
    /// Custom function for determining channel type based on endpoints and box dimensions
    #[allow(clippy::type_complexity)]
    Custom(fn(from: (f64, f64), to: (f64, f64), box_dims: (f64, f64)) -> ChannelType),
}

impl Default for ChannelTypeConfig {
    fn default() -> Self {
        Self::Adaptive {
            serpentine_config: SerpentineConfig::default(),
            arc_config: ArcConfig::default(),
            frustum_config: FrustumConfig::default(),
        }
    }
}

/// Builder for complex channel type configurations
pub struct ChannelTypeConfigBuilder {
    serpentine_config: SerpentineConfig,
    arc_config: ArcConfig,
    frustum_config: FrustumConfig,
    middle_zone_fraction: f64,
}

impl ChannelTypeConfigBuilder {
    /// Create a new builder with default configurations
    #[must_use]
    pub fn new() -> Self {
        Self {
            serpentine_config: SerpentineConfig::default(),
            arc_config: ArcConfig::default(),
            frustum_config: FrustumConfig::default(),
            middle_zone_fraction: constants::strategy_thresholds::DEFAULT_MIDDLE_ZONE_FRACTION,
        }
    }

    /// Set the serpentine configuration
    #[must_use]
    pub const fn with_serpentine_config(mut self, config: SerpentineConfig) -> Self {
        self.serpentine_config = config;
        self
    }

    /// Set the arc configuration
    #[must_use]
    pub const fn with_arc_config(mut self, config: ArcConfig) -> Self {
        self.arc_config = config;
        self
    }

    /// Set the frustum configuration
    #[must_use]
    pub const fn with_frustum_config(mut self, config: FrustumConfig) -> Self {
        self.frustum_config = config;
        self
    }

    /// Set the middle zone fraction for mixed by position
    #[must_use]
    pub const fn with_middle_zone_fraction(mut self, fraction: f64) -> Self {
        self.middle_zone_fraction = fraction;
        self
    }

    /// Build an adaptive channel type configuration
    #[must_use]
    pub const fn build_adaptive(self) -> ChannelTypeConfig {
        ChannelTypeConfig::Adaptive {
            serpentine_config: self.serpentine_config,
            arc_config: self.arc_config,
            frustum_config: self.frustum_config,
        }
    }

    /// Build a mixed by position channel type configuration
    #[must_use]
    pub const fn build_mixed_by_position(self) -> ChannelTypeConfig {
        ChannelTypeConfig::MixedByPosition {
            middle_zone_fraction: self.middle_zone_fraction,
            serpentine_config: self.serpentine_config,
            arc_config: self.arc_config,
        }
    }
}

impl Default for ChannelTypeConfigBuilder {
    fn default() -> Self {
        Self::new()
    }
}
