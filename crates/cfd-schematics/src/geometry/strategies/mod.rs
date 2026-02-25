//! Channel type generation strategies
//!
//! This module implements the Strategy pattern for channel type generation,
//! providing a clean separation of concerns and enabling easy extension
//! with new channel types while adhering to SOLID principles.

mod arc;
mod custom;
mod envelope;
mod factory;
mod frustum;
mod serpentine;
mod straight;

pub use self::arc::ArcChannelStrategy;
pub use self::custom::CustomChannelStrategy;
pub use self::envelope::{
    AdaptiveGaussianEnvelopeCalculator, EnvelopeCalculator, EnvelopeContext,
    SmoothEndpointEnvelopeCalculator,
};
pub use self::factory::ChannelTypeFactory;
pub use self::frustum::FrustumChannelStrategy;
pub use self::serpentine::SerpentineChannelStrategy;
pub use self::straight::{
    SmoothStraightChannelStrategy, SmoothTransitionConfig, StraightChannelStrategy,
};

use crate::config::GeometryConfig;
use crate::geometry::{ChannelType, Point2D};

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
