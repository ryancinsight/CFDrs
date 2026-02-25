//! Custom channel type strategy wrapper.

use crate::config::GeometryConfig;
use crate::geometry::{ChannelType, Point2D};

use super::ChannelTypeStrategy;

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
    /// use cfd_schematics::geometry::strategies::CustomChannelStrategy;
    /// use cfd_schematics::geometry::ChannelType;
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
