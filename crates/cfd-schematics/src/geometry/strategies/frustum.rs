//! Frustum (tapered) channel generation strategy.

use crate::config::{FrustumConfig, GeometryConfig};
use crate::geometry::{ChannelType, Point2D};

use super::ChannelTypeStrategy;

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
/// use cfd_schematics::geometry::strategies::FrustumChannelStrategy;
/// use cfd_schematics::config::FrustumConfig;
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
    /// use cfd_schematics::geometry::strategies::FrustumChannelStrategy;
    /// use cfd_schematics::config::FrustumConfig;
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
