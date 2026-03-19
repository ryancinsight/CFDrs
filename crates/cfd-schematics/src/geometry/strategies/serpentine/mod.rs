//! Serpentine channel generation strategy.

mod amplitude;
mod path;
mod phase;
mod wavelength;

use crate::config::{GeometryConfig, SerpentineConfig};
use crate::geometry::{ChannelType, Point2D};

use super::{ChannelGenerationContext, ChannelTypeStrategy};

/// Space metrics for amplitude calculation
#[derive(Debug, Clone)]
pub(super) struct SpaceMetrics {
    /// Available space for amplitude expansion
    pub available_space: f64,
}

/// Strategy for creating serpentine channels
#[derive(Debug, Clone)]
pub struct SerpentineChannelStrategy {
    pub(super) config: SerpentineConfig,
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
    /// use cfd_schematics::geometry::strategies::SerpentineChannelStrategy;
    /// use cfd_schematics::config::SerpentineConfig;
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
            WaveShape::Sine => (wave_phase + phase_offset).sin(),
            WaveShape::Square => {
                let sine_value = (wave_phase + phase_offset).sin();
                (square_sharpness * sine_value).tanh()
            }
            WaveShape::Triangular => {
                // Triangle wave: 2/pi * asin(sin(phase)) gives a linear
                // ramp between -1 and +1 with sharp apices at the peaks.
                let phase = wave_phase + phase_offset;
                (2.0 / std::f64::consts::PI) * phase.sin().asin()
            }
        }
    }
}
