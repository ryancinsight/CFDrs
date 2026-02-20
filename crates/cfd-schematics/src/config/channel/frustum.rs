use crate::error::{ConfigurationError, ConfigurationResult};
use serde::{Deserialize, Serialize};

/// Taper profile types for frustum channels
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum TaperProfile {
    /// Linear taper - constant rate of width change
    #[default]
    Linear,
    /// Exponential taper - smooth exponential transition
    Exponential,
    /// Smooth taper - uses cosine function for very smooth transitions
    Smooth,
}

/// Configuration for frustum (tapered) channels with venturi throat functionality
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct FrustumConfig {
    /// Inlet width (starting width) - must be positive (0.1 to 50.0)
    pub inlet_width: f64,
    /// Throat width (minimum width at center) - must be less than inlet/outlet widths (0.05 to 25.0)
    pub throat_width: f64,
    /// Outlet width (ending width) - must be positive (0.1 to 50.0)
    pub outlet_width: f64,
    /// Taper profile type - controls how width changes along the channel
    pub taper_profile: TaperProfile,
    /// Number of points to generate along the path - higher = smoother (10 to 1000)
    pub smoothness: usize,
    /// Throat position factor - 0.5 = center, 0.0 = at inlet, 1.0 = at outlet (0.1 to 0.9)
    pub throat_position: f64,
}

impl Default for FrustumConfig {
    fn default() -> Self {
        Self {
            inlet_width: 2.0,
            throat_width: 0.5,
            outlet_width: 2.0,
            taper_profile: TaperProfile::Linear,
            smoothness: 50,
            throat_position: 0.5, // Center
        }
    }
}

impl FrustumConfig {
    /// Create a new `FrustumConfig` with validation
    ///
    /// # Arguments
    /// * `inlet_width` - Starting width of the channel (must be positive)
    /// * `throat_width` - Minimum width at the throat (must be less than inlet/outlet)
    /// * `outlet_width` - Ending width of the channel (must be positive)
    /// * `taper_profile` - Type of taper profile to use
    /// * `smoothness` - Number of points for path generation
    /// * `throat_position` - Position of throat along channel (0.5 = center)
    ///
    /// # Returns
    /// * `ConfigurationResult<Self>` - The validated configuration or an error
    pub fn new(
        inlet_width: f64,
        throat_width: f64,
        outlet_width: f64,
        taper_profile: TaperProfile,
        smoothness: usize,
        throat_position: f64,
    ) -> ConfigurationResult<Self> {
        let config = Self {
            inlet_width,
            throat_width,
            outlet_width,
            taper_profile,
            smoothness,
            throat_position,
        };
        config.validate()?;
        Ok(config)
    }

    /// Validate the frustum configuration
    pub fn validate(&self) -> ConfigurationResult<()> {
        // Validate inlet width
        if self.inlet_width <= 0.0 || self.inlet_width > 50.0 {
            return Err(ConfigurationError::invalid_frustum_config(
                "inlet_width",
                self.inlet_width,
                "Inlet width must be positive and <= 50.0",
            ));
        }

        // Validate outlet width
        if self.outlet_width <= 0.0 || self.outlet_width > 50.0 {
            return Err(ConfigurationError::invalid_frustum_config(
                "outlet_width",
                self.outlet_width,
                "Outlet width must be positive and <= 50.0",
            ));
        }

        // Validate throat width
        if self.throat_width <= 0.0 || self.throat_width > 25.0 {
            return Err(ConfigurationError::invalid_frustum_config(
                "throat_width",
                self.throat_width,
                "Throat width must be positive and <= 25.0",
            ));
        }

        // Throat must be narrower than both inlet and outlet
        if self.throat_width >= self.inlet_width {
            return Err(ConfigurationError::invalid_frustum_config(
                "throat_width",
                self.throat_width,
                "Throat width must be less than inlet width",
            ));
        }

        if self.throat_width >= self.outlet_width {
            return Err(ConfigurationError::invalid_frustum_config(
                "throat_width",
                self.throat_width,
                "Throat width must be less than outlet width",
            ));
        }

        // Validate smoothness
        if self.smoothness < 10 || self.smoothness > 1000 {
            return Err(ConfigurationError::invalid_frustum_config(
                "smoothness",
                self.smoothness as f64,
                "Smoothness must be between 10 and 1000",
            ));
        }

        // Validate throat position
        if self.throat_position < 0.1 || self.throat_position > 0.9 {
            return Err(ConfigurationError::invalid_frustum_config(
                "throat_position",
                self.throat_position,
                "Throat position must be between 0.1 and 0.9",
            ));
        }

        Ok(())
    }

    /// Calculate width at a given position along the channel (0.0 to 1.0)
    #[must_use]
    pub fn width_at_position(&self, t: f64) -> f64 {
        let t = t.clamp(0.0, 1.0);

        match self.taper_profile {
            TaperProfile::Linear => self.linear_taper(t),
            TaperProfile::Exponential => self.exponential_taper(t),
            TaperProfile::Smooth => self.smooth_taper(t),
        }
    }

    /// Linear taper profile
    fn linear_taper(&self, t: f64) -> f64 {
        if t <= self.throat_position {
            // Inlet to throat
            let local_t = t / self.throat_position;
            (self.throat_width - self.inlet_width).mul_add(local_t, self.inlet_width)
        } else {
            // Throat to outlet
            let local_t = (t - self.throat_position) / (1.0 - self.throat_position);
            (self.outlet_width - self.throat_width).mul_add(local_t, self.throat_width)
        }
    }

    /// Exponential taper profile
    fn exponential_taper(&self, t: f64) -> f64 {
        if t <= self.throat_position {
            // Inlet to throat - exponential decay
            let local_t = t / self.throat_position;
            // Use a normalized exponential that goes from 1 to 0
            let exp_factor = 2.0;
            let factor =
                ((-exp_factor * local_t).exp() - (-exp_factor).exp()) / (1.0 - (-exp_factor).exp());
            (self.throat_width - self.inlet_width).mul_add(1.0 - factor, self.inlet_width)
        } else {
            // Throat to outlet - exponential growth
            let local_t = (t - self.throat_position) / (1.0 - self.throat_position);
            // Use a normalized exponential that goes from 0 to 1
            let exp_factor = 2.0;
            let factor = (1.0 - (-exp_factor * local_t).exp()) / (1.0 - (-exp_factor).exp());
            (self.outlet_width - self.throat_width).mul_add(factor, self.throat_width)
        }
    }

    /// Smooth taper profile using cosine function
    fn smooth_taper(&self, t: f64) -> f64 {
        if t <= self.throat_position {
            // Inlet to throat - smooth cosine transition
            let local_t = t / self.throat_position;
            let factor = 0.5 * (1.0 + (std::f64::consts::PI * local_t).cos());
            (self.inlet_width - self.throat_width).mul_add(factor, self.throat_width)
        } else {
            // Throat to outlet - smooth cosine transition
            let local_t = (t - self.throat_position) / (1.0 - self.throat_position);
            let factor = 0.5 * (1.0 - (std::f64::consts::PI * local_t).cos());
            (self.outlet_width - self.throat_width).mul_add(factor, self.throat_width)
        }
    }
}
