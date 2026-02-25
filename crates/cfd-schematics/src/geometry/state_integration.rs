//! Integration layer between existing strategies and new state management
//!
//! This module provides integration utilities to gradually migrate from
//! hardcoded parameters to the centralized state management system while
//! maintaining backward compatibility.

use crate::{
    config::{ArcConfig, GeometryConfig, OptimizationProfile, SerpentineConfig, WaveShape},
    error::{ConfigurationError, SchemeError, SchemeResult},
    geometry::Point2D,
    state_management::{
        adaptive::ChannelGenerationContext as StateChannelContext, ParameterRegistry,
    },
};
use std::collections::HashMap;

/// Integration helper for serpentine channel parameters
pub struct SerpentineParameterIntegration {
    /// Parameter registry for state management
    registry: ParameterRegistry,

    /// Whether to use adaptive parameters
    use_adaptive: bool,
}

impl SerpentineParameterIntegration {
    /// Create a new integration helper
    pub fn new() -> SchemeResult<Self> {
        let registry = ParameterRegistry::with_defaults().map_err(|e| {
            SchemeError::Configuration(ConfigurationError::MissingConfiguration {
                field: e.to_string(),
            })
        })?;

        Ok(Self {
            registry,
            use_adaptive: true,
        })
    }

    /// Create integration helper with custom registry
    #[must_use]
    pub const fn with_registry(registry: ParameterRegistry) -> Self {
        Self {
            registry,
            use_adaptive: true,
        }
    }

    /// Enable or disable adaptive parameter behavior
    pub const fn set_adaptive(&mut self, adaptive: bool) {
        self.use_adaptive = adaptive;
    }

    /// Apply `SerpentineConfig` to state-managed parameters
    pub fn apply_serpentine_config(&mut self, config: &SerpentineConfig) -> SchemeResult<()> {
        // Get mutable access to serpentine manager
        let serpentine_manager = self.registry.serpentine_mut().map_err(|e| {
            SchemeError::Configuration(ConfigurationError::MissingConfiguration {
                field: e.to_string(),
            })
        })?;

        // Map legacy parameters to new parameter system
        // Note: This is a simplified mapping - in practice, you might want more sophisticated conversion

        // Convert fill_factor to amplitude (rough approximation)
        let estimated_amplitude = config.fill_factor * 8.0; // Reasonable default scaling
        serpentine_manager
            .set_amplitude(estimated_amplitude, "legacy_config_conversion")
            .map_err(|e| {
                SchemeError::Configuration(ConfigurationError::ConflictingValues {
                    conflict: e.to_string(),
                })
            })?;

        // Apply wavelength factor directly
        serpentine_manager
            .set_wavelength_factor(config.wavelength_factor, "legacy_config_conversion")
            .map_err(|e| {
                SchemeError::Configuration(ConfigurationError::ConflictingValues {
                    conflict: e.to_string(),
                })
            })?;

        // Apply other parameters through the parameter system
        // This would be expanded to cover all SerpentineConfig fields

        Ok(())
    }

    /// Get parameters for serpentine generation with optional context adaptation
    pub fn get_serpentine_parameters(
        &self,
        from: Point2D,
        to: Point2D,
        geometry_config: &GeometryConfig,
        box_dims: (f64, f64),
        total_branches: usize,
        neighbor_info: Option<&[f64]>,
    ) -> SchemeResult<SerpentineParameters> {
        let serpentine_manager = self.registry.serpentine();

        if self.use_adaptive {
            // Create state context for adaptive behavior
            let state_context =
                StateChannelContext::new(*geometry_config, box_dims, total_branches, neighbor_info)
                    .with_endpoints(from, to);

            // Get adaptive parameters
            let wave_params = serpentine_manager.get_wave_parameters(Some(&state_context));

            Ok(SerpentineParameters::from_wave_params(wave_params))
        } else {
            // Get base parameters without adaptation
            let wave_params = serpentine_manager.get_wave_parameters(None);

            Ok(SerpentineParameters::from_wave_params(wave_params))
        }
    }

    /// Validate current parameter configuration
    pub fn validate(&self) -> SchemeResult<()> {
        self.registry.validate_all().map_err(|e| {
            SchemeError::Configuration(ConfigurationError::MissingConfiguration {
                field: e.to_string(),
            })
        })
    }

    /// Get parameter registry for advanced usage
    #[must_use]
    pub const fn registry(&self) -> &ParameterRegistry {
        &self.registry
    }

    /// Get mutable parameter registry for advanced usage
    pub const fn registry_mut(&mut self) -> SchemeResult<&mut ParameterRegistry> {
        Ok(&mut self.registry)
    }
}

impl Default for SerpentineParameterIntegration {
    fn default() -> Self {
        Self::new().unwrap_or_else(|_| {
            // Fallback to minimal configuration if registry creation fails
            Self {
                registry: ParameterRegistry::new().unwrap_or_default(),
                use_adaptive: true,
            }
        })
    }
}

/// Integration helper for arc channel parameters
pub struct ArcParameterIntegration {
    /// Parameter registry for state management
    registry: ParameterRegistry,

    /// Whether to use adaptive parameters
    use_adaptive: bool,
}

impl ArcParameterIntegration {
    /// Create a new integration helper
    pub fn new() -> SchemeResult<Self> {
        let registry = ParameterRegistry::with_defaults().map_err(|e| {
            SchemeError::Configuration(ConfigurationError::MissingConfiguration {
                field: e.to_string(),
            })
        })?;

        Ok(Self {
            registry,
            use_adaptive: true,
        })
    }

    /// Create integration helper with custom registry
    #[must_use]
    pub const fn with_registry(registry: ParameterRegistry) -> Self {
        Self {
            registry,
            use_adaptive: true,
        }
    }

    /// Enable or disable adaptive parameter behavior
    pub const fn set_adaptive(&mut self, adaptive: bool) {
        self.use_adaptive = adaptive;
    }

    /// Apply `ArcConfig` to state-managed parameters
    pub fn apply_arc_config(&mut self, _config: &ArcConfig) -> SchemeResult<()> {
        // Get mutable access to arc manager
        let _arc_manager = self.registry.arc_mut().map_err(|e| {
            SchemeError::Configuration(ConfigurationError::MissingConfiguration {
                field: e.to_string(),
            })
        })?;

        // Map parameters to new parameter system
        // For now, we'll store the config directly since the arc manager methods are stubs
        // In a full implementation, these would use proper parameter management

        // This is a simplified approach for the integration layer
        // The actual parameter setting would be implemented when the arc manager is fully developed

        Ok(())
    }

    /// Get parameters for arc generation with optional context adaptation
    pub fn get_arc_parameters(
        &self,
        from: Point2D,
        to: Point2D,
        geometry_config: &GeometryConfig,
        box_dims: (f64, f64),
        total_branches: usize,
        neighbor_info: Option<&[f64]>,
    ) -> SchemeResult<ArcParameters> {
        let arc_manager = self.registry.arc();

        if self.use_adaptive {
            // Create state context for adaptive behavior
            let state_context =
                StateChannelContext::new(*geometry_config, box_dims, total_branches, neighbor_info)
                    .with_endpoints(from, to);

            // Get adaptive parameters (for now, arc manager doesn't have adaptive behavior)
            // but we can still use the context for future enhancements
            Ok(ArcParameters::from_arc_manager(
                arc_manager,
                Some(&state_context),
            ))
        } else {
            // Get base parameters without adaptation
            Ok(ArcParameters::from_arc_manager(arc_manager, None))
        }
    }

    /// Validate current parameter configuration
    pub fn validate(&self) -> SchemeResult<()> {
        self.registry.validate_all().map_err(|e| {
            SchemeError::Configuration(ConfigurationError::MissingConfiguration {
                field: e.to_string(),
            })
        })
    }

    /// Get parameter registry for advanced usage
    #[must_use]
    pub const fn registry(&self) -> &ParameterRegistry {
        &self.registry
    }

    /// Get mutable parameter registry for advanced usage
    pub const fn registry_mut(&mut self) -> SchemeResult<&mut ParameterRegistry> {
        Ok(&mut self.registry)
    }
}

impl Default for ArcParameterIntegration {
    fn default() -> Self {
        Self::new().unwrap_or_else(|_| {
            // Fallback to minimal configuration if registry creation fails
            Self {
                registry: ParameterRegistry::new().unwrap_or_default(),
                use_adaptive: true,
            }
        })
    }
}

/// Structured parameters for serpentine generation
#[derive(Debug, Clone)]
pub struct SerpentineParameters {
    /// Wave amplitude
    pub amplitude: f64,

    /// Wavelength factor
    pub wavelength_factor: f64,

    /// Frequency multiplier
    pub frequency_multiplier: f64,

    /// Phase offset
    pub phase_offset: f64,

    /// Gaussian width factor
    pub gaussian_width_factor: f64,

    /// Wave density factor
    pub wave_density_factor: f64,

    /// Fill factor
    pub fill_factor: f64,

    /// Target fill ratio
    pub target_fill_ratio: f64,
}

impl SerpentineParameters {
    /// Create parameters from wave parameters map
    #[must_use]
    pub fn from_wave_params(wave_params: HashMap<String, f64>) -> Self {
        Self {
            amplitude: wave_params.get("amplitude").copied().unwrap_or(5.0),
            wavelength_factor: wave_params.get("wavelength_factor").copied().unwrap_or(2.0),
            frequency_multiplier: wave_params
                .get("frequency_multiplier")
                .copied()
                .unwrap_or(1.0),
            phase_offset: wave_params.get("phase_offset").copied().unwrap_or(0.0),
            gaussian_width_factor: wave_params
                .get("gaussian_width_factor")
                .copied()
                .unwrap_or(0.3),
            wave_density_factor: wave_params
                .get("wave_density_factor")
                .copied()
                .unwrap_or(2.0),
            fill_factor: wave_params.get("fill_factor").copied().unwrap_or(0.8),
            target_fill_ratio: wave_params.get("target_fill_ratio").copied().unwrap_or(0.9),
        }
    }

    /// Convert to legacy `SerpentineConfig` for backward compatibility
    #[must_use]
    pub fn to_legacy_config(&self) -> SerpentineConfig {
        SerpentineConfig {
            fill_factor: self.fill_factor,
            wavelength_factor: self.wavelength_factor,
            gaussian_width_factor: self.gaussian_width_factor,
            wave_density_factor: self.wave_density_factor,
            wave_phase_direction: 0.0,   // Auto-determine
            wave_shape: WaveShape::Sine, // Default
            optimization_enabled: false, // Default
            target_fill_ratio: self.target_fill_ratio,
            optimization_profile: OptimizationProfile::Balanced, // Default
            adaptive_config: crate::config::AdaptiveSerpentineConfig::default(),
        }
    }

    /// Validate parameter values
    pub fn validate(&self) -> SchemeResult<()> {
        if self.amplitude <= 0.0 {
            return Err(SchemeError::Configuration(
                ConfigurationError::ConflictingValues {
                    conflict: "Amplitude must be positive".to_string(),
                },
            ));
        }

        if self.wavelength_factor <= 0.0 {
            return Err(SchemeError::Configuration(
                ConfigurationError::ConflictingValues {
                    conflict: "Wavelength factor must be positive".to_string(),
                },
            ));
        }

        if self.frequency_multiplier <= 0.0 {
            return Err(SchemeError::Configuration(
                ConfigurationError::ConflictingValues {
                    conflict: "Frequency multiplier must be positive".to_string(),
                },
            ));
        }

        if self.gaussian_width_factor <= 0.0 {
            return Err(SchemeError::Configuration(
                ConfigurationError::ConflictingValues {
                    conflict: "Gaussian width factor must be positive".to_string(),
                },
            ));
        }

        if self.wave_density_factor <= 0.0 {
            return Err(SchemeError::Configuration(
                ConfigurationError::ConflictingValues {
                    conflict: "Wave density factor must be positive".to_string(),
                },
            ));
        }

        if !(0.0..=1.0).contains(&self.fill_factor) {
            return Err(SchemeError::Configuration(
                ConfigurationError::ConflictingValues {
                    conflict: "Fill factor must be between 0.0 and 1.0".to_string(),
                },
            ));
        }

        if !(0.0..=1.0).contains(&self.target_fill_ratio) {
            return Err(SchemeError::Configuration(
                ConfigurationError::ConflictingValues {
                    conflict: "Target fill ratio must be between 0.0 and 1.0".to_string(),
                },
            ));
        }

        Ok(())
    }
}

/// Structured parameters for arc generation
#[derive(Debug, Clone)]
pub struct ArcParameters {
    /// Curvature factor
    pub curvature_factor: f64,

    /// Number of smoothness points
    pub smoothness: usize,

    /// Curvature direction
    pub curvature_direction: f64,

    /// Minimum separation distance
    pub min_separation_distance: f64,

    /// Maximum curvature reduction factor
    pub max_curvature_reduction: f64,

    /// Enable collision prevention
    pub enable_collision_prevention: bool,

    /// Enable adaptive curvature
    pub enable_adaptive_curvature: bool,
}

impl ArcParameters {
    /// Create parameters from arc manager
    #[must_use]
    pub const fn from_arc_manager(
        arc_manager: &crate::state_management::ArcParameterManager,
        _context: Option<&StateChannelContext>,
    ) -> Self {
        Self {
            curvature_factor: arc_manager.get_curvature_factor(),
            smoothness: arc_manager.get_smoothness(),
            curvature_direction: 0.0,     // Default for now
            min_separation_distance: 2.0, // Default for now
            max_curvature_reduction: 0.8, // Default for now
            enable_collision_prevention: arc_manager.is_collision_prevention_enabled(),
            enable_adaptive_curvature: arc_manager.is_adaptive_curvature_enabled(),
        }
    }

    /// Convert to legacy `ArcConfig` for backward compatibility
    #[must_use]
    pub const fn to_legacy_config(&self) -> ArcConfig {
        ArcConfig {
            curvature_factor: self.curvature_factor,
            smoothness: self.smoothness,
            curvature_direction: self.curvature_direction,
            min_separation_distance: self.min_separation_distance,
            enable_collision_prevention: self.enable_collision_prevention,
            max_curvature_reduction: self.max_curvature_reduction,
            enable_adaptive_curvature: self.enable_adaptive_curvature,
        }
    }

    /// Validate parameter values
    pub fn validate(&self) -> SchemeResult<()> {
        if self.curvature_factor < 0.0 {
            return Err(SchemeError::Configuration(
                ConfigurationError::ConflictingValues {
                    conflict: "Curvature factor must be non-negative".to_string(),
                },
            ));
        }

        if self.smoothness < 3 {
            return Err(SchemeError::Configuration(
                ConfigurationError::ConflictingValues {
                    conflict: "Smoothness must be at least 3".to_string(),
                },
            ));
        }

        if self.min_separation_distance <= 0.0 {
            return Err(SchemeError::Configuration(
                ConfigurationError::ConflictingValues {
                    conflict: "Minimum separation distance must be positive".to_string(),
                },
            ));
        }

        if !(0.0..=1.0).contains(&self.max_curvature_reduction) {
            return Err(SchemeError::Configuration(
                ConfigurationError::ConflictingValues {
                    conflict: "Maximum curvature reduction must be between 0.0 and 1.0".to_string(),
                },
            ));
        }

        Ok(())
    }
}

/// Extension trait for existing `SerpentineChannelStrategy` to add state management
pub trait SerpentineStrategyStateExtension {
    /// Generate path using state-managed parameters
    #[allow(clippy::too_many_arguments)]
    fn generate_path_with_state_management(
        &self,
        from: Point2D,
        to: Point2D,
        geometry_config: &GeometryConfig,
        box_dims: (f64, f64),
        total_branches: usize,
        neighbor_info: Option<&[f64]>,
        integration: &SerpentineParameterIntegration,
    ) -> SchemeResult<Vec<Point2D>>;
}

/// Extension trait for existing `ArcChannelStrategy` to add state management
pub trait ArcStrategyStateExtension {
    /// Generate path using state-managed parameters
    #[allow(clippy::too_many_arguments)]
    fn generate_path_with_state_management(
        &self,
        from: Point2D,
        to: Point2D,
        geometry_config: &GeometryConfig,
        box_dims: (f64, f64),
        total_branches: usize,
        neighbor_info: Option<&[f64]>,
        integration: &ArcParameterIntegration,
    ) -> SchemeResult<Vec<Point2D>>;
}

/// Helper function to create a state-managed arc path
pub fn generate_state_managed_arc_path(
    from: Point2D,
    to: Point2D,
    geometry_config: &GeometryConfig,
    box_dims: (f64, f64),
    total_branches: usize,
    neighbor_info: Option<&[f64]>,
    integration: &ArcParameterIntegration,
) -> SchemeResult<Vec<Point2D>> {
    // Get state-managed parameters
    let params = integration.get_arc_parameters(
        from,
        to,
        geometry_config,
        box_dims,
        total_branches,
        neighbor_info,
    )?;

    // Validate parameters
    params.validate()?;

    // Generate path using the parameters
    generate_arc_path_with_params(from, to, &params, geometry_config)
}

/// Generate arc path with specific parameters
#[allow(clippy::unnecessary_wraps)]
fn generate_arc_path_with_params(
    p1: Point2D,
    p2: Point2D,
    params: &ArcParameters,
    geometry_config: &GeometryConfig,
) -> SchemeResult<Vec<Point2D>> {
    let n_points = geometry_config.generation.serpentine_points; // Use serpentine_points for now
    let mut path = Vec::with_capacity(n_points);

    let dx = p2.0 - p1.0;
    let dy = p2.1 - p1.1;
    let channel_length = dx.hypot(dy);

    // Calculate control points for Bezier curve
    let control_point_offset = params.curvature_factor * channel_length * 0.3; // Configurable factor

    // Calculate perpendicular direction for control points
    let angle = dy.atan2(dx);
    let perp_x = -angle.sin();
    let perp_y = angle.cos();

    // Determine curvature direction (simplified logic for now)
    let direction_factor = if params.curvature_direction == 0.0 {
        // Auto-determine based on position (could be enhanced with context)
        1.0
    } else {
        params.curvature_direction.signum()
    };

    // Calculate control points
    let mid_x = f64::midpoint(p1.0, p2.0);
    let mid_y = f64::midpoint(p1.1, p2.1);
    let control_x = (perp_x * control_point_offset).mul_add(direction_factor, mid_x);
    let control_y = (perp_y * control_point_offset).mul_add(direction_factor, mid_y);

    // Generate quadratic Bezier curve points
    for i in 0..n_points {
        let t = i as f64 / (n_points - 1) as f64;

        // Quadratic Bezier formula: B(t) = (1-t)²P₀ + 2(1-t)tP₁ + t²P₂
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

    // Apply collision prevention if enabled
    if params.enable_collision_prevention {
        apply_collision_prevention(&mut path, params);
    }

    Ok(path)
}

/// Apply collision prevention to arc path
fn apply_collision_prevention(path: &mut [Point2D], params: &ArcParameters) {
    // This is a simplified collision prevention implementation
    // In practice, this would use neighbor information and wall constraints

    if !params.enable_collision_prevention {
        return;
    }

    // Apply curvature reduction if needed (simplified logic)
    let reduction_factor = 1.0 - params.max_curvature_reduction;

    if reduction_factor < 1.0 {
        // Reduce curvature by moving points closer to the straight line
        let start = path[0];
        let end = path[path.len() - 1];

        let path_len = path.len();
        for (i, point) in path.iter_mut().enumerate() {
            let t = i as f64 / (path_len - 1) as f64;
            let straight_x = t.mul_add(end.0 - start.0, start.0);
            let straight_y = t.mul_add(end.1 - start.1, start.1);

            // Interpolate between curved and straight path
            point.0 = point
                .0
                .mul_add(reduction_factor, straight_x * (1.0 - reduction_factor));
            point.1 = point
                .1
                .mul_add(reduction_factor, straight_y * (1.0 - reduction_factor));
        }
    }
}

/// Helper function to create a state-managed serpentine path
pub fn generate_state_managed_serpentine_path(
    from: Point2D,
    to: Point2D,
    geometry_config: &GeometryConfig,
    box_dims: (f64, f64),
    total_branches: usize,
    neighbor_info: Option<&[f64]>,
    integration: &SerpentineParameterIntegration,
) -> SchemeResult<Vec<Point2D>> {
    // Get state-managed parameters
    let params = integration.get_serpentine_parameters(
        from,
        to,
        geometry_config,
        box_dims,
        total_branches,
        neighbor_info,
    )?;

    // Validate parameters
    params.validate()?;

    // Generate path using the parameters
    // This is a simplified implementation - the full version would use all the
    // sophisticated logic from the original SerpentineChannelStrategy
    generate_serpentine_path_with_params(from, to, &params, geometry_config)
}

/// Generate serpentine path with specific parameters
#[allow(clippy::unnecessary_wraps)]
fn generate_serpentine_path_with_params(
    p1: Point2D,
    p2: Point2D,
    params: &SerpentineParameters,
    geometry_config: &GeometryConfig,
) -> SchemeResult<Vec<Point2D>> {
    let n_points = geometry_config.generation.serpentine_points;
    let mut path = Vec::with_capacity(n_points);

    let dx = p2.0 - p1.0;
    let dy = p2.1 - p1.1;
    let channel_length = dx.hypot(dy);

    // Calculate wavelength and periods
    let base_wavelength = params.wavelength_factor * geometry_config.channel_width;
    let length_based_periods = (channel_length / base_wavelength) * params.wave_density_factor;
    let base_periods = length_based_periods.max(1.0);
    let half_periods = (base_periods * 2.0).round().max(1.0);

    // Generate path points
    for i in 0..n_points {
        let t = i as f64 / (n_points - 1) as f64;

        // Linear interpolation for base position
        let x = t.mul_add(dx, p1.0);
        let y = t.mul_add(dy, p1.1);

        // Calculate envelopes
        let smooth_envelope = calculate_smooth_envelope(t);
        let gaussian_envelope = calculate_gaussian_envelope(t, params.gaussian_width_factor);
        let envelope = smooth_envelope * gaussian_envelope;

        // Calculate wave phase
        let wave_phase = std::f64::consts::PI * half_periods * t * params.frequency_multiplier;

        // Calculate wave amplitude (sine wave for now)
        let wave_amplitude = (wave_phase + params.phase_offset).sin();

        // Calculate perpendicular offset
        let perpendicular_amplitude = params.amplitude * envelope * wave_amplitude;
        let angle = dy.atan2(dx);
        let perp_x = -angle.sin() * perpendicular_amplitude;
        let perp_y = angle.cos() * perpendicular_amplitude;

        path.push((x + perp_x, y + perp_y));
    }

    Ok(path)
}

/// Calculate smooth envelope for endpoints
fn calculate_smooth_envelope(t: f64) -> f64 {
    let constants = crate::config::ConstantsRegistry::new();
    let transition_zone = constants.get_transition_zone_factor();
    if t < transition_zone {
        0.5 * (1.0 - (std::f64::consts::PI * t / transition_zone).cos())
    } else if t > 1.0 - transition_zone {
        0.5 * (1.0 - (std::f64::consts::PI * (1.0 - t) / transition_zone).cos())
    } else {
        1.0
    }
}

/// Calculate Gaussian envelope
fn calculate_gaussian_envelope(t: f64, gaussian_width_factor: f64) -> f64 {
    let sigma = 1.0 / gaussian_width_factor;
    let center = 0.5;
    let exponent = -0.5 * ((t - center) / sigma).powi(2);
    exponent.exp()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_serpentine_parameter_integration() {
        let integration = SerpentineParameterIntegration::new()
            .expect("Failed to create SerpentineParameterIntegration for test");

        // Test parameter validation
        assert!(integration.validate().is_ok());
    }

    #[test]
    fn test_legacy_config_conversion() {
        let mut integration = SerpentineParameterIntegration::new()
            .expect("Failed to create SerpentineParameterIntegration for test");
        let legacy_config = SerpentineConfig::default();

        // Apply legacy config
        assert!(integration.apply_serpentine_config(&legacy_config).is_ok());

        // Validate after conversion
        assert!(integration.validate().is_ok());
    }

    #[test]
    fn test_parameter_retrieval() {
        let integration = SerpentineParameterIntegration::new()
            .expect("Failed to create SerpentineParameterIntegration for test");
        let geometry_config = GeometryConfig::default();

        let params = integration
            .get_serpentine_parameters(
                (0.0, 0.0),
                (100.0, 0.0),
                &geometry_config,
                (200.0, 100.0),
                4,
                None,
            )
            .expect("Failed to get serpentine parameters for test");

        // Validate retrieved parameters
        assert!(params.validate().is_ok());
        assert!(params.amplitude > 0.0);
        assert!(params.wavelength_factor > 0.0);
    }

    #[test]
    fn test_arc_parameter_integration() {
        let integration = ArcParameterIntegration::new()
            .expect("Failed to create ArcParameterIntegration for test");

        // Test parameter validation
        assert!(integration.validate().is_ok());
    }

    #[test]
    fn test_arc_legacy_config_conversion() {
        let mut integration = ArcParameterIntegration::new()
            .expect("Failed to create ArcParameterIntegration for test");
        let legacy_config = ArcConfig::default();

        // Apply legacy config
        assert!(integration.apply_arc_config(&legacy_config).is_ok());

        // Validate after conversion
        assert!(integration.validate().is_ok());
    }

    #[test]
    fn test_arc_parameter_retrieval() {
        let integration = ArcParameterIntegration::new()
            .expect("Failed to create ArcParameterIntegration for test");
        let geometry_config = GeometryConfig::default();

        let params = integration
            .get_arc_parameters(
                (0.0, 0.0),
                (100.0, 0.0),
                &geometry_config,
                (200.0, 100.0),
                4,
                None,
            )
            .expect("Failed to get arc parameters for test");

        // Validate retrieved parameters
        assert!(params.validate().is_ok());
        assert!(params.curvature_factor >= 0.0);
        assert!(params.smoothness >= 3);
    }

    #[test]
    fn test_state_managed_arc_path_generation() {
        let integration = ArcParameterIntegration::new().unwrap();
        let geometry_config = GeometryConfig::default();

        let path = generate_state_managed_arc_path(
            (0.0, 50.0),
            (100.0, 50.0),
            &geometry_config,
            (200.0, 100.0),
            4,
            None,
            &integration,
        )
        .unwrap();

        // Validate generated path
        assert!(!path.is_empty());
        assert_eq!(path.len(), geometry_config.generation.serpentine_points); // Using serpentine_points for now

        // Check that path starts and ends at correct points
        let first_point = path.first().unwrap();
        let last_point = path.last().unwrap();

        assert!((first_point.0 - 0.0).abs() < 1e-6);
        assert!((last_point.0 - 100.0).abs() < 1e-6);
    }

    #[test]
    fn test_state_managed_path_generation() {
        let integration = SerpentineParameterIntegration::new().unwrap();
        let geometry_config = GeometryConfig::default();

        let path = generate_state_managed_serpentine_path(
            (0.0, 50.0),
            (100.0, 50.0),
            &geometry_config,
            (200.0, 100.0),
            4,
            None,
            &integration,
        )
        .unwrap();

        // Validate generated path
        assert!(!path.is_empty());
        assert_eq!(path.len(), geometry_config.generation.serpentine_points);

        // Check that path starts and ends at correct points
        let first_point = path.first().unwrap();
        let last_point = path.last().unwrap();

        assert!((first_point.0 - 0.0).abs() < 1e-6);
        assert!((last_point.0 - 100.0).abs() < 1e-6);
    }
}
