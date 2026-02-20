use crate::config::channel::arc::ArcConfig;
use crate::config::channel::serpentine::{OptimizationProfile, SerpentineConfig, WaveShape};
use crate::config::constants::primitives as constants;
use crate::config::adaptive::AdaptiveSerpentineConfig;
use crate::config::geometry::{GeometryConfig, GeometryGenerationConfig};

/// Preset for microfluidic devices with fine features
#[must_use]
pub const fn fine_features() -> GeometryConfig {
    GeometryConfig {
        wall_clearance: 0.2,
        channel_width: 0.5,
        channel_height: 0.3,
        generation: GeometryGenerationConfig::high_quality(),
    }
}

/// Preset for standard microfluidic devices
#[must_use]
pub fn standard() -> GeometryConfig {
    GeometryConfig::default()
}

/// Preset for large-scale microfluidic devices
#[must_use]
pub const fn large_scale() -> GeometryConfig {
    GeometryConfig {
        wall_clearance: 2.0,
        channel_width: 5.0,
        channel_height: 2.0,
        generation: GeometryGenerationConfig::fast(),
    }
}

/// Preset for high-density serpentine channels
#[must_use]
pub fn high_density_serpentine() -> SerpentineConfig {
    SerpentineConfig {
        fill_factor: 0.9,
        wavelength_factor: 2.0,
        gaussian_width_factor: 8.0,
        wave_density_factor: 4.0,
        wave_phase_direction: 0.0,        // Auto-symmetric
        wave_shape: WaveShape::Sine,      // Default to sine wave
        optimization_enabled: false,
        target_fill_ratio: 0.9,
        optimization_profile: OptimizationProfile::Balanced,
        adaptive_config: AdaptiveSerpentineConfig::aggressive(), // High-density needs aggressive adaptation
    }
}

/// Preset for smooth serpentine channels
#[must_use]
pub fn smooth_serpentine() -> SerpentineConfig {
    SerpentineConfig {
        fill_factor: 0.6,
        wavelength_factor: 4.0,
        gaussian_width_factor: 10.0,
        wave_density_factor: 1.0, // Lower than default (1.5) for smoother, less dense waves
        wave_phase_direction: 0.0, // Auto-symmetric
        wave_shape: WaveShape::Sine, // Default to sine wave
        optimization_enabled: false,
        target_fill_ratio: 0.9,
        optimization_profile: OptimizationProfile::Balanced,
        adaptive_config: AdaptiveSerpentineConfig::conservative(), // Smooth channels need conservative adaptation
    }
}

/// Preset for inward-phase serpentine channels (all waves start inward)
#[must_use]
pub fn inward_serpentines() -> SerpentineConfig {
    SerpentineConfig {
        fill_factor: 0.8,
        wavelength_factor: 3.0,
        gaussian_width_factor: 6.0,
        wave_density_factor: 2.0,
        wave_phase_direction: -1.0,       // Force inward phase
        wave_shape: WaveShape::Sine,      // Default to sine wave
        optimization_enabled: false,
        target_fill_ratio: 0.9,
        optimization_profile: OptimizationProfile::Balanced,
        adaptive_config: AdaptiveSerpentineConfig::default(), // Default adaptive behavior
    }
}

/// Preset for outward-phase serpentine channels (all waves start outward)
#[must_use]
pub fn outward_serpentines() -> SerpentineConfig {
    SerpentineConfig {
        fill_factor: 0.8,
        wavelength_factor: 3.0,
        gaussian_width_factor: 6.0,
        wave_density_factor: 2.0,
        wave_phase_direction: 1.0,        // Force outward phase
        wave_shape: WaveShape::Sine,      // Default to sine wave
        optimization_enabled: false,
        target_fill_ratio: 0.9,
        optimization_profile: OptimizationProfile::Balanced,
        adaptive_config: AdaptiveSerpentineConfig::default(), // Default adaptive behavior
    }
}

/// Preset for length-optimized serpentine channels
#[must_use]
pub fn optimized_serpentine() -> SerpentineConfig {
    SerpentineConfig {
        fill_factor: 0.8,
        wavelength_factor: 3.0,
        gaussian_width_factor: 6.0,
        wave_density_factor: 2.0,
        wave_phase_direction: 0.0,        // Auto-symmetric
        wave_shape: WaveShape::Sine,      // Default to sine wave
        optimization_enabled: true,
        target_fill_ratio: 0.95, // Aggressive optimization target
        optimization_profile: OptimizationProfile::Balanced,
        adaptive_config: AdaptiveSerpentineConfig::default(), // Default adaptive behavior
    }
}

/// Preset for fast-optimized serpentine channels
#[must_use]
pub fn fast_optimized_serpentine() -> SerpentineConfig {
    SerpentineConfig {
        fill_factor: 0.8,
        wavelength_factor: 3.0,
        gaussian_width_factor: 6.0,
        wave_density_factor: 2.0,
        wave_phase_direction: 0.0,        // Auto-symmetric
        wave_shape: WaveShape::Sine,      // Default to sine wave
        optimization_enabled: true,
        target_fill_ratio: 0.9, // Moderate optimization target
        optimization_profile: OptimizationProfile::Fast,
        adaptive_config: AdaptiveSerpentineConfig::default(), // Default adaptive behavior
    }
}

/// Preset for thorough-optimized serpentine channels
#[must_use]
pub fn thorough_optimized_serpentine() -> SerpentineConfig {
    SerpentineConfig {
        fill_factor: 0.8,
        wavelength_factor: 3.0,
        gaussian_width_factor: 6.0,
        wave_density_factor: 2.0,
        wave_phase_direction: 0.0,        // Auto-symmetric
        wave_shape: WaveShape::Sine,      // Default to sine wave
        optimization_enabled: true,
        target_fill_ratio: 0.98, // Very aggressive optimization target
        optimization_profile: OptimizationProfile::Thorough,
        adaptive_config: AdaptiveSerpentineConfig::default(), // Default adaptive behavior
    }
}

/// Preset for square wave serpentine channels (angular transitions)
#[must_use]
pub fn square_wave_serpentine() -> SerpentineConfig {
    SerpentineConfig {
        fill_factor: 0.8,
        wavelength_factor: 3.0,
        gaussian_width_factor: 6.0,
        wave_density_factor: 2.0,
        wave_phase_direction: 0.0,     // Auto-symmetric
        wave_shape: WaveShape::Square, // Sharp square wave transitions
        optimization_enabled: false,
        target_fill_ratio: 0.9,
        optimization_profile: OptimizationProfile::Balanced,
        adaptive_config: AdaptiveSerpentineConfig::default(),
    }
}

/// Preset for subtle arc channels
#[must_use]
pub const fn subtle_arcs() -> ArcConfig {
    ArcConfig {
        curvature_factor: 0.2,
        smoothness: 30,
        curvature_direction: 0.0, // Auto-determine
        min_separation_distance: constants::DEFAULT_MIN_SEPARATION_DISTANCE,
        enable_collision_prevention: true,
        max_curvature_reduction: constants::DEFAULT_MAX_CURVATURE_REDUCTION,
        enable_adaptive_curvature: true,
    }
}

/// Preset for high-quality geometry generation
#[must_use]
pub const fn high_quality_generation() -> GeometryConfig {
    GeometryConfig {
        wall_clearance: constants::DEFAULT_WALL_CLEARANCE,
        channel_width: constants::DEFAULT_CHANNEL_WIDTH,
        channel_height: constants::DEFAULT_CHANNEL_HEIGHT,
        generation: GeometryGenerationConfig::high_quality(),
    }
}

/// Preset for fast geometry generation
#[must_use]
pub const fn fast_generation() -> GeometryConfig {
    GeometryConfig {
        wall_clearance: constants::DEFAULT_WALL_CLEARANCE,
        channel_width: constants::DEFAULT_CHANNEL_WIDTH,
        channel_height: constants::DEFAULT_CHANNEL_HEIGHT,
        generation: GeometryGenerationConfig::fast(),
    }
}

/// Preset for microfluidic research applications
#[must_use]
pub const fn research_grade() -> GeometryConfig {
    GeometryConfig {
        wall_clearance: 0.3,
        channel_width: 0.8,
        channel_height: 0.5,
        generation: GeometryGenerationConfig::high_quality(),
    }
}

/// Preset for industrial manufacturing applications
#[must_use]
pub const fn manufacturing_grade() -> GeometryConfig {
    GeometryConfig {
        wall_clearance: 1.0,
        channel_width: 2.0,
        channel_height: 1.5,
        generation: GeometryGenerationConfig::fast(),
    }
}

/// Preset for pronounced arc channels
#[must_use]
pub const fn pronounced_arcs() -> ArcConfig {
    ArcConfig {
        curvature_factor: 0.8,
        smoothness: 50,
        curvature_direction: 0.0, // Auto-determine
        min_separation_distance: constants::DEFAULT_MIN_SEPARATION_DISTANCE,
        enable_collision_prevention: true,
        max_curvature_reduction: constants::DEFAULT_MAX_CURVATURE_REDUCTION,
        enable_adaptive_curvature: true,
    }
}

/// Preset for inward-curving arc channels (all arcs curve toward center)
#[must_use]
pub const fn inward_arcs() -> ArcConfig {
    ArcConfig {
        curvature_factor: 0.5,
        smoothness: 30,
        curvature_direction: -1.0, // Force inward curvature
        min_separation_distance: constants::DEFAULT_MIN_SEPARATION_DISTANCE,
        enable_collision_prevention: true,
        max_curvature_reduction: constants::DEFAULT_MAX_CURVATURE_REDUCTION,
        enable_adaptive_curvature: true,
    }
}

/// Preset for outward-curving arc channels (all arcs curve away from center)
#[must_use]
pub const fn outward_arcs() -> ArcConfig {
    ArcConfig {
        curvature_factor: 0.5,
        smoothness: 30,
        curvature_direction: 1.0, // Force outward curvature
        min_separation_distance: constants::DEFAULT_MIN_SEPARATION_DISTANCE,
        enable_collision_prevention: true,
        max_curvature_reduction: constants::DEFAULT_MAX_CURVATURE_REDUCTION,
        enable_adaptive_curvature: true,
    }
}

/// Preset for safe high-curvature arcs with collision prevention
#[must_use]
pub const fn safe_high_curvature_arcs() -> ArcConfig {
    ArcConfig {
        curvature_factor: 1.5, // High curvature but with safety measures
        smoothness: 50,
        curvature_direction: 0.0,     // Auto-determine
        min_separation_distance: 2.0, // Increased separation for safety
        enable_collision_prevention: true,
        max_curvature_reduction: 0.3, // Allow significant reduction if needed
        enable_adaptive_curvature: true,
    }
}

/// Preset for maximum curvature with aggressive collision prevention
#[must_use]
pub const fn maximum_safe_arcs() -> ArcConfig {
    ArcConfig {
        curvature_factor: 2.0, // Maximum curvature
        smoothness: 100,
        curvature_direction: 0.0,     // Auto-determine
        min_separation_distance: 3.0, // Maximum separation for safety
        enable_collision_prevention: true,
        max_curvature_reduction: 0.1, // Allow maximum reduction if needed
        enable_adaptive_curvature: true,
    }
}

/// Preset for dense layouts with conservative curvature
#[must_use]
pub const fn dense_layout_arcs() -> ArcConfig {
    ArcConfig {
        curvature_factor: 0.3, // Conservative curvature
        smoothness: 25,
        curvature_direction: 0.0,     // Auto-determine
        min_separation_distance: 0.5, // Tight separation for dense layouts
        enable_collision_prevention: true,
        max_curvature_reduction: 0.7, // Allow significant reduction
        enable_adaptive_curvature: true,
    }
}
