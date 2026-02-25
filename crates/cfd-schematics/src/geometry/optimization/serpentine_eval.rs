//! Simplified serpentine path generation for optimization evaluation
//!
//! Provides a fast serpentine path generator used by the optimization algorithms
//! to evaluate candidate parameter sets. This is a streamlined version focused on
//! speed rather than the full feature set of `SerpentineChannelStrategy`.

use crate::config::{GeometryConfig, SerpentineConfig};
use crate::geometry::types::Point2D;

use super::constants;

/// Optimized serpentine path generation with aggressive amplitude calculation
#[must_use]
pub(super) fn generate_simplified_serpentine_path(
    p1: Point2D,
    p2: Point2D,
    geometry_config: &GeometryConfig,
    serpentine_config: &SerpentineConfig,
    box_dims: (f64, f64),
    neighbor_info: Option<&[f64]>,
) -> Vec<Point2D> {
    // Calculate amplitude first to check threshold
    let amplitude = calculate_optimized_amplitude(
        p1,
        p2,
        geometry_config,
        serpentine_config,
        box_dims,
        neighbor_info,
    );

    // If amplitude is too small, return straight line
    if amplitude <= 0.0 {
        let n_points = geometry_config.generation.optimization_points;
        let dx = p2.0 - p1.0;
        let dy = p2.1 - p1.1;

        return (0..n_points)
            .map(|i| {
                let t = i as f64 / (n_points - 1) as f64;
                (t.mul_add(dx, p1.0), t.mul_add(dy, p1.1))
            })
            .collect();
    }

    let n_points = geometry_config.generation.optimization_points;
    let mut path = Vec::with_capacity(n_points);

    let dx = p2.0 - p1.0;
    let dy = p2.1 - p1.1;
    let channel_length = dx.hypot(dy);

    // Calculate amplitude based on available space
    let channel_center_y = f64::midpoint(p1.1, p2.1);
    let box_height = box_dims.1;

    // Calculate available space considering neighbors
    let mut available_space_above = box_height - channel_center_y;
    let mut available_space_below = channel_center_y;

    if let Some(neighbors) = neighbor_info {
        for &neighbor_y in neighbors {
            if neighbor_y > channel_center_y {
                available_space_above = available_space_above.min(neighbor_y - channel_center_y);
            } else {
                available_space_below = available_space_below.min(channel_center_y - neighbor_y);
            }
        }
    }

    let available_space = available_space_above.min(available_space_below);

    // Calculate and validate wavelength for manufacturing constraints
    let initial_wavelength = serpentine_config.wavelength_factor * geometry_config.channel_width;
    let min_separation = geometry_config.channel_width;
    let min_wavelength = min_separation * 3.0; // Conservative constraint for proper spacing
    let base_wavelength = initial_wavelength.max(min_wavelength);

    // Calculate space-constrained amplitude
    let space_constrained_amplitude = available_space.max(0.0);

    // Calculate wavelength-constrained amplitude for manufacturing
    let min_separation_distance = geometry_config.channel_width;
    let wavelength_constrained_amplitude = if base_wavelength > min_separation_distance * 2.0 {
        // If wavelength is reasonable, allow full amplitude utilization
        space_constrained_amplitude // Use full available space
    } else if base_wavelength > min_separation_distance * 1.2 {
        // If wavelength is moderate, allow substantial amplitude
        space_constrained_amplitude * 0.85 // Use 85% of available space
    } else {
        // If wavelength is small, be more conservative but still generous
        space_constrained_amplitude * 0.6 // Use 60% of available space
    };

    // Apply neighbor constraints if present
    let neighbor_constrained_amplitude = if let Some(neighbor_info) = neighbor_info {
        let channel_center_y = f64::midpoint(p1.1, p2.1);
        let mut min_neighbor_amplitude = f64::INFINITY;

        for &neighbor_y in neighbor_info {
            let distance_to_neighbor = (channel_center_y - neighbor_y).abs();
            let safe_amplitude =
                ((distance_to_neighbor - geometry_config.channel_width) * 0.5).max(0.0);
            min_neighbor_amplitude = min_neighbor_amplitude.min(safe_amplitude);
        }

        min_neighbor_amplitude
    } else {
        f64::INFINITY
    };

    // Take the most restrictive constraint and apply fill factor
    let _max_safe_amplitude = space_constrained_amplitude
        .min(wavelength_constrained_amplitude)
        .min(neighbor_constrained_amplitude)
        .max(0.0);

    // Use the amplitude already calculated above

    // Use optimized wavelength calculation
    let base_wavelength = calculate_optimized_wavelength(geometry_config, serpentine_config);

    // Generate simplified serpentine path with smooth endpoints
    let base_periods = (channel_length / base_wavelength) * serpentine_config.wave_density_factor;
    // Round to nearest integer number of half-periods to ensure sin(pi*n) = 0 at endpoints
    let half_periods = (base_periods * 2.0).round().max(1.0);

    for i in 0..n_points {
        let t = i as f64 / (n_points - 1) as f64;

        let base_x = t.mul_add(dx, p1.0);
        let base_y = t.mul_add(dy, p1.1);

        // Apply smooth endpoint envelope combined with improved Gaussian envelope
        let smooth_envelope = calculate_smooth_endpoint_envelope_for_optimization(t);
        let gaussian_envelope = calculate_improved_envelope_for_optimization(
            t,
            channel_length,
            dx,
            dy,
            serpentine_config,
        );
        let envelope = smooth_envelope * gaussian_envelope;

        let wave_phase = std::f64::consts::PI * half_periods * t;
        // Calculate wave amplitude based on wave shape
        let wave_value = match serpentine_config.wave_shape {
            crate::config::WaveShape::Sine => wave_phase.sin(),
            crate::config::WaveShape::Square => {
                let sine_value = wave_phase.sin();
                let sharpness = constants::SQUARE_WAVE_SHARPNESS;
                (sharpness * sine_value).tanh()
            }
        };
        let wave_amplitude = amplitude * envelope * wave_value;

        let perp_x = -dy / channel_length;
        let perp_y = dx / channel_length;

        let x = wave_amplitude.mul_add(perp_x, base_x);
        let y = wave_amplitude.mul_add(perp_y, base_y);

        // Ensure exact endpoint matching for first and last points to maintain precision
        // The smooth envelope should make wave_amplitude ~= 0 at endpoints, but we ensure exactness
        if i == 0 {
            path.push(p1);
        } else if i == n_points - 1 {
            path.push(p2);
        } else {
            path.push((x, y));
        }
    }

    path
}

/// Calculate optimized amplitude with aggressive space utilization
fn calculate_optimized_amplitude(
    p1: Point2D,
    p2: Point2D,
    geometry_config: &GeometryConfig,
    serpentine_config: &SerpentineConfig,
    box_dims: (f64, f64),
    neighbor_info: Option<&[f64]>,
) -> f64 {
    let channel_center_y = f64::midpoint(p1.1, p2.1);

    // Dynamic space analysis with aggressive space utilization
    let available_space =
        calculate_available_space(p1, p2, geometry_config, box_dims, neighbor_info);

    // Use the full available space - safety margins already included in calculation
    let base_amplitude = available_space;

    // Position-based enhancement (more conservative)
    let center_distance = (channel_center_y - box_dims.1 / 2.0).abs();
    let position_factor = (1.0 - center_distance / (box_dims.1 / 2.0)).mul_add(0.08, 1.0); // Up to 8% boost

    // Density-based enhancement (more conservative)
    let box_area = box_dims.0 * box_dims.1;
    let density_factor = if box_area > 10000.0 { 1.05 } else { 1.02 }; // Modest boost for larger areas

    // Apply enhancements and fill factor
    let enhanced_amplitude = base_amplitude * position_factor * density_factor;
    let final_amplitude = enhanced_amplitude * serpentine_config.fill_factor;

    // Ensure minimum visibility
    let min_amplitude = geometry_config.channel_width * 0.08;
    final_amplitude.max(min_amplitude)
}

/// Calculate available space considering neighbors and walls with dynamic space utilization
fn calculate_available_space(
    p1: Point2D,
    p2: Point2D,
    geometry_config: &GeometryConfig,
    box_dims: (f64, f64),
    neighbor_info: Option<&[f64]>,
) -> f64 {
    let channel_center_y = f64::midpoint(p1.1, p2.1);
    let box_height = box_dims.1;
    let channel_width = geometry_config.channel_width;
    let wall_margin = geometry_config.wall_clearance + channel_width * 0.5;

    if let Some(neighbors) = neighbor_info {
        // Find the closest neighbors above and below
        let mut closest_above = box_height;
        let mut closest_below: f64 = 0.0;

        for &neighbor_y in neighbors {
            if (neighbor_y - channel_center_y).abs() > f64::EPSILON {
                // Exclude self
                if neighbor_y > channel_center_y {
                    closest_above = closest_above.min(neighbor_y);
                } else {
                    closest_below = closest_below.max(neighbor_y);
                }
            }
        }

        // Calculate available space above and below
        let space_above = if closest_above < box_height {
            ((closest_above - channel_center_y) - channel_width) * 0.5
        } else {
            // No neighbor above - use distance to wall
            box_height - channel_center_y - wall_margin
        };

        let space_below = if closest_below > 0.0 {
            ((channel_center_y - closest_below) - channel_width) * 0.5
        } else {
            // No neighbor below - use distance to wall
            channel_center_y - wall_margin
        };

        // Use the minimum of available spaces, but ensure it's positive
        space_above.min(space_below).max(0.0)
    } else {
        // Single channel - use full box space minus wall clearances
        let space_above = box_height - channel_center_y - wall_margin;
        let space_below = channel_center_y - wall_margin;
        space_above.min(space_below).max(0.0)
    }
}

/// Calculate optimized wavelength with manufacturing constraints
fn calculate_optimized_wavelength(
    geometry_config: &GeometryConfig,
    serpentine_config: &SerpentineConfig,
) -> f64 {
    let initial_wavelength = serpentine_config.wavelength_factor * geometry_config.channel_width;
    let min_separation = geometry_config.channel_width;
    let min_wavelength = min_separation * 3.0; // Conservative minimum for proper spacing

    initial_wavelength.max(min_wavelength)
}

/// Calculate smooth endpoint envelope for optimization
///
/// Uses smoothstep function for C1 continuity at endpoints
#[must_use]
fn calculate_smooth_endpoint_envelope_for_optimization(t: f64) -> f64 {
    // Smoothstep function: t^2(3-2t)
    t * t * constants::SMOOTHSTEP_COEFFICIENT_2.mul_add(-t, constants::SMOOTHSTEP_COEFFICIENT_1)
}

/// Calculate improved Gaussian envelope for optimization (helper function)
///
/// This mirrors the logic from `SerpentineChannelStrategy` but is available
/// for use in the optimization module.
#[must_use]
fn calculate_improved_envelope_for_optimization(
    t: f64,
    channel_length: f64,
    dx: f64,
    dy: f64,
    serpentine_config: &SerpentineConfig,
) -> f64 {
    // Calculate the actual distance between nodes
    let node_distance = dx.hypot(dy);

    // Determine if this is primarily a horizontal channel (middle section logic)
    let is_horizontal = dx.abs() > dy.abs();
    let horizontal_ratio = dx.abs() / node_distance;

    // For horizontal channels (middle sections), we want less aggressive tapering
    let middle_section_factor = if is_horizontal
        && horizontal_ratio > serpentine_config.adaptive_config.horizontal_ratio_threshold
    {
        (1.0 - serpentine_config
            .adaptive_config
            .middle_section_amplitude_factor)
            .mul_add(
                horizontal_ratio,
                serpentine_config
                    .adaptive_config
                    .middle_section_amplitude_factor,
            )
    } else {
        1.0
    };

    // Distance-based normalization
    let distance_normalization = if serpentine_config
        .adaptive_config
        .enable_distance_based_scaling
    {
        (node_distance
            / serpentine_config
                .adaptive_config
                .node_distance_normalization)
            .clamp(constants::MIN_DISTANCE_NORMALIZATION, constants::MAX_DISTANCE_NORMALIZATION)
    } else {
        constants::MAX_DISTANCE_NORMALIZATION // No distance-based scaling when disabled
    };

    // Calculate effective sigma
    let base_sigma = channel_length / serpentine_config.gaussian_width_factor;
    let effective_sigma = base_sigma * distance_normalization * middle_section_factor;

    // Center the envelope
    let center = constants::GAUSSIAN_CENTER;

    // Create smooth dome-shaped envelope instead of sharp Gaussian peaks
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
            let smoothstep =
                (edge_distance * edge_distance).mul_add(-2.0f64.mul_add(-edge_distance, 3.0), 1.0);
            smoothstep * 0.05 // Very small amplitude at edges
        } else {
            0.0
        }
    };

    // For horizontal middle sections, enhance the dome but keep it smooth
    if is_horizontal
        && horizontal_ratio > serpentine_config.adaptive_config.horizontal_ratio_threshold
    {
        let plateau_width = serpentine_config
            .adaptive_config
            .plateau_width_factor
            .min(0.3); // Limit plateau width
        let plateau_start = 0.5 - plateau_width / 2.0;
        let plateau_end = 0.5 + plateau_width / 2.0;

        if t >= plateau_start && t <= plateau_end {
            // In the plateau region, enhance the dome but keep it smooth
            let plateau_factor = 1.0 - ((t - 0.5).abs() / (plateau_width / 2.0));
            let enhanced_amplitude =
                (1.0 - serpentine_config.adaptive_config.plateau_amplitude_factor).mul_add(
                    plateau_factor,
                    serpentine_config.adaptive_config.plateau_amplitude_factor,
                );
            dome_envelope.max(enhanced_amplitude * dome_envelope)
        } else {
            dome_envelope
        }
    } else {
        dome_envelope
    }
}
