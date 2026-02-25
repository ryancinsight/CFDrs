//! Optimization algorithms for serpentine channel parameters
//!
//! Provides grid-search (fast), Nelder-Mead simplex (balanced), and
//! multi-start Nelder-Mead (thorough) optimization strategies.

use crate::config::{ConstantsRegistry, GeometryConfig, OptimizationProfile, SerpentineConfig};
use crate::geometry::types::Point2D;

use super::constants;
use super::serpentine_eval::generate_simplified_serpentine_path;
use super::{
    calculate_constraint_penalty, calculate_min_neighbor_distance, calculate_min_wall_distance,
    calculate_path_length, OptimizationParams, OptimizationResult,
};

/// Optimize serpentine parameters to maximize channel length using advanced algorithms
///
/// # Arguments
/// * `p1` - Start point of the channel
/// * `p2` - End point of the channel
/// * `geometry_config` - Geometry configuration
/// * `serpentine_config` - Serpentine configuration
/// * `box_dims` - Box dimensions
/// * `neighbor_info` - Optional neighbor channel positions
///
/// # Returns
/// Optimized parameters that maximize channel length while maintaining constraints
#[must_use]
pub fn optimize_serpentine_parameters(
    p1: Point2D,
    p2: Point2D,
    geometry_config: &GeometryConfig,
    serpentine_config: &SerpentineConfig,
    box_dims: (f64, f64),
    neighbor_info: Option<&[f64]>,
) -> OptimizationResult {
    let start_time = std::time::Instant::now();

    match serpentine_config.optimization_profile {
        OptimizationProfile::Fast => optimize_fast(
            p1,
            p2,
            geometry_config,
            serpentine_config,
            box_dims,
            neighbor_info,
            start_time,
        ),
        OptimizationProfile::Balanced => optimize_nelder_mead(
            p1,
            p2,
            geometry_config,
            serpentine_config,
            box_dims,
            neighbor_info,
            start_time,
        ),
        OptimizationProfile::Thorough => optimize_thorough(
            p1,
            p2,
            geometry_config,
            serpentine_config,
            box_dims,
            neighbor_info,
            start_time,
        ),
    }
}

/// Fast optimization using limited grid search
fn optimize_fast(
    p1: Point2D,
    p2: Point2D,
    geometry_config: &GeometryConfig,
    serpentine_config: &SerpentineConfig,
    box_dims: (f64, f64),
    neighbor_info: Option<&[f64]>,
    start_time: std::time::Instant,
) -> OptimizationResult {
    // Fast optimization with limited parameter exploration
    let min_clearance = geometry_config.wall_clearance;
    let channel_width = geometry_config.channel_width;

    // Get configurable parameter search ranges from constants registry
    let constants_reg = ConstantsRegistry::new();
    let wavelength_factors = constants_reg.get_fast_wavelength_factors();
    let wave_density_factors = constants_reg.get_fast_wave_density_factors();
    let fill_factors = constants_reg.get_fast_fill_factors();

    let mut best_result = OptimizationResult {
        params: OptimizationParams {
            wavelength_factor: serpentine_config.wavelength_factor,
            wave_density_factor: serpentine_config.wave_density_factor,
            fill_factor: serpentine_config.fill_factor,
        },
        path_length: 0.0,
        min_wall_distance: 0.0,
        min_neighbor_distance: 0.0,
        is_valid: false,
        iterations: 0,
        optimization_time: std::time::Duration::from_secs(0),
    };

    let mut iterations = 0;

    // Grid search over parameter combinations
    for &wavelength_factor in &wavelength_factors {
        for &wave_density_factor in &wave_density_factors {
            for &fill_factor in &fill_factors {
                iterations += 1;

                // Create test configuration without cloning the entire config
                let test_config = SerpentineConfig {
                    wavelength_factor,
                    wave_density_factor,
                    fill_factor,
                    gaussian_width_factor: serpentine_config.gaussian_width_factor,
                    wave_phase_direction: serpentine_config.wave_phase_direction,
                    wave_shape: serpentine_config.wave_shape,
                    optimization_enabled: false, // Disable nested optimization
                    target_fill_ratio: serpentine_config.target_fill_ratio,
                    optimization_profile: serpentine_config.optimization_profile,
                    adaptive_config: serpentine_config.adaptive_config,
                };

                // Generate test path using simplified serpentine generation logic
                let test_path = generate_simplified_serpentine_path(
                    p1,
                    p2,
                    geometry_config,
                    &test_config,
                    box_dims,
                    neighbor_info,
                );

                // Calculate metrics
                let path_length = calculate_path_length(&test_path);
                let min_wall_distance =
                    calculate_min_wall_distance(&test_path, box_dims, channel_width);
                let min_neighbor_distance = neighbor_info.map_or(f64::INFINITY, |neighbors| {
                    calculate_min_neighbor_distance(&test_path, neighbors, channel_width)
                });

                // Use penalty-based constraint handling for better optimization
                let penalty = calculate_constraint_penalty(
                    min_wall_distance,
                    min_neighbor_distance,
                    min_clearance,
                );
                let objective_score = path_length - penalty;

                // Update best result if this is better
                if objective_score > best_result.path_length {
                    let is_valid = penalty < constants::PENALTY_TOLERANCE;
                    best_result = OptimizationResult {
                        params: OptimizationParams {
                            wavelength_factor,
                            wave_density_factor,
                            fill_factor,
                        },
                        path_length: objective_score,
                        min_wall_distance,
                        min_neighbor_distance,
                        is_valid,
                        iterations,
                        optimization_time: start_time.elapsed(),
                    };
                }
            }
        }
    }

    // If no improvement found, return original parameters
    if best_result.path_length <= constants::MIN_PATH_LENGTH_THRESHOLD {
        let original_path = generate_simplified_serpentine_path(
            p1,
            p2,
            geometry_config,
            serpentine_config,
            box_dims,
            neighbor_info,
        );

        best_result = OptimizationResult {
            params: OptimizationParams {
                wavelength_factor: serpentine_config.wavelength_factor,
                wave_density_factor: serpentine_config.wave_density_factor,
                fill_factor: serpentine_config.fill_factor,
            },
            path_length: calculate_path_length(&original_path),
            min_wall_distance: calculate_min_wall_distance(&original_path, box_dims, channel_width),
            min_neighbor_distance: neighbor_info.map_or(f64::INFINITY, |neighbors| {
                calculate_min_neighbor_distance(&original_path, neighbors, channel_width)
            }),
            is_valid: true, // Assume original parameters are valid
            iterations,
            optimization_time: start_time.elapsed(),
        };
    }

    best_result
}

/// Balanced optimization using Nelder-Mead simplex algorithm
#[allow(clippy::needless_range_loop)]
fn optimize_nelder_mead(
    p1: Point2D,
    p2: Point2D,
    geometry_config: &GeometryConfig,
    serpentine_config: &SerpentineConfig,
    box_dims: (f64, f64),
    neighbor_info: Option<&[f64]>,
    start_time: std::time::Instant,
) -> OptimizationResult {
    // Start with current parameters as initial guess
    let initial_params = [
        serpentine_config.wavelength_factor,
        serpentine_config.wave_density_factor,
        serpentine_config.fill_factor,
    ];

    // Create initial simplex (triangle in 3D parameter space)
    let mut simplex = [
        initial_params,
        [
            initial_params[0] * constants::WAVELENGTH_PERTURBATION,
            initial_params[1],
            initial_params[2],
        ],
        [
            initial_params[0],
            initial_params[1] * constants::WAVE_DENSITY_PERTURBATION,
            initial_params[2],
        ],
        [
            initial_params[0],
            initial_params[1],
            initial_params[2] * constants::FILL_FACTOR_PERTURBATION,
        ],
    ];

    // Evaluate initial simplex
    let mut scores: Vec<f64> = simplex
        .iter()
        .map(|params| {
            evaluate_objective_function(
                *params,
                p1,
                p2,
                geometry_config,
                serpentine_config,
                box_dims,
                neighbor_info,
            )
        })
        .collect();

    let constants_reg = crate::config::ConstantsRegistry::new();
    let max_iterations = constants_reg.get_max_optimization_iterations();
    let tolerance = constants_reg.get_optimization_tolerance();
    let mut iterations = 0;

    // Nelder-Mead algorithm parameters
    let alpha = constants::REFLECTION_COEFFICIENT;
    let gamma = constants::EXPANSION_COEFFICIENT;
    let rho = constants::CONTRACTION_COEFFICIENT;
    let sigma = constants::SHRINK_COEFFICIENT;

    for _ in 0..max_iterations {
        iterations += 1;

        // Sort simplex by scores (best to worst)
        let mut indices: Vec<usize> = (0..simplex.len()).collect();
        indices.sort_by(|&a, &b| {
            scores[b]
                .partial_cmp(&scores[a])
                .unwrap_or(std::cmp::Ordering::Equal) // Handle NaN values gracefully
        });

        let best_idx = indices[0];
        let worst_idx = indices[indices.len() - 1];
        let second_worst_idx = indices[indices.len() - 2];

        // Check for convergence
        let score_range = scores[best_idx] - scores[worst_idx];
        if score_range < tolerance {
            break;
        }

        // Calculate centroid of all points except worst
        let mut centroid = [0.0; 3];
        for &idx in &indices[..indices.len() - 1] {
            for i in 0..3 {
                centroid[i] += simplex[idx][i];
            }
        }
        for val in &mut centroid {
            *val /= (simplex.len() - 1) as f64;
        }

        // Reflection
        let mut reflected = [0.0; 3];
        for i in 0..3 {
            reflected[i] = alpha.mul_add(centroid[i] - simplex[worst_idx][i], centroid[i]);
        }
        let reflected_score = evaluate_objective_function(
            reflected,
            p1,
            p2,
            geometry_config,
            serpentine_config,
            box_dims,
            neighbor_info,
        );

        if reflected_score > scores[second_worst_idx] && reflected_score <= scores[best_idx] {
            // Accept reflection
            simplex[worst_idx] = reflected;
            scores[worst_idx] = reflected_score;
        } else if reflected_score > scores[best_idx] {
            // Try expansion
            let mut expanded = [0.0; 3];
            for i in 0..3 {
                expanded[i] = gamma.mul_add(reflected[i] - centroid[i], centroid[i]);
            }
            let expanded_score = evaluate_objective_function(
                expanded,
                p1,
                p2,
                geometry_config,
                serpentine_config,
                box_dims,
                neighbor_info,
            );

            if expanded_score > reflected_score {
                simplex[worst_idx] = expanded;
                scores[worst_idx] = expanded_score;
            } else {
                simplex[worst_idx] = reflected;
                scores[worst_idx] = reflected_score;
            }
        } else {
            // Try contraction
            let mut contracted = [0.0; 3];
            for i in 0..3 {
                contracted[i] = rho.mul_add(simplex[worst_idx][i] - centroid[i], centroid[i]);
            }
            let contracted_score = evaluate_objective_function(
                contracted,
                p1,
                p2,
                geometry_config,
                serpentine_config,
                box_dims,
                neighbor_info,
            );

            if contracted_score > scores[worst_idx] {
                simplex[worst_idx] = contracted;
                scores[worst_idx] = contracted_score;
            } else {
                // Shrink simplex
                for i in 1..simplex.len() {
                    for j in 0..3 {
                        simplex[i][j] = sigma
                            .mul_add(simplex[i][j] - simplex[best_idx][j], simplex[best_idx][j]);
                    }
                    scores[i] = evaluate_objective_function(
                        simplex[i],
                        p1,
                        p2,
                        geometry_config,
                        serpentine_config,
                        box_dims,
                        neighbor_info,
                    );
                }
            }
        }
    }

    // Find best result
    let best_idx = scores
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| {
            a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal) // Handle NaN values gracefully
        })
        .map_or(0, |(idx, _)| idx); // Default to first element if no maximum found

    let best_params = simplex[best_idx];
    let best_config = SerpentineConfig {
        wavelength_factor: best_params[0],
        wave_density_factor: best_params[1],
        fill_factor: best_params[2],
        ..*serpentine_config
    };

    // Generate final path and calculate metrics
    let final_path = generate_simplified_serpentine_path(
        p1,
        p2,
        geometry_config,
        &best_config,
        box_dims,
        neighbor_info,
    );

    let path_length = calculate_path_length(&final_path);
    let min_wall_distance =
        calculate_min_wall_distance(&final_path, box_dims, geometry_config.channel_width);
    let min_neighbor_distance = if let Some(neighbors) = neighbor_info {
        calculate_min_neighbor_distance(&final_path, neighbors, geometry_config.channel_width)
    } else {
        f64::INFINITY
    };

    OptimizationResult {
        params: OptimizationParams {
            wavelength_factor: best_params[0],
            wave_density_factor: best_params[1],
            fill_factor: best_params[2],
        },
        path_length,
        min_wall_distance,
        min_neighbor_distance,
        is_valid: min_wall_distance >= geometry_config.wall_clearance
            && min_neighbor_distance >= geometry_config.wall_clearance,
        iterations,
        optimization_time: start_time.elapsed(),
    }
}

/// Evaluate objective function for optimization (length - penalties)
#[must_use]
fn evaluate_objective_function(
    params: [f64; 3],
    p1: Point2D,
    p2: Point2D,
    geometry_config: &GeometryConfig,
    serpentine_config: &SerpentineConfig,
    box_dims: (f64, f64),
    neighbor_info: Option<&[f64]>,
) -> f64 {
    // Clamp parameters to valid ranges
    let wavelength_factor = params[0].clamp(
        constants::MIN_WAVELENGTH_FACTOR,
        constants::MAX_WAVELENGTH_FACTOR,
    );
    let wave_density_factor = params[1].clamp(
        constants::MIN_WAVE_DENSITY_FACTOR,
        constants::MAX_WAVE_DENSITY_FACTOR,
    );
    let fill_factor = params[2].clamp(constants::MIN_FILL_FACTOR, constants::MAX_FILL_FACTOR);

    let test_config = SerpentineConfig {
        fill_factor,
        wavelength_factor,
        wave_density_factor,
        ..*serpentine_config
    };

    // Generate test path
    let test_path = generate_simplified_serpentine_path(
        p1,
        p2,
        geometry_config,
        &test_config,
        box_dims,
        neighbor_info,
    );

    // Calculate metrics
    let path_length = calculate_path_length(&test_path);
    let min_wall_distance =
        calculate_min_wall_distance(&test_path, box_dims, geometry_config.channel_width);
    let min_neighbor_distance = if let Some(neighbors) = neighbor_info {
        calculate_min_neighbor_distance(&test_path, neighbors, geometry_config.channel_width)
    } else {
        f64::INFINITY
    };

    // Calculate penalty
    let penalty = calculate_constraint_penalty(
        min_wall_distance,
        min_neighbor_distance,
        geometry_config.wall_clearance,
    );

    // Return objective score (maximize length, minimize penalty)
    path_length - penalty
}

/// Thorough optimization using multi-start Nelder-Mead with extensive exploration
fn optimize_thorough(
    p1: Point2D,
    p2: Point2D,
    geometry_config: &GeometryConfig,
    serpentine_config: &SerpentineConfig,
    box_dims: (f64, f64),
    neighbor_info: Option<&[f64]>,
    start_time: std::time::Instant,
) -> OptimizationResult {
    let mut best_result = OptimizationResult {
        params: OptimizationParams {
            wavelength_factor: serpentine_config.wavelength_factor,
            wave_density_factor: serpentine_config.wave_density_factor,
            fill_factor: serpentine_config.fill_factor,
        },
        path_length: constants::MIN_PATH_LENGTH_THRESHOLD,
        min_wall_distance: constants::MIN_PATH_LENGTH_THRESHOLD,
        min_neighbor_distance: constants::MIN_PATH_LENGTH_THRESHOLD,
        is_valid: false,
        iterations: 0,
        optimization_time: std::time::Duration::from_secs(0),
    };

    // Multiple starting points for thorough exploration
    let starting_points = vec![
        [
            serpentine_config.wavelength_factor,
            serpentine_config.wave_density_factor,
            serpentine_config.fill_factor,
        ],
        [1.0, 1.0, 0.7],
        [2.0, 2.0, 0.8],
        [3.0, 3.0, 0.9],
        [4.0, 1.5, 0.85],
        [1.5, 4.0, 0.75],
    ];

    let mut total_iterations = 0;

    for start_point in starting_points {
        // Create temporary config for this starting point
        let temp_config = SerpentineConfig {
            wavelength_factor: start_point[0],
            wave_density_factor: start_point[1],
            fill_factor: start_point[2],
            ..*serpentine_config
        };

        // Run Nelder-Mead from this starting point
        let result = optimize_nelder_mead(
            p1,
            p2,
            geometry_config,
            &temp_config,
            box_dims,
            neighbor_info,
            start_time,
        );

        total_iterations += result.iterations;

        // Keep the best result
        if result.path_length > best_result.path_length {
            best_result = result;
        }
    }

    // Update total iterations
    best_result.iterations = total_iterations;
    best_result.optimization_time = start_time.elapsed();

    best_result
}
