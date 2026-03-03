//! Nelder-Mead simplex optimization for serpentine channel parameters
//!
//! Implements the Nelder-Mead simplex algorithm for balanced optimization
//! of serpentine channel parameters (wavelength factor, wave density factor,
//! fill factor). Used directly for the Balanced profile and as a subroutine
//! for the Thorough multi-start profile.

use crate::config::GeometryConfig;
use crate::config::SerpentineConfig;
use crate::geometry::types::Point2D;

use super::super::constants;
use super::super::serpentine_eval::generate_simplified_serpentine_path;
use super::super::{
    calculate_min_neighbor_distance, calculate_min_wall_distance, calculate_path_length,
    OptimizationParams, OptimizationResult,
};
use super::evaluate_objective_function;

/// Balanced optimization using Nelder-Mead simplex algorithm
#[allow(clippy::needless_range_loop)]
pub(super) fn optimize_nelder_mead(
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
