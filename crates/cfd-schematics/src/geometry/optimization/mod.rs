//! Optimization utilities for serpentine channel generation
//!
//! This module provides utilities for optimizing serpentine channel parameters
//! to maximize channel length while maintaining proper wall clearance and
//! multi-channel compatibility.
//!
//! # Submodules
//!
//! - [`algorithms`]: Grid-search, Nelder-Mead, and multi-start optimizers
//! - [`serpentine_eval`]: Simplified serpentine path generation for evaluation

mod algorithms;
mod serpentine_eval;

use crate::geometry::types::Point2D;

// Re-export public API
pub use algorithms::optimize_serpentine_parameters;

/// Optimization algorithm constants
pub(super) mod constants {
    /// Minimum path length threshold for valid optimization results
    pub const MIN_PATH_LENGTH_THRESHOLD: f64 = 0.0;

    /// Penalty multiplier for constraint violations
    pub const CONSTRAINT_PENALTY_MULTIPLIER: f64 = 1000.0;

    /// Small penalty tolerance for validity checking
    pub const PENALTY_TOLERANCE: f64 = 1.0;

    /// Nelder-Mead algorithm coefficients
    pub const REFLECTION_COEFFICIENT: f64 = 1.0;
    pub const EXPANSION_COEFFICIENT: f64 = 2.0;
    pub const CONTRACTION_COEFFICIENT: f64 = 0.5;
    pub const SHRINK_COEFFICIENT: f64 = 0.5;

    /// Parameter bounds for optimization
    pub const MIN_WAVELENGTH_FACTOR: f64 = 0.5;
    pub const MAX_WAVELENGTH_FACTOR: f64 = 5.0;
    pub const MIN_WAVE_DENSITY_FACTOR: f64 = 0.5;
    pub const MAX_WAVE_DENSITY_FACTOR: f64 = 5.0;
    pub const MIN_FILL_FACTOR: f64 = 0.1;
    pub const MAX_FILL_FACTOR: f64 = 0.95;

    /// Simplex initialization perturbation factors
    pub const WAVELENGTH_PERTURBATION: f64 = 1.1;
    pub const WAVE_DENSITY_PERTURBATION: f64 = 1.1;
    pub const FILL_FACTOR_PERTURBATION: f64 = 1.05;

    /// Wave shape parameters
    pub const SQUARE_WAVE_SHARPNESS: f64 = 5.0;

    /// Envelope calculation constants
    pub const SMOOTHSTEP_COEFFICIENT_1: f64 = 3.0;
    pub const SMOOTHSTEP_COEFFICIENT_2: f64 = 2.0;
    pub const GAUSSIAN_CENTER: f64 = 0.5;

    /// Distance normalization bounds
    pub const MIN_DISTANCE_NORMALIZATION: f64 = 0.1;
    pub const MAX_DISTANCE_NORMALIZATION: f64 = 1.0;
}

/// Calculate the total path length of a serpentine channel
///
/// # Arguments
/// * `path` - Vector of points defining the serpentine path
///
/// # Returns
/// Total length of the path by summing Euclidean distances between consecutive points
#[must_use]
pub fn calculate_path_length(path: &[Point2D]) -> f64 {
    if path.len() < 2 {
        return 0.0;
    }

    path.windows(2)
        .map(|window| {
            let (p1, p2) = (window[0], window[1]);
            let dx = p2.0 - p1.0;
            let dy = p2.1 - p1.1;
            dx.hypot(dy)
        })
        .sum()
}

/// Calculate the minimum distance from a path to the box boundaries
///
/// # Arguments
/// * `path` - Vector of points defining the channel path
/// * `box_dims` - Box dimensions (width, height)
/// * `channel_width` - Width of the channel (for clearance calculation)
///
/// # Returns
/// Minimum distance from any path point to the nearest wall, considering channel width
#[must_use]
pub fn calculate_min_wall_distance(
    path: &[Point2D],
    box_dims: (f64, f64),
    channel_width: f64,
) -> f64 {
    let (box_width, box_height) = box_dims;
    let half_channel_width = channel_width / 2.0;

    path.iter()
        .map(|&(x, y)| {
            // Distance to each wall, accounting for channel width
            let dist_to_left = x - half_channel_width;
            let dist_to_right = box_width - x - half_channel_width;
            let dist_to_bottom = y - half_channel_width;
            let dist_to_top = box_height - y - half_channel_width;

            // Return minimum distance to any wall
            dist_to_left
                .min(dist_to_right)
                .min(dist_to_bottom)
                .min(dist_to_top)
        })
        .fold(f64::INFINITY, f64::min)
}

/// Calculate the minimum distance between a path and neighboring channels
///
/// # Arguments
/// * `path` - Vector of points defining the channel path
/// * `neighbor_positions` - Y-coordinates of neighboring channels
/// * `channel_width` - Width of the channel
///
/// # Returns
/// Minimum distance to any neighboring channel, considering channel width
#[must_use]
pub fn calculate_min_neighbor_distance(
    path: &[Point2D],
    neighbor_positions: &[f64],
    channel_width: f64,
) -> f64 {
    if neighbor_positions.is_empty() {
        return f64::INFINITY;
    }

    path.iter()
        .map(|&(_, y)| {
            neighbor_positions
                .iter()
                .map(|&neighbor_y| (y - neighbor_y).abs() - channel_width)
                .fold(f64::INFINITY, f64::min)
        })
        .fold(f64::INFINITY, f64::min)
}

/// Optimization parameters for serpentine channel generation
///
/// These parameters control the shape and density of serpentine channels
/// during the optimization process.
#[derive(Debug, Clone)]
pub struct OptimizationParams {
    /// Multiplier for channel width to determine wavelength (1.0 to 10.0)
    pub wavelength_factor: f64,
    /// Controls how many waves appear relative to channel length (0.5 to 5.0)
    pub wave_density_factor: f64,
    /// Fraction of available vertical space to fill (0.1 to 0.95)
    pub fill_factor: f64,
}

/// Result of serpentine optimization
///
/// Contains the optimized parameters and performance metrics from
/// the optimization process.
#[derive(Debug, Clone)]
pub struct OptimizationResult {
    /// The optimized parameters that produced the best result
    pub params: OptimizationParams,
    /// Total length of the optimized serpentine path
    pub path_length: f64,
    /// Minimum distance to any wall boundary
    pub min_wall_distance: f64,
    /// Minimum distance to any neighboring channel
    pub min_neighbor_distance: f64,
    /// Whether the optimization result meets all constraints
    pub is_valid: bool,
    /// Number of optimization iterations performed
    pub iterations: usize,
    /// Total time spent on optimization
    pub optimization_time: std::time::Duration,
}

/// Optimization statistics for monitoring performance
///
/// Provides detailed metrics about the optimization process for
/// performance analysis and debugging.
#[derive(Debug, Clone)]
pub struct OptimizationStats {
    /// Total number of parameter evaluations performed
    pub total_evaluations: usize,
    /// Number of cache hits during optimization
    pub cache_hits: usize,
    /// Number of cache misses during optimization
    pub cache_misses: usize,
    /// Best objective function score achieved
    pub best_score: f64,
    /// Number of iterations until convergence
    pub convergence_iterations: usize,
}

/// Calculate penalty for constraint violations
#[must_use]
pub(super) fn calculate_constraint_penalty(
    wall_distance: f64,
    neighbor_distance: f64,
    min_clearance: f64,
) -> f64 {
    let mut penalty = 0.0;

    // Heavy penalty for wall clearance violations
    if wall_distance < min_clearance {
        penalty += (min_clearance - wall_distance) * constants::CONSTRAINT_PENALTY_MULTIPLIER;
    }

    // Heavy penalty for neighbor clearance violations
    if neighbor_distance < min_clearance {
        penalty += (min_clearance - neighbor_distance) * constants::CONSTRAINT_PENALTY_MULTIPLIER;
    }

    penalty
}
