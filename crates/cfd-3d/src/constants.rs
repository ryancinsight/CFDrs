//! Named constants for 3D CFD algorithms
//!
//! This module provides a single source of truth for all numerical constants
//! used in 3D CFD simulations, following the SSOT principle.

use nalgebra::RealField;

// Numerical tolerances
/// Machine epsilon for floating point comparisons
pub const EPSILON_F64: f64 = 1e-10;
/// Small value for avoiding division by zero
pub const SMALL_NUMBER: f64 = 1e-14;
/// Default convergence tolerance
pub const DEFAULT_TOLERANCE: f64 = 1e-6;

// Mesh quality thresholds
/// Minimum acceptable mesh quality (0-1 scale)
pub const MIN_MESH_QUALITY: f64 = 0.1;
/// Warning threshold for mesh quality
pub const MESH_QUALITY_WARNING: f64 = 0.3;
/// Good mesh quality threshold
pub const GOOD_MESH_QUALITY: f64 = 0.7;

// Geometric constants
/// Factor for tetrahedron volume calculation (1/6)
pub const TETRAHEDRON_VOLUME_FACTOR: f64 = 6.0;
/// Number of vertices in a tetrahedron
pub const TETRAHEDRON_VERTICES: usize = 4;
/// Number of faces in a tetrahedron
pub const TETRAHEDRON_FACES: usize = 4;
/// Number of edges in a tetrahedron
pub const TETRAHEDRON_EDGES: usize = 6;

// Golden ratio for icosahedron generation
/// Golden ratio value
pub const GOLDEN_RATIO: f64 = 1.618033988749895; // (1 + sqrt(5)) / 2

// Format precision for mesh vertex hashing
/// Number of decimal places for vertex coordinate hashing
pub const VERTEX_HASH_PRECISION: usize = 6;

// Ray-triangle intersection tolerance
/// Minimum value for valid ray-triangle intersection
pub const RAY_INTERSECTION_EPSILON: f64 = 1e-10;

/// Get epsilon value for type T
pub fn epsilon<T: RealField>() -> T {
    T::from_f64(EPSILON_F64).unwrap_or_else(T::default_epsilon)
}

/// Get small number for type T
pub fn small_number<T: RealField>() -> T {
    T::from_f64(SMALL_NUMBER).unwrap_or_else(T::default_epsilon)
}

/// Get default tolerance for type T
pub fn default_tolerance<T: RealField>() -> T {
    T::from_f64(DEFAULT_TOLERANCE).unwrap_or_else(|| T::default_epsilon() * T::from_f64(100.0).unwrap())
}

/// Get minimum mesh quality threshold for type T
pub fn min_mesh_quality<T: RealField>() -> T {
    T::from_f64(MIN_MESH_QUALITY).unwrap()
}

/// Get tetrahedron volume factor for type T
pub fn tetrahedron_volume_factor<T: RealField>() -> T {
    T::from_f64(TETRAHEDRON_VOLUME_FACTOR).unwrap()
}

/// Get golden ratio for type T
pub fn golden_ratio<T: RealField>() -> T {
    T::from_f64(GOLDEN_RATIO).unwrap()
}

/// Get ray intersection epsilon for type T
pub fn ray_intersection_epsilon<T: RealField>() -> T {
    T::from_f64(RAY_INTERSECTION_EPSILON).unwrap_or_else(T::default_epsilon)
}