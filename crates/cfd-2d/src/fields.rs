//! Simulation field structures for 2D CFD calculations.
//!
//! This module implements efficient field storage using flattened vectors
//! for better cache locality and performance.

use nalgebra::{RealField, Vector2};
use num_traits::FromPrimitive;

/// Efficient 2D field storage using flattened vector for cache locality
#[derive(Debug, Clone)]
pub struct Field2D<T> {
    data: Vec<T>,
    nx: usize,
    ny: usize,
}

impl<T: Clone> Field2D<T> {
    /// Create new field with default value
    pub fn new(nx: usize, ny: usize, default_value: T) -> Self {
        Self {
            data: vec![default_value; nx * ny],
            nx,
            ny,
        }
    }

    /// Get immutable reference to element at (i, j)
    #[inline]
    pub fn at(&self, i: usize, j: usize) -> &T {
        debug_assert!(i < self.nx && j < self.ny, "Index out of bounds");
        &self.data[j * self.nx + i]
    }

    /// Get mutable reference to element at (i, j)
    #[inline]
    pub fn at_mut(&mut self, i: usize, j: usize) -> &mut T {
        debug_assert!(i < self.nx && j < self.ny, "Index out of bounds");
        &mut self.data[j * self.nx + i]
    }

    /// Get grid dimensions
    pub fn dimensions(&self) -> (usize, usize) {
        (self.nx, self.ny)
    }

    /// Get raw data slice for efficient operations
    pub fn data(&self) -> &[T] {
        &self.data
    }

    /// Get mutable raw data slice
    pub fn data_mut(&mut self) -> &mut [T] {
        &mut self.data
    }
}

/// Complete simulation field state for SIMPLE algorithm
#[derive(Debug, Clone)]
pub struct SimulationFields<T: RealField> {
    /// Current velocity field
    pub u: Field2D<Vector2<T>>,
    /// Predicted velocity field (momentum solution)
    pub u_star: Field2D<Vector2<T>>,
    /// Pressure field
    pub p: Field2D<T>,
    /// Pressure correction field
    pub p_prime: Field2D<T>,
    /// Grid dimensions
    pub nx: usize,
    pub ny: usize,
}

impl<T: RealField + FromPrimitive> SimulationFields<T> {
    /// Create new simulation fields with zero initialization
    pub fn new(nx: usize, ny: usize) -> Self {
        Self {
            u: Field2D::new(nx, ny, Vector2::zeros()),
            u_star: Field2D::new(nx, ny, Vector2::zeros()),
            p: Field2D::new(nx, ny, T::zero()),
            p_prime: Field2D::new(nx, ny, T::zero()),
            nx,
            ny,
        }
    }

    /// Reset all fields to zero
    pub fn reset(&mut self) {
        for u_val in self.u.data_mut() {
            *u_val = Vector2::zeros();
        }
        for u_star_val in self.u_star.data_mut() {
            *u_star_val = Vector2::zeros();
        }
        for p_val in self.p.data_mut() {
            *p_val = T::zero();
        }
        for p_prime_val in self.p_prime.data_mut() {
            *p_prime_val = T::zero();
        }
    }

    /// Get maximum velocity magnitude for stability analysis
    pub fn max_velocity_magnitude(&self) -> T {
        self.u.data()
            .iter()
            .map(|v| v.magnitude())
            .fold(T::zero(), |acc, mag| if mag > acc { mag } else { acc })
    }
}