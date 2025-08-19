//! Field abstraction for CFD grid data
//!
//! This module provides a clean abstraction for grid-based field data,
//! separating data storage from solver logic.

use nalgebra::{RealField, Vector2};
use std::ops::{Index, IndexMut};

/// A 2D field representing scalar or vector data on a grid
#[derive(Debug, Clone)]
pub struct Field<T> {
    /// Field data stored in row-major order for cache efficiency
    data: Vec<T>,
    /// Number of cells in x-direction
    nx: usize,
    /// Number of cells in y-direction
    ny: usize,
}

impl<T: Clone> Field<T> {
    /// Create a new field with given dimensions and initial value
    pub fn new(nx: usize, ny: usize, initial: T) -> Self {
        Self {
            data: vec![initial; nx * ny],
            nx,
            ny,
        }
    }

    /// Get field dimensions
    pub fn dimensions(&self) -> (usize, usize) {
        (self.nx, self.ny)
    }

    /// Get value at grid point (i, j)
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> &T {
        debug_assert!(i < self.nx && j < self.ny);
        &self.data[i * self.ny + j]
    }

    /// Get mutable value at grid point (i, j)
    #[inline]
    pub fn get_mut(&mut self, i: usize, j: usize) -> &mut T {
        debug_assert!(i < self.nx && j < self.ny);
        &mut self.data[i * self.ny + j]
    }

    /// Set value at grid point (i, j)
    #[inline]
    pub fn set(&mut self, i: usize, j: usize, value: T) {
        debug_assert!(i < self.nx && j < self.ny);
        self.data[i * self.ny + j] = value;
    }

    /// Iterate over all field values
    pub fn iter(&self) -> impl Iterator<Item = &T> {
        self.data.iter()
    }

    /// Iterate mutably over all field values
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut T> {
        self.data.iter_mut()
    }

    /// Apply a function to all field values
    pub fn apply<F>(&mut self, f: F)
    where
        F: Fn(&mut T),
    {
        for value in &mut self.data {
            f(value);
        }
    }

    /// Create a new field by mapping this field's values
    pub fn map<U, F>(&self, f: F) -> Field<U>
    where
        U: Clone,
        F: Fn(&T) -> U,
    {
        Field {
            data: self.data.iter().map(f).collect(),
            nx: self.nx,
            ny: self.ny,
        }
    }
}

/// Implement indexing for convenient access
impl<T> Index<(usize, usize)> for Field<T> {
    type Output = T;

    fn index(&self, (i, j): (usize, usize)) -> &Self::Output {
        self.get(i, j)
    }
}

impl<T> IndexMut<(usize, usize)> for Field<T> {
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut Self::Output {
        self.get_mut(i, j)
    }
}

/// Scalar field type alias
pub type ScalarField<T> = Field<T>;

/// Vector field type alias
pub type VectorField<T> = Field<Vector2<T>>;

/// Persistent solver state containing primary variables
#[derive(Debug, Clone)]
pub struct SolverState<T: RealField> {
    /// Velocity field
    pub velocity: VectorField<T>,
    /// Pressure field
    pub pressure: ScalarField<T>,
    /// Previous timestep velocity (for transient simulations)
    pub velocity_old: Option<VectorField<T>>,
}

impl<T: RealField> SolverState<T> {
    /// Create new solver state
    pub fn new(nx: usize, ny: usize) -> Self {
        Self {
            velocity: Field::new(nx, ny, Vector2::zeros()),
            pressure: Field::new(nx, ny, T::zero()),
            velocity_old: None,
        }
    }

    /// Initialize for transient simulation
    pub fn init_transient(&mut self) {
        if self.velocity_old.is_none() {
            self.velocity_old = Some(self.velocity.clone());
        }
    }

    /// Swap velocity fields for next timestep
    pub fn advance_timestep(&mut self) {
        if let Some(ref mut old) = self.velocity_old {
            std::mem::swap(&mut self.velocity, old);
        }
    }
}

/// Temporary workspace for PISO algorithm
pub struct PisoWorkspace<T: RealField> {
    /// First pressure correction
    pub p_prime: ScalarField<T>,
    /// Second pressure correction
    pub p_double_prime: ScalarField<T>,
    /// First velocity correction
    pub u_prime: VectorField<T>,
    /// Second velocity correction
    pub u_double_prime: VectorField<T>,
    /// Momentum equation coefficients
    pub coefficients: MomentumCoefficients<T>,
}

/// Momentum equation discretization coefficients
pub struct MomentumCoefficients<T: RealField> {
    /// Central coefficient
    pub a_p: ScalarField<T>,
    /// East neighbor coefficient
    pub a_e: ScalarField<T>,
    /// West neighbor coefficient
    pub a_w: ScalarField<T>,
    /// North neighbor coefficient
    pub a_n: ScalarField<T>,
    /// South neighbor coefficient
    pub a_s: ScalarField<T>,
}

impl<T: RealField> MomentumCoefficients<T> {
    /// Create new coefficient fields
    pub fn new(nx: usize, ny: usize) -> Self {
        Self {
            a_p: Field::new(nx, ny, T::zero()),
            a_e: Field::new(nx, ny, T::zero()),
            a_w: Field::new(nx, ny, T::zero()),
            a_n: Field::new(nx, ny, T::zero()),
            a_s: Field::new(nx, ny, T::zero()),
        }
    }
}

impl<T: RealField> PisoWorkspace<T> {
    /// Create workspace for given grid dimensions
    pub fn new(nx: usize, ny: usize) -> Self {
        Self {
            p_prime: Field::new(nx, ny, T::zero()),
            p_double_prime: Field::new(nx, ny, T::zero()),
            u_prime: Field::new(nx, ny, Vector2::zeros()),
            u_double_prime: Field::new(nx, ny, Vector2::zeros()),
            coefficients: MomentumCoefficients::new(nx, ny),
        }
    }

    /// Clear workspace for next timestep
    pub fn clear(&mut self) {
        self.p_prime.apply(|v| *v = T::zero());
        self.p_double_prime.apply(|v| *v = T::zero());
        self.u_prime.apply(|v| *v = Vector2::zeros());
        self.u_double_prime.apply(|v| *v = Vector2::zeros());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_field_creation() {
        let field: Field<f64> = Field::new(10, 10, 0.0);
        assert_eq!(field.dimensions(), (10, 10));
    }

    #[test]
    fn test_field_indexing() {
        let mut field = Field::new(5, 5, 0.0);
        field[(2, 3)] = 1.5;
        assert_eq!(field[(2, 3)], 1.5);
    }

    #[test]
    fn test_solver_state() {
        let mut state: SolverState<f64> = SolverState::new(10, 10);
        state.init_transient();
        assert!(state.velocity_old.is_some());
    }
}