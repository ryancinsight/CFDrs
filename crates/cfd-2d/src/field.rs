//! Field abstraction for grid-based data
//!
//! This module provides a clean abstraction for field data on structured grids,
//! separating data storage from solver logic.

use nalgebra::{RealField, Vector2};
use serde::{Deserialize, Serialize};
use std::ops::{Index, IndexMut};

/// A scalar field on a 2D structured grid
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScalarField<T: RealField> {
    /// Field data stored as [nx][ny]
    data: Vec<Vec<T>>,
    /// Grid dimensions
    nx: usize,
    ny: usize,
}

impl<T: RealField> ScalarField<T> {
    /// Create a new scalar field with given dimensions
    pub fn new(nx: usize, ny: usize) -> Self {
        Self {
            data: vec![vec![T::zero(); ny]; nx],
            nx,
            ny,
        }
    }

    /// Create a scalar field initialized with a constant value
    pub fn from_value(nx: usize, ny: usize, value: T) -> Self {
        Self {
            data: vec![vec![value; ny]; nx],
            nx,
            ny,
        }
    }

    /// Get field dimensions
    pub fn dimensions(&self) -> (usize, usize) {
        (self.nx, self.ny)
    }

    /// Get immutable reference to underlying data
    pub fn data(&self) -> &Vec<Vec<T>> {
        &self.data
    }

    /// Get mutable reference to underlying data
    pub fn data_mut(&mut self) -> &mut Vec<Vec<T>> {
        &mut self.data
    }

    /// Fill field with a constant value
    pub fn fill(&mut self, value: T) {
        for i in 0..self.nx {
            for j in 0..self.ny {
                self.data[i][j] = value;
            }
        }
    }

    /// Copy data from another field
    pub fn copy_from(&mut self, other: &Self) {
        assert_eq!(self.nx, other.nx);
        assert_eq!(self.ny, other.ny);
        for i in 0..self.nx {
            for j in 0..self.ny {
                self.data[i][j] = other.data[i][j];
            }
        }
    }
}

impl<T: RealField> Index<(usize, usize)> for ScalarField<T> {
    type Output = T;

    fn index(&self, (i, j): (usize, usize)) -> &Self::Output {
        &self.data[i][j]
    }
}

impl<T: RealField> IndexMut<(usize, usize)> for ScalarField<T> {
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut Self::Output {
        &mut self.data[i][j]
    }
}

/// A vector field on a 2D structured grid
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VectorField<T: RealField> {
    /// Field data stored as [nx][ny]
    data: Vec<Vec<Vector2<T>>>,
    /// Grid dimensions
    nx: usize,
    ny: usize,
}

impl<T: RealField> VectorField<T> {
    /// Create a new vector field with given dimensions
    pub fn new(nx: usize, ny: usize) -> Self {
        Self {
            data: vec![vec![Vector2::zeros(); ny]; nx],
            nx,
            ny,
        }
    }

    /// Create a vector field initialized with a constant value
    pub fn from_value(nx: usize, ny: usize, value: Vector2<T>) -> Self {
        Self {
            data: vec![vec![value; ny]; nx],
            nx,
            ny,
        }
    }

    /// Get field dimensions
    pub fn dimensions(&self) -> (usize, usize) {
        (self.nx, self.ny)
    }

    /// Get immutable reference to underlying data
    pub fn data(&self) -> &Vec<Vec<Vector2<T>>> {
        &self.data
    }

    /// Get mutable reference to underlying data
    pub fn data_mut(&mut self) -> &mut Vec<Vec<Vector2<T>>> {
        &mut self.data
    }

    /// Fill field with a constant value
    pub fn fill(&mut self, value: Vector2<T>) {
        for i in 0..self.nx {
            for j in 0..self.ny {
                self.data[i][j] = value;
            }
        }
    }

    /// Copy data from another field
    pub fn copy_from(&mut self, other: &Self) {
        assert_eq!(self.nx, other.nx);
        assert_eq!(self.ny, other.ny);
        for i in 0..self.nx {
            for j in 0..self.ny {
                self.data[i][j] = other.data[i][j];
            }
        }
    }
}

impl<T: RealField> Index<(usize, usize)> for VectorField<T> {
    type Output = Vector2<T>;

    fn index(&self, (i, j): (usize, usize)) -> &Self::Output {
        &self.data[i][j]
    }
}

impl<T: RealField> IndexMut<(usize, usize)> for VectorField<T> {
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut Self::Output {
        &mut self.data[i][j]
    }
}

/// Persistent solver state containing velocity and pressure fields
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SolverState<T: RealField> {
    /// Velocity field
    pub velocity: VectorField<T>,
    /// Pressure field
    pub pressure: ScalarField<T>,
    /// Previous time step velocity (for transient simulations)
    pub velocity_old: VectorField<T>,
}

impl<T: RealField> SolverState<T> {
    /// Create a new solver state with given dimensions
    pub fn new(nx: usize, ny: usize) -> Self {
        Self {
            velocity: VectorField::new(nx, ny),
            pressure: ScalarField::new(nx, ny),
            velocity_old: VectorField::new(nx, ny),
        }
    }

    /// Initialize state with uniform values
    pub fn initialize(&mut self, velocity: Vector2<T>, pressure: T) {
        self.velocity.fill(velocity);
        self.pressure.fill(pressure);
        self.velocity_old.fill(velocity);
    }

    /// Swap old and current velocity fields (for time stepping)
    pub fn swap_velocities(&mut self) {
        std::mem::swap(&mut self.velocity, &mut self.velocity_old);
    }
}

/// Temporary workspace for PISO solver computations
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
    pub a_p: ScalarField<T>,
    pub a_e: ScalarField<T>,
    pub a_w: ScalarField<T>,
    pub a_n: ScalarField<T>,
    pub a_s: ScalarField<T>,
    /// Source terms
    pub b_x: ScalarField<T>,
    pub b_y: ScalarField<T>,
}

impl<T: RealField> PisoWorkspace<T> {
    /// Create a new workspace with given dimensions
    pub fn new(nx: usize, ny: usize) -> Self {
        Self {
            p_prime: ScalarField::new(nx, ny),
            p_double_prime: ScalarField::new(nx, ny),
            u_prime: VectorField::new(nx, ny),
            u_double_prime: VectorField::new(nx, ny),
            a_p: ScalarField::new(nx, ny),
            a_e: ScalarField::new(nx, ny),
            a_w: ScalarField::new(nx, ny),
            a_n: ScalarField::new(nx, ny),
            a_s: ScalarField::new(nx, ny),
            b_x: ScalarField::new(nx, ny),
            b_y: ScalarField::new(nx, ny),
        }
    }

    /// Clear all workspace fields
    pub fn clear(&mut self) {
        self.p_prime.fill(T::zero());
        self.p_double_prime.fill(T::zero());
        self.u_prime.fill(Vector2::zeros());
        self.u_double_prime.fill(Vector2::zeros());
        // Note: Coefficients are recalculated each step, no need to clear
    }
}