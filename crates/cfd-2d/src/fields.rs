//! Simulation field structures for 2D CFD calculations.
//!
//! This module implements efficient field storage using flattened vectors
//! for better cache locality and performance, following SSOT and zero-copy principles.

use nalgebra::{RealField, Vector2};
use num_traits::FromPrimitive;
use cfd_core::Fluid;

/// Constants for field operations
pub mod constants {
    /// Default Reynolds number for laminar flow
    pub const DEFAULT_REYNOLDS: f64 = 100.0;
    
    /// Minimum allowed time step for stability
    pub const MIN_TIME_STEP: f64 = 1e-6;
    
    /// Maximum allowed time step
    pub const MAX_TIME_STEP: f64 = 1.0;
    
    /// Convergence tolerance for iterative methods
    pub const CONVERGENCE_TOLERANCE: f64 = 1e-6;
}

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

    /// Get value at element (i, j) - returns value for Copy types
    #[inline]
    pub fn at(&self, i: usize, j: usize) -> T 
    where 
        T: Copy
    {
        debug_assert!(i < self.nx && j < self.ny, "Index out of bounds");
        self.data[j * self.nx + i]
    }
    
    /// Get immutable reference to element at (i, j)
    #[inline]
    pub fn at_ref(&self, i: usize, j: usize) -> &T {
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
    
    /// Apply a function to all elements using iterator
    pub fn map_inplace<F>(&mut self, f: F) 
    where 
        F: Fn(&mut T)
    {
        self.data.iter_mut().for_each(f);
    }
    
    /// Create a new field by mapping a function over elements
    pub fn map<U, F>(&self, f: F) -> Field2D<U>
    where
        U: Clone,
        F: Fn(&T) -> U
    {
        Field2D {
            data: self.data.iter().map(f).collect(),
            nx: self.nx,
            ny: self.ny,
        }
    }
}

/// Complete simulation field state with proper abstraction for pressure-velocity coupling
/// 
/// This structure provides both vector and component access to velocity fields,
/// supporting the requirements of SIMPLE, PISO, and other algorithms.
#[derive(Debug, Clone)]
pub struct SimulationFields<T: RealField> {
    /// X-component of velocity field
    pub u: Field2D<T>,
    /// Y-component of velocity field  
    pub v: Field2D<T>,
    /// Predicted X-velocity (momentum solution)
    pub u_star: Field2D<T>,
    /// Predicted Y-velocity (momentum solution)
    pub v_star: Field2D<T>,
    /// Pressure field
    pub p: Field2D<T>,
    /// Pressure correction field
    pub p_prime: Field2D<T>,
    /// Density field (can be constant or variable)
    pub density: Field2D<T>,
    /// Dynamic viscosity field
    pub viscosity: Field2D<T>,
    /// Grid dimensions
    pub nx: usize,
    pub ny: usize,
}

impl<T: RealField + FromPrimitive> SimulationFields<T> {
    /// Create new simulation fields with zero initialization
    pub fn new(nx: usize, ny: usize) -> Self {
        Self {
            u: Field2D::new(nx, ny, T::zero()),
            v: Field2D::new(nx, ny, T::zero()),
            u_star: Field2D::new(nx, ny, T::zero()),
            v_star: Field2D::new(nx, ny, T::zero()),
            p: Field2D::new(nx, ny, T::zero()),
            p_prime: Field2D::new(nx, ny, T::zero()),
            density: Field2D::new(nx, ny, T::one()),
            viscosity: Field2D::new(nx, ny, T::from_f64(0.001).unwrap_or_else(T::one)),
            nx,
            ny,
        }
    }
    
    /// Create fields with specified fluid properties
    pub fn with_fluid(nx: usize, ny: usize, fluid: &Fluid<T>) -> Self {
        let mut fields = Self::new(nx, ny);
        fields.density.map_inplace(|d| *d = fluid.density.clone());
        // For initialization, use zero shear rate (Newtonian behavior)
        let viscosity = fluid.dynamic_viscosity(T::zero())
            .unwrap_or_else(|_| fluid.characteristic_viscosity());
        fields.viscosity.map_inplace(|v| *v = viscosity.clone());
        fields
    }

    /// Reset all fields to zero
    pub fn reset(&mut self) {
        self.u.map_inplace(|val| *val = T::zero());
        self.v.map_inplace(|val| *val = T::zero());
        self.u_star.map_inplace(|val| *val = T::zero());
        self.v_star.map_inplace(|val| *val = T::zero());
        self.p.map_inplace(|val| *val = T::zero());
        self.p_prime.map_inplace(|val| *val = T::zero());
    }
    
    /// Get velocity as Vector2 at point (i, j)
    #[inline]
    pub fn velocity_at(&self, i: usize, j: usize) -> Vector2<T> {
        Vector2::new(
            self.u.at(i, j).clone(),
            self.v.at(i, j).clone()
        )
    }
    
    /// Set velocity from Vector2 at point (i, j)
    #[inline]
    pub fn set_velocity_at(&mut self, i: usize, j: usize, vel: &Vector2<T>) {
        *self.u.at_mut(i, j) = vel.x.clone();
        *self.v.at_mut(i, j) = vel.y.clone();
    }
    
    /// Get predicted velocity as Vector2
    #[inline]
    pub fn velocity_star_at(&self, i: usize, j: usize) -> Vector2<T> {
        Vector2::new(
            self.u_star.at(i, j).clone(),
            self.v_star.at(i, j).clone()
        )
    }

    /// Get maximum velocity magnitude for stability analysis using iterators
    pub fn max_velocity_magnitude(&self) -> T {
        self.u.data()
            .iter()
            .zip(self.v.data().iter())
            .map(|(u, v)| {
                let u2 = u.clone() * u.clone();
                let v2 = v.clone() * v.clone();
                (u2 + v2).sqrt()
            })
            .fold(T::zero(), |acc, mag| if mag > acc { mag } else { acc })
    }
    
    /// Calculate kinematic viscosity field (nu = mu/rho)
    pub fn kinematic_viscosity(&self) -> Field2D<T> {
        Field2D {
            data: self.viscosity.data()
                .iter()
                .zip(self.density.data().iter())
                .map(|(mu, rho)| mu.clone() / rho.clone())
                .collect(),
            nx: self.nx,
            ny: self.ny,
        }
    }
    
    /// Calculate Reynolds number based on characteristic length and velocity
    pub fn reynolds_number(&self, characteristic_length: T, characteristic_velocity: T) -> T {
        let avg_density = self.density.data()
            .iter()
            .fold(T::zero(), |acc, d| acc + d.clone()) 
            / T::from_usize(self.density.data.len()).unwrap_or_else(T::one);
            
        let avg_viscosity = self.viscosity.data()
            .iter()
            .fold(T::zero(), |acc, v| acc + v.clone())
            / T::from_usize(self.viscosity.data.len()).unwrap_or_else(T::one);
            
        avg_density * characteristic_velocity * characteristic_length / avg_viscosity
    }
}