//! Flow field representations with zero-copy operations
//!
//! Provides core data structures for velocity, pressure, and scalar fields
//! following SSOT and zero-copy principles.

use nalgebra::{RealField, Vector3};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Flow field abstraction representing velocity, pressure, and scalar fields
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlowField<T: RealField + Copy> {
    /// Velocity field components
    pub velocity: VelocityField<T>,
    /// Pressure field
    pub pressure: PressureField<T>,
    /// Scalar fields (temperature, concentration, etc.)
    pub scalars: HashMap<String, ScalarField<T>>,
}

/// Velocity field representation with zero-copy operations
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VelocityField<T: RealField + Copy> {
    /// Velocity components (u, v, w)
    pub components: Vec<Vector3<T>>,
    /// Field dimensions
    pub dimensions: (usize, usize, usize),
}

/// Pressure field representation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PressureField<T: RealField + Copy> {
    /// Pressure values
    pub values: Vec<T>,
    /// Field dimensions
    pub dimensions: (usize, usize, usize),
}

/// Generic scalar field for temperature, concentration, etc.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScalarField<T: RealField + Copy> {
    /// Scalar values
    pub values: Vec<T>,
    /// Field dimensions
    pub dimensions: (usize, usize, usize),
    /// Field name/type
    pub name: String,
}

impl<T: RealField + Copy> FlowField<T> {
    /// Create a new flow field with specified dimensions
    #[must_use]
    pub fn new(nx: usize, ny: usize, nz: usize) -> Self {
        let total_points = nx * ny * nz;
        Self {
            velocity: VelocityField {
                components: vec![Vector3::zeros(); total_points],
                dimensions: (nx, ny, nz),
            },
            pressure: PressureField {
                values: vec![T::zero(); total_points],
                dimensions: (nx, ny, nz),
            },
            scalars: HashMap::new(),
        }
    }

    /// Add a scalar field
    pub fn add_scalar(&mut self, name: impl Into<String>, values: Vec<T>) {
        let name_str = name.into();
        let scalar = ScalarField {
            values,
            dimensions: self.velocity.dimensions,
            name: name_str.clone(),
        };
        self.scalars.insert(name_str, scalar);
    }

    /// Get total number of points
    #[must_use]
    pub fn total_points(&self) -> usize {
        let (nx, ny, nz) = self.velocity.dimensions;
        nx * ny * nz
    }
}

impl<T: RealField + Copy> VelocityField<T> {
    /// Get velocity at a specific grid point
    #[must_use]
    pub fn get(&self, i: usize, j: usize, k: usize) -> Option<&Vector3<T>> {
        let (nx, ny, _) = self.dimensions;
        let idx = k * nx * ny + j * nx + i;
        self.components.get(idx)
    }

    /// Get mutable velocity at a specific grid point
    pub fn get_mut(&mut self, i: usize, j: usize, k: usize) -> Option<&mut Vector3<T>> {
        let (nx, ny, _) = self.dimensions;
        let idx = k * nx * ny + j * nx + i;
        self.components.get_mut(idx)
    }

    /// Create a view of the velocity field as a slice
    #[must_use]
    pub fn as_slice(&self) -> &[Vector3<T>] {
        &self.components
    }

    /// Create a mutable view of the velocity field
    pub fn as_mut_slice(&mut self) -> &mut [Vector3<T>] {
        &mut self.components
    }
}

impl<T: RealField + Copy> PressureField<T> {
    /// Get pressure at a specific grid point
    #[must_use]
    pub fn get(&self, i: usize, j: usize, k: usize) -> Option<T> {
        let (nx, ny, _) = self.dimensions;
        let idx = k * nx * ny + j * nx + i;
        self.values.get(idx).copied()
    }

    /// Set pressure at a specific grid point
    pub fn set(&mut self, i: usize, j: usize, k: usize, value: T) {
        let (nx, ny, _) = self.dimensions;
        let idx = k * nx * ny + j * nx + i;
        if let Some(p) = self.values.get_mut(idx) {
            *p = value;
        }
    }

    /// Create a view of the pressure field as a slice
    #[must_use]
    pub fn as_slice(&self) -> &[T] {
        &self.values
    }
}
