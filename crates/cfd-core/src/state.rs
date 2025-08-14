//! Simulation state representations.

use indexmap::IndexMap;
use nalgebra::{DVector, RealField, Vector3};
use num_traits::cast::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Trait for simulation states
pub trait SimulationState: Send + Sync {
    /// Scalar type
    type Scalar: RealField;

    /// Get the current time
    fn time(&self) -> Self::Scalar;

    /// Set the current time
    fn set_time(&mut self, time: Self::Scalar);

    /// Get the time step
    fn time_step(&self) -> Self::Scalar;

    /// Set the time step
    fn set_time_step(&mut self, dt: Self::Scalar);

    /// Get the iteration count
    fn iteration(&self) -> usize;

    /// Increment the iteration count
    fn increment_iteration(&mut self);

    /// Reset the state
    fn reset(&mut self);
}

/// Field variable types
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum FieldVariable {
    /// Pressure
    Pressure,
    /// Velocity (vector)
    Velocity,
    /// Temperature
    Temperature,
    /// Density
    Density,
    /// Custom field
    Custom(u32),
}

/// Generic field data
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum FieldData<T: RealField> {
    /// Scalar field
    Scalar(DVector<T>),
    /// Vector field
    Vector(Vec<Vector3<T>>),
}

impl<T: RealField> FieldData<T> {
    /// Create a scalar field with given size
    pub fn scalar(size: usize) -> Self {
        Self::Scalar(DVector::zeros(size))
    }

    /// Create a vector field with given size
    pub fn vector(size: usize) -> Self {
        Self::Vector(vec![Vector3::zeros(); size])
    }

    /// Get the size of the field
    pub fn len(&self) -> usize {
        match self {
            Self::Scalar(v) => v.len(),
            Self::Vector(v) => v.len(),
        }
    }

    /// Check if the field is empty
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// Simulation state implementation for field variables and time integration
#[derive(Debug, Clone)]
pub struct FieldState<T: RealField> {
    /// Current simulation time
    pub time: T,
    /// Time step
    pub time_step: T,
    /// Iteration count
    pub iteration: usize,
    /// Field variables
    pub fields: IndexMap<FieldVariable, FieldData<T>>,
}

impl<T: RealField + FromPrimitive> FieldState<T> {
    /// Create a new field state
    pub fn new() -> Self {
        Self {
            time: T::zero(),
            time_step: T::from_f64(1e-3).unwrap(),
            iteration: 0,
            fields: IndexMap::new(),
        }
    }

    /// Add a scalar field
    pub fn add_scalar_field(&mut self, var: FieldVariable, size: usize) {
        self.fields.insert(var, FieldData::scalar(size));
    }

    /// Add a vector field
    pub fn add_vector_field(&mut self, var: FieldVariable, size: usize) {
        self.fields.insert(var, FieldData::vector(size));
    }

    /// Get a scalar field
    pub fn scalar_field(&self, var: FieldVariable) -> Option<&DVector<T>> {
        match self.fields.get(&var) {
            Some(FieldData::Scalar(v)) => Some(v),
            _ => None,
        }
    }

    /// Get a mutable scalar field
    pub fn scalar_field_mut(&mut self, var: FieldVariable) -> Option<&mut DVector<T>> {
        match self.fields.get_mut(&var) {
            Some(FieldData::Scalar(v)) => Some(v),
            _ => None,
        }
    }

    /// Get a vector field
    pub fn vector_field(&self, var: FieldVariable) -> Option<&[Vector3<T>]> {
        match self.fields.get(&var) {
            Some(FieldData::Vector(v)) => Some(v),
            _ => None,
        }
    }

    /// Get a mutable vector field
    pub fn vector_field_mut(&mut self, var: FieldVariable) -> Option<&mut Vec<Vector3<T>>> {
        match self.fields.get_mut(&var) {
            Some(FieldData::Vector(v)) => Some(v),
            _ => None,
        }
    }
}

impl<T: RealField + FromPrimitive> Default for FieldState<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy> SimulationState for FieldState<T> {
    type Scalar = T;

    fn time(&self) -> Self::Scalar {
        self.time
    }

    fn set_time(&mut self, time: Self::Scalar) {
        self.time = time;
    }

    fn time_step(&self) -> Self::Scalar {
        self.time_step
    }

    fn set_time_step(&mut self, dt: Self::Scalar) {
        self.time_step = dt;
    }

    fn iteration(&self) -> usize {
        self.iteration
    }

    fn increment_iteration(&mut self) {
        self.iteration += 1;
    }

    fn reset(&mut self) {
        self.time = T::zero();
        self.iteration = 0;
        // Reset field data to zeros
        for field in self.fields.values_mut() {
            match field {
                FieldData::Scalar(v) => v.fill(T::zero()),
                FieldData::Vector(v) => v.iter_mut().for_each(|x| *x = Vector3::zeros()),
            }
        }
    }
}

/// State snapshot for checkpointing
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StateSnapshot<T: RealField> {
    /// Simulation time
    pub time: T,
    /// Iteration number
    pub iteration: usize,
    /// Serialized state data
    pub data: Vec<u8>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_state() {
        let mut state = BasicState::<f64>::new();
        state.add_scalar_field(FieldVariable::Pressure, 10);
        state.add_vector_field(FieldVariable::Velocity, 10);

        assert_eq!(state.fields.len(), 2);
        assert!(state.scalar_field(FieldVariable::Pressure).is_some());
        assert!(state.vector_field(FieldVariable::Velocity).is_some());

        state.set_time(1.0);
        state.increment_iteration();
        assert_eq!(state.time(), 1.0);
        assert_eq!(state.iteration(), 1);
    }
}