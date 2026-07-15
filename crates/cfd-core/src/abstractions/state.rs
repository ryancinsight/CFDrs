//! Simulation state representations.

use eunomia::{FloatElement, RealField};
use indexmap::IndexMap;
use leto::{geometry::Vector3, Array1};
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
pub enum FieldData<T: RealField + Copy> {
    /// Scalar field
    Scalar(Array1<T>),
    /// Vector field
    Vector(Vec<Vector3<T>>),
}

impl<T: RealField + Copy> FieldData<T> {
    /// Create a scalar field with given size
    #[must_use]
    pub fn scalar(size: usize) -> Self {
        Self::Scalar(Array1::from_elem([size], T::ZERO))
    }

    /// Create a vector field with given size
    #[must_use]
    pub fn vector(size: usize) -> Self {
        Self::Vector(vec![Vector3::zeros(); size])
    }

    /// Get the size of the field
    #[must_use]
    pub fn len(&self) -> usize {
        match self {
            Self::Scalar(v) => v.size(),
            Self::Vector(v) => v.len(),
        }
    }

    /// Check if the field is empty
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// Simulation state implementation for field variables and time integration
#[derive(Debug, Clone)]
pub struct FieldState<T: RealField + Copy> {
    /// Current simulation time
    pub time: T,
    /// Time step
    pub time_step: T,
    /// Iteration count
    pub iteration: usize,
    /// Field variables
    pub fields: IndexMap<FieldVariable, FieldData<T>>,
}

impl<T: RealField + FloatElement + Copy> FieldState<T> {
    /// Create a new field state
    #[must_use]
    pub fn new() -> Self {
        Self {
            time: T::ZERO,
            time_step: <T as FloatElement>::from_f64(1e-3),
            iteration: 0,
            fields: IndexMap::new(),
        }
    }
}

impl<T: RealField + Copy> FieldState<T> {
    /// Add a scalar field
    pub fn add_scalar_field(&mut self, var: FieldVariable, size: usize) {
        self.fields.insert(var, FieldData::scalar(size));
    }

    /// Add a vector field
    pub fn add_vector_field(&mut self, var: FieldVariable, size: usize) {
        self.fields.insert(var, FieldData::vector(size));
    }

    /// Get a scalar field
    pub fn scalar_field(&self, var: FieldVariable) -> Option<&Array1<T>> {
        match self.fields.get(&var) {
            Some(FieldData::Scalar(v)) => Some(v),
            _ => None,
        }
    }

    /// Get a mutable scalar field
    pub fn scalar_field_mut(&mut self, var: FieldVariable) -> Option<&mut Array1<T>> {
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

impl<T: RealField + FloatElement + Copy> Default for FieldState<T> {
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
        self.time = T::ZERO;
        self.iteration = 0;
        // Reset field data to zeros
        for field in self.fields.values_mut() {
            match field {
                FieldData::Scalar(v) => v.fill(T::ZERO),
                FieldData::Vector(v) => v.iter_mut().for_each(|x| *x = Vector3::zeros()),
            }
        }
    }
}

/// State snapshot for checkpointing
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StateSnapshot<T: RealField + Copy> {
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
    use leto::Storage;

    #[test]
    fn test_core_state() {
        let mut state = FieldState::<f64>::new();
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

    #[test]
    fn field_state_reset_clears_leto_fields() {
        let mut state = FieldState::<f64>::new();
        state.add_scalar_field(FieldVariable::Pressure, 3);
        state.add_vector_field(FieldVariable::Velocity, 2);

        state.scalar_field_mut(FieldVariable::Pressure).unwrap()[1] = 5.0;
        state.vector_field_mut(FieldVariable::Velocity).unwrap()[0] = Vector3::new(1.0, 2.0, 3.0);
        state.set_time(2.0);
        state.increment_iteration();

        state.reset();

        assert_eq!(state.time(), 0.0);
        assert_eq!(state.iteration(), 0);
        assert_eq!(
            state
                .scalar_field(FieldVariable::Pressure)
                .unwrap()
                .storage()
                .as_slice(),
            &[0.0, 0.0, 0.0]
        );
        assert_eq!(
            state.vector_field(FieldVariable::Velocity).unwrap(),
            &[Vector3::zeros(), Vector3::zeros()]
        );
    }

    #[test]
    fn field_data_scalar_round_trips_as_leto_array() {
        let mut state = FieldState::<f64>::new();
        state.add_scalar_field(FieldVariable::Pressure, 2);
        state.scalar_field_mut(FieldVariable::Pressure).unwrap()[0] = 4.0;

        let encoded = serde_json::to_string(&state.fields[&FieldVariable::Pressure]).unwrap();
        let decoded: FieldData<f64> = serde_json::from_str(&encoded).unwrap();

        match decoded {
            FieldData::Scalar(values) => {
                assert_eq!(values.shape(), [2]);
                assert_eq!(values.storage().as_slice(), &[4.0, 0.0]);
            }
            FieldData::Vector(_) => panic!("expected scalar field data"),
        }
    }
}
