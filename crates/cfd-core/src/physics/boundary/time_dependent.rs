//! Time-dependent boundary condition functionality

use crate::physics::boundary::BoundaryCondition;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

// Mathematical constants
const TWO_PI: f64 = 2.0 * std::f64::consts::PI;

/// Time-dependent boundary condition specification
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TimeDependentSpec<T: RealField + Copy> {
    /// Time function type
    pub function_type: TimeFunctionType,
    /// Function parameters
    pub parameters: Vec<T>,
}

/// Time function types
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum TimeFunctionType {
    /// Constant value
    Constant,
    /// Linear ramp: value = initial + rate * time
    Linear,
    /// Sinusoidal variation: value = amplitude * sin(frequency * time + phase) + offset
    Sinusoidal,
    /// Exponential decay/growth: value = initial * exp(rate * time)
    Exponential,
    /// Polynomial: value = sum(coefficients[i] * time^i)
    Polynomial,
    /// Custom function (user-defined)
    Custom(String),
}

impl<T: RealField + Copy + FromPrimitive> TimeDependentSpec<T> {
    /// Create a constant time function
    #[must_use]
    pub fn constant() -> Self {
        Self {
            function_type: TimeFunctionType::Constant,
            parameters: vec![],
        }
    }

    /// Create a linear ramp function
    pub fn linear(initial: T, rate: T) -> Self {
        Self {
            function_type: TimeFunctionType::Linear,
            parameters: vec![initial, rate],
        }
    }

    /// Create a sinusoidal function
    pub fn sinusoidal(amplitude: T, frequency: T, phase: T, offset: T) -> Self {
        Self {
            function_type: TimeFunctionType::Sinusoidal,
            parameters: vec![amplitude, frequency, phase, offset],
        }
    }

    /// Create an exponential function
    pub fn exponential(initial: T, rate: T) -> Self {
        Self {
            function_type: TimeFunctionType::Exponential,
            parameters: vec![initial, rate],
        }
    }

    /// Create a polynomial function
    #[must_use]
    pub fn polynomial(coefficients: Vec<T>) -> Self {
        Self {
            function_type: TimeFunctionType::Polynomial,
            parameters: coefficients,
        }
    }

    /// Evaluate the time function at a given time
    pub fn evaluate(&self, time: T) -> T {
        match self.function_type {
            TimeFunctionType::Constant | TimeFunctionType::Custom(_) => T::one(),

            TimeFunctionType::Linear => {
                if self.parameters.len() >= 2 {
                    self.parameters[0] + self.parameters[1] * time
                } else {
                    T::one()
                }
            }

            TimeFunctionType::Sinusoidal => {
                if self.parameters.len() >= 4 {
                    let two_pi = T::from_f64(TWO_PI).unwrap_or_else(T::one);
                    let amplitude = self.parameters[0];
                    let frequency = self.parameters[1];
                    let phase = self.parameters[2];
                    let offset = self.parameters[3];
                    amplitude * (two_pi * frequency * time + phase).sin() + offset
                } else {
                    T::one()
                }
            }

            TimeFunctionType::Exponential => {
                if self.parameters.len() >= 2 {
                    self.parameters[0] * (self.parameters[1] * time).exp()
                } else {
                    T::one()
                }
            }

            TimeFunctionType::Polynomial => {
                let mut result = T::zero();
                let mut time_power = T::one();
                for coeff in &self.parameters {
                    result += *coeff * time_power;
                    time_power *= time;
                }
                result
            }
        }
    }

    /// Apply time function to a boundary condition
    pub fn apply_to_condition(
        &self,
        condition: &BoundaryCondition<T>,
        time: T,
    ) -> BoundaryCondition<T> {
        let factor = self.evaluate(time);

        match condition {
            BoundaryCondition::Dirichlet {
                value,
                component_values,
            } => BoundaryCondition::Dirichlet {
                value: *value * factor,
                component_values: component_values.as_ref().map(|comps| {
                    comps.iter().map(|opt| opt.map(|v| v * factor)).collect()
                }),
            },
            BoundaryCondition::Neumann { gradient } => BoundaryCondition::Neumann {
                gradient: *gradient * factor,
            },
            BoundaryCondition::Robin { alpha, beta, gamma } => BoundaryCondition::Robin {
                alpha: *alpha * factor,
                beta: *beta * factor,
                gamma: *gamma * factor,
            },
            other => other.clone(),
        }
    }
}
