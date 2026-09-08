//! Time-dependent boundary condition functionality

use crate::physics::boundary::BoundaryCondition;
use eunomia::FloatElement;
use eunomia::RealField;
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
    /// Polynomial: value = sum(coefficients\[i\] * time^i)
    Polynomial,
    /// Custom function (user-defined)
    Custom(String),
}

impl<T: RealField + FloatElement + Copy> TimeDependentSpec<T> {
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
            TimeFunctionType::Constant | TimeFunctionType::Custom(_) => T::ONE,

            TimeFunctionType::Linear => {
                if self.parameters.len() >= 2 {
                    self.parameters[0] + self.parameters[1] * time
                } else {
                    T::ONE
                }
            }

            TimeFunctionType::Sinusoidal => {
                if self.parameters.len() >= 4 {
                    let two_pi = <T as FloatElement>::from_f64(TWO_PI);
                    let amplitude = self.parameters[0];
                    let frequency = self.parameters[1];
                    let phase = self.parameters[2];
                    let offset = self.parameters[3];
                    amplitude * FloatElement::sin(two_pi * frequency * time + phase) + offset
                } else {
                    T::ONE
                }
            }

            TimeFunctionType::Exponential => {
                if self.parameters.len() >= 2 {
                    self.parameters[0] * FloatElement::exp(self.parameters[1] * time)
                } else {
                    T::ONE
                }
            }

            TimeFunctionType::Polynomial => {
                let mut result = T::ZERO;
                let mut time_power = T::ONE;
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
                component_values: component_values
                    .as_ref()
                    .map(|comps| comps.iter().map(|opt| opt.map(|v| v * factor)).collect()),
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

#[cfg(test)]
mod tests {
    use super::TimeDependentSpec;
    use crate::physics::boundary::BoundaryCondition;

    #[test]
    fn sinusoidal_time_function_matches_closed_form() {
        let spec = TimeDependentSpec::sinusoidal(2.0_f64, 0.25, 0.5, 1.0);

        let value = spec.evaluate(0.5);
        let expected = 2.0 * f64::sin(2.0 * std::f64::consts::PI * 0.25 * 0.5 + 0.5) + 1.0;

        assert!((value - expected).abs() <= 8.0 * f64::EPSILON * expected.abs().max(1.0));
    }

    #[test]
    fn exponential_time_function_matches_closed_form() {
        let spec = TimeDependentSpec::exponential(3.0_f64, -0.25);

        let value = spec.evaluate(4.0);
        let expected = 3.0 * f64::exp(-0.25 * 4.0);

        assert!((value - expected).abs() <= 8.0 * f64::EPSILON * expected);
    }

    #[test]
    fn time_function_scales_boundary_condition_components() {
        let spec = TimeDependentSpec::linear(2.0_f64, 0.5);
        let condition = BoundaryCondition::Dirichlet {
            value: 4.0,
            component_values: Some(vec![Some(1.0), None, Some(3.0)]),
        };

        let scaled = spec.apply_to_condition(&condition, 2.0);

        assert_eq!(
            scaled,
            BoundaryCondition::Dirichlet {
                value: 12.0,
                component_values: Some(vec![Some(3.0), None, Some(9.0)]),
            }
        );
    }
}
