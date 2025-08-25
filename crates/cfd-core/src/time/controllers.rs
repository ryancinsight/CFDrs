//! Time step controllers.

use nalgebra::RealField;
use num_traits::cast::FromPrimitive;
use num_traits::Float;

/// Adaptive time stepping controller for error-based time step adjustment
pub struct AdaptiveTimeStepController<T: RealField + Copy> {
    /// Target error tolerance
    pub target_error: T,
    /// Safety factor for time step adjustment
    pub safety_factor: T,
    /// Maximum time step increase factor
    pub max_increase: T,
    /// Minimum time step decrease factor
    pub min_decrease: T,
}

impl<T: RealField + FromPrimitive + Copy> Default for AdaptiveTimeStepController<T> {
    fn default() -> Self {
        Self {
            target_error: T::from_f64(1e-6).unwrap_or_else(T::zero),
            safety_factor: T::from_f64(0.9).unwrap_or_else(T::zero),
            max_increase: T::from_f64(2.0).unwrap_or_else(T::zero),
            min_decrease: T::from_f64(0.1).unwrap_or_else(T::zero),
        }
    }
}

impl<T: RealField + FromPrimitive + Copy> AdaptiveTimeStepController<T> {
    /// Calculate new time step based on error estimate
    pub fn calculate_dt(&self, current_dt: T, error: T, order: usize) -> T {
        let factor = (self.target_error / error)
            .powf(T::one() / T::from_usize(order + 1).unwrap_or_else(T::one));
        let factor = factor * self.safety_factor;
        let factor = factor.min(self.max_increase).max(self.min_decrease);
        current_dt * factor
    }
}

/// Variable time step controller with min/max constraints
pub struct VariableTimeStep<T: RealField + Copy> {
    /// Minimum allowed time step
    pub dt_min: T,
    /// Maximum allowed time step
    pub dt_max: T,
    /// Safety factor for step size adjustment
    pub safety_factor: T,
    /// Target error tolerance
    pub target_error: T,
}

impl<T: RealField + FromPrimitive + Copy> Default for VariableTimeStep<T> {
    fn default() -> Self {
        Self {
            dt_min: T::from_f64(1e-10).unwrap_or_else(|| T::one()),
            dt_max: T::from_f64(0.1).unwrap_or_else(|| T::one()),
            safety_factor: T::from_f64(0.9).unwrap_or_else(|| T::one()),
            target_error: T::from_f64(1e-6).unwrap_or_else(|| T::one()),
        }
    }
}

impl<T: RealField + Copy + FromPrimitive + Float> VariableTimeStep<T> {
    /// Calculate new time step based on error estimate
    pub fn calculate_dt(&self, current_dt: T, error: T, order: usize) -> T {
        if error < T::epsilon() {
            return num_traits::Float::min(
                self.dt_max,
                current_dt * T::from_f64(2.0).unwrap_or_else(|| T::one()),
            );
        }

        let exponent = T::one() / T::from_usize(order).unwrap_or_else(|| T::one());
        let factor =
            self.safety_factor * num_traits::Float::powf(self.target_error / error, exponent);

        let current_dt = current_dt * factor;
        let max_dt = num_traits::Float::max(current_dt, self.dt_min);
        num_traits::Float::min(max_dt, self.dt_max)
    }
}

#[cfg(test)]
mod tests {
    use crate::time::integrators::ForwardEuler;
    use crate::time::integrators::TimeIntegrator;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_forward_euler() {
        let integrator = ForwardEuler;

        // Test with a simple ODE: dy/dt = -y
        // Solution: y(t) = y0 * exp(-t)
        let mut state = nalgebra::DVector::from_element(1, 1.0);
        let dt = 0.1;

        // Define derivative function
        let derivative = |_t: f64, y: &nalgebra::DVector<f64>| -> nalgebra::DVector<f64> { -y };

        // Take one step
        integrator
            .step(&mut state, 0.0, dt, derivative)
            .expect("CRITICAL: Add proper error handling");

        // After one step: y ≈ y0 * (1 - dt) = 1.0 * (1 - 0.1) = 0.9
        assert_abs_diff_eq!(state[0], 0.9, epsilon = 1e-10);
    }
}
