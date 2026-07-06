//! Time step controllers.

use crate::error::{Error, Result};
use eunomia::{FloatElement, NumericElement, RealField};

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

impl<T: RealField + FloatElement + Copy> Default for AdaptiveTimeStepController<T> {
    fn default() -> Self {
        Self {
            target_error: <T as FloatElement>::from_f64(1e-6),
            safety_factor: <T as FloatElement>::from_f64(0.9),
            max_increase: <T as FloatElement>::from_f64(2.0),
            min_decrease: <T as FloatElement>::from_f64(0.1),
        }
    }
}

impl<T: RealField + FloatElement + Copy> AdaptiveTimeStepController<T> {
    /// Calculate new time step based on error estimate
    pub fn calculate_dt(&self, current_dt: T, error: T, order: usize) -> Result<T> {
        let exponent_denominator = integration_order_denominator::<T>(order, 1)?;
        let factor = <T as FloatElement>::powf(
            self.target_error / error,
            <T as NumericElement>::ONE / exponent_denominator,
        );
        let factor = factor * self.safety_factor;
        let factor = clamp(factor, self.min_decrease, self.max_increase);
        Ok(current_dt * factor)
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

impl<T: RealField + FloatElement + Copy> Default for VariableTimeStep<T> {
    fn default() -> Self {
        Self {
            dt_min: <T as FloatElement>::from_f64(1e-10),
            dt_max: <T as FloatElement>::from_f64(0.1),
            safety_factor: <T as FloatElement>::from_f64(0.9),
            target_error: <T as FloatElement>::from_f64(1e-6),
        }
    }
}

impl<T: RealField + Copy + FloatElement> VariableTimeStep<T> {
    /// Calculate new time step based on error estimate
    pub fn calculate_dt(&self, current_dt: T, error: T, order: usize) -> Result<T> {
        if error < <T as RealField>::EPSILON {
            return Ok(min_value(
                self.dt_max,
                current_dt * <T as FloatElement>::from_f64(2.0),
            ));
        }

        let exponent = <T as NumericElement>::ONE / integration_order_denominator::<T>(order, 0)?;
        let factor =
            self.safety_factor * <T as FloatElement>::powf(self.target_error / error, exponent);

        let current_dt = current_dt * factor;
        Ok(clamp(current_dt, self.dt_min, self.dt_max))
    }
}

fn integration_order_denominator<T: FloatElement>(order: usize, offset: usize) -> Result<T> {
    if order == 0 {
        return Err(Error::InvalidConfiguration(
            "time-integration order must be positive".to_string(),
        ));
    }

    let denominator = order.checked_add(offset).ok_or_else(|| {
        Error::InvalidConfiguration("time-integration order denominator overflowed".to_string())
    })?;
    let denominator = u32::try_from(denominator).map_err(|_| {
        Error::InvalidConfiguration(
            "time-integration order exceeds the supported u32 range".to_string(),
        )
    })?;

    Ok(<T as FloatElement>::from_f64(f64::from(denominator)))
}

fn min_value<T: PartialOrd + Copy>(left: T, right: T) -> T {
    if left <= right {
        left
    } else {
        right
    }
}

fn clamp<T: PartialOrd + Copy>(value: T, min: T, max: T) -> T {
    if value < min {
        min
    } else if value > max {
        max
    } else {
        value
    }
}

#[cfg(test)]
mod tests {
    use super::{AdaptiveTimeStepController, VariableTimeStep};
    use crate::compute::time::integrators::ForwardEuler;
    use crate::compute::time::integrators::TimeIntegrator;
    use approx::assert_abs_diff_eq;
    use leto::Array1;

    fn state_from(values: Vec<f64>) -> Array1<f64> {
        Array1::from_vec([values.len()], values).expect("test state shape matches supplied values")
    }

    #[test]
    fn test_forward_euler() {
        let integrator = ForwardEuler;

        // Test with a simple ODE: dy/dt = -y
        // Solution: y(t) = y0 * exp(-t)
        let mut state = Array1::from_elem([1], 1.0);
        let dt = 0.1;

        // Define derivative function
        let derivative = |_t: f64, y: &Array1<f64>| -> Array1<f64> { state_from(vec![-y[[0]]]) };

        // Take one step
        integrator
            .step(&mut state, 0.0, dt, derivative)
            .expect("CRITICAL: Add proper error handling");

        // After one step: y ≈ y0 * (1 - dt) = 1.0 * (1 - 0.1) = 0.9
        assert_abs_diff_eq!(state[[0]], 0.9, epsilon = 1e-10);
    }

    #[test]
    fn adaptive_controller_defaults_are_value_semantic() {
        let controller = AdaptiveTimeStepController::<f64>::default();

        assert_eq!(controller.target_error, 1e-6);
        assert_eq!(controller.safety_factor, 0.9);
        assert_eq!(controller.max_increase, 2.0);
        assert_eq!(controller.min_decrease, 0.1);
    }

    #[test]
    fn adaptive_controller_rejects_zero_order() {
        let controller = AdaptiveTimeStepController::<f64>::default();

        let err = controller.calculate_dt(0.1, 1e-5, 0).unwrap_err();

        assert_eq!(
            err.to_string(),
            "Invalid configuration: time-integration order must be positive"
        );
    }

    #[test]
    fn variable_controller_clamps_to_bounds() {
        let controller = VariableTimeStep::<f64>::default();

        let doubled = controller
            .calculate_dt(0.08, 0.0, 2)
            .expect("zero error should use bounded doubling");
        let reduced = controller
            .calculate_dt(1e-12, 1e-3, 2)
            .expect("finite error should compute bounded reduction");

        assert_eq!(doubled, 0.1);
        assert_eq!(reduced, 1e-10);
    }
}
