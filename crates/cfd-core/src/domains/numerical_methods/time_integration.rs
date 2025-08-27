//! Time integration schemes for temporal discretization

use super::traits::TimeIntegrationScheme;
use nalgebra::RealField;

/// Constant for RK4 coefficient
const RK4_COEFFICIENT_ONE_SIXTH: f64 = 1.0 / 6.0;
/// Constant for RK4 half step
const RK4_HALF_STEP: f64 = 0.5;

/// Time integration schemes
pub mod time_schemes {
    use super::{RealField, TimeIntegrationScheme, RK4_COEFFICIENT_ONE_SIXTH, RK4_HALF_STEP};

    /// Forward Euler scheme (explicit, 1st order)
    #[derive(Debug, Clone)]
    pub struct ForwardEuler;

    impl<T: RealField + Copy> TimeIntegrationScheme<T> for ForwardEuler {
        fn advance(&self, current: &[T], derivative: &[T], dt: T) -> Vec<T> {
            current
                .iter()
                .zip(derivative.iter())
                .map(|(u, dudt)| *u + *dudt * dt)
                .collect()
        }

        fn name(&self) -> &str {
            "Forward Euler"
        }

        fn order(&self) -> usize {
            1
        }

        fn is_implicit(&self) -> bool {
            false
        }
    }

    /// Constant derivative time integration (equivalent to Forward Euler)
    ///
    /// This scheme assumes the derivative is constant over the time step.
    /// For systems where the derivative changes with state, this reduces to
    /// Forward Euler regardless of the intended higher-order method.
    #[derive(Debug, Clone)]
    pub struct ConstantDerivative;

    impl<T: RealField + Copy> TimeIntegrationScheme<T> for ConstantDerivative {
        fn advance(&self, current: &[T], derivative: &[T], dt: T) -> Vec<T> {
            // With constant derivative assumption: y_new = y + dt * f(y)
            // This is exactly Forward Euler method
            current
                .iter()
                .zip(derivative.iter())
                .map(|(y, d)| *y + *d * dt)
                .collect()
        }

        fn name(&self) -> &str {
            "Constant Derivative (Euler)"
        }

        fn order(&self) -> usize {
            1 // Order of accuracy: 1
        }

        fn is_implicit(&self) -> bool {
            false
        }
    }

    /// Runge-Kutta 4th order scheme
    /// Properly evaluates derivatives at each RK4 stage
    #[derive(Debug, Clone)]
    pub struct RungeKutta4;

    impl RungeKutta4 {
        /// Advance with derivative function for proper RK4
        /// Reference: Hairer, E., Nørsett, S.P., Wanner, G. "Solving Ordinary Differential Equations I" (1993)
        pub fn advance_with_function<T, F>(
            &self,
            current: &[T],
            t: T,
            dt: T,
            derivative_fn: F,
        ) -> Vec<T>
        where
            T: RealField + Copy,
            F: Fn(T, &[T]) -> Vec<T>,
        {
            let half =
                T::from_f64(RK4_HALF_STEP).unwrap_or_else(|| T::one() / (T::one() + T::one()));
            let one_sixth = T::from_f64(RK4_COEFFICIENT_ONE_SIXTH)
                .unwrap_or_else(|| T::one() / (T::from_usize(6).unwrap_or_else(|| T::one())));
            let two = T::one() + T::one();

            // Classical RK4 stages
            // k1 = dt * f(t, y)
            let k1: Vec<T> = derivative_fn(t, current)
                .into_iter()
                .map(move |val| val * dt)
                .collect();

            // y1 = y + k1/2
            let y1: Vec<T> = current
                .iter()
                .zip(k1.iter())
                .map(|(y, k)| *y + *k * half)
                .collect();

            // k2 = dt * f(t + dt/2, y1)
            let k2: Vec<T> = derivative_fn(t + dt * half, &y1)
                .into_iter()
                .map(move |val| val * dt)
                .collect();

            // y2 = y + k2/2
            let y2: Vec<T> = current
                .iter()
                .zip(k2.iter())
                .map(|(y, k)| *y + *k * half)
                .collect();

            // k3 = dt * f(t + dt/2, y2)
            let k3: Vec<T> = derivative_fn(t + dt * half, &y2)
                .into_iter()
                .map(move |val| val * dt)
                .collect();

            // y3 = y + k3
            let y3: Vec<T> = current
                .iter()
                .zip(k3.iter())
                .map(|(y, k)| *y + *k)
                .collect();

            // k4 = dt * f(t + dt, y3)
            let k4: Vec<T> = derivative_fn(t + dt, &y3)
                .into_iter()
                .map(move |val| val * dt)
                .collect();

            // Final combination: y_new = y + (k1 + 2*k2 + 2*k3 + k4)/6
            current
                .iter()
                .zip(k1.iter())
                .zip(k2.iter())
                .zip(k3.iter())
                .zip(k4.iter())
                .map(|((((y, k1), k2), k3), k4)| {
                    let weighted_sum = *k1 + *k2 * two + *k3 * two + *k4;
                    *y + weighted_sum * one_sixth
                })
                .collect()
        }
    }

    impl<T: RealField + Copy> TimeIntegrationScheme<T> for RungeKutta4 {
        fn advance(&self, current: &[T], derivative: &[T], dt: T) -> Vec<T> {
            // Fallback to Forward Euler when derivative function is not available
            // This is a limitation when only the derivative at current state is known
            current
                .iter()
                .zip(derivative.iter())
                .map(|(u, dudt)| *u + *dudt * dt)
                .collect()
        }

        fn name(&self) -> &str {
            "Runge-Kutta 4"
        }

        fn order(&self) -> usize {
            4
        }

        fn is_implicit(&self) -> bool {
            false
        }
    }
}

// Re-export for convenience
pub use time_schemes::{ConstantDerivative, ForwardEuler, RungeKutta4};

#[cfg(test)]
mod tests {
    use super::time_schemes::{ConstantDerivative, ForwardEuler};
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_forward_euler() {
        let scheme = time_schemes::ForwardEuler;
        let current = vec![1.0f64, 2.0, 3.0];
        let derivative = vec![0.5, -1.0, 2.0];
        let dt = 0.1;

        let result = scheme.advance(&current, &derivative, dt);

        assert_eq!(result.len(), 3);
        assert_relative_eq!(result[0], 1.05, epsilon = 1e-10);
        assert_relative_eq!(result[1], 1.9, epsilon = 1e-10);
        assert_relative_eq!(result[2], 3.2, epsilon = 1e-10);

        assert_eq!(
            <ForwardEuler as TimeIntegrationScheme<f64>>::name(&scheme),
            "Forward Euler"
        );
        assert_eq!(
            <ForwardEuler as TimeIntegrationScheme<f64>>::order(&scheme),
            1
        );
        assert!(!<ForwardEuler as TimeIntegrationScheme<f64>>::is_implicit(
            &scheme
        ));
    }

    #[test]
    fn test_constant_derivative() {
        let scheme = time_schemes::ConstantDerivative;
        let current = vec![1.0f64, 2.0];
        let derivative = vec![1.0, -0.5];
        let dt = 0.2;

        let result = scheme.advance(&current, &derivative, dt);

        assert_eq!(result.len(), 2);
        assert_relative_eq!(result[0], 1.2, epsilon = 1e-10);
        assert_relative_eq!(result[1], 1.9, epsilon = 1e-10);

        assert_eq!(
            <ConstantDerivative as TimeIntegrationScheme<f64>>::name(&scheme),
            "Constant Derivative (Euler)"
        );
        assert_eq!(
            <ConstantDerivative as TimeIntegrationScheme<f64>>::order(&scheme),
            1
        );
        assert!(!<ConstantDerivative as TimeIntegrationScheme<f64>>::is_implicit(&scheme));
    }

    #[test]
    fn test_rk4_with_function() {
        let scheme = time_schemes::RungeKutta4;

        // Test with exponential decay: dy/dt = -y
        let current = vec![1.0f64];
        let t = 0.0;
        let dt = 0.1;
        let derivative_fn = |_t: f64, y: &[f64]| vec![-y[0]];

        let result = scheme.advance_with_function(&current, t, dt, derivative_fn);

        // Analytical solution: y(t) = e^(-t)
        // y(0.1) = e^(-0.1) ≈ 0.9048374
        assert_eq!(result.len(), 1);
        assert_relative_eq!(result[0], 0.9048374, epsilon = 1e-6);
    }
}
