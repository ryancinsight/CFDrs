//! Numerical differentiation methods for CFD applications.
//!
//! This module provides finite difference schemes optimized for CFD simulations
//! with support for various boundary conditions and grid types.

use nalgebra::{RealField, Vector3, DVector};
use cfd_core::{constants};
use cfd_core::error::{Error, Result};
use num_traits::FromPrimitive;

/// Finite difference schemes
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FiniteDifferenceScheme {
    /// Forward difference: f'(x) ≈ (f(x+h) - f(x)) / h
    Forward,
    /// Backward difference: f'(x) ≈ (f(x) - f(x-h)) / h
    Backward,
    /// Central difference: f'(x) ≈ (f(x+h) - f(x-h)) / (2h)
    Central,
    /// Second-order forward: f'(x) ≈ (-3f(x) + 4f(x+h) - f(x+2h)) / (2h)
    ForwardSecondOrder,
    /// Second-order backward: f'(x) ≈ (f(x-2h) - 4f(x-h) + 3f(x)) / (2h)
    BackwardSecondOrder,
}

/// Finite difference operator for 1D problems
pub struct FiniteDifference<T: RealField + Copy> {
    scheme: FiniteDifferenceScheme,
    spacing: T,
}

impl<T: RealField + From<f64> + FromPrimitive + Copy> FiniteDifference<T> {
    /// Create a new finite difference operator
    pub fn new(scheme: FiniteDifferenceScheme, spacing: T) -> Self {
        Self { scheme, spacing }
    }

    /// Create central difference operator
    pub fn central(spacing: T) -> Self {
        Self::new(FiniteDifferenceScheme::Central, spacing)
    }

    /// Create forward difference operator
    pub fn forward(spacing: T) -> Self {
        Self::new(FiniteDifferenceScheme::Forward, spacing)
    }

    /// Create backward difference operator
    pub fn backward(spacing: T) -> Self {
        Self::new(FiniteDifferenceScheme::Backward, spacing)
    }

    /// Compute first derivative using iterator combinators and zero-copy operations
    pub fn first_derivative(&self, values: &[T]) -> Result<DVector<T>> {
        if values.len() < 2 {
            return Err(Error::InvalidConfiguration(
                "Need at least 2 points for differentiation".to_string()
            ));
        }

        let n = values.len();
        let mut result = DVector::zeros(n);
        let inv_spacing = T::one() / self.spacing;

        match self.scheme {
            FiniteDifferenceScheme::Forward => {
                // Use windows() for efficient forward differences
                values.windows(2)
                    .enumerate()
                    .for_each(|(i, window)| {
                        result[i] = (window[1] - window[0]) * inv_spacing;
                    });

                // Backward difference for last point
                if n > 1 {
                    result[n-1] = (values[n-1] - values[n-2]) * inv_spacing;
                }
            },
            FiniteDifferenceScheme::Backward => {
                // Forward difference for first point
                result[0] = (values[1] - values[0]) * inv_spacing;

                // Use windows() for backward differences: (values[i] - values[i-1])
                values.windows(2)
                    .enumerate()
                    .for_each(|(i, window)| {
                        result[i + 1] = (window[1] - window[0]) * inv_spacing;
                    });
            },
            FiniteDifferenceScheme::Central => {
                // Forward difference for first point
                result[0] = (values[1] - values[0]) * inv_spacing;

                // Central difference using windows(3) for interior points
                let two_inv_spacing = inv_spacing / T::from_f64(constants::TWO).unwrap_or_else(|| T::zero());
                values.windows(3)
                    .enumerate()
                    .for_each(|(i, window)| {
                        result[i + 1] = (window[2] - window[0]) * two_inv_spacing;
                    });

                // Backward difference for last point
                if n > 1 {
                    result[n-1] = (values[n-1] - values[n-2]) * inv_spacing;
                }
            },
            FiniteDifferenceScheme::ForwardSecondOrder => {
                if n < 3 {
                    return Err(Error::InvalidConfiguration(
                        "Need at least 3 points for second-order forward difference".to_string()
                    ));
                }

                let two = T::from_f64(constants::TWO).unwrap_or_else(|| T::zero());
                let three = T::from_f64(constants::THREE).unwrap_or_else(|| T::zero());
                let four = T::from_f64(constants::FOUR).unwrap_or_else(|| T::zero());

                // Use forward difference for first n-2 points
                result.iter_mut()
                    .take(n.saturating_sub(2))
                    .enumerate()
                    .for_each(|(i, r)| {
                        *r = (-three * values[i] +
                              four * values[i+1] -
                              values[i+2]) / (two * self.spacing);
                    });

                // Use central difference for remaining points
                result.iter_mut()
                    .skip(n.saturating_sub(2))
                    .enumerate()
                    .for_each(|(idx, r)| {
                        let i = idx + n.saturating_sub(2);
                        if i > 0 && i < n-1 {
                            *r = (values[i+1] - values[i-1]) /
                                 (two * self.spacing);
                        }
                    });
            },
            FiniteDifferenceScheme::BackwardSecondOrder => {
                if n < 3 {
                    return Err(Error::InvalidConfiguration(
                        "Need at least 3 points for second-order backward difference".to_string()
                    ));
                }

                let two = T::from_f64(constants::TWO).unwrap_or_else(|| T::zero());
                let three = T::from_f64(constants::THREE).unwrap_or_else(|| T::zero());
                let four = T::from_f64(constants::FOUR).unwrap_or_else(|| T::zero());

                // Use central difference for first points
                result.iter_mut()
                    .take(2)
                    .enumerate()
                    .for_each(|(i, r)| {
                        if i > 0 && i < n-1 {
                            *r = (values[i+1] - values[i-1]) /
                                 (two * self.spacing);
                        }
                    });

                // Use backward difference for remaining points
                result.iter_mut()
                    .skip(2)
                    .enumerate()
                    .for_each(|(idx, r)| {
                        let i = idx + 2;
                        *r = (values[i-2] -
                              four * values[i-1] +
                              three * values[i]) / (two * self.spacing);
                    });
            },
        }

        Ok(result)
    }

    /// Compute second derivative using central differences
    pub fn second_derivative(&self, values: &[T]) -> Result<DVector<T>> {
        if values.len() < 3 {
            return Err(Error::InvalidConfiguration(
                "Need at least 3 points for second derivative".to_string()
            ));
        }

        let n = values.len();
        let mut result = DVector::zeros(n);
        let h_squared = self.spacing * self.spacing;

        // Use forward difference for first point
        result[0] = (values[2] - T::from_f64(constants::TWO).unwrap_or_else(|| T::zero()) * values[1] + values[0]) / h_squared;

        // Central difference for interior points
        for i in 1..n-1 {
            result[i] = (values[i+1] - T::from_f64(constants::TWO).unwrap_or_else(|| T::zero()) * values[i] + values[i-1]) / h_squared;
        }

        // Use backward difference for last point
        result[n-1] = (values[n-1] - T::from_f64(constants::TWO).unwrap_or_else(|| T::zero()) * values[n-2] + values[n-3]) / h_squared;

        Ok(result)
    }

    /// Get the finite difference scheme
    pub fn scheme(&self) -> FiniteDifferenceScheme {
        self.scheme
    }

    /// Get the grid spacing
    pub fn spacing(&self) -> T {
        self.spacing
    }
}

impl<T: RealField + From<f64> + FromPrimitive + Copy> Default for FiniteDifference<T> {
    fn default() -> Self {
        Self::central(T::one())
    }
}

/// Gradient computation for multi-dimensional fields
pub struct Gradient<T: RealField + Copy> {
    dx: T,
    dy: T,
    dz: T,
}

impl<T: RealField + From<f64> + FromPrimitive + Copy> Gradient<T> {
    /// Create new gradient operator
    pub fn new(dx: T, dy: T, dz: T) -> Self {
        Self { dx, dy, dz }
    }

    /// Create uniform spacing gradient operator
    pub fn uniform(spacing: T) -> Self {
        Self::new(spacing, spacing, spacing)
    }

    /// Compute gradient of a 1D field
    pub fn gradient_1d(&self, field: &[T]) -> Result<Vec<T>> {
        let fd = FiniteDifference::central(self.dx);
        let grad = fd.first_derivative(field)?;
        Ok(grad.data.as_vec())
    }

    /// Compute gradient of a 2D field (stored row-major)
    pub fn gradient_2d(&self, field: &[T], nx: usize, ny: usize) -> Result<Vec<Vector3<T>>> {
        if field.len() != nx * ny {
            return Err(Error::InvalidConfiguration(
                "Field size doesn't match grid dimensions".to_string()
            ));
        }

        let mut gradients = Vec::with_capacity(nx * ny);
        let two = T::from_f64(constants::TWO).unwrap_or_else(|| T::zero());

        // Use iterator combinators instead of nested loops
        gradients.extend((0..ny).flat_map(|j| {
            let two = two;
            let dx = self.dx;
            let dy = self.dy;
            (0..nx).map(move |i| {
                let idx = j * nx + i;

                // Compute x-derivative
                let dfdx = if i == 0 {
                    // Forward difference
                    (field[idx + 1] - field[idx]) / dx
                } else if i == nx - 1 {
                    // Backward difference
                    (field[idx] - field[idx - 1]) / dx
                } else {
                    // Central difference
                    (field[idx + 1] - field[idx - 1]) / (two * dx)
                };

                // Compute y-derivative
                let dfdy = if j == 0 {
                    // Forward difference
                    (field[idx + nx] - field[idx]) / dy
                } else if j == ny - 1 {
                    // Backward difference
                    (field[idx] - field[idx - nx]) / dy
                } else {
                    // Central difference
                    (field[idx + nx] - field[idx - nx]) / (two * dy)
                };

                Vector3::new(dfdx, dfdy, T::zero())
            })
        }));

        Ok(gradients)
    }

    /// Compute gradient of a 3D field (stored in z-y-x order)
    pub fn gradient_3d(&self, field: &[T], nx: usize, ny: usize, nz: usize) -> Result<Vec<Vector3<T>>> {
        if field.len() != nx * ny * nz {
            return Err(Error::InvalidConfiguration(
                "Field size doesn't match grid dimensions".to_string()
            ));
        }

        let mut gradients = Vec::with_capacity(nx * ny * nz);
        let two = T::from_f64(constants::TWO).unwrap_or_else(|| T::zero());

        // Use iterator combinators for better performance
        gradients.extend((0..nz).flat_map(|k| {
            let two = two;
            let dx = self.dx;
            let dy = self.dy;
            let dz = self.dz;
            (0..ny).flat_map(move |j| {
                let two = two;
                let dx = dx;
                let dy = dy;
                let dz = dz;
                (0..nx).map(move |i| {
                    let idx = k * nx * ny + j * nx + i;

                    // Compute x-derivative
                    let dfdx = if i == 0 {
                        (field[idx + 1] - field[idx]) / dx
                    } else if i == nx - 1 {
                        (field[idx] - field[idx - 1]) / dx
                    } else {
                        (field[idx + 1] - field[idx - 1]) / (two * dx)
                    };

                    // Compute y-derivative
                    let dfdy = if j == 0 {
                        (field[idx + nx] - field[idx]) / dy
                    } else if j == ny - 1 {
                        (field[idx] - field[idx - nx]) / dy
                    } else {
                        (field[idx + nx] - field[idx - nx]) / (two * dy)
                    };

                    // Compute z-derivative
                    let dfdz = if k == 0 {
                        (field[idx + nx * ny] - field[idx]) / dz
                    } else if k == nz - 1 {
                        (field[idx] - field[idx - nx * ny]) / dz
                    } else {
                        (field[idx + nx * ny] - field[idx - nx * ny]) / (two * dz)
                    };

                    Vector3::new(dfdx, dfdy, dfdz)
                })
            })
        }));

        Ok(gradients)
    }

    /// Compute divergence of a vector field in 2D
    pub fn divergence_2d(&self, field: &[Vector3<T>], nx: usize, ny: usize) -> Result<Vec<T>> {
        if field.len() != nx * ny {
            return Err(Error::InvalidConfiguration(
                "Field size doesn't match grid dimensions".to_string()
            ));
        }

        let two = T::from_f64(constants::TWO).unwrap_or_else(|| T::zero());

        let divergence: Vec<T> = (0..ny)
            .flat_map(|j| (0..nx).map(move |i| (i, j)))
            .map(|(i, j)| {
                let idx = j * nx + i;

                // Compute ∂u/∂x
                let dudx = if i == 0 {
                    (field[idx + 1].x - field[idx].x) / self.dx
                } else if i == nx - 1 {
                    (field[idx].x - field[idx - 1].x) / self.dx
                } else {
                    (field[idx + 1].x - field[idx - 1].x) / (two * self.dx)
                };

                // Compute ∂v/∂y
                let dvdy = if j == 0 {
                    (field[idx + nx].y - field[idx].y) / self.dy
                } else if j == ny - 1 {
                    (field[idx].y - field[idx - nx].y) / self.dy
                } else {
                    (field[idx + nx].y - field[idx - nx].y) / (two * self.dy)
                };

                dudx + dvdy
            })
            .collect();

        Ok(divergence)
    }

    /// Compute curl of a vector field in 2D (returns z-component only)
    pub fn curl_2d(&self, field: &[Vector3<T>], nx: usize, ny: usize) -> Result<Vec<T>> {
        if field.len() != nx * ny {
            return Err(Error::InvalidConfiguration(
                "Field size doesn't match grid dimensions".to_string()
            ));
        }

        let mut curl = Vec::with_capacity(nx * ny);
        let two = T::from_f64(constants::TWO).unwrap_or_else(|| T::zero());

        for j in 0..ny {
            for i in 0..nx {
                let idx = j * nx + i;

                // Compute ∂v/∂x
                let dvdx = if i == 0 {
                    (field[idx + 1].y - field[idx].y) / self.dx
                } else if i == nx - 1 {
                    (field[idx].y - field[idx - 1].y) / self.dx
                } else {
                    (field[idx + 1].y - field[idx - 1].y) / (two * self.dx)
                };

                // Compute ∂u/∂y
                let dudy = if j == 0 {
                    (field[idx + nx].x - field[idx].x) / self.dy
                } else if j == ny - 1 {
                    (field[idx].x - field[idx - nx].x) / self.dy
                } else {
                    (field[idx + nx].x - field[idx - nx].x) / (two * self.dy)
                };

                // curl_z = ∂v/∂x - ∂u/∂y
                curl.push(dvdx - dudy);
            }
        }

        Ok(curl)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use nalgebra::Vector3;

    #[test]
    fn test_finite_difference_central() {
        let fd = FiniteDifference::central(1.0);

        // Test on quadratic function: f(x) = x^2, f'(x) = 2x
        let x_values = vec![0.0, 1.0, 4.0, 9.0, 16.0]; // x^2 for x = 0,1,2,3,4
        let derivatives = fd.first_derivative(&x_values).expect("CRITICAL: Add proper error handling");

        // Expected derivatives: [1, 2, 4, 6, 7] (approximate due to boundary conditions)
        assert_relative_eq!(derivatives[1], 2.0, epsilon = 1e-10); // Interior point
        assert_relative_eq!(derivatives[2], 4.0, epsilon = 1e-10); // Interior point
        assert_relative_eq!(derivatives[3], 6.0, epsilon = 1e-10); // Interior point
    }

    #[test]
    fn test_finite_difference_accuracy_polynomial() {
        // Test accuracy on polynomial functions
        // Literature: LeVeque (2007), "Finite Difference Methods for ODEs and PDEs"

        let h = 0.1;
        let fd = FiniteDifference::central(h);

        // Test on quadratic polynomial: f(x) = x², f'(x) = 2x
        let x_points: Vec<f64> = (0..11).map(|i| i as f64 * h).collect();
        let f_values: Vec<f64> = x_points.iter().map(|&x| x.powi(2)).collect();
        let derivatives = fd.first_derivative(&f_values).expect("CRITICAL: Add proper error handling");

        // Check interior points (central difference should be exact for linear derivative)
        for i in 1..derivatives.len()-1 {
            let x = x_points[i];
            let expected = 2.0 * x;
            let computed = derivatives[i];
            assert_relative_eq!(computed, expected, epsilon = 1e-12);
        }
    }

    #[test]
    fn test_finite_difference_convergence_rate() {
        // Test convergence rate of finite difference schemes
        // Literature: Fornberg (1988), "Calculation of Weights in Finite Difference Formulas"

        let test_function = |x: f64| x.sin();
        let test_derivative = |x: f64| x.cos();
        let x_test = 1.0;

        let grid_sizes = vec![0.1, 0.05, 0.025, 0.0125];
        let mut errors = Vec::new();

        for &h in &grid_sizes {
            let fd = FiniteDifference::central(h);

            // Create local grid around test point
            let x_values = vec![x_test - h, x_test, x_test + h];
            let f_values: Vec<f64> = x_values.iter().map(|&x| test_function(x)).collect();

            let derivatives = fd.first_derivative(&f_values).expect("CRITICAL: Add proper error handling");
            let computed = derivatives[1]; // Central point
            let expected = test_derivative(x_test);
            let error = (computed - expected).abs();
            errors.push(error);
        }

        // Check that error decreases quadratically (central difference is O(h²))
        for i in 1..errors.len() {
            let ratio = errors[i-1] / errors[i];
            // Should be approximately 4 (since h is halved each time, error should decrease by factor of 4)
            assert!(ratio > 3.5 && ratio < 4.5, "Convergence rate not quadratic: ratio = {}", ratio);
        }
    }

    #[test]
    fn test_finite_difference_forward() {
        let fd = FiniteDifference::forward(1.0);

        // Test on linear function: f(x) = 2x, f'(x) = 2
        let x_values = vec![0.0, 2.0, 4.0, 6.0, 8.0];
        let derivatives = fd.first_derivative(&x_values).expect("CRITICAL: Add proper error handling");

        // Should be exactly 2.0 for linear function
        for &deriv in derivatives.iter() {
            assert_relative_eq!(deriv, 2.0, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_finite_difference_backward() {
        let fd = FiniteDifference::backward(1.0);

        // Test on linear function: f(x) = 3x, f'(x) = 3
        let x_values = vec![0.0, 3.0, 6.0, 9.0, 12.0];
        let derivatives = fd.first_derivative(&x_values).expect("CRITICAL: Add proper error handling");

        // Should be exactly 3.0 for linear function
        for &deriv in derivatives.iter() {
            assert_relative_eq!(deriv, 3.0, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_second_derivative() {
        let fd = FiniteDifference::central(1.0);

        // Test on quadratic function: f(x) = x^2, f''(x) = 2
        let x_values = vec![0.0, 1.0, 4.0, 9.0, 16.0, 25.0]; // x^2 for x = 0,1,2,3,4,5
        let second_derivatives = fd.second_derivative(&x_values).expect("CRITICAL: Add proper error handling");

        // Should be exactly 2.0 for quadratic function (interior points)
        for i in 1..second_derivatives.len()-1 {
            assert_relative_eq!(second_derivatives[i], 2.0, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_gradient_1d() {
        let grad = Gradient::uniform(1.0);

        // Test on quadratic function: f(x) = x^2, f'(x) = 2x
        let field = vec![0.0, 1.0, 4.0, 9.0, 16.0]; // x^2 for x = 0,1,2,3,4
        let gradient = grad.gradient_1d(&field).expect("CRITICAL: Add proper error handling");

        // Check interior points
        assert_relative_eq!(gradient[1], 2.0, epsilon = 1e-10); // f'(1) = 2
        assert_relative_eq!(gradient[2], 4.0, epsilon = 1e-10); // f'(2) = 4
        assert_relative_eq!(gradient[3], 6.0, epsilon = 1e-10); // f'(3) = 6
    }

    #[test]
    fn test_gradient_2d() {
        let grad = Gradient::uniform(1.0);

        // Test on function f(x,y) = x^2 + y^2
        // ∂f/∂x = 2x, ∂f/∂y = 2y
        let nx = 3;
        let ny = 3;
        let mut field = vec![0.0; nx * ny];

        for j in 0..ny {
            for i in 0..nx {
                let x = i as f64;
                let y = j as f64;
                field[j * nx + i] = x * x + y * y;
            }
        }

        let gradients = grad.gradient_2d(&field, nx, ny).expect("CRITICAL: Add proper error handling");

        // Check center point (1,1): should have gradient (2, 2, 0)
        let center_grad = &gradients[1 * nx + 1];
        assert_relative_eq!(center_grad.x, 2.0, epsilon = 1e-10);
        assert_relative_eq!(center_grad.y, 2.0, epsilon = 1e-10);
        assert_relative_eq!(center_grad.z, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_divergence_2d() {
        let grad = Gradient::uniform(1.0);

        // Test on vector field (x, y, 0) which has divergence = 2
        let nx = 3;
        let ny = 3;
        let mut field = vec![Vector3::zeros(); nx * ny];

        for j in 0..ny {
            for i in 0..nx {
                let x = i as f64;
                let y = j as f64;
                field[j * nx + i] = Vector3::new(x, y, 0.0);
            }
        }

        let divergence = grad.divergence_2d(&field, nx, ny).expect("CRITICAL: Add proper error handling");

        // Divergence should be 2.0 everywhere for this field
        for &div in &divergence {
            assert_relative_eq!(div, 2.0, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_curl_2d() {
        let grad = Gradient::uniform(1.0);

        // Test on vector field (-y, x, 0) which has curl_z = 2
        let nx = 3;
        let ny = 3;
        let mut field = vec![Vector3::zeros(); nx * ny];

        for j in 0..ny {
            for i in 0..nx {
                let x = i as f64;
                let y = j as f64;
                field[j * nx + i] = Vector3::new(-y, x, 0.0);
            }
        }

        let curl = grad.curl_2d(&field, nx, ny).expect("CRITICAL: Add proper error handling");

        // Curl should be 2.0 everywhere for this field (interior points)
        assert_relative_eq!(curl[1 * nx + 1], 2.0, epsilon = 1e-10); // Center point
    }

    #[test]
    fn test_gradient_3d() {
        let grad = Gradient::uniform(1.0);

        // Test on function f(x,y,z) = x^2 + y^2 + z^2
        // ∂f/∂x = 2x, ∂f/∂y = 2y, ∂f/∂z = 2z
        let nx = 3;
        let ny = 3;
        let nz = 3;
        let mut field = vec![0.0; nx * ny * nz];

        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let x = i as f64;
                    let y = j as f64;
                    let z = k as f64;
                    field[k * nx * ny + j * nx + i] = x * x + y * y + z * z;
                }
            }
        }

        let gradients = grad.gradient_3d(&field, nx, ny, nz).expect("CRITICAL: Add proper error handling");

        // Check center point (1,1,1): should have gradient (2, 2, 2)
        let center_grad = &gradients[1 * nx * ny + 1 * nx + 1];
        assert_relative_eq!(center_grad.x, 2.0, epsilon = 1e-10);
        assert_relative_eq!(center_grad.y, 2.0, epsilon = 1e-10);
        assert_relative_eq!(center_grad.z, 2.0, epsilon = 1e-10);
    }

    #[test]
    fn test_finite_difference_schemes() {
        let values = vec![1.0, 4.0, 9.0, 16.0, 25.0]; // x^2 for x = 1,2,3,4,5
        let spacing = 1.0;

        // Test different schemes
        let schemes = [
            FiniteDifferenceScheme::Forward,
            FiniteDifferenceScheme::Backward,
            FiniteDifferenceScheme::Central,
            FiniteDifferenceScheme::ForwardSecondOrder,
            FiniteDifferenceScheme::BackwardSecondOrder,
        ];

        for scheme in &schemes {
            let fd = FiniteDifference::new(*scheme, spacing);
            let result = fd.first_derivative(&values);
            assert!(result.is_ok(), "Scheme {:?} failed", scheme);
        }
    }

    #[test]
    fn test_error_conditions() {
        let grad = Gradient::uniform(1.0);

        // Test insufficient data
        let empty_field = vec![];
        assert!(grad.gradient_1d(&empty_field).is_err());

        let single_point = vec![1.0];
        assert!(grad.gradient_1d(&single_point).is_err());

        // Test mismatched dimensions
        let field_2d = vec![1.0, 2.0, 3.0];
        assert!(grad.gradient_2d(&field_2d, 2, 2).is_err()); // 3 elements for 2x2 grid

        let field_3d = vec![1.0; 8];
        assert!(grad.gradient_3d(&field_3d, 3, 3, 3).is_err()); // 8 elements for 3x3x3 grid
    }
}