//! Numerical differentiation methods for CFD applications.
//!
//! This module provides finite difference schemes optimized for CFD simulations
//! with support for various boundary conditions and grid types.

use cfd_core::{Error, Result};
use nalgebra::{RealField, Vector3, DVector};
use num_traits::cast::FromPrimitive;

/// Finite difference scheme types
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
pub struct FiniteDifference<T: RealField> {
    scheme: FiniteDifferenceScheme,
    spacing: T,
}

impl<T: RealField + FromPrimitive> FiniteDifference<T> {
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

    /// Compute first derivative at interior points
    pub fn first_derivative(&self, values: &[T]) -> Result<DVector<T>> {
        if values.len() < 2 {
            return Err(Error::InvalidConfiguration(
                "Need at least 2 points for differentiation".to_string()
            ));
        }

        let n = values.len();
        let mut result = DVector::zeros(n);

        match self.scheme {
            FiniteDifferenceScheme::Forward => {
                for i in 0..n-1 {
                    result[i] = (values[i+1].clone() - values[i].clone()) / self.spacing.clone();
                }
                // Use backward difference for last point
                result[n-1] = (values[n-1].clone() - values[n-2].clone()) / self.spacing.clone();
            },
            FiniteDifferenceScheme::Backward => {
                // Use forward difference for first point
                result[0] = (values[1].clone() - values[0].clone()) / self.spacing.clone();
                for i in 1..n {
                    result[i] = (values[i].clone() - values[i-1].clone()) / self.spacing.clone();
                }
            },
            FiniteDifferenceScheme::Central => {
                // Forward difference for first point
                result[0] = (values[1].clone() - values[0].clone()) / self.spacing.clone();

                // Central difference for interior points
                let two = T::from_f64(2.0).unwrap();
                for i in 1..n-1 {
                    result[i] = (values[i+1].clone() - values[i-1].clone()) /
                               (two.clone() * self.spacing.clone());
                }

                // Backward difference for last point
                result[n-1] = (values[n-1].clone() - values[n-2].clone()) / self.spacing.clone();
            },
            FiniteDifferenceScheme::ForwardSecondOrder => {
                if n < 3 {
                    return Err(Error::InvalidConfiguration(
                        "Need at least 3 points for second-order forward difference".to_string()
                    ));
                }

                let two = T::from_f64(2.0).unwrap();
                let three = T::from_f64(3.0).unwrap();
                let four = T::from_f64(4.0).unwrap();

                for i in 0..n-2 {
                    result[i] = (-three.clone() * values[i].clone() +
                                four.clone() * values[i+1].clone() -
                                values[i+2].clone()) / (two.clone() * self.spacing.clone());
                }

                // Use central difference for remaining points
                for i in n-2..n {
                    if i > 0 && i < n-1 {
                        result[i] = (values[i+1].clone() - values[i-1].clone()) /
                                   (two.clone() * self.spacing.clone());
                    }
                }
            },
            FiniteDifferenceScheme::BackwardSecondOrder => {
                if n < 3 {
                    return Err(Error::InvalidConfiguration(
                        "Need at least 3 points for second-order backward difference".to_string()
                    ));
                }

                let two = T::from_f64(2.0).unwrap();
                let three = T::from_f64(3.0).unwrap();
                let four = T::from_f64(4.0).unwrap();

                // Use central difference for first points
                for i in 0..2 {
                    if i > 0 && i < n-1 {
                        result[i] = (values[i+1].clone() - values[i-1].clone()) /
                                   (two.clone() * self.spacing.clone());
                    }
                }

                for i in 2..n {
                    result[i] = (values[i-2].clone() -
                                four.clone() * values[i-1].clone() +
                                three.clone() * values[i].clone()) / (two.clone() * self.spacing.clone());
                }
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
        let h_squared = self.spacing.clone() * self.spacing.clone();

        // Use forward difference for first point
        result[0] = (values[2].clone() - T::from_f64(2.0).unwrap() * values[1].clone() + values[0].clone()) / h_squared.clone();

        // Central difference for interior points
        for i in 1..n-1 {
            result[i] = (values[i+1].clone() - T::from_f64(2.0).unwrap() * values[i].clone() + values[i-1].clone()) / h_squared.clone();
        }

        // Use backward difference for last point
        result[n-1] = (values[n-1].clone() - T::from_f64(2.0).unwrap() * values[n-2].clone() + values[n-3].clone()) / h_squared.clone();

        Ok(result)
    }

    /// Get the finite difference scheme
    pub fn scheme(&self) -> FiniteDifferenceScheme {
        self.scheme
    }

    /// Get the grid spacing
    pub fn spacing(&self) -> T {
        self.spacing.clone()
    }
}

impl<T: RealField + FromPrimitive> Default for FiniteDifference<T> {
    fn default() -> Self {
        Self::central(T::one())
    }
}

/// Gradient computation for multi-dimensional fields
pub struct Gradient<T: RealField> {
    dx: T,
    dy: T,
    dz: T,
}

impl<T: RealField + FromPrimitive> Gradient<T> {
    /// Create new gradient operator
    pub fn new(dx: T, dy: T, dz: T) -> Self {
        Self { dx, dy, dz }
    }

    /// Create uniform spacing gradient operator
    pub fn uniform(spacing: T) -> Self {
        Self::new(spacing.clone(), spacing.clone(), spacing)
    }

    /// Compute gradient of a 1D field
    pub fn gradient_1d(&self, field: &[T]) -> Result<Vec<T>> {
        let fd = FiniteDifference::central(self.dx.clone());
        let grad = fd.first_derivative(field)?;
        Ok(grad.data.as_vec().clone())
    }

    /// Compute gradient of a 2D field (stored row-major)
    pub fn gradient_2d(&self, field: &[T], nx: usize, ny: usize) -> Result<Vec<Vector3<T>>> {
        if field.len() != nx * ny {
            return Err(Error::InvalidConfiguration(
                "Field size doesn't match grid dimensions".to_string()
            ));
        }

        let mut gradients = Vec::with_capacity(nx * ny);
        let two = T::from_f64(2.0).unwrap();

        for j in 0..ny {
            for i in 0..nx {
                let idx = j * nx + i;

                // Compute x-derivative
                let dfdx = if i == 0 {
                    // Forward difference
                    (field[idx + 1].clone() - field[idx].clone()) / self.dx.clone()
                } else if i == nx - 1 {
                    // Backward difference
                    (field[idx].clone() - field[idx - 1].clone()) / self.dx.clone()
                } else {
                    // Central difference
                    (field[idx + 1].clone() - field[idx - 1].clone()) / (two.clone() * self.dx.clone())
                };

                // Compute y-derivative
                let dfdy = if j == 0 {
                    // Forward difference
                    (field[idx + nx].clone() - field[idx].clone()) / self.dy.clone()
                } else if j == ny - 1 {
                    // Backward difference
                    (field[idx].clone() - field[idx - nx].clone()) / self.dy.clone()
                } else {
                    // Central difference
                    (field[idx + nx].clone() - field[idx - nx].clone()) / (two.clone() * self.dy.clone())
                };

                gradients.push(Vector3::new(dfdx, dfdy, T::zero()));
            }
        }

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
        let two = T::from_f64(2.0).unwrap();

        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let idx = k * nx * ny + j * nx + i;

                    // Compute x-derivative
                    let dfdx = if i == 0 {
                        (field[idx + 1].clone() - field[idx].clone()) / self.dx.clone()
                    } else if i == nx - 1 {
                        (field[idx].clone() - field[idx - 1].clone()) / self.dx.clone()
                    } else {
                        (field[idx + 1].clone() - field[idx - 1].clone()) / (two.clone() * self.dx.clone())
                    };

                    // Compute y-derivative
                    let dfdy = if j == 0 {
                        (field[idx + nx].clone() - field[idx].clone()) / self.dy.clone()
                    } else if j == ny - 1 {
                        (field[idx].clone() - field[idx - nx].clone()) / self.dy.clone()
                    } else {
                        (field[idx + nx].clone() - field[idx - nx].clone()) / (two.clone() * self.dy.clone())
                    };

                    // Compute z-derivative
                    let dfdz = if k == 0 {
                        (field[idx + nx * ny].clone() - field[idx].clone()) / self.dz.clone()
                    } else if k == nz - 1 {
                        (field[idx].clone() - field[idx - nx * ny].clone()) / self.dz.clone()
                    } else {
                        (field[idx + nx * ny].clone() - field[idx - nx * ny].clone()) / (two.clone() * self.dz.clone())
                    };

                    gradients.push(Vector3::new(dfdx, dfdy, dfdz));
                }
            }
        }

        Ok(gradients)
    }

    /// Compute divergence of a vector field in 2D
    pub fn divergence_2d(&self, field: &[Vector3<T>], nx: usize, ny: usize) -> Result<Vec<T>> {
        if field.len() != nx * ny {
            return Err(Error::InvalidConfiguration(
                "Field size doesn't match grid dimensions".to_string()
            ));
        }

        let mut divergence = Vec::with_capacity(nx * ny);
        let two = T::from_f64(2.0).unwrap();

        for j in 0..ny {
            for i in 0..nx {
                let idx = j * nx + i;

                // Compute ∂u/∂x
                let dudx = if i == 0 {
                    (field[idx + 1].x.clone() - field[idx].x.clone()) / self.dx.clone()
                } else if i == nx - 1 {
                    (field[idx].x.clone() - field[idx - 1].x.clone()) / self.dx.clone()
                } else {
                    (field[idx + 1].x.clone() - field[idx - 1].x.clone()) / (two.clone() * self.dx.clone())
                };

                // Compute ∂v/∂y
                let dvdy = if j == 0 {
                    (field[idx + nx].y.clone() - field[idx].y.clone()) / self.dy.clone()
                } else if j == ny - 1 {
                    (field[idx].y.clone() - field[idx - nx].y.clone()) / self.dy.clone()
                } else {
                    (field[idx + nx].y.clone() - field[idx - nx].y.clone()) / (two.clone() * self.dy.clone())
                };

                divergence.push(dudx + dvdy);
            }
        }

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
        let two = T::from_f64(2.0).unwrap();

        for j in 0..ny {
            for i in 0..nx {
                let idx = j * nx + i;

                // Compute ∂v/∂x
                let dvdx = if i == 0 {
                    (field[idx + 1].y.clone() - field[idx].y.clone()) / self.dx.clone()
                } else if i == nx - 1 {
                    (field[idx].y.clone() - field[idx - 1].y.clone()) / self.dx.clone()
                } else {
                    (field[idx + 1].y.clone() - field[idx - 1].y.clone()) / (two.clone() * self.dx.clone())
                };

                // Compute ∂u/∂y
                let dudy = if j == 0 {
                    (field[idx + nx].x.clone() - field[idx].x.clone()) / self.dy.clone()
                } else if j == ny - 1 {
                    (field[idx].x.clone() - field[idx - nx].x.clone()) / self.dy.clone()
                } else {
                    (field[idx + nx].x.clone() - field[idx - nx].x.clone()) / (two.clone() * self.dy.clone())
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
        let derivatives = fd.first_derivative(&x_values).unwrap();

        // Expected derivatives: [1, 2, 4, 6, 7] (approximate due to boundary conditions)
        assert_relative_eq!(derivatives[1], 2.0, epsilon = 1e-10); // Interior point
        assert_relative_eq!(derivatives[2], 4.0, epsilon = 1e-10); // Interior point
        assert_relative_eq!(derivatives[3], 6.0, epsilon = 1e-10); // Interior point
    }

    #[test]
    fn test_finite_difference_forward() {
        let fd = FiniteDifference::forward(1.0);

        // Test on linear function: f(x) = 2x, f'(x) = 2
        let x_values = vec![0.0, 2.0, 4.0, 6.0, 8.0];
        let derivatives = fd.first_derivative(&x_values).unwrap();

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
        let derivatives = fd.first_derivative(&x_values).unwrap();

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
        let second_derivatives = fd.second_derivative(&x_values).unwrap();

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
        let gradient = grad.gradient_1d(&field).unwrap();

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

        let gradients = grad.gradient_2d(&field, nx, ny).unwrap();

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

        let divergence = grad.divergence_2d(&field, nx, ny).unwrap();

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

        let curl = grad.curl_2d(&field, nx, ny).unwrap();

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

        let gradients = grad.gradient_3d(&field, nx, ny, nz).unwrap();

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