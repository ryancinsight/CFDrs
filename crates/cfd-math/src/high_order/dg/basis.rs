//! Basis functions for Discontinuous Galerkin methods.
//!
//! This module provides orthogonal polynomial basis functions and their derivatives
//! for use in DG methods, including both modal and nodal representations.

use super::{DGError, Result};
use nalgebra::{DMatrix, DVector};
use std::f64::consts::PI;

/// Type of basis functions
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BasisType {
    /// Orthogonal basis (Legendre polynomials)
    Orthogonal,
    /// Nodal basis (Lagrange polynomials at Gauss-Lobatto points)
    Nodal,
}

/// Represents a set of basis functions for a DG element
#[derive(Debug, Clone)]
pub struct DGBasis {
    /// Polynomial order
    pub order: usize,
    /// Number of basis functions (order + 1)
    pub num_basis: usize,
    /// Type of basis functions
    pub basis_type: BasisType,
    /// Quadrature points for integration
    pub quad_points: DVector<f64>,
    /// Quadrature weights
    pub quad_weights: DVector<f64>,
    /// Basis function values at quadrature points (num_basis × num_quad_points)
    pub phi: DMatrix<f64>,
    /// Basis function derivatives at quadrature points (num_basis × num_quad_points)
    pub dphi_dx: DMatrix<f64>,
    /// Mass matrix (num_basis × num_basis)
    pub mass_matrix: DMatrix<f64>,
    /// Stiffness matrix (num_basis × num_basis)
    pub stiffness_matrix: DMatrix<f64>,
    /// Differentiation matrix (num_basis × num_basis)
    pub diff_matrix: DMatrix<f64>,
}

impl DGBasis {
    /// Create a new DG basis of given order and type
    pub fn new(order: usize, basis_type: BasisType) -> Result<Self> {
        if order == 0 {
            return Err(DGError::InvalidOrder(order));
        }

        let num_basis = order + 1;

        // Use enough quadrature points to integrate the mass matrix exactly
        // For Orthogonal basis (Legendre), degree is 2*order, so we need at least order+1 Gauss-Legendre
        // or order+2 Gauss-Lobatto points.
        // For Nodal basis, we must use order+1 points to maintain the interpolation property.
        let num_quad = match basis_type {
            BasisType::Orthogonal => order + 2,
            BasisType::Nodal => order + 1,
        };

        // Compute quadrature points and weights
        let (quad_points, quad_weights) = gauss_lobatto_quadrature(num_quad)?;

        // Initialize basis function values and derivatives
        let mut phi = DMatrix::zeros(num_basis, num_quad);
        let mut dphi_dx = DMatrix::zeros(num_basis, num_quad);

        // Compute basis functions and derivatives at quadrature points
        for i in 0..num_basis {
            for (q, &xq) in quad_points.iter().enumerate() {
                match basis_type {
                    BasisType::Orthogonal => {
                        phi[(i, q)] = legendre_poly(i, xq);
                        dphi_dx[(i, q)] = legendre_poly_deriv(i, xq);
                    }
                    BasisType::Nodal => {
                        phi[(i, q)] = lagrange_basis(i, xq, quad_points.as_slice());
                        dphi_dx[(i, q)] = lagrange_basis_deriv(i, xq, quad_points.as_slice());
                    }
                }
            }
        }

        // Compute mass matrix M_ij = ∫ φ_i(x) φ_j(x) dx
        let mut mass_matrix = DMatrix::zeros(num_basis, num_basis);
        for i in 0..num_basis {
            for j in 0..num_basis {
                let mut m_ij = 0.0;
                for q in 0..num_quad {
                    m_ij += quad_weights[q] * phi[(i, q)] * phi[(j, q)];
                }
                mass_matrix[(i, j)] = m_ij;
            }
        }

        // Compute stiffness matrix K_ij = ∫ φ_i'(x) φ_j(x) dx
        let mut stiffness_matrix = DMatrix::zeros(num_basis, num_basis);
        for i in 0..num_basis {
            for j in 0..num_basis {
                let mut k_ij = 0.0;
                for q in 0..num_quad {
                    k_ij += quad_weights[q] * dphi_dx[(i, q)] * phi[(j, q)];
                }
                stiffness_matrix[(i, j)] = k_ij;
            }
        }

        // Compute differentiation matrix D_ij = φ_j'(x_i)
        let mut diff_matrix = DMatrix::zeros(num_basis, num_basis);
        for i in 0..num_basis {
            for j in 0..num_basis {
                diff_matrix[(i, j)] = match basis_type {
                    BasisType::Orthogonal => legendre_poly_deriv(j, quad_points[i]),
                    BasisType::Nodal => {
                        lagrange_basis_deriv(j, quad_points[i], quad_points.as_slice())
                    }
                };
            }
        }

        Ok(Self {
            order,
            num_basis,
            basis_type,
            quad_points,
            quad_weights,
            phi,
            dphi_dx,
            mass_matrix,
            stiffness_matrix,
            diff_matrix,
        })
    }

    /// Evaluate the i-th basis function at point x
    pub fn evaluate_basis(&self, i: usize, x: f64) -> f64 {
        match self.basis_type {
            BasisType::Orthogonal => legendre_poly(i, x),
            BasisType::Nodal => {
                // For nodal basis, we need the original nodes
                lagrange_basis(i, x, self.quad_points.as_slice())
            }
        }
    }

    /// Evaluate the derivative of the i-th basis function at point x
    pub fn evaluate_basis_deriv(&self, i: usize, x: f64) -> f64 {
        match self.basis_type {
            BasisType::Orthogonal => legendre_poly_deriv(i, x),
            BasisType::Nodal => {
                // For nodal basis, we need the original nodes
                lagrange_basis_deriv(i, x, self.quad_points.as_slice())
            }
        }
    }

    /// Project a function onto the DG basis
    pub fn project<F>(&self, f: F) -> DVector<f64>
    where
        F: Fn(f64) -> f64,
    {
        let mut rhs = DVector::zeros(self.num_basis);

        for i in 0..self.num_basis {
            for q in 0..self.quad_points.len() {
                let x = self.quad_points[q];
                let w = self.quad_weights[q];
                rhs[i] += w * f(x) * self.phi[(i, q)];
            }
        }

        // Solve M c = rhs, where M is the mass matrix
        self.mass_matrix
            .clone()
            .lu()
            .solve(&rhs)
            .unwrap_or_else(|| {
                // Fall back to interpolation if the mass matrix is singular
                let mut c = DVector::zeros(self.num_basis);
                for i in 0..self.num_basis {
                    c[i] = f(self.quad_points[i]);
                }
                c
            })
    }
}

/// Compute Legendre polynomial of degree n at point x
pub fn legendre_poly(n: usize, x: f64) -> f64 {
    match n {
        0 => 1.0,
        1 => x,
        _ => {
            let mut p0 = 1.0;
            let mut p1 = x;

            for k in 2..=n {
                let p2 = ((2 * k - 1) as f64 * x * p1 - (k - 1) as f64 * p0) / k as f64;
                p0 = p1;
                p1 = p2;
            }

            p1
        }
    }
}

/// Compute derivative of Legendre polynomial of degree n at point x
pub fn legendre_poly_deriv(n: usize, x: f64) -> f64 {
    if n == 0 {
        return 0.0;
    }

    // Handle endpoints to avoid division by zero
    if (x - 1.0).abs() < 1e-15 {
        return (n * (n + 1)) as f64 / 2.0;
    }
    if (x + 1.0).abs() < 1e-15 {
        let val = (n * (n + 1)) as f64 / 2.0;
        return if (n - 1).is_multiple_of(2) { val } else { -val };
    }

    n as f64 / (1.0 - x * x) * (legendre_poly(n - 1, x) - x * legendre_poly(n, x))
}

/// Compute Lagrange basis function i at point x
pub fn lagrange_basis(i: usize, x: f64, nodes: &[f64]) -> f64 {
    let mut result = 1.0;
    let xi = nodes[i];

    for (j, &xj) in nodes.iter().enumerate() {
        if i != j {
            result *= (x - xj) / (xi - xj);
        }
    }

    result
}

/// Compute derivative of Lagrange basis function i at point x
pub fn lagrange_basis_deriv(i: usize, x: f64, nodes: &[f64]) -> f64 {
    let mut result = 0.0;
    let xi = nodes[i];

    for (j, &xj) in nodes.iter().enumerate() {
        if i == j {
            continue;
        }

        let mut term = 1.0 / (xi - xj);

        for (k, &xk) in nodes.iter().enumerate() {
            if k != i && k != j {
                term *= (x - xk) / (xi - xk);
            }
        }

        result += term;
    }

    result
}

/// Compute Gauss-Lobatto quadrature points and weights
pub fn gauss_lobatto_quadrature(n: usize) -> Result<(DVector<f64>, DVector<f64>)> {
    if n < 2 {
        return Err(DGError::InvalidQuadrature(n));
    }

    let mut nodes = DVector::zeros(n);
    let mut weights = DVector::zeros(n);

    // Endpoints are always included in Gauss-Lobatto quadrature
    nodes[0] = -1.0;
    nodes[n - 1] = 1.0;

    // Initial guess for interior points (Chebyshev points)
    for i in 1..n - 1 {
        nodes[i] = -((i as f64 * PI) / (n - 1) as f64).cos();
    }

    // Newton iteration to find the roots of P_{n-1}'
    let max_iter = 100;
    let tol = 1e-15;

    for i in 1..n - 1 {
        let mut x = nodes[i];
        let mut iter = 0;
        let mut delta = 1.0;

        while delta > tol && iter < max_iter {
            let (p, dp) = legendre_poly_deriv_with_prev(n - 1, x);

            // P''_{n-1} = (2x P'_{n-1} - (n-1)n P_{n-1}) / (1-x^2)
            let denominator = 1.0 - x * x;
            if denominator.abs() < 1e-15 {
                return Err(DGError::NumericalError(format!(
                    "Newton iteration hit endpoint at x={x}"
                )));
            }
            let d2p = (2.0 * x * dp - ((n - 1) * n) as f64 * p) / denominator;

            if d2p.abs() < f64::EPSILON {
                return Err(DGError::NumericalError(format!(
                    "Zero second derivative at x = {x} for n = {}",
                    n - 1
                )));
            }

            let dx = dp / d2p;
            x -= dx;
            delta = dx.abs();
            iter += 1;
        }

        if iter >= max_iter {
            return Err(DGError::NumericalError(format!(
                "Failed to converge for node {i} of {n}"
            )));
        }

        nodes[i] = x;
    }

    // Compute weights
    for i in 0..n {
        let x = nodes[i];
        let p = legendre_poly(n - 1, x);
        if p.abs() < 1e-15 {
            return Err(DGError::NumericalError(format!("P_{}({x}) is zero", n - 1)));
        }
        weights[i] = 2.0 / (n as f64 * (n - 1) as f64) / (p * p);
        if weights[i].is_nan() {
            return Err(DGError::NumericalError(format!(
                "Weight {i} is NaN (n={n}, x={x}, p={p})"
            )));
        }
    }

    Ok((nodes, weights))
}

/// Compute Legendre polynomial and its derivative using recurrence
fn legendre_poly_deriv_with_prev(n: usize, x: f64) -> (f64, f64) {
    if n == 0 {
        return (1.0, 0.0);
    }

    let mut p_prev = 1.0;
    let mut p_curr = x;
    let mut dp_prev = 0.0;
    let mut dp_curr = 1.0;

    for k in 1..n {
        let p_next = ((2 * k + 1) as f64 * x * p_curr - k as f64 * p_prev) / (k + 1) as f64;
        let dp_next = dp_prev + (2 * k + 1) as f64 * p_curr;

        p_prev = p_curr;
        p_curr = p_next;
        dp_prev = dp_curr;
        dp_curr = dp_next;
    }

    (p_curr, dp_curr)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_legendre_polynomials() {
        // Test Legendre polynomial values
        assert_eq!(legendre_poly(0, 0.5), 1.0);
        assert_eq!(legendre_poly(1, 0.5), 0.5);
        assert_relative_eq!(legendre_poly(2, 0.5), -0.125);
        assert_relative_eq!(legendre_poly(3, 0.5), -0.4375);

        // Test orthogonality
        let n = 4;
        // Use more points for better accuracy of P_n^2 (degree 2n)
        let quad = gauss_lobatto_quadrature(2 * n + 1).unwrap();

        for i in 0..=n {
            for j in 0..=n {
                let mut integral = 0.0;

                for k in 0..quad.0.len() {
                    let x = quad.0[k];
                    let w = quad.1[k];
                    integral += w * legendre_poly(i, x) * legendre_poly(j, x);
                }

                // Orthogonality: ∫ P_i(x) P_j(x) dx = 2/(2i+1) δ_ij
                if i == j {
                    let expected = 2.0 / (2.0 * i as f64 + 1.0);
                    assert_relative_eq!(integral, expected, epsilon = 1e-10);
                } else {
                    assert_relative_eq!(integral, 0.0, epsilon = 1e-10);
                }
            }
        }
    }

    #[test]
    fn test_lagrange_basis() {
        let nodes = vec![-1.0, 0.0, 1.0]; // Quadratic elements

        // Test cardinality: φ_i(x_j) = δ_ij
        for (i, &_xi) in nodes.iter().enumerate() {
            for (j, &xj) in nodes.iter().enumerate() {
                let phi = lagrange_basis(i, xj, &nodes);
                let expected = if i == j { 1.0 } else { 0.0 };
                assert_relative_eq!(phi, expected, epsilon = 1e-10);
            }
        }

        // Test interpolation property
        let f = |x: f64| x * x + 2.0 * x + 1.0;
        let x_test = 0.5;
        let mut interpolant = 0.0;

        for (i, &xi) in nodes.iter().enumerate() {
            interpolant += f(xi) * lagrange_basis(i, x_test, &nodes);
        }

        // For quadratic elements, the interpolation should be exact
        assert_relative_eq!(interpolant, f(x_test), epsilon = 1e-10);
    }

    #[test]
    fn test_dg_basis() {
        let order = 3;
        let basis = DGBasis::new(order, BasisType::Orthogonal).unwrap();

        // Test mass matrix properties
        assert_eq!(basis.mass_matrix.nrows(), order + 1);
        assert_eq!(basis.mass_matrix.ncols(), order + 1);

        // Mass matrix should be diagonal for orthogonal bases
        for i in 0..=order {
            for j in 0..=order {
                if i == j {
                    assert!(basis.mass_matrix[(i, i)] > 0.0);
                } else {
                    assert_relative_eq!(basis.mass_matrix[(i, j)], 0.0, epsilon = 1e-10);
                }
            }
        }

        // Test projection
        let f = |x: f64| x.powi(3) - 2.0 * x + 1.0;
        let coeffs = basis.project(f);

        // For order >= 3, the projection should be exact
        for &x in basis.quad_points.iter() {
            let mut approx = 0.0;

            for i in 0..=order {
                approx += coeffs[i] * basis.evaluate_basis(i, x);
            }

            assert_relative_eq!(approx, f(x), epsilon = 1e-10);
        }
    }

    #[test]
    fn test_gauss_lobatto_quadrature() {
        // Test exactness of quadrature rule
        for n in 2..=5 {
            let (points, weights) = gauss_lobatto_quadrature(n).unwrap();

            // Test exactness for polynomials up to degree 2n-3
            for d in 0..=(2 * n - 3) {
                let exact = if d % 2 == 0 {
                    // Integral of x^d from -1 to 1 is 2/(d+1) for even d
                    2.0 / (d as f64 + 1.0)
                } else {
                    // Integral of odd function over symmetric interval is zero
                    0.0
                };

                let mut integral = 0.0;
                for i in 0..n {
                    integral += weights[i] * points[i].powi(d as i32);
                }

                assert!(
                    (integral - exact).abs() < 1e-10,
                    "Failed for n={n}, d={d}: integral={integral}, exact={exact}"
                );
            }
        }
    }
}
