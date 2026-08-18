//! Basis functions for Discontinuous Galerkin methods.
//!
//! This module provides orthogonal polynomial basis functions and their derivatives
//! for use in DG methods, including both modal and nodal representations.

use crate::error::Result;
use cfd_core::error::{Error, ErrorContext};
use leto::{Array1, Array2};
use leto_ops::{legendre_poly, legendre_poly_and_deriv};

use super::{matrix_solve, matrix_zeros, vector_len, vector_zeros};
use crate::high_order::lgl;

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
    pub quad_points: Array1<f64>,
    /// Quadrature weights
    pub quad_weights: Array1<f64>,
    /// Basis function values at quadrature points (num_basis × num_quad_points)
    pub phi: Array2<f64>,
    /// Basis function derivatives at quadrature points (num_basis × num_quad_points)
    pub dphi_dx: Array2<f64>,
    /// Mass matrix (num_basis × num_basis)
    pub mass_matrix: Array2<f64>,
    /// Stiffness matrix (num_basis × num_basis)
    pub stiffness_matrix: Array2<f64>,
    /// Differentiation matrix (num_basis × num_basis)
    pub diff_matrix: Array2<f64>,
}

impl DGBasis {
    /// Create a new DG basis of given order and type
    pub fn new(order: usize, basis_type: BasisType) -> Result<Self> {
        if order == 0 {
            return Err(Error::InvalidInput(format!(
                "Polynomial order must be at least 1, got {order}"
            )));
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
        let (quad_points, quad_weights) = gauss_lobatto_quadrature(num_quad)
            .context("computing Gauss-Lobatto quadrature for DG basis")?;

        // Initialize basis function values and derivatives
        let mut phi = matrix_zeros(num_basis, num_quad);
        let mut dphi_dx = matrix_zeros(num_basis, num_quad);
        let nodes = quad_points.iter().copied().collect::<Vec<_>>();

        // Compute basis functions and derivatives at quadrature points
        for i in 0..num_basis {
            for (q, &xq) in quad_points.iter().enumerate() {
                match basis_type {
                    BasisType::Orthogonal => {
                        phi[[i, q]] = legendre_poly(i, xq);
                        dphi_dx[[i, q]] = legendre_poly_and_deriv(i, xq).1;
                    }
                    BasisType::Nodal => {
                        phi[[i, q]] = lagrange_basis(i, xq, &nodes);
                        dphi_dx[[i, q]] = lagrange_basis_deriv(i, xq, &nodes);
                    }
                }
            }
        }

        // Compute mass matrix M_ij = ∫ φ_i(x) φ_j(x) dx
        let mut mass_matrix = matrix_zeros(num_basis, num_basis);
        for i in 0..num_basis {
            for j in 0..num_basis {
                let mut m_ij = 0.0;
                for q in 0..num_quad {
                    m_ij += quad_weights[q] * phi[[i, q]] * phi[[j, q]];
                }
                mass_matrix[[i, j]] = m_ij;
            }
        }

        // Compute stiffness matrix K_ij = ∫ φ_i'(x) φ_j(x) dx
        let mut stiffness_matrix = matrix_zeros(num_basis, num_basis);
        for i in 0..num_basis {
            for j in 0..num_basis {
                let mut k_ij = 0.0;
                for q in 0..num_quad {
                    k_ij += quad_weights[q] * dphi_dx[[i, q]] * phi[[j, q]];
                }
                stiffness_matrix[[i, j]] = k_ij;
            }
        }

        // Compute differentiation matrix D_ij = φ_j'(x_i)
        let mut diff_matrix = matrix_zeros(num_basis, num_basis);
        for i in 0..num_basis {
            for j in 0..num_basis {
                diff_matrix[[i, j]] = match basis_type {
                    BasisType::Orthogonal => legendre_poly_and_deriv(j, quad_points[i]).1,
                    BasisType::Nodal => lagrange_basis_deriv(j, quad_points[i], &nodes),
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
                let nodes = self.quad_points.iter().copied().collect::<Vec<_>>();
                lagrange_basis(i, x, &nodes)
            }
        }
    }

    /// Evaluate the derivative of the i-th basis function at point x
    pub fn evaluate_basis_deriv(&self, i: usize, x: f64) -> f64 {
        match self.basis_type {
            BasisType::Orthogonal => legendre_poly_and_deriv(i, x).1,
            BasisType::Nodal => {
                // For nodal basis, we need the original nodes
                let nodes = self.quad_points.iter().copied().collect::<Vec<_>>();
                lagrange_basis_deriv(i, x, &nodes)
            }
        }
    }

    /// Project a function onto the DG basis
    pub fn project<F>(&self, f: F) -> Result<Array1<f64>>
    where
        F: Fn(f64) -> f64,
    {
        let mut rhs = vector_zeros(self.num_basis);

        for i in 0..self.num_basis {
            for q in 0..vector_len(&self.quad_points) {
                let x = self.quad_points[q];
                let w = self.quad_weights[q];
                rhs[i] += w * f(x) * self.phi[[i, q]];
            }
        }

        // Solve M c = rhs, where M is the mass matrix
        matrix_solve(&self.mass_matrix, &rhs)
    }
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
///
/// This is a thin compatibility wrapper around `lgl::lgl_nodes_and_weights`.
/// `n` is the **number of points** (= order + 1).
pub fn gauss_lobatto_quadrature(n: usize) -> Result<(Array1<f64>, Array1<f64>)> {
    if n < 2 {
        return Err(Error::InvalidInput(format!(
            "Gauss-Lobatto quadrature requires at least 2 points, got {n}"
        )));
    }
    let order = n - 1;
    let (nodes_vec, weights_vec) = lgl::lgl_nodes_and_weights(order)
        .context("computing LGL nodes/weights for Gauss-Lobatto quadrature")?;
    let nodes =
        Array1::from_shape_vec([n], nodes_vec).expect("invariant: LGL returns exactly n nodes");
    let weights =
        Array1::from_shape_vec([n], weights_vec).expect("invariant: LGL returns exactly n weights");
    Ok((nodes, weights))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::high_order::dg::{matrix_cols, matrix_rows};
    use eunomia::assert_relative_eq;

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
        let quad = gauss_lobatto_quadrature(2 * n + 1)
            .expect("invariant: quadrature order supports the requested basis");

        for i in 0..=n {
            for j in 0..=n {
                let mut integral = 0.0;

                for k in 0..vector_len(&quad.0) {
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
        let basis =
            DGBasis::new(order, BasisType::Orthogonal).expect("invariant: valid DG basis order");

        // Test mass matrix properties
        assert_eq!(matrix_rows(&basis.mass_matrix), order + 1);
        assert_eq!(matrix_cols(&basis.mass_matrix), order + 1);

        // Mass matrix should be diagonal for orthogonal bases
        for i in 0..=order {
            for j in 0..=order {
                if i == j {
                    assert!(basis.mass_matrix[[i, i]] > 0.0);
                } else {
                    assert_relative_eq!(basis.mass_matrix[[i, j]], 0.0, epsilon = 1e-10);
                }
            }
        }

        // Test projection
        let f = |x: f64| x.powi(3) - 2.0 * x + 1.0;
        let coeffs = basis
            .project(f)
            .expect("invariant: projection function returns finite values");

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
            let (points, weights) =
                gauss_lobatto_quadrature(n).expect("invariant: valid quadrature order");

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
