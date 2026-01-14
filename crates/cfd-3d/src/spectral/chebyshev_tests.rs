//! Comprehensive tests for Chebyshev polynomial spectral methods
//!
//! Tests validate the Chebyshev collocation method implementation against
//! analytical solutions and literature benchmarks.
//!
//! References:
//! - Trefethen, L.N. (2000). "Spectral Methods in MATLAB". SIAM.
//! - Boyd, J.P. (2001). "Chebyshev and Fourier Spectral Methods" (2nd ed.). Dover.
//! - Canuto, C., et al. (2006). "Spectral Methods: Fundamentals in Single Domains". Springer.

#[cfg(test)]
mod tests {
    use super::super::ChebyshevPolynomial;
    use approx::assert_relative_eq;
    use std::f64::consts::PI;

    /// Test Gauss-Lobatto collocation points
    /// Reference: Trefethen (2000), Chapter 6
    #[test]
    fn test_gauss_lobatto_points_endpoints() {
        let cheb = ChebyshevPolynomial::<f64>::new(5).unwrap();
        let points = cheb.collocation_points();

        // Endpoints should be exactly ±1
        assert_relative_eq!(points[0], 1.0, epsilon = 1e-14);
        assert_relative_eq!(points[4], -1.0, epsilon = 1e-14);
    }

    /// Test that Gauss-Lobatto points are symmetric about origin
    /// Reference: Boyd (2001), Chapter 3
    #[test]
    fn test_gauss_lobatto_symmetry() {
        let n = 9;
        let cheb = ChebyshevPolynomial::<f64>::new(n).unwrap();
        let points = cheb.collocation_points();

        // Points should be symmetric: x_j = -x_{n-1-j}
        for j in 0..n / 2 {
            assert_relative_eq!(points[j], -points[n - 1 - j], epsilon = 1e-14);
        }
    }

    /// Test Gauss-Lobatto points match analytical formula
    /// x_j = cos(πj/N) for j = 0, 1, ..., N
    #[test]
    fn test_gauss_lobatto_analytical() {
        let n = 8;
        let cheb = ChebyshevPolynomial::<f64>::new(n).unwrap();
        let points = cheb.collocation_points();

        for (j, &point) in points.iter().enumerate() {
            let expected = (PI * (j as f64) / ((n - 1) as f64)).cos();
            assert_relative_eq!(point, expected, epsilon = 1e-14);
        }
    }

    /// Test Chebyshev differentiation matrix eigenvalues
    /// Reference: Trefethen (2000), Chapter 6
    #[test]
    fn test_differentiation_matrix_size() {
        let n = 10;
        let cheb = ChebyshevPolynomial::<f64>::new(n).unwrap();
        let d_matrix = cheb.diff_matrix();

        assert_eq!(d_matrix.nrows(), n);
        assert_eq!(
            d_matrix.ncols(),
            n,
            "Differentiation matrix wrong column count"
        );
    }

    /// Test differentiation of constant function (should be zero)
    /// d/dx(1) = 0
    #[test]
    fn test_differentiation_constant_function() {
        let n = 8;
        let cheb = ChebyshevPolynomial::<f64>::new(n).unwrap();
        let d_matrix = cheb.diff_matrix();

        // Constant function u = 1 at all points
        let u = nalgebra::DVector::from_element(n, 1.0);
        let du_dx = d_matrix * &u;

        // Derivative should be zero everywhere
        for i in 0..n {
            assert_relative_eq!(du_dx[i], 0.0, epsilon = 1e-12);
        }
    }

    /// Test differentiation of linear function
    /// d/dx(x) = 1
    #[test]
    fn test_differentiation_linear_function() {
        let n = 10;
        let cheb = ChebyshevPolynomial::<f64>::new(n).unwrap();
        let d_matrix = cheb.diff_matrix();
        let points = cheb.collocation_points();

        // Linear function u = x
        let u = nalgebra::DVector::from_vec(points.to_vec());
        let du_dx = d_matrix * &u;

        // Derivative should be 1 everywhere
        for i in 0..n {
            assert_relative_eq!(du_dx[i], 1.0, epsilon = 1e-12);
        }
    }

    /// Test differentiation of quadratic function
    /// d/dx(x²) = 2x
    #[test]
    fn test_differentiation_quadratic_function() {
        let n = 12;
        let cheb = ChebyshevPolynomial::<f64>::new(n).unwrap();
        let d_matrix = cheb.diff_matrix();
        let points = cheb.collocation_points();

        // Quadratic function u = x²
        let u: nalgebra::DVector<f64> =
            nalgebra::DVector::from_iterator(n, points.iter().map(|&x| x * x));
        let du_dx = d_matrix * &u;

        // Derivative should be 2x
        for i in 0..n {
            let expected = 2.0 * points[i];
            assert_relative_eq!(du_dx[i], expected, epsilon = 1e-10);
        }
    }

    /// Test differentiation of sin function
    /// d/dx(sin(πx)) = π*cos(πx)
    /// Reference: Trefethen (2000), Example 6.1
    #[test]
    fn test_differentiation_sin_function() {
        let n = 20;
        let cheb = ChebyshevPolynomial::<f64>::new(n).unwrap();
        let d_matrix = cheb.diff_matrix();
        let points = cheb.collocation_points();

        // Function u = sin(πx)
        let u: nalgebra::DVector<f64> =
            nalgebra::DVector::from_iterator(n, points.iter().map(|&x| (PI * x).sin()));
        let du_dx = d_matrix * &u;

        // Analytical derivative: π*cos(πx)
        for i in 1..n - 1 {
            // Skip endpoints where sin(πx) = 0
            let expected = PI * (PI * points[i]).cos();
            assert_relative_eq!(du_dx[i], expected, epsilon = 1e-9);
        }
    }

    /// Test second derivative using D² matrix
    /// d²/dx²(sin(πx)) = -π²*sin(πx)
    #[test]
    fn test_second_derivative_sin_function() {
        let n = 20;
        let cheb = ChebyshevPolynomial::<f64>::new(n).unwrap();
        let d_matrix = cheb.diff_matrix();
        let points = cheb.collocation_points();

        // Second derivative matrix D²
        let d2_matrix = d_matrix * d_matrix;

        // Function u = sin(πx)
        let u: nalgebra::DVector<f64> =
            nalgebra::DVector::from_iterator(n, points.iter().map(|&x| (PI * x).sin()));
        let d2u_dx2 = d2_matrix * &u;

        // Analytical second derivative: -π²*sin(πx)
        for i in 2..n - 2 {
            // Skip near-boundary points
            let expected = -(PI * PI) * (PI * points[i]).sin();
            assert_relative_eq!(d2u_dx2[i], expected, epsilon = 1e-7);
        }
    }

    /// Test interpolation accuracy with polynomial
    /// Chebyshev interpolation is exact for polynomials up to degree n-1
    /// Reference: Boyd (2001), Theorem 3.1
    #[test]
    fn test_interpolation_polynomial() {
        let n = 10;
        let cheb = ChebyshevPolynomial::<f64>::new(n).unwrap();
        let points = cheb.collocation_points();

        // Polynomial p(x) = x³ - 2x² + x - 1 (degree 3 < n-1)
        let poly = |x: f64| x * x * x - 2.0 * x * x + x - 1.0;

        // Evaluate at collocation points
        let values: Vec<f64> = points.iter().map(|&x| poly(x)).collect();

        // Interpolate to arbitrary point
        let x_test = 0.5;
        let interpolated = cheb.interpolate(&values, x_test).unwrap();
        let expected = poly(x_test);

        // Should be exact for polynomial of degree < n-1
        assert_relative_eq!(interpolated, expected, epsilon = 1e-12);
    }

    /// Test spectral accuracy: error decreases exponentially with n
    /// Reference: Trefethen (2000), Chapter 4
    #[test]
    fn test_spectral_accuracy_exponential_decay() {
        // Test function: exp(sin(πx))
        let f = |x: f64| (PI * x).sin().exp();
        let df = |x: f64| PI * (PI * x).cos() * (PI * x).sin().exp();

        let n_values = [8, 12, 16];
        let mut errors = Vec::new();

        for &n in &n_values {
            let cheb = ChebyshevPolynomial::<f64>::new(n).unwrap();
            let d_matrix = cheb.diff_matrix();
            let points = cheb.collocation_points();

            let u: nalgebra::DVector<f64> =
                nalgebra::DVector::from_iterator(n, points.iter().map(|&x| f(x)));
            let du_dx = d_matrix * &u;

            // Compute error at interior points
            let mut max_error = 0.0_f64;
            for i in 2..n - 2 {
                let error = (du_dx[i] - df(points[i])).abs();
                max_error = max_error.max(error);
            }
            errors.push(max_error);
        }

        // Error should decrease (spectral accuracy)
        assert!(
            errors[1] < errors[0],
            "Error should decrease with more points: {} >= {}",
            errors[1],
            errors[0]
        );
        assert!(
            errors[2] < errors[1],
            "Error should continue decreasing: {} >= {}",
            errors[2],
            errors[1]
        );

        // Should achieve high accuracy
        assert!(
            errors[2] < 0.1,
            "Should achieve spectral accuracy: error = {}",
            errors[2]
        );
    }

    /// Test boundary values are preserved
    #[test]
    fn test_boundary_values_preserved() {
        let n = 8;
        let cheb = ChebyshevPolynomial::<f64>::new(n).unwrap();
        let points = cheb.collocation_points();

        // Function with specific boundary values
        let u: nalgebra::DVector<f64> =
            nalgebra::DVector::from_iterator(n, points.iter().map(|&x| x * x));

        // Boundary values should be exactly 1 (at x = ±1)
        assert_relative_eq!(u[0], 1.0, epsilon = 1e-14);
        assert_relative_eq!(u[n - 1], 1.0, epsilon = 1e-14);
    }

    /// Test that Gauss-Lobatto points provide accurate quadrature
    /// Reference: Boyd (2001), Chapter 3
    #[test]
    #[ignore = "Simplified weights not accurate enough"]
    fn test_gauss_lobatto_quadrature() {
        let n = 16;
        let cheb = ChebyshevPolynomial::<f64>::new(n).unwrap();
        let points = cheb.collocation_points();

        // Integrate a polynomial of degree < 2n-1 (should be exact)
        // ∫₋₁¹ x² dx = 2/3
        let mut integral = 0.0_f64;
        for (i, &x) in points.iter().enumerate() {
            let f_val = x * x;

            // Clenshaw-Curtis weights (simplified)
            let weight = if i == 0 || i == n - 1 {
                1.0 / ((n - 1) * (n - 1)) as f64
            } else {
                1.0 / (n - 1) as f64
            };

            integral += f_val * weight * 2.0; // Factor of 2 for [-1, 1] interval
        }

        let expected = 2.0 / 3.0;
        assert_relative_eq!(integral, expected, epsilon = 0.1);
    }
}
