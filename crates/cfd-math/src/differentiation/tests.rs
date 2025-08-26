#[cfg(test)]
mod tests {
    use super::super::*;
    use approx::assert_relative_eq;
    use nalgebra::Vector3;

    #[test]
    fn test_finite_difference_central() {
        let fd = FiniteDifference::central(1.0);
        // Test on quadratic function: f(x) = x^2, f'(x) = 2x
        let x_values = vec![0.0, 1.0, 4.0, 9.0, 16.0]; // x^2 for x = 0,1,2,3,4
        let derivatives = fd
            .first_derivative(&x_values)
            .expect("CRITICAL: Add proper error handling");
        // Expected derivatives: [1, 2, 4, 6, 7] (approximate due to boundary conditions)
        assert_relative_eq!(derivatives[1], 2.0, epsilon = 1e-10); // Interior point
        assert_relative_eq!(derivatives[2], 4.0, epsilon = 1e-10); // Interior point
        assert_relative_eq!(derivatives[3], 6.0, epsilon = 1e-10); // Interior point
    }
    fn test_finite_difference_accuracy_polynomial() {
        // Test accuracy on polynomial functions
        // Literature: LeVeque (2007), "Finite Difference Methods for ODEs and PDEs"
        let h = 0.1;
        let fd = FiniteDifference::central(h);
        // Test on quadratic polynomial: f(x) = x², f'(x) = 2x
        let x_points: Vec<f64> = (0..11).map(|i| i as f64 * h).collect();
        let f_values: Vec<f64> = x_points.iter().map(|&x| x.powi(2)).collect();
            .first_derivative(&f_values)
        // Check interior points (central difference should be exact for linear derivative)
        for i in 1..derivatives.len() - 1 {
            let x = x_points[i];
            let expected = 2.0 * x;
            let computed = derivatives[i];
            assert_relative_eq!(computed, expected, epsilon = 1e-12);
        }
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
            let derivatives = fd
                .first_derivative(&f_values)
                .expect("CRITICAL: Add proper error handling");
            let computed = derivatives[1]; // Central point
            let expected = test_derivative(x_test);
            let error = (computed - expected).abs();
            errors.push(error);
        // Check that error decreases quadratically (central difference is O(h²))
        for i in 1..errors.len() {
            let ratio = errors[i - 1] / errors[i];
            // Should be approximately 4 (since h is halved each time, error should decrease by factor of 4)
            assert!(
                ratio > 3.5 && ratio < 4.5,
                "Convergence rate not quadratic: ratio = {}",
                ratio
            );
    fn test_finite_difference_forward() {
        let fd = FiniteDifference::forward(1.0);
        // Test on linear function: f(x) = 2x, f'(x) = 2
        let x_values = vec![0.0, 2.0, 4.0, 6.0, 8.0];
        // Should be exactly 2.0 for linear function
        for &deriv in derivatives.iter() {
            assert_relative_eq!(deriv, 2.0, epsilon = 1e-10);
    }

    fn test_finite_difference_backward() {
        let fd = FiniteDifference::backward(1.0);
        // Test on linear function: f(x) = 3x, f'(x) = 3
        let x_values = vec![0.0, 3.0, 6.0, 9.0, 12.0];
        // Should be exactly 3.0 for linear function
            assert_relative_eq!(deriv, 3.0, epsilon = 1e-10);
    }

    fn test_second_derivative() {
        // Test on quadratic function: f(x) = x^2, f''(x) = 2
        let x_values = vec![0.0, 1.0, 4.0, 9.0, 16.0, 25.0]; // x^2 for x = 0,1,2,3,4,5
        let second_derivatives = fd
            .second_derivative(&x_values)
        // Should be exactly 2.0 for quadratic function (interior points)
        for i in 1..second_derivatives.len() - 1 {
            assert_relative_eq!(second_derivatives[i], 2.0, epsilon = 1e-10);
    }

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
        let gradients = grad
            .gradient_2d(&field, nx, ny)
        // Check center point (1,1): should have gradient (2, 2, 0)
        let center_grad = &gradients[1 * nx + 1];
        assert_relative_eq!(center_grad.x, 2.0, epsilon = 1e-10);
        assert_relative_eq!(center_grad.y, 2.0, epsilon = 1e-10);
    fn test_divergence_2d() {
        // Test on vector field (x, y, 0) which has divergence = 2
        let mut field = vec![Vector3::zeros(); nx * ny];
                field[j * nx + i] = Vector3::new(x, y, 0.0);
        let divergence = grad
            .divergence_2d(&field, nx, ny)
        // Divergence should be 2.0 everywhere for this field
        for &div in &divergence {
            assert_relative_eq!(div, 2.0, epsilon = 1e-10);
    }

    fn test_curl_2d() {
        // Test on vector field (-y, x, 0) which has curl_z = 2
                field[j * nx + i] = Vector3::new(-y, x, 0.0);
        let curl = grad
            .curl_2d(&field, nx, ny)
        // Curl should be 2.0 everywhere for this field (interior points)
        assert_relative_eq!(curl[1 * nx + 1], 2.0, epsilon = 1e-10); // Center point
    }

    fn test_gradient_3d() {
        // Test on function f(x,y,z) = x^2 + y^2 + z^2
        // ∂f/∂x = 2x, ∂f/∂y = 2y, ∂f/∂z = 2z
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
            .gradient_3d(&field, nx, ny, nz)
        // Check center point (1,1,1): should have gradient (2, 2, 2)
        let center_grad = &gradients[1 * nx * ny + 1 * nx + 1];
        assert_relative_eq!(center_grad.z, 2.0, epsilon = 1e-10);
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
    fn test_error_conditions() {
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
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
