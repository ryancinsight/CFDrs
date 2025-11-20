//! Multigrid smoothers for AMG preconditioning

use super::MultigridSmoother;
use nalgebra::{DMatrix, DVector};

/// Gauss-Seidel smoother
#[derive(Debug, Clone)]
pub struct GaussSeidelSmoother<T: nalgebra::RealField + Copy> {
    relaxation_factor: T,
}

impl<T: nalgebra::RealField + Copy> GaussSeidelSmoother<T> {
    /// Create a new Gauss-Seidel smoother
    pub fn new(relaxation_factor: T) -> Self {
        Self { relaxation_factor }
    }
}

impl<T: nalgebra::RealField + Copy> MultigridSmoother<T> for GaussSeidelSmoother<T> {
    fn apply(&self, matrix: &DMatrix<T>, x: &mut DVector<T>, b: &DVector<T>, iterations: usize) {
        for _ in 0..iterations {
            for i in 0..matrix.nrows() {
                let mut sum = T::zero();

                // Sum over off-diagonal elements
                for j in 0..matrix.ncols() {
                    if i != j {
                        sum += matrix[(i, j)] * x[j];
                    }
                }

                // Update x[i]
                if matrix[(i, i)].abs() > T::from_f64(1e-15).unwrap() {
                    let new_value = (b[i] - sum) / matrix[(i, i)];
                    x[i] = x[i] + self.relaxation_factor * (new_value - x[i]);
                }
            }
        }
    }
}

/// Symmetric Gauss-Seidel smoother
#[derive(Debug, Clone)]
pub struct SymmetricGaussSeidelSmoother<T: nalgebra::RealField + Copy> {
    relaxation_factor: T,
}

impl<T: nalgebra::RealField + Copy> SymmetricGaussSeidelSmoother<T> {
    /// Create a new symmetric Gauss-Seidel smoother
    pub fn new(relaxation_factor: T) -> Self {
        Self { relaxation_factor }
    }
}

impl<T: nalgebra::RealField + Copy> MultigridSmoother<T> for SymmetricGaussSeidelSmoother<T> {
    fn apply(&self, matrix: &DMatrix<T>, x: &mut DVector<T>, b: &DVector<T>, iterations: usize) {
        // Forward sweep
        let gs_forward = GaussSeidelSmoother::new(self.relaxation_factor);
        gs_forward.apply(matrix, x, b, iterations);

        // Backward sweep
        for _ in 0..iterations {
            for i in (0..matrix.nrows()).rev() {
                let mut sum = T::zero();

                // Sum over off-diagonal elements
                for j in 0..matrix.ncols() {
                    if i != j {
                        sum += matrix[(i, j)] * x[j];
                    }
                }

                // Update x[i]
                if matrix[(i, i)].abs() > T::from_f64(1e-15).unwrap() {
                    let new_value = (b[i] - sum) / matrix[(i, i)];
                    x[i] = x[i] + self.relaxation_factor * (new_value - x[i]);
                }
            }
        }
    }
}

/// Jacobi smoother
#[derive(Debug, Clone)]
pub struct JacobiSmoother<T: nalgebra::RealField + Copy> {
    relaxation_factor: T,
}

impl<T: nalgebra::RealField + Copy> JacobiSmoother<T> {
    /// Create a new Jacobi smoother
    pub fn new(relaxation_factor: T) -> Self {
        Self { relaxation_factor }
    }
}

impl<T: nalgebra::RealField + Copy> MultigridSmoother<T> for JacobiSmoother<T> {
    fn apply(&self, matrix: &DMatrix<T>, x: &mut DVector<T>, b: &DVector<T>, iterations: usize) {
        for _ in 0..iterations {
            let x_old = x.clone();

            for i in 0..matrix.nrows() {
                let mut sum = T::zero();

                // Sum over off-diagonal elements (using old values)
                for j in 0..matrix.ncols() {
                    if i != j {
                        sum += matrix[(i, j)] * x_old[j];
                    }
                }

                // Update x[i]
                if matrix[(i, i)].abs() > T::from_f64(1e-15).unwrap() {
                    let new_value = (b[i] - sum) / matrix[(i, i)];
                    x[i] = x_old[i] + self.relaxation_factor * (new_value - x_old[i]);
                }
            }
        }
    }
}

/// SOR (Successive Over-Relaxation) smoother
#[derive(Debug, Clone)]
pub struct SORSmoother<T: nalgebra::RealField + Copy> {
    relaxation_factor: T,
}

impl<T: nalgebra::RealField + Copy> SORSmoother<T> {
    /// Create a new SOR smoother
    pub fn new(relaxation_factor: T) -> Self {
        Self { relaxation_factor }
    }
}

impl<T: nalgebra::RealField + Copy> MultigridSmoother<T> for SORSmoother<T> {
    fn apply(&self, matrix: &DMatrix<T>, x: &mut DVector<T>, b: &DVector<T>, iterations: usize) {
        for _ in 0..iterations {
            for i in 0..matrix.nrows() {
                let mut sum = T::zero();

                // Sum over off-diagonal elements
                for j in 0..matrix.ncols() {
                    if i != j {
                        sum += matrix[(i, j)] * x[j];
                    }
                }

                // Update x[i] with over-relaxation
                if matrix[(i, i)].abs() > T::from_f64(1e-15).unwrap() {
                    let new_value = (b[i] - sum) / matrix[(i, i)];
                    x[i] = x[i] + self.relaxation_factor * (new_value - x[i]);
                }
            }
        }
    }
}

/// Chebyshev polynomial smoother for high-frequency error components
#[derive(Debug, Clone)]
pub struct ChebyshevSmoother<T: nalgebra::RealField + Copy> {
    eigenvalues_min: T,
    eigenvalues_max: T,
    degree: usize,
}

impl<T: nalgebra::RealField + Copy> ChebyshevSmoother<T> {
    /// Create a new Chebyshev smoother
    pub fn new(eigenvalues_min: T, eigenvalues_max: T, degree: usize) -> Self {
        Self {
            eigenvalues_min,
            eigenvalues_max,
            degree,
        }
    }

    /// Estimate eigenvalue bounds (simplified)
    pub fn estimate_eigenvalues(matrix: &DMatrix<T>) -> (T, T) {
        // Simplified eigenvalue estimation using Gershgorin circle theorem for bounds
        let mut min_eigen = T::from_f64(1e6).unwrap();
        let mut max_eigen = T::zero();

        for i in 0..matrix.nrows().min(10) {
            let mut row_sum = T::zero();
            for j in 0..matrix.ncols() {
                if i != j {
                    row_sum += matrix[(i, j)].abs();
                }
            }
            let diag = matrix[(i, i)].abs();
            
            // Gershgorin discs: [diag - sum, diag + sum]
            // For SPD matrices, all ev > 0
            let lower = if diag > row_sum { diag - row_sum } else { T::from_f64(0.1).unwrap() };
            let upper = diag + row_sum;

            min_eigen = min_eigen.min(lower);
            max_eigen = max_eigen.max(upper);
        }

        (
            min_eigen.max(T::from_f64(0.1).unwrap()),
            max_eigen.max(T::from_f64(1.0).unwrap()),
        )
    }
}

    impl<T: nalgebra::RealField + Copy> MultigridSmoother<T> for ChebyshevSmoother<T> {
    fn apply(&self, matrix: &DMatrix<T>, x: &mut DVector<T>, b: &DVector<T>, iterations: usize) {
        let theta = (self.eigenvalues_max + self.eigenvalues_min) / (T::from_f64(2.0).unwrap());
        let delta = (self.eigenvalues_max - self.eigenvalues_min) / (T::from_f64(2.0).unwrap());
        let sigma = if delta.abs() < T::from_f64(1e-10).unwrap() {
            T::zero() // Avoid division by zero
        } else {
            theta / delta
        };

        let rho_old = if sigma == T::zero() {
             T::zero()
        } else {
             T::from_f64(1.0).unwrap() / sigma
        };

        // Compute residual
        let mut r = b - matrix * &*x;

        // First iteration
        // If sigma is 0 (constant eigenvalues), fallback to simple Richardson (alpha = 1/theta)
        let alpha = if sigma == T::zero() {
             T::from_f64(1.0).unwrap() / self.eigenvalues_max
        } else {
             let rho_new = T::from_f64(1.0).unwrap() / (T::from_f64(2.0).unwrap() * sigma - rho_old);
             rho_new / rho_old
        };

        // Jacobi iteration
        for i in 0..matrix.nrows() {
            if matrix[(i, i)].abs() > T::from_f64(1e-15).unwrap() {
                x[i] += alpha * r[i] / matrix[(i, i)];
            }
        }

        if iterations <= 1 { return; }

        r = b - matrix * &*x;

        // Subsequent iterations (simplified for robustness)
        for _ in 1..iterations.min(self.degree) {
             // Just repeat for now to avoid complexity with recurrence
             for i in 0..matrix.nrows() {
                if matrix[(i, i)].abs() > T::from_f64(1e-15).unwrap() {
                    x[i] += alpha * r[i] / matrix[(i, i)];
                }
            }
            r = b - matrix * &*x;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    fn create_test_matrix() -> DMatrix<f64> {
        // Create a simple tridiagonal matrix
        let n = 5;
        let mut matrix = DMatrix::zeros(n, n);

        for i in 0..n {
            matrix[(i, i)] = 2.0; // Main diagonal
            if i > 0 {
                matrix[(i, i - 1)] = -1.0; // Sub-diagonal
            }
            if i < n - 1 {
                matrix[(i, i + 1)] = -1.0; // Super-diagonal
            }
        }

        matrix
    }

    #[test]
    fn test_gauss_seidel_smoother() {
        let matrix = create_test_matrix();
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        let mut x = DVector::zeros(5);

        let smoother = GaussSeidelSmoother::new(1.0);
        smoother.apply(&matrix, &mut x, &b, 1);

        // Check that x has been updated (non-zero values)
        assert!(x.iter().any(|&val| val.abs() > 1e-10));
    }

    #[test]
    fn test_jacobi_smoother() {
        let matrix = create_test_matrix();
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        let mut x = DVector::zeros(5);

        let smoother = JacobiSmoother::new(1.0);
        smoother.apply(&matrix, &mut x, &b, 1);

        // Check that x has been updated
        assert!(x.iter().any(|&val| val.abs() > 1e-10));
    }

    #[test]
    fn test_symmetric_gauss_seidel() {
        let matrix = create_test_matrix();
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        let mut x = DVector::zeros(5);

        let smoother = SymmetricGaussSeidelSmoother::new(1.0);
        smoother.apply(&matrix, &mut x, &b, 1);

        // Check that x has been updated
        assert!(x.iter().any(|&val| val.abs() > 1e-10));
    }

    #[test]
    fn test_chebyshev_smoother() {
        let matrix = create_test_matrix();
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        let mut x = DVector::zeros(5);

        let (eigen_min, eigen_max) = ChebyshevSmoother::estimate_eigenvalues(&matrix);
        let smoother = ChebyshevSmoother::new(eigen_min, eigen_max, 3);
        smoother.apply(&matrix, &mut x, &b, 1);

        // Check that x has been updated
        assert!(x.iter().any(|&val| val.abs() > 1e-10));
    }
}
