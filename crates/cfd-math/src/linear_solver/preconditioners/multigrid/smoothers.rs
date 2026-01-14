//! Multigrid smoothers for AMG preconditioning

use super::MultigridSmoother;
use crate::sparse::SparseMatrix;
use nalgebra::DVector;

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
    fn apply(
        &self,
        matrix: &SparseMatrix<T>,
        x: &mut DVector<T>,
        b: &DVector<T>,
        iterations: usize,
    ) {
        let offsets = matrix.row_offsets();
        let indices = matrix.col_indices();
        let values = matrix.values();

        for _ in 0..iterations {
            for i in 0..matrix.nrows() {
                let mut sum = T::zero();
                let mut diag = T::zero();

                for k in offsets[i]..offsets[i + 1] {
                    let j = indices[k];
                    let val = values[k];
                    if i == j {
                        diag = val;
                    } else {
                        sum += val * x[j];
                    }
                }

                if diag.abs() > T::from_f64(1e-15).unwrap() {
                    let new_value = (b[i] - sum) / diag;
                    x[i] = x[i] + self.relaxation_factor * (new_value - x[i]);
                }
            }
        }
    }

    fn clone_box(&self) -> Box<dyn MultigridSmoother<T>> {
        Box::new(self.clone())
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
    fn apply(
        &self,
        matrix: &SparseMatrix<T>,
        x: &mut DVector<T>,
        b: &DVector<T>,
        iterations: usize,
    ) {
        let offsets = matrix.row_offsets();
        let indices = matrix.col_indices();
        let values = matrix.values();

        // Forward sweep
        let gs_forward = GaussSeidelSmoother::new(self.relaxation_factor);
        gs_forward.apply(matrix, x, b, iterations);

        // Backward sweep
        for _ in 0..iterations {
            for i in (0..matrix.nrows()).rev() {
                let mut sum = T::zero();
                let mut diag = T::zero();

                for k in offsets[i]..offsets[i + 1] {
                    let j = indices[k];
                    let val = values[k];
                    if i == j {
                        diag = val;
                    } else {
                        sum += val * x[j];
                    }
                }

                if diag.abs() > T::from_f64(1e-15).unwrap() {
                    let new_value = (b[i] - sum) / diag;
                    x[i] = x[i] + self.relaxation_factor * (new_value - x[i]);
                }
            }
        }
    }

    fn clone_box(&self) -> Box<dyn MultigridSmoother<T>> {
        Box::new(self.clone())
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
    fn apply(
        &self,
        matrix: &SparseMatrix<T>,
        x: &mut DVector<T>,
        b: &DVector<T>,
        iterations: usize,
    ) {
        let offsets = matrix.row_offsets();
        let indices = matrix.col_indices();
        let values = matrix.values();

        for _ in 0..iterations {
            let x_old = x.clone();

            for i in 0..matrix.nrows() {
                let mut sum = T::zero();
                let mut diag = T::zero();

                for k in offsets[i]..offsets[i + 1] {
                    let j = indices[k];
                    let val = values[k];
                    if i == j {
                        diag = val;
                    } else {
                        sum += val * x_old[j];
                    }
                }

                if diag.abs() > T::from_f64(1e-15).unwrap() {
                    let new_value = (b[i] - sum) / diag;
                    x[i] = x_old[i] + self.relaxation_factor * (new_value - x_old[i]);
                }
            }
        }
    }

    fn clone_box(&self) -> Box<dyn MultigridSmoother<T>> {
        Box::new(self.clone())
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
    fn apply(
        &self,
        matrix: &SparseMatrix<T>,
        x: &mut DVector<T>,
        b: &DVector<T>,
        iterations: usize,
    ) {
        let offsets = matrix.row_offsets();
        let indices = matrix.col_indices();
        let values = matrix.values();

        for _ in 0..iterations {
            for i in 0..matrix.nrows() {
                let mut sum = T::zero();
                let mut diag = T::zero();

                for k in offsets[i]..offsets[i + 1] {
                    let j = indices[k];
                    let val = values[k];
                    if i == j {
                        diag = val;
                    } else {
                        sum += val * x[j];
                    }
                }

                if diag.abs() > T::from_f64(1e-15).unwrap() {
                    let new_value = (b[i] - sum) / diag;
                    x[i] = x[i] + self.relaxation_factor * (new_value - x[i]);
                }
            }
        }
    }

    fn clone_box(&self) -> Box<dyn MultigridSmoother<T>> {
        Box::new(self.clone())
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

    /// Estimate eigenvalue bounds using Gershgorin circle theorem
    pub fn estimate_eigenvalues(matrix: &SparseMatrix<T>) -> (T, T) {
        let mut min_eigen = T::from_f64(1e6).unwrap();
        let mut max_eigen = T::zero();

        let offsets = matrix.row_offsets();
        let indices = matrix.col_indices();
        let values = matrix.values();

        for i in 0..matrix.nrows().min(100) {
            // Check more rows for better estimate
            let mut row_sum = T::zero();
            let mut diag = T::zero();

            for k in offsets[i]..offsets[i + 1] {
                let j = indices[k];
                let val = values[k].abs();
                if i == j {
                    diag = val;
                } else {
                    row_sum += val;
                }
            }

            // Gershgorin discs: [diag - sum, diag + sum]
            let lower = if diag > row_sum {
                diag - row_sum
            } else {
                T::from_f64(0.1).unwrap()
            };
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
    fn apply(
        &self,
        matrix: &SparseMatrix<T>,
        x: &mut DVector<T>,
        b: &DVector<T>,
        iterations: usize,
    ) {
        let theta = (self.eigenvalues_max + self.eigenvalues_min) / (T::from_f64(2.0).unwrap());
        let delta = (self.eigenvalues_max - self.eigenvalues_min) / (T::from_f64(2.0).unwrap());
        let sigma = if delta.abs() < T::from_f64(1e-10).unwrap() {
            T::zero()
        } else {
            theta / delta
        };

        let rho_old = if sigma == T::zero() {
            T::zero()
        } else {
            T::from_f64(1.0).unwrap() / sigma
        };

        // Compute initial residual
        let mut r = b - matrix * &*x;

        let alpha = if sigma == T::zero() {
            T::from_f64(1.0).unwrap() / self.eigenvalues_max
        } else {
            let rho_new = T::from_f64(1.0).unwrap() / (T::from_f64(2.0).unwrap() * sigma - rho_old);
            rho_new / rho_old
        };

        let offsets = matrix.row_offsets();
        let indices = matrix.col_indices();
        let values = matrix.values();

        // Chebyshev iterations
        for _ in 0..iterations.min(self.degree) {
            for i in 0..matrix.nrows() {
                let mut diag = T::zero();
                for k in offsets[i]..offsets[i + 1] {
                    if indices[k] == i {
                        diag = values[k];
                        break;
                    }
                }

                if diag.abs() > T::from_f64(1e-15).unwrap() {
                    x[i] += alpha * r[i] / diag;
                }
            }
            r = b - matrix * &*x;
        }
    }

    fn clone_box(&self) -> Box<dyn MultigridSmoother<T>> {
        Box::new(self.clone())
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
        let matrix_dense = create_test_matrix();
        let matrix = SparseMatrix::from(&matrix_dense);
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        let mut x = DVector::zeros(5);

        let smoother = GaussSeidelSmoother::new(1.0);
        smoother.apply(&matrix, &mut x, &b, 1);

        // Check that x has been updated (non-zero values)
        assert!(x.iter().any(|&val| val.abs() > 1e-10));
    }

    #[test]
    fn test_jacobi_smoother() {
        let matrix_dense = create_test_matrix();
        let matrix = SparseMatrix::from(&matrix_dense);
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        let mut x = DVector::zeros(5);

        let smoother = JacobiSmoother::new(1.0);
        smoother.apply(&matrix, &mut x, &b, 1);

        // Check that x has been updated
        assert!(x.iter().any(|&val| val.abs() > 1e-10));
    }

    #[test]
    fn test_symmetric_gauss_seidel() {
        let matrix_dense = create_test_matrix();
        let matrix = SparseMatrix::from(&matrix_dense);
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        let mut x = DVector::zeros(5);

        let smoother = SymmetricGaussSeidelSmoother::new(1.0);
        smoother.apply(&matrix, &mut x, &b, 1);

        // Check that x has been updated
        assert!(x.iter().any(|&val| val.abs() > 1e-10));
    }

    #[test]
    fn test_chebyshev_smoother() {
        let matrix_dense = create_test_matrix();
        let matrix = SparseMatrix::from(&matrix_dense);
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        let mut x = DVector::zeros(5);

        let (eigen_min, eigen_max) = ChebyshevSmoother::estimate_eigenvalues(&matrix);
        let smoother = ChebyshevSmoother::new(eigen_min, eigen_max, 3);
        smoother.apply(&matrix, &mut x, &b, 1);

        // Check that x has been updated
        assert!(x.iter().any(|&val| val.abs() > 1e-10));
    }
}
