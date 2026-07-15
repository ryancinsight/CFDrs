//! Multigrid smoothers for AMG preconditioning

use super::{MultigridSmoother, MultigridVector, SparseMatrix};
use eunomia::{FloatElement, NumericElement, RealField as EunomiaRealField};
use leto_ops::{spmv as leto_spmv, Scalar as LetoScalar};

#[inline]
fn from_f64<T: FloatElement>(value: f64) -> T {
    <T as FloatElement>::from_f64(value)
}

#[inline]
fn diagonal_epsilon<T: FloatElement>() -> T {
    from_f64(1e-15)
}

fn residual<T>(
    matrix: &SparseMatrix<T>,
    b: &MultigridVector<T>,
    x: &MultigridVector<T>,
) -> MultigridVector<T>
where
    T: EunomiaRealField + Copy + LetoScalar,
{
    let applied =
        leto_spmv(matrix, &x.view()).expect("invariant: smoother SpMV dimensions are valid");
    let mut out = MultigridVector::zeros([b.shape()[0]]);
    for i in 0..b.shape()[0] {
        out[i] = b[i] - applied[i];
    }
    out
}

/// Gauss-Seidel smoother
#[derive(Debug, Clone)]
pub struct GaussSeidelSmoother<T: EunomiaRealField + Copy> {
    relaxation_factor: T,
}

impl<T: EunomiaRealField + Copy> GaussSeidelSmoother<T> {
    /// Create a new Gauss-Seidel smoother
    pub fn new(relaxation_factor: T) -> Self {
        Self { relaxation_factor }
    }
}

impl<T: EunomiaRealField + Copy + FloatElement + LetoScalar> MultigridSmoother<T>
    for GaussSeidelSmoother<T>
{
    fn apply(
        &self,
        matrix: &SparseMatrix<T>,
        x: &mut MultigridVector<T>,
        b: &MultigridVector<T>,
        iterations: usize,
    ) {
        let offsets = matrix.row_ptr();
        let indices = matrix.col_indices();
        let values = matrix.values();

        for _ in 0..iterations {
            for i in 0..matrix.nrows() {
                let mut sum = <T as NumericElement>::ZERO;
                let mut diag = <T as NumericElement>::ZERO;

                for k in offsets[i]..offsets[i + 1] {
                    let j = indices[k];
                    let val = values[k];
                    if i == j {
                        diag = val;
                    } else {
                        sum += val * x[j];
                    }
                }

                if NumericElement::abs(diag) > diagonal_epsilon() {
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
pub struct SymmetricGaussSeidelSmoother<T: EunomiaRealField + Copy> {
    relaxation_factor: T,
}

impl<T: EunomiaRealField + Copy> SymmetricGaussSeidelSmoother<T> {
    /// Create a new symmetric Gauss-Seidel smoother
    pub fn new(relaxation_factor: T) -> Self {
        Self { relaxation_factor }
    }
}

impl<T: EunomiaRealField + Copy + FloatElement + LetoScalar> MultigridSmoother<T>
    for SymmetricGaussSeidelSmoother<T>
{
    fn apply(
        &self,
        matrix: &SparseMatrix<T>,
        x: &mut MultigridVector<T>,
        b: &MultigridVector<T>,
        iterations: usize,
    ) {
        let offsets = matrix.row_ptr();
        let indices = matrix.col_indices();
        let values = matrix.values();

        // Forward sweep
        let gs_forward = GaussSeidelSmoother::new(self.relaxation_factor);
        gs_forward.apply(matrix, x, b, iterations);

        // Backward sweep
        for _ in 0..iterations {
            for i in (0..matrix.nrows()).rev() {
                let mut sum = <T as NumericElement>::ZERO;
                let mut diag = <T as NumericElement>::ZERO;

                for k in offsets[i]..offsets[i + 1] {
                    let j = indices[k];
                    let val = values[k];
                    if i == j {
                        diag = val;
                    } else {
                        sum += val * x[j];
                    }
                }

                if NumericElement::abs(diag) > diagonal_epsilon() {
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
pub struct JacobiSmoother<T: EunomiaRealField + Copy> {
    relaxation_factor: T,
}

impl<T: EunomiaRealField + Copy> JacobiSmoother<T> {
    /// Create a new Jacobi smoother
    pub fn new(relaxation_factor: T) -> Self {
        Self { relaxation_factor }
    }
}

impl<T: EunomiaRealField + Copy + FloatElement + LetoScalar> MultigridSmoother<T>
    for JacobiSmoother<T>
{
    fn apply(
        &self,
        matrix: &SparseMatrix<T>,
        x: &mut MultigridVector<T>,
        b: &MultigridVector<T>,
        iterations: usize,
    ) {
        let offsets = matrix.row_ptr();
        let indices = matrix.col_indices();
        let values = matrix.values();

        for _ in 0..iterations {
            let x_old = x.clone();

            for i in 0..matrix.nrows() {
                let mut sum = <T as NumericElement>::ZERO;
                let mut diag = <T as NumericElement>::ZERO;

                for k in offsets[i]..offsets[i + 1] {
                    let j = indices[k];
                    let val = values[k];
                    if i == j {
                        diag = val;
                    } else {
                        sum += val * x_old[j];
                    }
                }

                if NumericElement::abs(diag) > diagonal_epsilon() {
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
pub struct SORSmoother<T: EunomiaRealField + Copy> {
    relaxation_factor: T,
}

impl<T: EunomiaRealField + Copy> SORSmoother<T> {
    /// Create a new SOR smoother
    pub fn new(relaxation_factor: T) -> Self {
        Self { relaxation_factor }
    }
}

impl<T: EunomiaRealField + Copy + FloatElement + LetoScalar> MultigridSmoother<T>
    for SORSmoother<T>
{
    fn apply(
        &self,
        matrix: &SparseMatrix<T>,
        x: &mut MultigridVector<T>,
        b: &MultigridVector<T>,
        iterations: usize,
    ) {
        let offsets = matrix.row_ptr();
        let indices = matrix.col_indices();
        let values = matrix.values();

        for _ in 0..iterations {
            for i in 0..matrix.nrows() {
                let mut sum = <T as NumericElement>::ZERO;
                let mut diag = <T as NumericElement>::ZERO;

                for k in offsets[i]..offsets[i + 1] {
                    let j = indices[k];
                    let val = values[k];
                    if i == j {
                        diag = val;
                    } else {
                        sum += val * x[j];
                    }
                }

                if NumericElement::abs(diag) > diagonal_epsilon() {
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
pub struct ChebyshevSmoother<T: EunomiaRealField + Copy> {
    eigenvalues_min: T,
    eigenvalues_max: T,
    degree: usize,
}

impl<T: EunomiaRealField + Copy + FloatElement + LetoScalar> ChebyshevSmoother<T> {
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
        let mut min_eigen: T = from_f64(1e6);
        let mut max_eigen = <T as NumericElement>::ZERO;

        let offsets = matrix.row_ptr();
        let indices = matrix.col_indices();
        let values = matrix.values();

        for i in 0..matrix.nrows().min(100) {
            // Check more rows for better estimate
            let mut row_sum = <T as NumericElement>::ZERO;
            let mut diag = <T as NumericElement>::ZERO;

            for k in offsets[i]..offsets[i + 1] {
                let j = indices[k];
                let val = NumericElement::abs(values[k]);
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
                from_f64::<T>(0.1)
            };
            let upper = diag + row_sum;

            min_eigen = min_eigen.min_scalar(lower);
            max_eigen = max_eigen.max_scalar(upper);
        }

        (
            min_eigen.max_scalar(from_f64::<T>(0.1)),
            max_eigen.max_scalar(from_f64::<T>(1.0)),
        )
    }
}

impl<T: EunomiaRealField + Copy + FloatElement + LetoScalar> MultigridSmoother<T>
    for ChebyshevSmoother<T>
{
    fn apply(
        &self,
        matrix: &SparseMatrix<T>,
        x: &mut MultigridVector<T>,
        b: &MultigridVector<T>,
        iterations: usize,
    ) {
        let theta = (self.eigenvalues_max + self.eigenvalues_min) / from_f64::<T>(2.0);
        let delta = (self.eigenvalues_max - self.eigenvalues_min) / from_f64::<T>(2.0);
        let sigma = if NumericElement::abs(delta) < from_f64::<T>(1e-10) {
            <T as NumericElement>::ZERO
        } else {
            theta / delta
        };

        let rho_old = if sigma == <T as NumericElement>::ZERO {
            <T as NumericElement>::ZERO
        } else {
            <T as NumericElement>::ONE / sigma
        };

        // Compute initial residual
        let mut r = residual(matrix, b, x);

        let alpha = if sigma == <T as NumericElement>::ZERO {
            <T as NumericElement>::ONE / self.eigenvalues_max
        } else {
            let rho_new = <T as NumericElement>::ONE / (from_f64::<T>(2.0) * sigma - rho_old);
            rho_new / rho_old
        };

        let offsets = matrix.row_ptr();
        let indices = matrix.col_indices();
        let values = matrix.values();

        // Chebyshev iterations
        for _ in 0..iterations.min(self.degree) {
            for i in 0..matrix.nrows() {
                let mut diag = <T as NumericElement>::ZERO;
                for k in offsets[i]..offsets[i + 1] {
                    if indices[k] == i {
                        diag = values[k];
                        break;
                    }
                }

                if NumericElement::abs(diag) > diagonal_epsilon() {
                    x[i] += alpha * r[i] / diag;
                }
            }
            r = residual(matrix, b, x);
        }
    }

    fn clone_box(&self) -> Box<dyn MultigridSmoother<T>> {
        Box::new(self.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::super::csr_from_parts;
    use super::*;

    fn assert_vector_eq<const N: usize>(vector: &MultigridVector<f64>, expected: [f64; N]) {
        assert_eq!(vector.shape(), [N]);
        for (idx, expected_value) in expected.into_iter().enumerate() {
            assert_eq!(vector[idx], expected_value);
        }
    }

    fn create_test_matrix() -> SparseMatrix<f64> {
        let n = 5;
        let mut row_ptr = Vec::with_capacity(n + 1);
        let mut col_indices = Vec::new();
        let mut values = Vec::new();
        row_ptr.push(0);
        for i in 0..n {
            if i > 0 {
                col_indices.push(i - 1);
                values.push(-1.0);
            }
            col_indices.push(i);
            values.push(2.0);
            if i < n - 1 {
                col_indices.push(i + 1);
                values.push(-1.0);
            }
            row_ptr.push(col_indices.len());
        }
        csr_from_parts(n, n, row_ptr, col_indices, values, "smoother test matrix").unwrap()
    }

    #[test]
    fn test_gauss_seidel_smoother() {
        let matrix = create_test_matrix();
        let b = MultigridVector::from_shape_vec([5], vec![1.0, 2.0, 3.0, 4.0, 5.0]).unwrap();
        let mut x = MultigridVector::zeros([5]);

        let smoother = GaussSeidelSmoother::new(1.0);
        smoother.apply(&matrix, &mut x, &b, 1);

        assert_vector_eq(&x, [0.5, 1.25, 2.125, 3.0625, 4.03125]);
    }

    #[test]
    fn test_jacobi_smoother() {
        let matrix = create_test_matrix();
        let b = MultigridVector::from_shape_vec([5], vec![1.0, 2.0, 3.0, 4.0, 5.0]).unwrap();
        let mut x = MultigridVector::zeros([5]);

        let smoother = JacobiSmoother::new(1.0);
        smoother.apply(&matrix, &mut x, &b, 1);

        assert_vector_eq(&x, [0.5, 1.0, 1.5, 2.0, 2.5]);
    }

    #[test]
    fn test_symmetric_gauss_seidel() {
        let matrix = create_test_matrix();
        let b = MultigridVector::from_shape_vec([5], vec![1.0, 2.0, 3.0, 4.0, 5.0]).unwrap();
        let mut x = MultigridVector::zeros([5]);

        let smoother = SymmetricGaussSeidelSmoother::new(1.0);
        smoother.apply(&matrix, &mut x, &b, 1);

        assert_vector_eq(&x, [2.291015625, 3.58203125, 4.6640625, 5.078125, 4.03125]);
    }

    #[test]
    fn test_sor_smoother() {
        let matrix = create_test_matrix();
        let b = MultigridVector::from_shape_vec([5], vec![1.0, 2.0, 3.0, 4.0, 5.0]).unwrap();
        let mut x = MultigridVector::zeros([5]);

        let smoother = SORSmoother::new(1.0);
        smoother.apply(&matrix, &mut x, &b, 1);

        assert_vector_eq(&x, [0.5, 1.25, 2.125, 3.0625, 4.03125]);
    }

    #[test]
    fn test_chebyshev_smoother() {
        let matrix = create_test_matrix();
        let b = MultigridVector::from_shape_vec([5], vec![1.0, 2.0, 3.0, 4.0, 5.0]).unwrap();
        let mut x = MultigridVector::zeros([5]);

        let (eigen_min, eigen_max) = ChebyshevSmoother::estimate_eigenvalues(&matrix);
        assert_eq!(eigen_min, 0.1);
        assert_eq!(eigen_max, 4.0);

        let smoother = ChebyshevSmoother::new(eigen_min, eigen_max, 3);
        smoother.apply(&matrix, &mut x, &b, 1);

        let theta = f64::midpoint(eigen_max, eigen_min);
        let delta = (eigen_max - eigen_min) / 2.0;
        let sigma = theta / delta;
        let rho_old = 1.0 / sigma;
        let rho_new = 1.0 / (2.0 * sigma - rho_old);
        let alpha = rho_new / rho_old;
        let expected_first = alpha * b[0] / 2.0;
        assert!((x[0] - expected_first).abs() <= f64::EPSILON * 16.0);
    }
}
