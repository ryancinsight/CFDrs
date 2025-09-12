//! Linear system solvers for algebraic equations

use super::traits::LinearSystemSolver;
use nalgebra::{DMatrix, DVector, RealField};

/// Tolerance for detecting breakdown in iterative solvers (machine epsilon level)
const BREAKDOWN_TOLERANCE: f64 = 1e-14;
/// Default maximum iterations for iterative solvers
const DEFAULT_MAX_ITERATIONS: usize = 1000;
/// Default convergence tolerance for iterative methods
const DEFAULT_TOLERANCE: f64 = 1e-6;

/// Conjugate Gradient solver for symmetric positive definite systems
#[derive(Debug, Clone)]
pub struct ConjugateGradient<T: RealField + Copy> {
    /// Maximum iterations
    pub max_iterations: usize,
    /// Convergence tolerance
    pub tolerance: T,
}

impl<T: RealField + Copy> LinearSystemSolver<T> for ConjugateGradient<T> {
    fn solve(&self, matrix: &DMatrix<T>, rhs: &DVector<T>) -> Option<DVector<T>> {
        // Conjugate Gradient implementation with preconditioning
        // Reference: Saad, Y. "Iterative Methods for Sparse Linear Systems" (2003)

        let n = rhs.len();
        if matrix.nrows() != n || matrix.ncols() != n {
            return None;
        }

        // Initialize solution vector
        let mut x = DVector::zeros(n);

        // Initial residual r = b - Ax
        let mut r = rhs.clone();
        let mut p = r.clone();
        let mut rsold = r.dot(&r);

        // Check for zero right-hand side
        if rhs.norm() < self.tolerance {
            return Some(x); // Solution is zero vector
        }

        for _iter in 0..self.max_iterations {
            // Matrix-vector product: Ap
            let ap = matrix * &p;

            // Step length: α = (r^T * r)/(p^T * Ap)
            let pap = p.dot(&ap);

            // Check for breakdown
            if pap.abs() < T::from_f64(BREAKDOWN_TOLERANCE).unwrap_or_else(T::zero) {
                // Try to return current solution if residual is small enough
                if r.norm() < self.tolerance {
                    return Some(x);
                }
                return None;
            }

            let alpha = rsold / pap;

            // Update solution: x = x + α * p
            x += &p * alpha;

            // Update residual: r = r - α * Ap
            r = &r - &ap * alpha;

            // Check convergence
            let rsnew = r.dot(&r);
            if rsnew.sqrt() < self.tolerance {
                return Some(x);
            }

            // Update search direction
            let beta = rsnew / rsold;
            p = &r + &p * beta;
            rsold = rsnew;
        }

        // Did not converge within max iterations
        None
    }

    fn name(&self) -> &str {
        "Conjugate Gradient"
    }

    fn is_iterative(&self) -> bool {
        true
    }
}

impl<T: RealField + Copy> Default for ConjugateGradient<T> {
    fn default() -> Self {
        Self {
            max_iterations: DEFAULT_MAX_ITERATIONS,
            tolerance: T::from_f64(DEFAULT_TOLERANCE).unwrap_or_else(|| T::one()),
        }
    }
}

/// Jacobi iterative solver
#[derive(Debug, Clone)]
pub struct Jacobi<T: RealField + Copy> {
    /// Maximum iterations
    pub max_iterations: usize,
    /// Convergence tolerance
    pub tolerance: T,
}

impl<T: RealField + Copy> LinearSystemSolver<T> for Jacobi<T> {
    fn solve(&self, matrix: &DMatrix<T>, rhs: &DVector<T>) -> Option<DVector<T>> {
        let n = rhs.len();
        if matrix.nrows() != n || matrix.ncols() != n {
            return None;
        }

        let mut x = DVector::zeros(n);
        let mut x_next = DVector::zeros(n);

        for _iter in 0..self.max_iterations {
            for i in 0..n {
                let mut sum = T::zero();
                for j in 0..n {
                    if i != j {
                        sum += matrix[(i, j)] * x[j];
                    }
                }
                let diag = matrix[(i, i)];
                if diag.abs() < T::from_f64(BREAKDOWN_TOLERANCE).unwrap_or_else(T::zero) {
                    return None; // Zero diagonal element
                }
                x_next[i] = (rhs[i] - sum) / diag;
            }

            // Check convergence
            let residual = matrix * &x_next - rhs;
            if residual.norm() < self.tolerance {
                return Some(x_next);
            }

            std::mem::swap(&mut x, &mut x_next);
        }

        None
    }

    fn name(&self) -> &str {
        "Jacobi"
    }

    fn is_iterative(&self) -> bool {
        true
    }
}

impl<T: RealField + Copy> Default for Jacobi<T> {
    fn default() -> Self {
        Self {
            max_iterations: DEFAULT_MAX_ITERATIONS,
            tolerance: T::from_f64(DEFAULT_TOLERANCE).unwrap_or_else(|| T::one()),
        }
    }
}

/// Gauss-Seidel iterative solver
#[derive(Debug, Clone)]
pub struct GaussSeidel<T: RealField + Copy> {
    /// Maximum iterations
    pub max_iterations: usize,
    /// Convergence tolerance
    pub tolerance: T,
}

impl<T: RealField + Copy> LinearSystemSolver<T> for GaussSeidel<T> {
    fn solve(&self, matrix: &DMatrix<T>, rhs: &DVector<T>) -> Option<DVector<T>> {
        let n = rhs.len();
        if matrix.nrows() != n || matrix.ncols() != n {
            return None;
        }

        let mut x = DVector::zeros(n);

        for _iter in 0..self.max_iterations {
            for i in 0..n {
                let mut sum = T::zero();
                for j in 0..i {
                    sum += matrix[(i, j)] * x[j]; // Use updated values
                }
                for j in (i + 1)..n {
                    sum += matrix[(i, j)] * x[j]; // Use old values
                }
                let diag = matrix[(i, i)];
                if diag.abs() < T::from_f64(BREAKDOWN_TOLERANCE).unwrap_or_else(T::zero) {
                    return None; // Zero diagonal element
                }
                x[i] = (rhs[i] - sum) / diag;
            }

            // Check convergence
            let residual = matrix * &x - rhs;
            if residual.norm() < self.tolerance {
                return Some(x);
            }
        }

        None
    }

    fn name(&self) -> &str {
        "Gauss-Seidel"
    }

    fn is_iterative(&self) -> bool {
        true
    }
}

impl<T: RealField + Copy> Default for GaussSeidel<T> {
    fn default() -> Self {
        Self {
            max_iterations: DEFAULT_MAX_ITERATIONS,
            tolerance: T::from_f64(DEFAULT_TOLERANCE).unwrap_or_else(|| T::one()),
        }
    }
}

/// Direct solver using LU decomposition
#[derive(Debug, Clone)]
pub struct DirectSolver;

impl<T: RealField + Copy> LinearSystemSolver<T> for DirectSolver {
    fn solve(&self, matrix: &DMatrix<T>, rhs: &DVector<T>) -> Option<DVector<T>> {
        let n = rhs.len();
        if matrix.nrows() != n || matrix.ncols() != n {
            return None;
        }

        // Use nalgebra's built-in LU decomposition
        matrix.clone().lu().solve(rhs)
    }

    fn name(&self) -> &str {
        "Direct (LU)"
    }

    fn is_iterative(&self) -> bool {
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_conjugate_gradient() {
        // Create a symmetric positive definite matrix
        let matrix = DMatrix::from_row_slice(3, 3, &[4.0, 1.0, 0.0, 1.0, 3.0, 1.0, 0.0, 1.0, 2.0]);
        let rhs = DVector::from_vec(vec![1.0, 2.0, 3.0]);

        let solver = ConjugateGradient::<f64>::default();
        let solution = solver.solve(&matrix, &rhs).unwrap();

        // Verify the solution
        let residual = &matrix * &solution - &rhs;
        assert!(residual.norm() < 1e-6);
    }

    #[test]
    fn test_jacobi() {
        // Create a diagonally dominant matrix
        let matrix =
            DMatrix::from_row_slice(3, 3, &[4.0, -1.0, 0.0, -1.0, 4.0, -1.0, 0.0, -1.0, 3.0]);
        let rhs = DVector::from_vec(vec![15.0, 10.0, 10.0]);

        let solver = Jacobi::<f64>::default();
        let solution = solver.solve(&matrix, &rhs).unwrap();

        // Verify the solution
        let residual = &matrix * &solution - &rhs;
        assert!(residual.norm() < 1e-6);
    }

    #[test]
    fn test_direct_solver() {
        // Test system with known exact solution:
        // 2x + y + z = 5
        // x + 3y + 2z = 8
        // x + 2y + 3z = 9
        // Solution: x = 1, y = 1, z = 2 (verified by substitution)
        let matrix = DMatrix::from_row_slice(3, 3, &[2.0, 1.0, 1.0, 1.0, 3.0, 2.0, 1.0, 2.0, 3.0]);
        let rhs = DVector::from_vec(vec![5.0, 8.0, 9.0]);

        let solver = DirectSolver;
        let solution = solver.solve(&matrix, &rhs).unwrap();

        // Validate exact analytical solution (not just residual)
        let expected = DVector::from_vec(vec![1.0_f64, 1.0_f64, 2.0_f64]);
        assert!(
            (solution[0] - expected[0]).abs() < 1e-14_f64,
            "x: expected {}, got {}",
            expected[0],
            solution[0]
        );
        assert!(
            (solution[1] - expected[1]).abs() < 1e-14_f64,
            "y: expected {}, got {}",
            expected[1],
            solution[1]
        );
        assert!(
            (solution[2] - expected[2]).abs() < 1e-14_f64,
            "z: expected {}, got {}",
            expected[2],
            solution[2]
        );

        // Also verify residual is near machine precision
        let residual = &matrix * &solution - &rhs;
        assert!(
            residual.norm() < 1e-14_f64,
            "Residual norm {} should be near machine precision",
            residual.norm()
        );
    }
}
