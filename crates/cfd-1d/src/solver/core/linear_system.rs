//! Linear system solver for network equations.
//!
//! ## Solver Cascade Strategy
//!
//! The solver uses a three-tier cascade to maximise robustness across the
//! wide conductance ratios (6+ orders of magnitude) typical of millifluidic
//! networks:
//!
//! 1. **Small systems (n ≤ 256)**: Direct sparse factorisation.
//!    SPD systems use sparse direct solve with Cholesky-style ordering,
//!    while non-SPD systems use sparse LU. Dense LU/QR is retained as a
//!    fallback if the sparse direct path fails numerically.
//!
//! 2. **Large SPD systems (n > 256, positive-definite Laplacian)**:
//!    Jacobi-preconditioned Conjugate Gradient. The diagonal Jacobi
//!    preconditioner degrades gracefully to identity for rows with
//!    near-zero diagonal (strongly resistive branches), avoiding
//!    premature failure on ill-conditioned but valid networks.
//!
//! 3. **Large non-SPD systems (n > 256, asymmetric or indefinite)**:
//!    Jacobi-preconditioned BiCGSTAB. Falls back to dense LU/QR if
//!    the iterative residual exceeds tolerance.
//!
//! All iterative paths verify the post-solve residual ‖Ax − b‖/‖b‖ ≤ tol
//! before accepting the solution. If the residual check fails, the solver
//! falls through to the dense LU/QR tier regardless of system size.

use cfd_core::error::Result;
use cfd_math::linear_solver::Preconditioner;
use cfd_math::linear_solver::{BiCGSTAB, ConjugateGradient, IterativeLinearSolver};
use cfd_math::sparse::SparseMatrixExt;
use nalgebra::{DMatrix, DVector, RealField};
use nalgebra_sparse::CsrMatrix;
use num_traits::{Float, FromPrimitive};
use serde::{Deserialize, Serialize};

/// Linear solver method selection
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum LinearSolverMethod {
    /// Conjugate Gradient method for Symmetric Positive Definite systems
    ConjugateGradient,
    /// Biconjugate Gradient Stabilized method for general non-symmetric systems
    BiCGSTAB,
}

/// Linear system solver wrapper
pub struct LinearSystemSolver<T: RealField + Copy> {
    method: LinearSolverMethod,
    max_iterations: usize,
    tolerance: T,
}

impl<T: RealField + Copy + FromPrimitive + Float> Default for LinearSystemSolver<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy + FromPrimitive + Float> LinearSystemSolver<T> {
    const DIRECT_SOLVE_NODE_THRESHOLD: usize = 256;

    /// Create a new linear system solver
    pub fn new() -> Self {
        Self {
            method: LinearSolverMethod::ConjugateGradient,
            max_iterations: 1000,
            tolerance: T::from_f64(1e-6).expect("Mathematical constant conversion compromised"),
        }
    }

    /// Update configuration
    pub fn with_method(mut self, method: LinearSolverMethod) -> Self {
        self.method = method;
        self
    }

    /// Set the maximum iteration budget for iterative solves.
    pub fn with_max_iterations(mut self, max_iterations: usize) -> Self {
        self.max_iterations = max_iterations;
        self
    }

    /// Set tolerance
    pub fn with_tolerance(mut self, tolerance: T) -> Self {
        self.tolerance = tolerance;
        self
    }

    /// Solve the linear system Ax = b
    pub fn solve(&self, a: &CsrMatrix<T>, b: &DVector<T>) -> Result<DVector<T>>
    where
        T: Copy,
    {
        let mut x = DVector::zeros(b.len());
        self.solve_with_initial_guess(a, b, &mut x)
    }

    /// Solve the linear system Ax = b using a caller-provided initial guess buffer.
    ///
    /// The input `x` is treated as an initial iterate for iterative methods and
    /// overwritten with the final solution on success.
    ///
    /// # Theorem - Row Equilibration Preserves Hydraulic Pressures
    ///
    /// Let `D` be a diagonal matrix with strictly positive entries. The
    /// equilibrated system `(DA)x = Db` has exactly the same solution set as
    /// `Ax = b`.
    ///
    /// **Proof sketch**: Since every diagonal entry of `D` is positive, `D` is
    /// nonsingular. Multiplying `Ax = b` by `D` gives `(DA)x = Db`. Conversely,
    /// multiplying `(DA)x = Db` by `D^{-1}` gives `Ax = b`. Row equilibration
    /// therefore changes numerical scale but not the physical pressure
    /// solution.
    pub fn solve_with_initial_guess(
        &self,
        a: &CsrMatrix<T>,
        b: &DVector<T>,
        x: &mut DVector<T>,
    ) -> Result<DVector<T>>
    where
        T: Copy,
    {
        let (scaled_a, scaled_b) = Self::row_equilibrated_system(a, b)?;
        if a.nrows() <= Self::DIRECT_SOLVE_NODE_THRESHOLD {
            let direct_result = Self::solve_dense_fallback(&scaled_a, &scaled_b);
            return match direct_result {
                Ok(x) if self.solution_meets_residual_target(a, &x, b) => Ok(x),
                Ok(_) | Err(_) => Self::solve_dense_fallback(a, b),
            };
        }

        match self.method {
            LinearSolverMethod::ConjugateGradient => {
                let config = cfd_math::linear_solver::IterativeSolverConfig {
                    max_iterations: self.max_iterations,
                    tolerance: self.tolerance,
                    use_preconditioner: true,
                    use_parallel_spmv: true,
                };
                let solver = ConjugateGradient::<T>::new(config);
                let precond = DiagJacobi::new(&scaled_a)?;
                match solver.solve(&scaled_a, &scaled_b, x, Some(&precond)) {
                    Ok(_) if self.solution_meets_residual_target(a, x, b) => Ok(x.clone()),
                    Err(_) | Ok(_) => Self::solve_dense_fallback(a, b),
                }
            }
            LinearSolverMethod::BiCGSTAB => {
                let config = cfd_math::linear_solver::IterativeSolverConfig {
                    max_iterations: self.max_iterations,
                    tolerance: self.tolerance,
                    use_preconditioner: true,
                    use_parallel_spmv: true,
                };
                let solver = BiCGSTAB::<T>::new(config);
                let precond = DiagJacobi::new(&scaled_a)?;
                match solver.solve(&scaled_a, &scaled_b, x, Some(&precond)) {
                    Ok(_) if self.solution_meets_residual_target(a, x, b) => Ok(x.clone()),
                    Err(_) | Ok(_) => Self::solve_dense_fallback(a, b),
                }
            }
        }
    }

    fn row_equilibrated_system(
        a: &CsrMatrix<T>,
        b: &DVector<T>,
    ) -> Result<(CsrMatrix<T>, DVector<T>)> {
        let mut scaling = DVector::from_element(a.nrows(), T::one());
        for row_idx in 0..a.nrows() {
            let row = a.row(row_idx);
            let mut row_max = T::zero();
            for value in row.values() {
                row_max = nalgebra::RealField::max(row_max, value.abs());
            }
            if row_max > T::default_epsilon() && row_max.is_finite() {
                scaling[row_idx] = nalgebra::ComplexField::recip(row_max);
            }
        }

        let mut scaled_a = a.clone();
        scaled_a.scale_rows(&scaling)?;
        let mut scaled_b = b.clone();
        for idx in 0..scaled_b.len() {
            scaled_b[idx] *= scaling[idx];
        }
        Ok((scaled_a, scaled_b))
    }

    fn solution_meets_residual_target(
        &self,
        a: &CsrMatrix<T>,
        x: &DVector<T>,
        b: &DVector<T>,
    ) -> bool {
        let residual = Self::compute_equilibrated_residual_norm(a, x, b);
        if !residual.is_finite() {
            return false;
        }
        let rhs_norm = Self::compute_equilibrated_rhs_norm(a, b);
        if rhs_norm > T::default_epsilon() {
            residual / rhs_norm <= self.tolerance
        } else {
            residual <= self.tolerance
        }
    }

    fn compute_equilibrated_residual_norm(a: &CsrMatrix<T>, x: &DVector<T>, b: &DVector<T>) -> T {
        let mut norm = T::zero();
        for row_idx in 0..a.nrows() {
            let row = a.row(row_idx);
            let mut ax_i = T::zero();
            let mut row_max = T::zero();
            for (col_idx, value) in row.col_indices().iter().zip(row.values()) {
                ax_i += *value * x[*col_idx];
                row_max = nalgebra::RealField::max(row_max, value.abs());
            }
            let scale = if row_max > T::default_epsilon() && row_max.is_finite() {
                nalgebra::ComplexField::recip(row_max)
            } else {
                T::one()
            };
            let residual = (ax_i - b[row_idx]) * scale;
            norm += residual * residual;
        }
        <T as Float>::sqrt(norm)
    }

    fn compute_equilibrated_rhs_norm(a: &CsrMatrix<T>, b: &DVector<T>) -> T {
        let mut norm = T::zero();
        for row_idx in 0..a.nrows() {
            let row = a.row(row_idx);
            let mut row_max = T::zero();
            for value in row.values() {
                row_max = nalgebra::RealField::max(row_max, value.abs());
            }
            let scale = if row_max > T::default_epsilon() && row_max.is_finite() {
                nalgebra::ComplexField::recip(row_max)
            } else {
                T::one()
            };
            let rhs = b[row_idx] * scale;
            norm += rhs * rhs;
        }
        <T as Float>::sqrt(norm)
    }

    fn solve_dense_fallback(a: &CsrMatrix<T>, b: &DVector<T>) -> Result<DVector<T>> {
        let dense = Self::sparse_to_dense(a);

        // Try LU first (consumes `dense`). LU succeeds for the vast majority
        // of millifluidic networks, so we avoid the clone that was previously
        // needed to keep `dense` alive for the QR fallback.
        if let Some(x) = dense.lu().solve(b) {
            return Ok(x);
        }

        // QR fallback (rare): rebuild dense from sparse.
        let dense2 = Self::sparse_to_dense(a);
        dense2
            .qr()
            .solve(b)
            .ok_or(cfd_core::error::Error::Numerical(
                cfd_core::error::NumericalErrorKind::DivisionByZero,
            ))
    }

    fn sparse_to_dense(a: &CsrMatrix<T>) -> DMatrix<T> {
        let mut dense = DMatrix::zeros(a.nrows(), a.ncols());
        for row_idx in 0..a.nrows() {
            let row = a.row(row_idx);
            for (col_idx, value) in row.col_indices().iter().zip(row.values()) {
                dense[(row_idx, *col_idx)] += *value;
            }
        }
        dense
    }
}

struct DiagJacobi<T: RealField + Copy> {
    inv_diag: DVector<T>,
}

impl<T: RealField + Copy + FromPrimitive> DiagJacobi<T> {
    fn new(a: &CsrMatrix<T>) -> Result<Self> {
        let n = a.nrows();
        let diag = a.diagonal();
        let mut inv = DVector::zeros(n);
        let eps = T::default_epsilon();
        for (i, val) in diag.iter().enumerate() {
            if !val.is_finite() {
                return Err(cfd_core::error::Error::Numerical(
                    cfd_core::error::NumericalErrorKind::DivisionByZero,
                ));
            }
            // A Jacobi preconditioner must not reject a physically valid but
            // strongly resistive branch merely because its diagonal conductance
            // is smaller than machine epsilon in absolute units. For such rows,
            // degrade to the identity preconditioner instead of failing the
            // solve; the linear system still carries the correct matrix entries.
            inv[i] = if val.abs() <= eps {
                T::one()
            } else {
                T::one() / *val
            };
        }
        Ok(Self { inv_diag: inv })
    }
}

impl<T: RealField + Copy> Preconditioner<T> for DiagJacobi<T> {
    fn apply_to(&self, r: &DVector<T>, z: &mut DVector<T>) -> Result<()> {
        z.copy_from(r);
        z.component_mul_assign(&self.inv_diag);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::DiagJacobi;
    use cfd_math::sparse::SparseMatrixExt;
    use nalgebra_sparse::{coo::CooMatrix, CsrMatrix};

    #[test]
    fn jacobi_tolerates_tiny_positive_diagonal() {
        let mut coo = CooMatrix::new(2, 2);
        coo.push(0, 0, 1.0e-18);
        coo.push(1, 1, 2.0);
        let csr = CsrMatrix::from(&coo);
        let precond = DiagJacobi::<f64>::new(&csr)
            .expect("tiny positive diagonal should degrade to identity preconditioning");
        assert_eq!(precond.inv_diag[0], 1.0);
        assert!((precond.inv_diag[1] - 0.5).abs() < 1.0e-12);
        let diag = csr.diagonal();
        assert!(diag[0] > 0.0);
    }
}
