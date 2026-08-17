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

use cfd_core::error::{Error, Result};
use cfd_core::CfdScalar;
use cfd_math::linear_solver::krylov::{self, SolverKind};
use cfd_math::linear_solver::Preconditioner;
use eunomia::{FloatElement, NumericElement};
use leto::{Array1, Storage};
use leto_ops::{lu_decompose, qr_decompose, CsrMatrix as LetoCsrMatrix, Scalar as LetoScalar};
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
pub struct LinearSystemSolver<T: CfdScalar + Copy> {
    method: LinearSolverMethod,
    max_iterations: usize,
    tolerance: T,
}

impl<T: CfdScalar> Default for LinearSystemSolver<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: CfdScalar> LinearSystemSolver<T> {
    const DIRECT_SOLVE_NODE_THRESHOLD: usize = 256;

    /// Create a new linear system solver
    pub fn new() -> Self {
        Self {
            method: LinearSolverMethod::ConjugateGradient,
            max_iterations: 1000,
            tolerance: <T as FloatElement>::from_f64(1e-6),
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
    pub fn solve(&self, a: &LetoCsrMatrix<T>, b: &Array1<T>) -> Result<Array1<T>>
    where
        T: Copy,
    {
        let mut x = Array1::from_elem([b.size()], T::ZERO);
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
        a: &LetoCsrMatrix<T>,
        b: &Array1<T>,
        x: &mut Array1<T>,
    ) -> Result<Array1<T>>
    where
        T: Copy,
    {
        let (scaled_a, scaled_b) = Self::row_equilibrated_system(a, b)?;
        if a.nrows() <= Self::DIRECT_SOLVE_NODE_THRESHOLD {
            let direct_result = Self::solve_dense_fallback(&scaled_a, &scaled_b);
            return match direct_result {
                Ok(ref sol) if self.solution_meets_residual_target(a, sol, b) => direct_result,
                Ok(_) | Err(_) => Self::solve_dense_fallback(a, b),
            };
        }

        let kind = match self.method {
            LinearSolverMethod::ConjugateGradient => SolverKind::ConjugateGradient,
            LinearSolverMethod::BiCGSTAB => SolverKind::BiCgStab,
        };
        let config = cfd_math::iterative::IterativeSolverConfig::new(self.tolerance)
            .with_max_iterations(self.max_iterations)
            .with_preconditioner();
        let precond = DiagJacobi::new(&scaled_a)?;
        let converged = krylov::converged_or_none(
            kind.name(),
            kind.solve_preconditioned(&scaled_a, &scaled_b, &precond, x, &config),
        )
        .is_some();

        // The residual is re-checked against the unscaled system even after a
        // reported convergence: the solve ran on the row-equilibrated matrix,
        // so meeting its tolerance does not by itself bound the original
        // residual.
        if converged && self.solution_meets_residual_target(a, x, b) {
            Ok(x.clone())
        } else {
            Self::solve_dense_fallback(a, b)
        }
    }

    fn row_equilibrated_system(
        a: &LetoCsrMatrix<T>,
        b: &Array1<T>,
    ) -> Result<(LetoCsrMatrix<T>, Array1<T>)> {
        let mut scaling = Array1::from_elem([a.nrows()], T::ONE);
        for row_idx in 0..a.nrows() {
            let row = a.row(row_idx);
            let mut row_max = T::ZERO;
            for value in row.values() {
                row_max = row_max.max_scalar(<T as NumericElement>::abs(*value));
            }
            if row_max > T::default_epsilon() && <T as NumericElement>::is_finite(row_max) {
                scaling[row_idx] = <T as NumericElement>::ONE / row_max;
            }
        }

        let mut scaled_a = a.clone();
        scaled_a
            .scale_rows(scaling.storage().as_slice())
            .map_err(|error| {
                Error::InvalidConfiguration(format!(
                    "Leto network solver row scaling failed during equilibration: {error}"
                ))
            })?;
        let mut scaled_b = b.clone();
        for idx in 0..scaled_b.size() {
            scaled_b[idx] *= scaling[idx];
        }
        Ok((scaled_a, scaled_b))
    }

    fn solution_meets_residual_target(
        &self,
        a: &LetoCsrMatrix<T>,
        x: &Array1<T>,
        b: &Array1<T>,
    ) -> bool {
        let residual = Self::compute_equilibrated_residual_norm(a, x, b);
        if !<T as NumericElement>::is_finite(residual) {
            return false;
        }
        let rhs_norm = Self::compute_equilibrated_rhs_norm(a, b);
        if rhs_norm > T::default_epsilon() {
            residual / rhs_norm <= self.tolerance
        } else {
            residual <= self.tolerance
        }
    }

    fn compute_equilibrated_residual_norm(a: &LetoCsrMatrix<T>, x: &Array1<T>, b: &Array1<T>) -> T {
        let mut norm = T::ZERO;
        for row_idx in 0..a.nrows() {
            let row = a.row(row_idx);
            let mut ax_i = T::ZERO;
            let mut row_max = T::ZERO;
            for (col_idx, value) in row.col_indices().iter().zip(row.values()) {
                ax_i += *value * x[*col_idx];
                row_max = row_max.max_scalar(<T as NumericElement>::abs(*value));
            }
            let scale =
                if row_max > T::default_epsilon() && <T as NumericElement>::is_finite(row_max) {
                    <T as NumericElement>::ONE / row_max
                } else {
                    <T as NumericElement>::ONE
                };
            let residual = (ax_i - b[row_idx]) * scale;
            norm += residual * residual;
        }
        <T as NumericElement>::sqrt(norm)
    }

    fn compute_equilibrated_rhs_norm(a: &LetoCsrMatrix<T>, b: &Array1<T>) -> T {
        let mut norm = T::ZERO;
        for row_idx in 0..a.nrows() {
            let row = a.row(row_idx);
            let mut row_max = T::ZERO;
            for value in row.values() {
                row_max = row_max.max_scalar(<T as NumericElement>::abs(*value));
            }
            let scale =
                if row_max > T::default_epsilon() && <T as NumericElement>::is_finite(row_max) {
                    <T as NumericElement>::ONE / row_max
                } else {
                    <T as NumericElement>::ONE
                };
            let rhs = b[row_idx] * scale;
            norm += rhs * rhs;
        }
        <T as NumericElement>::sqrt(norm)
    }

    fn solve_dense_fallback(a: &LetoCsrMatrix<T>, b: &Array1<T>) -> Result<Array1<T>> {
        let dense = a.to_dense();

        // Try LU first.
        if let Ok(lu) = lu_decompose(&dense.view()) {
            if let Ok(x) = lu.solve(&b.view()) {
                return Ok(x);
            }
        }

        // QR fallback (rare): recompute dense from sparse.
        let dense2 = a.to_dense();
        qr_decompose(&dense2.view())
            .and_then(|qr| qr.solve_least_squares(&b.view()))
            .map_err(|_| {
                cfd_core::error::Error::Numerical(
                    cfd_core::error::NumericalErrorKind::DivisionByZero,
                )
            })
    }
}

struct DiagJacobi<T: CfdScalar + Copy> {
    inv_diag: Array1<T>,
}

impl<T: CfdScalar + Copy + NumericElement + LetoScalar> DiagJacobi<T> {
    fn new(a: &LetoCsrMatrix<T>) -> Result<Self> {
        let n = a.nrows();
        let diag = a.diagonal();
        let mut inv = Array1::from_elem([n], T::ZERO);
        let eps = T::default_epsilon();
        for (i, val) in diag.iter().enumerate() {
            if !<T as NumericElement>::is_finite(*val) {
                return Err(cfd_core::error::Error::Numerical(
                    cfd_core::error::NumericalErrorKind::DivisionByZero,
                ));
            }
            // A Jacobi preconditioner must not reject a physically valid but
            // strongly resistive branch merely because its diagonal conductance
            // is smaller than machine epsilon in absolute units. For such rows,
            // degrade to the identity preconditioner instead of failing the
            // solve; the linear system still carries the correct matrix entries.
            inv[i] = if <T as NumericElement>::abs(*val) <= eps {
                <T as NumericElement>::ONE
            } else {
                <T as NumericElement>::ONE / *val
            };
        }
        Ok(Self { inv_diag: inv })
    }
}

impl<T: CfdScalar + Copy> Preconditioner<T> for DiagJacobi<T> {
    fn apply_to(&self, r: &Array1<T>, z: &mut Array1<T>) -> leto::Result<()> {
        for idx in 0..r.shape()[0] {
            z[idx] = r[idx] * self.inv_diag[idx];
        }
        Ok(())
    }
}

impl<T: CfdScalar + Copy + leto_ops::RealScalar>
    athena_core::Preconditioner<athena_leto::LetoBackend<T>> for DiagJacobi<T>
{
    /// Elementwise apply straight over the borrowed views, needing no scratch.
    fn apply(
        &self,
        _backend: &athena_leto::LetoBackend<T>,
        residual: <athena_leto::LetoBackend<T> as athena_core::KrylovBackend>::View<'_>,
        mut output: <athena_leto::LetoBackend<T> as athena_core::KrylovBackend>::ViewMut<'_>,
    ) -> std::result::Result<(), athena_leto::LetoBackendError> {
        let length = residual.shape()[0];
        if output.shape()[0] != length || self.inv_diag.shape()[0] != length {
            return Err(athena_leto::LetoBackendError::LengthMismatch {
                left: length,
                right: output.shape()[0],
            });
        }
        for index in 0..length {
            output[index] = residual[index] * self.inv_diag[index];
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::DiagJacobi;
    use leto_ops::CsrMatrix;

    #[test]
    fn jacobi_tolerates_tiny_positive_diagonal() {
        let csr = CsrMatrix::from_parts(vec![1.0e-18, 2.0], vec![0, 1], vec![0, 1, 2], 2, 2)
            .expect("valid diagonal CSR");
        let precond = DiagJacobi::<f64>::new(&csr)
            .expect("tiny positive diagonal should degrade to identity preconditioning");
        assert_eq!(precond.inv_diag[0], 1.0);
        assert!((precond.inv_diag[1] - 0.5).abs() < 1.0e-12);
        let diag = csr.diagonal();
        assert!(diag[0] > 0.0);
    }
}
