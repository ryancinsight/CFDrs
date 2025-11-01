//! Linear system solver for network equations

use cfd_core::error::Result;
use cfd_math::linear_solver::preconditioners::IdentityPreconditioner;
use cfd_math::linear_solver::{BiCGSTAB, ConjugateGradient, IterativeLinearSolver};
use nalgebra::{DVector, RealField};
use nalgebra_sparse::CsrMatrix;
use num_traits::FromPrimitive;

/// Linear solver method selection
#[derive(Debug, Clone, Copy)]
pub enum LinearSolverMethod {
    ConjugateGradient,
    BiCGSTAB,
}

/// Linear system solver wrapper
pub struct LinearSystemSolver<T: RealField + Copy> {
    method: LinearSolverMethod,
    max_iterations: usize,
    tolerance: T,
}

impl<T: RealField + Copy + FromPrimitive + Copy> Default for LinearSystemSolver<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy + FromPrimitive + Copy> LinearSystemSolver<T> {
    /// Create a new linear system solver
    pub fn new() -> Self {
        Self {
            method: LinearSolverMethod::BiCGSTAB,
            max_iterations: 1000,
            tolerance: T::from_f64(1e-6).unwrap_or_else(T::one),
        }
    }

    /// Update configuration
    pub fn with_method(mut self, method: LinearSolverMethod) -> Self {
        self.method = method;
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
        // Initial guess
        let x0 = DVector::zeros(b.len());

        match self.method {
            LinearSolverMethod::ConjugateGradient => {
                let config = cfd_math::linear_solver::IterativeSolverConfig {
                    max_iterations: self.max_iterations,
                    tolerance: self.tolerance,
                    use_preconditioner: false,
                    use_parallel_spmv: false,
                };
                let solver = ConjugateGradient::<T>::new(config);
                let mut x = x0.clone();
                solver.solve(a, b, &mut x, None::<&IdentityPreconditioner>)?;
                Ok(x)
            }
            LinearSolverMethod::BiCGSTAB => {
                let config = cfd_math::linear_solver::IterativeSolverConfig {
                    max_iterations: self.max_iterations,
                    tolerance: self.tolerance,
                    use_preconditioner: false,
                    use_parallel_spmv: false,
                };
                let solver = BiCGSTAB::<T>::new(config);
                let mut x = x0.clone();
                solver.solve(a, b, &mut x, None::<&IdentityPreconditioner>)?;
                Ok(x)
            }
        }
    }
}
