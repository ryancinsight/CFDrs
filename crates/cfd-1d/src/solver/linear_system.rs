//! Linear system solver for network equations

use nalgebra::{RealField, DVector};
use nalgebra_sparse::CsrMatrix;
use cfd_core::{Result, solver::LinearSolverConfig};
use cfd_math::linear_solver::{LinearSolver, ConjugateGradient, BiCGSTAB};
use num_traits::FromPrimitive;

/// Linear system solver wrapper
pub struct LinearSystemSolver<T: RealField> {
    config: LinearSolverConfig,
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField> LinearSystemSolver<T> {
    /// Create a new linear system solver
    pub fn new(config: LinearSolverConfig) -> Self {
        Self {
            config,
            _phantom: std::marker::PhantomData,
        }
    }

    /// Update solver configuration
    pub fn update_config(&mut self, config: &LinearSolverConfig) -> Result<()> {
        self.config = config.clone();
        Ok(())
    }
}

impl<T: RealField + FromPrimitive + std::fmt::Debug> LinearSystemSolver<T> {
    /// Solve the linear system Ax = b
    pub fn solve(&self, a: CsrMatrix<T>, b: DVector<T>) -> Result<DVector<T>> {
        let x0 = DVector::zeros(b.len());
        
        match self.config.method {
            cfd_core::solver::LinearSolverMethod::ConjugateGradient => {
                let solver = ConjugateGradient::<T>::new(
                    self.config.max_iterations,
                    T::from_f64(self.config.tolerance).unwrap_or_else(T::zero),
                );
                solver.solve(&a, &b, Some(&x0))
            }
            cfd_core::solver::LinearSolverMethod::BiCGSTAB => {
                let solver = BiCGSTAB::<T>::new(
                    self.config.max_iterations,
                    T::from_f64(self.config.tolerance).unwrap_or_else(T::zero),
                );
                solver.solve(&a, &b, Some(&x0))
            }
            _ => {
                // Fallback to CG for other methods
                let solver = ConjugateGradient::<T>::new(
                    self.config.max_iterations,
                    T::from_f64(self.config.tolerance).unwrap_or_else(T::zero),
                );
                solver.solve(&a, &b, Some(&x0))
            }
        }
    }
}