//! Linear system solver for network equations

use nalgebra::{RealField, DVector};
use nalgebra_sparse::CsrMatrix;
use cfd_core::Result;
use cfd_math::linear_solver::{LinearSolver, ConjugateGradient, BiCGSTAB};
use num_traits::FromPrimitive;

/// Linear solver method selection
#[derive(Debug, Clone, Copy)]
pub enum LinearSolverMethod {
    ConjugateGradient,
    BiCGSTAB,
}

/// Linear system solver wrapper
pub struct LinearSystemSolver<T: RealField> {
    method: LinearSolverMethod,
    max_iterations: usize,
    tolerance: T,
}

impl<T: RealField + FromPrimitive + Copy> LinearSystemSolver<T> {
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
                let config = cfd_math::linear_solver::LinearSolverConfig {
                    max_iterations: self.max_iterations,
                    tolerance: self.tolerance,
                    preconditioning: false,
                };
                let solver = ConjugateGradient::<T>::new(config);
                solver.solve(a, b, Some(&x0))
            }
            LinearSolverMethod::BiCGSTAB => {
                let config = cfd_math::linear_solver::LinearSolverConfig {
                    max_iterations: self.max_iterations,
                    tolerance: self.tolerance,
                    preconditioning: false,
                };
                let solver = BiCGSTAB::<T>::new(config);
                solver.solve(a, b, Some(&x0))
            }
        }
    }
}