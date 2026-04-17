use nalgebra::{DVector, RealField};
use std::collections::VecDeque;

/// Pre-allocated workspace buffers for the non-linear network solver.
///
/// Eliminates heap allocations during the Picard iteration hot path by
/// reusing vectors, deques, and preserving constant boundary conditions.
#[derive(Debug, Clone)]
pub struct SolverWorkspace<T: RealField + Copy> {
    /// Right-hand side vector `b` in `Ax = b`
    pub rhs: DVector<T>,
    /// Solution from the previous iteration
    pub last_solution: DVector<T>,
    /// Reusable initial guess / output buffer for the linear solve.
    pub linear_solution: DVector<T>,
    /// Anderson acceleration residual histories
    pub anderson_residuals: VecDeque<DVector<T>>,
    /// Anderson acceleration iterate histories
    pub anderson_iterates: VecDeque<DVector<T>>,
    /// Pre-classified Dirichlet boundary conditions (constant)
    pub dirichlet_values: Vec<Option<T>>,
    /// Pre-classified Neumann boundary conditions (constant)
    pub neumann_sources: Vec<Option<T>>,
}

impl<T: RealField + Copy> SolverWorkspace<T> {
    /// Create a new pre-allocated workspace
    pub fn new(n: usize, anderson_depth: usize) -> Self {
        Self {
            rhs: DVector::zeros(n),
            last_solution: DVector::zeros(n),
            linear_solution: DVector::zeros(n),
            anderson_residuals: VecDeque::with_capacity(anderson_depth),
            anderson_iterates: VecDeque::with_capacity(anderson_depth),
            dirichlet_values: vec![None; n],
            neumann_sources: vec![None; n],
        }
    }

    /// Resize or clear buffers to match a given node count `n`
    pub fn resize_and_clear(&mut self, n: usize) {
        if self.rhs.len() == n {
            self.rhs.fill(T::zero());
            self.last_solution.fill(T::zero());
            self.linear_solution.fill(T::zero());
        } else {
            self.rhs = DVector::zeros(n);
            self.last_solution = DVector::zeros(n);
            self.linear_solution = DVector::zeros(n);
            self.dirichlet_values.resize(n, None);
            self.neumann_sources.resize(n, None);
        }
        self.dirichlet_values.fill(None);
        self.neumann_sources.fill(None);
        self.anderson_residuals.clear();
        self.anderson_iterates.clear();
    }
}
