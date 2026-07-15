use crate::scalar::Cfd1dScalar;
use leto::Array1;

/// Pre-allocated workspace buffers for the non-linear network solver.
///
/// Eliminates heap allocations during the Picard iteration hot path by
/// reusing vectors and preserving constant boundary conditions.
#[derive(Debug, Clone)]
pub struct SolverWorkspace<T: Cfd1dScalar + Copy> {
    /// Right-hand side vector `b` in `Ax = b`
    pub rhs: Array1<T>,
    /// Solution from the previous iteration
    pub last_solution: Array1<T>,
    /// Reusable initial guess / output buffer for the linear solve.
    pub linear_solution: Array1<T>,
    /// Pre-classified Dirichlet boundary conditions (constant)
    pub dirichlet_values: Vec<Option<T>>,
    /// Pre-classified Neumann boundary conditions (constant)
    pub neumann_sources: Vec<Option<T>>,
}

impl<T: Cfd1dScalar + Copy> SolverWorkspace<T> {
    /// Create a new pre-allocated workspace
    pub fn new(n: usize) -> Self {
        Self {
            rhs: Array1::from_elem([n], T::zero()),
            last_solution: Array1::from_elem([n], T::zero()),
            linear_solution: Array1::from_elem([n], T::zero()),
            dirichlet_values: vec![None; n],
            neumann_sources: vec![None; n],
        }
    }

    /// Resize or clear buffers to match a given node count `n`
    pub fn resize_and_clear(&mut self, n: usize) {
        if self.rhs.size() == n {
            self.rhs.fill(T::zero());
            self.last_solution.fill(T::zero());
            self.linear_solution.fill(T::zero());
        } else {
            self.rhs = Array1::from_elem([n], T::zero());
            self.last_solution = Array1::from_elem([n], T::zero());
            self.linear_solution = Array1::from_elem([n], T::zero());
            self.dirichlet_values.resize(n, None);
            self.neumann_sources.resize(n, None);
        }
        self.dirichlet_values.fill(None);
        self.neumann_sources.fill(None);
    }
}
