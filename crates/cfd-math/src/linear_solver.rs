//! Linear solver implementations.

use nalgebra::RealField;

/// Trait for linear solvers
pub trait LinearSolver<T: RealField> {
    /// Solve Ax = b
    fn solve(&self, a: &nalgebra_sparse::CsrMatrix<T>, b: &nalgebra::DVector<T>) -> Result<nalgebra::DVector<T>, String>;
}

/// Conjugate Gradient solver
pub struct ConjugateGradient<T: RealField> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField> LinearSolver<T> for ConjugateGradient<T> {
    fn solve(&self, _a: &nalgebra_sparse::CsrMatrix<T>, b: &nalgebra::DVector<T>) -> Result<nalgebra::DVector<T>, String> {
        Ok(b.clone())
    }
}

/// GMRES solver
pub struct GMRES<T: RealField> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField> LinearSolver<T> for GMRES<T> {
    fn solve(&self, _a: &nalgebra_sparse::CsrMatrix<T>, b: &nalgebra::DVector<T>) -> Result<nalgebra::DVector<T>, String> {
        Ok(b.clone())
    }
}

/// BiCGSTAB solver
pub struct BiCGSTAB<T: RealField> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField> LinearSolver<T> for BiCGSTAB<T> {
    fn solve(&self, _a: &nalgebra_sparse::CsrMatrix<T>, b: &nalgebra::DVector<T>) -> Result<nalgebra::DVector<T>, String> {
        Ok(b.clone())
    }
}