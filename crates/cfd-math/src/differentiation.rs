//! Numerical differentiation methods.

use nalgebra::RealField;

/// Finite difference schemes
pub struct FiniteDifference<T: RealField> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField> FiniteDifference<T> {
    /// Create a new finite difference operator
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField> Default for FiniteDifference<T> {
    fn default() -> Self {
        Self::new()
    }
}

/// Gradient computation
pub struct Gradient<T: RealField> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField> Gradient<T> {
    /// Compute gradient
    pub fn compute(&self, _field: &[T]) -> Vec<nalgebra::Vector3<T>> {
        Vec::new()
    }
}