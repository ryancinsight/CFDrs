//! Matrix-free linear operators for CFD applications.

/// Poisson operators for pressure and potential fields
pub mod poisson;
/// Momentum and energy operators for fluid flow
pub mod momentum;
#[cfg(feature = "gpu")]
/// GPU-accelerated operators
pub mod gpu;

pub use poisson::{LaplacianOperator2D, PoissonOperator3D};
pub use momentum::{MomentumOperator1D, MomentumOperator2D, EnergyOperator2D};
#[cfg(feature = "gpu")]
pub use gpu::{BoundaryType, GpuLaplacianOperator2D};

use crate::linear_solver::traits::LinearOperator;
use cfd_core::error::Result;
use nalgebra::{DVector, RealField};

/// Identity operator: y = x
pub struct IdentityOperator;

impl<T: RealField + Copy> LinearOperator<T> for IdentityOperator {
    fn apply(&self, x: &DVector<T>, y: &mut DVector<T>) -> Result<()> {
        y.copy_from(x);
        Ok(())
    }

    fn size(&self) -> usize {
        0 // Size is determined by input vector (0 indicates size-agnostic)
    }

    fn is_symmetric(&self) -> bool {
        true
    }
}

/// Scaled operator: y = alpha * A * x
pub struct ScaledOperator<'a, T: RealField + Copy> {
    inner: &'a dyn LinearOperator<T>,
    alpha: T,
}

impl<'a, T: RealField + Copy> ScaledOperator<'a, T> {
    /// Create a new scaled operator
    pub fn new(inner: &'a dyn LinearOperator<T>, alpha: T) -> Self {
        Self { inner, alpha }
    }
}

impl<'a, T: RealField + Copy> LinearOperator<T> for ScaledOperator<'a, T> {
    fn apply(&self, x: &DVector<T>, y: &mut DVector<T>) -> Result<()> {
        self.inner.apply(x, y)?;
        y.scale_mut(self.alpha);
        Ok(())
    }

    fn size(&self) -> usize {
        self.inner.size()
    }

    fn is_symmetric(&self) -> bool {
        self.inner.is_symmetric()
    }
}
