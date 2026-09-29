//! Geometric Multigrid (GMG) for structured grids
//!
//! Geometric multigrid exploits structured grid hierarchies for optimal
//! convergence properties and computational efficiency.
//!
//! ## Mathematical Foundation
//!
//! For a structured grid with mesh size h, the geometric multigrid hierarchy is:
//! ```math
//! Ω¹ ⊂ Ω² ⊂ ⋯ ⊂ Ω^J
//! ```
//!
//! where Ω¹ is the finest grid and Ω^J is the coarsest grid.
//!
//! ## Algorithm Overview
//!
//! 1. **Relaxation**: Apply smoother on fine grid
//! 2. **Restriction**: Transfer residual to coarse grid
//! 3. **Coarsest Solve**: Direct or iterative solution on coarsest grid
//! 4. **Prolongation**: Transfer correction back to fine grid
//! 5. **Post-relaxation**: Apply smoother on fine grid
//!
//! ## Literature Compliance
//!
//! - Briggs, W. L., et al. (2000). *A multigrid tutorial*. SIAM. Chapter 3.
//! - Trottenberg, U., et al. (2001). *Multigrid*. Academic Press. Chapter 4.
//! - Wesseling, P. (1992). *An introduction to multigrid methods*. Wiley.

mod multigrid;
mod ops;
mod transfer;

pub use multigrid::GeometricMultigrid;
use ops::vector_sub;
#[cfg(test)]
use ops::{from_usize, l2_norm, matrix_vector_product};

use crate::error::Result;
use eunomia::RealField;
use leto::{Array1, Array2};

type GmgVector<T> = Array1<T>;
type GmgMatrix<T> = Array2<T>;

/// Trait for nonlinear operators in multigrid methods
///
/// This trait defines the interface for nonlinear operators that can be used
/// with the Full Approximation Scheme (FAS) multigrid method.
pub trait NonlinearOperator<T: RealField> {
    /// Compute the nonlinear residual: r = f - F(u)
    ///
    /// The default implementation computes f - apply(u).
    fn residual(&self, u: &GmgVector<T>, rhs: &GmgVector<T>, level: usize) -> GmgVector<T> {
        vector_sub(rhs, &self.apply(u, level))
    }

    /// Apply the nonlinear operator: F(u)
    fn apply(&self, u: &GmgVector<T>, level: usize) -> GmgVector<T>;

    /// Solve the nonlinear system on the coarsest level
    fn coarsest_solve(&self, u: &mut GmgVector<T>, rhs: &GmgVector<T>, level: usize) -> Result<()>;

    /// Restrict the residual to the next coarser level
    fn restrict_residual(&self, fine: &GmgVector<T>, level: usize) -> GmgVector<T>;

    /// Restrict the solution to the next coarser level
    fn restrict_solution(&self, fine: &GmgVector<T>, level: usize) -> GmgVector<T>;

    /// Prolongate a vector to the next finer level
    fn prolongate(&self, coarse: &GmgVector<T>, level: usize) -> GmgVector<T>;

    /// Apply smoothing to the solution u given RHS f
    fn smooth(&self, u: &mut GmgVector<T>, rhs: &GmgVector<T>, level: usize, iterations: usize);
}

#[cfg(test)]
mod tests;
