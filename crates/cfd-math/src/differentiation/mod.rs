//! Numerical differentiation methods for CFD applications.
//!
//! This module provides finite difference schemes optimized for CFD simulations
//! with support for various boundary conditions and grid types.

mod finite_difference;
mod gradient;
mod schemes;
pub use finite_difference::{differentiate_1d, differentiate_2d, FiniteDifference};
pub use gradient::{compute_gradient_2d, compute_gradient_3d, Gradient};
pub use schemes::FiniteDifferenceScheme;
// Re-export commonly used items at module level
pub use finite_difference::laplacian_2d;
#[cfg(test)]
mod tests;
