//! Numerical differentiation methods for CFD applications.
//!
//! This module provides finite difference schemes optimized for CFD simulations
//! with support for various boundary conditions and grid types.
//!
//! ## Theorem — Taylor Remainder Bound (Finite Difference Truncation)
//!
//! **Theorem**: For a function f ∈ C^{p+1}([a,b]) and a p-th order finite difference
//! stencil Dₕf(x), the truncation error satisfies:
//!
//! ```text
//! |Dₕf(x) - f'(x)| ≤ C · hᵖ · max_{ξ∈[a,b]} |f^{(p+1)}(ξ)|
//! ```
//!
//! where h is the grid spacing and C depends only on the stencil coefficients.
//!
//! **Corollary — Central Difference Superiority**: The 2nd-order central difference
//! `(f(x+h) - f(x-h)) / (2h)` eliminates the O(h) term via cancellation by symmetry,
//! giving O(h²) accuracy at no extra function evaluations compared to 1st-order upwind.
//!
//! ## Theorem — Spectral Differentiation Accuracy
//!
//! **Theorem (Spectral Convergence)**: For analytic periodic functions, the spectral
//! derivative computed via DFT achieves exponential convergence:
//!
//! ```text
//! ‖Dₛf - f'‖_∞ = O(e^{-c·N})
//! ```
//!
//! where N is the number of modes and c depends on the width of the analyticity strip.
//! This is the foundation of DNS (Direct Numerical Simulation) accuracy.
//!
//! ## Theorem — Discrete Conservation of Divergence
//!
//! **Theorem (Summation-by-Parts)**: The SBP finite difference operator D satisfies:
//!
//! ```text
//! vᵀ · H · D · u = -uᵀ · H · D · v + [boundary terms]
//! ```
//!
//! where H is a positive-definite diagonal norm matrix, enabling discrete integration
//! by parts and guaranteeing energy stability of the discretized PDE.
//!
//! ## Invariants
//! - All 1D stencil coefficient arrays sum to 0 (preservation of constants).
//! - Grid spacing h must be strictly positive; enforced at construction.
//! - Spectral modes N must be a power of 2 for FFT efficiency.

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
