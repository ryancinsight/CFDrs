//! Laplacian operator GPU implementation
//!
//! # Mathematical Foundation
//!
//! Discrete Laplacian ∇²u for 2D scalar fields on a uniform Cartesian grid using the
//! second-order 5-point finite difference stencil, with rigorously defined boundary stencils.
//!
//! ## Interior Stencil (Second Order)
//!
//! For u(x,y) sampled on indices (i,j) with spacings (Δx, Δy):
//!
//! ```text
//! ∇²u(i,j) ≈ (u_{i-1,j} - 2u_{i,j} + u_{i+1,j})/Δx² + (u_{i,j-1} - 2u_{i,j} + u_{i,j+1})/Δy²
//! ```
//!
//! Truncation error is O(Δx² + Δy²) via Taylor expansion (LeVeque 2007; Strikwerda 2004).
//!
//! ## Boundary Conditions (Endpoint-Inclusive Grids)
//!
//! Let the physical endpoints be included (`i=0` and `i=nx-1` for X; `j=0` and `j=ny-1` for Y).
//! Boundary stencils are chosen to preserve second-order accuracy and avoid bias.
//!
//! - Dirichlet (u=0): ghost points via odd reflection.
//!   - X-left (`i=0`): `u_{-1,j} = -u_{1,j}` ⇒ `d²u/dx²|_{i=0} = (-2 u_{0,j})/Δx²`
//!   - X-right (`i=nx-1`): `u_{nx,j} = -u_{nx-2,j}` ⇒ `d²u/dx²|_{i=nx-1} = (-2 u_{nx-1,j})/Δx²`
//!   - Y-bottom/top analogous in Y.
//!
//! - Neumann (∂u/∂n = 0): second-order one-sided second derivatives (Trefethen 1996-style one-sided stencils).
//!   - X-left (`i=0`): `d²u/dx² ≈ (2u₀ − 5u₁ + 4u₂ − u₃)/Δx²` (requires `nx ≥ 4`)
//!   - X-right (`i=nx-1`): `d²u/dx² ≈ (2u₀ − 5u₁ + 4u₂ − u₃)/Δx²` with `u₁=u_{nx-2}`, `u₂=u_{nx-3}`, `u₃=u_{nx-4}`
//!   - Y-bottom/top analogous in Y with `ny ≥ 4`.
//!   - Fallback for small axes (`< 4`) uses symmetric ghosting to avoid out-of-bounds, maintaining robustness.
//!
//! - Periodic: endpoint-inclusive wrapping to inner indices for second-order interior stencil consistency.
//!   - X-left wraps left neighbor to `nx-2`; X-right wraps right neighbor to `1`.
//!   - Y-bottom wraps bottom neighbor to `ny-2`; Y-top wraps top neighbor to `1`.
//!
//! These boundary treatments align CPU and GPU implementations exactly and prevent first-order bias at endpoints.
//!
//! ## Convergence, Stability, Spectral Properties
//!
//! - Consistency: O(Δx² + Δy²) convergence for smooth solutions.
//! - Stability: discrete Laplacian is negative semi-definite on uniform grids; eigenvalues satisfy bounds consistent with
//!   classical results (see Strikwerda 2004). Together, consistency + stability ⇒ convergence (Lax Equivalence Theorem).
//! - Spectral (Dirichlet, n×n): `λ_{k,l} = -4/h² [sin²(πk/2(n+1)) + sin²(πl/2(n+1))]`, `k,l=1..n`.
//!
//! ## Literature References
//!
//! - LeVeque, R. J. (2007). *Finite Difference Methods for Ordinary and Partial Differential Equations*. SIAM.
//! - Strikwerda, J. C. (2004). *Finite Difference Schemes and Partial Differential Equations* (2nd ed.). SIAM.
//! - Trefethen, L. N. (1996). *Finite Difference and Spectral Methods for ODEs/PDEs* (notes), stencils and accuracy analysis.

mod cpu_reference;
mod kernel;
mod types;

#[cfg(test)]
mod tests;

pub use kernel::Laplacian2DKernel;
pub use types::BoundaryType;
