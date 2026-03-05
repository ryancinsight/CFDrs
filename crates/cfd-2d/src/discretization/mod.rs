//! Discretization schemes for 2D CFD simulations
//!
//! This module contains various discretization schemes for convection and diffusion terms.
//!
//! # Theorem (Truncation Error Order — Taylor Expansion)
//!
//! For a smooth solution $\phi$, the central difference approximation
//! $\phi''(x) \approx (\phi_{i+1} - 2\phi_i + \phi_{i-1}) / \Delta x^2$
//! has truncation error $O(\Delta x^2)$, while first-order upwind
//! $\phi'(x) \approx (\phi_i - \phi_{i-1}) / \Delta x$ is $O(\Delta x)$.
//!
//! **Proof sketch**:
//! Taylor-expanding $\phi_{i\pm 1} = \phi_i \pm \phi'\Delta x + \phi''\Delta x^2/2
//! \pm \phi'''\Delta x^3/6 + \ldots$ and substituting into each stencil, the leading
//! cancelled terms determine the order. For central differences on the second
//! derivative, the $O(\Delta x)$ and $O(\Delta x^3)$ terms cancel by symmetry,
//! leaving $O(\Delta x^2)$.

pub mod convection;
pub mod extended_stencil;

// Re-export main discretization types
pub use convection::{
    CentralDifference, ConvectionScheme, ConvectionSchemeFactory, FirstOrderUpwind, HybridScheme,
    PowerLawScheme, QuadraticUpstreamInterpolationScheme,
};
