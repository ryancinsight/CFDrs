//! Explicit central-difference diffusion through the Hephaestus WGPU provider.
//!
//! The kernel advances a scalar field by
//! `u' = u + dt * diffusivity * Laplacian(u)` on a three-dimensional Cartesian
//! grid. Boundary values are copied unchanged. [`DiffusionConfig`] validates
//! the grid, physical coefficients, and the forward-Euler stability condition
//! `diffusivity * dt * sum(1 / spacing_i^2) <= 1/2` before dispatch.

mod kernel;

pub use kernel::{DiffusionConfig, GpuDiffusionKernel};

#[cfg(test)]
mod tests;
