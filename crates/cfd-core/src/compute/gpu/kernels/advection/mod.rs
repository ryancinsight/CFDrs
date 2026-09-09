//! First-order upwind advection through the Hephaestus WGPU provider.
//!
//! The kernel advances a scalar field independently on each z-plane:
//! `u' = u - dt (v_x D_x^upwind u + v_y D_y^upwind u)`.
//! X/Y boundary values are copied unchanged. [`AdvectionConfig`] validates the
//! grid and timestep; dispatch additionally enforces the unsplit stability
//! condition `dt (|v_x|/dx + |v_y|/dy) <= 1` at every voxel.

mod kernel;

pub use kernel::{AdvectionConfig, GpuAdvectionKernel};

#[cfg(test)]
mod tests;
