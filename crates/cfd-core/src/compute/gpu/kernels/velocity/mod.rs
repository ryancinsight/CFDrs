//! SIMPLE velocity correction and pressure-source assembly through Hephaestus.
//!
//! [`GpuVelocityKernel::correct`] applies
//! `velocity = velocity_star - (dt / density) * gradient(pressure)` while
//! [`GpuVelocityKernel::divergence_source`] evaluates
//! `(density / dt) * divergence(velocity_star)`. Both use centered differences
//! on a three-dimensional Cartesian grid and impose zero values at outer cells.

mod kernel;

pub use kernel::{GpuVelocityKernel, VelocityConfig};

#[cfg(test)]
mod tests;
