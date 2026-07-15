//! Weighted-Jacobi pressure iteration and residual evaluation through Hephaestus.
//!
//! The operation family discretizes `Laplacian(pressure) = source` on a
//! three-dimensional Cartesian grid. [`GpuPressureKernel::iterate`] performs
//! one weighted-Jacobi step and applies homogeneous Neumann boundary values;
//! [`GpuPressureKernel::residual`] evaluates the absolute pointwise residual.

mod kernel;

pub use kernel::{GpuPressureKernel, PressureConfig};

#[cfg(test)]
mod tests;
