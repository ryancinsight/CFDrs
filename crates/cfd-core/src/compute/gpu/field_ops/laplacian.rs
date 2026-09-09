//! Laplacian field-operation facade.

use super::GpuFieldOps;
use crate::error::Result;
use aequitas::systems::si::quantities::Length;

impl GpuFieldOps {
    /// Compute the two-dimensional Laplacian of `field` into `result`.
    ///
    /// # Errors
    /// Returns a typed configuration, dimension, or provider error.
    pub fn laplacian_2d(
        &self,
        field: &[f32],
        nx: usize,
        ny: usize,
        dx: Length<f32>,
        dy: Length<f32>,
        result: &mut [f32],
    ) -> Result<()> {
        self.laplacian_kernel.execute(field, nx, ny, dx, dy, result)
    }
}
