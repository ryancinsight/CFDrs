//! Laplacian field-operation facade.

use super::GpuFieldOps;

impl GpuFieldOps {
    /// Compute the two-dimensional Laplacian of `field` into `result`.
    pub fn laplacian_2d(
        &self,
        field: &[f32],
        nx: usize,
        ny: usize,
        dx: f32,
        dy: f32,
        result: &mut [f32],
    ) {
        self.laplacian_kernel.execute(field, nx, ny, dx, dy, result);
    }
}
