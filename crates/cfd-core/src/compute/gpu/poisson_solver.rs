//! GPU-accelerated Poisson equation solver.
//!
//! The solver owns Hephaestus WGPU provider kernels and buffers.

use super::GpuContext;
use crate::error::{Error, Result};
use hephaestus_wgpu::{
    ComputeDevice, DispatchGrid, MultiStorageKernel, WgpuBuffer, WgpuDevice,
    WgslMultiStorageKernel, WgslStorageBinding, WgslStorageBindingLayout,
};

/// GPU Poisson solver parameters.
#[repr(C)]
#[derive(Debug, Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
pub struct PoissonParams {
    /// Number of grid points in x direction.
    pub nx: u32,
    /// Number of grid points in y direction.
    pub ny: u32,
    /// Grid spacing in x direction.
    pub dx: f32,
    /// Grid spacing in y direction.
    pub dy: f32,
    /// SOR relaxation parameter.
    pub omega: f32,
}

/// GPU-accelerated Poisson equation solver.
pub struct GpuPoissonSolver {
    provider: WgpuDevice,
    jacobi_kernel: WgslMultiStorageKernel,
    red_black_kernel: WgslMultiStorageKernel,
    residual_kernel: WgslMultiStorageKernel,
    nx: usize,
    ny: usize,
    dx: f32,
    dy: f32,
}

impl GpuPoissonSolver {
    /// Create a new GPU Poisson solver from an acquired Hephaestus context.
    ///
    /// # Errors
    ///
    /// Returns an error if the Hephaestus WGSL kernels cannot be compiled for
    /// the context's provider device.
    pub fn from_context(
        context: &GpuContext,
        nx: usize,
        ny: usize,
        dx: f32,
        dy: f32,
    ) -> Result<Self> {
        let provider = context.provider().clone();
        let storage_layouts = [
            WgslStorageBindingLayout::read_only(1),
            WgslStorageBindingLayout::read_only(2),
            WgslStorageBindingLayout::read_write(3),
        ];

        let jacobi_kernel = Self::compile_kernel(
            &provider,
            "CFD Poisson Jacobi",
            "jacobi_iteration",
            &storage_layouts,
        )?;
        let red_black_kernel = Self::compile_kernel(
            &provider,
            "CFD Poisson Red-Black",
            "red_black_iteration",
            &storage_layouts,
        )?;
        let residual_kernel = Self::compile_kernel(
            &provider,
            "CFD Poisson Residual",
            "calculate_residual",
            &storage_layouts,
        )?;

        Ok(Self {
            provider,
            jacobi_kernel,
            red_black_kernel,
            residual_kernel,
            nx,
            ny,
            dx,
            dy,
        })
    }

    fn compile_kernel(
        provider: &WgpuDevice,
        label: &'static str,
        entry_point: &'static str,
        storage_layouts: &[WgslStorageBindingLayout],
    ) -> Result<WgslMultiStorageKernel> {
        WgslMultiStorageKernel::new(
            provider,
            label,
            include_str!("shaders/poisson.wgsl"),
            entry_point,
            storage_layouts,
            0,
        )
        .map_err(|error| Error::GpuCompute(format!("Hephaestus WGSL compile failed: {error}")))
    }

    fn params(&self, omega: f32) -> Result<PoissonParams> {
        let nx = u32::try_from(self.nx).map_err(|_| {
            Error::InvalidConfiguration(format!("Poisson nx {} exceeds u32::MAX", self.nx))
        })?;
        let ny = u32::try_from(self.ny).map_err(|_| {
            Error::InvalidConfiguration(format!("Poisson ny {} exceeds u32::MAX", self.ny))
        })?;
        Ok(PoissonParams {
            nx,
            ny,
            dx: self.dx,
            dy: self.dy,
            omega,
        })
    }

    fn validate_field_lengths(&self, phi: &[f32], source: &[f32]) -> Result<()> {
        let expected = self.nx.checked_mul(self.ny).ok_or_else(|| {
            Error::InvalidConfiguration(format!(
                "Poisson grid size overflows usize: nx={}, ny={}",
                self.nx, self.ny
            ))
        })?;
        if phi.len() != expected {
            return Err(Error::DimensionMismatch {
                expected,
                actual: phi.len(),
            });
        }
        if source.len() != expected {
            return Err(Error::DimensionMismatch {
                expected,
                actual: source.len(),
            });
        }
        Ok(())
    }

    fn dispatch_grid(&self) -> Result<DispatchGrid> {
        DispatchGrid::covering_domain([self.nx, self.ny, 1], [8, 8, 1])
            .map_err(|error| Error::GpuCompute(format!("Hephaestus dispatch grid failed: {error}")))
    }

    fn storage_bindings<'a>(
        input: &'a WgpuBuffer<f32>,
        source: &'a WgpuBuffer<f32>,
        output: &'a WgpuBuffer<f32>,
    ) -> [WgslStorageBinding<'a>; 3] {
        [
            WgslStorageBinding::new(1, input),
            WgslStorageBinding::new(2, source),
            WgslStorageBinding::new(3, output),
        ]
    }

    fn upload_pair(
        &self,
        phi: &[f32],
        source: &[f32],
    ) -> Result<(WgpuBuffer<f32>, WgpuBuffer<f32>, WgpuBuffer<f32>)> {
        let phi_a = self
            .provider
            .upload(phi)
            .map_err(|error| Error::GpuCompute(format!("Hephaestus phi upload failed: {error}")))?;
        let phi_b = self.provider.upload(phi).map_err(|error| {
            Error::GpuCompute(format!("Hephaestus work-buffer upload failed: {error}"))
        })?;
        let source = self.provider.upload(source).map_err(|error| {
            Error::GpuCompute(format!("Hephaestus source upload failed: {error}"))
        })?;
        Ok((phi_a, phi_b, source))
    }

    fn download_into(&self, buffer: &WgpuBuffer<f32>, out: &mut [f32]) -> Result<()> {
        self.provider
            .download(buffer, out)
            .map_err(|error| Error::GpuCompute(format!("Hephaestus download failed: {error}")))
    }

    /// Solve Poisson equation using Jacobi iteration.
    ///
    /// # Errors
    ///
    /// Returns an error if field dimensions do not match the configured grid or
    /// if Hephaestus buffer transfer/dispatch fails.
    pub fn solve_jacobi(
        &self,
        phi: &mut [f32],
        source: &[f32],
        iterations: usize,
        omega: f32,
    ) -> Result<()> {
        self.validate_field_lengths(phi, source)?;
        let params = self.params(omega)?;
        let grid = self.dispatch_grid()?;
        let (phi_a, phi_b, source_buffer) = self.upload_pair(phi, source)?;

        for iter in 0..iterations {
            let (input, output) = if iter % 2 == 0 {
                (&phi_a, &phi_b)
            } else {
                (&phi_b, &phi_a)
            };
            self.jacobi_kernel
                .dispatch(
                    &self.provider,
                    Self::storage_bindings(input, &source_buffer, output),
                    &params,
                    grid,
                )
                .map_err(|error| {
                    Error::GpuCompute(format!("Hephaestus Jacobi dispatch failed: {error}"))
                })?;
        }

        let final_buffer = if iterations.is_multiple_of(2) {
            &phi_a
        } else {
            &phi_b
        };
        self.download_into(final_buffer, phi)
    }

    /// Solve Poisson equation using Red-Black Gauss-Seidel iteration.
    ///
    /// # Errors
    ///
    /// Returns an error if field dimensions do not match the configured grid or
    /// if Hephaestus buffer transfer/dispatch fails.
    pub fn solve_red_black(
        &self,
        phi: &mut [f32],
        source: &[f32],
        iterations: usize,
        omega: f32,
    ) -> Result<()> {
        self.validate_field_lengths(phi, source)?;
        let params = self.params(omega)?;
        let grid = self.dispatch_grid()?;
        let (phi_a, phi_b, source_buffer) = self.upload_pair(phi, source)?;

        for _ in 0..iterations {
            self.red_black_kernel
                .dispatch(
                    &self.provider,
                    Self::storage_bindings(&phi_a, &source_buffer, &phi_b),
                    &params,
                    grid,
                )
                .map_err(|error| {
                    Error::GpuCompute(format!("Hephaestus red sweep dispatch failed: {error}"))
                })?;
            self.red_black_kernel
                .dispatch(
                    &self.provider,
                    Self::storage_bindings(&phi_b, &source_buffer, &phi_a),
                    &params,
                    grid,
                )
                .map_err(|error| {
                    Error::GpuCompute(format!("Hephaestus black sweep dispatch failed: {error}"))
                })?;
        }

        self.download_into(&phi_a, phi)
    }

    /// Calculate residual for convergence checking.
    ///
    /// Computes L2 norm of residual: `||laplacian(phi) - source||`.
    ///
    /// # Errors
    ///
    /// Returns an error if field dimensions do not match the configured grid or
    /// if Hephaestus buffer transfer/dispatch fails.
    pub fn calculate_residual(&self, phi: &[f32], source: &[f32]) -> Result<f32> {
        self.validate_field_lengths(phi, source)?;
        let params = self.params(1.0)?;
        let grid = self.dispatch_grid()?;
        let phi_buffer = self
            .provider
            .upload(phi)
            .map_err(|error| Error::GpuCompute(format!("Hephaestus phi upload failed: {error}")))?;
        let source_buffer = self.provider.upload(source).map_err(|error| {
            Error::GpuCompute(format!("Hephaestus source upload failed: {error}"))
        })?;
        let residual_buffer = self.provider.alloc_zeroed(phi.len()).map_err(|error| {
            Error::GpuCompute(format!("Hephaestus residual allocation failed: {error}"))
        })?;

        self.residual_kernel
            .dispatch(
                &self.provider,
                Self::storage_bindings(&phi_buffer, &source_buffer, &residual_buffer),
                &params,
                grid,
            )
            .map_err(|error| {
                Error::GpuCompute(format!("Hephaestus residual dispatch failed: {error}"))
            })?;

        let mut residuals = vec![0.0_f32; phi.len()];
        self.download_into(&residual_buffer, &mut residuals)?;
        let l2_norm_squared: f32 = residuals.iter().map(|&r| r * r).sum();
        Ok(l2_norm_squared.sqrt())
    }
}
