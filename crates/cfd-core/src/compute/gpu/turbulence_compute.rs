//! GPU-accelerated turbulence model computations
//!
//! This module provides GPU implementations for LES and DES turbulence models,
//! integrating with the CFD solver pipeline for high-performance turbulence calculations.

use super::{GpuBuffer, GpuContext};
use crate::compute::gpu::kernels::turbulence::{GpuDesKernel, GpuSmagorinskyKernel};
use crate::compute::traits::ComputeBuffer;
use crate::error::Result;
use std::sync::Arc;

/// GPU turbulence compute manager
pub struct GpuTurbulenceCompute {
    context: Arc<GpuContext>,
    smagorinsky_kernel: GpuSmagorinskyKernel<f32>,
    des_kernel: GpuDesKernel<f32>,
    // Reusable buffers
    velocity_u_buffer: Option<GpuBuffer<f32>>,
    velocity_v_buffer: Option<GpuBuffer<f32>>,
    sgs_viscosity_buffer: Option<GpuBuffer<f32>>,
    des_length_buffer: Option<GpuBuffer<f32>>,
}

impl std::fmt::Debug for GpuTurbulenceCompute {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("GpuTurbulenceCompute")
            .field("smagorinsky_kernel", &self.smagorinsky_kernel)
            .field("des_kernel", &self.des_kernel)
            .finish_non_exhaustive()
    }
}

impl GpuTurbulenceCompute {
    /// Create a new GPU turbulence compute manager
    ///
    /// # Errors
    /// Returns error if GPU context creation fails
    pub fn new() -> Result<Self> {
        let context = Arc::new(GpuContext::create()?);
        let smagorinsky_kernel = GpuSmagorinskyKernel::new();
        let des_kernel = GpuDesKernel::new();

        Ok(Self {
            context,
            smagorinsky_kernel,
            des_kernel,
            velocity_u_buffer: None,
            velocity_v_buffer: None,
            sgs_viscosity_buffer: None,
            des_length_buffer: None,
        })
    }

    /// Get reference to GPU context Arc
    #[must_use]
    pub fn context_arc(&self) -> Arc<GpuContext> {
        self.context.clone()
    }

    /// Get mutable reference to Smagorinsky kernel
    pub fn smagorinsky_kernel_mut(&mut self) -> &mut GpuSmagorinskyKernel<f32> {
        &mut self.smagorinsky_kernel
    }

    /// Get mutable reference to DES kernel
    pub fn des_kernel_mut(&mut self) -> &mut GpuDesKernel<f32> {
        &mut self.des_kernel
    }

    /// Get reference to GPU context
    #[must_use]
    pub fn context(&self) -> &GpuContext {
        &self.context
    }

    /// Compute Smagorinsky LES SGS viscosity on GPU
    ///
    /// # Arguments
    /// * `velocity_u` - U-velocity field
    /// * `velocity_v` - V-velocity field
    /// * `nx`, `ny` - Grid dimensions
    /// * `dx`, `dy` - Grid spacing
    /// * `c_s` - Smagorinsky constant
    ///
    /// # Returns
    /// SGS viscosity field as GPU buffer
    ///
    /// # Errors
    /// Returns error if GPU computation fails
    pub fn compute_smagorinsky_sgs(
        &mut self,
        velocity_u: &[f32],
        velocity_v: &[f32],
        nx: usize,
        ny: usize,
        dx: f32,
        dy: f32,
        c_s: f32,
    ) -> Result<GpuBuffer<f32>> {
        let size = nx * ny;

        // Manage velocity_u_buffer
        let velocity_u_buffer = if self
            .velocity_u_buffer
            .as_ref()
            .is_some_and(|b| b.size() == size)
        {
            let mut buffer = self.velocity_u_buffer.as_mut().unwrap().clone();
            buffer.write(velocity_u)?;
            buffer
        } else {
            let buffer = GpuBuffer::from_data(self.context.clone(), velocity_u)?;
            self.velocity_u_buffer = Some(buffer.clone());
            buffer
        };

        // Manage velocity_v_buffer
        let velocity_v_buffer = if self
            .velocity_v_buffer
            .as_ref()
            .is_some_and(|b| b.size() == size)
        {
            let mut buffer = self.velocity_v_buffer.as_mut().unwrap().clone();
            buffer.write(velocity_v)?;
            buffer
        } else {
            let buffer = GpuBuffer::from_data(self.context.clone(), velocity_v)?;
            self.velocity_v_buffer = Some(buffer.clone());
            buffer
        };

        // Manage sgs_viscosity_buffer
        let sgs_viscosity_buffer = if self
            .sgs_viscosity_buffer
            .as_ref()
            .is_some_and(|b| b.size() == size)
        {
            self.sgs_viscosity_buffer.as_ref().unwrap().clone()
        } else {
            let buffer = GpuBuffer::new(self.context.clone(), size)?;
            self.sgs_viscosity_buffer = Some(buffer.clone());
            buffer
        };

        // Execute SGS computation
        self.smagorinsky_kernel.compute_sgs_viscosity(
            &self.context.device,
            &self.context.queue,
            velocity_u_buffer.buffer(),
            velocity_v_buffer.buffer(),
            sgs_viscosity_buffer.buffer(),
            nx as u32,
            ny as u32,
            dx,
            dy,
            c_s,
        )?;

        Ok(sgs_viscosity_buffer)
    }

    /// Compute DES length scale on GPU
    ///
    /// # Arguments
    /// * `velocity_u` - U-velocity field
    /// * `velocity_v` - V-velocity field
    /// * `nx`, `ny` - Grid dimensions
    /// * `dx`, `dy` - Grid spacing
    /// * `des_constant` - DES constant (typically 0.65)
    ///
    /// # Returns
    /// DES length scale field as GPU buffer
    ///
    /// # Errors
    /// Returns error if GPU computation fails
    pub fn compute_des_length_scale(
        &mut self,
        velocity_u: &[f32],
        velocity_v: &[f32],
        nx: usize,
        ny: usize,
        dx: f32,
        dy: f32,
        des_constant: f32,
    ) -> Result<GpuBuffer<f32>> {
        let size = nx * ny;

        // Manage velocity_u_buffer
        let velocity_u_buffer = if self
            .velocity_u_buffer
            .as_ref()
            .is_some_and(|b| b.size() == size)
        {
            let mut buffer = self.velocity_u_buffer.as_mut().unwrap().clone();
            buffer.write(velocity_u)?;
            buffer
        } else {
            let buffer = GpuBuffer::from_data(self.context.clone(), velocity_u)?;
            self.velocity_u_buffer = Some(buffer.clone());
            buffer
        };

        // Manage velocity_v_buffer
        let velocity_v_buffer = if self
            .velocity_v_buffer
            .as_ref()
            .is_some_and(|b| b.size() == size)
        {
            let mut buffer = self.velocity_v_buffer.as_mut().unwrap().clone();
            buffer.write(velocity_v)?;
            buffer
        } else {
            let buffer = GpuBuffer::from_data(self.context.clone(), velocity_v)?;
            self.velocity_v_buffer = Some(buffer.clone());
            buffer
        };

        // Manage des_length_buffer
        let des_length_buffer = if self
            .des_length_buffer
            .as_ref()
            .is_some_and(|b| b.size() == size)
        {
            self.des_length_buffer.as_ref().unwrap().clone()
        } else {
            let buffer = GpuBuffer::new(self.context.clone(), size)?;
            self.des_length_buffer = Some(buffer.clone());
            buffer
        };

        // Execute DES length scale computation
        self.des_kernel.compute_des_length_scale(
            &self.context.device,
            &self.context.queue,
            velocity_u_buffer.buffer(),
            velocity_v_buffer.buffer(),
            des_length_buffer.buffer(),
            nx as u32,
            ny as u32,
            dx,
            dy,
            des_constant,
        )?;

        Ok(des_length_buffer)
    }

    /// Read back data from GPU buffer to CPU
    ///
    /// # Arguments
    /// * `buffer` - GPU buffer to read from
    ///
    /// # Returns
    /// Vector containing buffer data
    ///
    /// # Errors
    /// Returns error if GPU readback fails
    pub fn read_buffer(&self, buffer: &GpuBuffer<f32>) -> Result<Vec<f32>> {
        ComputeBuffer::read(buffer)
    }

    /// Get GPU performance information
    #[must_use]
    pub fn performance_info(&self) -> TurbulencePerformanceInfo {
        TurbulencePerformanceInfo {
            max_work_group_size: self.context.max_work_group_size(),
            max_buffer_size: self.context.max_buffer_size(),
            adapter_info: self.context.adapter_info(),
        }
    }
}

/// GPU performance information
#[derive(Debug, Clone)]
pub struct TurbulencePerformanceInfo {
    /// Maximum work group size
    pub max_work_group_size: usize,
    /// Maximum buffer size
    pub max_buffer_size: usize,
    /// Adapter information
    pub adapter_info: wgpu::AdapterInfo,
}

impl TurbulencePerformanceInfo {
    /// Estimate performance scaling for turbulence computations
    #[must_use]
    pub fn estimated_speedup(&self, problem_size: usize) -> f64 {
        let gpu_factor = if self.adapter_info.device_type == wgpu::DeviceType::DiscreteGpu {
            10.0 // Discrete GPU: ~10x speedup
        } else if self.adapter_info.device_type == wgpu::DeviceType::IntegratedGpu {
            3.0 // Integrated GPU: ~3x speedup
        } else {
            1.0 // CPU fallback
        };

        // Scale by problem size - GPU becomes more effective for larger problems
        let size_factor = (problem_size as f64 / 10000.0).clamp(1.0, 5.0);

        gpu_factor * size_factor
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(feature = "gpu")]
    #[test]
    fn test_gpu_turbulence_compute_creation() {
        let compute = GpuTurbulenceCompute::new();
        assert!(
            compute.is_ok(),
            "GPU turbulence compute should create successfully"
        );
    }

    #[cfg(feature = "gpu")]
    #[tokio::test]
    async fn test_smagorinsky_sgs_computation() {
        let mut compute = GpuTurbulenceCompute::new().unwrap();

        // Simple test data
        let nx = 16;
        let ny = 16;
        let velocity_u = vec![1.0; nx * ny];
        let velocity_v = vec![0.0; nx * ny];

        let result = compute.compute_smagorinsky_sgs(
            &velocity_u,
            &velocity_v,
            nx,
            ny,
            0.1,
            0.1,
            0.1, // C_S
        );

        assert!(result.is_ok(), "Smagorinsky SGS computation should succeed");
    }

    #[cfg(feature = "gpu")]
    #[test]
    fn test_performance_info() {
        let compute = GpuTurbulenceCompute::new().unwrap();
        let info = compute.performance_info();

        assert!(info.max_work_group_size > 0);
        assert!(info.max_buffer_size > 0);
        assert!(!info.adapter_info.name.is_empty());
    }

    #[cfg(feature = "gpu")]
    #[test]
    fn test_estimated_speedup() {
        let compute = GpuTurbulenceCompute::new().unwrap();
        let info = compute.performance_info();

        let speedup = info.estimated_speedup(10000);
        assert!(speedup > 1.0, "GPU should provide speedup over CPU");
    }
}
