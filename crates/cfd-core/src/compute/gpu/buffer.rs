//! GPU buffer management

use super::GpuContext;
use crate::compute::traits::ComputeBuffer;
use crate::error::{Error, Result};
use bytemuck::{Pod, Zeroable};
use eunomia::RealField;
use hephaestus_wgpu::{wgpu, ComputeDevice, WgpuBuffer as HephaestusWgpuBuffer};
use std::marker::PhantomData;
use std::sync::Arc;

/// GPU buffer wrapper
pub struct GpuBuffer<T: RealField + Pod + Zeroable> {
    /// Hephaestus-owned WGPU buffer.
    pub buffer: HephaestusWgpuBuffer<T>,
    /// Buffer size in elements
    size: usize,
    /// GPU context
    context: Arc<GpuContext>,
    /// Phantom data for type safety
    _phantom: PhantomData<T>,
}

impl<T: RealField + Pod + Zeroable + Copy> Clone for GpuBuffer<T> {
    fn clone(&self) -> Self {
        Self {
            buffer: self.buffer.clone(),
            size: self.size,
            context: self.context.clone(),
            _phantom: PhantomData,
        }
    }
}

impl<T: RealField + Pod + Zeroable + Copy> GpuBuffer<T> {
    /// Create a new GPU buffer
    ///
    /// # Errors
    ///
    /// Returns an error if the GPU buffer creation fails due to insufficient memory,
    /// invalid size parameters, or GPU device issues.
    pub fn new(context: Arc<GpuContext>, size: usize) -> Result<Self> {
        let buffer = context.provider().alloc_zeroed(size).map_err(|error| {
            Error::GpuCompute(format!("Hephaestus WGPU buffer allocation failed: {error}"))
        })?;

        Ok(Self {
            buffer,
            size,
            context,
            _phantom: PhantomData,
        })
    }

    /// Create buffer with initial data
    ///
    /// # Errors
    ///
    /// Returns an error if buffer creation fails due to data size issues,
    /// memory allocation failures, or GPU device limitations.
    pub fn from_data(context: Arc<GpuContext>, data: &[T]) -> Result<Self> {
        let buffer = context.provider().upload(data).map_err(|error| {
            Error::GpuCompute(format!("Hephaestus WGPU buffer upload failed: {error}"))
        })?;

        Ok(Self {
            buffer,
            size: data.len(),
            context,
            _phantom: PhantomData,
        })
    }

    /// Get underlying wgpu buffer
    pub fn buffer(&self) -> &wgpu::Buffer {
        self.buffer.raw()
    }
}

impl<T: RealField + Pod + Zeroable + Copy> ComputeBuffer<T> for GpuBuffer<T> {
    fn size(&self) -> usize {
        self.size
    }

    fn read(&self) -> Result<Vec<T>> {
        let mut result = vec![T::zeroed(); self.size];
        self.context
            .provider()
            .download(&self.buffer, &mut result)
            .map_err(|error| {
                Error::GpuCompute(format!("Hephaestus WGPU buffer download failed: {error}"))
            })?;
        Ok(result)
    }

    fn write(&mut self, data: &[T]) -> Result<()> {
        if data.len() != self.size {
            return Err(Error::InvalidConfiguration(format!(
                "Data size {} doesn't match buffer size {}",
                data.len(),
                self.size
            )));
        }

        self.context
            .provider()
            .write_buffer(&self.buffer, data)
            .map_err(|error| {
                Error::GpuCompute(format!("Hephaestus WGPU buffer write failed: {error}"))
            })?;

        Ok(())
    }

    fn map(&self) -> Option<&[T]> {
        // GPU buffers cannot be directly mapped for reading
        None
    }

    fn map_mut(&mut self) -> Option<&mut [T]> {
        // GPU buffers cannot be directly mapped for writing
        None
    }
}

impl<T: RealField + Pod + Zeroable> std::fmt::Debug for GpuBuffer<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("GpuBuffer")
            .field("size", &self.size)
            .field("buffer", &"<wgpu::Buffer>")
            .field("context", &"<GpuContext>") // Avoid recursive Debug
            .finish()
    }
}
