//! GPU buffer management

use super::GpuContext;
use crate::compute::traits::ComputeBuffer;
use crate::error::{Error, Result};
use bytemuck::{Pod, Zeroable};
use nalgebra::RealField;
use std::marker::PhantomData;
use std::sync::Arc;

/// GPU buffer wrapper
pub struct GpuBuffer<T: RealField + Pod + Zeroable> {
    /// Underlying wgpu buffer
    pub buffer: wgpu::Buffer,
    /// Buffer size in elements
    size: usize,
    /// GPU context
    context: Arc<GpuContext>,
    /// Phantom data for type safety
    _phantom: PhantomData<T>,
}

impl<T: RealField + Pod + Zeroable + Copy> GpuBuffer<T> {
    /// Create a new GPU buffer
    pub fn new(context: Arc<GpuContext>, size: usize) -> Result<Self> {
        let buffer = context.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("CFD Compute Buffer"),
            size: (size * std::mem::size_of::<T>()) as u64,
            usage: wgpu::BufferUsages::STORAGE
                | wgpu::BufferUsages::COPY_DST
                | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });

        Ok(Self {
            buffer,
            size,
            context,
            _phantom: PhantomData,
        })
    }

    /// Create buffer with initial data
    pub fn from_data(context: Arc<GpuContext>, data: &[T]) -> Result<Self> {
        use wgpu::util::DeviceExt;

        let buffer = context
            .device
            .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("CFD Compute Buffer"),
                contents: bytemuck::cast_slice(data),
                usage: wgpu::BufferUsages::STORAGE
                    | wgpu::BufferUsages::COPY_DST
                    | wgpu::BufferUsages::COPY_SRC,
            });

        Ok(Self {
            buffer,
            size: data.len(),
            context,
            _phantom: PhantomData,
        })
    }

    /// Get underlying wgpu buffer
    pub fn buffer(&self) -> &wgpu::Buffer {
        &self.buffer
    }
}

impl<T: RealField + Pod + Zeroable + Copy> ComputeBuffer<T> for GpuBuffer<T> {
    fn size(&self) -> usize {
        self.size
    }

    fn read(&self) -> Result<Vec<T>> {
        // Create staging buffer for reading
        let staging_buffer = self.context.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Staging Buffer"),
            size: (self.size * std::mem::size_of::<T>()) as u64,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        // Copy from GPU buffer to staging buffer
        let mut encoder =
            self.context
                .device
                .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                    label: Some("Buffer Read Encoder"),
                });

        encoder.copy_buffer_to_buffer(
            &self.buffer,
            0,
            &staging_buffer,
            0,
            (self.size * std::mem::size_of::<T>()) as u64,
        );

        self.context.queue.submit(Some(encoder.finish()));

        // Map staging buffer and read data
        let buffer_slice = staging_buffer.slice(..);
        let (tx, rx) = std::sync::mpsc::channel();

        buffer_slice.map_async(wgpu::MapMode::Read, move |result| {
            tx.send(result).unwrap();
        });

        self.context.device.poll(wgpu::Maintain::Wait);

        rx.recv()
            .map_err(|_| {
                Error::InvalidConfiguration("Failed to receive buffer mapping".to_string())
            })?
            .map_err(|_| Error::InvalidConfiguration("Failed to map buffer".to_string()))?;

        let data = buffer_slice.get_mapped_range();
        let result: Vec<T> = bytemuck::cast_slice(&data).to_vec();

        drop(data);
        staging_buffer.unmap();

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
            .queue
            .write_buffer(&self.buffer, 0, bytemuck::cast_slice(data));

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
            .field("buffer_id", &format!("{:?}", self.buffer.global_id()))
            .finish()
    }
}
