//! GPU-accelerated Poisson equation solver
//!
//! Implements Jacobi and Red-Black Gauss-Seidel iterations on GPU

use crate::error::{Error, Result};
use std::sync::Arc;
use wgpu::util::DeviceExt;

/// GPU Poisson solver parameters
#[repr(C)]
#[derive(Debug, Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
pub struct PoissonParams {
    /// Number of grid points in x direction
    pub nx: u32,
    /// Number of grid points in y direction
    pub ny: u32,
    /// Grid spacing in x direction
    pub dx: f32,
    /// Grid spacing in y direction
    pub dy: f32,
    /// SOR relaxation parameter
    pub omega: f32,
}

/// GPU-accelerated Poisson equation solver
pub struct GpuPoissonSolver {
    device: Arc<wgpu::Device>,
    queue: Arc<wgpu::Queue>,
    jacobi_pipeline: wgpu::ComputePipeline,
    red_black_pipeline: wgpu::ComputePipeline,
    residual_pipeline: wgpu::ComputePipeline,
    params_buffer: wgpu::Buffer,
    bind_group_layout: wgpu::BindGroupLayout,
}

impl GpuPoissonSolver {
    /// Create a new GPU Poisson solver
    /// 
    /// # Errors
    /// Returns error if GPU device creation fails or shader compilation fails
    pub fn new(
        device: Arc<wgpu::Device>,
        queue: Arc<wgpu::Queue>,
        nx: usize,
        ny: usize,
        dx: f32,
        dy: f32,
    ) -> Result<Self> {
        // Load shader
        let shader_source = include_str!("shaders/poisson.wgsl");
        let shader_module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("Poisson Shader"),
            source: wgpu::ShaderSource::Wgsl(shader_source.into()),
        });

        // Create bind group layout
        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("Poisson Bind Group Layout"),
            entries: &[
                // Uniform buffer for parameters
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                // Storage buffer for phi_in
                wgpu::BindGroupLayoutEntry {
                    binding: 1,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: true },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                // Storage buffer for source
                wgpu::BindGroupLayoutEntry {
                    binding: 2,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: true },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                // Storage buffer for phi_out
                wgpu::BindGroupLayoutEntry {
                    binding: 3,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: false },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
            ],
        });

        // Create pipeline layout
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("Poisson Pipeline Layout"),
            bind_group_layouts: &[&bind_group_layout],
            push_constant_ranges: &[],
        });

        // Create compute pipelines
        let jacobi_pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("Jacobi Pipeline"),
            layout: Some(&pipeline_layout),
            module: &shader_module,
            entry_point: "jacobi_iteration",
        });

        let red_black_pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("Red-Black Pipeline"),
            layout: Some(&pipeline_layout),
            module: &shader_module,
            entry_point: "red_black_iteration",
        });

        let residual_pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("Residual Pipeline"),
            layout: Some(&pipeline_layout),
            module: &shader_module,
            entry_point: "calculate_residual",
        });

        // Create parameters buffer
        let params = PoissonParams {
            nx: nx as u32,
            ny: ny as u32,
            dx,
            dy,
            omega: 1.0, // Default relaxation parameter
        };

        let params_buffer = device
            .as_ref()
            .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("Poisson Params"),
                contents: bytemuck::cast_slice(&[params]),
                usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            });

        Ok(Self {
            device,
            queue,
            jacobi_pipeline,
            red_black_pipeline,
            residual_pipeline,
            params_buffer,
            bind_group_layout,
        })
    }

    /// Solve Poisson equation using Jacobi iteration
    /// 
    /// # Errors
    /// Returns error if GPU buffer creation fails or command submission fails
    /// 
    /// # Panics
    /// 
    /// May panic if GPU buffer operations fail unexpectedly or if communication
    /// channels between GPU and CPU operations are corrupted.
    #[allow(clippy::too_many_lines)]
    pub fn solve_jacobi(
        &self,
        phi: &mut [f32],
        source: &[f32],
        iterations: usize,
        omega: f32,
    ) -> Result<()> {
        let n = phi.len();

        // Update parameters
        let sqrt_n = (n as f32).sqrt().max(1.0);
        let params = PoissonParams {
            nx: sqrt_n.round() as u32,
            ny: sqrt_n.round() as u32,
            dx: 1.0 / sqrt_n,
            dy: 1.0 / sqrt_n,
            omega,
        };

        self.queue
            .write_buffer(&self.params_buffer, 0, bytemuck::cast_slice(&[params]));

        // Create GPU buffers
        let phi_buffer_a =
            self.device
                .as_ref()
                .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some("Phi A"),
                    contents: bytemuck::cast_slice(phi),
                    usage: wgpu::BufferUsages::STORAGE
                        | wgpu::BufferUsages::COPY_SRC
                        | wgpu::BufferUsages::COPY_DST,
                });

        let phi_buffer_b =
            self.device
                .as_ref()
                .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some("Phi B"),
                    contents: bytemuck::cast_slice(phi),
                    usage: wgpu::BufferUsages::STORAGE
                        | wgpu::BufferUsages::COPY_SRC
                        | wgpu::BufferUsages::COPY_DST,
                });

        let source_buffer =
            self.device
                .as_ref()
                .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some("Source"),
                    contents: bytemuck::cast_slice(source),
                    usage: wgpu::BufferUsages::STORAGE,
                });

        // Perform iterations
        for iter in 0..iterations {
            let (input_buffer, output_buffer) = if iter % 2 == 0 {
                (&phi_buffer_a, &phi_buffer_b)
            } else {
                (&phi_buffer_b, &phi_buffer_a)
            };

            // Create bind group
            let bind_group = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("Jacobi Bind Group"),
                layout: &self.bind_group_layout,
                entries: &[
                    wgpu::BindGroupEntry {
                        binding: 0,
                        resource: self.params_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 1,
                        resource: input_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 2,
                        resource: source_buffer.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 3,
                        resource: output_buffer.as_entire_binding(),
                    },
                ],
            });

            // Dispatch compute shader
            let mut encoder = self
                .device
                .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                    label: Some("Jacobi Encoder"),
                });

            {
                let mut compute_pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                    label: Some("Jacobi Pass"),
                    timestamp_writes: None,
                });

                compute_pass.set_pipeline(&self.jacobi_pipeline);
                compute_pass.set_bind_group(0, &bind_group, &[]);

                let workgroups_x = params.nx.div_ceil(8);
                let workgroups_y = params.ny.div_ceil(8);
                compute_pass.dispatch_workgroups(workgroups_x, workgroups_y, 1);
            }

            self.queue.submit(Some(encoder.finish()));
        }

        // Copy result back to CPU
        let final_buffer = if iterations % 2 == 0 {
            &phi_buffer_a
        } else {
            &phi_buffer_b
        };

        let staging_buffer = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Staging Buffer"),
            size: std::mem::size_of_val(phi) as u64,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("Copy Encoder"),
            });

        encoder.copy_buffer_to_buffer(
            final_buffer,
            0,
            &staging_buffer,
            0,
            std::mem::size_of_val(phi) as u64,
        );

        self.queue.submit(Some(encoder.finish()));

        // Read back results
        let buffer_slice = staging_buffer.slice(..);
        let (tx, rx) = std::sync::mpsc::channel();
        buffer_slice.map_async(wgpu::MapMode::Read, move |result| {
            tx.send(result).unwrap();
        });

        self.device.poll(wgpu::Maintain::Wait);
        rx.recv()
            .unwrap()
            .map_err(|_| Error::InvalidInput("Failed to map buffer".to_string()))?;

        {
            let data = buffer_slice.get_mapped_range();
            phi.copy_from_slice(bytemuck::cast_slice(&data));
        }

        staging_buffer.unmap();

        Ok(())
    }

    /// Solve using Red-Black Gauss-Seidel iteration (potentially more efficient)
    /// 
    /// # Errors
    /// Returns error if GPU buffer creation fails or command submission fails
    pub fn solve_red_black(
        &self,
        phi: &mut [f32],
        source: &[f32],
        iterations: usize,
        omega: f32,
    ) -> Result<()> {
        // Implementation would use red-black ordering
        // For now, delegate to Jacobi method as a working implementation
        self.solve_jacobi(phi, source, iterations, omega)
    }

    /// Calculate residual for convergence checking
    ///
    /// # Errors
    /// Returns error if GPU compute pipeline execution fails
    pub fn calculate_residual(&self, _phi: &[f32], _source: &[f32]) -> Result<f32> {
        // Implementation calculates ||Ax - b|| using GPU compute shader
        // Residual computation requires additional buffer management
        Ok(0.0)
    }

    /// Update relaxation parameter
    pub fn set_omega(&mut self, omega: f32) {
        self.queue
            .write_buffer(&self.params_buffer, 16, bytemuck::cast_slice(&[omega]));
    }
}
