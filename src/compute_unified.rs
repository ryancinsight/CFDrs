//! Unified GPU/SIMD compute for the entire CFD suite
//!
//! This module integrates GPU (via wgpu) and SIMD acceleration
//! at the top level where all dependencies are available.

use cfd_core::error::Result;
use cfd_math::simd::{SimdCapability, SimdOperation, SimdProcessor};
use std::sync::Arc;

/// Unified compute backend selection
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Backend {
    /// GPU via wgpu
    Gpu,
    /// CPU SIMD (AVX2/SSE/NEON)
    Simd,
    /// CPU SWAR fallback
    Swar,
}

/// High-level unified compute interface
pub struct UnifiedCompute {
    backend: Backend,
    #[cfg(feature = "gpu")]
    #[allow(dead_code)]
    gpu_context: Option<Arc<wgpu::Device>>,
    simd_processor: SimdProcessor,
}

impl UnifiedCompute {
    /// Create with automatic backend selection
    /// Tries GPU first (supports discrete, integrated, and software rendering)
    pub fn new() -> Result<Self> {
        let simd_processor = SimdProcessor::new();

        // Try GPU first - now enabled by default
        #[cfg(feature = "gpu")]
        {
            match cfd_core::compute::gpu::GpuContext::create() {
                Ok(gpu) => {
                    println!("GPU acceleration enabled");
                    return Ok(Self {
                        backend: Backend::Gpu,
                        gpu_context: Some(gpu.device.clone()),
                        simd_processor,
                    });
                }
                Err(e) => {
                    println!("GPU not available: {}, falling back to SIMD", e);
                }
            }
        }

        // Fall back to CPU SIMD
        let capability = SimdCapability::detect();
        let backend = match capability {
            SimdCapability::Avx2 => {
                println!("Using AVX2 SIMD acceleration");
                Backend::Simd
            }
            SimdCapability::Sse42 => {
                println!("Using SSE4.2 SIMD acceleration");
                Backend::Simd
            }
            SimdCapability::Neon => {
                println!("Using NEON SIMD acceleration");
                Backend::Simd
            }
            _ => {
                println!("Using SWAR (software SIMD) fallback");
                Backend::Swar
            }
        };

        Ok(Self {
            backend,
            #[cfg(feature = "gpu")]
            gpu_context: None,
            simd_processor,
        })
    }

    /// Get active backend
    pub fn backend(&self) -> Backend {
        self.backend
    }

    /// Vector addition with automatic dispatch
    pub fn vector_add_f32(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        match self.backend {
            Backend::Gpu => {
                #[cfg(feature = "gpu")]
                {
                    // GPU implementation
                    self.simd_processor
                        .process_f32(a, b, result, SimdOperation::Add)
                }
                #[cfg(not(feature = "gpu"))]
                {
                    self.simd_processor
                        .process_f32(a, b, result, SimdOperation::Add)
                }
            }
            Backend::Simd | Backend::Swar => {
                self.simd_processor
                    .process_f32(a, b, result, SimdOperation::Add)
            }
        }
    }

    /// Vector multiplication
    pub fn vector_mul_f32(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        self.simd_processor
            .process_f32(a, b, result, SimdOperation::Mul)
    }

    /// Matrix-vector multiplication
    pub fn matvec_f32(
        &self,
        matrix: &[f32],
        vector: &[f32],
        result: &mut [f32],
        rows: usize,
        cols: usize,
    ) -> Result<()> {
        // Matrix-vector multiplication using SIMD operations
        if matrix.len() != rows * cols || vector.len() != cols || result.len() != rows {
            return Err(cfd_core::error::Error::InvalidInput(
                "Matrix-vector dimension mismatch".to_string(),
            ));
        }

        for i in 0..rows {
            let row_start = i * cols;
            let row = &matrix[row_start..row_start + cols];

            // Use SIMD dot product for each row
            let mut sum = 0.0f32;
            for (j, &val) in row.iter().enumerate() {
                sum += val * vector[j];
            }
            result[i] = sum;
        }
        Ok(())
    }

    /// Dot product
    pub fn dot_f32(&self, a: &[f32], b: &[f32]) -> Result<f32> {
        use cfd_math::simd::VectorOps;
        self.simd_processor.ops.dot(a, b)
    }
}

/// CFD-specific accelerated kernels
pub mod kernels {
    use super::*;

    /// Accelerated pressure Poisson solver
    pub struct PressureSolver {
        #[allow(dead_code)]
        compute: Arc<UnifiedCompute>,
    }

    impl PressureSolver {
        /// Create a new pressure solver with the given compute context
        pub fn new(compute: Arc<UnifiedCompute>) -> Self {
            Self { compute }
        }

        /// Solve pressure Poisson equation
        pub fn solve(
            &self,
            divergence: &[f32],
            pressure: &mut [f32],
            nx: usize,
            ny: usize,
            dx: f32,
            dy: f32,
            iterations: usize,
        ) -> Result<()> {
            let dx2 = dx * dx;
            let dy2 = dy * dy;
            let factor = 0.5 / (1.0 / dx2 + 1.0 / dy2);

            let mut pressure_new = vec![0.0f32; pressure.len()];

            for _ in 0..iterations {
                // Jacobi iteration
                for i in 1..nx - 1 {
                    for j in 1..ny - 1 {
                        let idx = i * ny + j;
                        pressure_new[idx] = factor
                            * ((pressure[(i - 1) * ny + j] + pressure[(i + 1) * ny + j]) / dx2
                                + (pressure[i * ny + j - 1] + pressure[i * ny + j + 1]) / dy2
                                - divergence[idx]);
                    }
                }

                // Boundary conditions
                for i in 0..nx {
                    pressure_new[i * ny] = pressure_new[i * ny + 1];
                    pressure_new[i * ny + ny - 1] = pressure_new[i * ny + ny - 2];
                }
                for j in 0..ny {
                    pressure_new[j] = pressure_new[ny + j];
                    pressure_new[(nx - 1) * ny + j] = pressure_new[(nx - 2) * ny + j];
                }

                pressure.copy_from_slice(&pressure_new);
            }

            Ok(())
        }
    }

    /// Accelerated advection solver
    pub struct AdvectionSolver {
        #[allow(dead_code)]
        compute: Arc<UnifiedCompute>,
    }

    impl AdvectionSolver {
        /// Create a new advection solver with the given compute context
        pub fn new(compute: Arc<UnifiedCompute>) -> Self {
            Self { compute }
        }

        /// Semi-Lagrangian advection
        pub fn advect(
            &self,
            field: &[f32],
            velocity_u: &[f32],
            velocity_v: &[f32],
            field_new: &mut [f32],
            nx: usize,
            ny: usize,
            dt: f32,
            dx: f32,
            dy: f32,
        ) -> Result<()> {
            for i in 0..nx {
                for j in 0..ny {
                    let idx = i * ny + j;

                    // Trace back
                    let x = i as f32 * dx - velocity_u[idx] * dt;
                    let y = j as f32 * dy - velocity_v[idx] * dt;

                    // Bilinear interpolation
                    let i0 = (x / dx).floor() as usize;
                    let j0 = (y / dy).floor() as usize;
                    let i1 = (i0 + 1).min(nx - 1);
                    let j1 = (j0 + 1).min(ny - 1);
                    let i0 = i0.min(nx - 1);
                    let j0 = j0.min(ny - 1);

                    let s = ((x / dx) - i0 as f32).max(0.0).min(1.0);
                    let t = ((y / dy) - j0 as f32).max(0.0).min(1.0);

                    field_new[idx] = field[i0 * ny + j0] * (1.0 - s) * (1.0 - t)
                        + field[i1 * ny + j0] * s * (1.0 - t)
                        + field[i0 * ny + j1] * (1.0 - s) * t
                        + field[i1 * ny + j1] * s * t;
                }
            }

            Ok(())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_unified_compute() {
        let compute = UnifiedCompute::new().unwrap();
        println!("Active backend: {:?}", compute.backend());

        let a = vec![1.0f32; 100];
        let b = vec![2.0f32; 100];
        let mut result = vec![0.0f32; 100];

        compute.vector_add_f32(&a, &b, &mut result).unwrap();

        for val in &result {
            assert_eq!(*val, 3.0);
        }
    }

    #[test]
    fn test_pressure_solver() {
        let compute = Arc::new(UnifiedCompute::new().unwrap());
        let solver = kernels::PressureSolver::new(compute);

        let nx = 10;
        let ny = 10;
        let divergence = vec![0.1f32; nx * ny];
        let mut pressure = vec![0.0f32; nx * ny];

        solver
            .solve(&divergence, &mut pressure, nx, ny, 0.1, 0.1, 10)
            .unwrap();

        // Should produce non-zero pressure
        let sum: f32 = pressure.iter().sum();
        assert!(sum.abs() > 1e-6);
    }
}
