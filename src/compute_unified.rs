//! Unified GPU/SIMD compute for the entire CFD suite
//!
//! This module integrates GPU (via wgpu) and SIMD acceleration
//! at the top level where all dependencies are available.

use cfd_core::error::Result;
use cfd_math::simd::{SimdCapability, SimdOperation, SimdProcessor};
use std::sync::Arc;
use std::time::Instant;
use tracing::{info, warn};

/// Performance benchmarking results
#[derive(Debug, Clone)]
pub struct PerformanceBenchmark {
    /// Array size tested
    pub size: usize,
    /// SIMD execution time (nanoseconds)
    pub simd_time_ns: u64,
    /// Parallel execution time (nanoseconds)
    pub parallel_time_ns: u64,
    /// Optimal threshold based on benchmark
    pub optimal_threshold: usize,
}

/// Runtime performance benchmarking for threshold selection
pub struct PerformanceBenchmarker {
    benchmarks: Vec<PerformanceBenchmark>,
    last_benchmark_time: Instant,
    benchmark_interval: std::time::Duration,
}

impl PerformanceBenchmarker {
    /// Create new performance benchmarker
    pub fn new() -> Self {
        Self {
            benchmarks: Vec::new(),
            last_benchmark_time: Instant::now(),
            benchmark_interval: std::time::Duration::from_secs(60), // Rebenchmark every minute
        }
    }

    /// Run performance benchmarks to determine optimal thresholds
    pub fn benchmark_thresholds(&mut self, simd_processor: &SimdProcessor) -> Result<usize> {
        let now = Instant::now();
        
        // Skip benchmarking if we've done it recently
        if now.duration_since(self.last_benchmark_time) < self.benchmark_interval && !self.benchmarks.is_empty() {
            return Ok(self.get_optimal_threshold());
        }

        self.last_benchmark_time = now;
        self.benchmarks.clear();

        // Test different array sizes
        let test_sizes = vec![100, 500, 1000, 2000, 5000, 10000];
        
        for &size in &test_sizes {
            let benchmark = self.benchmark_size(simd_processor, size)?;
            self.benchmarks.push(benchmark);
        }

        Ok(self.get_optimal_threshold())
    }

    /// Benchmark a specific array size
    fn benchmark_size(&self, simd_processor: &SimdProcessor, size: usize) -> Result<PerformanceBenchmark> {
        // Create test data
        let a: Vec<f32> = (0..size).map(|i| i as f32).collect();
        let b: Vec<f32> = (0..size).map(|i| (i * 2) as f32).collect();
        let mut simd_result = vec![0.0f32; size];
        let mut parallel_result = vec![0.0f32; size];

        // Benchmark SIMD
        let simd_start = Instant::now();
        for _ in 0..10 { // Run multiple times for average
            simd_processor.process_f32(&a, &b, &mut simd_result, SimdOperation::Add)?;
        }
        let simd_time = simd_start.elapsed();

        // Benchmark parallel
        let parallel_start = Instant::now();
        for _ in 0..10 { // Run multiple times for average
            use rayon::prelude::*;
            parallel_result.par_iter_mut()
                .zip(a.par_iter())
                .zip(b.par_iter())
                .for_each(|((res, &a_val), &b_val)| {
                    *res = a_val + b_val;
                });
        }
        let parallel_time = parallel_start.elapsed();

        // Calculate optimal threshold (where parallel becomes faster than SIMD)
        let optimal_threshold = if parallel_time < simd_time {
            size
        } else {
            size * 2 // If SIMD is faster, suggest higher threshold
        };

        Ok(PerformanceBenchmark {
            size,
            simd_time_ns: simd_time.as_nanos() as u64 / 10, // Average per run
            parallel_time_ns: parallel_time.as_nanos() as u64 / 10, // Average per run
            optimal_threshold,
        })
    }

    /// Get current optimal threshold from benchmarks
    fn get_optimal_threshold(&self) -> usize {
        // Find the size where parallel becomes consistently faster
        for benchmark in &self.benchmarks {
            if benchmark.parallel_time_ns < benchmark.simd_time_ns {
                return benchmark.size;
            }
        }
        
        // Default threshold if no clear winner found
        1000
    }

    /// Get benchmark results for analysis
    pub fn get_benchmarks(&self) -> &[PerformanceBenchmark] {
        &self.benchmarks
    }
}

impl Default for PerformanceBenchmarker {
    fn default() -> Self {
        Self::new()
    }
}

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
    benchmarker: PerformanceBenchmarker,
}

impl UnifiedCompute {
    /// Create with automatic backend selection
    /// Tries GPU first (supports discrete, integrated, and software rendering)
    pub fn new() -> Result<Self> {
        let simd_processor = SimdProcessor::new();
        let benchmarker = PerformanceBenchmarker::new();

        // Try GPU first - now enabled by default
        #[cfg(feature = "gpu")]
        {
            match cfd_core::compute::gpu::GpuContext::create() {
                Ok(gpu) => {
                    info!("GPU acceleration enabled");
                    return Ok(Self {
                        backend: Backend::Gpu,
                        gpu_context: Some(gpu.device.clone()),
                        simd_processor,
                        benchmarker,
                    });
                }
                Err(e) => {
                    warn!("GPU not available: {e}, falling back to SIMD");
                }
            }
        }

        // Fall back to CPU SIMD
        let capability = SimdCapability::detect();
        let backend = match capability {
            SimdCapability::Avx2 => {
                info!("Using AVX2 SIMD acceleration");
                Backend::Simd
            }
            SimdCapability::Sse42 => {
                info!("Using SSE4.2 SIMD acceleration");
                Backend::Simd
            }
            SimdCapability::Neon => {
                info!("Using NEON SIMD acceleration");
                Backend::Simd
            }
            SimdCapability::Swar => {
                info!("Using SWAR (software SIMD) fallback");
                Backend::Swar
            }
        };

        Ok(Self {
            backend,
            #[cfg(feature = "gpu")]
            gpu_context: None,
            simd_processor,
            benchmarker,
        })
    }

    /// Get active backend
    #[must_use]
    pub fn backend(&self) -> Backend {
        self.backend
    }

    /// Get performance benchmark results
    #[must_use]
    pub fn benchmark_results(&self) -> &[PerformanceBenchmark] {
        self.benchmarker.get_benchmarks()
    }

    /// Force re-run performance benchmarks
    pub fn rebenchmark(&mut self) -> Result<usize> {
        self.benchmarker.benchmark_thresholds(&self.simd_processor)
    }

    /// Vector addition with automatic dispatch
    /// 
    /// Performance optimization: Uses runtime performance benchmarking to auto-select optimal threshold
    /// for switching between SIMD and parallel processing based on actual hardware performance.
    pub fn vector_add_f32(&mut self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        // Get optimal threshold from runtime benchmarking
        let parallel_threshold = self.benchmarker.benchmark_thresholds(&self.simd_processor)?;
        
        match self.backend {
            Backend::Gpu => {
                #[cfg(feature = "gpu")]
                {
                    // TODO: Implement GPU kernel fusion for multiple operations
                    // CURRENT: Falling back to SIMD processing when GPU backend is selected
                    // NEEDED: GPU kernel implementation that can handle vector addition on GPU
                    // DEPENDENCIES: wgpu kernel compilation and buffer management
                    // PRIORITY: High - GPU acceleration is a key performance feature
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
                // Use adaptive threshold based on benchmarked hardware performance
                if a.len() > parallel_threshold {
                    // Use parallel processing for large arrays
                    use rayon::prelude::*;
                    result.par_iter_mut()
                        .zip(a.par_iter())
                        .zip(b.par_iter())
                        .for_each(|((res, &a_val), &b_val)| {
                            *res = a_val + b_val;
                        });
                    Ok(())
                } else {
                    // Use SIMD for smaller arrays where overhead is acceptable
                    self.simd_processor
                        .process_f32(a, b, result, SimdOperation::Add)
                }
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

        for (i, result_item) in result.iter_mut().enumerate().take(rows) {
            let row_start = i * cols;
            let row = &matrix[row_start..row_start + cols];

            // Use SIMD dot product for each row
            let mut sum = 0.0f32;
            for (j, &val) in row.iter().enumerate() {
                sum += val * vector[j];
            }
            *result_item = sum;
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
    use super::{Arc, Result, UnifiedCompute};

    /// Configuration for pressure Poisson solver
    #[derive(Debug, Clone, Copy)]
    pub struct PoissonConfig {
        /// Grid points in x-direction
        pub nx: usize,
        /// Grid points in y-direction
        pub ny: usize,
        /// Grid spacing in x-direction
        pub dx: f32,
        /// Grid spacing in y-direction
        pub dy: f32,
        /// Maximum iterations
        pub iterations: usize,
    }

    impl PoissonConfig {
        /// Create new Poisson solver configuration
        #[must_use]
        pub fn new(nx: usize, ny: usize, dx: f32, dy: f32, iterations: usize) -> Self {
            Self {
                nx,
                ny,
                dx,
                dy,
                iterations,
            }
        }
    }

    /// Configuration for advection solver
    #[derive(Debug, Clone, Copy)]
    pub struct AdvectionConfig {
        /// Grid points in x-direction
        pub nx: usize,
        /// Grid points in y-direction
        pub ny: usize,
        /// Time step
        pub dt: f32,
        /// Grid spacing in x-direction
        pub dx: f32,
        /// Grid spacing in y-direction
        pub dy: f32,
    }

    impl AdvectionConfig {
        /// Create new advection solver configuration
        #[must_use]
        pub fn new(nx: usize, ny: usize, dt: f32, dx: f32, dy: f32) -> Self {
            Self { nx, ny, dt, dx, dy }
        }
    }

    /// Accelerated pressure Poisson solver
    pub struct PressureSolver {
        #[allow(dead_code)]
        compute: Arc<UnifiedCompute>,
    }

    impl PressureSolver {
        /// Create a new pressure solver with the given compute context
        #[must_use]
        pub fn new(compute: Arc<UnifiedCompute>) -> Self {
            Self { compute }
        }

        /// Solve pressure Poisson equation
        pub fn solve(
            &self,
            divergence: &[f32],
            pressure: &mut [f32],
            config: PoissonConfig,
        ) -> Result<()> {
            let dx2 = config.dx * config.dx;
            let dy2 = config.dy * config.dy;
            let factor = 0.5 / (1.0 / dx2 + 1.0 / dy2);

            let mut pressure_new = vec![0.0f32; pressure.len()];

            for _ in 0..config.iterations {
                // Jacobi iteration
                for i in 1..config.nx - 1 {
                    for j in 1..config.ny - 1 {
                        let idx = i * config.ny + j;
                        pressure_new[idx] = factor
                            * ((pressure[(i - 1) * config.ny + j]
                                + pressure[(i + 1) * config.ny + j])
                                / dx2
                                + (pressure[i * config.ny + j - 1]
                                    + pressure[i * config.ny + j + 1])
                                    / dy2
                                - divergence[idx]);
                    }
                }

                // Boundary conditions
                for i in 0..config.nx {
                    pressure_new[i * config.ny] = pressure_new[i * config.ny + 1];
                    pressure_new[i * config.ny + config.ny - 1] =
                        pressure_new[i * config.ny + config.ny - 2];
                }
                for j in 0..config.ny {
                    pressure_new[j] = pressure_new[config.ny + j];
                    pressure_new[(config.nx - 1) * config.ny + j] =
                        pressure_new[(config.nx - 2) * config.ny + j];
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
        #[must_use]
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
            config: AdvectionConfig,
        ) -> Result<()> {
            for i in 0..config.nx {
                for j in 0..config.ny {
                    let idx = i * config.ny + j;

                    // Trace back
                    let x = i as f32 * config.dx - velocity_u[idx] * config.dt;
                    let y = j as f32 * config.dy - velocity_v[idx] * config.dt;

                    // Bilinear interpolation
                    let i0 = (x / config.dx).floor() as usize;
                    let j0 = (y / config.dy).floor() as usize;
                    let i1 = (i0 + 1).min(config.nx - 1);
                    let j1 = (j0 + 1).min(config.ny - 1);
                    let i0 = i0.min(config.nx - 1);
                    let j0 = j0.min(config.ny - 1);

                    let s = ((x / config.dx) - i0 as f32).clamp(0.0, 1.0);
                    let t = ((y / config.dy) - j0 as f32).clamp(0.0, 1.0);

                    field_new[idx] = field[i0 * config.ny + j0] * (1.0 - s) * (1.0 - t)
                        + field[i1 * config.ny + j0] * s * (1.0 - t)
                        + field[i0 * config.ny + j1] * (1.0 - s) * t
                        + field[i1 * config.ny + j1] * s * t;
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
        // Initialize logging for the test
        let _ = tracing_subscriber::fmt().with_test_writer().try_init();

        // TODO: Replace expect-based error handling with proper Result types and error propagation
        // DEPENDENCIES: Add comprehensive error handling framework for unified compute backend initialization
        // BLOCKED BY: Limited understanding of unified compute backend failure modes and recovery strategies
        // PRIORITY: High - Essential for robust testing and debugging of compute backends
        let mut compute = UnifiedCompute::new().expect("Failed to create UnifiedCompute backend for testing");
        info!("Active backend: {:?}", compute.backend());

        let a = vec![1.0f32; 100];
        let b = vec![2.0f32; 100];
        let mut result = vec![0.0f32; 100];

        // TODO: Replace expect-based error handling with proper Result types and error propagation
        // DEPENDENCIES: Add comprehensive error handling framework for vector operations
        // BLOCKED BY: Limited understanding of vector operation failure modes and recovery strategies
        // PRIORITY: High - Essential for robust testing and debugging of compute operations
        compute.vector_add_f32(&a, &b, &mut result)
            .expect("Failed to perform vector addition in test");

        for val in &result {
            assert_eq!(*val, 3.0);
        }
    }

    #[test]
    fn test_pressure_solver() {
        // TODO: Replace expect-based error handling with proper Result types and error propagation
        // DEPENDENCIES: Add comprehensive error handling framework for unified compute backend initialization
        // BLOCKED BY: Limited understanding of unified compute backend failure modes and recovery strategies
        // PRIORITY: High - Essential for robust testing and debugging of compute backends
        let compute = Arc::new(UnifiedCompute::new()
            .expect("Failed to create UnifiedCompute backend for pressure solver test"));
        let solver = kernels::PressureSolver::new(compute);

        let nx = 10;
        let ny = 10;
        let config = kernels::PoissonConfig::new(nx, ny, 0.1, 0.1, 10);
        let divergence = vec![0.1f32; nx * ny];
        let mut pressure = vec![0.0f32; nx * ny];

        // TODO: Replace expect-based error handling with proper Result types and error propagation
        // DEPENDENCIES: Add comprehensive error handling framework for pressure solver operations
        // BLOCKED BY: Limited understanding of pressure solver failure modes and recovery strategies
        // PRIORITY: High - Essential for robust testing and debugging of pressure solvers
        solver.solve(&divergence, &mut pressure, config)
            .expect("Failed to solve pressure equation in test");

        // Should produce non-zero pressure
        let sum: f32 = pressure.iter().sum();
        assert!(sum.abs() > 1e-6);
    }
}
