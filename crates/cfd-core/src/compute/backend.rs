//! Compute backend management and capability detection

use super::traits::ComputeBackend;
use std::sync::Arc;

/// Compute capability information
#[derive(Debug, Clone)]
pub struct ComputeCapability {
    /// Available backends
    pub backends: Vec<ComputeBackend>,
    /// Maximum work group size for GPU
    pub max_work_group_size: usize,
    /// Maximum buffer size (bytes)
    pub max_buffer_size: usize,
    /// Number of compute units (GPU cores or CPU threads)
    pub compute_units: usize,
    /// Available memory (bytes)
    pub available_memory: usize,
    /// Preferred backend based on problem size
    pub preferred_backend: ComputeBackend,
}

impl ComputeCapability {
    /// Detect system compute capabilities
    pub fn detect() -> Self {
        let mut backends = vec![ComputeBackend::Cpu];

        if ComputeBackend::Simd.is_available() {
            backends.push(ComputeBackend::Simd);
        }

        if ComputeBackend::Gpu.is_available() {
            backends.push(ComputeBackend::Gpu);
        }

        if backends.len() > 1 {
            backends.push(ComputeBackend::Hybrid);
        }

        let compute_units = num_cpus::get();
        let available_memory = Self::estimate_available_memory();

        // Determine preferred backend based on available options
        let preferred_backend = if backends.contains(&ComputeBackend::Gpu) {
            ComputeBackend::Gpu
        } else if backends.contains(&ComputeBackend::Simd) {
            ComputeBackend::Simd
        } else {
            ComputeBackend::Cpu
        };

        Self {
            backends,
            max_work_group_size: 256,              // Conservative default
            max_buffer_size: available_memory / 4, // Use up to 25% of memory
            compute_units,
            available_memory,
            preferred_backend,
        }
    }

    /// Estimate available system memory
    fn estimate_available_memory() -> usize {
        // Conservative estimate: 2GB
        2 * 1024 * 1024 * 1024
    }

    /// Select best backend for given problem size
    pub fn select_backend(&self, problem_size: usize) -> ComputeBackend {
        const GPU_THRESHOLD: usize = 100_000; // Minimum size for GPU efficiency
        const SIMD_THRESHOLD: usize = 1_000; // Minimum size for SIMD efficiency

        if problem_size >= GPU_THRESHOLD && self.backends.contains(&ComputeBackend::Gpu) {
            ComputeBackend::Gpu
        } else if problem_size >= SIMD_THRESHOLD && self.backends.contains(&ComputeBackend::Simd) {
            ComputeBackend::Simd
        } else {
            ComputeBackend::Cpu
        }
    }
}

/// Backend context for managing compute resources
pub struct BackendContext {
    /// Current backend
    pub backend: ComputeBackend,
    /// System capabilities
    pub capabilities: Arc<ComputeCapability>,
    #[cfg(feature = "gpu")]
    /// GPU context if available
    pub gpu_context: Option<Arc<super::gpu::GpuContext>>,
}

impl BackendContext {
    /// Create a new backend context
    pub fn new(backend: ComputeBackend) -> crate::error::Result<Self> {
        let capabilities = Arc::new(ComputeCapability::detect());

        if !capabilities.backends.contains(&backend) {
            return Err(crate::error::Error::InvalidConfiguration(format!(
                "Backend {:?} not available on this system",
                backend
            )));
        }

        Ok(Self {
            backend,
            capabilities,
            #[cfg(feature = "gpu")]
            gpu_context: if backend == ComputeBackend::Gpu || backend == ComputeBackend::Hybrid {
                Some(Arc::new(super::gpu::GpuContext::create()?))
            } else {
                None
            },
        })
    }

    /// Switch to a different backend
    ///
    /// # Errors
    /// Returns an error if:
    /// - The requested backend is not available on this system
    /// - Backend initialization fails
    /// - Hardware capabilities are insufficient
    pub fn switch_backend(&mut self, backend: ComputeBackend) -> crate::error::Result<()> {
        if !self.capabilities.backends.contains(&backend) {
            return Err(crate::error::Error::InvalidConfiguration(format!(
                "Backend {:?} not available on this system",
                backend
            )));
        }

        self.backend = backend;

        #[cfg(feature = "gpu")]
        if (backend == ComputeBackend::Gpu || backend == ComputeBackend::Hybrid)
            && self.gpu_context.is_none()
        {
            self.gpu_context = Some(Arc::new(super::gpu::GpuContext::create()?));
        }

        Ok(())
    }
}
