//! GPU compute backend using Hephaestus' WGPU provider ABI.

use crate::error::{Error, Result};
use hephaestus_wgpu::{ComputeDevice, WgpuDevice, wgpu};

const REQUIRED_STORAGE_BUFFERS_PER_SHADER_STAGE: u32 = 7;

pub mod buffer;
pub mod constants;
pub mod field_ops;
pub mod kernels;
pub mod poisson_solver;
pub mod shaders;
pub mod turbulence_compute;
pub mod validation_tests;

pub use buffer::GpuBuffer;
pub use field_ops::GpuFieldOps;
pub use kernels::turbulence::TurbulenceGrid;
pub use poisson_solver::{GpuPoissonSolver, PoissonParams};
pub use turbulence_compute::GpuTurbulenceCompute;

/// GPU context for managing device and queue
pub struct GpuContext {
    /// Hephaestus-owned WGPU provider and capability metadata.
    provider: WgpuDevice,
    /// GPU adapter metadata captured during provider acquisition.
    adapter_info: wgpu::AdapterInfo,
    /// Enabled WGPU features for this provider device.
    features: wgpu::Features,
    /// Device limits
    limits: wgpu::Limits,
}

impl GpuContext {
    /// Create a new GPU context
    ///
    /// # Errors
    /// Returns error if GPU device initialization fails or no suitable adapter found
    pub fn create() -> Result<Self> {
        let required_limits = wgpu::Limits {
            // Three velocity inputs, pressure, and three corrected outputs are
            // bound in one operation; requesting the derived limit makes the
            // provider reject unsupported adapters during acquisition.
            max_storage_buffers_per_shader_stage: REQUIRED_STORAGE_BUFFERS_PER_SHADER_STAGE,
            ..wgpu::Limits::downlevel_defaults()
        };
        let provider = WgpuDevice::try_default_with_limits("CFD GPU Device", required_limits)
            .map_err(|error| {
                Error::InvalidConfiguration(format!(
                    "Failed to acquire Hephaestus WGPU provider: {error}"
                ))
            })?;
        let adapter_info = provider.adapter_info().cloned().ok_or_else(|| {
            Error::InvalidConfiguration(
                "Hephaestus WGPU provider did not report adapter metadata".to_string(),
            )
        })?;

        // Log adapter info
        tracing::info!(
            "GPU adapter selected: {} ({:?}) - Backend: {:?}",
            adapter_info.name,
            adapter_info.device_type,
            adapter_info.backend
        );

        let features = provider.features();
        let limits = provider.limits();
        Ok(Self {
            provider,
            adapter_info,
            features,
            limits,
        })
    }

    /// Get the Hephaestus WGPU provider.
    #[must_use]
    pub(crate) fn provider(&self) -> &WgpuDevice {
        &self.provider
    }

    /// Get adapter info
    #[must_use]
    pub fn adapter_info(&self) -> wgpu::AdapterInfo {
        self.adapter_info.clone()
    }

    /// Check if GPU supports required features
    #[must_use]
    pub fn supports_features(&self, features: wgpu::Features) -> bool {
        self.features.contains(features)
    }

    /// Check whether timestamp query metrics are available.
    #[must_use]
    pub fn supports_timestamp_queries(&self) -> bool {
        self.supports_features(wgpu::Features::TIMESTAMP_QUERY)
    }

    /// Block until submitted GPU work visible to this context has completed.
    ///
    /// # Errors
    /// Returns an error if the Hephaestus provider cannot observe completion.
    pub fn synchronize(&self) -> Result<()> {
        self.provider.synchronize().map_err(Error::from)
    }

    /// Get maximum work group size
    #[must_use]
    pub fn max_work_group_size(&self) -> usize {
        self.limits.max_compute_invocations_per_workgroup as usize
    }

    /// Get maximum buffer size
    #[must_use]
    pub fn max_buffer_size(&self) -> usize {
        self.limits.max_buffer_size as usize
    }
}
