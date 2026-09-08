//! GPU compute backend using Hephaestus' WGPU provider ABI.

use crate::error::{Error, Result};
use hephaestus_wgpu::{
    ComputeDevice, ComputeDeviceAcquisition, ComputeDeviceCapabilities, DeviceFeature,
    DeviceLimits, DevicePreference, WgpuDevice,
};

const REQUIRED_STORAGE_BUFFERS_PER_SHADER_STAGE: u32 = 7;

/// The velocity-correction kernel binds seven storage buffers (three velocity
/// inputs, pressure, three corrected outputs) plus one uniform `Params` buffer.
/// WGPU 30 counts all eight against the combined per-stage budget, so this must
/// exceed `REQUIRED_STORAGE_BUFFERS_PER_SHADER_STAGE` by the uniform binding.
const REQUIRED_BUFFERS_AND_ACCELERATION_STRUCTURES_PER_SHADER_STAGE: u32 = 8;

pub mod buffer;
pub mod constants;
pub mod field_ops;
pub mod kernels;
pub mod poisson_solver;
pub mod turbulence_compute;
pub mod validation_tests;

pub use buffer::GpuBuffer;
pub use field_ops::GpuFieldOps;
pub use kernels::turbulence::TurbulenceGrid;
pub use poisson_solver::{GpuPoissonSolver, PoissonParams};
pub use turbulence_compute::GpuTurbulenceCompute;

/// GPU context for managing device and queue
pub struct GpuContext {
    /// Hephaestus-owned device provider.
    provider: WgpuDevice,
}

impl GpuContext {
    /// Create a new GPU context
    ///
    /// # Errors
    /// Returns error if GPU device initialization fails or no suitable adapter found
    pub fn create() -> Result<Self> {
        let provider = WgpuDevice::try_acquire_device(
            "CFD GPU Device",
            DevicePreference::HighPerformance,
            &[],
            required_device_limits(),
        )
        .map_err(|error| {
            Error::InvalidConfiguration(format!(
                "Failed to acquire Hephaestus GPU provider: {error}"
            ))
        })?;

        tracing::info!(backend = provider.backend_name(), "GPU provider selected");

        Ok(Self { provider })
    }

    /// Get the Hephaestus WGPU provider.
    #[must_use]
    pub(crate) fn provider(&self) -> &WgpuDevice {
        &self.provider
    }

    /// Return the selected provider's backend identifier.
    #[must_use]
    pub fn backend_name(&self) -> &'static str {
        self.provider.backend_name()
    }

    /// Check whether timestamp query metrics are available.
    #[must_use]
    pub fn supports_timestamp_queries(&self) -> bool {
        ComputeDeviceCapabilities::supports_device_feature(
            &self.provider,
            DeviceFeature::TimestampQuery,
        )
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
        ComputeDeviceCapabilities::device_limits(&self.provider)
            .max_compute_invocations_per_workgroup as usize
    }

    /// Get maximum buffer size
    #[must_use]
    pub fn max_buffer_size(&self) -> usize {
        ComputeDeviceCapabilities::device_limits(&self.provider).max_buffer_size as usize
    }
}

fn required_device_limits() -> DeviceLimits {
    let mut limits = WgpuDevice::downlevel_device_limits();
    // Three velocity inputs, pressure, and three corrected outputs are bound
    // in one operation; provider acquisition rejects unsupported adapters.
    limits.max_storage_buffers_per_shader_stage = Some(REQUIRED_STORAGE_BUFFERS_PER_SHADER_STAGE);
    // Those seven storage buffers plus the uniform Params buffer count against
    // WGPU 30's combined per-stage budget; request all eight explicitly so the
    // provider does not derive a seven-slot ceiling from the storage limit.
    limits.max_buffers_and_acceleration_structures_per_shader_stage =
        Some(REQUIRED_BUFFERS_AND_ACCELERATION_STRUCTURES_PER_SHADER_STAGE);
    limits
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn required_device_limits_preserve_downlevel_baseline() {
        let limits = required_device_limits();
        let mut expected = WgpuDevice::downlevel_device_limits();
        expected.max_storage_buffers_per_shader_stage =
            Some(REQUIRED_STORAGE_BUFFERS_PER_SHADER_STAGE);
        expected.max_buffers_and_acceleration_structures_per_shader_stage =
            Some(REQUIRED_BUFFERS_AND_ACCELERATION_STRUCTURES_PER_SHADER_STAGE);

        assert_eq!(limits, expected);
    }
}
