# GPU Configuration (wgpu)

## Overview
- Auto-detected backends: DX12 (Windows), Vulkan (Linux), Metal (macOS), GL as fallback.
- Power preference defaults to `HighPerformance`; configurable.
- Timestamp query feature used when supported; falls back to wall-clock timing otherwise.

## Configuration
- Type: `cfd_math::linear_solver::matrix_free::gpu_compute::Config`
  - `backends: wgpu::Backends` — override backend selection; default `PRIMARY`.
  - `power_preference: wgpu::PowerPreference` — `HighPerformance` or `LowPower`.
  - `enable_timestamps: bool` — request `TIMESTAMP_QUERY` when available.

## Initialization
- Create with defaults: `GpuComputeContext::new().await`.
- Explicit config: `GpuComputeContext::new_with_config(Config { backends, power_preference, enable_timestamps }).await`.
- Error handling: returns `std::io::Error` on adapter/device failure; operators fall back to CPU when GPU unavailable.

## Metrics
- `ComputeShader::execute_with_metrics` returns `DispatchMetrics { duration_ms, timestamp_supported }`.
- Use `device.poll(Maintain::Wait)` to ensure completion before timing capture.

## Platform Notes
- Windows: Prefer DX12; ensure recent drivers.
- Linux: Prefer Vulkan; headless compute supported.
- macOS: Prefer Metal; feature availability varies with OS version.

## Requirements
- `feature="gpu"` must be enabled for GPU acceleration.
- WGSL shaders compiled by wgpu; ensure compatible shader code.

## Upgrade Path
- CPU operators remain; GPU operators provide `apply_gpu` and `apply_gpu_with_metrics` variants.
- Enable `feature="gpu"` and create/retain a `GpuComputeContext` shared across operators.
