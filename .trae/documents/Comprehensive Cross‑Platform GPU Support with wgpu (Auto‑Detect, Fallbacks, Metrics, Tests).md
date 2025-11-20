## Overview
- Add robust GPU compute context built on wgpu with auto‑detected adapter/device, cross‑platform backends (DX12/Vulkan/Metal/GL), and CPU fallbacks.
- Preserve CPU paths for backward compatibility; provide clear upgrade hooks for GPU acceleration.

## Architecture
- GPU Context: Singleton/RAII `GpuComputeContext` with:
  - Instance/backends: auto‑select backend per OS; configurable override.
  - Adapter selection: `HighPerformance` then `LowPower`; respect `required_features` and limits.
  - Device/queue creation; timestamp query support for metrics when available.
  - Error handling via `thiserror`; `Result` throughout.
- Compute Shader Abstraction:
  - `ComputeShader` builder: module, pipeline, bind‑group layout(s), dispatch utilities.
  - Buffer wrappers: typed storage/uniform buffers; map/copy helpers; strict usage flags.
- Performance Metrics:
  - Enable timestamp queries (if supported) to measure GPU pass durations.
  - Fallback metrics when queries unsupported: wall‑clock timing around dispatch.
  - Metric structs collected per operator run; report via `tracing` and return opt data.
- Fallbacks:
  - No adapter/device → automatic CPU fallback with logged reason.
  - Feature unavailability (e.g., timestamp queries) → graceful downgrade.

## Backward Compatibility & Upgrade Path
- Keep existing CPU `LinearOperator` behavior.
- Feature‑gated GPU (`feature="gpu"`); when enabled, prefer GPU if adapter found, else CPU.
- Operators expose unified API; add optional `apply_gpu` for explicit control.

## Configuration & Docs
- Config options:
  - Backend override (DX12/Vulkan/Metal/GL/WebGPU), power preference, required features, limits, workgroup sizes.
  - Metrics toggle, debug validation layers, shader caching.
- Document options/requirements in crate docs and a GPU configuration guide.

## Tests
- Feature‑gated GPU tests:
  - Adapter detection: skip tests if no adapter; assert informative messages.
  - Numerical parity: GPU vs CPU outputs within tolerance for Laplacian (Dirichlet/Periodic/mixed), Momentum/Energy.
  - Resource lifecycle: create/use/drop buffers and pipelines without leaks; map/copy correctness.
  - Metrics presence: timestamp queries where supported; fallback timing recorded otherwise.
- Property tests: invariants (linearity, symmetry) on GPU path with tolerances.

## Cross‑Platform Optimizations
- Windows: prefer DX12 backend; enable shader compiler settings as needed.
- Linux: prefer Vulkan; validate driver limits; handle headless compute.
- macOS: prefer Metal; respect feature availability.

## Quality & Cleanup
- RAII cleanup for buffers/pipelines; queue flush and `device.poll` where necessary.
- Structured error handling and tracing spans around initialization/dispatch.

## Implementation Plan (Repo)
1. Implement `GpuComputeContext` in `cfd-math` (or `cfd-core` if shared) with auto‑detect init, error types, and metrics support.
2. Update `ComputeShader` and buffer wrappers to real wgpu implementations (load WGSL, create pipeline/layouts/bind‑groups).
3. Integrate context into existing GPU operators (Laplacian/Momentum/Poisson) with metrics and CPU fallback.
4. Add docs: GPU config guide and API rustdoc detailing options/invariants.
5. Add tests: adapter detection, parity, metrics, lifecycle; gate by `feature="gpu"` and skip when unavailable.
6. Add Criterion benches comparing CPU vs GPU where possible (opt‑in).

## Success Criteria
- Auto‑detect adapter/device; no crashes on unsupported systems; CPU fallback works.
- GPU operators produce results within documented tolerances vs CPU.
- Metrics available when supported; benches show speedups on capable hardware.
- Docs clearly describe configuration and platform nuances; API remains compatible.

## Risks & Contingencies
- Adapter absence or driver limitations → documented skip with CPU fallback.
- Timestamp queries unsupported → fallback timers.
- CI hardware variability → mark heavy GPU tests as skipped without adapter; rely on local GPU rigs for performance benches.
