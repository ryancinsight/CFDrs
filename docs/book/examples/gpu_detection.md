# Example: gpu_detection

**Crate**: `cfd-suite` (workspace root)  
**Run**: `cargo run --example gpu_detection`  
**Source**: [`examples/gpu_detection.rs`](../../examples/gpu_detection.rs)

## What This Example Demonstrates

Automatic compute-backend selection: discrete GPU → integrated GPU → SIMD → scalar,
using the unified `ComputeDispatcher` backed by `hephaestus-wgpu`.

| API | Purpose |
|---|---|
| `ComputeCapability::detect()` | Query available backends, memory, compute units |
| `ComputeDispatcher::new()` | Instantiate preferred backend |
| `dispatcher.current_backend()` | Inspect selected backend at runtime |

## Key Code Snippet

```rust
use cfd_core::compute::backend::ComputeCapability;
use cfd_core::compute::dispatch::ComputeDispatcher;
use cfd_core::compute::traits::ComputeBackend;

let caps = ComputeCapability::detect();
println!("Preferred: {:?}", caps.preferred_backend);
println!("Compute units: {}", caps.compute_units);
println!("Memory: {} MB", caps.available_memory / (1024 * 1024));

let dispatcher = ComputeDispatcher::new()?;
match dispatcher.current_backend() {
    ComputeBackend::Gpu  => println!("GPU active"),
    ComputeBackend::Simd => println!("SIMD active"),
    ComputeBackend::Cpu  => println!("Scalar fallback"),
}
```

## Atlas Integration

`ComputeDispatcher` wraps `hephaestus-wgpu` for the GPU path and `hermes-simd`
for the SIMD path. CFDrs kernels call through this dispatcher so they run on
whichever backend is available without compile-time feature flags.

## Book Chapter

[← SIMD, GPU, and Backend Migration](../performance_and_atlas.md)

