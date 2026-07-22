# Example: simd_performance_benchmark

**Crate**: `cfd-suite` (workspace root)  
**Run**: `cargo run --example simd_performance_benchmark`  
**Source**: [`examples/simd_performance_benchmark.rs`](../../../examples/simd_performance_benchmark.rs)

## What This Example Demonstrates

Runtime SIMD capability detection and throughput benchmark for CFD gradient and
flux kernels using `cfd-math::simd::CfdSimdOps` backed by `hermes`.

| Grid | Cells | Expected speedup |
|---|---|---|
| 16×16 | 256 | 2–4× on AVX2 |
| 32×32 | 1 024 | 2–4× |
| 64×64 | 4 096 | 2–4× |
| 128×64 | 8 192 | 2–4× |

## Key Code Snippet

```rust
use cfd_math::simd::cfd::CfdSimdOps;
use std::time::Instant;

let ops = CfdSimdOps::<f64>::new();
// ops dispatches to AVX2, NEON, or scalar depending on the host CPU

let t0 = Instant::now();
// ... gradient kernel call ...
println!("{}×{}: {:.3} ms", nx, ny, t0.elapsed().as_secs_f64() * 1e3);
```

## Atlas Integration

`CfdSimdOps` wraps the `hermes-simd` SIMD abstraction layer with runtime
dispatch so CFDrs code achieves SIMD acceleration without `unsafe` and
without committing to a specific ISA at compile time.

Expected fallback on hardware without AVX2/NEON: scalar path, 1× throughput.

## Book Chapter

[← SIMD, GPU, and Backend Migration](../performance_and_atlas.md)

