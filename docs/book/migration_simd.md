# Hermes: SIMD Lanes and Vectorized Kernels

CFDrs's hot paths — Stencil sweeps, sparse kernels, element-wise arrays —
migrate from hand-written SIMD intrinsics and `packed_simd` to **Hermes**'
portable SIMD lanes.  Hermes detects the host CPU once at startup and routes
every kernel through the best lane available.

## Lane Hierarchy

```rust
pub trait SimdLane: Copy + Send + Sync + 'static {
    const LANES: usize;
    type Scalar: RealField;
    fn splat(v: Self::Scalar) -> Self;
    fn add(self, rhs: Self) -> Self;
    fn mul(self, rhs: Self) -> Self;
    ...
}

pub struct SseLane<F>(PhantomData<F>);     // 4 lanes × f32 / 2 lanes × f64
pub struct Avx2Lane<F>(PhantomData<F>);    // 8 lanes × f32 / 4 lanes × f64
pub struct Avx512Lane<F>(PhantomData<F>); // 16 lanes × f32 / 8 lanes × f64
pub struct NeonLane<F>(PhantomData<F>);    // 4 lanes × f32 / 2 lanes × f64
pub struct WasmLane<F>(PhantomData<F>);    // 4 lanes × f32 / 2 lanes × f64
```

CFDrs CFG-time dispatch (in `cfd-math`'s build script):

```rust
#[cfg(target_arch = "x86_64")]
use hermes_simd::Avx2Lane as ActiveLane;
#[cfg(target_arch = "aarch64")]
use hermes_simd::NeonLane as ActiveLane;
```

The same source compiles to SSE on a Haswell server, AVX-512 on a Sapphire
Rapids node, and NEON on an Apple-Silicon workstation.  No `#ifdef` blocks
inside solver code.

## Migration From packed_simd

| Legacy | Atlas |
|---|---|
| `packed_simd::f64x4` | `hermes_simd::Avx2Lane<f64>` |
| `packed_simd_2::f32x8` | `hermes_simd::Avx2Lane<f32>` |
| manual `#[target_feature(enable = "avx2")]` | `Avx2Lane::add(...)` |
| runtime `is_x86_feature_detected!` | compile-time `#[cfg]` |

The shape of CFDrs SIMD kernels usually shrinks from “instruction-led“ to
“shape-led“:

```rust
use hermes_simd::SimdLane;

#[inline]
pub fn weighted_sum<F: FloatElement>(w: &[F], x: &[F]) -> F {
    let lane = F::Lane::load(w);
    let vec = F::Lane::load(x);
    let product = lane.mul(vec);
    product.reduce_add()
}
```

## Kernel Catalogue (CFDrs)

| Kernel | Hermes spelling |
|---|---|
| Cell-volume sum | `weighted_sum(weights, volumes)` |
| Spectral filter | `F::Lane::convolve(src, kernel)` |
| Stencil sweep | `F::Lane::stencil(grid, weights)` |
| Sparse SpMV | `F::Lane::scatter_mul(indices, vals, accum)` |
| Cavitation source term | `F::Lane::conditional(cond, then_, else_)` |

## Validation Examples

- [`simd_performance_benchmark`](examples/simd_performance_benchmark.md) —
  compares Hermes-vectorized SpMV against the legacy scalar loop.
- [`spectral_performance`](examples/spectral_performance.md) — FFT-driven
  spectral Poisson solve on `NdArray` slices.
- [`gpu_detection`](examples/gpu_detection.md) — when GPU is enabled this
  routes through `hephaestus` lanes; otherwise Hermes handles the same
  problem on CPU.

## Further Reading

- [`hermes-simd` source](../../../hermes/crates/hermes-simd/)
- [Leto: Arrays](migration_arrays.md)
- [Leto: GAT Tiling](migration_gat_tiles.md)
