# Migration Overview: ndarray/nalgebra → Atlas Stack

CFDrs is mid-flight in a bulk-then-cleanup migration from third-party mathematics
crates to the unified **Atlas** stack.  This Part documents the destination
crates, the principles they encode, and the per-crate migration surface that
each solver touches.

## Why Atlas

The Atlas stack replaces four categories of third-party dependency with one
coherent set of first-party crates.  The motivation is structural, not
feature-driven:

| Legacy dependency | Atlas replacement | Reason |
|---|---|---|
| `ndarray`, `nalgebra::Matrix` | `leto::NdArray<T,D>` (CPU) | one dense type across all crates |
| `nalgebra::DMatrix`, sparse | `leto::DMatrix`/`CsrMatrix` | unified linalg surface |
| `nalgebra::Point3/Vector3/Isometry` | `leto::Point3/Vector3/Isometry3` | same geometry across CPU/GPU |
| `tokio`, `rayon` | `moirai` executor + task graph | one runtime, async + parallel |
| `packed_simd`, manual unroll | `hermes-simd` lanes | portable SSE/AVX/NEON/WASM |
| `jemalloc`, `mimalloc` | `mnemosyne` allocator + `themis` placement | NUMA-aware, zero-fragmentation |
| `std::thread`, raw channels | `moirai::Executor` channels | typed task graph |
| `rustfft`, `realfft` | `apollo` forward FFT + autodiff | spectral methods w/ gradients |
| `num-traits`, `num-complex` | `eunomia::RealField/ComplexField` | single trait frontier |
| `burn::Tensor` (planned) | `coeus::Tensor<T>` + autograd | shared autodiff backends |

Each Atlas crate enforces the same design principles — **SRP** (one well-defined
surface per crate), **SSOT** (no overlap between crates), **DIP** (depend on
[eunomia] traits, not impls), **DRY** (no duplicated kernels across crates),
and zero-cost abstractions (ZSTs, `PhantomData`, const generics, GATs).

## Migration Status (CFDrs)

| Crate | Status | Notes |
|---|---|---|
| `cfd-core` | in-progress | scalar seams via [`eunomia::RealField`] |
| `cfd-math` | in-progress | dense paths via [`leto::NdArray`] |
| `cfd-1d` | legacy | targeted for `leto` + `moirai` first |
| `cfd-2d` | in-progress | spectral matrices ported to `leto::CsrMatrix` |
| `cfd-3d` | in-progress | FEM assembly kernel using `hermes-simd` |
| `cfd-validation` | partial | orthogonal oracle harness ported to eunomia traits |
| `cfd-schematics` | partial | CSG runtime uses `hephaestus` where available |
| `cfd-optim` | legacy | optimization audit behind stable facade |
| `cfd-io`, `cfd-python` | legacy | PyO3 boundary tied to ndarray interop |

## How To Read This Part

The sub-chapters are organized **from the lowest-level trait layer outward**.

1. [Eunomia: Numeric Traits](migration_eunomia.md) — the trait frontier every other crate depends on.
2. [Leto: Arrays and Linalg](migration_arrays.md) — CPU dense and sparse matrices, linear solvers.
3. [Leto: Geometry](migration_geometry.md) — points, vectors, isometries.
4. [Hermes: SIMD Lanes](migration_simd.md) — vectorized solvers.
5. [Mnemosyne and Themis: Memory](migration_memory.md) — arena allocation and NUMA placement.
6. [Moirai: Concurrency](migration_concurrency.md) — task graph, executors, channels.
7. [Apollo: FFT](migration_fft.md) — spectral transforms with autodiff integration.
8. [Leto: GAT Tiling](migration_gat_tiles.md) — lending iterators over data tiles.
9. [Migration Validation](migration_validation.md) — legacy-vs-atlas parity testing.

## Performance Contract

The Atlas stack promises the same per-flop throughput as the legacy crates,
but with strictly better **constants**:

- Zero heap allocations on hot paths (typed arenas via `mnemosyne`).
- Zero virtual dispatch on numeric kernels (trait specialization via `eunomia`).
- Zero copy for read-only views (`Cow` over `leto::NdArray`).
- Zero abstraction-tax polyfills (GATs encode lifetime variance).

Where Atlas does not yet match legacy throughput, this Part's
[Migration Validation](migration_validation.md) chapter records the gap and
governs the cleanup pass that closes it.
