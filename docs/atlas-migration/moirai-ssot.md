# Moirai parallel-surface SSOT formalization (kwavers-Atlas-migration-push)

This artifact declares the **Atlas-typed `moirai-parallel` macro-thread
/ task-rail** as the Single Source of Truth (SSOT) for all parallel
fan-out across the cfdrs workspace, replacing the legacy `rayon` +
`tokio` + `ndarray::Zip::indexed().par_*` fan-out surfaces. This
sub-batch formalizes the SSOT alignment under the
kwavers-Atlas-migration-push ceremony naming. cfdrs HEAD =
`1d768895 refactor(cfdrs): Atlas-provider migration push (Leto CSR +
Eunomia scalar + Hephaestus GPU + cfd-math / cfd-2d / cfd-3d /
cfd-1d / cfd-validation consumer cones)` is the trailing Atlas-typed
sweep.

## Surface inventory

### Legacy fan-out (NOT present in production)

- **0** call-sites of `rayon::` / `tokio::` / `par_iter` /
  `ndarray::Zip::indexed().par_*` in `repos/cfdrs/crates/`.
- The tokens appear only as string literals in
  `repos/cfdrs/xtask/src/migration_audit.rs`'s `LEGACY_SOURCE_TOKENS`
  const list (banned patterns), where the audit *detects* their
  presence or absence.

### Atlas-typed moirai adoption (29 call-sites)

| Surface pattern | Crates / files |
|---|---|
| `use moirai::prelude::{ParallelSlice, ParallelSliceMut}` | `cfd-core/src/physics/fluid_dynamics/operations.rs:8`, `cfd-math/src/simd/cfd.rs:13` |
| `use moirai::{fold_reduce_with, Adaptive}` | `cfd-3d/src/fem/solver.rs:43`, `cfd-optim/src/application/search/pool.rs:32`, `cfd-math/src/sparse/builder.rs:46`, `cfd-math/src/sparse/assembly.rs:7` |
| `use moirai::{map_collect_mut_with, Adaptive}` | `cfd-2d/src/network/solve.rs:15`, `cfd-2d/src/network/coupled.rs:15` |
| `use moirai::{reduce_index_with, Adaptive}` | `cfd-math/src/simd/vectorization.rs:19` |
| `use moirai::{map_collect_index_with, reduce_index_with, Adaptive}` | `cfd-math/src/simd/vector.rs:9` |
| `use moirai::ParallelSlice` | `cfd-optim/src/application/orchestration/milestone12/option2.rs:5` |
| `use moirai::prelude::ParallelSliceMut` | `cfd-math/src/simd/vectorization.rs:18` |

### Atlas co-resident dependencies

| Crate | Role |
|---|---|
| `hephaestus-wgpu` | GPU kernels (replaces direct `wgpu` in cfd-core/compute/gpu/{kernels,mod,pipeline}.rs) |
| `eunomia` | Typed scalar SSOT (replaces `num-traits`) |
| `leto` + `leto-ops` | CPU tensors + SVD/QR via Leto |
| `apollo-fft` | FFT (replaces `rustfft`) |
| `mnemosyne` | High-perf global allocator via `cfd-core/mnemosyne` feature |
| `hermes-simd` | SIMD abstraction |

## Validation

```
cargo run -p xtask -- migration-audit
```

Expects `LEGACY_FOUND: 0` and `LEGACY_DEP_TOKENS: 0` across all
production crates in the cfdrs workspace.

## Why moirai (cfdrs-specific rationale)

`cfd-core::mnemosyne` is the workspace's process `#[global_allocator]`.
mirroring `moirai`'s `no-global-alloc` feature avoids installing a
second allocator. moirai's `parallel` surface provides work-stealing
across workers and the `Adaptive` runtime-adaptive scheduler drives
most CFD fan-out sites. `ParallelSlice` / `ParallelSliceMut` extension
traits provide ndarray-class parallel iterators over `&[T]` and
`&mut [T]` slices with no extra coupling. The `fold_reduce_with`,
`map_collect_index_with`, and `reduce_index_with` primitives match
the CFD kernel shapes:

- FEM assemblers (`fold_reduce_with` over `row_offsets`)
- Sparse CSR builders (prefix-scan-based row counting)
- Womersley pulsatile 1D pipe sweeps
- Milestone 12 parameter-lattice optimization
- IBM delta-function kernels
- VOF scalar transport

Refs: kwavers-Atlas-migration-push ceremony; cfdrs sub-batch #1.
