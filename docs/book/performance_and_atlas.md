# Chapter 7 — Performance and Atlas Integration Overview

This Part is the **Atlas stack migration reference** for CFDrs.  It documents
how the legacy third-party crates have been — or are being — replaced
by the Atlas monolith: the unified, purpose-built crates that share one
trait frontier ([`eunomia::RealField`]) and one set of zero-cost
abstractions.

## Performance Mandate

Atlas is held to a strict performance contract:

- **Zero heap allocations on hot paths.** All transient storage lives in
  [`mnemosyne::Arena`] sub-arenas.
- **Zero virtual dispatch on numeric kernels.** Trait specialization
  through [`eunomia`] gives the compiler monomorphized code per
  scalar type and backend.
- **Zero copy for read-only views.** [`Cow`] over [`leto::NdArray`]
  eliminates `clone()` calls in adjoint passes.
- **Zero abstraction tax.** Const generics, ZSTs, and GATs encode
  capability and lifetime variance at the type level — there is no
  runtime machinery to charge.

Where Atlas does not yet match legacy throughput, the gap is recorded in
the per-crate migration chapters of this Part and governed by the parity
test harness described in [Migration Validation](migration_validation.md).

## Atlas Crate Layers

| Layer | Crate(s) | CFDrs consumer |
|---|---|---|
| Numeric traits | `eunomia::RealField` | every type |
| Dense/sparse storage | `leto::NdArray`, `leto::CsrMatrix`, `leto::DMatrix` | `cfd-math`, `cfd-1d`, `cfd-2d`, `cfd-3d` |
| Geometry | `leto::Point3`, `Vector3`, `Isometry3` | `cfd-schematics`, `cfd-3d` |
| SIMD | `hermes-simd` lanes | `cfd-2d::stencil`, `cfd-3d::fem`, `cfd-1d::kernel` |
| Memory | `mnemosyne::Arena`, `themis::Placement` | every crate |
| Concurrency | `moirai::Executor`, task graph | `cfd-2d::parallel`, `cfd-3d::parallel` |
| FFT | `apollo::FftPlan` | `cfd-3d::spectral`, `cfd-2d::spectral` |
| Tiling | `leto` GAT-based `TileStreaming` | `cfd-3d::stencil`, `cfd-2d::stencil` |

(Tensors, autodiff, GPU backends, image I/O — included in the full Atlas
stack — are used by `helios` and `kwavers` but not by CFDrs; see
[Atlas Dependency Map](appendix_dependencies.md) for surface compatibility.)

## Migration Table

| CFDrs crate | Legacy → Atlas | Status |
|---|---|---|
| `cfd-core` | `ndarray::ArrayBase` → `leto` traits | in-progress |
| `cfd-math` | `ndarray::Array2`/`nalgebra::DMatrix` → `leto::NdArray` | in-progress |
| `cfd-math` SpMV | `csr` raw → `leto::CsrMatrix` + `hermes` lanes | in-progress |
| `cfd-math` FFT | `rustfft` → `apollo::FftPlan` | graduating |
| `cfd-1d` | bespoke `Vec<f64>` → `Arena` + gutters | mid-flight |
| `cfd-2d` | `Array3` → `NdArray<f64, Ix3>` | in-progress |
| `cfd-3d` | `Array3` + explicit SIMD → `NdArray` + `hermes` | mid-flight |
| `cfd-validation` | orthogonal oracle harness | complete (eunomia traits) |
| `cfd-schematics` | CSG primitive catalogue | in-progress |
| `cfd-optim` | optimization audit harness | legacy (audit-only) |
| `cfd-io`, `cfd-python` | via `ndarray` interop | legacy |

See [BOOK_ORGANIZATION.md](BOOK_ORGANIZATION.md) for the four-phase
forward roadmap.

## Examples Referenced by This Chapter

- [`simd_performance_benchmark`](examples/simd_performance_benchmark.md) —
  Hermes SIMD vs legacy scalar Stencil benchmarks.
- [`gpu_detection`](examples/gpu_detection.md) — Atlas GPU detection
  (lets CFDrs flow into `hephaestus` once that's wired for CFD solvers).
- [`spectral_performance`](examples/spectral_performance.md) —
  Apollo FFT throughput vs `rustfft`.
- [`venturi_validated`](examples/venturi_validated.md) — venturi case
  validated against reference venturi phantom under both legacy and Atlas paths.

## How To Use This Part

Read [Migration Overview](migration_overview.md) first; then go to the
per-crate sub-chapters in the order that matches the layer your code
touches.  When closing a migration gap, the
[Migration Validation](migration_validation.md) chapter governs the
parity test that promotes the migration from *bulk* to *cleanup*.

## Further Reading

- [Atlas Dependency Map](appendix_dependencies.md)
- [Atlas Book (cross-reference)](../../../../README.md)
- [`BOOK_ORGANIZATION.md`](BOOK_ORGANIZATION.md) — forward roadmap
