# Appendix A — Atlas Crate Dependency Map

## Overview

This appendix is the **authoritative per-crate Atlas dependency map** for
CFDrs. It documents the Atlas stack surface, which CFDrs crate consumes which
Atlas crate, the Atlas-api boundary per consumer, the migration phase of each
consumer, the recommended migration order, and the cleanup gate that marks
completion. For the type-level migration (import rewrite table), see
[Migration Notes — ndarray/nalgebra/burn → leto/hephaestus/coeus](appendix_migration.md).

---

## Atlas Stack Summary

The Atlas monolith is a family of crates that share one numeric trait frontier
(`eunomia::RealField`) and one set of zero-cost memory abstractions
(`mnemosyne::Arena`, `themis::Placement`, `moirai::Executor`). The full family
(from the Atlas top-level book):

| Atlas crate | Surface | Formal role |
|---|---|---|
| `eunomia` | `RealField`, `ComplexField`, `FloatElement`, `IntElement` | Numeric trait unification for generics |
| `leto` (core) | `NdArray<T,D>`, `DMatrix`, `CsrMatrix`, `Point3`, `Vector3`, `Isometry3` | Dense/sparse arrays + geometry types |
| `leto` (GAT) | `LendingIterator`, `TileStreaming` | Zero-copy tiling with lifetime-parameterized iterators |
| `hermes-simd` | `SimdLane`, `Avx2Lane`, `NeonLane`, `WasmSimdLane` | Portable SIMD vectorization (SSE/AVX/NEON/WASM) |
| `mnemosyne` | `Arena`, `ScratchArena` | Bump/growing arenas, region-scoped allocation |
| `themis` | `Placement`, NUMA pinning | NUMA-aware memory placement, thread affinity |
| `moirai` | `Executor`, `TaskGraph`, channels | Unified async+parallel runtime replacing rayon/tokio |
| `apollo` | `FftPlan`, spectral autodiff bridge | FFT + autodiff dispatch (planned coeus tensor coupling) |
| `coeus` | `Tensor<T,D>` autograd | Automatic differentiation backend (not yet used by CFDrs) |
| `hephaestus` | `Backend`, GPU kernels, `DeviceBuffer` | GPU backend detection + storage (CUDA/ROCm/Metal) |
| `ritk` | Image I/O | not used by CFDrs |
| `consus` | Persistent storage / artifact store | planned for `cfd-io` |
| `tyche-core` | RngCore, Sobol', LHS | design-space sampling (via `cfd-optim`) |

---

## CFDrs-Facing Atlas Consumption Map

In CFDrs, seam-replacement of legacy crates ($\texttt{ndarray}$,
$\texttt{nalgebra}$, $\texttt{rayon}$, $\texttt{tokio}$, $\texttt{packed\_simd}$,
$\texttt{rustfft}$, $\texttt{num-traits}$) by Atlas is tracked per surface.

| Atlas crate | CFDrs surface consuming it | Legacy replaced | Compute vs I/O | Notes |
|---|---|---|---|---|
| `eunomia::RealField` | every CFDrs crate (generic bound on every solver/kernel) | `num-traits::Float`, `num-complex::Float` | compute | **COMPLETE** in all solver crates |
| `eunomia::FloatElement` | `cfd-core`, `cfd-1d`, `cfd-2d`, `cfd-3d`, `cfd-math`, `cfd-validation`, `cfd-schematics` | `num-traits` + ad hoc `Float+One+Zero+Sqrt` combos | compute | Convenience trait that composes `RealField+Default+Debug` |
| `leto::NdArray<T,D>` | `cfd-math` stencils + sparse ops, `cfd-1d`, `cfd-2d`, `cfd-3d`, `cfd-schematics` (sdf) | `ndarray::Array<T,D>`, `nalgebra::DMatrix` as dense | compute | Single dense type shared across all CFDrs |
| `leto::CsrMatrix` | `cfd-math::sparse`, `cfd-2d::pressure`, `cfd-3d::fem` | `nalgebra-sparse::CsrMatrix`, hand-rolled CSR | compute | Atlas CsrMatrix has builtin `apply` via `hermes-simd` |
| `leto::DMatrix` | `cfd-math::linalg`, `cfd-3d::spectral::diff_mat` (Chebyshev D) | `nalgebra::DMatrix` | compute | Alias over `NdArray<T, Ix2>` with BLAS dispatch |
| `leto::Point3/Vector3` | `cfd-schematics::primitives`, `cfd-3d::mesh` | `nalgebra::Point3/Vector3` | compute | Same storage but different trait frontier: `Copy`, `RealField`-generic |
| `leto::Isometry3` | `cfd-schematics` (mesh transforms) | `nalgebra::Isometry3` | compute | Transforms point+vector via dedicated API (no nalgebra left-multiply) |
| `leto` GAT (`LendingIterator`) | `cfd-2d::stencil`, `cfd-3d::stencil` (stencil tiles) | `Iterator<Item=&Tile>` workaround (GAT via nightly) | compute | Streaming zero-copy stencil tile views with GAT lifetime |
| `hermes-simd::Avx2Lane` | FEM quad loops, stencil inner loops, VOF compression | `packed_simd::f64x4`, `std::arch::x86_64` intrinsics | compute | AVX2/AVX-512/NEON/WASM per-target; scalar fallback on remaining |
| `mnemosyne::Arena` | every crate's per-step transient storage | `Vec<T>::with_capacity` + repeated realloc | compute | Bump-only; lifetime-scoped region — no dealloc churn |
| `themis::Placement` | `cfd-3d::mesh`, domain decomposition weights | raw `Vec<T>` + affinity guesswork | compute | NUMA first-touch placement hint for mesh blocks |
| `moirai::Executor` | `cfd-2d::parallel`, `cfd-3d::parallel`, `cfd-optim::batch` | `rayon::par_iter`, `tokio::spawn`, `std::thread::spawn` | compute | Scope-based tasks; cooperative async tasks for io-bound sim pipeline steps |
| `apollo::FftPlan` | `cfd-3d::spectral`, `cfd-2d::spectral` | `rustfft::FftPlanner`, `realfft::RealFftPlanner` | compute | **graduating**: Atlas FFT path performance within 1% of rustfft; switch is type-level flag |
| `hephaestus::Backend` | `cfd-schematics` GPU SDF + `cfd-3d` opt-in GPU transport | raw wgpu/cuda checks in examples | compute | Opt-in behind `gpu` feature — not yet mandatory path |
| `tyche-core` (Atlas rng + sampling) | `cfd-optim::design_space::sampling::*` | `rand`, `rand_chacha`, hand-rolled LHS | compute | Basis for LHS/Sobol sampling in audit/optim |
| `coeus::Tensor` (planned) | n/a | `burn::Tensor` (planned, not yet used) | compute | No autodiff surface in CFDrs yet |
| `consus` (planned) | `cfd-io` future HDF5/NetCDF | `hdf5`, `netcdf` crate deps | I/O | Planned — not yet started |
| `ritk` | n/a | n/a | I/O | Image I/O not needed |

---

## Per-Crate Migration Phase

Per-crate phase and still-open items:

| CFDrs crate | Phase | Completed surfaces | Remaining work | Atlas deps |
|---|---|---|---|---|
| `cfd-core` | **bulk** — finalize eunomia frontier | `eunomia::FloatElement` seam complete; `BoundaryKind<F>`, `StateVec<F,D>` generic over `F` | Finalize `leto::NdArray` conversion for fields that are still `ndarray::Array` via nightly GAT gating | `eunomia`, `mnemosyne` |
| `cfd-math` | **graduation** — NdArray migration final | `leto::NdArray` in new stencils; `apollo::FftPlan` path within 1% bandwidth of `rustfft` | Prune `ndarray`/`nalgebra` from `Cargo.toml`; delete nalgebra-based `DMatrix` wrapper | `eunomia`, `leto`, `hermes-simd`, `apollo` |
| `cfd-1d` | **bulk** — mid-migration | `eunomia` complete; rheology dispatch complete; Womersley exact path complete | Krylov laplacian graph `CsrMatrix` still ndarray-templated pending leto `CsrMatrix::from_graph` path; FL correction complete | `eunomia`, `leto`, `mnemosyne`, `moirai` |
| `cfd-2d` | **bulk** — mid-migration | Stencil ops trait-migrated (backward compat kept via feature); IBM masking complete | Finish `leto::NdArray<F,Ix2>` migration for ghost-cell buffers; `hermes-simd` lane in advection kernel; `moirai::Executor::scope` for parallel assembly | `eunomia`, `leto`, `mnemosyne`, `moirai`, `hermes-simd`, `apollo` |
| `cfd-3d` | **bulk** — mid-migration | FEM P1 tet assembly partially migrated; spectral Poisson with FFT decoupling via Atlas traits; `leto::Point3` behind flag | `hermes-simd` for quadrature inner loops; `moirai` parallel domain decomp; `hephaestus` opt-in transport | `eunomia`, `leto`, `mnemosyne`, `themis`, `moirai`, `hermes-simd`, `apollo`, opt. `hephaestus` |
| `cfd-validation` | **graduation** | eunomia trait harness complete; oracle tables independent of solver internals | Prune remaining legacy deps after solvers finish; add `cfd-3d` new oracle modules as graduated | `eunomia`, `leto` (output) |
| `cfd-schematics` | **bulk** | CSG is already mesh-agnostic; primitive SDF path uses leto types where available | Complete `leto::Point3/Vector3` conversion; `hephaestus` mesh extraction (GPU SDF→mesh) opt-in behind `gpu` feature | `leto`, `hephaestus` (GPU opt-in) |
| `cfd-io` | **legacy** — planned `consus` | reads/writes via HDF5/NetCDF/CSV libraries | Strongly coupled to `ndarray` view; wait until let o IO feature stabilizes; then wrap | (planned `consus`) |
| `cfd-optim` | **legacy / audit** | audit and GA logic legacy-complete | Retrofit batch evaluator `rayon::par_iter → moirai::Executor::scope` per `performance_and_atlas` | (planned `moirai`) |
| `cfd-python` | **legacy** — PyO3 thin wrapper | PyO3 boundary (thin `#[pyclass]` wrappers) | `ndarray`→`leto` inside PyO3 glue becomes transparent only after upstream crates fully switch; binding shim needed | PyO3 boundary |

---

## Migration Order (Recommended — Dependency-Ordered)

Nominates the sequence in which crates should graduate to minimize re-work:

1. **Complete `cfd-math` NdArray migration first** — it is the transitive
   dependency of all solver crates for dense and sparse storage. Until `cfd-math`
   exports `leto::NdArray` + `leto::CsrMatrix` as its primary API (not as a
   secondary impl behind a feature flag), all other crates still carry a legacy
   path for interop.

2. **Migrate `cfd-1d` next** — smallest surface among solvers; exercises the
   most-used Atlas APIs (`leto` arrays, `eunomia` traits, `mnemosyne` arena).
   No GPU path, simplifies first Atlas graduation. Validates the
   `CsrMatrix::from_graph` pattern for graph Laplacian construction.

3. **Migrate `cfd-2d` and `cfd-3d` together** — they share solver patterns
   (stencil, Krylov, SIMPLEC/PIMPLE, adaptive time-stepping). `cfd-2d` as the
   cheaper validation loop (shorter test runs); `cfd-3d` follows.

4. **Migrate `cfd-schematics` last** (after solver migrations are stable) — it
   depends on the new geometry types becoming consistent across all consumers.
   Its migration is shallow (point/vector type rename) but must synchronize with
   `cfd-3d::mesh` module's type declarations.

5. **After all solvers graduate**, run the `comprehensive_validation_suite`
   parity harness against legacy reference outputs, then prune the legacy crates
   (`ndarray`, `nalgebra`, `rayon`, `tokio`, `rustfft`) from every `Cargo.toml`.
   `cargo tree | grep -E '(ndarray|nalgebra|rayon|tokio|rustfft)'` must be
   empty for completion.

6. **`cfd-io` + `cfd-python`** migrate after solvers — their legacy is coupled
   to `ndarray` view types; wait for consistent `leto::NdArray` output from
   `cfd-2d`/`cfd-3d` so adapters become trivial re-writes.

7. **`cfd-optim`** migrates in parallel but with lower priority — it is already
   functionally complete, only its parallelism backend changes.

---

## Cleanup Gate — Markers That Migration Is Complete

Before legacy dependencies can be removed and a major version tag
incremented:

- [ ] `cargo tree | grep -E '(ndarray|nalgebra|rayon|tokio|rustfft|num-traits|packed_simd)'` is **empty** on the workspace `Cargo.lock`.
- [ ] `cfd-validation` parity passes on every canonical scenario: Ghia R100/400/1000,
      Hagen-Poiseuille, Richardson extrapolation $p\approx2$, turbulent channel RANS
      $U^+$ vs Moser DNS, venturi inception $\sigma_i$ vs Franc&Michel.
- [ ] `validation_reports/` updated with `migration_<crate>.md` (per-crate side-by-side
      before/after diff of validation results for that crate).
- [ ] No `unsafe { core::arch::x86_64::_mm...` } intrinsic outside `hermes-simd`
      crate — verified by grepping `crates/*/src/` (maintenance guard).
- [ ] No `unsafe { Vec::from_raw_parts ... }` outside `mnemosyne` — same grep guard.
- [ ] `cfd-io` and `cfd-python` consuming `leto::NdArray` without `ndarray` shim.
- [ ] `cfd-optim` batch evaluator via `moirai::Executor::scope`.
- [ ] Dedicated bench harness (`criterion`) shows Atlas path no slower than legacy
      path on canonical cases (cavity $Re=400$, pipe parabolic, spectral Poisson).
- [ ] No `#[cfg(legacy)]` / `#[cfg(atlas)]` feature gates left in shipped code —
      every dual-path has collapsed to Atlas-only with legacy path removed.

---

## Further Reading

- [Migration Notes — ndarray/nalgebra/burn → leto/hephaestus/coeus](appendix_migration.md) — type alias mapping + recipe.
- [Migration Overview](migration_overview.md) — Atlas motivation and per-layer forward roadmap.
- [Performance and Atlas Integration Overview](performance_and_atlas.md) — performance contract and migration table.
- [Atlas Book — top-level book](../../../../README.md) — Atlas crate inventory and design rationales.
- [`BOOK_ORGANIZATION.md`](BOOK_ORGANIZATION.md) — forward roadmap phases.
