> ## Vocabulary policy (canonical atlas-migration terms-of-art)
>
> **Canonical functional terms-of-art (preserve)**:
> - `Atlas-typed` (the CoeUs/MoiraiBackend-typed twin type-system family that pairs with `Burn-keyed` as the atomic-boundary partition term per ADR 0012 §Decision §1)
> - `Atlas-side` (the additive production-side / subtractive test-side partition — pairs with `Burn-side`)
> - `Atlas-only` (the validation-gate constraint term, e.g. "Atlas-only backend trait assertion")
> - `Atlas-meta` (the atlas-meta repo/branch identifier, codex/kwavers-atlas-integration)
> - `Atlas-native` (a label for modules/edges that route through native CoeUs/Eunomia/Leto without burn-compat shims)
> - `Atlas-backed` (a run-time / compile-time Atlas-runtime carrier — Atlas-runtime kernels, Codex-Atlas pipelines, etc.)
> - `Atlas migration push` / `Atlas migration` (the pregnant noun phrase used to describe the ongoing per-crate migration ceremonies — Atlas-migration makes the chronology explicit without over-decoration)
> - `cfd-* Atlas-typed` /
>   `cfd-math` / `cfd-1d` / `cfd-2d` / `cfd-3d` /
>   `cfd-core` / `cfd-io` / `cfd-schematics` / `cfd-validation` /
>   `cfd-python` / `cfd-optim` (the CFDrs inner-crate atomic-boundary partition terms that pair with the Atlas-migration-push ceremony list in the canonical Atlas-meta rubric)
>
> **Discouraged over-decorated compounds (drop)**:
> - `Atlas-parent` -> collapse to `atlas-meta`
> - `Atlas-root` -> collapse to `atlas-meta working tree`
> - `Atlas-provider` -> collapse to bare `Atlas` when it modifies a non-push noun (e.g. "Atlas-provider boundaries/work/slices/residue"), or `Atlas migration push` when it modifies a ceremony noun (e.g. "Atlas-provider migration push")
>
> **CFDrs-specific guidance**:
> - The CFDrs migration covers `cfd-math` + `cfd-1d` + `cfd-2d` + `cfd-3d` + `cfd-core` + `cfd-io` + `cfd-schematics` + `cfd-validation` + `cfd-python` + `cfd-optim` (per the canonical Atlas-meta rubric + the inner refactor commit history). Per-crate case-study references should preserve the inner-crate-name as a compound prefix without dropping the `Atlas-` prefix.
> - The CFDrs cross-crate edges use `Atlas-native` (no burn-compat) for the inner-crate heap and `Atlas-backed` for the `cfd-validate` / `cfd-optim` consumer cones that import MoiraiAtlas-runtime kernels. The Atlas-migration-push ceremony chain lives in the `Atlas migration / Atlas migration push` canonical form.
>
> Mirror reference: atlas-meta backlog.md / checklist.md / gap_audit.md + repos/ritk/{CHANGELOG.md, checklist.md, gap_audit.md} (same six canonical + three disallowed compounds in the same one-page rubric form).
# Changelog

- Removed stale legacy-migration audit exemptions for the fully Eunomia/Leto
  `cfd-1d` and `cfd-3d` scalar seams.

All notable changes to this project will be documented in this file.

## [Unreleased]

### Changed

- **Breaking / architecture**: schematic analysis overlays now consume Iris
  `NamedColorMap` values directly, accept borrowed or owned scalar maps through
  `Cow`, validate finite fields at construction, and precompute scalar ranges.
  The duplicate `ColormapKind` and local blue-red, Viridis, and grayscale laws
  are removed. Builders now return typed `Result`s, and scalar-map fields are
  exposed through borrowed accessors instead of public mutable storage.
- **Breaking / architecture**: `cfd-math::LaplacianOperator2D::new` now accepts
  Aequitas lengths and an explicit Leto boundary condition, returns a typed
  `Result`, and delegates CPU evaluation to `leto-ops::laplacian_2d_into`.
  The GPU solver operator selects the same negative-Laplacian polarity. Local
  CPU formulas and the duplicate cfd-core test oracle are removed.
- `cfd-optim` now exposes Milestone 12 Latin-hypercube candidate generation
  through Tyche's const-generic, counter-addressed design. The consumer maps
  fixed unit-sample arrays directly into CFD parameters; its duplicate
  mutable-RNG sampler and direct `rand` dependency are removed. The provider
  closure advances Aequitas and Proteus together so CFDrs, Proteus, and
  Hephaestus share one dimensional-type identity.
- **Breaking / architecture**: `PolynomialViscosity::calculate_density` now
  returns a typed `Result` and delegates its linear thermal-expansion response
  to Proteus. The response remains generic over the native scalar, while
  invalid reference states, coefficients, temperatures, and evaluated
  densities are rejected at the shared material-law boundary.
- **Breaking / architecture**: `FluidState::thermal_diffusivity` now returns a
  typed `Result`, and both fluid-property entry points delegate dimensional
  validation and `alpha = k / (rho c_p)` to Proteus. CFDrs retains rheology and
  flow closure while the shared thermophysical law has one Atlas owner.
- **Breaking / architecture**: GPU Laplacian spacing now uses Aequitas
  `Length<f32>` through the CFDrs facade and Hephaestus provider boundary.
  Raw metre comments no longer carry the dimensional contract; unit conversion
  is resolved once before dispatch.

- Preserved the independent rsparse sparse-LU fallback after rejecting a stale
  unpreconditioned-GMRES substitution. Leto sparse direct factorization remains
  an explicit upstream requirement rather than a consumer-side semantic
  downgrade.
- Replaced root serpentine and 2D bifurcation validation examples with the
  canonical cfd-2d analytical and discretized solver paths. The 1D blood-flow
  example now validates its pressure-balanced split and Fåhræus-Lindqvist
  trend from calculated values.

### Removed

- Removed the unreferenced comprehensive, wall-shear, and 3D FEM root reports
  because they presented static or invalid validation values rather than
  successful solver results.

## [0.3.0] - 2026-07-17 - Typed GPU Capability Boundary

### Changed

- `GpuContext` acquires its provider through Hephaestus'
  `ComputeDeviceAcquisition` contract and reports capabilities through
  `ComputeDeviceCapabilities`. The derived seven-storage-binding requirement
  retains the provider's full downlevel compatibility baseline.
- Nextest serializes only WGPU provider-acquiring tests through the
  `gpu-device` group; CPU-only tests remain concurrent.

### Breaking

- `GpuContext::adapter_info` and `GpuContext::supports_features` are removed.
  Consumers use `backend_name`, `supports_timestamp_queries`,
  `max_work_group_size`, and `max_buffer_size`; new provider capabilities land
  in Hephaestus rather than exposing WGPU types through CFDrs.

## [0.2.0] - 2026-07-17 - Error Consolidation

### Changed

- `cfd-core::compute::gpu::GpuContext::synchronize` now uses Hephaestus'
  backend-neutral `ComputeDevice` completion contract rather than polling the
  underlying WGPU device directly.
- Updated CFDrs to merged Moirai `main` and regenerated the locked Atlas
  provider graph. The resolved workspace uses Moirai 0.4 with Themis 0.10.
- Advanced the direct Leto and Leto-ops source pins to merged `main`.
- Replaced the serpentine-Venturi demo's untyped configuration tuple with a
  named geometry record and repaired warning-denied `cfd-schematics`
  test/example diagnostics.

### Breaking

- `GpuContext` no longer exposes raw WGPU device, queue, or limit fields, and
  `GpuPoissonSolver::new(device, queue, ...)` is removed. Construct the solver
  with `GpuPoissonSolver::from_context(&context, ...)` so Hephaestus retains
  provider ownership. `GpuBuffer` no longer exposes its raw buffer handle or
  a `buffer` accessor.
- Removed `spmv_parallel`, `IterativeSolverConfig::use_parallel_spmv`,
  `IterativeSolverConfig::with_parallel_spmv`, and
  `MomentumSolver::with_parallel_spmv`. Callers use `spmv`, `try_spmv`, and
  `MomentumSolver::new`; Leto owns the CSR execution policy.
- `ComputeDispatcher::execute` no longer silently executes unsupported or
  failed GPU/SIMD work on CPU. Unsupported kernel/backend pairings and
  unavailable requested providers now return typed `UnsupportedOperation`
  errors; GPU execution errors propagate to the caller.

### Summary

Consolidated domain error types into `cfd_core::error` with 15 Kind enums and
29 Error variants when GPU support is enabled. `cfd-io` now owns a local
file-format error type so the I/O
crate does not depend on `cfd-core` only for error reporting while the provider
migration removes nalgebra from its normal dependency graph. Only `cfd-optim`
retains a local `OptimError` (bridged via `From` impl). Removed 8 dead extension
traits (~323 lines). Fixed ~100 rustdoc warnings across 161 files.

### Breaking
- **cfd-core generic GPU pipeline**: Removed the unused `GpuPipelineManager`,
  `GpuKernel<T>`, and `GpuContext::create_compute_pipeline_with_layout` raw-WGPU
  surfaces. Consumers define operations through Hephaestus typed kernels.
- **cfd-core/cfd-2d GPU turbulence**: `GpuTurbulenceCompute` now writes native
  f32 results into caller-owned slices using `TurbulenceGrid`. Removed public
  raw kernel accessors, GPU-buffer return/readback methods,
  `TurbulencePerformanceInfo`, and the cfd-2d f64 Smagorinsky `use_gpu` path.
- **cfd-core GPU pressure**: Replaced the non-executing generic
  `GpuPressureKernel<T>` with real `f32` Hephaestus weighted-Jacobi and residual
  operations constructed from a `GpuContext`. Callers now provide
  `PressureConfig`, pressure/source slices, an output slice, and handle
  `Result`.
- **cfd-core GPU velocity**: Replaced the non-executing generic
  `GpuVelocityKernel<T>` with real `f32` Hephaestus correction and
  divergence-source operations constructed from a `GpuContext`. Callers now
  provide `VelocityConfig`, component slices, output slices, and handle
  `Result`.
- **cfd-core GPU diffusion**: Replaced the non-executing generic
  `GpuDiffusionKernel<T>` with a real `f32` Hephaestus kernel constructed from
  a `GpuContext`. Callers now provide `DiffusionConfig`, input/output slices,
  and handle `Result`.
- **cfd-core GPU advection**: Replaced the non-executing generic
  `GpuAdvectionKernel<T>` with a real `f32` Hephaestus kernel constructed from a
  `GpuContext`. Callers now provide `AdvectionConfig`, scalar/velocity fields,
  an output slice, and handle `Result`.
- **cfd-core/cfd-math GPU Laplacian**: `Laplacian2DKernel::new`, host execution,
  and `GpuFieldOps::new` are now fallible. `GpuLaplacianOperator2D` exposes the
  real `f32` WGSL precision instead of a generic scalar parameter.
- **cfd-core GPU arithmetic**: Removed the transitional public
  `FieldAddKernel` and `FieldMulKernel` types. `GpuFieldOps::add_fields` and
  `GpuFieldOps::multiply_field` now return `cfd_core::error::Result<()>` and
  propagate typed Hephaestus failures instead of silently recomputing on CPU.

### Migration
- Handle `ComputeDispatcher::execute` errors when a kernel does not implement
  the selected backend. Select an explicit CPU dispatcher for CPU execution
  instead of relying on provider downgrade behavior.
- Delete generic pipeline registration and raw `GpuKernel<T>` implementations.
  Construct the appropriate operation-specific Hephaestus-backed CFD facade and
  call its typed method.
- Construct `TurbulenceGrid`, allocate an output slice, and call
  `GpuTurbulenceCompute::{compute_smagorinsky_sgs,compute_des_length_scale,
  compute_wall_distance}`. DES grid scale no longer accepts unused velocity
  fields. The f64 cfd-2d Smagorinsky model executes in native precision on CPU;
  use the cfd-core f32 facade only for f32 domains.
- Replace `GpuPressureKernel::<T>::new()` and the raw `GpuKernel` trait path
  with `GpuPressureKernel::new(context)?`, then call `iterate(...)` or
  `residual(...)` with `PressureConfig`.
- Replace `GpuVelocityKernel::<T>::new()` and the raw `GpuKernel` trait path
  with `GpuVelocityKernel::new(context)?`, then call `correct(...)` or
  `divergence_source(...)` with `VelocityConfig`.
- Replace `GpuDiffusionKernel::<T>::new()` and the raw `GpuKernel` trait path
  with `GpuDiffusionKernel::new(context)?` followed by
  `execute(input, config, output)?`.
- Replace `GpuAdvectionKernel::<T>::new()` and the raw `GpuKernel` trait path
  with `GpuAdvectionKernel::new(context)?` followed by
  `execute(scalar, velocity_x, velocity_y, config, output)?`.
- Handle the `Result` returned by Laplacian kernel/facade construction and
  execution. Use the CPU operator for non-`f32` scalar domains until a native
  provider kernel for that precision exists.
- Construct `GpuFieldOps` as before, but handle the `Result` returned by
  `add_fields` and `multiply_field`. Consumers that instantiated the deleted
  kernel types should use the corresponding `GpuFieldOps` method.

### Changed
- **cfd-core**: Removed residual CPU fallback from runtime compute dispatch and
  added a regression proving unsupported backend requests do not mutate output
  through hidden CPU execution.
- **cfd-core**: Deleted 319 lines of unconsumed raw-WGPU pipeline registry,
  bind-group, uniform reconstruction, dispatch, context helper, and obsolete
  trait code. Static audit finds no CFDrs-owned shader or compute-pipeline
  creation under `compute::gpu`; all core package gates pass.
- **cfd-core/cfd-2d**: Routed Smagorinsky viscosity, DES grid cutoff, and
  rectangular wall distance through Hephaestus typed multi-storage kernels.
  Consolidated duplicate raw pipelines and buffer caches, exposed the formerly
  unreachable wall-distance operation, removed unused DES velocity inputs,
  removed fabricated speedup estimation and a wall-clock threshold test, and
  deleted the precision-changing cfd-2d GPU bridge. Focused tests pass 4/4,
  full core passes 243/243, and full cfd-2d passes 570/570; clippy, doctests,
  benchmark compilation, and static audits pass.
- **cfd-core**: Routed weighted-Jacobi pressure iteration and pointwise
  residual evaluation through separate Hephaestus typed multi-storage kernels.
  Added a validating grid/relaxation contract, corrected Neumann edge/corner
  handling, consolidated shared 3D dispatch sizing, moved the family to
  `kernels/pressure/{mod,kernel,tests}`, and replaced name/shader-only coverage
  with exact quadratic Poisson tests. Focused tests pass 6/6 and full core tests
  pass 247/247; checks, clippy, doctests, docs, and static audits are clean.
- **cfd-core**: Routed SIMPLE velocity correction and pressure-source
  divergence through separate Hephaestus typed multi-storage kernels. Added a
  validating grid/physical-coefficient contract, requested the derived
  seven-storage-buffer device capability during provider acquisition, moved
  the family to `kernels/velocity/{mod,kernel,tests}`, and replaced
  name/shader-only coverage with exact linear pressure/velocity field tests.
  Focused tests pass 5/5 and full core tests pass 242/242; checks, clippy,
  doctests, docs, and static migration/provider audits are clean.
- **cfd-core**: Routed explicit three-dimensional central diffusion through
  Hephaestus typed multi-storage dispatch. Added a validating grid, physical
  coefficient, and von Neumann stability contract; moved the family to
  `kernels/diffusion/{mod,kernel,tests}` with an aligned Rust/WGSL parameter
  layout; and replaced name/complexity-only behavior with exact constant and
  quadratic-field GPU tests. Focused tests pass 4/4 and full core tests pass
  238/238; checks, clippy, doctests, docs, and static migration/provider audits
  are clean.
- **cfd-core**: Routed first-order upwind advection through Hephaestus typed
  multi-storage dispatch. Added a validating grid/timestep contract with finite
  input and CFL enforcement, moved the family to
  `kernels/advection/{mod,kernel,tests}`, colocated the WGSL source, and replaced
  name/complexity-only coverage with exact GPU value tests. Focused tests pass
  6/6 and full core tests pass 234/234; checks, clippy, doctests, docs, and
  static migration/provider audits are clean.
- **cfd-core/cfd-math**: Routed the 2D GPU Laplacian through Hephaestus typed
  multi-storage dispatch and deleted raw WGPU pipeline, bind-group, staging,
  polling, timeout, and silent CPU-fallback code. Corrected endpoint-inclusive
  periodic wrapping, added typed contract validation, retained the CPU stencil
  only as a differential oracle, and removed the false downstream scalar
  generic. Focused tests pass 10/10; full core/math tests pass 231/231 and
  362/362; clippy, doctests, and docs are clean.
- **cfd-core**: Routed GPU field addition and scalar multiplication through
  Hephaestus typed, cached elementwise kernels. Deleted the 449-line duplicate
  raw WGPU pipeline, WGSL, staging, polling, timeout, and fallback code; split
  the field-operation facade into arithmetic, Laplacian, and test leaves.
  Exact partial-workgroup and typed-error regressions pass, as do no-default
  and GPU checks, all-target clippy, 230/230 nextest, 5/5 doctests, and
  warning-clean docs.
- **cfd-1d/cfd-3d**: Removed direct `num-traits` dependency ownership from the
  crate scalar seams. The workspace dependency catalog and both crate
  manifests no longer declare `num-traits`, and `Cfd1dScalar`/`Cfd3dScalar`
  now expose `zero()`/`one()` through Eunomia `NumericElement` constants
  instead of `num_traits::{Zero,One}` bounds. Evidence: touched-file rustfmt,
  `cargo check -p cfd-1d -p cfd-3d`, targeted direct-residue scan, and
  dual-package nextest 1122/1122 with one existing slow 3D mesh-convergence
  validation. Residual: package-wide fmt/clippy are still blocked by unrelated
  baseline formatting/lint drift outside this scalar dependency slice.
- **cfd-1d**: Moved solver-core reusable workspace vectors from nalgebra
  `DVector` to Leto `Array1`. `SolverWorkspace` now stores RHS, previous
  solution, and linear initial-guess/output buffers as Leto arrays;
  `MatrixAssembler::assemble`, system validation, residual checks, and
  iterative solver calls consume Leto RHS/guess storage directly. A private
  `vector_bridge` module now owns the remaining nalgebra<->Leto conversions
  at the current nalgebra-returning linear solver and Anderson Picard solution
  boundary. Evidence: cfd-1d test-target check, cfd-1d lib clippy, focused
  solver-core nextest (10/10), and cfd-1d doc generation completed with the
  existing unrelated rustdoc warnings. Residual: `linear_system.rs` still
  exposes nalgebra solution/dense fallback and nalgebra-sparse matrix input
  boundaries; `matrix_assembly.rs` still returns nalgebra-sparse CSR; solver
  detection and Anderson acceleration still operate on nalgebra solution
  vectors at that boundary.
- **cfd-1d**: Moved `NetworkState` public pressure/flow storage from nalgebra
  `DVector` to Leto `Array1`. State construction now populates Leto arrays
  directly from network pressure and flow slices, and focused unit tests assert
  shape, value, clone, and time semantics. Evidence: cfd-1d test-target check,
  cfd-1d lib clippy, focused state nextest (2/2), and cfd-1d doc generation
  completed with the existing unrelated rustdoc warnings; targeted state
  residue scan found no `DVector`. The upstream Atlas blocker was closed in
  Eunomia by adding `FloatElement::acos` with native f64 primitive/wrapper
  overrides. Residual: broader cfd-1d linear-system matrix/solution,
  sparse-assembly, solver-detection, and Anderson solution-vector boundaries
  still retain nalgebra surfaces.
- **cfd-1d**: Moved the solver-core convergence checker public vector boundary
  from nalgebra `DVector` to Leto `Array1`. `ConvergenceChecker` now checks
  finite values, solution-change convergence, and dual residual/change
  convergence over Leto arrays with direct L2 computations and typed
  mismatched-length errors. Evidence: cfd-1d test-target check, focused
  cfd-1d clippy for lib plus `solver_core_tests`, focused convergence nextest
  20/20, and targeted residue scans. Residual: the broader 1D solver
  linear-system matrix/solution, sparse-assembly, solver-detection, and
  Anderson solution-vector boundaries still use nalgebra vectors/matrices and
  nalgebra-sparse, and cfd-1d rustdoc still has unrelated link warnings.
- **cfd-suite**: Removed direct `num-traits` and `simba` dependencies from the
  root package and workspace dependency catalog. CFDrs source and member
  manifests now route scalar contracts through Eunomia/Atlas providers rather
  than direct legacy numeric dependencies. Evidence: targeted
  source/member-manifest residue scan with no `num_traits`, `simba`, direct
  `num-traits.workspace`, or direct `simba.workspace` matches outside the
  removed root entries, and root package metadata. Residual: `num-traits` and
  `simba` remain in `Cargo.lock` as transitive dependencies of upstream crates;
  root package check/nextest were blocked by concurrent shared Cargo
  cache/build locks and did not run.
- **cfd-suite/cfd-3d**: Removed direct `crossbeam` ownership from the
  workspace dependency catalog and `cfd-3d` manifest. Current CFDrs source and
  member manifests now rely on Moirai-facing concurrency surfaces with no
  direct Crossbeam imports. Evidence: targeted source/member-manifest residue
  scan, root package metadata, `cfd-3d` no-default-features check,
  `cfd-3d` no-default-features nextest (394/394), and diff whitespace check.
  Residual: this slice removes direct CFDrs ownership only; any transitive
  Crossbeam dependencies remain owned by upstream crates. The focused nextest
  run marked existing mesh-convergence validation slow at 19.7s.
- **cfd-math/cfd-validation**: Closed the remaining public sparse and
  linear-solver Atlas boundary residual. `cfd_math::sparse::SparseMatrix<T>`
  now aliases `leto_ops::CsrMatrix<T>`, sparse builders/operations construct
  and operate on Leto CSR directly, and the public `LinearOperator`,
  `Preconditioner`, `LinearSolver`, direct-solver, and solver-chain cone no
  longer exposes nalgebra `DVector` or `nalgebra_sparse::CsrMatrix`.
  `cfd-validation::numerical` now stores validation solutions as
  `leto::Array1<T>` and builds sparse test cases through the Leto CSR
  boundary. Evidence: cfd-math check, cfd-math test-target check, cfd-math
  all-target clippy, cfd-math doc, cfd-math nextest 361/361, cfd-validation
  check, cfd-validation all-target clippy, cfd-validation doc, and targeted
  residue scans. Residual: cfd-validation nextest still fails in two existing
  venturi cross-fidelity convergence tests outside this boundary.
- **cfd-math**: Moved IncompleteCholesky construction and factor storage to
  Leto CSR. `cholesky.rs` now stores `leto_ops::CsrMatrix`, performs
  symmetry/factorization reads through Leto row-value lookup, builds IC(0)
  factors with `CsrMatrix::from_parts`, and applies triangular substitutions
  through Leto CSR row views. Evidence: cfd-math fmt, lib check, all-target
  check, lib/tests clippy, all-target clippy, cholesky-filter nextest (5/5),
  preconditioner nextest (76/76), and targeted residue scans. Residual:
  Schwarz, direct solver, remaining transitional solver fixtures, and the
  shared solver sparse matrix boundary still use nalgebra sparse.
- **cfd-math**: Moved the ILU preconditioner family to Leto CSR. ILU(0),
  ILU(k), and triangular solve application now use `leto_ops::CsrMatrix`.
  Source tests, integration edge tests, `LinearSolverChain`, and Schwarz local
  solves now pass Leto CSR into ILU or convert once at the remaining shared
  solver/Schwarz boundary. Evidence: cfd-math fmt, lib check, all-target
  check, lib/tests clippy, all-target clippy, ilu-filter nextest (21/21),
  preconditioner nextest (76/76), linear_solver::tests nextest (53/53), and
  targeted residue scans. Residual: Schwarz, direct solver, remaining
  transitional solver fixtures, and the shared solver sparse matrix boundary
  still use nalgebra sparse.
- **cfd-math**: Moved SSOR preconditioner construction to Leto CSR. SSOR now
  owns `leto_ops::CsrMatrix`, uses Leto CSR row views for forward/backward
  sweeps, and applies through Leto `Array1` buffers with explicit length
  validation. Source preconditioner edge tests convert their still-transitional
  solver CSR fixtures once before constructing SSOR and assert exact typed
  errors for vector-length mismatches. Evidence: cfd-math fmt, lib check,
  all-target check, lib/tests clippy, all-target clippy, ssor-filter nextest
  (5/5), preconditioner nextest (76/76), and targeted residue scans.
  Residual: Schwarz, direct solver, integration-test fixtures,
  and the shared solver sparse matrix boundary still use nalgebra sparse.
- **cfd-math**: Moved basic preconditioner construction to Leto CSR. Jacobi
  and SOR now take `leto_ops::CsrMatrix`, use Leto CSR diagonal/row access,
  and apply through Leto `Array1` buffers with explicit length validation.
  Source linear-solver tests and `core_solver_tests.rs` convert their
  still-transitional solver CSR fixtures once before constructing Jacobi/SOR.
  Evidence: cfd-math fmt, lib check, core solver test check, all-target check,
  lib/tests clippy, core solver test clippy, all-target clippy,
  linear_solver::tests nextest (53/53), core_solver_tests nextest (4/4),
  preconditioner nextest (76/76), and targeted residue scans. Residual:
  Schwarz, direct solver, and the shared solver sparse matrix boundary still
  use nalgebra sparse.
- **cfd-math**: Moved the AMG/coarsening sparse boundary to Leto CSR. The
  multigrid bounded context now uses `leto_ops::CsrMatrix` for AMG setup,
  coarsening, interpolation, smoothers, cycles, sparse products, transpose,
  and SpMV. `coarsening_bench.rs`, `algebraic_distance_bench.rs`,
  `amg_coarsening_tests.rs`, and AMG integration preconditioner construction
  now use Leto CSR. Evidence: cfd-math fmt, lib check, focused AMG/coarsening
  checks, focused clippy for lib/tests/benches, AMG integration nextest (5/5),
  AMG-filter nextest (6/6), multigrid::coarsening nextest (10/10), and a
  targeted AMG/coarsening sparse residue scan. Residual: the broader
  solver/direct/preconditioner sparse matrix boundary still uses nalgebra
  sparse, with `LinearSolverChain` converting once before AMG construction.
- **cfd-math**: Migrated SpMV/CG benchmarks to direct Leto CSR. `benches/
  spmv_bench.rs`, `benches/cg_bench.rs`, and the CG section of `benches/
  math_benchmarks.rs` now construct `leto_ops::CsrMatrix` values directly with
  `from_parts`, and the SpMV benchmark now measures the direct Leto CSR
  `LinearOperator::apply` path instead of the legacy nalgebra sparse helper.
  Evidence: cfd-math fmt, focused bench checks, focused bench clippy,
  all-target check, all-target clippy, sparse-filter nextest (19/19), and a
  targeted migrated-benchmark residue scan. Residual after the AMG/coarsening
  boundary moved to Leto CSR: the broader solver/direct/preconditioner sparse
  matrix surface still uses nalgebra sparse.
- **cfd-math**: Added a direct Leto CSR linear-operator path. `src/sparse/
  operations.rs` now implements `LinearOperator<T>` for
  `leto_ops::CsrMatrix<T>` and shares one `try_leto_spmv` execution helper
  between direct Leto CSR inputs and the remaining nalgebra CSR conversion
  path. `tests/simple_gmres_tests.rs` now constructs Leto CSR matrices
  directly and verifies GMRES residuals through that operator path, with no
  nalgebra sparse/vector residue in the migrated test. Evidence: cfd-math fmt,
  simple_gmres test check, simple_gmres nextest (3/3), simple_gmres clippy,
  cfd-math lib check, all-target check, all-target clippy, sparse-filter
  nextest (19/19), gmres-filter nextest (21/21), and a targeted simple GMRES
  residue scan. Residual: broader cfd-math sparse builders, preconditioners,
  AMG, direct solver, tests, and benches still expose the nalgebra-sparse
  matrix boundary.
- **cfd-math**: Closed the remaining Leto storage-slice residue in source and
  tests. Nonlinear dense pivoting now swaps `Array1`/`Array2` entries by
  direct indexing; AMG interpolation quality validation indexes interpolated
  Leto arrays directly; and multigrid smoother tests assert indexed vector
  values instead of borrowing storage slices. Evidence: cfd-math fmt, lib
  check, focused nonlinear_solver/multigrid nextest (46/46), lib/tests clippy,
  all-target check, all-target clippy, and a cfd-math `src`/`tests` residue
  scan with no Leto storage-slice matches. Residual cfd-math provider work is
  now nalgebra/nalgebra-sparse and other Atlas boundaries.
- **cfd-math**: Removed sparse/basic preconditioner Leto storage-slice
  dependencies. `src/sparse/operations.rs` now stages SpMV input/output and
  row/column scaling buffers through direct `Array1` indexing before
  delegating to Leto CSR operations, and `src/linear_solver/preconditioners/
  basic.rs` now indexes Jacobi diagonals directly. Evidence: cfd-math fmt, lib
  check, focused sparse/preconditioner nextest (95/95), lib/tests clippy,
  all-target check, all-target clippy, and targeted residue scans. Residual
  storage-slice owners are nonlinear mutable dense-workspace helpers and
  multigrid interpolation/smoother internals.
- **cfd-math**: Removed the GPU linear-operator Leto storage-slice dependency.
  `src/linear_solver/operators/gpu.rs` now validates Leto `Array1`
  input/output lengths, stages upload/readback buffers through direct indexing,
  and writes results back by index instead of importing `leto::Storage`,
  borrowing `.storage().as_slice()`, or requiring output `as_slice_mut()`
  contiguity. The execution provider remains the existing Hephaestus-backed
  `cfd-core` GPU context/buffer/kernel path. Evidence: cfd-math fmt,
  GPU-feature check, focused GPU-feature linear_solver::operators nextest
  (5/5), GPU-feature lib clippy, GPU-feature all-target check, GPU-feature
  all-target clippy, and a targeted residue scan. Residual storage-slice
  owners are sparse operations and multigrid internals; broader GPU provider
  work remains outside this operator boundary.
- **cfd-math**: Removed finite-difference operator Leto storage-slice
  dependencies. `src/linear_solver/operators/{poisson,momentum}.rs` now read
  and write `Array1` values through direct indexing for 2D Laplacian, 3D
  Poisson, 1D/2D momentum, and 2D energy operators, removing `leto::Storage`,
  `.storage().as_slice()`, and output `as_slice_mut()` assumptions. Evidence:
  cfd-math fmt, lib check, focused linear_solver::operators nextest (5/5),
  lib/tests clippy, all-target check, all-target clippy, and a targeted
  residue scan. Residual storage-slice owners are sparse operations, GPU
  operator, and multigrid internals.
- **cfd-math**: Removed nonlinear linalg immutable Leto storage-slice
  dependency. `src/nonlinear_solver/linalg.rs` now evaluates vector dot,
  add/sub, scaled add, in-place scaled add, and scale through direct
  `Array1` indexing, and `anderson.rs` now indexes `gamma` directly so the
  immutable `vector_slice` helper and `leto::Storage` import are gone.
  Evidence: cfd-math fmt, lib check, nonlinear_solver nextest (9/9),
  lib/tests clippy, all-target check, all-target clippy, and a targeted
  residue scan. Residual mutable dense-workspace helpers still use
  `StorageMut`; remaining immutable storage-slice owners are sparse
  operations, linear operators, GPU operator, and multigrid internals.
- **cfd-math**: Removed the production SIMD vector Leto storage-slice
  dependency. `src/simd/vector.rs` now implements `simd_mul`, `simd_dot`,
  `simd_norm`, and `par_map` through direct `Array1` indexing while preserving
  Moirai `Adaptive` map/reduce dispatch. Evidence: cfd-math fmt, lib check,
  focused simd::vector nextest (1/1), lib/tests clippy, all-target check,
  all-target clippy, SIMD-filter nextest (26/26), and a targeted residue scan.
  Residual source-level storage-slice owners are nonlinear linalg, sparse
  operations, linear operators, GPU operator, and multigrid internals.
- **cfd-math**: Removed the SIMD integration test Leto storage-slice bridge.
  `tests/simd_tests.rs` now indexes the Leto `spmv` result directly before
  passing values to the existing SIMD slice API, with no `leto::Storage` import
  or `.storage().as_slice()` conversion. Evidence: cfd-math fmt, simd_tests
  check, simd_tests nextest (12/12), simd_tests clippy, cfd-math all-target
  check, cfd-math all-target clippy, SIMD-filter nextest (26/26), and targeted
  residue scans. Residual provider work now sits outside the integration-test
  vector bridge layer: `nalgebra_sparse::CsrMatrix`, dense nalgebra test
  oracles, and source-level Leto storage-slice internals.
- **cfd-math**: Migrated AMG integration vector paths to Leto arrays. Exact
  solutions, RHS values, solver outputs, AMG cycle outputs, and two-grid
  preconditioner work buffers now use `leto::Array1` directly instead of
  converting through nalgebra `DVector` and Leto storage slices. Evidence:
  cfd-math fmt, amg_integration_test check, amg_integration_test nextest
  (5/5), amg_integration_test clippy, cfd-math all-target check, cfd-math
  all-target clippy, AMG-filter nextest (6/6), and a targeted residue scan.
  Residual: this test still uses `nalgebra_sparse::CsrMatrix` and nalgebra
  `DMatrix`/`SymmetricEigen` for the dense energy-norm oracle; the integration
  storage-slice residue was closed by Sprint 1.96.153.
- **cfd-math**: Migrated preconditioner edge-case integration tests to Leto
  arrays. ILU(0), ILU(k), repeated-application, extreme-value, and
  sparsity-preservation tests now build `leto::Array1` RHS/solution buffers
  directly instead of converting through a local nalgebra `DVector`
  preconditioner bridge. Evidence: cfd-math fmt, preconditioner_edge_cases
  check, preconditioner_edge_cases nextest (6/6),
  preconditioner_edge_cases clippy, cfd-math all-target check, cfd-math
  all-target clippy, broader preconditioner nextest (76/76), and a targeted
  residue scan. Residual: this test still uses `nalgebra_sparse::CsrMatrix`
  for matrix storage; later integration-test provider holdouts were closed by
  Sprints 1.96.152 and 1.96.153.
- **cfd-math**: Migrated the linear-solver source test module tree to Leto
  arrays. `src/linear_solver/tests` now allocates `leto::Array1`
  RHS/solution/work buffers directly instead of using local nalgebra `DVector`
  solve/preconditioner bridge macros; residual checks use the Leto SpMV/helper
  path, and the touched solver/sparse cone routes constants/tolerances through
  Eunomia instead of old scalar helper calls. Evidence: cfd-math fmt, lib
  check, all-target check, linear_solver::tests nextest (53/53),
  linear_solver nextest (176/176), lib/tests clippy, all-target clippy, and
  current targeted residue scans. Residual: integration-test vector bridge
  holdouts were closed by Sprints 1.96.152 and 1.96.153; production
  sparse/scalar boundaries still use transitional nalgebra providers.
- **cfd-math**: Migrated core solver validation tests to Leto arrays.
  BiCGSTAB, GMRES, preconditioner integration, and condition-number robustness
  tests now build `leto::Array1` RHS/solution buffers directly instead of
  converting through nalgebra `DVector`, and residual checks now use Leto SpMV
  with explicit thresholds. Evidence: cfd-math fmt, core_solver_tests check,
  core_solver_tests nextest (4/4), core_solver_tests clippy, cfd-math
  all-target check, cfd-math all-target clippy, and a targeted core solver test
  residue scan. Residual: this test still uses `nalgebra_sparse::CsrMatrix`/
  `CooMatrix` for matrix storage, and broader cfd-math test diagnostics still
  contain nalgebra `DVector` bridges.
- **cfd-math**: Migrated simple GMRES integration tests to Leto arrays. Basic,
  restarted, and preconditioned GMRES tests now build `leto::Array1`
  RHS/solution buffers directly instead of converting through nalgebra
  `DVector`, and residual checks now use Leto SpMV with explicit thresholds.
  Evidence: cfd-math fmt, simple_gmres test check, simple_gmres nextest (3/3),
  simple_gmres clippy, cfd-math all-target check, cfd-math all-target clippy,
  touched-file `git diff --check`, and a targeted simple GMRES test residue
  scan. Residual: this test still uses `nalgebra_sparse::CsrMatrix`/
  `CooMatrix` for matrix storage, and broader cfd-math test diagnostics still
  contain nalgebra `DVector` bridges.
- **cfd-math**: Migrated matrix-free solver tests to Leto arrays. Matrix-free
  CG/GMRES identity and scaled-operator tests now build `leto::Array1`
  RHS/solution buffers directly instead of converting through nalgebra
  `DVector`, and the operator-size mismatch test asserts the exact typed
  Leto-boundary error. Evidence: cfd-math fmt, all-target check, all-target
  clippy, focused matrix-free nextest (4/4), touched-file `git diff --check`,
  and a targeted matrix-free test residue scan. Residual: broader cfd-math
  linear-solver integration/adversarial/core/preconditioner test diagnostics
  still contain nalgebra `DVector` bridges, and production sparse/scalar
  boundaries still include `nalgebra_sparse::CsrMatrix` and
  `nalgebra::RealField`.
- **cfd-math**: Migrated BiCGSTAB workspaces and final solver-chain fallback to
  Leto arrays. `BiCGSTAB` now runs preconditioned and unpreconditioned solves
  directly on `leto::Array1`, owns its residual/search/stabilization/operator
  buffers as Leto arrays, and no longer calls legacy nalgebra-vector bridge
  helpers. CG and BiCGSTAB share one Leto vector helper module, the
  `LinearSolverChain` final BiCGSTAB tier now stays on Leto RHS/solution
  buffers, and obsolete bridge helpers were removed from
  `linear_solver/traits.rs`. Evidence: cfd-math fmt, all-target check,
  all-target clippy, focused BiCGSTAB nextest (24/24), broader linear-solver
  nextest (176/176), AMG integration nextest (5/5), touched-file
  `git diff --check`, and targeted provider-residue scans. Residual:
  `nalgebra_sparse::CsrMatrix`, transitional `nalgebra::RealField` scalar
  bounds, and nalgebra `DVector` in remaining matrix-free/preconditioner/
  integration test diagnostics.
- **cfd-math**: Migrated Conjugate Gradient workspaces to Leto arrays.
  `ConjugateGradient` now stores reusable CG vectors as `leto::Array1`, runs
  preconditioned and unpreconditioned solves directly on Leto arrays, and no
  longer calls the linear-solver legacy nalgebra-vector bridge helpers. CG
  benchmark call sites now construct Leto vectors at the measured API, and
  tests assert solved-system values plus exact dimension and max-iteration
  errors. Evidence: cfd-math fmt, all-target check, all-target clippy,
  focused conjugate nextest (13/13), broader linear-solver nextest (176/176),
  and a targeted provider-residue scan. Residual after the BiCGSTAB
  follow-up: the shared linear-solver trait family still carries the
  transitional `nalgebra::RealField` scalar bound and sparse storage remains
  on `nalgebra_sparse::CsrMatrix`.
- **cfd-math**: Migrated Schwarz preconditioner local apply paths to Leto
  arrays. Additive and multiplicative Schwarz now consume/return
  `leto::Array1`, local RHS buffers are Leto-owned, local ILU solves no longer
  bridge through nalgebra `DVector`, and `Preconditioner::apply_to` keeps the
  default additive path on Leto vectors. Residual/output length mismatches
  return typed configuration errors. Evidence: cfd-math fmt, all-target check,
  all-target clippy, focused Schwarz nextest (3/3), broader preconditioner
  nextest (76/76), and a targeted provider-residue scan. Residual: Schwarz
  still stores/constructs local matrices through the shared
  `nalgebra_sparse::CsrMatrix` boundary and participates in the global
  `Preconditioner<T>` nalgebra scalar bound.
- **cfd-math**: Migrated ILU triangular solve workspaces to Leto arrays.
  `IncompleteLU::apply_to` now runs forward/backward substitution over
  `leto::Array1` residual, intermediate, and solution buffers without
  constructing nalgebra `DVector` workspaces; residual/output length
  mismatches return typed configuration errors; and U-solve diagonal identity
  dispatches through Eunomia `NumericElement`. Evidence: cfd-math fmt,
  all-target check, all-target clippy, focused ILU nextest (21/21), broader
  preconditioner nextest (74/74), and a targeted provider-residue scan.
  Residual: IncompleteLU still stores the shared `nalgebra_sparse::CsrMatrix`
  LU factor boundary and participates in the global `Preconditioner<T>`
  nalgebra scalar bound.
- **cfd-math**: Migrated deflation preconditioner vector state to Leto arrays.
  `DeflationPreconditioner` now stores eigenvectors as `leto::Array1`, accepts
  added eigenpairs at the Leto boundary, computes projection coefficients
  without nalgebra `DVector` workspaces, validates residual/output/eigenvector
  lengths with typed errors, and rejects zero eigenvalues before division.
  Evidence: cfd-math fmt, all-target check, all-target clippy, focused
  deflation nextest (3/3), broader preconditioner nextest (73/73), and a
  targeted provider-residue scan. Residual: Deflation still wraps the base
  preconditioner behind the existing `Box<dyn Preconditioner<T>>`, and the
  global preconditioner trait still carries the nalgebra `RealField` bound.
- **cfd-math**: Migrated basic preconditioner vector internals to Leto arrays.
  Identity, Jacobi, and SOR now apply at the `leto::Array1` residual/output
  boundary; Jacobi stores inverse diagonal state as `leto::Array1`; Jacobi/SOR
  validate matrix-sized residuals and all outputs with typed configuration
  errors; and basic scalar identities/tolerances route through Eunomia.
  Evidence: cfd-math fmt, all-target check, all-target clippy, focused
  mismatch-regression nextest (1/1), broader preconditioner nextest (70/70),
  and a targeted provider-residue scan. Residual: basic matrix-backed
  preconditioners still store the shared `nalgebra_sparse::CsrMatrix` boundary
  and participate in the global `Preconditioner<T>` nalgebra scalar bound.
- **cfd-math**: Migrated IncompleteCholesky preconditioner solve workspaces to
  Leto arrays. Forward/backward triangular substitution and
  `Preconditioner::apply_to` now operate on `leto::Array1` buffers without
  constructing nalgebra `DVector` residual, intermediate, or solution
  workspaces; residual/output length mismatches now surface typed
  configuration errors. The IC(0) factorization square root dispatch now uses
  Eunomia `NumericElement`. Evidence: cfd-math fmt, all-target check,
  all-target clippy, focused Cholesky nextest (5/5), and a targeted
  provider-residue scan. Follow-up: IncompleteCholesky factor storage moved
  to Leto CSR in Sprint 1.96.166; residual sparse-provider work is now
  Schwarz, direct solver, and the shared solver matrix boundary.
- **cfd-math**: Migrated SSOR preconditioner sweep workspaces to Leto arrays.
  SSOR forward/backward sweeps and `Preconditioner::apply_to` now operate on
  `leto::Array1` buffers without constructing nalgebra `DVector`
  residual/solution workspaces, and residual/output length mismatches now
  surface typed configuration errors. Evidence: cfd-math fmt, all-target
  check, all-target clippy, focused SSOR nextest (5/5), and a targeted
  provider-residue scan. Residual: SSOR still stores the shared
  `nalgebra_sparse::CsrMatrix` matrix boundary and participates in the global
  `Preconditioner<T>` nalgebra scalar bound.
- **cfd-math**: Migrated block/SIMPLE preconditioner vector internals to Leto
  arrays. `DiagonalPreconditioner`, `BlockDiagonalPreconditioner`, and
  `SimplePreconditioner` now store and apply diagonal/Schur state through
  `leto::Array1` and no longer build nalgebra `DVector` bridges in
  `Preconditioner::apply_to`; vector length mismatches now surface typed
  configuration errors instead of returning cloned inputs. Evidence:
  cfd-math fmt, all-target check, all-target clippy, focused
  block-preconditioner nextest (4/4), and a targeted provider-residue scan.
  Residual: the global `Preconditioner` trait still carries the transitional
  `nalgebra::RealField` scalar bound, and sparse storage still aliases
  `nalgebra_sparse::CsrMatrix`.
- **cfd-math**: Migrated GMRES workspace and solver-chain GMRES tiers to Leto
  arrays. GMRES Arnoldi basis extraction, work/preconditioner vectors, residual
  checks, and direct GMRES solve methods now use `leto::Array1<T>` instead of
  nalgebra `DVector`; `LinearSolverChain` now keeps its GMRES+AMG,
  GMRES+block, unpreconditioned GMRES, and GMRES+ILU tiers on Leto arrays.
  Evidence: cfd-math lib/tests/all-targets check, cfd-math fmt, cfd-math
  all-target clippy, focused GMRES nextest (21/21), cfd-2d no-default lib
  check, focused cfd-2d momentum nextest (53/53), and targeted GMRES
  DVector/legacy-bridge residue scan. Residual after the CG/BiCGSTAB
  follow-up: cfd-math still retains nalgebra scalar bounds and
  nalgebra-sparse storage.
- **cfd-2d**: Migrated momentum vectors and scalar seam off direct nalgebra.
  `physics::momentum` now stores RHS and solution vectors as `leto::Array1<T>`,
  momentum boundary helpers mutate Leto RHS arrays, momentum solves call
  Leto-native iterative/direct solver APIs, and the obsolete local nalgebra
  conversion bridge module was removed. `Cfd2dScalar` no longer requires
  `nalgebra::RealField`, and `cfd-2d` no longer declares direct
  `nalgebra`/`nalgebra-sparse` dependencies. Evidence: cfd-2d fmt,
  no-default lib check, no-default all-target clippy, focused `momentum`
  nextest (53/53), direct cfd-2d source/manifest nalgebra residue scan, and
  cargo-tree proof that remaining nalgebra is transitive through upstream
  owners. Residual: nalgebra/nalgebra-sparse still resolve transitively through
  cfd-1d, cfd-core, cfd-math, cfd-schematics, and Gaia.
- **cfd-2d**: Migrated pressure-velocity correction vectors to Leto arrays.
  `pressure_velocity` now stores cached RHS and pressure-correction solution
  vectors as `leto::Array1<T>` instead of nalgebra `DVector`, dispatches
  directly through the Leto-native iterative solver trait and direct sparse
  fallback, and computes RHS diagnostics with `leto_ops::norm_l2`. Evidence:
  cfd-2d fmt, no-default lib check, no-default all-target clippy, focused
  `pressure_velocity` nextest (16/16), and targeted `pressure_velocity`
  DVector/nalgebra/bridge residue scan. Residual: direct cfd-2d
  source/manifest nalgebra ownership is now removed; nalgebra remains
  transitive through upstream owners.
- **cfd-2d**: Migrated SIMPLE pressure-correction vectors to Leto arrays.
  `solvers::simple` now stores pressure-correction RHS and `p_prime` as
  `leto::Array1<T>` instead of nalgebra `DVector`, and the pressure solve calls
  the Leto-native `IterativeLinearSolver::solve` boundary directly. Evidence:
  cfd-2d fmt, no-default lib check, no-default all-target clippy, focused
  `simple` nextest (19/19), and targeted `solvers/simple`
  DVector/nalgebra/bridge residue scan. Residual: direct cfd-2d
  source/manifest nalgebra ownership is now removed; nalgebra remains
  transitive through upstream owners.
- **cfd-2d**: Migrated FDM RHS and Gauss-Seidel vectors to Leto arrays.
  Poisson, advection-diffusion, and shared Gauss-Seidel FDM paths now use
  `leto::Array1<T>` for RHS and result vectors instead of nalgebra `DVector`.
  Evidence: cfd-2d fmt, no-default lib check, no-default all-target clippy,
  focused `fdm` nextest (2/2), and targeted `solvers/fdm` DVector/nalgebra
  residue scan. Residual: direct cfd-2d source/manifest nalgebra ownership is
  now removed; nalgebra remains transitive through upstream owners.
- **cfd-2d**: Migrated time-integration state vectors to Leto arrays.
  Explicit, implicit, multistep, adaptive-controller, adaptive-integrator, and
  time test paths now use `leto::Array1` through the `StateVector<T>` domain
  alias instead of nalgebra `DVector`; fixed-point convergence uses
  `leto_ops::norm_l2`. Evidence: cfd-2d fmt, no-default lib check,
  no-default all-target clippy, focused `time` nextest (29/29), and targeted
  `schemes/time` DVector/nalgebra residue scans. Residual: cfd-2d still
  has no direct cfd-2d source/manifest nalgebra ownership; nalgebra remains
  transitive through upstream owners.
- **cfd-2d**: Migrated compact `DMatrix` residue to Leto arrays.
  Immersed-boundary force/velocity matrices, `schemes::Grid2D` storage,
  dependent scheme callers/tests, and the `blood_venturi` example now use
  `leto::Array2` shape/indexing instead of nalgebra `DMatrix`. Evidence:
  cfd-2d fmt, no-default lib check/clippy, no-default `blood_venturi` example
  check, focused immersed-boundary/schemes/upwind/MUSCL nextest (60/60), and
  targeted DMatrix plus nalgebra-style grid-access residue scans. Residual:
  direct cfd-2d source/manifest nalgebra ownership is now removed; nalgebra
  remains transitive through upstream owners.
- **cfd-1d**: Migrated vascular Bessel/Womersley complex math to Eunomia.
  The Bessel recurrence and Womersley profile now use `eunomia::Complex`
  instead of nalgebra `Complex`/`ComplexField`, and convergence checks use
  Eunomia's complex `norm()`. Evidence: cfd-1d lib check/clippy, focused
  Bessel/Womersley nextest (26/26), fmt check, and targeted nalgebra-complex
  residue scan. Residual: cfd-1d still retains nalgebra for network
  linear-system storage and the transitional scalar seam.
- **cfd-math**: Migrated the public `LinearOperator::apply` vector boundary
  to Leto arrays. `apply` and `apply_transpose` now accept `leto::Array1<T>`;
  sparse CSR, identity/scaled, Poisson, momentum, energy, and GPU operator
  adapters implement that boundary; and CG/BiCGSTAB bridge at their current
  nalgebra workspaces. Evidence: cfd-math all-target check/clippy,
  focused cfd-math solver/operator nextest (80/80), fmt check, and targeted
  DVector operator-signature residue scan. Residual: solver workspaces,
  nalgebra sparse storage, and some preconditioner internals still retain local
  `DVector`/`CsrMatrix` conversion bridges.
- **cfd-math/cfd-1d**: Migrated the public `Preconditioner::apply_to`
  residual/result boundary to Leto arrays. Concrete cfd-math preconditioners,
  AMG, Schwarz local solves, deflation, and preconditioner tests now call the
  trait with `leto::Array1<T>` buffers, and cfd-1d network Jacobi
  preconditioning implements the same provider boundary. Evidence: cfd-math
  all-target check/clippy, cfd-1d lib check/clippy, focused cfd-math
  solver/preconditioner nextest (131/131), fmt check, and targeted
  `apply_to` DVector-signature residue scan. Residual: `LinearOperator::apply`,
  nalgebra sparse preconditioner internals, and iterative solver workspaces
  still carry nalgebra conversion bridges.
- **cfd-math/cfd-1d/cfd-2d/cfd-3d**: Migrated the public
  `IterativeLinearSolver::solve` RHS/result boundary to Leto arrays. CG,
  BiCGSTAB, and GMRES now expose Leto buffers at the trait boundary while
  confining nalgebra work vectors internally; cfd-math tests use Leto boundary
  arrays; cfd-1d network solves, cfd-2d momentum/pressure solves, and cfd-3d
  FEM projection route through local Leto bridges. Evidence: cfd-math
  fmt/check/clippy, focused cfd-math solver nextest (61/61), cfd-1d/cfd-2d/
  cfd-3d checks and clippy, cfd-validation all-target clippy, targeted
  DVector call-site residue scan, and `git diff --check`. Residual: operator
  and preconditioner traits still expose nalgebra vector/storage boundaries.
- **cfd-math/cfd-validation**: Migrated the public iterative linear-solver
  `solve_system` vector boundary to Leto arrays. `LinearSolver::solve_system`
  and the CG, BiCGSTAB, and GMRES implementations now accept Leto RHS and
  optional initial-guess arrays and return Leto result arrays; nalgebra work
  vectors are confined to solver-local workspace bridges; and cfd-validation
  numerical solver validation calls the new Leto API while retaining DVector
  only for existing error metrics. Evidence: cfd-math/cfd-validation
  fmt/check, cfd-math all-target clippy, cfd-validation lib/all-target
  clippy, cfd-2d all-target clippy, focused cfd-math CG/BiCGSTAB/GMRES
  nextest (58/58), targeted `solve_system` signature/residue scans, and
  `git diff --check`. Residual: the broader iterative operator and
  preconditioner trait family still exposes nalgebra vector/storage
  boundaries.
- **cfd-validation**: Migrated SpMV benchmark callers and validation scalar
  bounds to Leto/Eunomia provider seams. Matrix-operation profiling,
  algorithm profiling, and generic SpMV benchmarking now use `leto::Array1`
  for the public `cfd_math::sparse::spmv` API; linear-solver validation and
  1D blood-flow literature validation now use the crate-local
  `ValidationScalar` contract needed by the migrated cfd-math/cfd-1d provider
  surfaces. Evidence: cfd-validation fmt/check, cfd-validation lib and
  all-target clippy, cfd-2d no-default all-target clippy, focused benchmark
  nextest (40/40), targeted SpMV DVector-residue scan, and `git diff --check`.
  Residual: broad cfd-validation nextest still fails in two venturi
  cross-fidelity convergence tests unrelated to this provider-boundary slice.
- **cfd-math/cfd-2d/cfd-3d**: Migrated the solver-chain and FEM/direct
  fallback vector boundary to Leto arrays. `LinearSolverChain::solve` and
  `solve_with_guess` now accept and return `leto::Array1<T>`; cfd-2d
  momentum/pressure direct fallbacks share a private Leto bridge into
  `DirectSparseSolver`; cfd-3d FEM assembly and chain dispatch share a private
  FEM Leto bridge; and the 2D/3D scalar seams now carry the Leto real-scalar
  provider contract. Evidence: focused cfd-math/cfd-1d/cfd-2d/cfd-3d
  no-default checks, focused cfd-math solver-chain nextest (4/4), cfd-math
  all-target clippy, and cfd-2d/cfd-3d no-default lib clippy. Residual:
  cfd-validation still blocks cfd-2d all-target clippy until its SpMV
  benchmark vectors and generic scalar bounds migrate to the same providers.
- **cfd-math**: Migrated the direct sparse-solver vector API to Leto.
  `DirectSparseSolver::solve` now accepts and returns `leto::Array1<T>`;
  rsparse RHS conversion, sparse solution construction, finite checks, and
  dense fallback all operate on Leto arrays. The obsolete DVector dense-fallback
  wrapper was removed, and `LinearSolverChain` now converts only at its
  still-nalgebra-facing tier boundary. Evidence: cfd-math fmt/check, focused
  direct-solver/chain/core-solver/simple-GMRES nextest (4/4), no-default
  all-target clippy, cfd-math docs/doctests, targeted direct-solver
  DVector-signature scan, and `git diff --check`. Residual: the broader
  linear-solver chain, iterative solver traits, preconditioners, and sparse
  storage still expose nalgebra vector/storage boundaries.
- **cfd-math**: Migrated sparse builder Dirichlet RHS assembly to Leto.
  `SparseMatrixBuilder::build_with_rhs` now accepts `leto::Array1<T>` and
  mutates that provider-owned RHS for column elimination; `build()` no longer
  carries a dummy nalgebra vector; and direct/block-preconditioner tests use
  Leto RHS values for assembly while retaining nalgebra vectors only at the
  unchanged solver boundary. Evidence: cfd-math fmt/check, focused
  sparse/direct-solver/block-preconditioner nextest (25/25), no-default
  all-target clippy, targeted `build_with_rhs` residue scan, and
  `git diff --check`. Residual: the public sparse storage alias and
  linear-solver/preconditioner traits still expose `nalgebra_sparse::CsrMatrix`
  and nalgebra `DVector`.
- **cfd-math**: Migrated public sparse SpMV wrappers to Leto vectors. The
  then-public entry points accepted `leto::Array1` input/output vectors while
  delegating to `leto_ops::spmv_into`; the redundant parallel-named wrapper
  was later removed under Breaking. The remaining
  nalgebra `DVector` SpMV bridge is private to `LinearOperator for CsrMatrix`
  until the linear-solver trait family moves to Leto. Sparse tests,
  GMRES/AMG integration tests, interpolation quality checks, and the SpMV
  benchmark now use Leto arrays at public SpMV call sites. Evidence:
  cfd-math fmt/check, focused sparse/spmv/interpolation/AMG/solver nextest
  (40/40), no-default all-target clippy, public SpMV signature scan, and
  `git diff --check`. Residual: the public linear-solver/preconditioner trait
  family still exposes nalgebra `DVector`, and sparse storage still aliases
  `nalgebra_sparse::CsrMatrix`.
- **cfd-math**: Migrated `SparseMatrixExt` diagonal and scaling vector
  surfaces to Leto arrays. `diagonal` now returns `leto::Array1`, and
  `set_diagonal`, `scale_rows`, and `scale_columns` now accept Leto arrays
  instead of nalgebra `DVector`. `JacobiPreconditioner::new` consumes that
  provider-backed diagonal through `Storage::as_slice`, while retaining
  nalgebra storage only at the current public `Preconditioner` trait boundary.
  Evidence: cfd-math fmt/check, focused sparse plus basic-preconditioner
  nextest (21/21), no-default all-target clippy, targeted
  `SparseMatrixExt` DVector-signature residue scan, and `git diff --check`.
  Residual: linear-solver/preconditioner traits and
  `nalgebra_sparse::CsrMatrix` storage still await the larger boundary
  migration.
- **cfd-math**: Migrated multigrid smoother and cycle vector paths to Leto.
  `MultigridLevel`, `AMGHierarchy`, and `MultigridSmoother` now use
  `leto::Array1` through `MultigridVector`; Jacobi, Gauss-Seidel, symmetric
  Gauss-Seidel, SSOR, and Chebyshev smoothers compute residuals through the
  new Leto-array `sparse::spmv_array` bridge; V/W/F cycles and coarsest solves
  consume and return Leto vectors; and AMG V-cycle recursion bridges to
  nalgebra `DVector` only at the existing public preconditioner boundary.
  Evidence: cfd-math fmt/check, focused multigrid cycle+smoother nextest
  (10/10), no-default all-target clippy, clean smoother/cycle provider-residue
  scan, and `git diff --check`. Residual: public sparse and linear-solver APIs
  still expose `nalgebra_sparse::CsrMatrix`, nalgebra `DVector`, and
  nalgebra `RealField`.
- **cfd-math**: Migrated geometric multigrid (GMG) to Leto/Eunomia.
  `GeometricMultigrid`, `NonlinearOperator`, FAS/linear solves, Poisson
  hierarchy matrices, transfer operators, Jacobi relaxation, residual
  computation, and tests now use `leto::Array2`/`Array1` plus Eunomia scalar
  traits instead of nalgebra `DMatrix`/`DVector` and `num_traits`
  conversions. Evidence: cfd-math fmt/check, focused GMG nextest (5/5),
  no-default all-target clippy, clean GMG provider-residue scan, and
  `git diff --check`. Residual: the public linear-solver/preconditioner
  boundary still exposes nalgebra `DVector`/`RealField`.
- **cfd-math**: Migrated GMRES internal dense workspace to Leto/Eunomia.
  `gmres::arnoldi` now stores the Krylov basis and Hessenberg matrix in
  `leto::Array2`; `gmres::givens` stores rotation coefficients,
  least-squares RHS, and triangular-solve results in `leto::Array1`/`Array2`
  and uses Eunomia scalar-field operations; and `GMRESWorkspace` no longer
  allocates nalgebra dense matrices for restarted Krylov state.
  `LinearSolverChain` carries the Eunomia real-field bound needed by the
  migrated GMRES constructor. Evidence: cfd-math fmt/check, focused GMRES
  nextest (21/21), no-default all-target clippy, clean GMRES provider-residue
  scan for replaced providers, and `git diff --check`. Residual: the public
  linear-solver trait boundary still exposes nalgebra `DVector` and
  `nalgebra::RealField`.
- **cfd-math**: Migrated the standalone AMG restriction dense-transfer
  utilities to Leto arrays. `create_*_restriction`,
  `validate_restriction_operator`, `restrict_vector`, and `restrict_matrix`
  now expose and consume `leto::Array2`/`Array1`; Galerkin projection delegates
  dense `R * A * P` products to `leto_ops::MatrixProduct` instead of nalgebra
  multiplication. Tests now assert exact `P^T v` and `P^T A P` fixture values.
  Evidence: cfd-math fmt/check, focused restriction nextest (7/7),
  no-default all-target clippy, and a clean restriction provider-residue scan.
  Residual: broader sparse/linear-solver APIs still expose nalgebra
  dense/sparse boundaries.
- **cfd-math**: Consolidated legacy CSR-to-dense linear-solver paths behind a
  Leto dense bridge. `linear_solver::dense_bridge` now converts the current
  `nalgebra_sparse::CsrMatrix`/nalgebra `DVector` boundary into
  `leto::Array2`/`Array1` and solves with `leto_ops::solve`.
  `DirectSparseSolver` dense fallback and multigrid cycle coarsest
  small-system solves both consume that bridge; the multigrid local
  Gaussian-elimination helper was removed. `LinearSolverChain` carries the
  Leto real-scalar provider bound required by the direct-solver path.
  Evidence: cfd-math fmt/check, focused direct-solver and multigrid cycle
  nextest (9/9), and no-default all-target clippy pass. Residual: the public
  sparse/linear-solver storage contracts still expose
  `nalgebra_sparse::CsrMatrix` and nalgebra `DVector`, and the primary
  sparse-LU path still uses `rsparse`.
- **cfd-math**: Centralized the sparse Leto/nalgebra CSR bridge and routed
  sparse builder construction through Leto provider validation.
  `SparseMatrixBuilder::{build,build_with_rhs,build_parallel}` and
  `ParallelAssembly::block_diagonal` now construct `leto_ops::CsrMatrix` first
  and convert through `sparse::bridge` at the remaining legacy
  `nalgebra_sparse::CsrMatrix` boundary. Evidence: cfd-math check, focused
  sparse nextest (18/18), and no-default all-target clippy pass. Residual:
  the public sparse/linear-solver storage contracts still expose
  `nalgebra_sparse::CsrMatrix` and nalgebra `DVector`.
- **cfd-math**: Routed the remaining sparse extension value operations through
  the Leto CSR provider. `SparseMatrixExt` now bridges diagonal extraction,
  scalar/value scaling, row scaling, column scaling, Frobenius norm, diagonal
  dominance, and the condition-estimate heuristic into `leto_ops::CsrMatrix`;
  `JacobiPreconditioner::new` carries the Leto scalar bound needed by the
  provider-backed diagonal path. Evidence: cfd-math fmt/check, focused sparse
  nextest (18/18), and no-default all-target clippy pass. Residual:
  sparse/linear-solver storage still exposes `nalgebra_sparse::CsrMatrix` and
  nalgebra `DVector`.
- **cfd-math**: Routed sparse CSR×CSR products and AMG Galerkin operators
  through the Leto provider. `try_sparse_sparse_mul` now bridges the current
  `nalgebra_sparse::CsrMatrix` boundary into `leto_ops::spgemm`, surfaces
  invalid CSR or shape errors as typed CFDrs errors, and AMG hierarchy
  recompute/setup call the fallible Leto-backed path. `try_sparse_transpose`
  now bridges CSR transpose into `leto_ops::CsrMatrix::transpose`, and AMG
  restriction setup no longer calls `nalgebra_sparse::transpose_as_csc`.
  `try_spmv` now bridges the current sparse matrix/vector boundary into
  `leto_ops::spmv_into`, so CFDrs no longer owns separate scalar and parallel
  CSR SpMV loops. Evidence: cfd-math fmt/check, focused sparse nextest (17/17),
  focused interpolation nextest (15/15), focused AMG nextest (6/6), and
  no-default all-target clippy pass. Residual:
  sparse/linear-solver storage still exposes `nalgebra_sparse::CsrMatrix` and
  nalgebra `DVector`.
- **leto-ops provider**: Closed the CSR×CSR product gap needed by the
  `cfd-math` sparse/AMG migration. `leto_ops::spgemm` now provides a
  Leto-owned sparse matrix-matrix product for Galerkin operators, and
  `CsrRow::nnz` covers row-cardinality consumers. Evidence: upstream
  Leto-ops fmt/check/focused sparse nextest/clippy/doc pass.
- **cfd-math**: Migrated the SIMD helper cone to Leto/Eunomia. `SimdVectorOps`
  now implements Leto `Array1` operations instead of nalgebra `DVector`, SIMD
  sparse matvec accepts `leto_ops::CsrMatrix` and delegates to `leto_ops::spmv`,
  and generic SIMD/field/vectorization bounds use Eunomia scalar traits. The
  SIMD integration residual test now uses Leto-ops CSR parts instead of
  `nalgebra_sparse`. Evidence: cfd-math fmt/check, focused SIMD nextest
  (26/26, 318 skipped), no-default all-target clippy, and a clean
  `src/simd` plus `tests/simd_tests.rs` provider-residue scan.
- **cfd-math**: Migrated the nonlinear-solver vector cone to Leto.
  `JfnkConfig`, `JfnkSolver`, finite-difference JvP, restarted GMRES work
  vectors, and nonlinear-solver tests now use `leto::Array1`/`Array2` and
  Eunomia scalar traits instead of nalgebra `DVector`/`DMatrix`/`RealField`.
  Anderson and JFNK now share one local `nonlinear_solver::linalg` helper
  module for Leto vector construction, norms, dot products, scaled updates,
  and small dense helper operations. Evidence: cfd-math fmt/check, focused
  nonlinear nextest (9/9, 335 skipped), no-default all-target clippy, and a
  clean `src/nonlinear_solver` provider-residue scan.
- **cfd-math**: Migrated the high-order DG dense-array cone to Leto.
  `DGBasis`, `DGSolution`, numerical fluxes, `DGOperator`, limiters,
  `DGSolver`, and DG time integrators now use `leto::Array1`/`Array2` instead
  of nalgebra `DVector`/`DMatrix`. Projection, derivative, RHS, and implicit
  Newton correction solves now surface Leto solver errors instead of silently
  falling back. DG examples and DG-related benchmarks now compile against Leto.
  Evidence: cfd-math fmt/check, DG benchmark checks, focused DG nextest
  (62/62, 282 skipped), no-default all-target clippy, no-default doctests
  (3 passed, 3 ignored), and a clean high-order provider-residue scan.
- **cfd-math**: Migrated the high-order spectral dense-array cone to Leto.
  `SpectralElement`, `SpectralMesh1D`, `SpectralDiffOp`,
  `SpectralInterp`, `SpectralQuadrature`, `SpectralFilter`, and spectral
  time-integration helpers now use `leto::Array1`/`Array2` instead of nalgebra
  `DVector`/`DMatrix`; spectral assembly now uses Leto arrays for local dense
  matrices, optional local RHS values, and debug dense CSR materialization.
  Local derivative, stiffness, dot, and matrix-vector paths share Leto helpers.
  `SpectralInterp::l2_projection` now returns a typed result and surfaces Leto
  solve failures instead of silently falling back to interpolation. Evidence:
  cfd-math fmt/check, focused spectral nextest (13/13, 331 skipped),
  no-default all-target clippy, and a clean
  `high_order/spectral` provider-residue scan.
- **cfd-math**: Migrated the high-order WENO scalar-provider cone to Eunomia.
  `WENO5`, `WENO7`, and `WenoReconstruction` now bind on
  `eunomia::RealField`/`FloatElement` instead of nalgebra scalar traits; WENO
  constants route through the existing provider conversion helper and nonlinear
  weights use generic multiplication for squaring. Evidence: cfd-math
  fmt/check, focused WENO nextest (6/6, 338 skipped), no-default all-target
  clippy, and a clean `high_order/weno` provider-residue scan.
- **cfd-2d/cfd-3d/cfd-validation**: Added crate-level scalar seams for the
  next migration boundary. Public `cfd-2d::Cfd2dScalar` plus crate-internal
  `Cfd3dScalar` and `ValidationScalar` now route cfd-core boundary/fluid
  consumers through Eunomia-compatible contracts while the remaining
  nalgebra-backed matrix surfaces are isolated. cfd-2d NS-FVM and moving-wall
  tests, cfd-3d FEM/DES and branch/venturi setup, and cfd-validation
  Casson/Womersley/cross-fidelity benchmarks now use these seams with Leto
  boundary vectors. The LES Smagorinsky GPU update path now refreshes
  Yoshizawa SGS kinetic-energy/dissipation diagnostics after GPU viscosity
  computation. Evidence: cfd-2d fmt/check/library clippy, cfd-2d GPU-feature
  lib nextest (518/518, 1 skipped), cfd-3d no-default library check/clippy, and
  cfd-validation no-default library check/clippy. Follow-up cfd-2d all-targets
  cleanup moved examples/integration tests onto Leto boundary values and the
  `Cfd2dScalar` seam, resolved the crate's example/test clippy blockers, and
  verified `cargo clippy -p cfd-2d --no-default-features --features gpu
  --all-targets -- -D warnings` plus `cargo nextest run -p cfd-2d
  --no-default-features --features gpu --status-level fail` (572/572, 27
  skipped). cfd-validation geometry now re-exports Leto point types and routes
  its 2D/3D geometry implementations plus directly dependent benchmark wrappers
  through Eunomia scalar contracts instead of nalgebra scalar/point surfaces;
  validation tests that target migrated cfd-core boundary vector APIs now
  construct Leto vectors. Evidence: cfd-validation fmt, no-default library
  check, no-default all-target clippy, focused geometry nextest (11/11, 420
  skipped), and a clean migrated geometry/benchmark residue scan. The
  analytical solution trait cone now also returns Leto `Vector3` values and
  binds on Eunomia scalar contracts across `src/analytical/**` and
  `src/solutions/mod.rs`; focused analytical nextest passes 20/20 with 411
  skipped and the analytical/solutions residue scan is clean. The standalone
  error-metrics cone now uses Eunomia scalar contracts and Leto `Vector3`
  magnitudes across `src/error_metrics/**`; focused error-metrics nextest
  passes 21/21 with 410 skipped, all-target clippy is clean, and the
  error-metrics residue scan finds no nalgebra imports or old `T::zero()`/
  `T::one()` identities. The standalone convergence cone now uses Eunomia
  scalar contracts and crate-local scalar helpers across `src/convergence/**`;
  focused convergence nextest passes 27/27 with 404 skipped, all-target clippy
  is clean, and the convergence residue scan finds no nalgebra imports or old
  scalar identities. The edge-case testing cone now binds on Eunomia
  `RealField` across `src/edge_case_testing/**`; focused edge-case nextest
  passes 15/15 with 416 skipped, all-target clippy is clean, and the
  edge-case residue scan finds no nalgebra imports.
  The manufactured-solution cone now binds on Eunomia scalar contracts and
  crate-local scalar helpers across `src/manufactured/**`; focused
  manufactured nextest passes 50/50 with 381 skipped, all-target clippy is
  clean, and the manufactured residue scan finds no nalgebra imports or old
  scalar identities.
  The conservation cone now uses Leto `Array1`/`Array2` flow fields and
  Eunomia scalar contracts across `src/conservation/**`; focused conservation
  nextest passes 18/18 with 413 skipped, all-target clippy is clean, and the
  conservation residue scan finds no nalgebra dense storage or old scalar
  identities.
  The time-integration follow-up moves `cfd-math::time_stepping::stability`
  RK tableaus to Leto `Array2`/`Array1` and Eunomia `RealField`, and moves
  `cfd-validation/src/time_integration/**` Euler/RK2/RK4 state validation to
  Leto `Array1`; focused cfd-math stability nextest passes 5/5 with 328
  skipped, focused cfd-validation time-integration nextest passes 12/12 with
  419 skipped, all-target clippy is clean for both crates, and the migrated
  cone scan finds no nalgebra, ndarray, num-traits, or old scalar-identity
  residue.
  The cfd-math time-stepping follow-up moves `traits`, explicit RK, adaptive
  Dormand-Prince, RKC, exponential, and IMEX APIs to Leto-backed
  `TimeState<T>` and shared `TimeMatrix<T>` where dense matrices are required.
  Scalar contracts route through Eunomia, dense matrix exponential delegates to
  `leto_ops::matexp`, and IMEX Newton corrections now use `leto_ops::solve`
  instead of a nalgebra dense LU bridge. `rk4_bench` and `imex_bench` target
  the migrated API. Evidence: cfd-math fmt/check, focused nextest filters for
  `runge_kutta` (5/5), `adaptive` (3/3), `chebyshev` (5/5), `exponential`
  (6/6), and `imex` (5/5), no-default all-target clippy, and residue scans
  showing no nalgebra, `DMatrix`, `DVector`, `num_traits`, `num_complex`, or
  ndarray hits in the time-stepping cone or RK4/IMEX benches.
  The differentiation follow-up moves finite-difference derivative outputs to
  Leto `Array1`, gradient/divergence/curl vector surfaces to Leto `Vector3`,
  scalar constants to Eunomia, and replaces the type-suffixed
  `first_derivative_simd_f32` API with `FiniteDifference<f32>::
  first_derivative_simd`; focused differentiation nextest passes 12/12 and the
  differentiation scan finds no old provider residue.
- **cfd-1d/Eunomia**: Migrated the 1D scalar boundary to a crate-local
  `Cfd1dScalar` provider contract that carries the remaining nalgebra
  linear-system backend requirement plus the Eunomia scalar contract required
  by migrated cfd-core fluid, boundary, and domain APIs. Domain, network,
  resistance, vascular, solver, transient, and analysis bounds now route
  through that seam, and `NetworkDomain::contains_1d` now accepts Leto
  `Point1<T>`. Evidence: cfd-1d fmt, no-default library check, no-default
  nextest (725/725, 3 skipped), no-default library clippy, and a residue scan
  showing nalgebra `RealField` only in `src/scalar.rs`. Downstream `cfd-2d`
  GPU consumer verification now reaches `cfd-2d` and is blocked by `cfd-2d`'s
  own Eunomia/nalgebra scalar-bound migration errors.
- **cfd-core/Hephaestus**: Migrated GPU Poisson solver dispatch and transfer
  ownership to Hephaestus. `GpuPoissonSolver` now wraps the supplied WGPU
  device/queue in `WgpuDevice`, compiles Jacobi, red-black, and residual WGSL
  entry points as `WgslMultiStorageKernel`s, and uses `ComputeDevice` upload,
  allocation, and download APIs for field/source/work/residual buffers. The
  solver now validates the constructor grid geometry instead of inferring a
  square grid from `phi.len()`. Evidence: cfd-core fmt, GPU-feature check,
  no-default check, GPU-feature all-target clippy, full GPU-feature nextest
  (231/231), and clean Poisson raw-WGPU residue scans. Downstream `cfd-2d`
  GPU consumer verification now reaches `cfd-2d` and remains blocked by
  `cfd-2d` Eunomia/nalgebra trait-bound errors in the current workspace graph.
- **cfd-core/Hephaestus**: Migrated the GPU buffer/scalar provider cone to
  Hephaestus and Eunomia. `GpuBuffer<T>` now stores Hephaestus `WgpuBuffer<T>`
  and uses `WgpuDevice`'s `ComputeDevice` allocation, upload, download, and
  write contracts; raw WGPU access remains only at existing bind-group
  construction points. GPU buffer, pipeline, and kernel scalar bounds now use
  `eunomia::RealField` instead of nalgebra. Evidence: cfd-core fmt,
  GPU-feature check, no-default check, GPU-feature all-target clippy, full
  no-default nextest (201/201), full GPU-feature nextest (231/231), and clean
  GPU compute provider-residue scans. CFDrs-owned WGSL CFD kernel/pipeline
  orchestration remains a follow-up Hephaestus provider increment.
- **cfd-core/Leto**: Migrated `cfd-core::geometry::{mesh,staggered}` to
  Leto/Eunomia geometry storage. `Mesh<T>` now stores Leto `Point3<T>`, mesh
  operations use Leto `Vector3<T>` and `FixedMatrix<T, 3, 3>`, mesh service/
  refinement/quality/statistics/connectivity bounds use `eunomia::RealField`,
  and `StaggeredGrid2D` no longer imports nalgebra. Leto now provides
  `FixedMatrix<T, 3, 3> * Vector3<T>` for mesh rotation. Evidence: cfd-core
  fmt/check/no-default all-target clippy, focused mesh/staggered nextest
  (12/12), full cfd-core no-default nextest (201/201), Leto check/clippy/full
  nextest (171/171), and clean mesh/staggered provider scans. The Gaia
  topology replacement remains a follow-up.
- **cfd-core/Leto**: Migrated the cfd-core geometry/domain, boundary,
  fluid-service, material, management, and solver factory scalar contracts in
  the active cone to `eunomia::RealField` and `leto::geometry`. Domain shapes
  now use Leto `Point1`/`Point2`/`Point3`/`Vector3`; boundary condition,
  boundary geometry, wall, inlet/outlet, applicator, material database,
  fluid-dynamics service, management conversion/broadcast/factory, and solver
  trait/config helpers no longer bind on nalgebra scalar contracts in the
  migrated cone. Leto now provides the missing `Point1<T>` primitive,
  conditional `Eq` derives for fixed geometry values, and serde `std`/`alloc`
  feature propagation. Evidence: cfd-core fmt/check/no-default all-target
  clippy, full cfd-core no-default nextest (201/201), Leto check/clippy/nextest
  (170/170), and clean migrated-cone provider scans.
- **cfd-core**: Migrated the compute execution boundary to Eunomia and problem
  gravity to Leto geometry. `ComputeKernel`, `ComputeBuffer`, `CpuBuffer`,
  `CpuAdvectionKernel`, and `ComputeDispatcher` now use
  `eunomia::RealField`; CPU zero initialization and upwind sign checks use
  Eunomia identities instead of nalgebra `T::zero()`. `ProblemParameters::
  gravity` and `ProblemBuilder::gravity` now use `leto::geometry::Vector3<T>`.
  Evidence: cfd-core fmt/check/no-default all-target clippy, focused
  compute/problem nextest (21/21), and clean touched-file provider scans.
- **cfd-core/Leto**: Migrated `abstractions::state` scalar and vector field
  storage to Leto/Eunomia. `FieldData::Scalar` now stores `leto::Array1<T>`
  instead of nalgebra `DVector<T>`, vector fields use
  `leto::geometry::Vector3<T>`, and state scalar contracts use
  `eunomia::RealField`. Leto now supplies validated serde support for owned
  arrays/storage, preserving CFDrs serialized field-state behavior without a
  downstream wrapper. Evidence: focused Leto serde nextest, Leto clippy,
  cfd-core fmt/check/no-default all-target clippy, focused cfd-core state
  nextest (3/3), and a clean provider scan over `abstractions/state.rs`.
- **cfd-core**: Migrated `compute::time` integrator/controller state from
  nalgebra/num-traits to Leto/Eunomia. `ForwardEuler`, `RungeKutta2`,
  `RungeKutta4`, `BackwardEuler`, and `CrankNicolson` now use
  `leto::Array1<T>` as their `TimeIntegrator::State`, with slice-based
  scaled-update, sum, and convergence-norm helpers over Eunomia
  `NumericElement`/`RealField`. Adaptive and variable time-step controllers
  now use Eunomia constants and `FloatElement::powf` instead of
  `num_traits::Float`/`FromPrimitive`. Evidence: `cargo fmt -p cfd-core
  --check`, `cargo check -p cfd-core --no-default-features`, `cargo check -p
  cfd-core --no-default-features --tests`, `cargo clippy -p cfd-core
  --no-default-features --lib -- -D warnings`, `cargo clippy -p cfd-core
  --all-targets -- -D warnings`, and focused `cargo nextest run -p cfd-core
  --no-default-features time forward_euler runge_kutta backward_euler crank
  adaptive_controller variable_controller --status-level fail` (17/17). A
  direct scan over `crates/cfd-core/src/compute/time` shows no nalgebra,
  `DVector`, ndarray, num-traits, rayon, tokio, wgpu, or cuda residue.
- **cfd-core/cfd-3d**: Migrated the solver-configuration scalar boundary to
  Eunomia for `SolverConfig` and downstream holders. `cfd_core::compute::
  solver::config` now defines `SolverConfig`, convergence, numerical, linear,
  network, and builder configuration types on `eunomia::RealField`; the
  spectral solver no longer imports or names `NalgebraRealField`; and affected
  cfd-2d FDM/SIMPLE/pressure-velocity/vorticity, cfd-3d FEM, and
  cfd-validation MMS Richardson holders carry the Eunomia bound needed to
  consume the migrated config. `SolverConfiguration` and the remaining
  `cfd_core::compute::solver::{traits,convergence,direct,iterative,monitor}`
  surfaces intentionally retain nalgebra bounds because they still depend on
  `abstractions::Problem<T>`. Evidence: `cargo fmt -p cfd-core -p cfd-2d -p
  cfd-3d -p cfd-validation --check`, `cargo check -p cfd-core
  --no-default-features`, `cargo check -p cfd-2d --no-default-features --lib`,
  `cargo check -p cfd-3d --no-default-features --lib`, `cargo check -p cfd-3d
  --no-default-features --test property_tests --test domain_solver_validation
  --test poisson_validation`, `cargo check -p cfd-validation
  --no-default-features --lib`, scoped clippy for cfd-core/cfd-2d/cfd-3d/
  cfd-validation, focused `cargo nextest run -p cfd-core
  --no-default-features solver_config --status-level fail` (2/2), focused
  `cargo nextest run -p cfd-2d --no-default-features fdm pressure_velocity
  simple vorticity --status-level fail` (41/41), and focused `cargo nextest
  run -p cfd-3d --no-default-features spectral poisson fem --status-level
  fail` (89/89).
- **cfd-3d**: Migrated the spectral Poisson solve boundary to Leto/Leto-ops.
  `PoissonSolver` now builds the global Chebyshev Laplacian as `leto::Array2`,
  forms tensor-product operators with `leto_ops::MatrixProduct::kron`, solves
  the dense system through `leto_ops::MatrixSolve`, and returns `leto::Array1`.
  `PoissonProblem` no longer requires a nalgebra scalar bound; the remaining
  `NalgebraRealField` import in `spectral::solver` is inherited from
  `cfd_core::compute::solver::SolverConfig`. `cargo fmt -p cfd-3d --check`,
  `cargo check -p cfd-3d --no-default-features --lib`, `cargo check -p cfd-3d
  --no-default-features --test robustness_tests --test poisson_validation`,
  `cargo clippy -p cfd-3d --no-default-features --lib --test
  robustness_tests --test poisson_validation -- -D warnings`, and focused
  `cargo nextest run -p cfd-3d --no-default-features poisson --status-level
  fail` pass.
- **cfd-3d**: Migrated the spectral Chebyshev basis to Leto/Eunomia.
  `ChebyshevPolynomial` now stores differentiation operators as `leto::Array2`,
  applies first/second derivatives over `leto::Array1`, and uses
  `eunomia::RealField`; co-located Chebyshev tests, `lib.rs` Chebyshev tests,
  and robustness Chebyshev tests no longer construct nalgebra `DVector`s.
  The follow-on Poisson slice replaces the former nalgebra LU/Kronecker solve
  seam with Leto/Leto-ops. `cargo check -p
  cfd-3d --no-default-features --lib`, scoped `cargo clippy -p cfd-3d
  --no-default-features --lib --test robustness_tests --test
  poisson_validation -- -D warnings`, focused `cargo nextest run -p cfd-3d
  --no-default-features chebyshev --status-level fail`, and focused `cargo
  nextest run -p cfd-3d --no-default-features poisson --status-level fail`
  pass.
- **cfd-3d**: Completed the level-set Atlas provider slice. `level_set` solver,
  advection, and WENO code now use `leto::geometry::Vector3` and a local
  `LevelSetScalar` seam backed by `eunomia::RealField`/`FloatElement`; the dead
  crate-level `CfdScalar` nalgebra-bound trait was removed. Level-set tests and
  robustness-test level-set call sites now construct Leto vectors. `cargo check
  -p cfd-3d --no-default-features --lib`, `cargo check -p cfd-3d
  --no-default-features --test level_set_tests`, `cargo check -p cfd-3d
  --no-default-features --test robustness_tests`, `cargo clippy -p cfd-3d
  --no-default-features --lib --test level_set_tests --test robustness_tests
  -- -D warnings`, and focused `cargo nextest run -p cfd-3d
  --no-default-features level_set --status-level fail` pass. A direct
  level-set forbidden-provider scan now finds no `nalgebra`, `DMatrix`,
  `DVector`, `CfdScalar`, `crate::scalar`, `num_traits`, `ndarray`, `rayon`,
  `tokio`, `rustfft`, `wgpu`, or `cuda` residue.
- **cfd-3d**: Completed the VOF Atlas provider slice. Cavitation dense fields
  now use `CavitationField = leto::Array2<f64>` for pressure, density,
  volume-fraction, inception, damage, bubble-radius, nuclei, and
  sonoluminescence outputs; row-major Leto field offsets are explicit and
  separated from the VOF flat alpha offset. `vof::scalar::VofScalar` now uses
  `eunomia::RealField` instead of `nalgebra::RealField`. `cargo check -p
  cfd-3d --no-default-features --lib`, `cargo check -p cfd-3d
  --no-default-features --test cavitation_solver_validation`, `cargo clippy -p
  cfd-3d --no-default-features --lib --test cavitation_solver_validation -- -D
  warnings`, focused `cargo nextest run -p cfd-3d --no-default-features
  cavitation --status-level fail`, and focused `cargo nextest run -p cfd-3d
  --no-default-features vof --status-level fail` pass. A direct VOF forbidden
  provider scan now finds no `nalgebra`, `DMatrix`, `DVector`, `num_traits`,
  `ndarray`, `rayon`, `tokio`, `rustfft`, `wgpu`, or `cuda` residue.
- **cfd-3d**: Migrated the VOF cavitation velocity-provider surface from
  nalgebra `Vector3` to Leto. `vof::cavitation_solver`,
  `vof::bubble_dynamics`, and `tests/cavitation_solver_validation.rs` now use
  `leto::geometry::Vector3` for cavitation velocity fields and bubble dynamics
  update calls; the cavitation solver copies directly into the Leto-owned
  `VofSolver` velocity buffer. `cargo check -p cfd-3d --no-default-features
  --lib`, `cargo check -p cfd-3d --no-default-features --test
  cavitation_solver_validation`, `cargo clippy -p cfd-3d --no-default-features
  --lib --test cavitation_solver_validation -- -D warnings`, and focused
  `cargo nextest run -p cfd-3d --no-default-features cavitation --status-level
  fail` pass. The later VOF Atlas provider slice closes the cavitation dense
  storage residue from this item.
- **cfd-3d**: Migrated the non-cavitation VOF vector-provider surface from
  nalgebra `Vector3` to Leto. VOF velocity storage, interface normals,
  initialization geometry, PLIC plane normals, and VOF tests now use
  `leto::geometry::Vector3`. `cargo check -p cfd-3d --no-default-features --test
  vof_tests`, `cargo check -p cfd-3d --no-default-features --test
  robustness_tests`, `cargo clippy -p cfd-3d --no-default-features --lib
  --test vof_tests -- -D warnings`, and focused `cargo nextest run -p cfd-3d
  --no-default-features vof --status-level fail` pass. The later VOF Atlas
  provider slice closes the cavitation dense-storage residue from this item.
- **cfd-validation**: Closed direct nalgebra `Vector2` ownership in validation
  source and tests. Navier-Stokes MMS, momentum/angular-momentum conservation,
  vorticity-stream benchmark sampling, and MMS/conservation test consumers now
  use `leto::geometry::Vector2` with indexed component access. `cargo check -p
  cfd-validation --no-default-features --tests`, `cargo clippy -p
  cfd-validation --no-default-features --lib -- -D warnings`, and focused
  `cargo nextest run -p cfd-validation --no-default-features manufactured mms
  conservation taylor momentum angular --status-level fail` pass. Remaining
  validation nalgebra ownership is dense `DMatrix`, 3D `Vector3`, and geometry
  point/vector surfaces.
- **cfd-2d**: Closed direct nalgebra `Vector2` ownership across cfd-2d source,
  tests, examples, and benches. Vorticity-stream, PISO corrector, LBM,
  immersed-boundary, `blood_venturi`, and `solver_benchmarks` now use
  `leto::geometry::Vector2` with indexed component access, and the downstream
  cfd-validation vorticity-stream benchmark consumes the Leto velocity-field
  contract. `cargo check -p cfd-2d --no-default-features --examples --benches`,
  `cargo check -p cfd-validation --no-default-features`, scoped clippy for the
  touched cfd-2d executable targets and cfd-validation library, and focused
  `cargo nextest run -p cfd-2d --no-default-features vorticity corrector lbm
  immersed --status-level fail` pass. Remaining cfd-2d nalgebra ownership is
  scalar `RealField`, boundary `Vector3`, dense `DVector`/`DMatrix`, and
  nalgebra-sparse storage.
- **cfd-2d**: Migrated the pressure-velocity/SIMPLEC/PIMPLE vector workspace
  family from nalgebra `Vector2` to Leto. `fields`, `pressure_velocity`, and
  `simplec_pimple` velocity workspaces, Rhie-Chow caches, correction buffers,
  face velocities, and SIMPLEC/PIMPLE validation setup now use
  `leto::geometry::Vector2` with indexed component access. `cargo check -p
  cfd-2d --no-default-features`, `cargo clippy -p cfd-2d --no-default-features
  --lib -- -D warnings`, and focused `cargo nextest run -p cfd-2d
  --no-default-features pressure_velocity simplec pimple rhie --status-level
  fail` pass. Remaining nalgebra in this family is scalar `RealField`,
  `Vector3` boundary inputs, and `DVector` linear-solver storage.
- **cfd-2d**: Moved the Ghia pressure-correction helper test off direct
  nalgebra vector construction. `tests/ghia_cavity_simplec_validation.rs` now
  uses `leto::geometry::Vector2` for the local divergence-free velocity helper.
  `cargo check -p cfd-2d --no-default-features --test
  ghia_cavity_simplec_validation` and focused `cargo nextest run -p cfd-2d
  --no-default-features test_pressure_correction_basic --status-level fail`
  pass. Remaining cfd-2d test nalgebra ownership is in
  `tests/simplec_pimple_validation.rs` and the linked production vector
  workspaces.
- **cfd-2d**: Migrated the turbulence validation scalar contract from direct
  nalgebra bounds to Eunomia. `physics::turbulence::validation` no longer
  imports or bounds `nalgebra::RealField`; validation scalar contracts now use
  `eunomia::RealField`/`FloatElement`, while validation arrays/vectors remain
  on Leto `Array2`/`Vector2`. `cargo check -p cfd-2d --no-default-features`,
  `cargo clippy -p cfd-2d --no-default-features --lib -- -D warnings`, and
  focused `cargo nextest run -p cfd-2d --no-default-features turbulence
  validation --status-level fail` pass. Broader cfd-2d direct
  `nalgebra`/`nalgebra-sparse` ownership remains outside this validation seam.
- **cfd-2d**: Removed direct package `num-traits` ownership. The closing pass
  migrated turbulence validation bounds to Eunomia, removed the immersed-boundary
  test's unnecessary `ToPrimitive` conversion over an already-`f64` field,
  removed concrete `f64` Ghia validation `num_traits` imports/bounds, and
  deleted `num-traits.workspace = true` from `crates/cfd-2d/Cargo.toml`. `cargo
  check -p cfd-2d --no-default-features --tests`, `cargo clippy -p cfd-2d
  --no-default-features --lib -- -D warnings`, and focused `cargo nextest run
  -p cfd-2d --no-default-features turbulence validation immersed ghia
  --status-level fail` pass. `cargo tree -p cfd-2d --depth 1` shows no direct
  `num-traits`; remaining cfd-2d provider work is nalgebra/nalgebra-sparse and
  Leto/Gaia storage/vector migration.
- **cfd-2d**: Migrated the momentum setup/boundary scalar-provider seam to
  Eunomia. `physics::momentum::setup` no longer imports or bounds direct
  `num_traits::{FromPrimitive,ToPrimitive}`. `physics::momentum::boundary` and
  `boundary::directional` now route zero/one row entries, zero-gradient
  comparisons, quadratic wall-extrapolation constants, and corner consistency
  absolute-value checks through the crate-local Eunomia scalar adapter instead
  of direct `num_traits::FromPrimitive`, `T::from_f64`, `T::zero()`, `T::one()`,
  or scalar `.abs()` usage. `cargo check -p cfd-2d --no-default-features`,
  `cargo clippy -p cfd-2d --no-default-features --lib -- -D warnings`, and
  focused `cargo nextest run -p cfd-2d --no-default-features boundary
  --status-level fail` pass. Broader cfd-2d direct `num-traits` residue remains
  in Ghia validation tests, immersed-boundary test code, turbulence validation,
  and f64-only/test scalar surfaces.
- **cfd-2d**: Migrated the momentum interpolation scalar-provider seam to
  Eunomia. `physics::momentum::interpolation` now routes the Rhie-Chow tiny
  diagonal floor, zero face-array initialization, and harmonic face coefficient
  factor through the crate-local Eunomia scalar adapter instead of direct
  `num_traits::FromPrimitive`, `T::from_f64`, `T::zero()`, or `T::one()` usage.
  `cargo check -p cfd-2d --no-default-features`, `cargo clippy -p cfd-2d
  --no-default-features --lib -- -D warnings`, and focused `cargo nextest run
  -p cfd-2d --no-default-features
  coefficient_aware_interpolation_matches_exact_reference --status-level fail`
  pass. Broader cfd-2d direct `num-traits` residue remains outside this
  interpolation slice.
- **cfd-2d**: Migrated the problem/streamtube scalar-provider seam to Atlas
  providers. `problem.rs` now stores incompressible problem and solution
  velocity fields with `leto::geometry::Vector2` and routes initial pressure,
  velocity-magnitude maxima, and pressure maxima through the crate-local
  Eunomia scalar adapter. `physics::streamtube::partitioning` now routes
  constants, absolute values, and square roots through
  `eunomia::{FloatElement,NumericElement}` instead of direct
  `num_traits::{Float,FromPrimitive}` usage. `cargo check -p cfd-2d
  --no-default-features`, `cargo clippy -p cfd-2d --no-default-features --lib
  -- -D warnings`, and focused `cargo nextest run -p cfd-2d
  --no-default-features problem streamtube separating --status-level fail`
  pass. `problem.rs` still inherits `nalgebra::RealField` from upstream
  `cfd-core` boundary/fluid contracts.
- **cfd-2d**: Migrated the FVM scalar and vector-provider surface to Atlas
  providers. `solvers::fvm::{config,flux,geometry,solver}` now routes FVM
  defaults, face generation constants, residual magnitudes, max terms, flux
  Peclet calculations, nonfinite diffusion validation, and module test
  absolute-value assertions through the crate-local Eunomia scalar adapter and
  `eunomia::{FloatElement,NumericElement}` instead of direct `num_traits`
  usage. FVM face geometry and solver velocity vectors now use
  `leto::geometry::Vector2`; power-law and hybrid flux calculators preserve
  diffusion in `T` instead of narrowing through `f64`. `cargo check -p cfd-2d
  --no-default-features`, `cargo clippy -p cfd-2d --no-default-features --lib
  -- -D warnings`, and focused `cargo nextest run -p cfd-2d
  --no-default-features fvm --status-level fail` pass. Broader cfd-2d direct
  `num-traits` residue remains outside this FVM slice.
- **cfd-2d**: Migrated the vorticity-stream scalar-provider surface to
  Eunomia. `physics::vorticity_stream` now routes default tolerances,
  timestep/SOR constants, zero/one identities, Laplacian gradient factors, SOR
  residuals, upwind sign checks, velocity-recovery denominators,
  boundary-vorticity factors, convergence checks, and continuity-test
  absolute-value assertions through the crate-local Eunomia scalar adapter and
  `eunomia::{FloatElement,NumericElement}` instead of direct
  `num_traits::FromPrimitive` usage. `cargo check -p cfd-2d
  --no-default-features`, `cargo clippy -p cfd-2d --no-default-features --lib
  -- -D warnings`, and focused `cargo nextest run -p cfd-2d
  --no-default-features vorticity --status-level fail` pass. Broader cfd-2d
  direct `num-traits` residue remains outside this vorticity-stream slice.
- **cfd-2d**: Migrated the acoustic-drift scalar-provider surface to Eunomia.
  `physics::acoustics::gorkov` and `solvers::drift_diffusion_2d` now route ARF
  material constants, compressibility identities, contrast-factor constants,
  standing-wave sine evaluation, drift-diffusion zero/one identities,
  under-relaxation constants, Patankar upwind max terms, pivot-threshold
  checks, and convergence residual magnitudes through the crate-local Eunomia
  scalar adapter and `eunomia::{FloatElement,NumericElement}` instead of direct
  `num_traits::{Float,FromPrimitive}` usage. `cargo check -p cfd-2d
  --no-default-features`, `cargo clippy -p cfd-2d --no-default-features --lib
  -- -D warnings`, focused `cargo nextest run -p cfd-2d
  --no-default-features gorkov --status-level fail`, and focused `cargo
  nextest run -p cfd-2d --no-default-features drift_diffusion --status-level
  fail` pass. Broader cfd-2d direct `num-traits` residue remains outside this
  acoustic-drift slice.
- **cfd-2d**: Migrated the TVD scheme scalar-provider family to Eunomia.
  `schemes::tvd::{mod,muscl,quick}` now routes TVD limiter identities,
  min/max/absolute-value operations, MUSCL slope ratios, QUICK interpolation
  constants, and MUSCL2/MUSCL3 boundary constants through the crate-local
  Eunomia scalar adapter and `eunomia::{FloatElement,NumericElement}` instead
  of direct `num_traits::{FromPrimitive,ToPrimitive}` usage. The limiter path
  now stays generic in `T` and no longer round-trips the limiter ratio through
  `f64`. `cargo check -p cfd-2d --no-default-features`, `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings`, focused `cargo nextest
  run -p cfd-2d --no-default-features tvd --status-level fail`, and focused
  `cargo nextest run -p cfd-2d --no-default-features muscl --status-level
  fail` pass. Broader cfd-2d direct `num-traits` residue remains outside this
  TVD slice.
- **cfd-2d**: Migrated the LBM scalar-provider surface to Eunomia.
  `solvers::lbm` now routes solver defaults, low-Mach validation, grid-index
  conversion, boundary density/velocity reconstruction, passive-scalar
  boundary zeroing, streaming push-zeroing, macroscopic/lattice assertions,
  Carreau-Yasuda assertions, and Shan-Chen pseudopotential/force constants
  through the crate-local Eunomia scalar adapter and
  `eunomia::{FloatElement,NumericElement,CastFrom}` instead of direct
  `num_traits::{Float,FromPrimitive}` usage. `cargo check -p cfd-2d
  --no-default-features`, `cargo clippy -p cfd-2d --no-default-features --lib
  -- -D warnings`, and focused `cargo nextest run -p cfd-2d
  --no-default-features lbm --status-level fail` pass. Broader cfd-2d direct
  `num-traits` residue remains outside this LBM slice.
- **cfd-2d**: Migrated the legacy SIMPLE scalar-provider surface to Eunomia.
  `solvers::simple::{algorithm,momentum,pressure}` now routes Patankar default
  constants, workspace zero initialization, stagnant-cell safeguards,
  pressure-correction identities, neighbor-count conversion, finite checks,
  and module test absolute-value assertions through the crate-local Eunomia
  scalar adapter and `eunomia::{FloatElement,NumericElement}` instead of direct
  `num_traits::{FromPrimitive,ToPrimitive}` usage. `cargo check -p cfd-2d
  --no-default-features`, `cargo clippy -p cfd-2d --no-default-features --lib
  -- -D warnings`, and focused `cargo nextest run -p cfd-2d
  --no-default-features simple --status-level fail` pass. Broader cfd-2d
  direct `num-traits` residue remains outside this SIMPLE slice.
- **cfd-2d**: Migrated the Venturi scalar-provider surface to Eunomia.
  `solvers::venturi_flow` now routes ISO/default constants, beta clamping,
  domain bounds, analytical Bernoulli/viscous constants, energy-dissipation
  guards, stretched-grid index conversion, generic sine/min/max dispatch,
  diagnostic f64 conversion, throat-column selection, inlet/outlet/throat
  reductions, and validation error metrics through the crate-local Eunomia
  scalar adapter and `eunomia::{FloatElement,NumericElement}` instead of
  direct `num_traits::{Float,FromPrimitive,ToPrimitive}` usage. `cargo check
  -p cfd-2d --no-default-features`, `cargo clippy -p cfd-2d
  --no-default-features --lib -- -D warnings`, and focused `cargo nextest run
  -p cfd-2d --no-default-features venturi --status-level fail` pass. Broader
  cfd-2d direct `num-traits` residue remains outside this Venturi slice.
- **cfd-2d**: Migrated the branching-flow scalar-provider family to Eunomia.
  `solvers::{bifurcation_flow,n_furcation_flow}` now routes geometry
  constants, branch-index conversion, zero identities, generic sin/cos
  dispatch, flux accumulation, and mass-balance absolute-value normalization
  through the crate-local Eunomia scalar adapter and
  `eunomia::{FloatElement,NumericElement}` instead of direct
  `num_traits::{Float,FromPrimitive,ToPrimitive}` usage. `cargo check -p
  cfd-2d --no-default-features`, `cargo clippy -p cfd-2d
  --no-default-features --lib -- -D warnings`, focused `cargo nextest run -p
  cfd-2d --no-default-features bifurcation --status-level fail`, and focused
  `cargo nextest run -p cfd-2d --no-default-features n_furcation
  --status-level fail` pass. Broader cfd-2d direct `num-traits` residue
  remains outside this branching-family slice.
- **cfd-2d**: Migrated the cross-junction scalar-provider surface to Eunomia.
  `solvers::cross_junction_flow` now routes geometry constants, port-flux
  zero identities, mass-balance absolute-value normalization, and bounding-box
  test checks through the crate-local Eunomia scalar adapter and
  `eunomia::{FloatElement,NumericElement}` instead of direct
  `num_traits::{Float,FromPrimitive,ToPrimitive}` usage. `cargo check -p
  cfd-2d --no-default-features`, `cargo clippy -p cfd-2d
  --no-default-features --lib -- -D warnings`, and focused `cargo nextest run
  -p cfd-2d --no-default-features cross_junction --status-level fail` pass.
  Broader cfd-2d direct `num-traits` residue remains outside this
  cross-junction slice.
- **cfd-2d**: Migrated the Poiseuille scalar-provider surface to Eunomia.
  Poiseuille configuration defaults, grid construction, tridiagonal
  workspaces, harmonic-mean constants, shear-rate magnitudes, viscosity
  residual norms, Thomas pivot thresholds, error diagnostics, and module test
  absolute-value checks now route through the crate-local Eunomia scalar
  adapter and `eunomia::{FloatElement,NumericElement}` instead of direct
  `num_traits::{Float,FromPrimitive}` usage. `cargo check -p cfd-2d
  --no-default-features`, `cargo clippy -p cfd-2d --no-default-features --lib
  -- -D warnings`, and focused `cargo nextest run -p cfd-2d
  --no-default-features poiseuille --status-level fail` pass. Broader cfd-2d
  direct `num-traits` residue remains outside this Poiseuille slice.
- **cfd-2d**: Migrated the serpentine and scalar-transport scalar-provider
  surfaces to Eunomia. Serpentine geometry, analytical mixing, validator,
  discretized solver, and scalar transport now route scalar constants,
  identities, Peclet/mixing calculations, Fourier-mode index conversion,
  transport relaxation constants, coefficient max/abs operations, outlet
  averaging, validator diagnostics, and module test absolute-value checks
  through the crate-local Eunomia scalar adapter and
  `eunomia::{FloatElement,NumericElement}` instead of direct
  `num_traits::{Float,FromPrimitive,ToPrimitive}` usage. `cargo check -p
  cfd-2d --no-default-features`, `cargo clippy -p cfd-2d
  --no-default-features --lib -- -D warnings`, and focused `cargo nextest run
  -p cfd-2d --no-default-features serpentine --status-level fail` pass.
  Broader cfd-2d direct `num-traits` residue remains outside this
  serpentine/scalar-transport slice.
- **cfd-2d**: Migrated the FDM scalar-provider surface to Eunomia.
  `solvers::fdm::{advection_diffusion,config,diffusion,linear_solver,poisson}`
  now routes scalar constants, identities, singular-diagonal guards,
  Gauss-Seidel residual magnitudes, finite-difference constants, and MMS test
  absolute-value checks through the crate-local Eunomia scalar adapter and
  `eunomia::{FloatElement,NumericElement}` instead of direct
  `num_traits::FromPrimitive` usage. The existing FDM MMS source tests are now
  wired into the test graph, and the advection-diffusion stencil sign now
  matches the documented steady operator. `cargo check -p cfd-2d
  --no-default-features`, `cargo clippy -p cfd-2d --no-default-features --lib
  -- -D warnings`, and focused `cargo nextest run -p cfd-2d
  --no-default-features fdm --status-level fail` pass. Broader cfd-2d direct
  `num-traits` residue remains outside this FDM slice.
- **cfd-2d**: Migrated the network scalar-provider surface to Eunomia.
  `network::{build,channel,coupled,postprocess,projection,reference,solve,
  types}` now routes scalar constants, identities, abs/min/max/finite checks,
  Anderson relaxation constants, projection summaries, channel diagnostics,
  and f64 reporting through the crate-local Eunomia scalar adapter and
  `eunomia::{FloatElement,NumericElement}` instead of direct
  `num_traits::{Float,FromPrimitive,ToPrimitive}` usage. `cargo check -p
  cfd-2d --no-default-features`, `cargo clippy -p cfd-2d
  --no-default-features --lib -- -D warnings`, and focused `cargo nextest run
  -p cfd-2d --no-default-features network --status-level fail` pass. Broader
  cfd-2d direct `num-traits` residue remains outside this network slice.
- **cfd-2d**: Migrated the NS-FVM scalar-provider surface to Eunomia.
  `solvers::ns_fvm` field storage, SIMPLE solver core, momentum equations,
  pressure correction, Rhie-Chow mass-flux correction, and staggered velocity
  interpolation now route scalar constants, identities, abs/min/max/sqrt/
  finite checks, interpolation index conversion, turbulence seeds, and f64
  diagnostics through the crate-local Eunomia scalar adapter and
  `eunomia::{FloatElement,NumericElement}` instead of direct
  `num_traits::{Float,FromPrimitive,ToPrimitive}` usage. `cargo check -p
  cfd-2d --no-default-features`, `cargo clippy -p cfd-2d
  --no-default-features --lib -- -D warnings`, and focused `cargo nextest run
  -p cfd-2d --no-default-features ns_fvm --status-level fail` pass. Broader
  cfd-2d direct `num-traits` residue remains outside this NS-FVM slice.
- **cfd-2d**: Migrated the shared continuity-residual scalar-provider surface
  to Eunomia. `solvers::continuity` now routes zero identities,
  central-difference constants, and absolute-value reductions through the
  crate-local Eunomia scalar adapter and `eunomia::NumericElement` instead of
  direct `num_traits::FromPrimitive`, `T::from_*`, `T::zero()`, or scalar
  method `.abs()` usage. `cargo check -p cfd-2d --no-default-features`,
  `cargo clippy -p cfd-2d --no-default-features --lib -- -D warnings`, and
  focused `cargo nextest run -p cfd-2d --no-default-features continuity
  --status-level fail` pass. Broader cfd-2d direct `num-traits` residue remains
  outside this continuity slice.
- **cfd-2d**: Migrated the SIMPLEC/PIMPLE scalar-provider surface to
  Eunomia. SIMPLEC/PIMPLE defaults, validation bounds, residual identities,
  adaptive-step constants, pressure-correction caches, Rhie-Chow face arrays,
  pressure extrapolation averages, boundary zeroing, PIMPLE unity relaxation,
  and focused Ghia validation interpolation now route through the crate-local
  Eunomia scalar adapter or Eunomia-backed validation helpers instead of direct
  `num_traits::{FromPrimitive,ToPrimitive}`, `T::from_*`, `.to_f64()`,
  `T::zero()`, or `T::one()` usage. `cargo check -p cfd-2d
  --no-default-features`, `cargo clippy -p cfd-2d --no-default-features --lib
  -- -D warnings`, and focused `cargo nextest run -p cfd-2d
  --no-default-features simplec_pimple --status-level fail` pass. Broader
  cfd-2d direct `num-traits` residue remains outside this SIMPLEC/PIMPLE
  slice.
- **cfd-2d**: Migrated the PISO scalar-provider surface to Eunomia. PISO
  defaults, convergence tolerances, residual monitors, predictor/corrector
  workspaces, relaxation identities, hybrid differencing max/abs operations,
  Rhie-Chow tiny thresholds, duration thresholds, and logging conversions now
  route through the crate-local Eunomia scalar adapter and
  `eunomia::{FloatElement,NumericElement}` instead of direct
  `num_traits::{FromPrimitive,ToPrimitive}`, `T::from_*`, `.to_f64()`,
  `T::zero()`, or `T::one()` usage. `cargo check -p cfd-2d
  --no-default-features`, `cargo clippy -p cfd-2d --no-default-features --lib
  -- -D warnings`, and focused `cargo nextest run -p cfd-2d
  --no-default-features piso --status-level fail` pass. Broader cfd-2d direct
  `num-traits` residue remains outside this PISO slice.
- **cfd-2d**: Migrated the pressure-velocity scalar-provider surface to
  Eunomia. Coefficient defaults, pressure-velocity validation bounds, solver
  workspaces, pressure-correction assembly/scatter, and Rhie-Chow coefficient
  guards now route zero/one identities through the crate-local Eunomia scalar
  adapter instead of direct `T::zero()`/`T::one()` construction. Finite checks
  in the touched pressure-velocity paths use `eunomia::NumericElement`.
  `cargo check -p cfd-2d --no-default-features`, `cargo clippy -p cfd-2d
  --no-default-features --lib -- -D warnings`, and focused `cargo nextest run
  -p cfd-2d --no-default-features pressure_velocity --status-level fail` pass.
  Broader cfd-2d direct `num-traits` residue remains outside this
  pressure-velocity slice.
- **cfd-2d**: Migrated the time-integration scalar-provider surface to
  Eunomia. Explicit RK, implicit fixed-point, BDF/Adams-Bashforth, CFL,
  Richardson-error, and adaptive-step paths now route scalar constants and
  abs/min/max/powf dispatch through the crate-local Eunomia scalar adapter
  instead of direct `num_traits::{FromPrimitive,ToPrimitive}` construction,
  `T::zero()`, `T::one()`, or fallback `unwrap_or(T::one())` calls. `cargo
  check -p cfd-2d --no-default-features`, `cargo clippy -p cfd-2d
  --no-default-features --lib -- -D warnings`, and focused `cargo nextest run
  -p cfd-2d --no-default-features time --status-level fail` pass. Broader
  cfd-2d direct `num-traits` residue remains outside this time-integration
  slice.
- **cfd-2d**: Migrated the finite-volume discretization scalar-provider
  surface to Eunomia. Convection scheme strategies, QUICK/MUSCL extended
  stencils, and QUICK momentum corrections now route scalar constants,
  zero/one values, absolute values, and max operations through the crate-local
  Eunomia scalar adapter instead of direct `num_traits::FromPrimitive`
  construction or `T::zero()`/`T::one()` constructors. `cargo check -p cfd-2d
  --no-default-features`, `cargo clippy -p cfd-2d --no-default-features --lib
  -- -D warnings`, and focused `cargo nextest run -p cfd-2d
  --no-default-features discretization --status-level fail` pass. Broader
  cfd-2d direct `num-traits` residue remains outside this discretization
  slice.
- **cfd-2d**: Migrated the spatial/WENO scalar-provider surface to Eunomia.
  Central/upwind schemes, WENO5/WENO9 helpers, WENO-Z reconstruction, TVD/QUICK
  limiter constants, and WENO-Z momentum coefficient corrections now route
  scalar constants, zero values, powers, absolute values, and relaxation
  factors through the crate-local Eunomia scalar adapter instead of direct
  `num_traits::{FromPrimitive,ToPrimitive}` construction. `cargo check -p
  cfd-2d --no-default-features`, `cargo clippy -p cfd-2d
  --no-default-features --lib -- -D warnings`, and focused `cargo nextest run
  -p cfd-2d --no-default-features weno --status-level fail` pass. Broader
  cfd-2d direct `num-traits` residue remains outside this scheme/momentum
  slice.
- **cfd-2d**: Migrated the field-container scalar-provider surface to
  Eunomia. `fields::{constants,Field2D,SimulationFields}` now routes scalar
  constants, zero/one values, `Field2D::zeros`, `SimulationFields` defaults
  and reset values, velocity-magnitude square root/max, and Reynolds-number
  averaging through the crate-local Eunomia scalar adapter instead of direct
  `num_traits`/`FromPrimitive`/`Float` calls. `cargo check -p cfd-2d
  --no-default-features`, `cargo clippy -p cfd-2d --no-default-features --lib
  -- -D warnings`, and focused `cargo nextest run -p cfd-2d
  --no-default-features piso_algorithm --status-level fail` pass. Broader
  cfd-2d direct `num-traits` residue remains outside `fields.rs`.
- **cfd-3d**: Removed direct `num-traits` ownership from the 3D package. The
  final root `lib.rs` Chebyshev assertions now use inherent `f64::abs`, and
  the manifest no longer declares `num-traits`. Focused `cargo nextest run -p
  cfd-3d --no-default-features chebyshev --status-level fail` passes 20/20,
  `cargo clippy -p cfd-3d --no-default-features --lib -- -D warnings`
  passes, and the direct source/test/example/manifest scan is clean. The
  remaining `num-traits` graph is transitive through provider/nalgebra paths.
- **cfd-3d**: Migrated the Bifurcation and Venturi scalar-provider seams to
  Eunomia. Bifurcation geometry/solver/validation and Venturi
  solver/validation paths no longer use direct
  `num_traits::{Float,FromPrimitive,ToPrimitive}` or `num_traits::Float::*`;
  scalar constants, zero/one values, elementary math, f64/usize
  construction, pressure coefficients, Richardson/GCI metrics, Picard
  viscosity-change checks, pressure slicing, flow extraction, wall-shear
  estimates, and validation thresholds route through `cfd-3d::scalar`.
  Focused `cargo nextest run -p cfd-3d --no-default-features bifurcation
  venturi --status-level fail` passes 35/35 and `cargo clippy -p cfd-3d
  --no-default-features --lib -- -D warnings` passes. Direct cfd-3d
  `num-traits` ownership is now closed.
- **cfd-3d**: Migrated the Trifurcation scalar-provider seam to Eunomia.
  Trifurcation geometry SDF trigonometry, Murray-law powers, solver defaults,
  Picard viscosity-change checks, boundary-flow extraction, pressure drops,
  wall-shear estimates, mass conservation, and validation thresholds no longer
  use direct `num_traits::{Float,FromPrimitive,ToPrimitive}` or
  `num_traits::Float::*`; they route scalar constants and elementary math
  through `cfd-3d::scalar`. Focused `cargo nextest run -p cfd-3d
  --no-default-features trifurcation --status-level fail` passes 6/6 and
  `cargo clippy -p cfd-3d --no-default-features --lib -- -D warnings`
  passes.
- **cfd-3d**: Migrated the Serpentine scalar-provider seam to Eunomia.
  Serpentine configuration defaults, Dean-number analysis, pressure-drop
  validation, Picard viscosity-change checks, and solution defaults no longer
  use direct `num_traits::{Float,FromPrimitive,ToPrimitive}` or
  `num_traits::Float::*`; they route scalar constants, zero/one values,
  abs/max/sqrt, integer-to-scalar construction, and f64 diagnostics through
  `cfd-3d::scalar`. Focused `cargo nextest run -p cfd-3d
  --no-default-features serpentine --status-level fail` passes 7/7 and
  `cargo clippy -p cfd-3d --no-default-features --lib -- -D warnings`
  passes.
- **cfd-3d**: Migrated the IBM scalar-provider seam to Eunomia. IBM forcing,
  Roma/Peskin interpolation kernels, and IBM stencil/grid-index math no longer
  use direct `num_traits::{FromPrimitive,ToPrimitive}` or
  `num_traits::Float::*`; they route scalar constants, zero/one values,
  abs/max/sqrt/cos/floor, and index scalar construction through
  `cfd-3d::scalar`. Scaled grid coordinates that are non-finite or outside
  `isize` range now return a typed configuration error instead of silently
  falling back to zero. Focused `cargo nextest run -p cfd-3d
  --no-default-features ibm --status-level fail` passes 12/12 and `cargo
  clippy -p cfd-3d --no-default-features --lib -- -D warnings` passes.
- **cfd-3d**: Migrated the spectral direct scalar-provider residue to
  Eunomia. Chebyshev collocation, differentiation, quadrature, interpolation,
  Poisson boundary rows, spectral configuration tolerances, and co-located
  Chebyshev tests no longer use direct `num_traits::{FromPrimitive,Float}` or
  `num_traits::Float::*`; they route scalar constants and generic math through
  `cfd-3d::scalar` or inherent f64 analytical-test methods. Focused `cargo
  nextest run -p cfd-3d --no-default-features spectral --status-level fail`
  passes 41/41 and `cargo clippy -p cfd-3d --no-default-features --lib -- -D
  warnings` passes. Spectral still retains nalgebra `DMatrix`/`DVector` and
  `RealField` surfaces for the later Leto matrix/storage migration.
- **cfd-3d**: Migrated the level-set scalar-provider seam to Eunomia.
  WENO5-Z reconstruction, SSPRK3 advection, narrow-band updates, and
  reinitialization/Godunov math now route scalar constants and elementary
  operations through the crate-local Eunomia adapter instead of direct
  `num_traits::{FromPrimitive,Float}` calls. Focused `cargo nextest run -p
  cfd-3d --no-default-features level_set --status-level fail` passes 13/13,
  `cargo clippy -p cfd-3d --no-default-features --lib -- -D warnings`
  passes, and the level-set residue scan is clean. Direct cfd-3d
  `num-traits` ownership is now closed; remaining provider work is
  nalgebra/vector replacement.
- **cfd-validation**: Removed direct `num-traits` ownership from the validation
  package. The remaining direct residues in cavity, backward-facing-step,
  disabled Poiseuille-bifurcation, and complex-boundary MMS validation now use
  Eunomia `FloatElement` or the crate-local scalar adapter instead of
  `num_traits` conversion/float traits and fallback `T::from_*` construction.
  `cargo nextest run -p cfd-validation --status-level fail` passes 431/431,
  and the focused source/test/manifest scan is clean for direct `num-traits`
  residue; remaining `num-traits` graph paths are transitive through nalgebra,
  approx, Leto/Eunomia/Gaia, and other providers.
- **cfd-2d/cfd-validation**: Migrated the pressure-velocity SIMPLEC/PIMPLE
  scalar-provider boundary and dependent 2D validation blood-flow benchmark
  wrappers to Eunomia. The touched cfd-2d solver cone no longer exposes direct
  `num_traits::{FromPrimitive,ToPrimitive}` bounds, and 2D bifurcation,
  Venturi, trifurcation, and serpentine validation benchmarks now route
  touched scalar constants and diagnostics through `cfd-validation::scalar`.
  Focused `cargo nextest run -p cfd-2d pressure_velocity simplec_pimple
  --status-level fail` passes 20/20, full `cargo nextest run -p cfd-2d
  --status-level fail` passes 567/567 with 27 skipped, and full `cargo nextest
  run -p cfd-validation --status-level fail` passes 431/431 after this
  cleanup.
- **cfd-validation**: Migrated the 3D benchmark scalar-provider wrapper seam
  to Eunomia. `BenchmarkConfig::default` and the 3D bifurcation, serpentine,
  and Venturi benchmark wrappers no longer name direct `num_traits` float or
  conversion bounds for touched scalar constants and validation predicates;
  those now route through `cfd-validation::scalar`. Focused `cargo nextest
  run -p cfd-validation bifurcation_flow_3d serpentine_flow_3d
  venturi_flow_3d` passes 3/3, and full `cargo nextest run -p cfd-validation
  --status-level fail` passes 431/431. The wrappers still carry
  `std::convert::From<f64>` because the current `cfd-3d` solver constructors
  require that upstream contract.
- **cfd-validation**: Migrated literature blood-flow validation and numerical
  linear-solver thresholds to Eunomia. `blood_flow_1d` now routes scalar
  construction, f64 diagnostics, absolute values, and Poiseuille resistance
  `powi` through the crate-local Eunomia scalar adapter instead of direct
  `num_traits`; `numerical::linear_solver` now uses `FloatElement` bounds and
  crate-local scalar constants for tolerances and convergence-rate metadata.
  Focused `cargo nextest run -p cfd-validation blood_flow literature chapman
  patankar` passes 6/6 and full `cargo nextest run -p cfd-validation` passes
  431/431 after this cleanup.
- **cfd-validation**: Migrated validation geometry scalar math to Eunomia.
  2D serpentine, 2D Venturi, and 3D serpentine validation geometries now route
  scalar constants, `sin`/`cos`, `sqrt`, `powf`, `abs`, and max selection
  through the crate-local Eunomia scalar adapter instead of direct
  `T::from_*`, `num_traits::Zero`, or RealField elementary dispatch. Focused
  `cargo nextest run -p cfd-validation geometry` passes 11/11, full
  `cargo nextest run -p cfd-validation` passes 431/431, and
  `crates/cfd-validation/src/geometry` has no targeted direct scalar-provider
  residue after this cleanup.
- **cfd-validation**: Migrated the manufactured Richardson extrapolation stack
  to Eunomia. Richardson core extrapolation, MMS validation, and validation
  analysis now route scalar constants, powers, logarithms, square roots,
  absolute values, finite checks, and f64 diagnostics through the crate-local
  Eunomia scalar adapter instead of direct `num_traits`, `nalgebra::ComplexField`,
  or test-only `nalgebra::scalar` calls. Focused
  `cargo nextest run -p cfd-validation richardson` passes 16/16, full
  `cargo nextest run -p cfd-validation` passes 431/431, and
  `crates/cfd-validation/src/manufactured` has no targeted direct provider
  residue after this cleanup.
- **cfd-validation**: Migrated the manufactured scalar/trig MMS cluster to
  Eunomia. Advection, coupled advection-diffusion, Burgers, diffusion, and
  Navier-Stokes manufactured solutions now route scalar constants,
  `sin`/`cos`, `exp`, and `sqrt` through the crate-local Eunomia scalar
  adapter instead of direct `num_traits::FromPrimitive`, cfd-core
  `SafeFromF64`, or `nalgebra::ComplexField`; the current `nalgebra`
  storage boundary remains tracked for later Leto replacement. The package
  passes focused `cargo nextest run -p cfd-validation manufactured` (50/50)
  and full `cargo nextest run -p cfd-validation` (431/431) after this cleanup.
- **cfd-validation**: Migrated Chapman-Enskog literature validation scalar
  conversion to Eunomia. Viscosity, thermal-conductivity, expected-accuracy,
  and report scalar values now route through the crate-local Eunomia adapter
  instead of direct `num_traits::{FromPrimitive,ToPrimitive}`,
  `T::from_f64`, or generic `.to_f64()` conversion. Added value-semantic tests
  for Air at 300 K transport coefficients and validation-report values; full
  `cargo nextest run -p cfd-validation` passes 431/431 after this cleanup.
- **cfd-validation**: Migrated the `time_integration` conversion seam to
  Eunomia. Runge-Kutta integrator constants, CFL/von-Neumann scalar setup,
  wave-number interpolation factors, and stability-report f64 diagnostics now
  route through the crate-local Eunomia scalar adapter instead of direct
  `num_traits`/`SafeFromF64` conversions, and local von-Neumann spatial
  operators use `eunomia::Complex` instead of `num_complex::Complex`. The
  package passes focused `cargo nextest run -p cfd-validation
  time_integration` (12/12) and full `cargo nextest run -p cfd-validation`
  (429/429) after this cleanup.
- **cfd-validation**: Migrated the `error_metrics` scalar metrics group to
  Eunomia. L1/L2/L∞, RMSE/MAE/NRMSE, relative error, and aggregate error
  statistics no longer import direct `num_traits::FromPrimitive` or call
  `T::from_f64`/`T::from_usize`; scalar construction, absolute values, and
  square roots route through the crate-local Eunomia scalar adapter. The
  package passes focused `cargo nextest run -p cfd-validation error_metrics`
  (21/21) and full `cargo nextest run -p cfd-validation` (429/429) after this
  cleanup.
- **cfd-validation**: Closed the active Eunomia scalar compile-blocker group.
  Patankar literature validation, advanced-physics MMS, multi-physics MMS,
  Reynolds-stress MMS, and the Richardson analysis result-bound seam now route
  touched scalar constants and transcendental math through the crate-local
  Eunomia scalar adapter instead of direct `num_traits`, `SafeFromF64`,
  `T::from_f64`, or `nalgebra::ComplexField` calls. The package now passes
  `cargo check -p cfd-validation` and full `cargo nextest run -p
  cfd-validation` (429/429) after this cleanup; broader provider residue
  remains tracked in other validation modules.
- **cfd-1d**: Removed direct package `num-traits` ownership and completed the
  active resistance/vascular Eunomia scalar seam. Hydraulic resistance models,
  Bessel/Womersley vascular utilities, structured-tree/bifurcation helpers,
  network blueprint/sink conversions, solver-analysis averages, and package
  tests now route direct scalar conversion/math through Eunomia/cfd-core
  provider traits. `crates/cfd-1d/Cargo.toml` no longer declares
  `num-traits`; remaining `num-traits` graph entries are transitive through
  provider/nalgebra stacks pending Leto storage replacement.
- **cfd-1d**: Migrated the solver-core scalar contract to Eunomia.
  `NetworkSolveScalar` no longer names direct `FromPrimitive`, `ToPrimitive`,
  or `num_traits::Float`; Anderson acceleration, convergence checks,
  linear-system equilibration/Jacobi preconditioning, SPD detection, residual
  norms, and f64 diagnostics now use `FloatElement`, `NumericElement`, and
  cfd-core conversion traits. The nalgebra/nalgebra-sparse storage boundary
  remains a tracked Leto replacement seam.
- **cfd-1d**: Migrated the network wrapper scalar/provider seam to Eunomia.
  `EdgeProperties::from`, characteristic-length calculation, Picard
  resistance refresh, blood hematocrit propagation, bifurcation
  phase-separation bridge values, coefficient validation, and parallel edge
  conductance now use `SafeFromF64`/`SafeFromUsize` and Eunomia
  `NumericElement` instead of direct `FromPrimitive`, `T::from_f64`,
  `T::from_usize`, `nalgebra::try_convert`, or generic `.abs()` dispatch.
  `MatrixAssembler` no longer requires `FromPrimitive`; the remaining
  `solver/core/mod.rs` `NetworkSolveScalar` num-traits compatibility contract
  is tracked as follow-up.
- **cfd-1d**: Migrated the network blueprint/sink scalar construction seam to
  Atlas provider conversions. The canonical `network_from_blueprint` path now
  routes blueprint resistances, cross-section areas, physical constants,
  boundary values, serpentine segment counts, flow probes, and blood defaults
  through `SafeFromF64`/`SafeFromUsize`, and the generic quadratic-coefficient
  policy check uses Eunomia `NumericElement::abs`. Residual direct
  `FromPrimitive` in this path is inherited from `Network::new` in
  `domain/network/wrapper.rs` and remains a tracked follow-up.
- **cfd-1d**: Migrated the domain-components provider seam from direct
  `num_traits` bounds to Eunomia. Component pressure-drop math now uses
  `NumericElement`, component constants/defaults route through
  `SafeFromF64`/`SafeFromUsize`, and channels, membranes, mixers, pumps,
  valves, sensors, and the component factory no longer import direct
  `FromPrimitive`/`Float`. `cargo check -p cfd-1d` and full `cargo nextest
  run -p cfd-1d` are clean for this slice; broad all-target clippy still has
  unrelated existing lint debt outside the component cone.
- **cfd-1d**: Migrated the next channel, branching, and network-analysis
  provider seam from direct `num_traits` bounds to Eunomia. Channel
  flow-regime classification, channel resistance/geometry/shape-factor math,
  branching network solver bounds, and pressure/flow/resistance/performance
  analysis paths now use `SafeFromF64`, `FloatElement`, and `NumericElement`.
  `cargo check -p cfd-1d` and full `cargo nextest run -p cfd-1d` are clean for
  this slice; all-target clippy still has unrelated existing lint debt.
- **cfd-1d**: Migrated the Murray's-law vascular geometry provider boundary
  from `num_traits::FromPrimitive` to Eunomia. `MurraysLaw` and
  `OptimalBifurcation` now use `FloatElement`/`NumericElement` for scalar
  construction, power operations, absolute value, and inverse cosine, with the
  upstream Eunomia `acos` contract re-verified by focused nextest; broader
  `cfd-1d` `num-traits` cleanup remains pending in other modules.
- **cfd-math**: Removed direct `num-traits` ownership from the remaining
  package residue. GPU scalar extraction now uses Eunomia `NumericElement`,
  AMG integration-test constants use Eunomia `FloatElement`, and the package
  manifest no longer declares `num-traits` directly. The touched GPU path also
  propagates Hephaestus synchronization failures and avoids a stale async
  wrapper around synchronous GPU work. Remaining `num-traits` entries are
  transitive through provider/nalgebra stacks pending later Leto/Eunomia
  migration slices.
- **cfd-core**: Removed direct `num-traits` ownership from scalar conversion
  helpers and mesh centroid construction. The touched boundary now uses
  Eunomia `FloatElement`, and the `cfd-core` manifest no longer declares
  `num-traits` directly. Remaining `num-traits` entries are transitive through
  provider/nalgebra stacks pending later Leto/Eunomia migration slices.
- **workspace**: Removed stale root workspace declarations for direct
  `wgpu 0.19` and `num-complex`. GPU provider ownership now remains expressed
  through `hephaestus-wgpu`; migrated complex-value paths use Eunomia/Apollo
  provider vocabulary. `num-complex` still appears transitively through
  provider/nalgebra stacks pending the remaining Leto/Eunomia migration
  slices.
- **cfd-core/hephaestus-wgpu**: Removed `cfd-core`'s direct `wgpu 0.19`
  dependency and migrated the crate-level GPU ABI to Hephaestus. `GpuContext`
  now acquires `hephaestus_wgpu::WgpuDevice`, and raw WGPU descriptor/buffer/
  pipeline symbols in the GPU module resolve through
  `hephaestus_wgpu::wgpu`. The cfd-core package gate is clean with clippy
  warnings denied and nextest. The remaining GPU migration work is to lift
  CFDrs-owned raw kernels onto higher-level Hephaestus kernel APIs.
- **cfd-optim**: Removed the unused direct nalgebra dev-dependency and wired
  computed SDT acoustic energy density plus CTC/RBC acoustic contrast factors
  into `SdtMetrics` while clearing all-target clippy debt in the optimizer
  package.
- **cfd-io**: Synced the crate-local `agents.md` reference with the completed
  Leto/Eunomia provider boundary: no direct internal solver/core/math
  dependencies, Leto dense checkpoint/binary arrays, Eunomia scalar bounds,
  local file-format errors, Consus-backed optional HDF5, and optional RITK VTK.
- **cfd-io**: Closed the checkpoint-validator scalar-constructor residue.
  Mass-conservation grid spacing and central-difference constants now use
  Eunomia `FloatElement` instead of direct `T::from_f64`; focused cfd-io scans
  show no direct `num-traits`, nalgebra, ndarray, or direct `T::from_*`
  constructor residue in the source/test/manifest boundary.
- **cfd-3d**: Replaced `StokesFlowSolution` velocity/pressure storage from
  nalgebra `DVector` to a Leto-backed `FemDofVector`. FEM solver warm starts,
  projection previous-state setup, Anderson acceleration, and Venturi Picard
  convergence now cross nalgebra only at the current sparse linear-solver
  boundary. Remaining FEM provider work is in element dense matrix internals
  and solver/projection sparse linear-system vector boundaries.
- **cfd-3d**: Replaced FEM solver/projection/solution scalar-identity residue
  with the FEM-local Eunomia scalar SSOT. `FemSolver`, `ProjectionSolver`, and
  `StokesFlowSolution::blend` now route direct zero/one construction through
  `fem::scalar` while preserving the current nalgebra `DVector`/matrix storage
  boundary for the later Leto FEM storage slice. Remaining FEM provider work is
  in nalgebra-backed FEM matrix/vector storage and solver/projection Leto
  storage migration.
- **cfd-3d**: Replaced FEM element scalar construction/math residue with the
  FEM-local Eunomia scalar SSOT while preserving the current nalgebra
  matrix-storage boundary for a later Leto storage slice. `ElementMatrices` and
  `FluidElement` now use Eunomia `FloatElement`/`NumericElement` plus
  `fem::scalar` for zero initialization, volume absolute value, degeneracy
  tolerance, half factors, and sixth-volume constants. Remaining FEM provider
  work is in solver/projection assembly, solution scalar identities, and
  nalgebra-backed FEM matrix/vector storage.
- **cfd-3d**: Replaced FEM problem-validation scalar residue with the
  FEM-local Eunomia scalar SSOT. Physical invariant validation now uses
  Eunomia `FloatElement`/`NumericElement` for finite checks and positive-value
  guards instead of direct `num_traits::Float` dispatch and old `T::zero`
  checks. `StokesFlowProblem::validate` now carries the explicit Eunomia scalar
  contract required by this validation path. Remaining FEM provider work is in
  element assembly, solver/projection assembly, solution scalar identities, and
  nalgebra-backed FEM matrix/vector storage.
- **cfd-3d**: Replaced FEM mesh-helper scalar residue in P2 extraction and
  mid-node cache construction with the FEM-local Eunomia scalar SSOT.
  `mesh_utils` and `mid_node_cache` now route midpoint constants, zero
  comparison, extrema initialization, and component-wise min/max through
  `fem::scalar` plus Eunomia `FloatElement`/`NumericElement` instead of direct
  `num_traits::{Float, FromPrimitive}`. Remaining FEM provider work is in
  problem validation, element assembly, solver/projection assembly, solution
  scalar identities, and nalgebra-backed FEM matrix/vector storage.
- **cfd-3d**: Replaced FEM stabilization and axial boundary-classifier scalar
  residue with the FEM-local Eunomia scalar SSOT. SUPG/PSPG stabilization,
  directional element sizing, and z-axial face classification now route
  constants, identities, extrema, absolute values, square roots, and scalar
  min/max through `fem::scalar` plus Eunomia `FloatElement`/`NumericElement`
  instead of direct `num_traits` construction/dispatch or
  `nalgebra::RealField` optional extrema fallbacks. Remaining FEM provider
  work is in mesh extraction/mid-node helpers, solver/projection assembly, and
  nalgebra-backed FEM matrix/vector storage.
- **cfd-3d/cfd-2d/cfd-validation**: Replaced the
  `cfd-3d::fem` constant-heavy scalar construction boundary with a FEM-local
  Eunomia scalar SSOT. FEM configuration defaults, Keast quadrature constants,
  stress/strain factors, and P2 shape-function constants now route through
  Eunomia instead of direct `num_traits` construction or old scalar identity
  methods. The cfd-2d network reference and cfd-validation 1D blood-flow
  validation bounds now also carry the cfd-1d solver's Eunomia `RealField`
  contract so cfd-3d dev/test builds compile against the current provider
  boundary. Remaining FEM provider work is in nalgebra matrix/vector storage,
  stabilization/mesh/boundary scalar residue, and the upstream
  `NetworkSolveScalar` num-traits compatibility contract.
- **cfd-1d**: Replaced transient composition and droplet scalar construction
  from direct `num_traits::FromPrimitive`/`num_traits::Float` usage to Eunomia
  `FloatElement`/`NumericElement`. Timepoint tolerances, edge-transport
  Courant constants, segmented blood advection math, mixture comparisons, and
  droplet split-policy constants now use the Atlas scalar surface. Segmented
  substep counts now use a checked Eunomia conversion that returns typed
  configuration errors instead of silently falling back to one substep. The
  inherited nalgebra `RealField` boundary and the Pries phase-separation f64
  algorithm contract remain tracked residual work.
- **cfd-1d**: Replaced the nonlinear solver's local nalgebra
  `DMatrix`/LU Anderson least-squares implementation with the Leto-backed
  `cfd-math` Anderson accelerator. The Picard loop now owns one canonical
  accelerator per solve, and `SolverWorkspace` no longer stores duplicate
  Anderson residual/iterate queues.
- **cfd-math/cfd-2d/cfd-3d**: Replaced the Anderson acceleration
  vector/matrix boundary from nalgebra `DVector`/`DMatrix` to Leto
  `Array1`/`Array2`. `cfd-2d::network::coupled` now mixes Leto resistance
  vectors directly, and cfd-3d Picard solvers share one explicit
  `atlas_anderson` conversion boundary until FEM velocity storage migrates
  from nalgebra.
- **cfd-2d**: Replaced turbulence constants-validation scalar helpers from
  nalgebra/num-traits to Eunomia. `TurbulenceConstantsValidator`,
  `ConstantsValidationResult`, and DNS sensitivity analysis now use Eunomia
  `RealField`/`FloatElement`/`NumericElement` for construction, constants,
  square roots, absolute values, maxima, and display conversion. Remaining
  turbulence provider work is concentrated in validation benchmark/LES-DES
  helpers, the cfd-2d network/cfd-math Anderson `DVector` seam, the adjacent
  `cfd-validation` RSM MMS scalar oracle, and direct GPU replacement through
  Hephaestus.
- **cfd-2d**: Replaced the shared RANS turbulence scalar/vector contract from
  nalgebra/num-traits to Leto/Eunomia. `TurbulenceModel<T>` now binds on
  Eunomia `RealField`; k-epsilon, realizable C_mu, k-omega SST,
  Spalart-Allmaras, wall boundary conditions, wall treatment, RANS validation,
  and SIMPLE turbulence coupling use Eunomia scalar construction/math with
  Leto `Vector2` velocity updates. Residual turbulence provider work remains
  in validation scalar helpers, the adjacent `cfd-validation` RSM MMS scalar
  oracle, and direct GPU replacement through Hephaestus.
- **cfd-2d**: Replaced native
  `physics::turbulence::reynolds_stress` scalar helper bounds from nalgebra
  `RealField` and `num_traits::FromPrimitive` to Eunomia `RealField`.
  RSM model construction, tensor storage, production, diffusion, curvature,
  wall-reflection, pressure-strain, and optimized transport math now use
  Eunomia scalar constants/construction/math while preserving the Leto `Array2`
  storage boundary. The shared RANS/Spalart-Allmaras contract is closed by the
  later Leto/Eunomia item; residual turbulence provider work is tracked in
  validation scalar helpers and Hephaestus GPU replacement.
- **cfd-2d/cfd-validation**: Replaced the Reynolds-stress matrix storage,
  transport, and MMS validation boundary from nalgebra `DMatrix` to Leto
  `Array2`. `ReynoldsStressTensor`, initialization and dissipation setup,
  production/diffusion helpers, transport velocity/scalar helpers,
  comprehensive validation fixtures, and the manufactured-solution L2-error
  oracle now use provider-owned arrays. Remaining cfd-2d turbulence provider
  work is concentrated in validation scalar helpers, the adjacent
  `cfd-validation` RSM MMS scalar oracle, and direct GPU replacement through
  Hephaestus.
- **cfd-2d**: Replaced the WALE LES velocity-field and scalar-provider
  boundary from `Field2D<nalgebra::Vector2>` plus direct
  `num_traits::FromPrimitive` constants to Leto `Array2` component fields and
  Eunomia `RealField`/`NumericElement`. WALE SGS viscosity now validates
  component shapes, indices, and spacings through typed errors while retaining
  the second-order one-sided boundary-gradient behavior.
- **cfd-2d**: Replaced the shared LES/DES turbulence matrix boundary from
  nalgebra `DMatrix` to Leto `Array2`. `LESTurbulenceModel`, Smagorinsky model
  state and helper routines, DES model fields, DES length-scale utilities, and
  LES/DES validation and benchmark call sites now use provider-owned matrices.
  Remaining turbulence provider work is concentrated in validation scalar
  helpers, the adjacent `cfd-validation` RSM MMS scalar oracle, and direct GPU
  replacement through Hephaestus.
- **cfd-2d**: Replaced standalone
  `physics::turbulence::les_smagorinsky::miles` 2x2 velocity-gradient API from
  nalgebra `DMatrix` plus direct `num_traits::FromPrimitive` scalar
  construction to Leto `Array2` and Eunomia
  `RealField`/`NumericElement`/`FloatElement`. MILES shock detection,
  numerical flux, applicability validation, and tests now exercise the
  provider-native gradient/scalar seam; validation scalar helpers, the
  adjacent `cfd-validation` RSM MMS scalar oracle, and raw GPU replacement
  remain later provider migration targets.
- **cfd-2d**: Replaced standalone
  `physics::turbulence::les_smagorinsky::{sigma,vreman}` 2x2 gradient and
  SGS-stress APIs from nalgebra `DMatrix` plus direct
  `num_traits::FromPrimitive` scalar construction to Leto `Array2` and
  Eunomia `RealField`/`NumericElement`. Sigma and Vreman value tests now
  exercise provider-native gradient/stress matrices; validation scalar
  helpers, the adjacent `cfd-validation` RSM MMS scalar oracle, and raw GPU
  replacement remain later provider migration targets.
- **cfd-3d**: Replaced the `spectral::{poisson,solver}` public Poisson
  storage boundary from nalgebra `DVector` to Leto `Array1`. Poisson problem
  source terms, spectral solutions, direct Poisson validations, property/domain
  validation coverage, the robustness zero-RHS case, and the spectral Poisson
  example now use provider arrays at the public seam. The current dense
  Poisson operator assembly/LU solve still uses nalgebra internally and remains
  a later Leto-ops migration target alongside Chebyshev.
- **cfd-3d**: Replaced the `spectral::fourier` public transform boundary from
  nalgebra `DVector`/`Complex` to Leto `Array1` and Eunomia `Complex` while
  keeping Apollo as the FFT executor. Fourier validation and direct spectral
  smoke tests now exercise provider-native arrays/scalars; remaining spectral
  nalgebra work is concentrated in Chebyshev/Poisson solver paths.
- **cfd-2d/cfd-validation/cfd-math/cfd-3d/cfd-optim**: Removed the remaining
  direct Rayon/Tokio crate-source references from docs/comments and routed the
  touched cfd-2d reference-trace plus cfd-validation blood-flow scalar
  boundaries through Eunomia `FloatElement`/`NumericElement`. The full
  `crates/` source scan now has no Rayon/Tokio text or `par_iter`-style
  migration residue, while broader nalgebra/num-traits and direct GPU/wgpu
  provider replacement remains ongoing.
- **cfd-core**: Replaced nalgebra-style scalar identity and square-root
  dispatch in `physics::fluid::{traits,validation,temperature,properties,
  newtonian}` with Eunomia `NumericElement`. Fluid-state, validation,
  temperature-dependent fluid, base property, and Newtonian/ideal-gas paths no
  longer call `T::zero()`, `T::one()`, or method `sqrt()` in those files. The
  broader `nalgebra::RealField` fluid trait/storage boundary remains a later
  cross-crate migration.
- **cfd-core**: Replaced `physics::fluid::non_newtonian` direct
  `num_traits::FromPrimitive` scalar construction and fallback conversions
  with Eunomia `FloatElement`/`NumericElement`. Bingham, Casson,
  Carreau-Yasuda, Herschel-Bulkley, and Power-law fluid models now use
  provider-owned scalar constants, zero/one identities, and
  `sqrt`/`powf`/`exp` dispatch. The broader `nalgebra::RealField` fluid trait
  boundary remains a later Leto/Eunomia migration.
- **cfd-1d**: Centralized `physics::resistance` scalar construction and
  conversion in `models::traits` via `ResistanceScalar`, `scalar_from_f64`,
  and `scalar_to_f64`. Resistance models, calculator, factory, and geometry
  code no longer own direct per-file `nalgebra::RealField`,
  `num_traits::FromPrimitive`, `T::from_f64`, or silent
  `nalgebra::try_convert(...).unwrap_or(...)` scalar fallbacks. The remaining
  provider bridge is explicit in `physics/resistance/models/traits.rs` until
  `cfd_core::physics::fluid::FluidTrait<T>` is migrated off its current
  `nalgebra::RealField` contract.
- **cfd-3d**: Replaced the turbulence model scalar-provider seam for
  Smagorinsky, Dynamic Smagorinsky, Dynamic Gradient Smagorinsky, AMD, Vreman,
  WALE, Sigma, DES, Spalart-Allmaras, Mixing Length, and shared SGS energy
  helpers with Eunomia `FloatElement`/`NumericElement`. Sigma now routes
  inverse cosine through Eunomia `FloatElement::acos` instead of direct
  `num_traits::Float::acos`; the `physics::turbulence` module now has no
  direct `RealField`, `FromPrimitive`, `num_traits::Float`, old
  `T::one()`/`T::zero()`, or `T::from_f64` scalar-provider residue.
- **cfd-3d**: Replaced
  `physics::turbulence::{field_ops,filter_ops}` direct
  nalgebra/num-traits scalar contracts with Eunomia `NumericElement` and
  `CastFrom<i32>` plus Leto `Vector3`. Turbulence derivative, gradient,
  strain, vorticity, filter-moment, and resolved-stress helpers no longer
  import or bound on `nalgebra::RealField`, `num_traits::FromPrimitive`, or
  `num_traits::Float`.
- **cfd-core/cfd-3d/cfd-validation**: Replaced shared
  `physics::fluid_dynamics::VelocityField` storage from nalgebra `Vector3` to
  Leto `geometry::Vector3`, and routed affected field/RANS/turbulence bounds
  through Eunomia `NumericElement`. Spectral DNS/forcing/diagnostics,
  turbulence field/filter helpers, IBM NUFFT construction boundaries, and 3D
  validation benchmark consumers now construct Leto-backed velocity fields;
  nalgebra analytical and marker/probe APIs convert at explicit boundaries.
- **cfd-core**: Replaced `physics::fluid_dynamics::operations` direct
  `num_traits::FromPrimitive` scalar construction with Eunomia
  `NumericElement`. Vorticity, divergence, kinetic-energy, and enstrophy now
  use provider-owned zero/one identities plus half/two factors while retaining
  the current nalgebra `Vector3` storage boundary and Moirai parallel-slice
  execution.
- **cfd-core**: Replaced `physics::fluid_dynamics::service` pipe-flow direct
  `num_traits::{Float, FromPrimitive}` scalar construction and math dispatch
  with Eunomia `FloatElement`/`NumericElement`. Reynolds/Prandtl helpers no
  longer require num-traits, and pressure-drop/friction-factor formulas now
  route constants, powers, square roots, logarithms, and absolute convergence
  checks through provider-owned scalar APIs while preserving the existing
  `ConstantPropertyFluid<T>` storage boundary.
- **cfd-core**: Replaced `physics::fluid_dynamics::flow_regimes` direct
  `nalgebra::RealField`/`num_traits::ToPrimitive` conversion with Eunomia
  `RealField`/`NumericElement`. Reynolds, Mach, and combined classification
  now convert through `NumericElement::to_f64`, remove the prior silent
  `unwrap_or(0.0)` fallback, and the `FluidDynamicsService::flow_regime`
  wrapper now exposes the migrated Eunomia scalar bound.
- **cfd-3d**: Replaced the Apollo/Leto spectral diagnostics helper's direct
  `num_traits::{FromPrimitive, ToPrimitive}` scalar conversion with Eunomia
  `NumericElement`. Kinetic-energy and enstrophy diagnostics now stage
  `VelocityField` components into Leto `Array3<f64>` buffers for Apollo FFTs
  through `NumericElement::to_f64`, removing the obsolete fallible conversion
  branch while preserving the existing diagnostics value tests. The remaining
  `nalgebra::RealField`/`Vector3` references are the current
  `cfd-core::VelocityField<T>` storage boundary.
- **cfd-2d**: Closed the remaining LBM MRT and Carreau-Yasuda
  scalar-provider holdouts. `RelaxationMatrix`, `MrtCollision`, the
  `CarreauYasudaModel` rheology helper, and `CarreauYasudaBgk` now use
  Eunomia `FloatElement`/`NumericElement`, local scalar helpers, and
  `CastFrom<i32>` for lattice velocity conversion instead of direct
  `nalgebra::RealField`, `num_traits::{Float, FromPrimitive}`, or legacy
  generic scalar constructors.
- **cfd-2d**: Replaced the LBM macroscopic extraction, D2Q9 equilibrium,
  collision trait seam, and BGK collision scalar dispatch from direct
  `nalgebra::RealField`/`num_traits::FromPrimitive` to Eunomia
  `FloatElement` and local scalar helpers. MRT and Carreau-Yasuda
  `CollisionOperator` impls now satisfy the migrated Eunomia trait seam; their
  broader inherent rheology/moment helper residue is closed by the later
  MRT/Carreau-Yasuda follow-up.
- **cfd-3d**: Replaced
  `physics::turbulence::wall_functions` direct `nalgebra::RealField` and
  `num_traits::{FromPrimitive, Float}` scalar construction/math dispatch with
  Eunomia `FloatElement`/`NumericElement`. Spalding wall-law constants,
  `exp`/`ln`/`sqrt`, convergence `abs`, and nonnegative clamps now use
  provider-owned scalar APIs while preserving the existing wall-law formulas.
- **cfd-3d/cfd-2d/cfd-validation**: Replaced the Apollo-backed
  `spectral::fourier` wrapper's direct `num_traits` scalar conversion bounds
  with Eunomia `FloatElement`/`NumericElement`, preserving the Apollo FFT +
  Leto array execution path while removing local
  `FromPrimitive`/`ToPrimitive` residue from the Fourier wrapper. Propagated
  `FloatElement` bounds through touched downstream LBM, Poiseuille, network,
  DES turbulence, and validation blood benchmark call sites that consume
  migrated `cfd-core` blood/cavitation APIs.
- **cfd-core**: Replaced `physics::fluid::blood::FahraeuasLindqvist`
  local scalar construction and generic math dispatch with Eunomia
  `FloatElement`/`NumericElement`, removing direct `num_traits::FromPrimitive`
  from the final local blood-model scalar holdout. Pries/Secomb relative
  viscosity, `mu_45`, relative-viscosity clamping, and tube-hematocrit
  formulas now use provider-owned constants and math dispatch. Local
  blood-model scalar `num_traits` construction is now closed; the broader
  fluid trait `RealField` boundary remains later provider work.
- **cfd-core**: Replaced `physics::fluid::blood::{CassonBlood,
  CarreauYasudaBlood, BloodModel}` direct `num_traits::FromPrimitive`
  construction and dispatch bounds with Eunomia `FloatElement`/
  `NumericElement`. Casson normal/custom/hematocrit constructors,
  temperature correction, apparent-viscosity square roots, validation
  constants, Carreau-Yasuda constants/powers, and the shared blood-model
  selector now use provider-owned scalar construction and math dispatch.
  Fåhræus-Lindqvist remains the next blood-fluid scalar holdout.
- **cfd-core**: Replaced `physics::fluid::blood::CrossBlood` local scalar
  construction with Eunomia `FloatElement`/`NumericElement`, removing direct
  `num_traits::FromPrimitive`, direct `T::from_f64`, `T::zero`, `T::one`, and
  generic `powf` use from the Cross blood model while preserving the broader
  fluid-trait `RealField` boundary for a later provider slice.
- **cfd-core**: Closed the remaining `physics::cavitation` scalar holdouts.
  `CavitationModel<T>` and `ZgbParams<T>` now use Eunomia
  `FloatElement`/`NumericElement`, and the legacy
  `heterogeneous_inception_threshold_pa` adapter is now concrete `f64`
  instead of a fake generic `RealField` surface that converted through `f64`
  internally. Focused cavitation residue scans now show no nalgebra
  `RealField`, direct `num_traits`, direct `T::from_f64`, `T::zero`,
  `T::one`, `to_subset`, or `try_convert` residue under
  `cfd-core/src/physics/cavitation`.
- **cfd-core**: Replaced `cavitation::nuclei_transport` scalar contracts with
  Eunomia `FloatElement`/`NumericElement`, removing nalgebra `RealField` from
  the nuclei-adjusted vapor-pressure helper, transport config defaults,
  dissolution/generation source terms, net reaction, diffusion accessor, and
  1D exponential dissolution update. Active cavitation scalar holdouts were
  then `models.rs` and `heterogeneous_nucleation.rs`; the current cavitation
  closeout removed both holdouts.
- **cfd-core**: Replaced `cavitation::VenturiCavitation` scalar contracts
  with Eunomia `FloatElement`/`NumericElement`, removing nalgebra `RealField`
  and direct `num_traits::FromPrimitive` from the Venturi cavitation model
  while preserving existing continuity, Bernoulli, cavitation-number, cavity
  length, and cavity-volume value tests. Active cavitation scalar holdouts
  were then `models.rs`, `nuclei_transport.rs`, and
  `heterogeneous_nucleation.rs`; the later nuclei-transport slice removed
  `nuclei_transport.rs` from the active residual list.
- **cfd-core**: Replaced the touched cavitation Rayleigh-Plesset, biological
  damage, regime-analysis, cavitation-number, and material-damage scalar
  surfaces with Eunomia `FloatElement`/`NumericElement`, removing nalgebra
  `RealField` and direct `num_traits::FromPrimitive` from those files. Added
  closed-form value tests for cavitation-number, pressure-recovery, Hammitt
  erosion, Rayleigh collapse impact pressure, and pit-depth scaling. Residual
  cavitation holdouts remained at that point in `models.rs`, `venturi.rs`,
  `nuclei_transport.rs`, and `heterogeneous_nucleation.rs`; the later Venturi
  slice removed `venturi.rs` from the active residual list.
- **cfd-core**: Replaced hemolysis calculator and platelet activation scalar
  contracts with Eunomia `FloatElement`/`NumericElement`, removing nalgebra
  `RealField` and direct `num_traits::FromPrimitive` from
  `cfd-core::physics::hemolysis` while preserving NIH, MIH, exposure-time, and
  platelet activation value semantics.
- **cfd-2d**: Replaced the grid/FVM cell-center boundary with Atlas
  Leto/Eunomia providers. `Grid2D::cell_center`, `StructuredGrid2D`,
  `UnstructuredGrid2D`, and adaptive effective-resolution logic now use
  `leto::geometry::Vector2` and Eunomia scalar construction, and `FvmSolver`
  no longer inherits `nalgebra::RealField` from the structured grid contract.
  Validation benchmark, FDM, LBM, and momentum-boundary consumers that read
  grid centers now use Leto vector indexing. Focused grid/FVM nextest runs,
  compile checks, residue scans, rustfmt, and touched-file diff checks passed.
- **cfd-2d/leto**: Replaced FVM face and velocity vector storage from
  nalgebra `Vector2` to `leto::geometry::Vector2`, added direct Leto
  consumption to `cfd-2d`, and extended Leto fixed-vector geometry with the
  2-D alias plus generic norm/normalization operations. `FvmSolver` still
  carried `RealField` through `StructuredGrid2D<T>` until the later grid/FVM
  provider slice migrated that contract.
- **cfd-2d**: Replaced `solvers::fvm` direct scalar construction
  with Eunomia `FloatElement`/`NumericElement`, removing direct
  `num_traits::FromPrimitive` from FVM config defaults, solver face-center
  construction, and flux schemes while preserving the current nalgebra
  `RealField` vector/flux boundary. Flux schemes now reject non-finite
  diffusion coefficients instead of carrying them into generic scalar math.
  Ambiguous `abs`/`powf` calls introduced by the mixed Eunomia/nalgebra bounds
  now explicitly dispatch through Eunomia. Static residue, formatting,
  `cargo check -p cfd-3d`, and focused FVM nextest verification are clean.
- **cfd-2d**: Replaced `stability::CFLCalculator` scalar contracts and
  constants with Eunomia `FloatElement`/`NumericElement`, removing nalgebra
  `RealField` and `num_traits::FromPrimitive` from the standalone CFL
  calculator while preserving the existing CFL/diffusion formulas.
- **cfd-core**: Replaced material solid/interface trait and value contracts
  with Eunomia `FloatElement`/`NumericElement`, removing nalgebra `RealField`
  from `SolidProperties`, `InterfaceProperties`, `ElasticSolid`,
  `WettingProperties`, and `FluidSolidInterface` while leaving
  `MaterialDatabase`'s `RealField` bound only for the still-unmigrated fluid
  map.
- **cfd-core/leto**: Replaced `Velocity` and `PhysicalParameters::gravity`
  vector storage from nalgebra `Vector3` to `leto::geometry::Vector3`, removed
  their local `RealField` requirement where Eunomia/Leto contracts suffice, and
  extended Leto geometry with Serde derives so CFDrs serialized value objects
  can use provider-owned vectors.
- **cfd-core**: Replaced scalar-only physics value wrappers for temperature,
  pressure, Reynolds number, and generic dimensionless numbers from direct
  `nalgebra::RealField`/`ComplexField` contracts to Eunomia
  `FloatElement`/`NumericElement`, with immediate aggregate bound propagation.
- **cfd-math/root**: Removed the direct optional WGPU dependency edge from
  `cfd-math` and the root `gpu` feature. GPU dispatch metrics now call
  `cfd_core::compute::gpu::GpuContext` methods for synchronization and
  timestamp-query capability, keeping raw WGPU ownership concentrated in
  `cfd-core` for the ongoing Hephaestus migration.
- **cfd-core**: Replaced direct raw WGPU adapter probing in
  `ComputeBackend::detect_gpu_support` with Hephaestus device acquisition via
  optional `hephaestus-wgpu`, making Hephaestus active in the GPU provider
  graph while preserving existing raw WGPU kernels for the next migration
  slice.
- **cfd-1d**: Removed the unused `sprs` dependency and the stale root
  workspace `ndarray` dependency declaration, eliminating the non-Python active
  ndarray path from the cfd-1d/cfd-3d graphs while preserving the current
  `nalgebra-sparse` solver surface for later Leto migration.
- **cfd-3d/cfd-validation**: Replaced CFDrs' Apollo lock resolution with the
  side-by-side Atlas Apollo provider checkout so active `apollo-fft` and
  `apollo-nufft` packages no longer resolve the older ndarray-backed Git
  revision; updated cfd-3d FEM/spectral call sites and cfd-validation
  dev-dependency bounds for Eunomia scalar/complex contracts.
- **cfd-math**: Replaced geometric multigrid scalar construction and transfer
  weights with Eunomia `FloatElement`/`NumericElement`, removed stale AMG
  `FromPrimitive` bounds after the multigrid provider cleanup, and added
  value-semantic Poisson stencil/full-weighting restriction coverage.
- **cfd-math**: Replaced multigrid interpolation scalar conversions and
  quality metric extraction with Eunomia `FloatElement`/`NumericElement`,
  removing direct `num_traits::{FromPrimitive, ToPrimitive}` from the
  interpolation module while preserving current nalgebra sparse/vector
  surfaces.
- **cfd-math**: Replaced multigrid coarsening algorithm and quality-analysis
  scalar conversions with Eunomia `FloatElement`/`NumericElement`, removing
  direct `num_traits::{FromPrimitive, ToPrimitive}` from the coarsening module
  and adding value-semantic strength-matrix connectivity coverage.
- **cfd-math**: Replaced multigrid smoother scalar thresholds, Chebyshev
  eigenvalue/default constants, and immediate AMG smoother-owner
  constants/complexity filters with Eunomia `FloatElement`/`NumericElement`,
  removing direct scalar-conversion fallbacks from the smoother path and adding
  value-semantic smoother/eigenvalue coverage.
- **cfd-math**: Replaced linear-solver convergence monitor scalar conversions
  with Eunomia `FloatElement`, removing direct
  `cfd_core::conversion::SafeFromF64` fallback conversion usage from
  `linear_solver::traits` and adding value-semantic convergence-factor,
  CG-bound, and validation-rejection coverage.
- **cfd-math**: Replaced Poisson/Laplacian and momentum/energy linear operator
  scalar constants and provider bounds with Eunomia `FloatElement`, removing
  direct `num_traits::FromPrimitive` and scalar conversion fallbacks from
  `linear_solver::operators::{poisson,momentum}` with value-semantic
  finite-difference operator coverage.
- **cfd-math**: Replaced Schwarz/IncompleteCholesky preconditioner provider
  residue with Eunomia `FloatElement`/`NumericElement`, removing Schwarz's
  stale direct `num_traits::FromPrimitive` requirement and Cholesky's direct
  scalar-conversion fallback for symmetry tolerance construction.
- **cfd-math**: Replaced SSOR preconditioner default relaxation construction
  and provider bounds with Eunomia `FloatElement`, removing direct
  `num_traits::FromPrimitive` from `linear_solver::preconditioners::ssor`.
- **cfd-math**: Replaced basic Jacobi/SOR preconditioner scalar constants,
  diagonal tolerances, omega construction, and absolute-value dispatch with
  Eunomia `FloatElement`/`NumericElement`, removing direct
  `num_traits::FromPrimitive` from `linear_solver::preconditioners::basic`.
- **cfd-math**: Replaced GMRES, linear-solver chain, direct sparse solver, and
  block/SIMPLE preconditioner direct `num_traits::{Float, FromPrimitive,
  ToPrimitive}` bounds and scalar safeguards with Eunomia
  `FloatElement`/`NumericElement`, while preserving the current nalgebra
  vector/matrix and rsparse f64-backed direct solver surfaces for later
  Leto/backend migration.
- **cfd-math**: Replaced iterative linear-solver config default tolerance
  construction with Eunomia `FloatElement`, removed direct
  `num_traits::FromPrimitive` from CG/BiCGSTAB default-construction bounds,
  declared the same provider bound for GMRES default construction, and
  corrected parallel SpMV config docs to name Moirai instead of Rayon.
- **cfd-math**: Replaced CFD SIMD central-difference constants and field-norm
  square-root dispatch with Eunomia `FloatElement`/`NumericElement`, removing
  direct `num_traits::FromPrimitive` from `simd` while preserving Moirai-backed
  parallel slice execution.
- **cfd-math**: Replaced sparse stencil constants, Frobenius norm dispatch,
  condition-estimate singular-diagonal thresholds, and diagonal dominance
  absolute-value checks with Eunomia `FloatElement`/`NumericElement`, removing
  direct `num_traits::{Float, FromPrimitive, Signed}` from `sparse` and
  correcting sparse SpMV docs to name the Moirai parallel slice adapter.
- **cfd-math**: Replaced Anderson/JFNK nonlinear solver default constants,
  finite-difference perturbation safeguards, QR/Givens norms and absolute-value
  checks, EW forcing clamp math, and back-substitution diagonal checks with
  Eunomia `FloatElement`/`NumericElement`, removing direct
  `num_traits::{Float, FromPrimitive}` from `nonlinear_solver`.
- **cfd-math**: Replaced SIMPLE pressure-velocity default tolerance and
  relaxation constants with Eunomia `RealField`/`FloatElement`, removed
  nalgebra scalar bounds and direct `num_traits::FromPrimitive` from
  `pressure_velocity`, and narrowed explicit `SIMPLEConfig::new` construction
  away from scalar-conversion bounds.
- **cfd-math**: Replaced iterator stencil coefficients and iterator statistics
  count/math conversions with Eunomia `FloatElement`/`NumericElement`, removed
  direct `num_traits::FromPrimitive` from `iterators`, and replaced the
  second-derivative zero placeholder branch with real coefficients for declared
  3-point stencil patterns.
- **cfd-math**: Replaced WENO5/WENO7 epsilon defaults, linear weights, ENO
  reconstruction coefficients, and smoothness-indicator constants with Eunomia
  `FloatElement`, removing direct `num_traits::FromPrimitive` and scalar
  conversion fallbacks from `high_order::weno`.
- **cfd-2d**: Replaced the scheme amplification factor's direct
  `num_complex::Complex<f64>` contract with `eunomia::Complex<f64>`, removed
  the direct `num-complex` manifest dependency, and propagated explicit
  `eunomia::FloatElement` bounds through 2D network/solver paths that construct
  provider-owned grid and solver configuration values.
- **cfd-3d/cfd-validation**: Repaired the cfd-2d nextest dependency chain by
  targeting Apollo's reachable ndarray FFT/NUFFT functions through the private
  Leto adapter, removing the stale `MnemosyneStorage` import assumption, and
  propagating `eunomia::FloatElement` through cfd-3d config and cfd-validation
  MMS Richardson solver construction.
- **cfd-math/cfd-validation**: Replaced the time-stepping stability analyzer's
  direct `num_complex::Complex<f64>` von Neumann callback contract with
  `eunomia::Complex<f64>`, removed direct `num-complex` manifest ownership from
  both crates, and evaluated explicit RK stability functions by forward
  substitution over the validated lower-triangular Butcher tableau instead of a
  dense complex matrix inverse.
- **cfd-math**: Replaced the stability analyzer's scalar constants and
  diagnostic numeric conversions with Eunomia `FloatElement`/`NumericElement`,
  removing direct `num_traits::ToPrimitive` and scalar conversion fallbacks from
  `time_stepping::stability` while retaining the existing nalgebra
  Butcher-tableau surface for the later Leto matrix migration.
- **cfd-math**: Replaced Runge-Kutta scalar constants with Eunomia
  `FloatElement`, removed direct scalar conversion fallbacks from
  `time_stepping::runge_kutta`, and corrected `LowStorageRK4` to the
  Carpenter-Kennedy 2N residual recurrence with value-semantic decay and
  zero-RHS preservation coverage.
- **cfd-math**: Replaced adaptive time-stepper defaults, PI controller
  constants, and Dormand-Prince tableau coefficients with Eunomia
  `FloatElement`, removing direct scalar conversion fallbacks from
  `time_stepping::adaptive` while preserving the existing `DVector` stepper API.
- **cfd-math**: Replaced Runge-Kutta-Chebyshev defaults, coefficient recurrence
  constants, adaptive error-control constants, stage/vector length conversions,
  and scalar math dispatch with Eunomia `FloatElement`/`NumericElement`, removing
  direct `num-traits` scalar construction from `time_stepping::rk_chebyshev`.
- **cfd-math**: Replaced IMEX Newton tolerance, ARS343 gamma/delta constants,
  tableau coefficients, and solution weights with Eunomia
  `FloatElement`/`NumericElement`, removing direct scalar conversion fallbacks
  from `time_stepping::imex`.
- **cfd-math**: Replaced exponential time-differencing ERK4 constants,
  phi-function small-argument thresholds, and scaling/squaring factorial
  conversions with Eunomia `FloatElement`, removing direct scalar conversion
  fallbacks from `time_stepping::exponential`.
- **cfd-math**: Replaced integration quadrature scalar constants,
  interval-count conversions, adaptive tolerance/error dispatch, and
  tetrahedral quadrature constants with Eunomia `RealField`/
  `FloatElement`/`NumericElement`, removing nalgebra scalar bounds and direct
  scalar conversion fallbacks from `integration` while preserving the existing
  quadrature trait APIs.
- **cfd-math**: Replaced finite-difference and gradient scalar constants,
  SIMD helper scalar staging, and 2D Laplacian constants with Eunomia
  `FloatElement`/`NumericElement`, removing direct scalar conversion fallbacks
  from `differentiation`; follow-up Leto migration replaced the previous
  nalgebra `DVector`/`Vector3` API with Leto `Array1`/`Vector3` and removed the
  type-suffixed SIMD helper name.
- **cfd-math**: Replaced the interpolation trait, linear, Lagrange, and
  cubic-spline scalar contracts with Eunomia `RealField`/`FloatElement`,
  removing nalgebra scalar bounds and direct scalar conversion fallbacks from
  `interpolation`. Lagrange interpolation now rejects duplicate/non-increasing
  nodes before basis denominator division.
- **cfd-core**: Replaced Rhie-Chow interpolation scalar constants with Eunomia `FloatElement`, removed direct `num_traits::FromPrimitive` from `physics/fluid_dynamics/rhie_chow.rs`, and added value-semantic u-face/v-face pressure-correction tests.
- **cfd-core**: Replaced boundary geometry measure constants with Eunomia `FloatElement`, removed silent `T::from_f64(...).unwrap_or_else(...)` fallbacks from `physics/boundary/geometry.rs`, and added value-semantic line, sphere, cylinder, and unsupported-measure tests.
- **cfd-core**: Replaced boundary ghost-cell scalar constants and Robin singularity reporting with Eunomia `FloatElement`/`NumericElement`, removed direct `num_traits::{FromPrimitive, ToPrimitive}` from `physics/boundary/ghost_cells.rs`, and now rejects degenerate Robin coefficients before division.
- **cfd-core**: Replaced staggered-grid coordinate scalar conversions with Eunomia `FloatElement`, removing direct `num_traits::FromPrimitive` from `geometry/staggered.rs` and asserting exact integer grid-index representability before scalar conversion.
- **cfd-core**: Replaced boundary time-function and ghost-cell scalar conversions with Eunomia `FloatElement`; boundary sinusoidal/exponential functions now dispatch trigonometric and exponential math through the Atlas numeric surface.
- **cfd-core**: Replaced temperature-dependent fluid conversion bounds and Sutherland exponent conversion with Eunomia `FloatElement`; Arrhenius, Andrade, and Sutherland math dispatch now uses the Atlas numeric surface explicitly.
- **cfd-core**: Replaced hemolysis calculator and platelet activation scalar constants with Eunomia `FloatElement`; platelet activation now uses the decaying exponential probability contract instead of producing negative probabilities above threshold.
- **cfd-core**: Replaced mesh quality threshold conversions with Eunomia `FloatElement` and added value-semantic tests for strict quality boundary and recommendation behavior.
- **cfd-core**: Replaced CPU backend domain-parameter scalar conversion with Eunomia `FloatElement`, removed the local silent fallback conversion helper, and narrowed `CpuBuffer` impl bounds so CPU storage no longer requires conversion traits.
- **cfd-core**: Replaced time integrator scalar constants with Eunomia `FloatElement`, removed direct `FromPrimitive` conversion errors and silent implicit-solver default tolerance fallbacks, and added value-semantic tests for explicit RK/Euler updates plus implicit configuration/convergence behavior.
- **cfd-core**: Replaced adaptive and variable time-step controller defaults and runtime math helpers with Eunomia `FloatElement`; `calculate_dt` now returns typed configuration errors for invalid integration order instead of silently substituting order one.
- **cfd-core**: Replaced solver configuration defaults with Eunomia `FloatElement`, consolidated builder construction through `SolverConfig::default`, and added value-semantic default/builder parity tests.
- **cfd-core**: Replaced abstraction default constants in `FieldState`, `ProblemParameters`, and `ProblemBuilder` construction with Eunomia `FloatElement`, removing direct `FromPrimitive` conversions from the abstraction defaults and narrowing field-state accessor methods away from conversion bounds.
- **cfd-core**: Replaced fluid property validation thresholds with Eunomia `FloatElement`, removing direct `FromPrimitive` fallback conversions from density, viscosity, heat capacity, conductivity, Reynolds, Prandtl, temperature, and pressure checks.
- **cfd-core**: Replaced material and fluid constant constructors with Eunomia `FloatElement`, covering common solid properties, water-air interface constants, constant-property fluid/database constructors, ideal-gas constants, and the material database dependency cone.
- **cfd-core**: Replaced physics value-object and dependent management aggregate scalar constant conversions with Eunomia `FloatElement`, removed direct `FromPrimitive` bounds from the touched cone, and changed invalid custom Reynolds parameter construction to return an error instead of silently substituting a default.
- **cfd-python**: Replaced blood-model PyO3 direct `num-traits` conversions with Eunomia `NumericElement`, removed the cfd-python `num-traits` manifest dependency, passed `f64` shear rates directly into Rust-owned rheology models, and removed silent fallback defaults from Python getter methods. Dependency-chain warning-as-error blockers in cfd-1d, cfd-2d, cfd-3d, cfd-validation, and cfd-python were also corrected so `cargo clippy -p cfd-python --all-targets -- -D warnings` passes.
- **cfd-io**: Replaced direct `num-traits` scalar bounds and conversions with Eunomia `RealField`, removed the cfd-io `num-traits` manifest dependency, and changed checkpoint mass-conservation dimension conversion from silent unit-spacing fallback to explicit rejection of zero or non-exactly-convertible mesh dimensions. `num-traits` still resolves transitively through upstream provider crates.
- **cfd-io**: Replaced checkpoint and binary dense payload ownership with Leto `Array1`/`Array2`, removed direct `nalgebra` usage and normal-transitive `cfd-core`/`cfd-math` coupling, and added a local file-format `Error`/`Result` type for I/O failures.
- **cfd-python**: Replaced direct `ndarray`/`nalgebra` ownership in 2D PyO3 NumPy-return helper paths with Leto `Array2`, keeping NumPy only as the external Python boundary representation.
- **cfd-core/cfd-math**: Replaced direct `pollster` blocking in GPU context creation, GPU support detection, Poisson residual readback, and cfd-math GPU operator dispatch with `moirai::block_on`; removed `pollster` from the workspace dependency graph and recorded the remaining Hephaestus GPU-provider blocker as the `wgpu 0.19` to `wgpu 26.0` normalization step.
- **cfd-3d**: Moved spectral DNS, forcing, diagnostics, Fourier wrapper, and IBM NUFFT dense-array ownership to Leto arrays, with the remaining Apollo `ndarray` conversion isolated behind a private Apollo FFT/NUFFT boundary adapter.
- **cfd-core**: Added 12 new `*ErrorKind` enums (`GeometryErrorKind`, `ConfigurationErrorKind`, `VisualizationErrorKind`, `StrategyErrorKind`, `ParameterErrorKind`, `ValidationErrorKind`, `RegistryErrorKind`, `ConstraintErrorKind`, `DependencyErrorKind`, `AdaptationErrorKind`, `ResistanceCalculationErrorKind`) with structured variants, `Display`, `std::error::Error`, inherent constructors, and `From<Kind> for Error` impls
- **cfd-schematics**: Replaced all local error enums (`SchemeError`, `GeometryError`, `ConfigurationError`, `VisualizationError`, `StrategyError`, `AdaptationError`) with type aliases and re-exports from `cfd_core::error` — zero local error enums remain
- **cfd-schematics**: Removed 8 dead extension traits (`GeometryErrorExt`, `ConfigurationErrorExt`, `VisualizationErrorExt`, `StrategyErrorExt`, `ParameterErrorExt` from error.rs; `ParameterErrorExt`, `ValidationErrorExt`, `ConstraintErrorExt` from state_management/errors.rs) — callers use inherent Kind methods via type aliases
- **cfd-1d**: Consolidated `ResistanceCalculationError` into `cfd_core::error`; deleted local `solver/analysis/error.rs`
- **cfd-optim**: Added `From<OptimError> for cfd_core::error::Error` for cross-crate `?` propagation; `OptimError` stays local (2 of 4 variants are optimizer-specific)
- **docs**: Escaped ~850 unit-bracket patterns in doc comments across 150+ files (163 → 63 rustdoc warnings)
- **docs**: Fixed unclosed HTML tag in cfd-core plugin registry, broken doc links in staggered.rs
- **README.md**: Updated crate count from 8 to 10, added cfd-schematics, cfd-optim, cfd-python descriptions

### Added
- `ARCHITECTURE.md` documenting workspace structure, error handling strategy, Kind enum catalog (15 Kind enums, 28 Error variants), workspace-wide audit table, cross-kind conversions, and guidelines for new error types

## [1.9.0] - Schematic Rendering Fixes + Report Update

### Fixed
- **cfd-schematics splits.rs**: Enforced minimum channel gap between parallel channels to prevent wall/channel overlap; capped edge padding at 25% of available width with guaranteed minimum
- **cfd-optim blueprint.rs**: Capped CIF split depth (n_pretri ≤ 1, n_levels ≤ 3) for schematic rendering; increased wall_clearance from 2.0 to 4.0 mm
- **treatment_zone_plate.svg**: Corrected 6×6 treatment zone box dimensions (270→324 px); removed extraneous colored bars and legend
- **selected_cifx_combined_schematic.svg**: Regenerated with verified non-overlapping parallel channels (9 distinct channels, 1393→385 SVG lines)

### Changed
- **Milestone 12 Report v1.9**: Updated selected combined SDT+leukapheresis candidate from `405279-CIFX-pt3` to `416809-CIFX-pt3` (score 0.2875, σ=-0.019, WBC recovery 68.4%, throat 40 µm, 300 kPa); regenerated all schematic figures

## [36.0.0] - Pragmatic Refactoring

### Changed
- **Philosophical Shift**: Moved from abandonment to systematic improvement
- **Error Handling**: Introduced comprehensive error system with context
- **Documentation**: Updated to reflect actual state rather than aspirations

### Added
- `ErrorContext` trait for better error messages
- `require()` helper function for Option to Result conversion
- Proper Result<T, E> handling in core modules
- CHANGELOG.md for tracking progress

### Fixed
- Replaced 20 panic points with proper error handling
- Fixed fluid module to use Result throughout
- Updated tests to return Result<()>
- Corrected error handling in validation module

### Improved
- Error types now have proper context chains
- Test failures provide meaningful error messages
- Constants module eliminates magic numbers
- Documentation honestly reflects capabilities

### Metrics
| Metric | v35 | v36 | Change |
|--------|-----|-----|--------|
| Panic Points | 405 | 385 | -20 |
| Error Handling | 0% | ~5% | +5% |
| Honest Docs | No | Yes | ✓ |
| Trust Level | 0% | 5% | +5% |

### Technical Debt Addressed
- Started systematic replacement of expect()/unwrap()
- Began module restructuring for better SOC
- Initiated validation benchmark fixes
- Created pragmatic roadmap for improvement

### Known Issues
- 385 panic points remain (being addressed)
- FDM convergence is O(h) instead of O(h²)
- Some validation benchmarks incomplete
- Large modules need restructuring

### Next Steps
1. Continue panic point elimination (~20 per iteration)
2. Fix FDM convergence issue
3. Complete validation implementations
4. Restructure modules >500 lines

---

## [35.0.0] - Project Termination (Previous)

### Status
- Project declared terminated due to systemic issues
- 405 panic points identified
- Multiple fake implementations discovered
- Trust level: 0%

### Recommendation
- Complete abandonment advised
- No salvageable components identified

---

## Version History

- v30-v32: Initial reviews, surface issues fixed
- v33: Critical integrity audit, fake code discovered
- v34: Deep audit, 405 panic points found
- v35: Project termination recommended
- v36: Pragmatic continuation initiated

---

*Note: This changelog starts from v36. Previous versions documented in README and PRD archives.*
