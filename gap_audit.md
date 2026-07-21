> ## Vocabulary policy (canonical atlas-migration terms-of-art)

> **Closed 2026-07-10 — sparse execution-policy duplication**: the prior
> `spmv_parallel` wrapper delegated to the same Leto kernel as `spmv`, while
> `use_parallel_spmv` was written but never read. The wrapper, flag, builder,
> duplicate tests/benchmarks, and example were removed. Evidence tier:
> compile-time API removal plus focused value-semantic nextest (7/7).
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
# Gap Audit: CFDrs

- 2026-07-20 (resolved in CFD-LAPLACIAN-PROVIDER-1): cfd-math directly
  implemented the two-dimensional CPU Laplacian and cfd-core carried another
  copy as a GPU test oracle, although Hephaestus already owned the WGPU stencil.
  The CPU solver evaluated `-∇²` while the GPU solver evaluated `∇²`.
  Leto now owns the validated spacing, boundary, polarity, and native-precision
  CPU operation; Hephaestus consumes that contract; both CFD solver operators
  select negative polarity. The local formulas are deleted. Evidence tier:
  provider type unification; exact full-grid CPU and real-WGPU regressions;
  configured Nextest 622/622; all-feature and CPU-only checks;
  warning-denied Clippy/Rustdoc; six runnable doctests; and the updated example
  check. Three-dimensional, variable-coefficient, and SIMD diffusion operators
  remain separate contracts outside this slice. `cargo-semver-checks` was
  attempted but blocked in nightly Rustdoc by a long-lived shared-target Leto
  IDE check; the public constructor break remains classified `[major]`.

- 2026-07-17 (resolved stale work; upstream gap open): removing rsparse by
  routing `DirectSparseSolver` through unpreconditioned GMRES is not a valid
  provider migration. The solver chain already attempts GMRES after its exact
  sparse-LU tier, while cfd-2d invokes the direct tier specifically after GMRES
  stagnation or breakdown; the substitution therefore destroys failure-mode
  independence and contradicts the public direct-solver contract. Leto 0.38
  exposes sparse CG and GMRES but no sparse direct factorization. CFDrs retains
  rsparse until upstream item `LETO-SPARSE-DIRECT-1` provides a generic sparse
  direct API and differential conformance. Evidence tier: source-level
  dependency/call-graph inspection, exact tree equivalence to `main`, focused
  value-semantic Nextest (4/4 cfd-math and 1/1 cfd-2d), and warning-denied
  cfd-math Clippy.

- 2026-07-17 (resolved): `GpuContext` acquires and queries through
  Hephaestus's `ComputeDeviceAcquisition` and `ComputeDeviceCapabilities`
  seams after provider release 0.16.1 repaired typed downlevel acquisition.
  The derived seven-storage-binding request preserves the full downlevel
  descriptor; raw adapter and feature methods are deleted. Nextest serializes
  only provider-acquiring tests through `gpu-device`, eliminating process-level
  WGPU device races while CPU tests remain concurrent. Evidence tier:
  compile-time API removal and empty source scan; value-semantic typed-limit
  regression; cfd-core GPU 245/245, cfd-math GPU 362/362, cfd-2d GPU 570/570
  (27 pre-existing skips), root integration 26/26; warning-denied touched
  targets; doctest/rustdoc; and SemVer's expected major-only classification.
  The root all-target example lint baseline is independently tracked by
  CFD-EXAMPLE-CLIPPY-1.

- 2026-07-17 (resolved): root `cfd-suite --all-targets` Clippy no longer
  reports the 29 diagnostics formerly distributed across seven validation
  examples. Four retained examples now execute provider-owned cfd-1d/cfd-2d
  calculations; three unreferenced static reports are deleted rather than
  presenting hardcoded validation outputs. Evidence tier: executable examples
  plus warning-denied all-target Clippy.

- 2026-07-17 (open): `BifurcationSolver3D` builds an unlabeled SDF volume mesh
  but integrates daughter flow only across `outlet_0` and `outlet_1` labels.
  The resulting zero daughter flows are deterministic, so the invalid root FEM
  example is deleted. CFD-3D-BIFURCATION-BOUNDARIES-1 owns the upstream mesh
  terminal-facet contract and cfd-3d flow regressions. Evidence tier: direct
  executable reproduction and source inspection of mesh construction and
  label-based integration.

- 2026-07-17: `cfd-core::compute::gpu::GpuContext::synchronize` now delegates
  completion to Hephaestus `ComputeDevice::synchronize`; `GpuContext` no longer
  exposes raw WGPU device, queue, or limit fields; and cfd-2d creates its
  Poisson solver through `GpuPoissonSolver::from_context`. Evidence tier:
  compile-time provider integration plus GPU-enabled value-semantic regression
  coverage (244/244 cfd-core and 2/2 accelerated cfd-2d nextest),
  warning-denied cfd-core all-target Clippy, and exact source audits with no
  `device.poll`, old Poisson constructor, context device/queue access, or
  public raw-buffer accessor. The remaining adapter/feature introspection risk
  is resolved by the 0.3.0 typed capability boundary above.

- **Verification closure (2026-07-17)**: `cargo doc -p cfd-core --no-deps
  --features gpu --locked` completes warning-clean after the final raw-buffer
  visibility change. The final source state is verified by cfd-core GPU
  nextest 244/244, cfd-2d accelerated nextest 2/2, warning-denied cfd-core/
  cfd-2d Clippy, and the package documentation gate.

- **SemVer classification (2026-07-17)**: Git-baseline semver checks identify
  the intentional removals as three breaking API classes under a minor-change
  assumption. The explicit major-change classification passes. CFDrs remains
  pre-1.0, so the workspace advances from `0.1.0` to `0.2.0` and records the
  migration in the `0.2.0` changelog section without retaining a compatibility
  surface.

- 2026-07-17: Preserved the stale peer's valid Leto source revision
  `6aedde0c7835238867d6f3cd17b030f7e69cb6f2`, which is merged on Leto `main`,
  and advanced its Moirai companion pin to merged `main`
  `5ead788c70c728d971237d7afa0b915ea7cf87e3`. Locked metadata, all-feature
  `cfd-schematics` check, and warning-denied Clippy pass; evidence tier is
  compile-time integration.

- 2026-07-16: Updated the workspace Moirai source pin to merged `main`
  `5ead788c70c728d971237d7afa0b915ea7cf87e3`. Locked metadata resolves Moirai
  0.4 and Themis 0.10; `cfd-schematics` all-feature check, warning-denied
  Clippy, focused nextest, doctests, and docs pass. The source- and
  test-level evidence is compile-time integration plus value-semantic test
  coverage.
- 2026-07-16: Removed the `cfd-schematics` strict-Clippy baseline in the
  touched test/example cone. Direct geometry/phase values use exact
  bit-pattern assertions, and the Venturi example exposes named physical
  fields instead of positional tuple entries. The workspace-wide formatter
  remains blocked by unrelated pre-existing formatting in
  `crates/cfd-schematics/src/error.rs`; touched files pass `rustfmt`.

- 2026-07-10: Removed stale allowlist entries for `cfd-1d/src/scalar.rs` and
  `cfd-3d/src/scalar.rs`. Both seams are provider-native and contain no legacy
  dependency tokens; the active `cfd-core` compute-dispatch diff is untouched.

## Sprint 2026-07-10: dead GPU pipeline surface removal

- **Resolved redundant lifecycle ownership**: The unconsumed
  `GpuPipelineManager` independently compiled shaders, created layouts,
  reconstructed untyped uniform vectors, stored name-keyed pipelines, built
  bind groups, and submitted raw WGPU commands. Every live CFD operation now
  uses a typed Hephaestus kernel, so the manager was deleted.
- **Resolved obsolete abstraction seam**: `GpuKernel<T>` exposed raw WGPU
  devices, pipelines, encoders, and the legacy `KernelParams` shape but had no
  implementors after the operation-family migrations. The trait and the
  unconsumed context pipeline constructor were deleted rather than retained as
  compatibility layers.
- **Evidence tier**: Compile-time integration and static consumer/provider
  audit. No `GpuPipelineManager`, `GpuKernel`, shader-module creation, or
  compute-pipeline creation remains in `cfd-core::compute::gpu`; cfd-core
  nextest passes 243/243; warning-denied all-target clippy, GPU/no-default and
  downstream cfd-2d checks, 3/3 doctests, warning-clean core docs, and migration
  audit pass.
- **Residual risk**: `GpuBuffer` and public raw device/queue fields still expose
  consumer-side WGPU storage infrastructure. Audit their live consumers and
  collapse them onto Hephaestus typed buffers next.

---

## Sprint 2026-07-10: GPU turbulence Hephaestus consolidation

- **Resolved duplicated infrastructure**: `GpuSmagorinskyKernel<T>`,
  `GpuDesKernel<T>`, and `GpuTurbulenceCompute` independently owned raw shader,
  pipeline, bind-group, buffer-cache, and dispatch state. One
  `GpuTurbulenceCompute` now owns three typed Hephaestus kernels and writes into
  caller-owned slices.
- **Resolved fake precision and degradation**: The cfd-2d Smagorinsky model
  converted its concrete f64 fields to f32, computed on GPU, converted back to
  f64, and silently used CPU when provider creation failed. That invalid bridge
  and its `use_gpu` configuration were removed. The provider facade now states
  its real f32 contract; the f64 domain model remains native-precision CPU.
- **Resolved false and unreachable contracts**: DES accepted two velocity
  buffers but never used them; its API now exposes the actual grid-cutoff
  quantity. Rectangular wall distance existed only as an unreachable shader
  entry point and now has a direct operation. The hardcoded device-type
  speedup estimator and wall-clock threshold test were deleted; Criterion is
  the only performance instrument.
- **Resolved topology drift**: Turbulence kernels now live in
  `kernels/turbulence/{mod,kernel,smagorinsky,des_grid_scale,wall_distance}`;
  the facade remains the single consumer-facing module. Shared uniform layout,
  field validation, provider transfers, and dispatch live once.
- **Evidence tier**: Type-level grid/provider contracts, exact linear-strain,
  constant-grid-scale, and rectangular-distance tests, typed negative tests,
  cross-crate compilation, and static audit. Focused nextest passes 4/4; full
  core passes 243/243; full cfd-2d passes 570/570 with 27 existing skips;
  warning-denied all-target clippy, legacy benchmark compilation, doctests,
  and migration audit pass. cfd-core docs are warning-clean; cfd-2d docs retain
  eight unrelated pre-existing intra-doc-link warnings.
- **Residual risk**: `GpuBuffer` and public raw device/queue fields still expose
  consumer-side WGPU storage infrastructure.

---

## Sprint 2026-07-10: GPU pressure Hephaestus dispatch

- **Resolved cosmetic execution surface**: The old
  `GpuPressureKernel<T>` returned `UnsupportedOperation`, stored a raw shader
  module, and advertised arbitrary precision over `f32` WGSL storage. The
  replacement executes weighted-Jacobi iteration and residual evaluation
  through Hephaestus typed multi-storage kernels.
- **Resolved unreachable computation**: The previous residual entry point had
  no Rust dispatch path. `GpuPressureKernel::residual` now exposes the absolute
  pointwise Poisson residual with explicit zero boundary semantics.
- **Resolved boundary defect**: Boundary selection previously copied along the
  first matching axis, so edges and corners could copy another boundary cell.
  Iteration now clamps all three coordinates to the nearest interior cell,
  enforcing the homogeneous Neumann contract at faces, edges, and corners.
- **Resolved numerical and structural drift**: `PressureConfig` validates 3D
  dimensions, finite positive spacing, representable inverse-square factors,
  and weighted-Jacobi relaxation in `(0, 1]`. The family now lives under
  `kernels/pressure/{mod,kernel,tests}` with one WGSL leaf per operation. Shared
  3D dispatch-grid construction now has one home in the kernel-family parent.
- **Evidence tier**: Type-level configuration boundary, exact quadratic
  Poisson identities, typed negative tests, compile-time provider integration,
  and static source audit. Focused pressure nextest passes 6/6; full
  `cfd-core` nextest passes 247/247; GPU/no-default checks, warning-denied
  clippy, 3/3 doctests, docs, and migration audit are clean.
- **Residual risk**: The generic pipeline manager still owns raw
  WGPU lifecycle code. Continue by complete operation family without reviving
  the deleted trait shape.

---

## Sprint 2026-07-10: GPU velocity Hephaestus dispatch

- **Resolved cosmetic execution surface**: The old
  `GpuVelocityKernel<T>` returned `UnsupportedOperation`, stored a raw shader
  module, and advertised arbitrary precision over `f32` WGSL storage. The
  replacement executes both SIMPLE correction and pressure-source divergence
  through Hephaestus typed multi-storage kernels.
- **Resolved unreachable computation**: The previous WGSL contained a
  divergence entry point that the Rust type never compiled or dispatched. The
  public `divergence_source` operation now assembles
  `(density / dt) * divergence(velocity)` as documented.
- **Resolved device-contract mismatch**: Real correction requires three input
  velocity buffers, pressure, and three output buffers. `GpuContext` now asks
  Hephaestus for the derived seven-storage-buffer limit during acquisition, so
  unsupported adapters fail at the provider boundary rather than during
  pipeline creation.
- **Resolved physical and topology drift**: `VelocityConfig` validates the 3D
  centered-difference grid, finite positive spacing, timestep, density, checked
  element count, and representable coefficients. The family now lives under
  `kernels/velocity/{mod,kernel,tests}` with correction and divergence WGSL
  leaves separated by responsibility.
- **Evidence tier**: Type-level configuration and device-capability boundaries,
  exact analytical linear-field tests, typed negative tests, compile-time
  provider integration, and static source audit. Focused velocity nextest
  passes 5/5; full `cfd-core` nextest passes 242/242; GPU and no-default checks,
  warning-denied clippy, 3/3 doctests, docs, and migration audit are clean.
- **Residual risk**: Turbulence and the generic pipeline manager
  still own raw WGPU lifecycle code. Continue by complete operation family
  without reviving the deleted trait shape.

---

## Sprint 2026-07-10: GPU diffusion Hephaestus dispatch

- **Resolved cosmetic execution surface**: The old
  `GpuDiffusionKernel<T>` never computed through `ComputeKernel`; it returned
  `UnsupportedOperation`, stored a raw shader module, and exposed arbitrary
  scalar precision over an `f32` WGSL buffer contract. The replacement is a
  real `GpuDiffusionKernel` compiled and dispatched through Hephaestus typed
  multi-storage bindings.
- **Resolved physical contract gap**: `DiffusionConfig` validates dimensions,
  checked element count, `u32` representability, positive finite spacing,
  nonnegative finite timestep/diffusivity, and the exact three-dimensional
  forward-Euler stability condition before dispatch. Execution validates field
  lengths and rejects non-finite scalar input.
- **Resolved topology and layout drift**: The operation family now lives in
  `kernels/diffusion/{mod,kernel,tests}` with its WGSL source colocated beside
  the kernel. A pair of aligned four-lane parameter blocks is the SSOT for the
  Rust/WGSL uniform layout. The old flat Rust/WGSL files were deleted.
- **Evidence tier**: Type-level configuration boundary, exact analytical value
  tests, typed negative tests, compile-time provider integration, and static
  source audit. Focused diffusion nextest passes 4/4; full `cfd-core` nextest
  passes 238/238; GPU and no-default checks pass; warning-denied all-target
  clippy passes; doctests pass 3/3; docs are warning-clean; the migration
  allowlist and provider/fake-generic audits are clean.
- **Residual risk**: Pressure, turbulence, and the generic pipeline
  manager still own raw WGPU lifecycle code. Continue by complete operation
  family without reviving the deleted trait shape.

---

## Sprint 2026-07-10: GPU advection Hephaestus dispatch

- **Resolved cosmetic execution surface**: The old
  `GpuAdvectionKernel<T>` never computed through `ComputeKernel`; it returned
  `UnsupportedOperation`, stored a raw shader module, and exposed arbitrary
  scalar precision over an `f32` WGSL buffer contract. The replacement is a
  real `GpuAdvectionKernel` compiled and dispatched through Hephaestus typed
  multi-storage bindings.
- **Resolved primitive configuration boundary**: `AdvectionConfig` validates
  dimensions, checked element count, `u32` representability, positive finite
  spacing, and nonnegative finite timestep once. Dispatch validates all field
  lengths and finite values and enforces the unsplit first-order upwind CFL
  condition before GPU work is submitted.
- **Resolved topology drift**: The operation family now lives in
  `kernels/advection/{mod,kernel,tests}` with its WGSL source colocated beside
  the kernel. The old flat Rust/WGSL sibling files and separate raw-pipeline
  integration path were deleted.
- **Evidence tier**: Type-level configuration boundary, exact analytical value
  tests, typed negative tests, compile-time provider integration, and static
  source audit. Focused advection nextest passes 6/6; full `cfd-core` nextest
  passes 234/234; GPU and no-default checks pass; warning-denied all-target
  clippy passes; doctests pass 3/3; docs are warning-clean; the migration
  allowlist and provider/fake-generic audits are clean.
- **Residual risk**: Velocity, pressure, turbulence, and the generic
  pipeline manager still own raw WGPU lifecycle code. Continue by complete
  operation family without reviving the deleted trait shape.

---

## Sprint 2026-07-10: 2D Laplacian Hephaestus dispatch

- **Resolved duplicate GPU infrastructure**: `Laplacian2DKernel` now compiles
  and dispatches the CFD stencil through Hephaestus
  `WgslMultiStorageKernel`, typed storage bindings, dispatch-grid validation,
  provider buffers, and provider-owned transfer synchronization. Consumer-owned
  pipelines, layouts, bind groups, uniform allocation, staging buffers,
  mapping channels, polling, and timeouts were deleted.
- **Resolved hidden degradation**: Kernel construction and execution return
  typed errors. Small inputs and provider failures no longer trigger silent CPU
  recomputation; the CPU stencil is compiled only for differential tests.
- **Resolved fake generic**: `GpuLaplacianOperator2D` now exposes its actual
  WGSL `f32` contract. The prior generic path converted spacing to `f32` while
  reinterpreting arbitrary `T` device storage as `f32`.
- **Resolved dimensional ambiguity**: The CFD facade, kernel wrapper, and
  matrix-free operator now carry Aequitas `Length<f32>` spacing. CFDrs
  validates finite positive metre values, and Hephaestus converts quantities
  once into its provider-owned POD parameter block; no unit wrapper enters
  device storage or the WGSL loop.
- **Resolved periodic-boundary drift**: Endpoint-inclusive periodic neighbors
  now wrap to the opposite inner point (`nx-2`/`1`, `ny-2`/`1`) as documented
  and as implemented by the independent CPU oracle.
- **Resolved diagnostic debt**: The clean Anderson test cone now uses
  `Option::map_or`, closing the warning surfaced by the full all-target gate.
- **Evidence tier**: Compile-time provider integration, exact and
  analytically-bounded GPU/CPU differential tests, typed negative tests, and
  static source audit. No-default and GPU checks pass; the typed-spacing
  focused Laplacian nextest passes 13/13; full `cfd-core` and `cfd-math`
  nextest for the provider migration pass 231/231 and
  362/362; warning-denied all-target clippy passes; doctests pass 6/6 with 3
  intentionally ignored; package docs are warning-clean.
- **SemVer evidence limit**: `cargo-semver-checks check-release -p cfd-core -p
  cfd-math` cannot retrieve a baseline because neither crate is published to
  crates.io. The fallible-constructor and scalar-contract changes are manually
  classified as breaking and documented with migration instructions.
- **Residual risk**: Advection, diffusion, velocity, pressure, turbulence, and
  generic pipeline modules still own raw WGPU orchestration. Migrate each
  operation family independently through the existing Hephaestus authored
  kernel seam.

---

## Sprint 2026-07-10: cfd-core GPU arithmetic Hephaestus elementwise

- **Resolved duplicated provider orchestration**: Removed `FieldAddKernel`,
  `FieldMulKernel`, their local WGSL, raw WGPU pipelines, staging buffers,
  polling channels, and timeout logic. `GpuFieldOps` now uploads typed buffers,
  dispatches Hephaestus `AddOp`/`MulOp`, and downloads through the provider seam.
- **Resolved hidden degradation**: GPU failures now propagate through the typed
  `Error::GpuProvider` variant. No small-input or failure-triggered CPU branch
  remains in the arithmetic path.
- **Resolved partial-workgroup defect**: Exact-value regressions cover lengths
  1, 63, 64, 65, and 257, closing the former aligned-readback panic and proving
  tail dispatch through the real provider.
- **Evidence tier**: Compile-time integration, exact empirical differential
  tests, typed negative tests, and static source audit. No-default and GPU
  checks pass; all-target GPU clippy passes; cfd-core nextest passes 230/230
  with no skips; doctests pass 5/5; docs are warning-clean.
- **Residual risk**: Other cfd-core raw WGPU kernels and their independent CPU
  fallback policies remain outside this arithmetic slice and require separate
  provider-owned migrations.

---

## Sprint 2026-07-07: cfd-1d/cfd-3d Eunomia identity seam

- **Resolved direct scalar dependency**:
  `Cargo.toml`, `crates/cfd-1d/Cargo.toml`, and
  `crates/cfd-3d/Cargo.toml` no longer declare direct `num-traits`
  dependencies for the 1D/3D solver scalar seams.
- **Resolved identity ownership**:
  `Cfd1dScalar` and `Cfd3dScalar` now expose `zero()` and `one()` through the
  Eunomia `NumericElement` constants already required by the crate-local
  scalar contracts, removing `num_traits::{Zero,One}` as a supertrait
  requirement.
- **Evidence tier**: compile-time integration, empirical nextest coverage,
  touched-file formatting, and static source/manifest audit. In
  `D:/atlas/repos/CFDrs`, touched-file rustfmt passed; `rustup run nightly
  cargo check -p cfd-1d -p cfd-3d` passed; direct residue scan found no
  `num_traits` or direct `num-traits` hits in the touched manifests/scalar
  cones; and `rustup run nightly cargo nextest run -p cfd-1d -p cfd-3d
  --status-level fail` passed 1122/1122 with one existing slow 3D
  mesh-convergence validation.
- **Residual risk**: package-wide fmt and all-targets clippy remain blocked by
  pre-existing unrelated formatting/lint debt outside this slice. Lockfile
  `num-traits` entries, if present, are transitive provider dependencies owned
  by upstream crates rather than direct CFDrs scalar-seam dependencies.

---

## Sprint 2026-07-05: cfd-1d solver workspace Leto vector storage

- **Resolved workspace vector storage boundary**:
  `crates/cfd-1d/src/solver/core/workspace.rs` no longer imports or stores
  nalgebra `DVector`. `SolverWorkspace::{rhs,last_solution,linear_solution}`
  are now `leto::Array1<T>` buffers, preserving the reusable allocation model
  while removing nalgebra from the workspace state.
- **Resolved duplicated conversion helpers**:
  `crates/cfd-1d/src/solver/core/vector_bridge.rs` is now the private SSOT for
  nalgebra<->Leto vector conversion, Leto vector norm calculation, and dense
  Leto workspace copies. The previous ad hoc conversion helpers in
  `mod.rs`, `linear_system.rs`, and `anderson_acceleration.rs` were removed.
- **Resolved RHS validation/storage path**:
  `MatrixAssembler::assemble` returns a Leto RHS array, and solver validation,
  residual calculation, convergence RHS norms, row equilibration, and
  iterative linear solver calls consume Leto RHS/initial-guess buffers directly.
- **Evidence tier**: compile-time provider integration, focused empirical
  tests, clippy, rustdoc execution, and static source audit. In
  `D:/atlas/repos/CFDrs`, `rustup run nightly cargo check -p cfd-1d --tests`
  passed; `rustup run nightly cargo clippy -p cfd-1d --lib -- -D warnings`
  passed; focused solver-core nextest passed 10/10; cfd-1d doc generation
  completed with the existing 11 rustdoc warnings; targeted residue scan found
  no nalgebra `DVector` hits in `workspace.rs`.
- **Residual risk**: `linear_system.rs` still exposes nalgebra `DVector`
  solutions, dense `DMatrix` fallback storage, and the nalgebra-sparse matrix
  input boundary. `matrix_assembly.rs` still returns
  `nalgebra_sparse::CsrMatrix`. `solver_detection.rs` and
  `anderson_acceleration.rs` still operate on nalgebra solution vectors at the
  current linear-solver boundary.

---

## Sprint 2026-07-05: cfd-1d NetworkState Leto vector boundary

- **Resolved public state vector boundary**:
  `crates/cfd-1d/src/solver/core/state.rs` no longer imports or exposes
  nalgebra `DVector`. `NetworkState::{pressures,flow_rates}` now use
  `leto::Array1<T>`, and `NetworkState::new`/`from_network` build Leto arrays
  directly from validated network dimensions and network pressure/flow slices.
- **Resolved value-semantic state coverage**:
  focused unit tests assert Leto array shapes, zero initialization, direct
  indexed values, clone preservation, and time mutation semantics for the
  migrated state type.
- **Evidence tier**: compile-time provider integration, focused unit tests,
  clippy, rustdoc execution, and static source audit. In
  `D:/atlas/repos/CFDrs`, cfd-1d test-target check passed; cfd-1d lib clippy
  passed; focused state nextest passed 2/2; cfd-1d doc generation completed
  with the existing rustdoc warnings; targeted residue scan found no `DVector`
  in `state.rs`. The upstream Atlas blocker was closed in
  `D:/atlas/repos/eunomia` by adding `FloatElement::acos` with native f64
  primitive/wrapper overrides.
- **Residual risk**: The broader cfd-1d linear-system matrix/solution,
  sparse-assembly, solver-detection, and Anderson solution-vector boundaries
  still retain nalgebra vector/matrix/sparse storage.

---

## Sprint 2026-07-05: cfd-1d convergence checker Leto vector boundary

- **Resolved convergence-checker vector boundary**:
  `crates/cfd-1d/src/solver/core/convergence.rs` no longer imports or exposes
  nalgebra `DVector`. `ConvergenceChecker::{check,has_converged,
  has_converged_dual}` now takes `leto::Array1<T>` and computes finite checks,
  L2 norms, and L2 delta norms directly from vector values.
- **Resolved error semantics for shape mismatch**:
  mismatched current/previous solution lengths now return a typed
  `Error::InvalidInput` with the offending lengths instead of relying on
  nalgebra vector subtraction shape behavior.
- **Evidence tier**: compile-time provider integration, focused empirical
  tests, clippy, rustdoc execution, and static source audit. In
  `D:/atlas/repos/CFDrs`, cfd-1d test-target check passed; focused cfd-1d
  clippy for lib plus `solver_core_tests` passed; focused convergence nextest
  passed 20/20; cfd-1d doc generation completed; targeted residue scan found
  no `DVector` in the checker implementation or migrated checker tests.
- **Residual risk**: This is not a full 1D solver storage migration. The
  solver loop still stores workspaces as nalgebra vectors and converts once at
  the convergence-check boundary; `linear_system`, `matrix_assembly`,
  `solver_detection`, and `anderson_acceleration` still retain nalgebra
  `DVector`/`DMatrix`/`nalgebra_sparse` storage. `cargo doc -p cfd-1d
  --no-deps` still emits 11 pre-existing broken/private intra-doc link
  warnings outside this convergence slice.

---

## Sprint 2026-07-05: cfd-suite direct numeric dependency removal

- **Resolved direct legacy numeric dependency**:
  the root `Cargo.toml` no longer declares `num-traits` or `simba` in the
  workspace dependency catalog or root package dependency list. CFDrs source
  and member manifests now rely on Eunomia/Atlas scalar contracts rather than
  direct `num_traits` or `simba` imports.
- **Evidence tier**: static source/member-manifest audit plus empirical
  package metadata verification. In
  `D:/atlas/repos/CFDrs`, the targeted residue scan found no `num_traits` or
  `simba` source imports and no direct `num-traits.workspace` or
  `simba.workspace` member dependencies after the removal. Root package
  metadata passed.
- **Residual risk**: `num-traits` and `simba` remain in `Cargo.lock` as
  transitive dependencies of upstream crates. This slice removes CFDrs direct
  ownership; it does not claim the full dependency graph is free of transitive
  numeric provider crates. `rustup run nightly cargo check -p cfd-suite
  --no-default-features` stayed queued behind concurrent shared Cargo
  cache/build locks and was terminated without compiler diagnostics; nextest
  did not run.

---

## Sprint 2026-07-05: cfd-suite direct Crossbeam removal

- **Resolved direct legacy concurrency dependency**:
  the root `Cargo.toml` no longer declares `crossbeam` in the workspace
  dependency catalog, and `crates/cfd-3d/Cargo.toml` no longer depends on
  `crossbeam.workspace`. CFDrs source and member manifests now rely on
  Moirai-facing concurrency surfaces without direct Crossbeam imports.
- **Evidence tier**: static source/member-manifest audit plus empirical
  package metadata, compile, nextest, and diff whitespace validation. In
  `D:/atlas/repos/CFDrs`, the targeted residue scan found no `crossbeam`
  source imports and no direct `crossbeam.workspace` member dependencies after
  the removal. Root package metadata passed. `rustup run nightly cargo check -p
  cfd-3d --no-default-features` passed. `rustup run nightly cargo nextest run
  -p cfd-3d --no-default-features` passed 394/394 tests, with existing
  mesh-convergence validation marked slow at 19.7s.
- **Residual risk**: this slice removes direct CFDrs ownership only. It does
  not claim the full dependency graph is Crossbeam-free; any transitive
  Crossbeam dependencies remain upstream-owned. The slow mesh-convergence test
  remains a performance follow-up outside this dependency-ownership slice.

---

## Sprint 2026-07-05: public sparse/linear-solver Leto boundary

- **Resolved public sparse storage boundary**:
  `crates/cfd-math/src/sparse/mod.rs` now exposes
  `leto_ops::CsrMatrix<T>` as the public `SparseMatrix<T>` alias. Sparse
  builders, assembly helpers, sparse operations, and sparse tests construct
  and operate on Leto CSR directly rather than converting through
  nalgebra-sparse.
- **Resolved public solver vector/matrix boundary**:
  `LinearOperator`, `Preconditioner`, `LinearSolver`, direct solver,
  solver-chain, and migrated solver fixtures use Leto CSR and Leto
  `Array1<T>` boundaries in the requested cone. `cfd-validation::numerical`
  now stores computed and analytical validation vectors as `leto::Array1<T>`
  and uses Leto CSR for linear-solver validation test cases.
- **Evidence tier**: compile-time provider integration, package empirical
  regression tests, clippy, rustdoc, and static source audit. In
  `D:/atlas/repos/CFDrs`, cfd-math check passed; cfd-math test-target check
  passed; cfd-math all-target clippy passed; cfd-math doc passed; cfd-math
  nextest passed 361/361; cfd-validation check passed; cfd-validation
  all-target clippy passed; cfd-validation doc passed after fixing a stale
  intra-doc link; the targeted residue scan found no
  `nalgebra_sparse::CsrMatrix`, public nalgebra-sparse re-export, `DVector`,
  `row_offsets()`, `try_from_csr_data`, `CooMatrix`, or `nalgebra::` matches
  under the migrated sparse/linear-solver/validation files.
- **Residual risk**: The requested sparse/linear-solver public boundary has no
  known nalgebra sparse/vector holdouts in the scanned cone. Full
  `cfd-validation` package nextest remains blocked by the existing venturi
  cross-fidelity convergence tests
  `option2_selected_45um_geometry_routes_to_fallback_and_converges` and
  `microventuri_35um_case_produces_converged_informative_2d_result`, which are
  outside this boundary. Broader CFDrs provider migration still has nalgebra
  residue in other crates and contexts outside this slice.

---

## Sprint 2026-07-05: cfd-math IncompleteCholesky Leto CSR

- **Resolved IncompleteCholesky sparse provider boundary**:
  `crates/cfd-math/src/linear_solver/preconditioners/cholesky.rs` now stores
  and factorizes `leto_ops::CsrMatrix`, uses Leto row-value lookup for
  symmetry/factorization reads, and performs triangular substitutions through
  Leto CSR row views.
- **Resolved direct nalgebra residue in the migrated boundary**:
  `cholesky.rs` no longer imports `nalgebra_sparse` or uses `row_offsets()`,
  `try_from_csr_data`, `get_entry()`, or `SparseEntry`. The integration
  preconditioner edge-case test now constructs the Cholesky negative case with
  Leto CSR directly.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest, clippy, static source audit, and format hygiene. In
  `D:/atlas/repos/CFDrs`, cfd-math fmt passed; cfd-math lib check passed;
  cfd-math all-target check passed; cfd-math lib/tests clippy passed;
  cfd-math all-target clippy passed; cholesky-filter nextest passed 5/5;
  preconditioner-filter nextest passed 76/76; targeted Cholesky residue scan
  found no nalgebra sparse/direct CSR/row-offset/get-entry residue in
  `cholesky.rs` or `tests/preconditioner_edge_cases.rs`.
- **Residual risk**: Schwarz, direct solver, remaining transitional solver
  fixtures, and the shared `crate::sparse::SparseMatrix` solver matrix
  boundary still use nalgebra-sparse. The next sparse-provider increment
  should move Schwarz or the shared solver matrix boundary.

---

## Sprint 2026-07-05: cfd-math ILU Leto CSR

- **Resolved ILU sparse provider boundary**:
  `crates/cfd-math/src/linear_solver/preconditioners/ilu` now stores and
  factorizes `leto_ops::CsrMatrix` for ILU(0), ILU(k), and triangular solve
  application.
- **Resolved direct nalgebra residue in the migrated boundary**:
  the ILU module no longer imports `nalgebra_sparse` or uses `row_offsets()`,
  `try_from_csr_data`, `get_entry()`, or `SparseEntry`. Source tests,
  integration edge tests, `LinearSolverChain`, and Schwarz local-solver setup
  pass Leto CSR matrices into ILU or convert once at the still-transitional
  boundary.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest, clippy, static source audit, and format hygiene. In
  `D:/atlas/repos/CFDrs`, cfd-math fmt passed; cfd-math lib check passed;
  cfd-math all-target check passed; cfd-math lib/tests clippy passed;
  cfd-math all-target clippy passed; ilu-filter nextest passed 21/21;
  preconditioner-filter nextest passed 76/76; linear_solver::tests nextest
  passed 53/53; targeted ILU residue scan found no nalgebra sparse/direct
  CSR/row-offset/get-entry residue in `preconditioners/ilu`.
- **Residual risk**: Schwarz, direct solver, remaining transitional solver
  fixtures, and the shared `crate::sparse::SparseMatrix` solver matrix
  boundary still use nalgebra-sparse. The next sparse-provider increment
  should move Schwarz or the shared solver matrix boundary.

---

## Sprint 2026-07-05: cfd-math SSOR Leto CSR

- **Resolved SSOR sparse provider boundary**:
  `crates/cfd-math/src/linear_solver/preconditioners/ssor.rs` now stores
  `leto_ops::CsrMatrix` and performs both SOR sweeps through Leto CSR row
  views.
- **Resolved direct nalgebra residue in the migrated boundary**:
  `ssor.rs` no longer imports `nalgebra_sparse` or uses `row_offsets()`,
  `try_from_csr_data`, or `get_entry()`. Source preconditioner edge tests
  convert legacy solver CSR fixtures once into Leto CSR before constructing
  SSOR.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest, clippy, static source audit, and format hygiene. In
  `D:/atlas/repos/CFDrs`, cfd-math fmt passed; cfd-math lib check passed;
  cfd-math all-target check passed; cfd-math lib/tests clippy passed;
  cfd-math all-target clippy passed; ssor-filter nextest passed 5/5;
  preconditioner-filter nextest passed 76/76; targeted SSOR residue scan found
  no nalgebra sparse/direct CSR/row-offset/get-entry residue in `ssor.rs`.
- **Residual risk**: Schwarz, direct solver, integration-test
  fixtures, and the shared `crate::sparse::SparseMatrix` solver matrix
  boundary still use nalgebra-sparse. The next sparse-provider increment
  should move another whole preconditioner family rather than adding SSOR
  compatibility constructors.

---

## Sprint 2026-07-05: cfd-math Basic Preconditioner Leto CSR

- **Resolved basic preconditioner sparse provider boundary**:
  `crates/cfd-math/src/linear_solver/preconditioners/basic.rs` now constructs
  Jacobi and SOR preconditioners from `leto_ops::CsrMatrix`, with Leto CSR
  diagonal extraction, row access, and Leto `Array1` application buffers.
- **Resolved direct nalgebra residue in the migrated boundary**:
  `basic.rs` no longer imports `nalgebra_sparse`, `SparseMatrixExt`, or uses
  `row_offsets()`, `try_from_csr_data`, or `get_entry()`. Tests that still
  exercise legacy solver matrices convert those fixtures once into Leto CSR
  before constructing Jacobi/SOR.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest, clippy, static source audit, and format hygiene. In
  `D:/atlas/repos/CFDrs`, cfd-math fmt passed; cfd-math lib check passed; core
  solver test check passed; cfd-math all-target check passed; cfd-math
  lib/tests clippy passed; core solver test clippy passed; cfd-math all-target
  clippy passed; linear_solver::tests nextest passed 53/53; core_solver_tests
  nextest passed 4/4; preconditioner-filter nextest passed 76/76; targeted
  basic-preconditioner residue scan found no nalgebra sparse/direct CSR/row
  offset/get-entry residue in `basic.rs`.
- **Residual risk**: Schwarz, direct solver, and the shared
  `crate::sparse::SparseMatrix` solver matrix boundary still use
  nalgebra-sparse. The next sparse-provider increment should move another
  whole preconditioner family, not add compatibility constructors to basic
  preconditioners.

---

## Sprint 2026-07-05: cfd-math AMG/Coarsening Leto CSR

- **Resolved AMG sparse provider boundary**:
  `crates/cfd-math/src/linear_solver/preconditioners/multigrid` now uses a
  local `leto_ops::CsrMatrix` sparse alias and construction/value helpers for
  AMG setup, coarsening, interpolation, smoothers, cycles, sparse products,
  transpose, row access, and SpMV.
- **Resolved benchmark/test holdouts**:
  `crates/cfd-math/benches/coarsening_bench.rs`,
  `crates/cfd-math/benches/algebraic_distance_bench.rs`,
  `crates/cfd-math/tests/amg_coarsening_tests.rs`, and AMG preconditioner
  construction in `crates/cfd-math/tests/amg_integration_test.rs` now use Leto
  CSR. The AMG integration test keeps the Krylov solver matrix on the existing
  solver sparse boundary but assembles it through `SparseMatrixBuilder` instead
  of direct nalgebra-sparse APIs.
- **Evidence tier**: compile-time provider integration, focused empirical AMG
  nextest, focused clippy, static source audit, and format hygiene. In
  `D:/atlas/repos/CFDrs`, cfd-math fmt passed; cfd-math lib check passed;
  focused checks passed for `amg_coarsening_tests`, `amg_integration_test`,
  `coarsening_bench`, and `algebraic_distance_bench`; focused clippy passed
  for cfd-math lib, both AMG tests, and both migrated benches; AMG integration
  nextest passed 5/5; AMG-filter nextest passed 6/6; multigrid::coarsening
  nextest passed 10/10; targeted AMG/coarsening sparse scan found no
  `nalgebra_sparse`, `CooMatrix`, `try_from_csr_data`, `get_entry`,
  `crate::sparse::spmv`, `spmv_array`, or `row_offsets()` residue in the
  migrated boundary.
- **Residual risk**: cfd-math still exposes nalgebra-sparse through the broader
  `crate::sparse::SparseMatrix` solver/direct/preconditioner boundary.
  `LinearSolverChain` still converts that solver matrix boundary once into
  Leto CSR before constructing AMG. The next sparse-provider increment should
  move a whole solver/preconditioner boundary to direct Leto CSR rather than
  adding benchmark-local conversions.

---

## Sprint 2026-07-05: cfd-math Leto CSR Benchmarks

- **Resolved benchmark nalgebra CSR construction**:
  `crates/cfd-math/benches/{spmv_bench,cg_bench,math_benchmarks}.rs` now
  construct `leto_ops::CsrMatrix` values directly with `from_parts` instead of
  constructing `nalgebra_sparse::CsrMatrix` values through
  `try_from_csr_data`.
- **Resolved SpMV benchmark dispatch**:
  `spmv_bench.rs` now measures the direct Leto CSR `LinearOperator::apply`
  path rather than the legacy `cfd_math::sparse::spmv` compatibility helper.
- **Evidence tier**: compile-time provider integration, focused benchmark
  clippy, sparse empirical nextest, static source audit, and format hygiene.
  In `D:/atlas/repos/CFDrs`, cfd-math fmt passed; focused checks passed for
  `spmv_bench`, `cg_bench`, and `math_benchmarks`; focused clippy passed for
  all three migrated benches; cfd-math all-target check passed; cfd-math
  all-target clippy passed; sparse-filter nextest passed 19/19; targeted
  migrated-benchmark scan found no `nalgebra_sparse`, `CooMatrix`, `DVector`,
  `DMatrix`, `cfd_math::sparse::spmv`, or `try_from_csr_data` residue.
- **Residual risk**: after the AMG/coarsening sparse boundary moved to Leto CSR
  in Sprint 2026-07-05, the remaining sparse-provider risk is the broader
  solver/direct/preconditioner `crate::sparse::SparseMatrix` boundary that
  still resolves to nalgebra-sparse.

---

## Sprint 2026-07-05: cfd-math Leto CSR LinearOperator

- **Resolved Atlas-native sparse operator path**:
  `crates/cfd-math/src/sparse/operations.rs` now implements
  `LinearOperator<T>` for `leto_ops::CsrMatrix<T>`, allowing iterative solvers
  to consume Leto CSR matrices directly.
- **Consolidated SpMV execution**:
  the remaining `nalgebra_sparse::CsrMatrix` SpMV boundary now converts once to
  Leto CSR and delegates to the same `try_leto_spmv` helper used by direct
  Leto CSR inputs.
- **Migrated proof path**:
  `crates/cfd-math/tests/simple_gmres_tests.rs` now constructs Leto CSR
  matrices directly and validates GMRES residuals through the direct Leto CSR
  operator path, with no nalgebra COO/CSR/vector residue in that integration
  test.
- **Evidence tier**: compile-time provider integration, focused empirical
  GMRES tests, broader sparse/GMRES empirical filters, clippy over the focused
  test and all cfd-math targets, static source audit, and format hygiene. In
  `D:/atlas/repos/CFDrs`, cfd-math fmt passed; simple_gmres test check
  passed; simple_gmres nextest passed 3/3; simple_gmres clippy passed;
  cfd-math lib check passed; cfd-math all-target check passed; cfd-math
  all-target clippy passed; sparse-filter nextest passed 19/19; gmres-filter
  nextest passed 21/21; targeted simple GMRES scan found no
  `nalgebra_sparse`, `CooMatrix`, `nalgebra::`, `DVector`, `DMatrix`, or
  `cfd_math::sparse` residue.
- **Residual risk**: the broad cfd-math sparse matrix boundary still exposes
  `nalgebra_sparse::CsrMatrix` in sparse builders, preconditioners, AMG,
  direct solver, source tests, integration tests, and benches. This slice
  opens and verifies the direct Atlas CSR solver path; the remaining work is
  replacing those public matrix surfaces with Leto CSR in owner-specific
  follow-up slices.

---

## Sprint 2026-07-05: cfd-math Storage-Slice Closure

- **Resolved nonlinear mutable storage bridge**:
  `crates/cfd-math/src/nonlinear_solver/linalg.rs` no longer exposes
  `StorageMut`-backed mutable slice helpers, and
  `crates/cfd-math/src/nonlinear_solver/anderson.rs` performs pivot swaps
  through direct Leto array indexing.
- **Resolved multigrid storage bridges**:
  `crates/cfd-math/src/linear_solver/preconditioners/multigrid/
  interpolation.rs` indexes interpolated Leto arrays directly for
  constant-preservation error, and `smoothers.rs` tests compare expected
  vectors through indexed assertions.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest over nonlinear/multigrid tests, clippy over lib/tests and all
  cfd-math targets, static source audit, and format hygiene. In
  `D:/atlas/repos/CFDrs`, cfd-math fmt passed; cfd-math lib check passed;
  nonlinear_solver/multigrid nextest passed 46/46; cfd-math lib/tests clippy
  passed; cfd-math all-target check passed; cfd-math all-target clippy passed;
  cfd-math `src`/`tests` scan found no `leto::Storage`, `StorageMut`,
  `.storage().as_slice()`, `as_slice_mut()`, `vector_slice_mut`, or
  `matrix_slice_mut` residue.
- **Residual risk**: the Leto storage-slice cleanup is closed for cfd-math
  source/tests. Remaining cfd-math Atlas work is now broader
  nalgebra/nalgebra-sparse replacement and other provider boundaries.

---

## Sprint 2026-07-05: cfd-math Sparse/Basic Leto Array1

- **Resolved sparse vector storage bridge**:
  `crates/cfd-math/src/sparse/operations.rs` now stages SpMV input/output and
  row/column scaling buffers through direct `Array1` indexing before
  delegating to Leto CSR operations.
- **Resolved Jacobi diagonal storage bridge**:
  `crates/cfd-math/src/linear_solver/preconditioners/basic.rs` now indexes the
  Leto diagonal directly when building the inverse diagonal.
- **Preserved provider ownership**:
  sparse matrix multiply, transpose, SpMV, scaling, diagonal extraction,
  Frobenius norm, condition estimate, and diagonal-dominance checks remain
  delegated to Leto CSR operations through the existing sparse bridge. This
  slice removes CFDrs raw Leto storage access rather than reimplementing those
  operations locally.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest over sparse/preconditioner tests, clippy over lib/tests and all
  cfd-math targets, static source audit, and format hygiene. In
  `D:/atlas/repos/CFDrs`, cfd-math fmt passed; cfd-math lib check passed;
  sparse/preconditioner nextest passed 95/95; cfd-math lib/tests clippy
  passed; cfd-math all-target check passed; cfd-math all-target clippy passed;
  targeted scans found no `Storage`, `.storage().as_slice()`, `as_slice_mut()`,
  or obsolete SpMV output-contiguity diagnostic residue in the migrated
  sparse/basic files.
- **Residual risk**: remaining cfd-math storage-slice owners are
  `crates/cfd-math/src/nonlinear_solver/linalg.rs` mutable `StorageMut`
  helpers and multigrid interpolation/smoother internals.

---

## Sprint 2026-07-05: cfd-math GPU Operator Leto Array1

- **Resolved GPU operator storage bridge**:
  `crates/cfd-math/src/linear_solver/operators/gpu.rs` now stages GPU
  upload/readback through direct `Array1` indexing rather than borrowing raw
  Leto storage slices.
- **Removed contiguity assumptions**:
  the GPU linear-operator `LinearOperator::apply` implementation validates
  input/output vector lengths and no longer requires output `as_slice_mut()`
  contiguity before writing results.
- **Preserved GPU provider routing**:
  the operator still executes through the existing Hephaestus-backed
  `cfd-core` `GpuContext`/`GpuBuffer`/`Laplacian2DKernel` path; this slice
  does not introduce a local GPU fallback or duplicate kernel provider.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest over GPU-feature linear-operator tests, GPU-feature clippy over lib
  and all cfd-math targets, static source audit, and format hygiene. In
  `D:/atlas/repos/CFDrs`, cfd-math fmt passed; cfd-math GPU-feature check
  passed; GPU-feature linear_solver::operators nextest passed 5/5; cfd-math
  GPU-feature lib clippy passed; cfd-math GPU-feature all-target check passed;
  cfd-math GPU-feature all-target clippy passed; targeted scan found no
  `Storage`, `.storage().as_slice()`, `as_slice_mut()`, or obsolete output
  contiguity diagnostic residue in `operators/gpu.rs`.
- **Residual risk**: remaining cfd-math storage-slice owners are
  `crates/cfd-math/src/sparse/operations.rs` and multigrid
  smoother/interpolation internals. Full GPU migration still requires
  replacing broader CFDrs raw WGPU context/kernel ownership with Hephaestus
  abstractions where Hephaestus exposes the needed contracts.

---

## Sprint 2026-07-05: cfd-math Finite-Difference Operators Leto Array1

- **Resolved CPU operator storage bridge**:
  `crates/cfd-math/src/linear_solver/operators/{poisson,momentum}.rs` now read
  and write `Array1` values through direct indexing for 2D Laplacian, 3D
  Poisson, 1D/2D momentum, and 2D energy operators.
- **Removed contiguity assumptions**:
  the migrated CPU finite-difference operators no longer import
  `leto::Storage`, borrow input storage slices, or require output
  `as_slice_mut()` contiguity before writing stencil values.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest over the linear-operator tests, clippy over lib/tests and all
  cfd-math targets, static source audit, and format hygiene. In
  `D:/atlas/repos/CFDrs`, cfd-math fmt passed; cfd-math lib check passed;
  linear_solver::operators nextest passed 5/5; cfd-math lib/tests clippy
  passed; cfd-math all-target check passed; cfd-math all-target clippy
  passed; targeted scan found no `Storage`, `.storage().as_slice()`, or
  `as_slice_mut()` residue in the migrated operator files.
- **Residual risk**: remaining storage-slice owners are
  `crates/cfd-math/src/sparse/operations.rs`,
  `crates/cfd-math/src/linear_solver/operators/gpu.rs`, and multigrid
  smoother/interpolation internals. The GPU operator should be handled in a
  Hephaestus-specific slice rather than folded into this CPU operator cleanup.

---

## Sprint 2026-07-05: cfd-math Nonlinear Linalg Leto Array1

- **Resolved nonlinear immutable storage bridge**:
  `crates/cfd-math/src/nonlinear_solver/linalg.rs` now evaluates dot,
  add/sub, scaled add, in-place scaled add, and scale through direct `Array1`
  indexing instead of an immutable `vector_slice` helper backed by
  `leto::Storage`.
- **Resolved Anderson coefficient bridge**:
  `crates/cfd-math/src/nonlinear_solver/anderson.rs` now indexes the Anderson
  least-squares coefficient vector directly, allowing the immutable
  `vector_slice` helper to be deleted.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest over nonlinear solver tests, clippy over lib/tests and all cfd-math
  targets, static source audit, and format hygiene. In
  `D:/atlas/repos/CFDrs`, cfd-math fmt passed; cfd-math lib check passed;
  nonlinear_solver nextest passed 9/9; cfd-math lib/tests clippy passed;
  cfd-math all-target check passed; cfd-math all-target clippy passed;
  targeted scan found no immutable `.storage().as_slice()`, `leto::Storage`,
  or `vector_slice` residue in nonlinear linalg/Anderson.
- **Residual risk**: mutable dense-workspace helpers in nonlinear linalg still
  use `StorageMut`; remaining immutable storage-slice owners are sparse
  operations, linear operators, GPU operator, and multigrid internals.

---

## Sprint 2026-07-05: cfd-math Production SIMD Vector Leto Array1

- **Resolved SIMD vector storage bridge**:
  `crates/cfd-math/src/simd/vector.rs` now evaluates `simd_mul`, `simd_dot`,
  `simd_norm`, and `par_map` through direct `Array1` indexing instead of
  importing `leto::Storage` and borrowing raw storage slices.
- **Preserved execution provider**:
  the implementation still routes work through Moirai `Adaptive`
  `map_collect_index_with` and `reduce_index_with`; this slice only removes
  the Leto storage trait coupling from the SIMD vector extension.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest over the SIMD vector module, broader empirical nextest over the SIMD
  filter, clippy over lib/tests and all cfd-math targets, static source audit,
  and format hygiene. In `D:/atlas/repos/CFDrs`, cfd-math fmt passed;
  cfd-math lib check passed; simd::vector nextest passed 1/1; cfd-math
  lib/tests clippy passed; cfd-math all-target check passed; cfd-math
  all-target clippy passed; SIMD-filter nextest passed 26/26; targeted scan
  found no `Storage` or `.storage().as_slice()` residue in
  `src/simd/vector.rs`.
- **Residual risk**: source-level storage-slice owners remain in
  `crates/cfd-math/src/nonlinear_solver/linalg.rs`,
  `crates/cfd-math/src/sparse/operations.rs`,
  `crates/cfd-math/src/linear_solver/operators/{poisson,momentum,gpu}.rs`,
  and multigrid smoother/interpolation internals; those need owner-specific
  provider slices because some are sparse/gpu backend boundaries.

---

## Sprint 2026-07-05: cfd-math SIMD Integration Test Leto Array1

- **Resolved SIMD test storage bridge**:
  `crates/cfd-math/tests/simd_tests.rs` now reads the `leto_ops::spmv` result
  through direct `Array1` indexing and passes those values to the existing
  SIMD slice API, removing the final cfd-math integration-test
  `leto::Storage`/`.storage().as_slice()` bridge.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest over SIMD integration tests, broader empirical nextest over the SIMD
  filter, clippy over the focused test and all cfd-math targets, static source
  audit, and format hygiene. In `D:/atlas/repos/CFDrs`, cfd-math fmt passed;
  simd_tests check passed; simd_tests nextest passed 12/12; simd_tests clippy
  passed; cfd-math all-target check passed; cfd-math all-target clippy passed;
  SIMD-filter nextest passed 26/26; targeted scans found no `DVector`,
  nalgebra vector import, local preconditioner bridge, `Storage`, or
  storage-slice conversion residue under `crates/cfd-math/tests`.
- **Residual risk**: remaining provider migration work is outside the
  cfd-math integration-test vector bridge layer and includes
  `nalgebra_sparse::CsrMatrix`, dense nalgebra test oracles, and source-level
  Leto storage-slice internals that should be handled in their owning
  production/provider slices.

---

## Sprint 2026-07-05: cfd-math AMG Integration Test Leto Array1

- **Resolved AMG vector bridge**:
  `crates/cfd-math/tests/amg_integration_test.rs` now constructs exact
  solutions, RHS vectors, solver outputs, AMG cycle outputs, and two-grid
  preconditioner work buffers as `leto::Array1` values directly instead of
  converting through nalgebra `DVector` and Leto storage slices.
- **Resolved preconditioner bridge**:
  AMG cycle and two-grid tests now call `Preconditioner::apply_to` directly at
  the Leto RHS/output boundary; solution-error and cycle-output assertions use
  Leto vector helper reductions.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest over AMG integration tests, broader empirical nextest over the AMG
  filter, clippy over the focused test and all cfd-math targets, static source
  audit, and format hygiene. In `D:/atlas/repos/CFDrs`, cfd-math fmt passed;
  amg_integration_test check passed; amg_integration_test nextest passed 5/5;
  amg_integration_test clippy passed; cfd-math all-target check passed;
  cfd-math all-target clippy passed; AMG-filter nextest passed 6/6; targeted
  scan found no `DVector`, nalgebra vector import, `Storage`, storage-slice
  conversion, or local preconditioner bridge residue in the migrated test.
- **Residual risk**: `amg_integration_test.rs` still uses
  `nalgebra_sparse::CsrMatrix` as the sparse storage boundary and nalgebra
  `DMatrix`/`SymmetricEigen` as the dense energy-norm oracle. The integration
  storage-slice holdout was closed by Sprint 1.96.153; production
  sparse/scalar boundaries still include transitional nalgebra providers.

---

## Sprint 2026-07-05: cfd-math Preconditioner Edge-Case Tests Leto Array1

- **Resolved preconditioner edge-case bridge**:
  `crates/cfd-math/tests/preconditioner_edge_cases.rs` now constructs
  `leto::Array1` RHS and solution vectors directly for ILU(0), ILU(k),
  repeated-application, extreme-value, and sparsity-preservation tests instead
  of converting through a local nalgebra `DVector` preconditioner bridge.
- **Resolved provider-boundary assertion path**:
  preconditioner application now calls `Preconditioner::apply_to` directly at
  the Leto boundary; the finite/bounded value assertions inspect the Leto
  output buffers without storage-slice conversion.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest over preconditioner edge-case tests, broader empirical nextest over
  preconditioner-filtered cfd-math tests, clippy over the focused test and
  cfd-math all targets, static source audit, and format hygiene. In
  `D:/atlas/repos/CFDrs`, cfd-math fmt passed; preconditioner_edge_cases check
  passed; preconditioner_edge_cases nextest passed 6/6;
  preconditioner_edge_cases clippy passed; cfd-math all-target check passed;
  cfd-math all-target clippy passed; preconditioner nextest passed 76/76;
  targeted scan found no `DVector`, nalgebra vector import, `Storage`,
  storage-slice conversion, or local preconditioner bridge residue in the
  migrated test.
- **Residual risk**: `preconditioner_edge_cases.rs` still uses the current
  `nalgebra_sparse::CsrMatrix` matrix storage boundary. Integration-test
  provider holdouts were closed by Sprints 1.96.152 and 1.96.153; production
  sparse/scalar boundaries still include `nalgebra_sparse::CsrMatrix` and
  transitional nalgebra scalar/test utilities.

---

## Sprint 2026-07-05: cfd-math Linear-Solver Test Module Leto Array1

- **Resolved source test bridge**:
  `crates/cfd-math/src/linear_solver/tests/{mod,edge_case_tests,
  adversarial_solver_tests,extended_edge_case_tests}.rs` now constructs
  `leto::Array1` RHS, solution, residual, and work buffers directly instead of
  converting through nalgebra `DVector` bridge macros.
- **Resolved scalar helper drift**:
  the touched linear-solver/sparse cone no longer uses the old
  `T::zero()`, `T::one()`, `default_epsilon()`, or `to_subset()` helper calls;
  constants and tolerances route through Eunomia trait surfaces.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest over source linear-solver tests, broader empirical nextest over the
  linear-solver filter, clippy over lib/tests and all cfd-math targets, static
  source audits, and format hygiene. In `D:/atlas/repos/CFDrs`, cfd-math fmt
  passed; cfd-math lib check passed; cfd-math all-target check passed;
  linear_solver::tests nextest passed 53/53; linear_solver nextest passed
  176/176; cfd-math lib/tests clippy passed; cfd-math all-target clippy
  passed; current targeted scans found no `DVector`, nalgebra vector import,
  `Storage`, storage-slice conversion, local solve/preconditioner bridge,
  matrix-vector `&a * &x`, or vector `.norm()` residue in
  `src/linear_solver/tests`, and no old scalar helper residue in the searched
  solver/sparse cone.
- **Residual risk**: integration-test provider holdouts were closed by Sprints
  1.96.152 and 1.96.153; production sparse/scalar boundaries still include
  `nalgebra_sparse::CsrMatrix` and transitional nalgebra scalar/test utilities.

---

## Sprint 2026-07-04: cfd-math Core Solver Tests Leto Array1

- **Resolved core solver test bridge**:
  `crates/cfd-math/tests/core_solver_tests.rs` now constructs `leto::Array1`
  RHS and solution vectors directly for BiCGSTAB, GMRES, preconditioner
  integration, and condition-number robustness tests instead of converting
  from nalgebra `DVector` through local solve/preconditioner bridge helpers.
- **Resolved residual assertions**:
  residual verification now uses `cfd_math::sparse::spmv` on Leto vectors and
  asserts explicit residual thresholds, replacing the remaining nalgebra
  matrix-vector multiplication/vector-norm checks in this integration module.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest over core solver tests, clippy over the focused test and cfd-math all
  targets, static source audit, and format hygiene. In
  `D:/atlas/repos/CFDrs`, cfd-math fmt passed; core_solver_tests check passed;
  core_solver_tests nextest passed 4/4; core_solver_tests clippy passed;
  cfd-math all-target check passed; cfd-math all-target clippy passed;
  targeted scan found no `DVector`, nalgebra vector, `Storage`, storage-slice
  conversion, local solve/preconditioner bridge, matrix-vector `&a * &x`, or
  vector `.norm()` residue in the core solver test module.
- **Residual risk**: `core_solver_tests.rs` still uses the current
  `nalgebra_sparse::CsrMatrix`/`CooMatrix` matrix storage boundary, and broader
  cfd-math integration/adversarial/preconditioner diagnostics still contain
  nalgebra `DVector` bridges.

---

## Sprint 2026-07-04: cfd-math Simple GMRES Tests Leto Array1

- **Resolved simple GMRES test bridge**:
  `crates/cfd-math/tests/simple_gmres_tests.rs` now constructs
  `leto::Array1` RHS and solution vectors directly for basic, restarted, and
  preconditioned GMRES integration tests instead of converting from nalgebra
  `DVector` through a local bridge macro.
- **Resolved residual assertions**:
  residual verification now uses `cfd_math::sparse::spmv` on Leto vectors and
  asserts explicit residual thresholds for all three tests, replacing the
  remaining nalgebra matrix-vector multiplication/vector-norm checks.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest over simple GMRES tests, clippy over the focused test and cfd-math
  all targets, static source audit, and format/diff hygiene. In
  `D:/atlas/repos/CFDrs`, cfd-math fmt passed; simple_gmres test check
  passed; simple_gmres nextest passed 3/3; simple_gmres clippy passed;
  cfd-math all-target check passed; cfd-math all-target clippy passed;
  touched-file `git diff --check` passed; targeted scan found no `DVector`,
  nalgebra vector, `Storage`, storage-slice conversion, matrix-vector
  `&a * &x`, or vector `.norm()` residue in the simple GMRES test module.
- **Residual risk**: `simple_gmres_tests.rs` still uses the current
  `nalgebra_sparse::CsrMatrix`/`CooMatrix` matrix storage boundary, and
  broader cfd-math integration/adversarial/core/preconditioner tests still
  contain nalgebra `DVector` bridge diagnostics.

---

## Sprint 2026-07-04: cfd-math Matrix-Free Tests Leto Array1

- **Resolved matrix-free test bridge**:
  `crates/cfd-math/src/linear_solver/matrix_free/tests.rs` now constructs
  `leto::Array1` RHS and solution vectors directly for CG/GMRES matrix-free
  solver tests instead of converting from nalgebra `DVector` through a local
  bridge macro.
- **Resolved mismatch assertion**:
  the operator-size mismatch test now keeps the solution buffer at the same
  Leto length as the RHS so it exercises the operator-size check and asserts
  the exact `InvalidConfiguration` message.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest over matrix-free tests, clippy over cfd-math all targets, static
  source audit, and format/diff hygiene. In `D:/atlas/repos/CFDrs`, cfd-math
  fmt passed; cfd-math all-target check passed; cfd-math matrix-free nextest
  passed 4/4; cfd-math all-target clippy passed; touched-file
  `git diff --check` passed; targeted scan found no `DVector`, nalgebra,
  `Storage`, storage-slice conversion, or local bridge macro residue in the
  matrix-free test module.
- **Residual risk**: broader cfd-math linear-solver integration/adversarial/
  core/preconditioner test diagnostics still contain nalgebra `DVector`
  bridges, and production sparse/scalar provider boundaries still include
  `nalgebra_sparse::CsrMatrix` and transitional `nalgebra::RealField`.

---

## Sprint 2026-07-04: cfd-math BiCGSTAB Leto Array1

- **Resolved BiCGSTAB vector provider**:
  `crates/cfd-math/src/linear_solver/bicgstab/mod.rs` now stores residual,
  shadow residual, search, operator-product, stabilization, preconditioned,
  and initial-product workspaces as `leto::Array1` instead of nalgebra
  `DVector`.
- **Resolved BiCGSTAB solve bridge**:
  preconditioned/unpreconditioned BiCGSTAB methods and the
  `IterativeLinearSolver` implementation now run directly on Leto
  RHS/solution arrays and call `LinearOperator::apply`/
  `Preconditioner::apply_to` without legacy nalgebra-vector bridge helpers.
- **Resolved chain fallback bridge**:
  `crates/cfd-math/src/linear_solver/chain.rs` now keeps the final
  last-resort BiCGSTAB tier on the existing Leto RHS/solution buffers instead
  of converting to nalgebra vectors.
- **Resolved duplicated vector helpers**:
  CG and BiCGSTAB now share `linear_solver/array_ops.rs` for Leto-array dot,
  norm, copy, residual, axpy, and scale-add operations; the obsolete legacy
  bridge helpers were removed from `linear_solver/traits.rs`.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest over BiCGSTAB tests, broader empirical nextest over linear-solver
  tests, empirical AMG integration coverage, clippy over cfd-math all targets,
  static source audit, and format/diff hygiene. In
  `D:/atlas/repos/CFDrs`, cfd-math fmt passed; cfd-math all-target check
  passed; focused cfd-math BiCGSTAB nextest passed 24/24; broader cfd-math
  linear-solver nextest passed 176/176; AMG integration nextest passed 5/5;
  cfd-math all-target clippy passed; touched-file `git diff --check` passed;
  targeted scans found no legacy bridge helper residue and no `DVector`/
  nalgebra vector-operation residue in the migrated BiCGSTAB/chain/traits
  source.
- **Residual risk**: cfd-math still owns the current
  `nalgebra_sparse::CsrMatrix` storage/provider boundary and transitional
  `nalgebra::RealField` scalar bounds. Remaining nalgebra `DVector` usage in
  this area is confined to other matrix-free/preconditioner/integration test
  diagnostics until those boundaries migrate.

---

## Sprint 2026-07-04: cfd-math Conjugate Gradient Leto Array1

- **Resolved CG vector provider**:
  `crates/cfd-math/src/linear_solver/conjugate_gradient/mod.rs` now stores CG
  residual, direction, preconditioned residual, operator product, and
  initial-product workspaces as `leto::Array1` instead of nalgebra `DVector`.
- **Resolved CG solve bridge**:
  preconditioned/unpreconditioned CG methods and the `IterativeLinearSolver`
  implementation now run directly on Leto RHS/solution arrays and call
  `LinearOperator::apply`/`Preconditioner::apply_to` without legacy
  nalgebra-vector bridge helpers.
- **Resolved benchmark boundary**:
  `crates/cfd-math/benches/{cg_bench,math_benchmarks}.rs` now constructs Leto
  vectors at the measured CG API instead of timing a nalgebra-vector boundary.
- **Resolved test semantics**:
  co-located CG tests assert solved-system values and exact typed error
  outcomes for dimension mismatch and max-iteration exhaustion.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest over CG tests, broader empirical nextest over linear-solver tests,
  clippy over cfd-math all targets, static source audit, and format hygiene.
  In `D:/atlas/repos/CFDrs`, cfd-math fmt passed; cfd-math all-target check
  passed; cfd-math all-target clippy passed; focused cfd-math conjugate
  nextest passed 13/13; broader cfd-math linear-solver nextest passed 176/176;
  targeted source scan found no `DVector`, legacy bridge helper calls,
  `num_traits`, Leto storage-slice conversion, nalgebra vector math calls, or
  nalgebra `copy_from` residue in the migrated CG source and CG benchmark call
  sites.
- **Residual risk after the BiCGSTAB follow-up**: the shared linear-solver
  trait family still carries the transitional `nalgebra::RealField` scalar
  bound and sparse storage remains on `nalgebra_sparse::CsrMatrix` until the
  remaining scalar-provider stack moves fully to Leto/Eunomia.

---

## Sprint 2026-07-04: cfd-math Schwarz Preconditioner Leto Array1

- **Resolved Schwarz vector provider**:
  `crates/cfd-math/src/linear_solver/preconditioners/schwarz.rs` now exposes
  additive and multiplicative Schwarz application over `leto::Array1` instead
  of nalgebra `DVector`.
- **Resolved local RHS bridge**:
  Schwarz local RHS extraction now writes directly into Leto arrays, local ILU
  solves consume those arrays, and the default additive
  `Preconditioner::apply_to` path copies from a Leto result without nalgebra
  vector conversion.
- **Resolved mismatch validation**:
  residual/output length mismatches now return typed configuration errors, and
  tests assert exact identity-preconditioned values plus exact mismatch
  messages.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest over Schwarz tests, broader empirical nextest over preconditioner
  tests, clippy over cfd-math all targets, static source audit, and format
  hygiene. In `D:/atlas/repos/CFDrs`, cfd-math fmt passed; cfd-math
  all-target check passed; cfd-math all-target clippy passed; focused
  cfd-math Schwarz nextest passed 3/3; broader cfd-math preconditioner
  nextest passed 76/76; targeted source scan found no `DVector`, `Storage`,
  `num_traits`, `FromPrimitive`, local-RHS conversion bridge, `DMatrix`, or
  `ndarray` residue in `schwarz.rs`.
- **Residual risk**: Schwarz still stores and constructs local sparse matrices
  through the current `nalgebra_sparse::CsrMatrix` boundary and remains
  constrained by the global `Preconditioner<T>` nalgebra scalar bound until the
  wider sparse/preconditioner stack migrates fully to Leto storage and Eunomia
  scalar traits.

---

## Sprint 2026-07-04: cfd-math ILU Triangular Solve Leto Array1

- **Resolved ILU vector provider**:
  `crates/cfd-math/src/linear_solver/preconditioners/ilu/{types,triangular}.rs`
  now performs ILU forward and backward substitution directly on
  `leto::Array1` residual, intermediate, and solution buffers instead of
  nalgebra `DVector`.
- **Resolved ILU apply bridge**:
  `IncompleteLU::apply_to` now validates the caller-provided Leto
  residual/output buffers and copies the Leto triangular-solve result directly
  into the output array without vector conversion.
- **Resolved mismatch and scalar dispatch**:
  residual/output length mismatches now return typed configuration errors, and
  the backward substitution diagonal identity dispatches through Eunomia
  `NumericElement`.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest over ILU tests, broader empirical nextest over preconditioner tests,
  clippy over cfd-math all targets, static source audit, and format hygiene.
  In `D:/atlas/repos/CFDrs`, cfd-math fmt passed; cfd-math all-target check
  passed; cfd-math all-target clippy passed; focused cfd-math ILU nextest
  passed 21/21; broader cfd-math preconditioner nextest passed 74/74; targeted
  source scan found no `DVector`, `DMatrix`, `Storage`, old scalar identities,
  fallback wording, `component_mul`, `rows_mut`, row-view residue,
  `num_traits`, or `num_complex` in `ilu/types.rs` or `ilu/triangular.rs`.
- **Residual risk**: IncompleteLU still stores the current
  `nalgebra_sparse::CsrMatrix` LU factor boundary and remains constrained by
  the global `Preconditioner<T>` nalgebra scalar bound until the wider
  sparse/preconditioner stack migrates fully to Leto and Eunomia.

---

## Sprint 2026-07-04: cfd-math Deflation Preconditioner Leto Array1

- **Resolved deflation eigenvector provider**:
  `crates/cfd-math/src/linear_solver/preconditioners/deflation.rs` now stores
  deflation eigenvectors as `leto::Array1` instead of nalgebra `DVector`.
- **Resolved deflation apply bridge**:
  `DeflationPreconditioner::apply_to` now keeps the base-preconditioner result
  and projection correction in Leto arrays without constructing a nalgebra
  deflated work vector.
- **Resolved mismatch and division validation**:
  output/eigenvector length mismatches now return typed configuration errors,
  added eigenpairs use the Leto vector boundary, and zero eigenvalues are
  rejected before projection division.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest over deflation tests, broader empirical nextest over preconditioner
  tests, clippy over cfd-math all targets, static source audit, and format
  hygiene. In `D:/atlas/repos/CFDrs`, cfd-math fmt passed; cfd-math
  all-target check passed; cfd-math all-target clippy passed; focused
  cfd-math deflation nextest passed 3/3; broader cfd-math preconditioner
  nextest passed 73/73; targeted source scan found no `DVector`, `DMatrix`,
  `Storage`, old scalar identities, fallback wording, `component_mul`,
  `rows_mut`, row-view residue, `num_traits`, or `num_complex` in
  `deflation.rs`.
- **Residual risk**: Deflation still wraps the base preconditioner behind the
  existing `Box<dyn Preconditioner<T>>`, and cfd-math's global
  `Preconditioner<T>` trait still requires the transitional nalgebra
  `RealField` scalar bound.

---

## Sprint 2026-07-04: cfd-math Basic Preconditioners Leto Array1

- **Resolved basic-preconditioner vector provider**:
  `crates/cfd-math/src/linear_solver/preconditioners/basic.rs` now applies
  Identity, Jacobi, and SOR directly over `leto::Array1` residual/output
  buffers instead of nalgebra vector state.
- **Resolved Jacobi state provider**:
  Jacobi inverse diagonal storage now uses `leto::Array1`, and diagonal
  extraction uses the Leto storage slice boundary after the sparse diagonal is
  produced.
- **Resolved mismatch and scalar dispatch**:
  Identity output mismatches and Jacobi/SOR residual/output mismatches now
  return typed configuration errors; Jacobi/SOR scalar identities and diagonal
  tolerance dispatch through Eunomia traits.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest over the new mismatch regression, broader empirical nextest over
  preconditioner tests, clippy over cfd-math all targets, static source audit,
  and format hygiene. In `D:/atlas/repos/CFDrs`, cfd-math fmt passed;
  cfd-math all-target check passed; cfd-math all-target clippy passed; focused
  cfd-math basic mismatch nextest passed 1/1; broader cfd-math preconditioner
  nextest passed 70/70; targeted source scan found no `DVector`, `DMatrix`,
  old scalar identities, fallback wording, `component_mul`, `rows_mut`, or
  row-view residue in `basic.rs`.
- **Residual risk**: Jacobi and SOR still store/consume the current
  `nalgebra_sparse::CsrMatrix` sparse matrix boundary and remain constrained by
  the global `Preconditioner<T>` nalgebra scalar bound until the wider
  preconditioner trait family migrates fully to Eunomia and Leto sparse
  storage.

---

## Sprint 2026-07-04: cfd-math IncompleteCholesky Preconditioner Leto Array1

- **Resolved Cholesky vector provider**:
  `crates/cfd-math/src/linear_solver/preconditioners/cholesky.rs` now performs
  forward and backward triangular substitution directly on `leto::Array1`
  residual, intermediate, and solution buffers instead of nalgebra `DVector`.
- **Resolved Cholesky apply bridge**:
  `IncompleteCholesky::apply_to` now copies the Leto substitution result into
  the caller-provided Leto output buffer without constructing `DVector`
  workspaces.
- **Resolved mismatch and scalar dispatch**:
  residual/output length mismatches now return typed configuration errors
  before substitution, and IC(0) factorization square-root dispatch is routed
  through Eunomia `NumericElement`.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest over Cholesky tests, clippy over cfd-math all targets, static source
  audit, and format hygiene. In `D:/atlas/repos/CFDrs`, cfd-math fmt passed;
  cfd-math all-target check passed; cfd-math all-target clippy passed;
  focused cfd-math Cholesky nextest passed 5/5; targeted scans found no
  `DVector`, `DMatrix`, `Storage`, old scalar identities, silent fallback
  wording, ambiguous `.sqrt()`, or residual `DVector` workspace construction
  in `cholesky.rs`.
- **Residual risk**: resolved for IncompleteCholesky by Sprint 1.96.166, which
  moved its factor boundary to Leto CSR. Remaining sparse-provider risk is
  Schwarz, direct solver, and the shared solver matrix boundary.

---

## Sprint 2026-07-04: cfd-math SSOR Preconditioner Leto Array1

- **Resolved SSOR vector provider**:
  `crates/cfd-math/src/linear_solver/preconditioners/ssor.rs` now performs
  forward and backward SOR sweeps directly on `leto::Array1` residual and
  solution buffers instead of nalgebra `DVector`.
- **Resolved SSOR apply bridge**:
  `SSOR::apply_to` now initializes and mutates the caller-provided Leto output
  buffer directly. The previous residual and solution `DVector` allocations
  are removed.
- **Resolved mismatch handling**:
  residual/output length mismatches now return typed configuration errors
  before sweep execution.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest over SSOR tests, clippy over cfd-math all targets, static source
  audit, and format hygiene. In `D:/atlas/repos/CFDrs`, cfd-math fmt passed;
  cfd-math all-target check passed; cfd-math all-target clippy passed;
  focused cfd-math SSOR nextest passed 5/5; targeted scans found no
  `DVector`, `DMatrix`, `Storage`, nalgebra row-view operations,
  `component_mul`, old scalar identities, silent fallback wording, or clone
  fallback in `ssor.rs`.
- **Residual risk**: SSOR still stores the current
  `nalgebra_sparse::CsrMatrix` sparse matrix boundary and remains constrained
  by the global `Preconditioner<T>` nalgebra scalar bound until the wider
  preconditioner trait family migrates fully to Eunomia.

---

## Sprint 2026-07-04: cfd-math Block/SIMPLE Preconditioner Leto Array1

- **Resolved block-preconditioner vector provider**:
  `crates/cfd-math/src/linear_solver/block_preconditioner.rs` now stores
  Jacobi diagonal inverses and SIMPLE Schur diagonal inverses as
  `leto::Array1` instead of nalgebra `DVector`.
- **Resolved preconditioner apply bridge**:
  Direct block/SIMPLE `apply` methods consume Leto arrays and
  `Preconditioner::apply_to` writes from those Leto results directly, without
  constructing residual/result `DVector` bridges.
- **Resolved silent mismatch fallback**:
  Block/SIMPLE/diagonal preconditioner length mismatches now return typed
  configuration errors instead of cloning the mismatched input vector.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest over block-preconditioner tests, clippy over cfd-math all targets,
  static source audit, and format hygiene. In `D:/atlas/repos/CFDrs`,
  cfd-math fmt passed; cfd-math all-target check passed; cfd-math all-target
  clippy passed; focused cfd-math block-preconditioner nextest passed 4/4;
  targeted scans found no `DVector`, `DMatrix`, nalgebra import, nalgebra
  row-view operations, `component_mul`, old scalar identities, silent fallback
  wording, or clone fallback in `block_preconditioner.rs`.
- **Residual risk**: the crate-level `Preconditioner<T>` trait still requires
  the transitional `nalgebra::RealField` scalar bound, so the trait impls keep
  that bound until the broader linear-solver trait family moves fully to
  Eunomia. Sparse storage still aliases `nalgebra_sparse::CsrMatrix`.

---

## Sprint 2026-07-04: cfd-math GMRES Leto Workspace

- **Resolved GMRES vector provider**:
  `crates/cfd-math/src/linear_solver/gmres/{solver,arnoldi}.rs` now stores
  GMRES work, preconditioned work, Arnoldi basis-column extraction, and
  residual verification vectors as `leto::Array1` instead of nalgebra
  `DVector`.
- **Resolved solver-chain GMRES bridge**:
  `crates/cfd-math/src/linear_solver/chain.rs` now passes Leto arrays directly
  through GMRES+AMG, GMRES+block, unpreconditioned GMRES, and GMRES+ILU tiers.
  The chain converts to nalgebra only for the still-legacy final BiCGSTAB
  fallback.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest over GMRES and the cfd-2d momentum consumer, clippy over cfd-math
  all targets, static residue scan, and format hygiene. In
  `D:/atlas/repos/CFDrs`, cfd-math lib/tests/all-targets checks passed;
  cfd-math fmt passed; cfd-math all-target clippy passed; focused cfd-math
  GMRES nextest passed 21/21; cfd-2d no-default lib check passed; focused
  cfd-2d momentum nextest passed 53/53; targeted scans found no `DVector` or
  legacy operator/preconditioner bridge calls under `linear_solver/gmres`.
- **Residual risk**: `cfd-math` still has nalgebra scalar bounds, CG and
  BiCGSTAB nalgebra workspaces, `LinearSolverChain`'s final BiCGSTAB bridge,
  and `nalgebra_sparse::CsrMatrix` storage pending later Leto/Eunomia slices.

---

## Sprint 2026-07-04: cfd-2d Momentum Leto Array1 and Direct nalgebra Removal

- **Resolved momentum vector provider**:
  `crates/cfd-2d/src/physics/momentum/{solver,solve,boundary/**}.rs` now uses
  `leto::Array1` for RHS and solution buffers instead of nalgebra `DVector`.
  Momentum boundary helpers mutate Leto arrays directly.
- **Resolved cfd-2d local solve bridge**:
  Momentum solve paths now call `IterativeLinearSolver::solve` and
  `DirectSparseSolver::solve` directly. The obsolete
  `crates/cfd-2d/src/linear_solver_bridge.rs` module and its `lib.rs` entry
  were removed.
- **Resolved direct cfd-2d nalgebra ownership**:
  `crates/cfd-2d/src/scalar.rs` no longer requires `nalgebra::RealField`, and
  `crates/cfd-2d/Cargo.toml` no longer declares direct `nalgebra` or
  `nalgebra-sparse` dependencies.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest over momentum tests, clippy over the affected crate targets, static
  source/manifest residue scan, dependency graph audit, and format hygiene. In
  `D:/atlas/repos/CFDrs`, cfd-2d fmt passed; no-default lib check passed;
  no-default all-target clippy passed; focused cfd-2d momentum nextest passed
  53/53; targeted cfd-2d source/manifest scans found no `DVector`, `DMatrix`,
  nalgebra, nalgebra-sparse, or `linear_solver_bridge`; `cargo tree -p cfd-2d
  -e normal -i nalgebra` and `-i nalgebra-sparse` show remaining graph edges
  are transitive through upstream crates.
- **Residual risk**: cfd-2d still resolves nalgebra/nalgebra-sparse
  transitively through upstream owners including cfd-1d, cfd-core, cfd-math,
  cfd-schematics, and Gaia while those crates continue their provider
  migrations.

---

## Sprint 2026-07-04: cfd-2d Pressure-Velocity Leto Array1

- **Resolved pressure-velocity vector provider**:
  `crates/cfd-2d/src/pressure_velocity/{pressure,correction,faces}.rs` now
  caches pressure-correction RHS and solution buffers as `leto::Array1`
  instead of nalgebra `DVector`.
- **Resolved pressure-velocity solve bridge**:
  Pressure-correction dispatch now calls the Leto-native
  `IterativeLinearSolver::solve` boundary and `DirectSparseSolver::solve`
  fallback directly instead of routing through the local nalgebra conversion
  bridge. RHS diagnostics now use `leto_ops::norm_l2`.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest over pressure-velocity tests, clippy over the affected crate targets,
  static residue scan, and format hygiene. In `D:/atlas/repos/CFDrs`, cfd-2d
  fmt passed; no-default lib check passed; no-default all-target clippy passed;
  focused cfd-2d pressure-velocity nextest passed 16/16; targeted scans found
  no `DVector`, nalgebra, or `linear_solver_bridge` residue under
  `crates/cfd-2d/src/pressure_velocity`.
- **Residual risk**: direct cfd-2d source/manifest nalgebra ownership is now
  removed; nalgebra/nalgebra-sparse still resolve transitively through upstream
  owners while the broader Atlas migration continues.

---

## Sprint 2026-07-04: cfd-2d SIMPLE Pressure Correction Leto Array1

- **Resolved SIMPLE pressure-correction vector provider**:
  `crates/cfd-2d/src/solvers/simple/algorithm.rs` now stores the
  pressure-correction RHS and `p_prime` buffers as `leto::Array1` instead of
  nalgebra `DVector`.
- **Resolved SIMPLE solve bridge**:
  `crates/cfd-2d/src/solvers/simple/pressure.rs` now calls the Leto-native
  `IterativeLinearSolver::solve` boundary directly for the SIMPLE
  pressure-correction solve instead of routing through
  `linear_solver_bridge::iterative_solve`.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest over SIMPLE tests, clippy over the affected crate targets, static
  residue scan, and format hygiene. In `D:/atlas/repos/CFDrs`, cfd-2d fmt
  passed; no-default lib check passed; no-default all-target clippy passed;
  focused cfd-2d SIMPLE nextest passed 19/19; targeted scans found no
  `DVector`, nalgebra, or `linear_solver_bridge` residue under
  `crates/cfd-2d/src/solvers/simple`.
- **Residual risk**: direct cfd-2d source/manifest nalgebra ownership is now
  removed; nalgebra/nalgebra-sparse still resolve transitively through upstream
  owners while the broader Atlas migration continues.

---

## Sprint 2026-07-04: cfd-2d FDM Leto Array1

- **Resolved FDM vector provider**:
  `crates/cfd-2d/src/solvers/fdm/linear_solver.rs`,
  `crates/cfd-2d/src/solvers/fdm/poisson.rs`, and
  `crates/cfd-2d/src/solvers/fdm/advection_diffusion.rs` now use
  `leto::Array1` for RHS and solution vectors instead of nalgebra `DVector`.
- **Resolved FDM solver boundary**:
  `solve_gauss_seidel` now accepts a Leto RHS and returns a Leto solution.
  Poisson and advection-diffusion stencil builders mutate Leto RHS arrays
  directly before solving.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest over migrated FDM tests, clippy over the affected crate targets,
  static residue scan, and format hygiene. In `D:/atlas/repos/CFDrs`, cfd-2d
  fmt passed; no-default lib check passed; no-default all-target clippy passed;
  focused cfd-2d FDM nextest passed 2/2; targeted scans found no `DVector` or
  nalgebra residue under `crates/cfd-2d/src/solvers/fdm`.
- **Residual risk**: direct cfd-2d source/manifest nalgebra ownership is now
  removed; nalgebra/nalgebra-sparse still resolve transitively through upstream
  owners while the broader Atlas migration continues.

---

## Sprint 2026-07-04: cfd-2d Time Integration Leto Array1

- **Resolved time-state vector provider**:
  `crates/cfd-2d/src/schemes/time/**` now uses `leto::Array1` for
  time-integration state vectors through `StateVector<T>` instead of nalgebra
  `DVector`.
- **Resolved time-kernel vector operations**:
  Explicit, implicit, multistep, adaptive-controller, and adaptive-integrator
  paths now use Leto vector arithmetic, indexing, and shape access. Fixed-point
  convergence checks use `leto_ops::norm_l2`.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest over migrated time tests, clippy over the affected crate library,
  static residue scans, and format hygiene. In `D:/atlas/repos/CFDrs`,
  cfd-2d fmt passed; no-default lib check passed; no-default all-target clippy
  passed; focused cfd-2d time nextest passed 29/29; targeted scans found no
  `DVector` or nalgebra residue under `crates/cfd-2d/src/schemes/time`.
- **Residual risk**: direct cfd-2d source/manifest nalgebra ownership is now
  removed; nalgebra/nalgebra-sparse still resolve transitively through upstream
  owners while the broader Atlas migration continues.

---

## Sprint 2026-07-04: cfd-2d DMatrix Residue Leto Array2

- **Resolved immersed-boundary matrix provider**:
  `crates/cfd-2d/src/physics/immersed_boundary.rs` now accepts
  `leto::Array2<f64>` force and velocity matrices with explicit `[nx * ny, 2]`
  shape semantics instead of nalgebra `DMatrix`.
- **Resolved scheme grid storage provider**:
  `crates/cfd-2d/src/schemes/grid.rs` now stores `Grid2D<T>::data` as
  `leto::Array2<T>`. Scheme callers in `upwind.rs`, `tvd/muscl.rs`, TVD tests,
  and WENO tests now use Leto `shape()` and `[[i, j]]` indexing.
- **Resolved example boundary**:
  `crates/cfd-2d/examples/blood_venturi.rs` now passes Leto arrays through the
  immersed-boundary coupling path instead of constructing `DMatrix` buffers.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest over the migrated scheme and immersed-boundary surface, clippy over
  the affected crate library, static residue scans, and format hygiene. In
  `D:/atlas/repos/CFDrs`, cfd-2d fmt passed; no-default lib check and clippy
  passed; no-default `blood_venturi` example check passed; focused cfd-2d
  nextest passed 60/60; targeted scans found no cfd-2d source/example/test/bench
  `DMatrix` residue and no nalgebra-style `Grid2D.data` tuple access.
- **Residual risk**: cfd-2d still retains nalgebra in non-DMatrix
  `DVector`/sparse linear-system boundaries and other provider seams while the
  broader Atlas migration continues.

---

## Sprint 2026-07-04: cfd-1d Vascular Complex Math Eunomia

- **Resolved vascular complex provider**:
  `crates/cfd-1d/src/physics/vascular/bessel.rs` and
  `crates/cfd-1d/src/physics/vascular/womersley/profile.rs` now use
  `eunomia::Complex` for the Bessel recurrence and Womersley analytical
  profile instead of nalgebra complex numbers.
- **Resolved nalgebra ComplexField dependency**:
  The Bessel convergence stop condition now uses Eunomia complex `norm()`, so
  the vascular path no longer imports `nalgebra::ComplexField` only for
  `modulus()`.
- **Evidence tier**: compile-time provider integration, focused empirical
  Bessel/Womersley nextest, clippy over the affected crate library, static
  residue scan, and format hygiene. In `D:/atlas/repos/CFDrs`,
  `rustup run nightly cargo check -p cfd-1d --no-default-features --lib`,
  `cargo fmt -p cfd-1d --check`, `cargo clippy -p cfd-1d
  --no-default-features --lib -- -D warnings`, and `cargo nextest run -p
  cfd-1d --no-default-features bessel womersley --status-level fail` passed
  26/26. A targeted vascular residue scan found no remaining nalgebra complex
  imports or `ComplexField::modulus()` calls.
- **Residual risk**: cfd-1d still depends on nalgebra for network sparse/dense
  linear-system storage and the transitional `Cfd1dScalar` seam while those
  boundaries continue migrating to Leto/Eunomia.

---

## Sprint 2026-07-04: cfd-math LinearOperator Trait Leto Vectors

- **Resolved public operator vector API**:
  `crates/cfd-math/src/linear_solver/traits.rs` now exposes
  `LinearOperator::apply` and `apply_transpose` with `leto::Array1<T>` input
  and output buffers instead of nalgebra `DVector`.
- **Resolved operator implementations**:
  sparse CSR, identity/scaled, Poisson, momentum, energy, and GPU operator
  adapters implement the Leto vector boundary. The Laplacian CPU benchmark now
  calls the Leto operator API.
- **Confined residual bridges**:
  CG, BiCGSTAB, GMRES, and Arnoldi call an internal
  `apply_operator_to_legacy` bridge only because their current workspaces still
  use nalgebra vectors.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest, clippy over the producer crate, static residue scan, and format
  hygiene. In `D:/atlas/repos/CFDrs`, `rustup run nightly cargo check -p
  cfd-math --no-default-features --all-targets`, `cargo fmt -p cfd-math
  --check`, `cargo clippy -p cfd-math --no-default-features --all-targets --
  -D warnings`, and `cargo nextest run -p cfd-math --no-default-features
  conjugate_gradient bicgstab gmres operator poisson momentum matrix_free
  core_solver simple_gmres --status-level fail` passed 80/80. A targeted
  DVector operator-signature residue scan found no old public signature.
- **Residual risk after the CG/BiCGSTAB follow-up**: nalgebra sparse storage
  and some preconditioner internals still retain local `DVector`/`CsrMatrix`
  conversion bridges. cfd-validation numerical result and error metric storage
  still use `DVector`. Broad cfd-validation nextest still fails in the
  existing venturi cross-fidelity convergence tests.

---

## Sprint 2026-07-04: cfd-math Preconditioner Trait Leto Vectors

- **Resolved public preconditioner vector API**:
  `crates/cfd-math/src/linear_solver/traits.rs` now exposes
  `Preconditioner::apply_to` with `leto::Array1<T>` residual and output
  buffers instead of nalgebra `DVector`.
- **Resolved producer implementations and tests**:
  basic, ILU, incomplete Cholesky, SSOR, Schwarz, deflation, block, and AMG
  preconditioner paths implement or call the Leto boundary; legacy nalgebra
  work vectors are confined to local conversion bridges where matrix internals
  still require them.
- **Resolved downstream caller**:
  `crates/cfd-1d/src/solver/core/linear_system.rs` now implements
  `DiagJacobi` as a Leto-array preconditioner for network iterative solves.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest, clippy over producer and affected consumer paths, static residue
  scan, and format hygiene. In `D:/atlas/repos/CFDrs`,
  `rustup run nightly cargo check -p cfd-math --no-default-features
  --all-targets`, `cargo check -p cfd-1d --no-default-features --lib`,
  `cargo fmt -p cfd-math -p cfd-1d --check`, `cargo clippy -p cfd-math
  --no-default-features --all-targets -- -D warnings`, `cargo clippy -p
  cfd-1d --no-default-features --lib -- -D warnings`, and `cargo nextest run
  -p cfd-math --no-default-features conjugate_gradient bicgstab gmres
  preconditioner ilu matrix_free core_solver simple_gmres --status-level fail`
  passed 131/131. A targeted `apply_to` DVector-signature residue scan found
  no old public signature.
- **Residual risk**: `LinearOperator::apply`, nalgebra sparse preconditioner
  internals, and iterative solver workspaces still retain
  `DVector`/`CsrMatrix` conversion bridges. cfd-validation numerical result
  and error metric storage still use `DVector`. Broad cfd-validation nextest
  still fails in the existing venturi cross-fidelity convergence tests.

---

## Sprint 2026-07-04: cfd-math Iterative Solver Trait Leto Vectors

- **Resolved public iterative solver vector API**:
  `crates/cfd-math/src/linear_solver/traits.rs` now exposes
  `IterativeLinearSolver::solve` with `leto::Array1<T>` RHS and mutable
  solution buffers instead of nalgebra `DVector`.
- **Resolved solver implementations and tests**:
  CG, BiCGSTAB, and GMRES implement the Leto public boundary and confine
  nalgebra vector conversion to solver-local workspace bridges. cfd-math
  solver tests call the public boundary with Leto arrays.
- **Resolved downstream callers**:
  cfd-1d network iterative solves, cfd-2d momentum/pressure iterative solves,
  and cfd-3d FEM projection iterative solves now convert through local
  Leto bridge modules instead of passing DVector buffers into the public
  iterative trait.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest, clippy over producer and affected consumers, static residue scan,
  and diff hygiene. In `D:/atlas/repos/CFDrs`, cfd-math fmt/check/clippy
  passed; focused `cargo nextest run -p cfd-math --no-default-features
  conjugate_gradient bicgstab gmres matrix_free core_solver simple_gmres
  --status-level fail` passed 61/61; cfd-1d/cfd-2d/cfd-3d focused checks and
  clippy passed; cfd-validation all-target clippy passed; a targeted
  DVector-call-site residue scan and `git diff --check` passed.
- **Residual risk**: `LinearOperator::apply`, `Preconditioner::apply_to`,
  concrete preconditioners, and internal iterative workspaces still expose
  nalgebra `DVector`/`CsrMatrix`. cfd-validation numerical solver result and
  error metric storage still use `DVector`. Broad cfd-validation nextest still
  fails in the existing venturi cross-fidelity convergence tests.

---

## Sprint 2026-07-04: cfd-math LinearSolver Trait Leto Vectors

- **Resolved public solver trait vector API**:
  `crates/cfd-math/src/linear_solver/traits.rs` now exposes
  `LinearSolver::solve_system` with `leto::Array1<T>` RHS, optional
  initial-guess, and result vectors instead of nalgebra `DVector`.
- **Resolved iterative solver implementations**:
  CG, BiCGSTAB, and GMRES now implement the Leto public boundary and confine
  nalgebra `DVector` conversion to solver-local workspace bridges.
- **Resolved validation caller boundary**:
  `crates/cfd-validation/src/numerical/linear_solver.rs` now calls
  `solve_system` with Leto arrays and converts back to nalgebra vectors only
  for its existing analytical/error-metric storage.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest, clippy over affected producer and consumer paths, static
  signature/residue scans, and diff hygiene. In `D:/atlas/repos/CFDrs`,
  `rustup run nightly cargo fmt -p cfd-math -p cfd-validation --check`,
  `cargo check -p cfd-math --no-default-features --lib`, `cargo check -p
  cfd-validation --no-default-features --lib`, `cargo clippy -p cfd-math
  --no-default-features --all-targets -- -D warnings`, `cargo clippy -p
  cfd-validation --no-default-features --lib -- -D warnings`, `cargo clippy
  -p cfd-validation --no-default-features --all-targets -- -D warnings`,
  `cargo clippy -p cfd-2d --no-default-features --all-targets -- -D
  warnings`, `cargo nextest run -p cfd-math --no-default-features
  conjugate_gradient bicgstab gmres --status-level fail` (58/58), targeted
  `solve_system` signature/residue scans, and `git diff --check` passed.
- **Residual risk**: `IterativeLinearSolver`, `LinearOperator`,
  `Preconditioner`, preconditioners, and internal iterative workspaces still
  expose nalgebra `DVector`/`CsrMatrix`. cfd-validation numerical solver
  result and error metric storage still use `DVector`. Broad cfd-validation
  nextest still fails in the existing venturi cross-fidelity convergence
  tests recorded in the previous sprint.

---

## Sprint 2026-07-04: cfd-validation Leto SpMV and Scalar Bounds

- **Resolved validation SpMV benchmark callers**:
  `crates/cfd-validation/src/benchmarking/{memory.rs,performance.rs}` now pass
  `leto::Array1` vectors into the public `cfd_math::sparse::spmv` wrappers.
- **Resolved validation scalar-bound propagation**:
  `crates/cfd-validation/src/numerical/linear_solver.rs` and
  `crates/cfd-validation/src/literature/blood_flow_1d.rs` now use the
  crate-local `ValidationScalar` seam so sparse linear-solver validation and
  cfd-1d network literature cases satisfy the Leto/Eunomia provider contracts
  reached through cfd-math and cfd-1d.
- **Evidence tier**: compile-time provider integration, clippy over the
  affected all-target consumer path, focused empirical nextest, static residue
  scan, and diff hygiene. In `D:/atlas/repos/CFDrs`, `rustup run nightly cargo
  fmt -p cfd-validation --check`, `cargo check -p cfd-validation
  --no-default-features --lib`, `cargo clippy -p cfd-validation
  --no-default-features --lib -- -D warnings`, `cargo clippy -p
  cfd-validation --no-default-features --all-targets -- -D warnings`, `cargo
  clippy -p cfd-2d --no-default-features --all-targets -- -D warnings`, a
  targeted cfd-validation SpMV DVector-residue scan, `cargo nextest run -p
  cfd-validation --no-default-features benchmark --status-level fail` (40/40),
  and `git diff --check` passed.
- **Residual risk**: broad `cargo nextest run -p cfd-validation
  --no-default-features --status-level fail` compiled but failed in
  `numerical::venturi_cross_fidelity` tests
  `microventuri_35um_case_produces_converged_informative_2d_result` and
  `option2_selected_45um_geometry_routes_to_fallback_and_converges`; both
  report non-converged 2D fallback results and remain the next validation
  baseline gap. The broader cfd-math solver/preconditioner trait family still
  exposes nalgebra `DVector` and `nalgebra_sparse::CsrMatrix`.

---

## Sprint 2026-07-04: Solver Chain and FEM Consumer Leto Vector Boundary

- **Resolved chain vector API**:
  `crates/cfd-math/src/linear_solver/chain.rs` now exposes
  `LinearSolverChain::solve` and `solve_with_guess` with `leto::Array1<T>`
  RHS/result vectors instead of nalgebra `DVector`.
- **Resolved 2D direct-fallback consumers**:
  `crates/cfd-2d/src/linear_solver_bridge.rs` is the single cfd-2d conversion
  boundary from nalgebra work vectors into the Leto-backed
  `DirectSparseSolver`; momentum and pressure fallback paths, plus the
  momentum regression test, route through it.
- **Resolved 3D FEM assembly consumers**:
  `crates/cfd-3d/src/fem/leto_bridge.rs` is the single FEM conversion boundary
  for nalgebra work vectors crossing `SparseMatrixBuilder::build_with_rhs` and
  `LinearSolverChain`. `FemSolver` and `ProjectionSolver` now use the
  crate-level `Cfd3dScalar` seam, which carries the Leto real-scalar provider
  bound.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest, clippy on touched library surfaces, static residue scan, and diff
  hygiene. In `D:/atlas/repos/CFDrs`, `rustup run nightly cargo fmt -p
  cfd-math -p cfd-1d -p cfd-2d -p cfd-3d --check`, `cargo check -p cfd-math
  --no-default-features --lib`, `cargo check -p cfd-1d --no-default-features
  --lib`, `cargo check -p cfd-2d --no-default-features --lib`, `cargo check
  -p cfd-3d --no-default-features --lib`, `cargo nextest run -p cfd-math
  --no-default-features chain direct_solver core_solver simple_gmres
  --status-level fail` (4/4), `cargo clippy -p cfd-math
  --no-default-features --all-targets -- -D warnings`, `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings`, and `cargo clippy -p
  cfd-3d --no-default-features --lib -- -D warnings` passed.
- **Residual risk**: `cargo clippy -p cfd-2d --no-default-features
  --all-targets -- -D warnings` is blocked in `cfd-validation`, not cfd-2d:
  validation benchmark modules still pass nalgebra `DVector` to public
  `cfd_math::sparse::spmv`, and generic validation/1D literature paths need
  the Leto scalar bound propagated after the cfd-1d network solver seam moved.
  The broader iterative solver/preconditioner traits still expose nalgebra
  `DVector` until that trait family moves to Leto arrays.

---

## Sprint 2026-07-04: cfd-math Direct Solver Leto Vectors

- **Resolved direct solver vector API**:
  `crates/cfd-math/src/linear_solver/direct_solver.rs` now exposes
  `DirectSparseSolver::solve` with `leto::Array1` RHS and result vectors
  instead of nalgebra `DVector`.
- **Resolved fallback bridge**:
  `crates/cfd-math/src/linear_solver/dense_bridge.rs` no longer carries the
  obsolete DVector dense-solve wrapper; the direct solver calls the Leto-array
  dense fallback directly.
- **Confined residual bridge**:
  `crates/cfd-math/src/linear_solver/chain.rs` converts between its current
  nalgebra `DVector` boundary and the Leto direct solver only at the direct
  tier.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest, all-target clippy, rustdoc/doctest, static provider-residue scan,
  and diff hygiene. In `D:/atlas/repos/CFDrs`, `rustup run nightly cargo fmt
  -p cfd-math --check`, `cargo check -p cfd-math --no-default-features --lib`,
  `cargo nextest run -p cfd-math --no-default-features direct_solver chain
  core_solver simple_gmres --status-level fail` (4/4), `cargo clippy -p
  cfd-math --no-default-features --all-targets -- -D warnings`, `cargo doc -p
  cfd-math --no-default-features --no-deps`, `cargo test --doc -p cfd-math
  --no-default-features` (3 passed, 3 ignored), targeted direct-solver
  DVector-signature residue scan, and `git diff --check` passed.
- **Residual risk**: `LinearSolverChain`, iterative solver traits,
  preconditioners, and sparse storage still expose nalgebra `DVector` and
  `nalgebra_sparse::CsrMatrix` boundaries.

---

## Sprint 2026-07-04: cfd-math Sparse Builder Leto RHS

- **Resolved builder RHS API**:
  `crates/cfd-math/src/sparse/builder.rs` now exposes
  `SparseMatrixBuilder::build_with_rhs` with a `leto::Array1` RHS instead of
  nalgebra `DVector`.
- **Resolved Dirichlet assembly path**: Column elimination mutates the Leto RHS
  directly, and `build()` no longer fabricates a dummy nalgebra vector for the
  no-RHS path.
- **Resolved call sites**:
  `crates/cfd-math/src/linear_solver/{direct_solver.rs,block_preconditioner.rs}`
  now use Leto RHS values for matrix assembly while keeping nalgebra `DVector`
  only where the current solver/preconditioner traits still require it.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest, all-target clippy, static provider-residue scan, and diff hygiene.
  In `D:/atlas/repos/CFDrs`, `rustup run nightly cargo fmt -p cfd-math
  --check`, `cargo check -p cfd-math --no-default-features --lib`, `cargo
  nextest run -p cfd-math --no-default-features sparse direct_solver
  block_preconditioner --status-level fail` (25/25), `cargo clippy -p
  cfd-math --no-default-features --all-targets -- -D warnings`, targeted
  `build_with_rhs` residue scan, and `git diff --check` passed.
- **Residual risk**: The public sparse storage alias and the
  linear-solver/preconditioner trait family still expose
  `nalgebra_sparse::CsrMatrix` and nalgebra `DVector`.

---

## Sprint 2026-07-04: cfd-math Public SpMV Leto Vector Boundary

- **Resolved public SpMV vector API**:
  `crates/cfd-math/src/sparse/operations.rs` now exposes `spmv`,
  `spmv_parallel`, and `try_spmv` with `leto::Array1` input/output vectors.
- **Resolved public call sites**: Sparse tests, GMRES/AMG integration tests,
  interpolation quality checks, and `crates/cfd-math/benches/spmv_bench.rs`
  now call the public sparse SpMV API with Leto arrays.
- **Confined residual bridge**: The nalgebra `DVector` SpMV path is now
  private (`try_spmv_dvector`) and exists only for the current
  `LinearOperator for CsrMatrix` implementation.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest, all-target clippy, static provider-residue scan, and diff hygiene.
  In `D:/atlas/repos/CFDrs`, `rustup run nightly cargo fmt -p cfd-math
  --check`, `cargo check -p cfd-math --no-default-features --lib`, `cargo
  nextest run -p cfd-math --no-default-features sparse spmv interpolation amg
  simple_gmres core_solver --status-level fail` (40/40), `cargo clippy -p
  cfd-math --no-default-features --all-targets -- -D warnings`, and `git diff
  --check` passed. The public signature scan shows `spmv`, `spmv_parallel`,
  and `try_spmv` accept `Array1`; only private `try_spmv_dvector` retains
  nalgebra `DVector`.
- **Residual risk**: The public `LinearOperator`/solver/preconditioner trait
  family still exposes nalgebra `DVector`, and sparse storage still aliases
  `nalgebra_sparse::CsrMatrix`.

---

## Sprint 2026-07-04: cfd-math SparseMatrixExt Leto Vector Boundary

- **Resolved sparse extension vector API**:
  `crates/cfd-math/src/sparse/operations.rs` no longer exposes nalgebra
  `DVector` in `SparseMatrixExt::diagonal`, `set_diagonal`, `scale_rows`, or
  `scale_columns`. Those methods now use `leto::Array1`.
- **Resolved Jacobi construction dependency**:
  `crates/cfd-math/src/linear_solver/preconditioners/basic.rs` consumes the
  provider-backed Leto diagonal through `Storage::as_slice` instead of
  iterating a nalgebra `DVector`.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest, all-target clippy, static targeted provider-residue scan, and diff
  hygiene. In `D:/atlas/repos/CFDrs`, `rustup run nightly cargo fmt -p
  cfd-math --check`, `cargo check -p cfd-math --no-default-features --lib`,
  `cargo nextest run -p cfd-math --no-default-features sparse basic
  --status-level fail` (21/21), `cargo clippy -p cfd-math
  --no-default-features --all-targets -- -D warnings`, and `git diff
  --check` passed. A targeted scan found no `DVector` signatures for
  `SparseMatrixExt` diagonal/scaling methods.
- **Residual risk**: This closes only the sparse extension vector API. The
  public linear-solver/preconditioner traits still carry nalgebra `DVector`;
  sparse storage still aliases `nalgebra_sparse::CsrMatrix`.

---

## Sprint 2026-07-04: cfd-math Multigrid Smoother/Cycle Leto Vectors

- **Resolved smoother vector boundary**:
  `crates/cfd-math/src/linear_solver/preconditioners/multigrid/smoothers.rs`
  no longer exposes nalgebra `DVector` or nalgebra scalar traits in smoother
  contracts. Jacobi, Gauss-Seidel, symmetric Gauss-Seidel, SSOR, and
  Chebyshev smoothers now operate on `leto::Array1` through
  `MultigridVector`.
- **Resolved cycle vector boundary**:
  `crates/cfd-math/src/linear_solver/preconditioners/multigrid/cycles.rs`
  now accepts and returns `MultigridVector<f64>` for V-cycle, W-cycle, and
  F-cycle correction paths, including coarsest-level solves through the
  Leto dense bridge.
- **Resolved local SpMV gap**: `crates/cfd-math/src/sparse/operations.rs`
  now exposes `spmv_array`, which delegates the remaining
  `nalgebra_sparse::CsrMatrix` storage boundary to `leto_ops::spmv_into`
  while accepting and returning Leto vectors.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest, all-target clippy, static provider-residue scan, and diff hygiene.
  In `D:/atlas/repos/CFDrs`, `rustup run nightly cargo fmt -p cfd-math
  --check`, `cargo check -p cfd-math --no-default-features --lib`, `cargo
  nextest run -p cfd-math --no-default-features multigrid::cycles smoothers
  --status-level fail` (10/10), `cargo clippy -p cfd-math
  --no-default-features --all-targets -- -D warnings`, and `git diff
  --check` passed. The provider-residue scan over `smoothers.rs` and
  `cycles.rs` found no `DMatrix`, `DVector`, `nalgebra::`, `num_traits`,
  `ndarray`, `rayon`, `tokio`, `rustfft`, or `wgpu` hits.
- **Residual risk**: This closes the multigrid smoother/cycle internal vector
  path only. The public sparse and linear-solver traits still expose
  `nalgebra_sparse::CsrMatrix`, nalgebra `DVector`, and nalgebra
  `RealField`; AMG still bridges at `Preconditioner::apply_to` until that
  boundary migrates to Leto/Eunomia.

---

## Sprint 2026-07-04: cfd-math GMG Leto/Eunomia Migration

- **Resolved GMG dense storage**:
  `crates/cfd-math/src/linear_solver/preconditioners/multigrid/gmg/{mod.rs,transfer.rs}`
  no longer import nalgebra or expose `DMatrix`/`DVector` in GMG hierarchy,
  transfer, FAS, or test signatures.
- **Resolved scalar conversion path**: GMG no longer uses `num_traits`
  `FromPrimitive`/`Zero` conversion paths; grid constants and relaxation
  factors now route through Eunomia `FloatElement`/`NumericElement`.
- **Resolved local dense operations**: GMG owns explicit Leto-array helpers
  for matrix-vector product, residual, vector add/subtract, vector update, and
  L2 norm instead of relying on nalgebra operator overloads.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest, all-target clippy, static provider-residue scan, and diff hygiene.
  In `D:/atlas/repos/CFDrs`, `rustup run nightly cargo fmt -p cfd-math
  --check`, `cargo check -p cfd-math --no-default-features --lib`, `cargo
  nextest run -p cfd-math --no-default-features gmg --status-level fail`
  (5/5), `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings`, and `git diff --check` passed.
- **Residual risk**: This closes the GMG cone only. The public
  linear-solver/preconditioner traits still expose nalgebra
  `DVector`/`RealField`; those remain separate provider migration slices.

---

## Sprint 2026-07-04: cfd-math GMRES Leto/Eunomia Workspace

- **Resolved Krylov dense workspace**:
  `crates/cfd-math/src/linear_solver/gmres/{arnoldi.rs,solver.rs}` now store
  the Krylov basis and Hessenberg matrix in `leto::Array2` instead of nalgebra
  `DMatrix`.
- **Resolved Givens dense workspace**:
  `crates/cfd-math/src/linear_solver/gmres/givens.rs` now stores rotation
  coefficients, least-squares RHS, and triangular-solve output in
  `leto::Array1`/`Array2`, with Eunomia `RealField`/`NumericElement` as the
  scalar operation surface.
- **Resolved chain contract**: `LinearSolverChain` now carries the Eunomia
  real-field bound required by the migrated GMRES constructor.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest, all-target clippy, static provider-residue scan, and diff hygiene.
  In `D:/atlas/repos/CFDrs`, `rustup run nightly cargo fmt -p cfd-math
  --check`, `cargo check -p cfd-math --no-default-features --lib`, `cargo
  nextest run -p cfd-math --no-default-features gmres --status-level fail`
  (21/21), `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings`, and `git diff --check` passed.
- **Residual risk**: This closes GMRES internal dense workspace only. The
  public `LinearOperator`, `Preconditioner`, `IterativeLinearSolver`, and
  `LinearSolver` boundary still exposes nalgebra `DVector` and
  `nalgebra::RealField`; replacing that contract is the next larger
  linear-solver API migration.

---

## Sprint 2026-07-04: cfd-math AMG Restriction Leto Arrays

- **Resolved restriction matrix boundary**:
  `crates/cfd-math/src/linear_solver/preconditioners/multigrid/restriction.rs`
  no longer imports nalgebra or exposes `DMatrix`/`DVector` in its transfer
  utility signatures.
- **Resolved Galerkin path**: `restrict_matrix` now delegates dense
  `R * A * P` products to `leto_ops::MatrixProduct`.
- **Resolved weak vector oracle**: Restriction tests now assert concrete
  `P^T v` values instead of finite-only output.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest, all-target clippy, and static provider-residue scan. In
  `D:/atlas/repos/CFDrs`, `rustup run nightly cargo fmt -p cfd-math --check`,
  `cargo check -p cfd-math --no-default-features --lib`, `cargo nextest run
  -p cfd-math --no-default-features restriction --status-level fail` (7/7),
  and `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings` passed.
- **Residual risk**: This closes the standalone restriction dense-transfer
  module only. Sparse/linear-solver public storage contracts still carry
  nalgebra `DVector` or `nalgebra_sparse::CsrMatrix` and remain separate
  migration slices.

---

## Sprint 2026-07-04: cfd-math Linear Solver Leto Dense Bridge

- **Resolved dense bridge**: `linear_solver::dense_bridge` now owns the
  conversion from the current `nalgebra_sparse::CsrMatrix`/nalgebra `DVector`
  boundary into row-major Leto dense storage and calls `leto_ops::solve`.
- **Resolved fallback path**: `DirectSparseSolver` no longer owns a local dense
  fallback implementation.
- **Resolved multigrid path**: Multigrid cycle coarsest small-system solves no
  longer use the local Gaussian-elimination helper; they consume the shared
  Leto dense bridge.
- **Resolved chain contract**: `LinearSolverChain` carries the Leto real-scalar
  bound required by the direct-solver fallback provider path.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest, and all-target clippy. In `D:/atlas/repos/CFDrs`, `rustup run
  nightly cargo fmt -p cfd-math --check`, `cargo check -p cfd-math
  --no-default-features --lib`, `cargo nextest run -p cfd-math
  --no-default-features multigrid::cycles direct_solver --status-level fail`
  (9/9), and `cargo clippy -p cfd-math --no-default-features --all-targets --
  -D warnings` passed.
- **Residual risk**: This is still a bridge slice. The direct solver public API
  remains `nalgebra_sparse::CsrMatrix`/nalgebra `DVector`, and the primary
  sparse direct solve still uses `rsparse`. A later provider slice should
  replace those contracts with `leto_ops::CsrMatrix`, `leto::Array1`, and an
  Atlas-owned sparse direct-solver path.

---

## Sprint 2026-07-04: cfd-math Sparse Builder Leto Construction Bridge

- **Resolved construction path**: `SparseMatrixBuilder` no longer constructs
  CSR storage by calling `nalgebra_sparse::CsrMatrix::try_from_csr_data`
  directly. It now builds `leto_ops::CsrMatrix`, uses the provider's CSR
  structural validation, and crosses the legacy storage boundary through the
  centralized `sparse::bridge` module.
- **Resolved bridge duplication**: Leto/nalgebra CSR conversion and DVector to
  Leto view construction now live in `crates/cfd-math/src/sparse/bridge.rs`
  instead of being local to sparse operations.
- **Resolved assembly edge**: `ParallelAssembly::block_diagonal` also routes
  its empty matrix case through Leto CSR construction.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest, and all-target clippy. In `D:/atlas/repos/CFDrs`, `rustup run
  nightly cargo check -p cfd-math --no-default-features --lib`, `cargo nextest
  run -p cfd-math --no-default-features sparse --status-level fail` (18/18),
  and `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings` passed.
- **Residual risk**: This still preserves the legacy public sparse type:
  `cfd-math::sparse::SparseMatrix` aliases `nalgebra_sparse::CsrMatrix`, and
  solver APIs still use nalgebra `DVector`. Replacing those public contracts is
  the next larger linear-solver storage migration.

---

## Sprint 2026-07-04: cfd-math Sparse Extension Leto Provider Consumption

- **Resolved consumer path**: `SparseMatrixExt` no longer implements diagonal
  extraction, scalar/value scaling, row scaling, column scaling, Frobenius
  norm, diagonal dominance, or the condition-estimate heuristic with
  CFDrs-local CSR traversal loops. Each operation bridges the current
  `nalgebra_sparse::CsrMatrix` storage boundary into `leto_ops::CsrMatrix` and
  delegates to the provider.
- **Resolved preconditioner contract**: `JacobiPreconditioner::new` now carries
  the Leto scalar bound required by provider-backed diagonal extraction.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest, and all-target clippy. In `D:/atlas/repos/CFDrs`, `rustup run
  nightly cargo fmt -p cfd-math --check`, `cargo check -p cfd-math
  --no-default-features --lib`, `cargo nextest run -p cfd-math
  --no-default-features sparse --status-level fail` (18/18), and `cargo
  clippy -p cfd-math --no-default-features --all-targets -- -D warnings`
  passed.
- **Residual risk**: This remains a bridge slice. `crates/cfd-math/src/sparse`
  and the linear-solver/preconditioner APIs still expose
  `nalgebra_sparse::CsrMatrix` and nalgebra `DVector`; the next larger slice
  should replace those public storage contracts with `leto_ops::CsrMatrix` and
  `leto::Array1`.

---

## Sprint 2026-07-04: cfd-math AMG Leto SpGEMM Consumption

- **Resolved consumer path**: `cfd-math::sparse::try_sparse_sparse_mul` now
  delegates CSR×CSR products to `leto_ops::spgemm` while the surrounding
  linear-solver storage boundary remains `nalgebra_sparse::CsrMatrix`.
- **Resolved AMG call sites**: AMG hierarchy recompute, cached-hierarchy
  setup, and standard hierarchy setup now use the fallible Leto-backed product
  path for `R * A * P` Galerkin operators.
- **Resolved restriction construction**: `cfd-math::sparse::try_sparse_transpose`
  now delegates CSR transpose to `leto_ops::CsrMatrix::transpose`, and AMG
  setup uses it for `R = P^T` instead of `nalgebra_sparse::transpose_as_csc`.
- **Resolved SpMV computation**: `cfd-math::sparse::try_spmv` now delegates
  matrix-vector products to `leto_ops::spmv_into`; the legacy `spmv` and
  `spmv_parallel` entry points no longer contain CFDrs-owned CSR traversal
  loops.
- **Evidence tier**: compile-time integration, focused empirical nextest, and
  all-target clippy. In `D:/atlas/repos/CFDrs`, `rustup run nightly cargo fmt
  -p cfd-math --check`, `cargo check -p cfd-math --no-default-features
  --lib`, `cargo nextest run -p cfd-math --no-default-features sparse
  --status-level fail` (17/17), `cargo nextest run -p cfd-math
  --no-default-features interpolation --status-level fail` (15/15), `cargo
  nextest run -p cfd-math --no-default-features amg --status-level fail`
  (6/6), and `cargo clippy -p
  cfd-math --no-default-features --all-targets -- -D warnings` passed.
- **Residual risk**: This is a consumer bridge, not full sparse/linear-solver
  migration. `crates/cfd-math/src/sparse` and AMG still expose
  `nalgebra_sparse::CsrMatrix` and nalgebra `DVector`; the next provider slice
  should move those boundaries to `leto_ops::CsrMatrix` and `leto::Array1`.

---

## Sprint 2026-07-04: Leto-ops CSR Product Provider Gap

- **Resolved upstream**: `leto_ops::spgemm` now computes CSR×CSR products with
  sorted output rows and exact-zero cancellation removal. `CsrRow::nnz` now
  exposes row cardinality directly.
- **Consumer driver**: `cfd-math` AMG still needs a Leto-owned replacement for
  its `nalgebra_sparse` Galerkin product path before the sparse/linear-solver
  cone can migrate without a CFDrs-local sparse multiply.
- **Evidence tier**: upstream provider compile/clippy/doc and empirical
  nextest. In `D:/atlas/repos/leto`, `rustup run nightly cargo fmt -p
  leto-ops --check`, `cargo check -p leto-ops`, `cargo nextest run -p
  leto-ops --test ops_tests sparse --status-level fail` (14/14),
  `cargo clippy -p leto-ops --all-targets -- -D warnings`, and `cargo doc -p
  leto-ops --no-deps` passed.
- **Residual risk**: CFDrs sparse and AMG still need the consumer-side
  migration from `nalgebra_sparse::CsrMatrix`/`DVector` to
  `leto_ops::CsrMatrix`/`leto::Array1`.

---

## Sprint 2026-07-04: cfd-math SIMD Leto/Eunomia Providers

- **Resolved**: `crates/cfd-math/src/simd` no longer imports nalgebra scalar
  traits or implements vector helpers for nalgebra `DVector`. `SimdVectorOps`
  now targets Leto `Array1`, and generic SIMD/field/vectorization methods use
  Eunomia scalar traits.
- **Resolved**: The SIMD sparse matvec helper now accepts `leto_ops::CsrMatrix`
  and delegates to `leto_ops::spmv`; the external SIMD integration residual
  test builds its CSR matrix from Leto-ops validated parts instead of
  `nalgebra_sparse`.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and all-target clippy. Passing commands: `rustup run
  nightly cargo fmt -p cfd-math --check`; `cargo check -p cfd-math
  --no-default-features --lib`; `cargo nextest run -p cfd-math
  --no-default-features simd --status-level fail` (26/26 passed, 318
  skipped); and `cargo clippy -p cfd-math --no-default-features
  --all-targets -- -D warnings`. The direct scan over
  `crates/cfd-math/src/simd` and `crates/cfd-math/tests/simd_tests.rs` found
  no `use nalgebra`, `nalgebra::`, `DVector`, `DMatrix`, `num_traits`,
  `num_complex`, `ndarray`, `rayon`, `tokio`, `rustfft`, `FromPrimitive`,
  `ToPrimitive`, `T::zero()`, `T::one()`, or `default_epsilon()` residue.
- **Residual risk**: cfd-math sparse and linear-solver provider replacement
  remain separate migration slices.

---

## Sprint 2026-07-04: cfd-math Nonlinear Solver Leto Vectors

- **Resolved**: `crates/cfd-math/src/nonlinear_solver` no longer exposes or
  computes through nalgebra `DVector`, `DMatrix`, or `RealField`. JFNK
  configuration, state vectors, finite-difference JvP, restarted GMRES work
  vectors, convergence history, and tests now use Leto `Array1`/`Array2` plus
  Eunomia scalar traits.
- **Resolved**: Anderson and JFNK now share one local
  `nonlinear_solver::linalg` helper module for Leto vector construction,
  norms, dot products, scaled updates, and the small dense helpers used by
  Anderson's least-squares path.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and all-target clippy. Passing commands: `rustup run
  nightly cargo fmt -p cfd-math --check`; `cargo check -p cfd-math
  --no-default-features --lib`; `cargo nextest run -p cfd-math
  --no-default-features nonlinear --status-level fail` (9/9 passed, 335
  skipped); and `cargo clippy -p cfd-math --no-default-features
  --all-targets -- -D warnings`. The direct scan over
  `crates/cfd-math/src/nonlinear_solver` found no `use nalgebra`,
  `nalgebra::`, `DVector`, `DMatrix`, `num_traits`, `num_complex`, `ndarray`,
  `rayon`, `tokio`, `rustfft`, `FromPrimitive`, `ToPrimitive`, `T::zero()`,
  `T::one()`, or `default_epsilon()` residue.
- **Residual risk**: cfd-math sparse and linear-solver provider replacement
  remain separate migration slices.

---

## Sprint 2026-07-04: cfd-math DG Leto Dense Arrays

- **Resolved**: `crates/cfd-math/src/high_order/dg` no longer exposes or
  computes through nalgebra `DVector` or `DMatrix`. DG basis, solution,
  numerical flux, operator, limiter, solver, and time-integration APIs now use
  Leto `Array1`/`Array2` dense storage.
- **Resolved**: DG mass-matrix projection, derivative, RHS, and implicit
  Newton correction solves now route through Leto dense solve helpers and
  return typed solver errors instead of silently falling back to interpolation
  or residual updates.
- **Resolved**: DG public examples and the DG-related Criterion benchmarks now
  construct and consume Leto arrays.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and all-target clippy. Passing commands: `rustup run
  nightly cargo fmt -p cfd-math --check`; `cargo check -p cfd-math
  --no-default-features --lib`; `cargo check -p cfd-math
  --no-default-features --bench dg_benchmarks`; `cargo check -p cfd-math
  --no-default-features --bench flux_alloc_bench`; `cargo nextest run -p
  cfd-math --no-default-features dg --status-level fail` (62/62 passed, 282
  skipped); and `cargo clippy -p cfd-math --no-default-features
  --all-targets -- -D warnings`; `cargo test --doc -p cfd-math
  --no-default-features` (3 passed, 3 ignored). The direct scan over
  `crates/cfd-math/src/high_order` plus the DG benches found no `use
  nalgebra`, `nalgebra::`, `DVector`, `DMatrix`, `num_traits`, `num_complex`,
  `ndarray`, `rayon`, `tokio`, `rustfft`, `FromPrimitive`, `ToPrimitive`,
  `T::zero()`, or `T::one()` residue.
- **Residual risk**: This closes the high-order provider cone. cfd-math sparse
  and linear-solver provider replacement remain separate migration slices.

---

## Sprint 2026-07-04: cfd-math Spectral Leto Dense Arrays

- **Resolved**: `crates/cfd-math/src/high_order/spectral` no longer exposes
  or computes through nalgebra `DVector` or `DMatrix`. Spectral elements, 1D
  spectral meshes, differential operators, interpolation, quadrature, filters,
  spectral time integration, and spectral assembly now use Leto
  `Array1`/`Array2` storage at dense vector/matrix boundaries.
- **Resolved**: The spectral module now owns shared Leto helpers for
  shape-checked vector construction, matrix-vector multiplication, dot
  products, and stiffness assembly instead of duplicating dense-array loops
  across element and operator code.
- **Resolved**: `SpectralInterp::l2_projection` now returns
  `Result<Array1<f64>>` and maps Leto linear-solve failures to
  `Error::Solver` instead of silently degrading to interpolation.
- **Resolved**: `GlobalAssembly::add_element_matrix` now accepts Leto local
  matrices and optional Leto local RHS values, while `SparseMatrixCSR::
  to_dense` materializes a Leto `Array2`.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and all-target clippy. Passing commands: `rustup run
  nightly cargo fmt -p cfd-math --check`; `cargo check -p cfd-math
  --no-default-features --lib`; `cargo nextest run -p cfd-math
  --no-default-features spectral --status-level fail` (13/13 passed, 331
  skipped); and `cargo clippy -p cfd-math --no-default-features
  --all-targets -- -D warnings`. The direct scan over
  `crates/cfd-math/src/high_order/spectral` found no `use nalgebra`,
  `nalgebra::`, `DVector`, `DMatrix`, `num_traits`, `num_complex`, `ndarray`,
  `rayon`, `tokio`, `rustfft`, `FromPrimitive`, `ToPrimitive`, `T::zero()`,
  or `T::one()` residue.
- **Residual risk**: This closes the high-order spectral dense-array cone.
  Sparse and linear-solver provider replacement remain separate migration
  slices; SIMD is closed by the later July 4 provider slice.

---

## Sprint 2026-07-04: cfd-math WENO Eunomia Scalar Follow-up

- **Resolved**: `crates/cfd-math/src/high_order/weno/{mod,weno7}.rs` no
  longer imports or binds on nalgebra scalar traits. `WENO5`, `WENO7`, and
  `WenoReconstruction` now use the Eunomia `RealField`/`FloatElement` scalar
  provider surface.
- **Resolved**: The WENO constant-conversion path remains a single
  provider-owned helper over `FloatElement`, and nonlinear-weight squaring uses
  the existing generic `squared` helper instead of requiring a nalgebra `powi`
  dispatch.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and all-target clippy. Passing commands: `rustup run
  nightly cargo fmt -p cfd-math --check`; `cargo check -p cfd-math
  --no-default-features --lib`; `cargo nextest run -p cfd-math
  --no-default-features weno --status-level fail` (6/6 passed, 338 skipped);
  and `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings`. The direct scan over `crates/cfd-math/src/high_order/weno` found
  no `use nalgebra`, `nalgebra::`, `num_traits`, `num_complex`, `ndarray`,
  `rayon`, `tokio`, `rustfft`, `FromPrimitive`, `ToPrimitive`, `T::zero()`, or
  `T::one()` residue.
- **Residual risk**: This closes only the high-order WENO scalar-provider cone.
  Sparse and linear-solver surfaces still carry broader nalgebra/Leto provider
  replacement work; SIMD is closed by the later July 4 provider slice.

---

## Sprint 2026-07-04: cfd-math Time-Stepping State Leto/Eunomia Follow-up

- **Resolved**: `crates/cfd-math/src/time_stepping/{traits,runge_kutta,
  adaptive,rk_chebyshev}.rs` now expose Leto-backed `TimeState<T>` state
  vectors and bind scalar operations on Eunomia `RealField`/`FloatElement`
  instead of nalgebra `DVector`/`RealField` or `num_traits` identities.
- **Resolved**: `crates/cfd-math/src/time_stepping/imex.rs` now uses
  `TimeState<T>` for ODE state, explicit and implicit stage vectors, error
  norms, and `TimeStepper` integration. IMEX jacobians now use the shared
  Leto-backed `TimeMatrix<T>` alias and Newton corrections delegate to
  `leto_ops::solve` instead of the former nalgebra dense LU bridge.
- **Resolved**: `crates/cfd-math/benches/{rk4_bench,imex_bench}.rs` now
  compile against the migrated `TimeState<f64>` API and construct benchmark
  vectors with Leto `Array1`; the IMEX benchmark now constructs its zero
  jacobian as `TimeMatrix<f64>`.
- **Resolved**: `crates/cfd-math/src/time_stepping/exponential.rs` now uses
  Leto-backed `TimeState<T>`/`TimeMatrix<T>` and delegates dense matrix
  exponential evaluation to `leto_ops::matexp`. The local nalgebra
  `DMatrix`/`DVector`/`ComplexField` surface and truncated matrix-exponential
  power series were removed.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and all-target clippy. Passing commands: `rustup run
  nightly cargo fmt -p cfd-math --check`; `cargo check -p cfd-math
  --no-default-features --lib`; `cargo nextest run -p cfd-math
  --no-default-features runge_kutta --status-level fail` (5/5 passed, 328
  skipped); `cargo nextest run -p cfd-math --no-default-features adaptive
  --status-level fail` (3/3 passed, 330 skipped); `cargo nextest run -p
  cfd-math --no-default-features chebyshev --status-level fail` (5/5 passed,
  328 skipped); `cargo nextest run -p cfd-math --no-default-features imex
  --status-level fail` (5/5 passed, 328 skipped); `cargo nextest run -p
   cfd-math --no-default-features exponential --status-level fail` (6/6 passed,
   329 skipped); and `cargo clippy -p cfd-math --no-default-features
  --all-targets -- -D warnings`. The post-IMEX direct scan over
  `crates/cfd-math/src/time_stepping`, `crates/cfd-math/benches/imex_bench.rs`,
  and `crates/cfd-math/benches/rk4_bench.rs` found no nalgebra, `DMatrix`,
  `DVector`, `num_traits`, `num_complex`, or ndarray hits.
- **Residual risk**: This closes the explicit/adaptive/RKC vector-state
  surface, exponential dense matrix-function surface, and IMEX dense-solve
  bridge, not the full `cfd-math` migration. Broader `cfd-math` still has
  nalgebra-owned sparse and linear-solver surfaces outside the
  time-stepping/differentiation/SIMD cones.

---

## Sprint 2026-07-04: cfd-math Differentiation Leto/Eunomia Follow-up

- **Resolved**: `crates/cfd-math/src/differentiation/{finite_difference,
  gradient,tests}.rs` now uses Leto `Array1` for 1D derivative outputs and
  Leto `Vector3` for gradient/divergence/curl vector surfaces instead of
  nalgebra `DVector`/`Vector3`.
- **Resolved**: Differentiation scalar constants and identities route through
  Eunomia `FloatElement`/`NumericElement`, and the type-suffixed
  `first_derivative_simd_f32` method was replaced by the
  `FiniteDifference<f32>::first_derivative_simd` inherent method.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and all-target clippy. Passing commands: `rustup run
  nightly cargo fmt -p cfd-math --check`; `cargo check -p cfd-math
  --no-default-features --lib`; `cargo nextest run -p cfd-math
  --no-default-features differentiation --status-level fail` (12/12 passed,
  323 skipped); and `cargo clippy -p cfd-math --no-default-features
  --all-targets -- -D warnings`. The direct scan over
  `crates/cfd-math/src/differentiation` found no nalgebra, `DVector`, old
  scalar identities, `to_subset`, `num_traits`, `num_complex`, ndarray, or
  removed type-suffixed SIMD name.
- **Residual risk**: This closes the differentiation cone, not the full
  `cfd-math` migration. Broader `cfd-math` still has nalgebra-owned sparse
  and linear-solver surfaces; high-order, nonlinear-solver, and SIMD cones are
  closed by later July 4 slices.

---

## Sprint 2026-07-04: cfd-math Stability and cfd-validation Time Integration Leto/Eunomia Follow-up

- **Resolved**: `crates/cfd-math/src/time_stepping/stability/**` no longer
  exposes nalgebra `DMatrix`/`DVector` or nalgebra scalar traits for explicit
  Runge-Kutta stability analysis. The RK stability-region, absolute-limit, and
  von Neumann RK APIs now accept Leto `Array2`/`Array1` Butcher tableaus and
  bind on Eunomia `RealField`; built-in RK3/RK4 tableaus use Leto constructors.
- **Resolved**: `crates/cfd-validation/src/time_integration/**` no longer owns
  nalgebra vectors or nalgebra scalar traits. Validation Euler/RK2/RK4 state is
  Leto `Array1`, generic identities route through crate-local Eunomia helpers,
  edge-case tests exercise the migrated Leto state API, and the harmonic
  oscillator validator no longer boxes integrator closures behind `dyn Fn`.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and all-target clippy. Passing commands: `rustup run
  nightly cargo fmt -p cfd-math -p cfd-validation --check`; `cargo check -p
  cfd-math --no-default-features --lib`; `cargo check -p cfd-validation
  --no-default-features --lib`; `cargo nextest run -p cfd-math
  --no-default-features stability --status-level fail` (5/5 passed, 328
  skipped); `cargo nextest run -p cfd-validation --no-default-features
  time_integration --status-level fail` (12/12 passed, 419 skipped); `cargo
  clippy -p cfd-math --no-default-features --all-targets -- -D warnings`; and
  `cargo clippy -p cfd-validation --no-default-features --all-targets -- -D
  warnings`. A direct scan over the migrated cone found no `use nalgebra`,
  `nalgebra::`, `DMatrix`, `DVector`, `T::zero()`, `T::one()`, `num_traits`,
  `num_complex`, or `ndarray` residue.
- **Residual risk**: This closes the cfd-math stability RK surface and the
  cfd-validation time-integration cone, not the full CFDrs provider migration.
  Remaining validation nalgebra surfaces are still present in
  `src/adaptive_mesh/**`, `src/literature/**`,
  `src/numerical/{error_metrics,linear_solver,test_cases,validation_result}.rs`,
  `src/analytical_benchmarks.rs`, selected benchmark/benchmarking modules, and
  validation tests such as `tests/analytical_poiseuille_3d.rs`,
  `tests/twelve_steps.rs`, and `tests/complex_boundary_mms_validation.rs`.
  Remaining cfd-math nalgebra surfaces outside this slice include other
  time-stepping modules and linear-solver/sparse storage paths.

---

## Sprint 2026-07-04: cfd-validation Conservation Leto/Eunomia Follow-up

- **Resolved**: `crates/cfd-validation/src/conservation/**` no longer owns
  nalgebra dense matrices/vectors or nalgebra scalar traits. Conservation
  reports, tolerances, history, mass, momentum, energy, angular-momentum,
  vorticity, and geometric-conservation checks now bind on Eunomia `RealField`
  and use Leto `Array1`/`Array2` for flow fields.
- **Resolved**: Generic identities and array indexing in the migrated
  conservation cone now route through crate-local Eunomia scalar helpers and
  Leto multidimensional indexing.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and all-target clippy. Passing commands: `rustup run
  nightly cargo fmt -p cfd-validation --check`; `cargo check -p
  cfd-validation --no-default-features --lib`; `cargo nextest run -p
  cfd-validation --no-default-features conservation --status-level fail`
  (18/18 passed, 413 skipped); and `cargo clippy -p cfd-validation
  --no-default-features --all-targets -- -D warnings`. A direct scan over
  `src/conservation` found no `use nalgebra`, `nalgebra::`, `DMatrix`,
  `DVector`, `T::zero()`, or `T::one()` residue.
- **Residual risk**: This closes the standalone conservation cone, not the
  entire validation crate. Remaining nalgebra surfaces are still present in
  `src/adaptive_mesh/**`, `src/time_integration/**`, `src/literature/**`,
  `src/numerical/{error_metrics,linear_solver,test_cases,validation_result}.rs`,
  `src/analytical_benchmarks.rs`, selected benchmark/benchmarking modules, and
  validation tests such as `tests/analytical_poiseuille_3d.rs`,
  `tests/twelve_steps.rs`, and `tests/complex_boundary_mms_validation.rs`.

---

## Sprint 2026-07-04: cfd-validation Manufactured Eunomia Follow-up

- **Resolved**: `crates/cfd-validation/src/manufactured/**` no longer imports
  nalgebra scalar traits or calls nalgebra `T::zero()`/`T::one()` identities.
  Manufactured-solution traits, MMS implementations, and Richardson MMS
  analysis now bind on Eunomia scalar contracts and use the crate-local
  Eunomia scalar helpers. Navier-Stokes MMS continues to use Leto `Vector2`.
- **Resolved**: The Richardson MMS Poisson bridge keeps the remaining
  cfd-2d solver scalar requirement localized as `cfd_2d::Cfd2dScalar`, rather
  than reintroducing a direct nalgebra bound into the manufactured module.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and all-target clippy. Passing commands: `rustup run
  nightly cargo fmt -p cfd-validation --check`; `cargo check -p
  cfd-validation --no-default-features --lib`; `cargo nextest run -p
  cfd-validation --no-default-features manufactured --status-level fail`
  (50/50 passed, 381 skipped); and `cargo clippy -p cfd-validation
  --no-default-features --all-targets -- -D warnings`. A direct scan over
  `src/manufactured` found no `use nalgebra`, `nalgebra::`, `T::zero()`, or
  `T::one()` residue.
- **Residual risk**: This closes the standalone manufactured-solution cone,
  not the entire validation crate. Remaining nalgebra surfaces are still
  present in `src/adaptive_mesh/**`, `src/time_integration/**`,
  `src/literature/**`,
  `src/numerical/{error_metrics,linear_solver,test_cases,validation_result}.rs`,
  `src/analytical_benchmarks.rs`, selected benchmark/benchmarking modules, and
  validation tests such as `tests/analytical_poiseuille_3d.rs`,
  `tests/twelve_steps.rs`, and `tests/complex_boundary_mms_validation.rs`.

---

## Sprint 2026-07-04: cfd-validation Edge-Case Testing Eunomia Follow-up

- **Resolved**: `crates/cfd-validation/src/edge_case_testing/**` no longer
  imports nalgebra scalar traits. `EdgeCaseTestSuite<T>` and its test-case
  implementation now bind on Eunomia `RealField`.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and all-target clippy. Passing commands: `rustup run
  nightly cargo fmt -p cfd-validation --check`; `cargo check -p
  cfd-validation --no-default-features --lib`; `cargo nextest run -p
  cfd-validation --no-default-features edge_case --status-level fail` (15/15
  passed, 416 skipped); and `cargo clippy -p cfd-validation
  --no-default-features --all-targets -- -D warnings`. A direct scan over
  `src/edge_case_testing`, `src/convergence`, and `src/error_metrics` found no
  `use nalgebra`, `nalgebra::`, `T::zero()`, or `T::one()` residue.
- **Residual risk**: This closes the standalone edge-case testing cone, not the
  entire validation crate. Remaining nalgebra surfaces are still present in
  `src/adaptive_mesh/**`, `src/time_integration/**`,
  `src/literature/**`,
  `src/numerical/{error_metrics,linear_solver,test_cases,validation_result}.rs`,
  `src/analytical_benchmarks.rs`, selected benchmark/benchmarking modules, and
  validation tests such as `tests/analytical_poiseuille_3d.rs`,
  `tests/twelve_steps.rs`, and `tests/complex_boundary_mms_validation.rs`.

---

## Sprint 2026-07-04: cfd-validation Convergence Leto/Eunomia Follow-up

- **Resolved**: `crates/cfd-validation/src/convergence/**` no longer imports
  nalgebra scalar traits or calls nalgebra `T::zero()`/`T::one()` identities.
  Convergence order classification, convergence monitors, grid convergence
  index, Richardson extrapolation, and convergence-study regression helpers now
  bind on Eunomia `RealField`/`FloatElement` and use the crate-local Eunomia
  scalar helpers.
- **Resolved**: The Richardson MMS result holders that embed
  `ConvergenceStudy<T>` now carry the Eunomia field bound required by the
  migrated convergence API; the rest of the manufactured Richardson module
  remains a separate provider-migration cone.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and all-target clippy. Passing commands: `rustup run
  nightly cargo fmt -p cfd-validation --check`; `cargo check -p
  cfd-validation --no-default-features --lib`; `cargo nextest run -p
  cfd-validation --no-default-features convergence --status-level fail`
  (27/27 passed, 404 skipped); and `cargo clippy -p cfd-validation
  --no-default-features --all-targets -- -D warnings`. A direct scan over
  `src/convergence` found no `use nalgebra`, `nalgebra::`, `T::zero()`,
  `T::one()`, or `f64::INFINITY` residue.
- **Residual risk**: This closes the standalone convergence cone, not the
  entire validation crate. Remaining nalgebra surfaces are still present in
  `src/adaptive_mesh/**`, `src/time_integration/**`,
  `src/literature/**`,
  `src/numerical/{error_metrics,linear_solver,test_cases,validation_result}.rs`,
  `src/analytical_benchmarks.rs`, selected benchmark/benchmarking modules, and
  validation tests such as `tests/analytical_poiseuille_3d.rs`,
  `tests/twelve_steps.rs`, and `tests/complex_boundary_mms_validation.rs`.

---

## Sprint 2026-07-04: cfd-validation Error Metrics Leto/Eunomia Follow-up

- **Resolved**: `crates/cfd-validation/src/error_metrics/**` no longer imports
  nalgebra scalar or vector types. Scalar error metrics, normalized metrics,
  statistics, convergence analysis, and vector error helpers now bind on
  Eunomia `RealField`/`FloatElement` and use Leto `Vector3` magnitudes.
- **Resolved**: Generic identities and reductions in the migrated error-metric
  cone now route through the crate-local Eunomia scalar helpers, with explicit
  generic type parameters where inference is ambiguous.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and all-target clippy. Passing commands: `rustup run
  nightly cargo fmt -p cfd-validation --check`; `cargo check -p
  cfd-validation --no-default-features --lib`; `cargo nextest run -p
  cfd-validation --no-default-features error_metrics --status-level fail`
  (21/21 passed, 410 skipped); and `cargo clippy -p cfd-validation
  --no-default-features --all-targets -- -D warnings`. A direct scan over
  `src/error_metrics` found no `use nalgebra`, `nalgebra::`, `T::zero()`, or
  `T::one()` residue.
- **Residual risk**: This closes the standalone error-metric cone, not the
  entire validation crate. Remaining nalgebra surfaces are still present in
  `src/adaptive_mesh/**`, `src/time_integration/**`,
  `src/literature/**`,
  `src/numerical/{error_metrics,linear_solver,test_cases,validation_result}.rs`,
  `src/analytical_benchmarks.rs`, selected benchmark/benchmarking modules, and
  validation tests such as `tests/analytical_poiseuille_3d.rs`,
  `tests/twelve_steps.rs`, and `tests/complex_boundary_mms_validation.rs`.

---

## Sprint 2026-07-04: cfd-validation Analytical Leto/Eunomia Follow-up

- **Resolved**: `crates/cfd-validation/src/analytical/**` and
  `crates/cfd-validation/src/solutions/mod.rs` no longer expose nalgebra
  `Vector3` or `RealField` in their analytical-solution traits and
  implementations. Blasius, Couette, Poiseuille, 2D rheological Poiseuille,
  Stokes, Taylor-Green, Womersley, and analytical utility surfaces now use
  Leto `Vector3` and Eunomia scalar contracts.
- **Resolved**: Generic zero/one identities in the migrated analytical cone
  now route through the crate-local Eunomia scalar helpers instead of nalgebra
  `T::zero()`/`T::one()` methods.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and all-target clippy. Passing commands: `rustup run
  nightly cargo fmt -p cfd-validation --check`; `cargo check -p
  cfd-validation --no-default-features --lib`; `cargo nextest run -p
  cfd-validation --no-default-features analytical --status-level fail` (20/20
  passed, 411 skipped); and `cargo clippy -p cfd-validation
  --no-default-features --all-targets -- -D warnings`. A direct scan over
  `src/analytical` and `src/solutions` found no nalgebra `RealField`/`Vector3`
  imports or old `T::zero()`/`T::one()` residue.
- **Residual risk**: This closes the analytical solution trait cone, not the
  entire validation crate. Remaining nalgebra surfaces are still present in
  adaptive mesh, conservation, time integration, manufactured/literature,
  `analytical_benchmarks.rs`, and selected benchmark/numerical modules.

---

## Sprint 2026-07-04: cfd-validation Geometry Leto/Eunomia Follow-up

- **Resolved**: `crates/cfd-validation/src/geometry/**` no longer exports
  nalgebra points or imports nalgebra `RealField` in the migrated geometry
  cone. The 2D/3D geometry traits and rectangular/circular/annular/
  bifurcation/serpentine/trifurcation/venturi implementations use Leto
  point/vector types and Eunomia scalar helpers.
- **Resolved**: The directly dependent bifurcation, serpentine, trifurcation,
  and venturi benchmark wrappers now bind their geometry fields on
  `eunomia::RealField`, while keeping `ValidationScalar` for solver-backed
  cross-crate calls. Validation tests that pass velocity vectors into migrated
  cfd-core boundary conditions now construct Leto vectors at that API boundary.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and all-target clippy. Passing commands: `rustup run
  nightly cargo fmt -p cfd-validation --check`; `cargo check -p
  cfd-validation --no-default-features --lib`; `cargo clippy -p
  cfd-validation --no-default-features --all-targets -- -D warnings`; and
  `cargo nextest run -p cfd-validation --no-default-features geometry
  --status-level fail` (11/11 passed, 420 skipped). A direct scan over the
  migrated geometry/benchmark cone found no `nalgebra::RealField`, nalgebra
  point re-export, `T::zero()`, or `T::one()` residue.
- **Residual risk**: This closes the validation geometry and immediate
  benchmark boundary, not the entire validation crate. Later sections close
  the error-metric and convergence cones; other validation modules outside
  those cones still contain nalgebra vector/scalar surfaces and need their own
  Leto/Eunomia slices.

---

## Sprint 2026-07-04: cfd-2d/cfd-3d/cfd-validation Scalar Boundary

- **Resolved**: `cfd-2d` now exposes `Cfd2dScalar` as the crate-level scalar
  seam for migrated cfd-core boundary/fluid contracts. NS-FVM construction and
  moving-wall tests use Leto `Vector3` for boundary velocities, and the direct
  nalgebra scalar residue scan over `crates/cfd-2d/src` finds nalgebra
  `RealField` only in `crates/cfd-2d/src/scalar.rs`.
- **Resolved**: `cfd-3d` now exposes `Cfd3dScalar`. FEM problem/stress/
  validation and DES cfd-core fluid/boundary integrations use that seam, while
  cfd-3d branch/cascade/serpentine/trifurcation/venturi boundary-condition
  construction emits Leto `Vector3` values instead of nalgebra matrices.
- **Resolved**: `cfd-validation` now exposes `ValidationScalar` for
  cross-fidelity cfd-1d/2d/3d validation callers. Casson Poiseuille,
  Womersley, and solver-backed 2D/3D benchmarks use the migrated scalar seam.
- **Correctness fix**: LES Smagorinsky GPU updates now refresh Yoshizawa SGS
  kinetic-energy and dissipation diagnostics after GPU viscosity computation;
  before this fix the GPU path could leave diagnostic fields at zero while
  viscosity was nonzero.
- **Evidence tier**: static source audit, compile-time integration, empirical
  nextest, and scoped clippy. Passing commands: `rustup run nightly cargo fmt
  -p cfd-2d --check`; `cargo check -p cfd-2d --no-default-features --features
  gpu --lib`; `cargo check -p cfd-2d --no-default-features --lib`; `cargo
  clippy -p cfd-2d --no-default-features --features gpu --lib -- -D
  warnings`; `cargo nextest run -p cfd-2d --no-default-features --features gpu
  --lib --status-level fail` (518/518 passed, 1 skipped); `cargo check -p
  cfd-3d --no-default-features --lib`; `cargo clippy -p cfd-3d
  --no-default-features --lib -- -D warnings`; `cargo check -p cfd-validation
  --no-default-features --lib`; and `cargo clippy -p cfd-validation
  --no-default-features --lib -- -D warnings`.
- **All-targets follow-up**: The former cfd-2d all-targets blocker is closed.
  `examples/blood_venturi.rs` and `tests/simplec_pimple_validation.rs` now use
  Leto boundary vectors and the `Cfd2dScalar` integration-test seam; the
  example/test lint blockers in TPMS blood, bifurcation schematic, Ghia cavity,
  grid/refinement, energy, streamtube, DES, LBM, Poiseuille, WENO, and related
  solver tests were resolved without relaxing assertions or solver workloads.
  Evidence: `cargo clippy -p cfd-2d --no-default-features --features gpu
  --all-targets -- -D warnings` and `cargo nextest run -p cfd-2d
  --no-default-features --features gpu --status-level fail` (572/572 passed,
  27 skipped).
- **Residual risk**: This closes the cfd-2d crate-level all-targets blocker,
  not the full CFDrs objective. Broader replacement of remaining ndarray/
  nalgebra, raw WGPU, rayon/tokio, rustfft, and local allocation surfaces must
  continue crate by crate with provider gaps pushed into the owning Atlas repos.

---

## Sprint 2026-07-04: cfd-1d Eunomia/Leto Scalar Boundary

- **Resolved**: `cfd-1d` no longer exposes scattered nalgebra `RealField`
  imports as its domain scalar contract. The crate now routes domain, network,
  component, resistance, vascular, solver, transient, and analysis generic
  bounds through `Cfd1dScalar`, which explicitly combines the remaining
  nalgebra linear-system backend requirement with the Eunomia scalar provider
  contract consumed by migrated `cfd-core` APIs.
- **Geometry contract corrected**: `NetworkDomain::contains_1d` now accepts
  Leto `Point1<T>`, matching the migrated `cfd-core::geometry::Domain`
  contract.
- **Evidence tier**: static source audit, compile-time integration, full
  empirical nextest, and scoped clippy. `rustup run nightly cargo fmt -p
  cfd-1d --check`, `cargo check -p cfd-1d --no-default-features --lib`, `cargo
  nextest run -p cfd-1d --no-default-features --status-level fail` (725/725, 3
  skipped), and `cargo clippy -p cfd-1d --no-default-features --lib -- -D
  warnings` passed. A direct scan found nalgebra `RealField` only in
  `crates/cfd-1d/src/scalar.rs`, where it documents the remaining matrix
  backend boundary.
- **Residual risk**: `cargo check -p cfd-2d --no-default-features --features
  gpu --lib` now reaches `cfd-2d` and fails on `cfd-2d`'s own nalgebra
  `RealField` bounds around migrated `cfd-core` boundary/fluid APIs. `cargo
  clippy -p cfd-1d --no-default-features --all-targets -- -D warnings` remains
  blocked by pre-existing lint debt in tests, examples, and cell-separation
  modules outside this scalar-boundary slice.

---

## Sprint 2026-07-04: cfd-core GPU Poisson Hephaestus Kernels

- **Resolved**: `crates/cfd-core/src/compute/gpu/poisson_solver.rs` no longer
  owns raw WGPU compute pipelines, bind-group layouts, parameter buffers,
  staging buffers, manual `map_async` readback, or `futures`/`mpsc` mapping
  channels. Jacobi, red-black, and residual entry points now dispatch through
  Hephaestus `WgslMultiStorageKernel`.
- **Provider boundary tightened**: Poisson field/source/residual storage now
  uses `WgpuDevice`'s `ComputeDevice` upload, allocation, and download
  contracts. The public constructor still accepts the existing WGPU device and
  queue handles, then immediately wraps them in a Hephaestus provider.
- **Shape contract corrected**: The solver now stores `nx`, `ny`, `dx`, and
  `dy` from construction and validates `phi.len() == source.len() == nx * ny`
  before dispatch. The previous implementation inferred a square grid from
  `phi.len()`, ignoring the constructor geometry.
- **Evidence tier**: static source audit, compile-time integration, full
  empirical nextest, and scoped clippy. `rustup run nightly cargo fmt -p
  cfd-core --check`, `cargo check -p cfd-core --features gpu`, `cargo check -p
  cfd-core --no-default-features`, `cargo clippy -p cfd-core --features gpu
  --all-targets -- -D warnings`, and full `cargo nextest run -p cfd-core
  --features gpu --status-level fail` (231/231) passed. A direct scan of
  `poisson_solver.rs` found no `create_buffer_init`, `create_buffer(`,
  `params_buffer`, `ComputePipeline`, `BindGroupLayout`, `map_async`,
  `futures::channel`, `poll(wgpu::PollType`, or `std::sync::mpsc` residue.
- **Residual risk**: `cfd-2d` accelerated Poisson consumer verification now
  reaches `cfd-2d` and remains blocked by `cfd-2d` Eunomia/nalgebra trait-bound
  errors before the consumer reaches the GPU Poisson path. Broader CFDrs GPU
  cleanup still has raw WGPU orchestration in non-Poisson kernels and tests.

---

## Sprint 2026-07-04: cfd-core GPU Hephaestus/Eunomia Provider

- **Resolved**: `crates/cfd-core/src/compute/gpu/**` no longer imports
  nalgebra scalar contracts. GPU buffer, pipeline, and kernel bounds now use
  `eunomia::RealField`.
- **Provider boundary tightened**: `GpuBuffer<T>` now stores a Hephaestus
  `WgpuBuffer<T>` and routes allocation, host upload, device download, and
  overwrite writes through `WgpuDevice`'s `ComputeDevice` implementation. Raw
  `wgpu::Buffer` access remains only through `GpuBuffer::buffer()` for the
  existing CFDrs bind groups.
- **Evidence tier**: static source audit, compile-time integration, full
  empirical nextest, and scoped clippy. `rustup run nightly cargo fmt -p
  cfd-core --check`, `cargo check -p cfd-core --features gpu`, `cargo check -p
  cfd-core --no-default-features`, `cargo clippy -p cfd-core --features gpu
  --all-targets -- -D warnings`, full `cargo nextest run -p cfd-core
  --no-default-features --status-level fail` (201/201), and full `cargo
  nextest run -p cfd-core --features gpu --status-level fail` (231/231)
  passed. A direct scan over `crates/cfd-core/src/compute/{gpu,traits.rs}`
  found no `nalgebra::RealField`, `use nalgebra`, `T::zero()`, `T::one()`, or
  `num_traits` residue.
- **Residual risk**: This closes the selected GPU scalar/buffer provider cone,
  not the final Hephaestus ownership model for all CFD GPU kernels. CFDrs still
  owns raw WGSL kernel strings and bind-group/pipeline assembly; reusable CFD
  kernel/pipeline orchestration should move into Hephaestus with contract tests
  before the downstream raw WGPU assembly can be deleted.

---

## Sprint 2026-07-04: cfd-core Mesh/Staggered Leto-Eunomia Geometry

- **Resolved**: `crates/cfd-core/src/geometry/mesh/**` and
  `crates/cfd-core/src/geometry/staggered.rs` no longer import nalgebra.
  `Mesh<T>` stores `leto::geometry::Point3<T>`, mesh translation/bounds/
  centroid paths use Leto point/vector arithmetic, mesh rotation accepts
  `leto::FixedMatrix<T, 3, 3>`, and mesh/staggered scalar contracts use
  `eunomia::RealField`. Old nalgebra `T::zero()`/`T::one()` identity calls are
  removed from this cone.
- **Provider gap closed**: Leto `FixedMatrix<T, 3, 3>` now multiplies
  `leto::geometry::Vector3<T>` directly, so CFDrs mesh rotation no longer needs
  nalgebra `Matrix3 * Vector3`.
- **Evidence tier**: static source audit, compile-time integration, focused and
  full empirical nextest, and scoped clippy. `rustup run nightly cargo fmt -p
  cfd-core --check`, `cargo check -p cfd-core --no-default-features`, focused
  `cargo nextest run -p cfd-core --no-default-features mesh staggered
  --status-level fail` (12/12), full `cargo nextest run -p cfd-core
  --no-default-features --status-level fail` (201/201), and `cargo clippy -p
  cfd-core --no-default-features --all-targets -- -D warnings` passed. Provider
  `cargo fmt -p leto --check`, `cargo check -p leto`, full `cargo nextest run
  -p leto --status-level fail` (171/171), and `cargo clippy -p leto
  --all-targets -- -D warnings` passed. A direct scan over the mesh/staggered
  cone found no `nalgebra`, `T::zero()`, `T::one()`, `default_epsilon`,
  `z_axis`, or `.scale(` residue.
- **Residual risk**: This is a geometry-storage migration, not the final Gaia
  mesh-topology replacement. CFDrs still owns `Mesh<T>`/`Element` topology; the
  next mesh increment should replace or bridge that topology with Gaia-owned
  `IndexedMesh`/topology primitives, with explicit conversion tests rather than
  a compatibility shim.

---

## Sprint 2026-07-04: cfd-core Geometry/Boundary Eunomia-Leto Contract

- **Resolved**: `crates/cfd-core/src/geometry/shapes.rs` no longer imports
  nalgebra point/vector/scalar contracts. Domain shape primitives now use
  `leto::geometry::{Point1, Point2, Point3, Vector3}` and
  `eunomia::RealField`. The scalar contract was propagated through the directly
  dependent problem/aggregate management, solver factory/trait/config helpers,
  material/fluid-dynamics services, and boundary condition/geometry/applicator
  family. Old nalgebra identities and epsilon helpers in that cone now use
  Eunomia constants.
- **Provider gap closed**: Leto now supplies `Point1<T>` for one-dimensional
  domains, conditional `Eq` derives on fixed geometry values for downstream
  generic enum derives, and serde `std`/`alloc` feature propagation so direct
  Leto checks compile the existing Vec-backed serde surfaces.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and scoped clippy. `rustup run nightly cargo fmt -p
  cfd-core --check`, `cargo check -p cfd-core --no-default-features`, full
  `cargo nextest run -p cfd-core --no-default-features --status-level fail`
  (201/201), and `cargo clippy -p cfd-core --no-default-features
  --all-targets -- -D warnings` passed. Provider `cargo check -p leto`,
  `cargo nextest run -p leto --status-level fail` (170/170), and `cargo
  clippy -p leto --all-targets -- -D warnings` passed. A direct scan over the
  migrated cone found no nalgebra scalar imports, nalgebra geometry imports,
  `T::zero()`, `T::one()`, `default_epsilon`, `.scale(`, or `z_axis` residue.
- **Residual risk**: This closes the selected geometry/boundary/fluid-management
  cone only. Mesh, GPU/compute backends, and other unmigrated CFDrs surfaces
  still need their own provider slices before nalgebra/raw GPU dependencies can
  be removed from the package graph.

---

## Sprint 2026-07-04: cfd-core Compute Execution Eunomia Boundary
- **Resolved**: `crates/cfd-core/src/compute/{traits,cpu,dispatch}.rs` no
  longer import or bind on nalgebra. `ComputeKernel`, `ComputeBuffer`,
  `CpuBuffer`, `CpuAdvectionKernel`, and `ComputeDispatcher` now use
  `eunomia::RealField`; CPU buffer initialization, Dirichlet-zero boundary
  values, and upwind sign comparisons use Eunomia identities instead of
  nalgebra `T::zero()`.
- **Adjacent vector cleanup**: `crates/cfd-core/src/abstractions/problem.rs`
  now stores problem gravity as `leto::geometry::Vector3<T>` and the builder
  test asserts the stored Leto vector value. Direct nalgebra vector ownership is
  removed from this file.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and scoped clippy. `cargo fmt -p cfd-core --check`
  passed. `cargo check -p cfd-core --no-default-features` passed. `cargo
  clippy -p cfd-core --no-default-features --all-targets -- -D warnings`
  passed. Focused `cargo nextest run -p cfd-core --no-default-features
  problem_builder cpu compute dispatcher backend --status-level fail` passed
  21/21 tests. A direct scan over the touched files found no nalgebra vector
  ownership, `DVector`, `num_traits`, `ndarray`, `rayon`, `tokio`, `rustfft`,
  or `T::zero()` residue.
- **Residual risk**: `Problem<T>` still imports `nalgebra::RealField` because
  the upstream `Domain<T>` and `FluidTrait<T>` contracts still own nalgebra
  scalar bounds. The solver trait/convergence family remains a separate
  provider slice after those contracts move to Eunomia.

---

## Sprint 2026-07-04: cfd-core State Leto/Eunomia Storage
- **Resolved**: `crates/cfd-core/src/abstractions/state.rs` no longer imports
  nalgebra or owns `DVector` scalar field state. `FieldData::Scalar` now stores
  `leto::Array1<T>`, vector fields use `leto::geometry::Vector3<T>`, and the
  state scalar vocabulary is `eunomia::RealField`. Reset and initialization use
  Eunomia identities while preserving value-semantic field mutation and reset
  behavior.
- **Provider gap closed**: Leto owned arrays/storage now implement serde with
  deserialization routed through `Array::new`, so decoded layout/storage
  combinations are validated before constructing an array. This preserves the
  serialized-state contract required by `FieldData` without adding a CFDrs
  compatibility wrapper.
- **Evidence tier**: provider source audit, compile-time integration, focused
  empirical nextest, and scoped clippy. `cargo fmt -p leto --check` passed.
  Focused `cargo nextest run -p leto
  owned_array_round_trips_shape_and_values_through_serde --status-level fail`
  passed 1/1 test. `cargo clippy -p leto --all-targets -- -D warnings`
  passed. `cargo fmt -p cfd-core --check` passed. `cargo check -p cfd-core
  --no-default-features` passed. `cargo clippy -p cfd-core
  --no-default-features --all-targets -- -D warnings` passed. Focused `cargo
  nextest run -p cfd-core --no-default-features field_state field_data
  test_core_state --status-level fail` passed 3/3 tests. A direct scan over
  `crates/cfd-core/src/abstractions/state.rs` found no `nalgebra`, `DVector`,
  `num_traits`, `ndarray`, `rayon`, `tokio`, `wgpu`, or `cuda` residue.
- **Residual risk**: This closes state storage only. `abstractions::problem`,
  `compute::{traits,cpu,gpu,solver}`, mesh, fluid, boundary, material, and GPU
  kernels still need separate Atlas slices before cfd-core can remove
  nalgebra/raw GPU providers from its manifest.

---

## Sprint 2026-07-04: cfd-core Time Integration Leto/Eunomia State
- **Resolved**: `crates/cfd-core/src/compute/time/integrators.rs` no longer
  imports nalgebra or owns `DVector` state. All cfd-core time integrators now
  use `leto::Array1<T>` as `TimeIntegrator::State`, with local contiguous
  slice helpers for scaled updates, sums, and convergence norms. `crates/
  cfd-core/src/compute/time/controllers.rs` now uses Eunomia `RealField`,
  `FloatElement`, and `NumericElement` for tolerances, powers, constants, and
  epsilon checks instead of nalgebra/num-traits.
- **Boundary**: This closes the cfd-core `compute::time` provider boundary. It
  does not claim full cfd-core closure: `compute::solver`, `abstractions::
  state/problem`, GPU kernels, mesh, fluid, and boundary modules still own
  nalgebra/raw GPU provider residue and need separate slices.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and scoped clippy. `cargo fmt -p cfd-core --check`
  passed. `cargo check -p cfd-core --no-default-features` passed. `cargo check
  -p cfd-core --no-default-features --tests` passed. `cargo clippy -p cfd-core
  --no-default-features --lib -- -D warnings` passed. `cargo clippy -p
  cfd-core --all-targets -- -D warnings` passed. Focused `cargo nextest run -p
  cfd-core --no-default-features time forward_euler runge_kutta backward_euler
  crank adaptive_controller variable_controller --status-level fail` passed
  17/17 tests. A direct scan over `crates/cfd-core/src/compute/time` found no
  `nalgebra`, `DVector`, `T::one()`, `default_epsilon`, `ndarray`,
  `num_traits`, `rayon`, `tokio`, `wgpu`, or `cuda` residue.
- **Follow-up closed**: The benchmark feature-hygiene issue is fixed.
  `crates/cfd-core/Cargo.toml` marks `turbulence_benchmark` with
  `required-features = ["gpu"]`, and `cargo clippy -p cfd-core
  --no-default-features --all-targets -- -D warnings` now passes.

---

## Sprint 2026-07-04: cfd-core SolverConfig Eunomia Scalar Boundary
- **Resolved**: `crates/cfd-core/src/compute/solver/config.rs` now defines
  `SolverConfig` and nested convergence, numerical, linear, network, and
  builder config types on `eunomia::RealField`. The default relaxation value
  uses Eunomia's numeric constant, and `crates/cfd-3d/src/spectral/solver.rs`
  no longer imports `NalgebraRealField`. The required bound was propagated
  through affected cfd-2d FDM/SIMPLE/pressure-velocity/vorticity config
  holders, cfd-3d FEM config/solver holders, and cfd-validation MMS Richardson
  validation.
- **Boundary**: This is a configuration-scalar slice, not a complete
  `compute::solver` trait migration. `SolverConfiguration` remains dual-bound
  so existing solver users compile, and `cfd_core::compute::solver::{traits,
  convergence,direct,iterative,monitor}` still import `nalgebra::RealField`
  because those APIs depend on `abstractions::Problem<T>`, which remains the
  next upstream scalar-contract migration boundary.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and scoped clippy. `cargo fmt -p cfd-core -p cfd-2d -p
  cfd-3d -p cfd-validation --check` passed. `cargo check -p cfd-core
  --no-default-features`, `cargo check -p cfd-2d --no-default-features --lib`,
  `cargo check -p cfd-3d --no-default-features --lib`, `cargo check -p cfd-3d
  --no-default-features --test property_tests --test domain_solver_validation
  --test poisson_validation`, and `cargo check -p cfd-validation
  --no-default-features --lib` passed. Scoped clippy passed for cfd-core,
  cfd-2d, cfd-3d, and cfd-validation. Focused `cargo nextest run -p cfd-core
  --no-default-features solver_config --status-level fail` passed 2/2 tests.
  Focused `cargo nextest run -p cfd-2d --no-default-features fdm
  pressure_velocity simple vorticity --status-level fail` passed 41/41 tests.
  Focused `cargo nextest run -p cfd-3d --no-default-features spectral poisson
  fem --status-level fail` passed 89/89 tests. A direct scan over
  `spectral::solver` found no `NalgebraRealField`, `nalgebra::RealField`,
  `DMatrix`, `DVector`, `ndarray`, `num_traits`, `rustfft`, `rayon`, `tokio`,
  `wgpu`, or `cuda` residue.
- **Residual risk**: The next solver scalar-provider increment is the
  `abstractions::Problem<T>` contract and the solver trait family that depends
  on it. Until that is migrated, cfd-core cannot remove nalgebra from the
  solver trait cone.

---

## Sprint 2026-07-04: cfd-3d Spectral Poisson Leto/Leto-ops Provider Slice
- **Resolved**: `crates/cfd-3d/src/spectral/poisson.rs` no longer imports or
  uses nalgebra `DMatrix`/`DVector` for the global Chebyshev Poisson solve.
  `PoissonSolver` now builds the Laplacian as `leto::Array2`, assembles
  Kronecker products through `leto_ops::MatrixProduct::kron`, applies boundary
  rows in Leto storage, solves through `leto_ops::MatrixSolve`, and returns
  `leto::Array1`. `PoissonProblem` in `spectral::solver` now carries the
  Leto-ops `RealScalar` bound required by the Poisson solve and no longer
  requires `nalgebra::RealField`.
- **Boundary**: This closes the direct Poisson matrix/vector provider boundary.
  `crates/cfd-3d/src/spectral/solver.rs` still imports `NalgebraRealField`
  because `cfd_core::compute::solver::SolverConfig` remains bound on
  `nalgebra::RealField`; that upstream solver-config scalar contract is a
  separate cfd-core migration slice.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and scoped clippy. `cargo fmt -p cfd-3d --check` passed.
  `cargo check -p cfd-3d --no-default-features --lib` passed. `cargo check -p
  cfd-3d --no-default-features --test robustness_tests --test
  poisson_validation` passed. `cargo clippy -p cfd-3d --no-default-features
  --lib --test robustness_tests --test poisson_validation -- -D warnings`
  passed. Focused `cargo nextest run -p cfd-3d --no-default-features poisson
  --status-level fail` passed 12/12 tests. A direct provider scan over
  `spectral::poisson` and `spectral::solver` shows no `DMatrix`, `DVector`,
  `ndarray`, `num_traits`, `rustfft`, `rayon`, `tokio`, `wgpu`, or `cuda`
  residue; the only remaining match is the `NalgebraRealField` alias in
  `spectral::solver` for `SolverConfig`.
- **Residual risk**: Full spectral closure still requires migrating
  `cfd_core::compute::solver::SolverConfig` off nalgebra scalar bounds and
  continuing the Fourier/Apollo and execution-provider cleanup outside this
  Poisson solve boundary.

---

## Sprint 2026-07-04: cfd-3d Spectral Chebyshev Leto/Eunomia Provider Slice
- **Resolved**: `crates/cfd-3d/src/spectral/chebyshev.rs` no longer imports or
  stores nalgebra `DMatrix`/`DVector` or nalgebra `RealField`. The Chebyshev
  basis stores differentiation operators as `leto::Array2`, applies first and
  second derivatives over `leto::Array1`, validates mismatched vector lengths
  through `cfd_core::Error::DimensionMismatch`, and uses `eunomia::RealField`
  for scalar vocabulary. `crates/cfd-3d/src/spectral/basis.rs` now uses
  Eunomia `RealField`. Chebyshev tests in `crates/cfd-3d/src/spectral/
  chebyshev_tests.rs`, `crates/cfd-3d/src/lib.rs`, and
  `crates/cfd-3d/tests/robustness_tests.rs` no longer construct nalgebra
  `DVector`s.
- **Boundary**: This closes the Chebyshev representation/API surface, not the
  whole spectral module. The follow-on Poisson section now closes the former
  nalgebra LU/Kronecker global solve boundary. `spectral::solver` carries an
  explicit `NalgebraRealField` alias only because upstream
  `cfd-core::SolverConfig` still requires nalgebra.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and scoped clippy. `cargo fmt -p cfd-3d --check` passed.
  `cargo check -p cfd-3d --no-default-features --lib` passed. `cargo check -p
  cfd-3d --no-default-features --test robustness_tests --test
  poisson_validation` passed. `cargo clippy -p cfd-3d --no-default-features
  --lib --test robustness_tests --test poisson_validation -- -D warnings`
  passed. Focused `cargo nextest run -p cfd-3d --no-default-features
  chebyshev --status-level fail` passed 20/20 tests. Focused `cargo nextest
  run -p cfd-3d --no-default-features poisson --status-level fail` passed
  12/12 tests. A spectral provider scan shows remaining nalgebra hits confined
  to `spectral::poisson` plus the explicit Poisson-boundary alias in
  `spectral::solver`.
- **Residual risk**: Full spectral provider closure still requires migrating
  `cfd-core::SolverConfig` off nalgebra scalar bounds and continuing the
  Fourier/Apollo and execution-provider cleanup outside the Chebyshev/Poisson
  storage boundary.

---

## Sprint 2026-07-04: cfd-3d Level-Set Atlas Provider Closure
- **Resolved**: `crates/cfd-3d/src/level_set/{solver.rs,advection.rs,weno.rs}`
  no longer imports nalgebra or the broad crate-level scalar seam. The
  level-set solver/advection path now uses `leto::geometry::Vector3` for
  velocity storage and a local `LevelSetScalar` trait backed by
  `eunomia::RealField`/`FloatElement`. The dead crate-level `CfdScalar` trait
  and its `nalgebra::RealField` import were removed from
  `crates/cfd-3d/src/scalar.rs`. `crates/cfd-3d/tests/level_set_tests.rs` and
  the level-set sections of `crates/cfd-3d/tests/robustness_tests.rs` now pass
  Leto vectors into `LevelSetSolver`.
- **Boundary**: This closes direct forbidden-provider ownership in the
  level-set source cone and dedicated level-set tests. It does not claim full
  `cfd-3d` provider closure: FEM, spectral Chebyshev/Poisson, IBM, turbulence,
  validation, and non-level-set robustness-test surfaces still own
  nalgebra/ndarray/raw GPU/execution-provider residue.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and scoped clippy. `cargo fmt -p cfd-3d --check` passed.
  `cargo check -p cfd-3d --no-default-features --lib` passed. `cargo check -p
  cfd-3d --no-default-features --test level_set_tests` passed. `cargo check -p
  cfd-3d --no-default-features --test robustness_tests` passed. `cargo clippy
  -p cfd-3d --no-default-features --lib --test level_set_tests --test
  robustness_tests -- -D warnings` passed. Focused `cargo nextest run -p
  cfd-3d --no-default-features level_set --status-level fail` passed 13/13
  tests. A forbidden-provider scan over `crates/cfd-3d/src/level_set` and
  `crates/cfd-3d/tests/level_set_tests.rs` found no `nalgebra`, `DMatrix`,
  `DVector`, `CfdScalar`, `crate::scalar`, `num_traits`, `num-traits`,
  `ndarray`, `rayon`, `tokio`, `rustfft`, `wgpu`, or `cuda` residue.
- **Residual risk**: Full cfd-3d provider closure remains incomplete outside
  level-set and VOF. Next high-value slices are FEM matrix/vector storage to
  Leto, spectral Chebyshev/Poisson storage/FFT provider cleanup through
  Leto/Apollo, raw GPU replacement through Hephaestus, and execution-provider
  cleanup through Moirai.

---

## Sprint 2026-07-04: cfd-3d VOF Atlas Provider Closure
- **Resolved**: `crates/cfd-3d/src/vof/cavitation_solver.rs` no longer exposes
  nalgebra `DMatrix` for cavitation pressure, density, volume-fraction,
  inception, damage, bubble-radius, nuclei, or sonoluminescence fields. These
  dense fields now use `CavitationField = leto::Array2<f64>` with explicit
  row-major Leto offsets and explicit mapping to the VOF flat alpha storage.
  `crates/cfd-3d/tests/cavitation_solver_validation.rs` now constructs and
  indexes cavitation dense fields through Leto `Array2`.
  `crates/cfd-3d/src/vof/scalar.rs` now binds `VofScalar` on
  `eunomia::RealField` instead of
  `nalgebra::RealField`.
- **Boundary**: This closes direct forbidden-provider ownership in the VOF
  source/test cone. It does not remove legacy provider dependencies from
  `cfd-3d`: FEM, spectral Chebyshev/Poisson, level-set, IBM, turbulence, and
  validation surfaces still own nalgebra/ndarray/raw GPU/execution-provider
  residue and need separate crate-level Leto/Gaia/Eunomia/Hephaestus/Moirai
  slices before the manifest can drop those dependencies.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and scoped clippy. `cargo fmt -p cfd-3d --check` passed.
  `cargo check -p cfd-3d --no-default-features --lib` passed. `cargo check -p
  cfd-3d --no-default-features --test cavitation_solver_validation` passed.
  `cargo clippy -p cfd-3d --no-default-features --lib --test
  cavitation_solver_validation -- -D warnings` passed. Focused `cargo nextest
  run -p cfd-3d --no-default-features cavitation --status-level fail` passed
  23/23 tests. Focused `cargo nextest run -p cfd-3d --no-default-features vof
  --status-level fail` passed 42/42 tests. A direct forbidden-provider scan
  over `crates/cfd-3d/src/vof`,
  `crates/cfd-3d/tests/cavitation_solver_validation.rs`, and
  `crates/cfd-3d/tests/vof_tests.rs`
  found no `nalgebra`, `DMatrix`, `DVector`, `num_traits`, `num-traits`,
  `ndarray`, `rayon`, `tokio`, `rustfft`, `wgpu`, or `cuda` residue.
- **Residual risk**: Full cfd-3d provider closure remains incomplete outside
  VOF. Next high-value crate-level slices are FEM matrix/vector storage to Leto,
  spectral Chebyshev/Poisson storage/FFT provider cleanup through Leto/Apollo,
  raw GPU replacement through Hephaestus, and execution-provider cleanup through
  Moirai.

---

## Sprint 2026-07-04: cfd-3d VOF Cavitation Leto Velocity Provider
- **Resolved**: `crates/cfd-3d/src/vof/{cavitation_solver.rs,bubble_dynamics.rs}` and `crates/cfd-3d/tests/cavitation_solver_validation.rs` no longer import or use nalgebra `Vector3`. Cavitation velocity fields and bubble dynamics update calls now use `leto::geometry::Vector3`, and `CavitationVofSolver::step` copies that Leto velocity slice directly into the Leto-owned `VofSolver` velocity buffer.
- **Boundary**: This closes the remaining direct nalgebra `Vector3` ownership inside the VOF/cavitation module family. The later VOF Atlas provider closure replaces the cavitation dense-field residue from this item. Other cfd-3d FEM/spectral/level-set/test surfaces still own nalgebra vectors, matrices, points, and scalar bounds.
- **Evidence tier**: static source audit, compile-time integration, focused empirical nextest, and scoped clippy. `cargo fmt -p cfd-3d --check` passed. `cargo check -p cfd-3d --no-default-features --lib` passed. `cargo check -p cfd-3d --no-default-features --test cavitation_solver_validation` passed. `cargo clippy -p cfd-3d --no-default-features --lib --test cavitation_solver_validation -- -D warnings` passed. Focused `cargo nextest run -p cfd-3d --no-default-features cavitation --status-level fail` passed 23/23 tests. A direct scan over the cavitation VOF cone found no `nalgebra::Vector3`, `use nalgebra::Vector3`, or `use nalgebra::{..., Vector3}` residue.
- **Residual risk**: Full cfd-3d nalgebra removal still requires continuing FEM, spectral Chebyshev/Poisson, level-set, and validation provider migrations.

---

## Sprint 2026-07-04: cfd-3d VOF Leto Vector Provider
- **Resolved**: `crates/cfd-3d/src/vof/{solver.rs,reconstruction.rs,initialization.rs,plic_geometry.rs,advection.rs}` and `crates/cfd-3d/tests/vof_tests.rs` no longer import or use nalgebra `Vector3`. VOF velocity fields, interface normals, initialization geometry, PLIC plane normals, and VOF tests now use `leto::geometry::Vector3`. `crates/cfd-3d/tests/robustness_tests.rs` passes Leto vectors at the VOF API call sites.
- **Boundary**: This closes direct nalgebra `Vector3` ownership for the non-cavitation VOF module family. The later VOF Atlas provider closure replaces the cavitation dense-field residue from this item. Other cfd-3d FEM/spectral/level-set/test surfaces still own nalgebra matrices, points, vectors, and scalar bounds.
- **Evidence tier**: static source audit, compile-time integration, focused empirical nextest, and scoped clippy. `cargo fmt -p cfd-3d --check` passed. `cargo check -p cfd-3d --no-default-features --test vof_tests` passed. `cargo check -p cfd-3d --no-default-features --test robustness_tests` passed. `cargo clippy -p cfd-3d --no-default-features --lib --test vof_tests -- -D warnings` passed. Focused `cargo nextest run -p cfd-3d --no-default-features vof --status-level fail` passed 42/42 tests. A direct scan over the migrated VOF files found no `nalgebra::Vector3`, `use nalgebra::Vector3`, `use nalgebra::{..., Vector3}`, stale `dot(&...)`, or stale vector `into_iter()` residue.
- **Residual risk**: Full cfd-3d nalgebra removal still requires FEM, spectral Chebyshev/Poisson, level-set, and test provider migrations outside VOF.

---

## Sprint 2026-07-04: cfd-validation Direct Vector2 Provider Closure
- **Resolved**: `crates/cfd-validation/src/manufactured/navier_stokes.rs`, `crates/cfd-validation/src/conservation/{momentum.rs,angular_momentum.rs,mod.rs}`, `crates/cfd-validation/src/benchmarks/vorticity_stream.rs`, and MMS/conservation validation tests no longer import or use nalgebra `Vector2`. These surfaces now use `leto::geometry::Vector2`; component access uses Leto `[0]`/`[1]` indexing.
- **Boundary**: This closes direct nalgebra `Vector2` ownership in cfd-validation source and tests. It does not remove nalgebra from cfd-validation because dense `DMatrix`, 3D `Vector3`, and geometry point/vector surfaces remain separate Leto/Gaia provider migrations.
- **Evidence tier**: static source audit, compile-time integration, focused empirical nextest, and scoped clippy. `cargo fmt -p cfd-validation --check` passed. `cargo check -p cfd-validation --no-default-features --tests` passed. `cargo clippy -p cfd-validation --no-default-features --lib -- -D warnings` passed. Focused `cargo nextest run -p cfd-validation --no-default-features manufactured mms conservation taylor momentum angular --status-level fail` passed 104/104 tests. A direct scan over `crates/cfd-validation/{src,tests}` found no `nalgebra::Vector2`, `use nalgebra::Vector2`, or `use nalgebra::{..., Vector2}` residue.
- **Residual risk**: Full `cargo clippy -p cfd-validation --no-default-features --lib --tests -- -D warnings` remains blocked by pre-existing unrelated lints in cross-fidelity, schematics, and other validation tests outside this provider slice. Remaining provider work should move cfd-validation dense matrices to Leto storage and geometry types to Gaia/Leto before nalgebra can be removed from its manifest.

---

## Sprint 2026-07-04: cfd-2d Direct Vector2 Provider Closure
- **Resolved**: `crates/cfd-2d/src/physics/vorticity_stream.rs`, `crates/cfd-2d/src/piso_algorithm/corrector.rs`, `crates/cfd-2d/src/solvers/lbm/solver.rs`, `crates/cfd-2d/src/physics/immersed_boundary.rs`, `crates/cfd-2d/examples/blood_venturi.rs`, and `crates/cfd-2d/benches/solver_benchmarks.rs` no longer import or use nalgebra `Vector2`. These surfaces now use `leto::geometry::Vector2`; component access uses Leto `[0]`/`[1]` indexing. `crates/cfd-validation/src/benchmarks/vorticity_stream.rs` was updated as the downstream public-API consumer of `VorticityStreamSolver::velocity_field`.
- **Boundary**: This closes direct nalgebra `Vector2` ownership for cfd-2d source/tests/examples/benches. It does not remove nalgebra from cfd-2d because scalar `RealField`, boundary `Vector3`, dense `DVector`/`DMatrix`, and nalgebra-sparse matrix surfaces remain separate provider migrations. `physics::immersed_boundary` still uses nalgebra `DMatrix` for force/velocity matrices pending a Leto dense-storage migration.
- **Evidence tier**: static source audit, compile-time integration, focused empirical nextest, and scoped clippy. `cargo fmt -p cfd-2d -p cfd-validation --check` passed. `cargo check -p cfd-2d --no-default-features --examples --benches` passed. `cargo check -p cfd-validation --no-default-features` passed. `cargo clippy -p cfd-2d --no-default-features --example blood_venturi --bench solver_benchmarks -- -D warnings` passed. `cargo clippy -p cfd-validation --no-default-features --lib -- -D warnings` passed. Focused `cargo nextest run -p cfd-2d --no-default-features vorticity corrector lbm immersed --status-level fail` passed 44/44 tests. A direct scan over cfd-2d source/tests/examples/benches plus the touched cfd-validation benchmark found no `nalgebra::Vector2`, `use nalgebra::Vector2`, or `use nalgebra::{..., Vector2}` residue.
- **Residual risk**: Full `cargo clippy -p cfd-2d --no-default-features --all-targets -- -D warnings` remains blocked by pre-existing unrelated lints in examples/tests/modules outside this provider slice. The provider migration still needs cfd-math Leto linear-solver storage, cfd-core boundary vector migration, and Eunomia scalar contract replacement before cfd-2d can remove direct nalgebra/nalgebra-sparse manifest ownership.

---

## Sprint 2026-07-04: cfd-2d Pressure-Velocity/SIMPLEC Leto Vector Workspaces
- **Resolved**: `crates/cfd-2d/src/fields.rs`, `crates/cfd-2d/src/pressure_velocity/{solver.rs,pressure.rs,rhie_chow.rs}`, `crates/cfd-2d/src/simplec_pimple/{algorithms.rs,diagnostics.rs,interpolation.rs,pimple.rs,simplec.rs,solver.rs}`, `crates/cfd-2d/src/physics/momentum/interpolation.rs` reference tests, and `crates/cfd-2d/tests/simplec_pimple_validation.rs` no longer use nalgebra `Vector2`. Velocity workspaces, Rhie-Chow caches, correction buffers, face velocities, and validation setup now use `leto::geometry::Vector2`; component access uses Leto `[0]`/`[1]` indexing.
- **Boundary**: This is a vector-provider migration, not a full nalgebra closure. The same family still uses nalgebra scalar `RealField` where upstream `cfd-core`/`cfd-math` contracts require it, nalgebra `Vector3` for boundary-condition inputs, and nalgebra `DVector` for pressure-correction linear-solver buffers.
- **Evidence tier**: static source audit, compile-time integration, focused empirical nextest, and lib clippy. `cargo fmt -p cfd-2d --check` passed. `cargo check -p cfd-2d --no-default-features` passed. `cargo clippy -p cfd-2d --no-default-features --lib -- -D warnings` passed. Focused `cargo nextest run -p cfd-2d --no-default-features pressure_velocity simplec pimple rhie --status-level fail` passed 29/29 tests. `git diff --check` passed for the touched family with only Git LF-to-CRLF warnings. A direct scan over the touched family found no `nalgebra::Vector2`, `use nalgebra::Vector2`, or `use nalgebra::{RealField, Vector2}` residue.
- **Residual risk**: Full cfd-2d nalgebra removal still requires Eunomia scalar replacement for `RealField`, a cfd-core boundary `Vector3` provider migration, and Leto-backed vector/matrix linear-solver storage replacing `DVector`/sparse nalgebra surfaces.

---

## Sprint 2026-07-04: cfd-2d Ghia Validation Leto Vector Helper
- **Resolved**: `crates/cfd-2d/tests/ghia_cavity_simplec_validation.rs` no
  longer imports `nalgebra::Vector2` for the local pressure-correction helper.
  The divergence-free helper velocity field uses `leto::geometry::Vector2`.
- **Boundary**: This is a test-helper provider cleanup. It does not migrate the
  larger SIMPLEC/PIMPLE validation test or the production pressure-velocity
  vector workspaces, which still require a coordinated Leto vector/storage
  family slice.
- **Evidence tier**: static source audit, compile-time integration, and focused
  empirical nextest. `cargo fmt -p cfd-2d --check` passed. `cargo check -p
  cfd-2d --no-default-features --test ghia_cavity_simplec_validation` passed.
  Focused `cargo nextest run -p cfd-2d --no-default-features
  test_pressure_correction_basic --status-level fail` passed 1/1 test. A direct
  scan over `crates/cfd-2d/tests/ghia_cavity_simplec_validation.rs` found no
  `nalgebra`, `nalgebra-sparse`, `nalgebra_sparse`, `DVector`, or `DMatrix`
  residue.
- **Residual risk**: Remaining cfd-2d test nalgebra ownership is concentrated in
  `crates/cfd-2d/tests/simplec_pimple_validation.rs`; production
  SIMPLEC/PIMPLE and pressure-velocity workspaces still use nalgebra vectors and
  should migrate as one provider-family slice.

---

## Sprint 2026-07-04: cfd-2d Turbulence Validation Eunomia Scalar Contract
- **Resolved**: `crates/cfd-2d/src/physics/turbulence/validation/{mod.rs,rans.rs,les_des.rs,benchmarks.rs}`
  no longer imports or bounds `nalgebra::RealField`. Turbulence validation now
  uses `eunomia::RealField`/`FloatElement` for scalar contracts while retaining
  Leto `Array2`/`Vector2` for validation storage and vectors.
- **Boundary**: This is a turbulence-validation provider seam, not a full
  cfd-2d nalgebra closure. `crates/cfd-2d/Cargo.toml` still directly declares
  `nalgebra` and `nalgebra-sparse` because other cfd-2d modules and tests still
  use nalgebra field, dense vector/matrix, geometric vector, and sparse matrix
  surfaces.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-2d --check` passed.
  `cargo check -p cfd-2d --no-default-features` passed. `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings` passed. Focused `cargo
  nextest run -p cfd-2d --no-default-features turbulence validation
  --status-level fail` passed 207/207 tests. A direct scan over
  `crates/cfd-2d/src/physics/turbulence/validation` found no `nalgebra`,
  `nalgebra-sparse`, `nalgebra_sparse`, `DVector`, or `DMatrix` residue.
- **Residual risk**: Remaining cfd-2d nalgebra/nalgebra-sparse work should move
  module families to Eunomia scalar contracts, Leto CPU arrays/vectors, Gaia mesh
  contracts, and Hephaestus GPU providers before the manifest dependencies can
  be removed.

---

## Sprint 2026-07-04: cfd-2d Direct num-traits Ownership Removal
- **Resolved**: `cfd-2d` no longer owns direct `num-traits` usage in source,
  tests, or its package manifest. The final pass removed direct
  `num_traits::{FromPrimitive,ToPrimitive}` bounds from
  `crates/cfd-2d/src/physics/turbulence/validation/{mod.rs,les_des.rs,benchmarks.rs}`,
  removed the immersed-boundary test's `ToPrimitive` conversion over an
  already-`f64` field, removed concrete `f64` Ghia validation `num_traits`
  imports/bounds, and deleted `num-traits.workspace = true` from
  `crates/cfd-2d/Cargo.toml`.
- **Boundary**: This closes direct `num-traits` ownership for the `cfd-2d`
  crate. It does not claim full provider graph removal because `num-traits`
  still enters the workspace through other crates/transitive dependencies, and
  `cfd-2d` still has direct nalgebra/nalgebra-sparse surfaces that require
  Leto/Eunomia/Gaia provider migration.
- **Evidence tier**: static source/dependency audit, compile-time integration,
  focused empirical nextest, and lib clippy. `cargo fmt -p cfd-2d --check`
  passed. `cargo check -p cfd-2d --no-default-features --tests` passed. `cargo
  clippy -p cfd-2d --no-default-features --lib -- -D warnings` passed. Focused
  `cargo nextest run -p cfd-2d --no-default-features turbulence validation
  immersed ghia --status-level fail` passed 213/213 tests. `cargo tree -p
  cfd-2d --depth 1` shows no direct `num-traits` dependency. `git diff
  --check` passed for touched files and PM artifacts. A direct scan over
  `crates/cfd-2d/{src,tests,Cargo.toml}` found no `num_traits`, `num-traits`,
  `FromPrimitive`, or `ToPrimitive` residue.
- **Residual risk**: Older cfd-2d slice-local audit entries that mention
  remaining direct `num-traits` are superseded by this crate-level closure.
  Remaining cfd-2d provider work is direct nalgebra/nalgebra-sparse removal,
  Leto array/vector migration, and transitive provider cleanup.

---

## Sprint 2026-07-04: cfd-2d Momentum Setup/Boundary Eunomia Scalar Seam
- **Resolved**: `crates/cfd-2d/src/physics/momentum/setup.rs` no longer uses
  direct `num_traits::{FromPrimitive,ToPrimitive}` in the boundary setup apply
  path. `crates/cfd-2d/src/physics/momentum/boundary/{mod.rs,directional.rs}`
  no longer uses direct `num_traits::FromPrimitive`, `T::from_f64(...)`,
  `T::zero()`, `T::one()`, or scalar `.abs()` in the touched boundary assembly
  paths. No-slip/slip/symmetry/periodic/outflow zero/one row entries,
  zero-gradient comparisons, quadratic wall-extrapolation constants, and corner
  consistency absolute-value checks now route through
  `crates/cfd-2d/src/scalar.rs` and `eunomia::FloatElement`.
- **Boundary**: This slice is limited to scalar-provider cleanup in the momentum
  setup and boundary assembly family. It deliberately preserves nalgebra
  `RealField`, `DVector`, and `Vector3` surfaces because pressure/momentum vector
  storage and upstream boundary types remain separate Leto/provider migrations.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-2d --check` passed.
  `cargo check -p cfd-2d --no-default-features` passed. `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings` passed. Focused `cargo
  nextest run -p cfd-2d --no-default-features boundary --status-level fail`
  passed 42/42 tests. `git diff --check` passed for the touched code files. A
  direct-provider scan over the touched files found no direct scalar-provider
  residue.
- **Residual risk**: Direct `num-traits` remains in cfd-2d Ghia validation tests,
  immersed-boundary test code, turbulence validation, and f64-only/test scalar
  surfaces; full cfd-2d direct `num-traits` removal remains a larger crate-level
  Eunomia migration before the manifest dependency can be dropped.

---

## Sprint 2026-07-04: cfd-2d Momentum Interpolation Eunomia Scalar Seam
- **Resolved**: `crates/cfd-2d/src/physics/momentum/interpolation.rs` no longer
  uses direct `num_traits::FromPrimitive`, `T::from_f64(...).unwrap_or_else`,
  `T::zero()`, or `T::one()` in the Rhie-Chow interpolation path. The tiny
  diagonal floor, zero face-array initialization, and harmonic face coefficient
  factor route through `crates/cfd-2d/src/scalar.rs` and
  `eunomia::FloatElement`.
- **Boundary**: This slice is limited to scalar-provider cleanup in momentum
  Rhie-Chow interpolation. It deliberately preserves the test's nalgebra vector
  reference path because `pressure_velocity::RhieChowInterpolation` remains a
  nalgebra-vector boundary and is a separate provider migration.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-2d --check` passed.
  `cargo check -p cfd-2d --no-default-features` passed. `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings` passed. Focused `cargo
  nextest run -p cfd-2d --no-default-features
  coefficient_aware_interpolation_matches_exact_reference --status-level fail`
  passed 1/1 test. `git diff --check` passed for the touched interpolation and
  PM artifact files. A direct-provider scan over the touched file found no direct
  scalar-provider residue.
- **Residual risk**: Direct `num-traits` remains in cfd-2d Ghia validation tests,
  immersed-boundary test code, turbulence validation, and f64-only/test scalar
  surfaces; full cfd-2d direct `num-traits` removal remains a larger crate-level
  Eunomia migration before the manifest dependency can be dropped.

---

## Sprint 2026-07-04: cfd-2d Problem/Streamtube Atlas Provider Seam
- **Resolved**: `crates/cfd-2d/src/problem.rs` now stores incompressible
  problem and solution velocity fields with `leto::geometry::Vector2` and
  routes initial pressure, velocity-magnitude maxima, and pressure maxima
  through `crates/cfd-2d/src/scalar.rs`/Eunomia instead of local nalgebra
  vector storage or direct `T::zero()` folds.
  `crates/cfd-2d/src/physics/streamtube/partitioning.rs` no longer uses direct
  `num_traits::{Float,FromPrimitive}`, `T::from_f64(...).unwrap()`,
  `T::zero()`, `T::one()`, `Float::abs`, `Float::sqrt`, or scalar `.abs()` in
  the touched APIs and tests; constants, absolute values, and square roots now
  route through `eunomia::{FloatElement,NumericElement}` and the crate-local
  scalar adapter.
- **Boundary**: This slice is limited to the problem setup and streamtube
  partitioning scalar/vector-provider seam. `problem.rs` still carries
  `nalgebra::RealField` because `cfd_core::physics::{boundary,fluid}` types
  are still nalgebra-bound upstream; removing that bound requires an upstream
  cfd-core provider migration. It does not remove the cfd-2d manifest's direct
  `num-traits` dependency because direct residues remain outside this slice.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-2d --check` passed.
  `cargo check -p cfd-2d --no-default-features` passed. `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings` passed. Focused `cargo
  nextest run -p cfd-2d --no-default-features problem streamtube separating
  --status-level fail` passed 4/4 tests. `git diff --check` passed for the
  touched problem/streamtube and PM artifact files. A direct-provider scan over
  both touched files found no `num_traits`, `FromPrimitive`, `Float::`,
  `T::from_*`, `.to_f64()`, `T::zero()`, `T::one()`, scalar `.abs()`, or
  local nalgebra `Vector2` residue.
- **Residual risk**: Direct `num-traits` remains in cfd-2d immersed-boundary
  tests, momentum setup/interpolation/boundary, turbulence validation, and
  f64-only/test scalar surfaces; full cfd-2d direct `num-traits` removal
  remains a larger crate-level Eunomia migration before the manifest
  dependency can be dropped. Full `problem.rs` nalgebra-bound removal is blocked
  by upstream `cfd-core` boundary/fluid contracts.

---

## Sprint 2026-07-04: cfd-2d FVM Atlas Provider Seam
- **Resolved**: `crates/cfd-2d/src/solvers/fvm/{config.rs,flux.rs,geometry.rs,solver.rs}`
  no longer use direct `num_traits`, `T::from_*`, `.to_f64()`, `T::zero()`,
  `T::one()`, scalar `.abs()`, or `nalgebra` imports in the touched APIs and
  tests. FVM defaults, face generation constants, residual magnitudes, max
  terms, flux Peclet calculations, nonfinite diffusion validation, and module
  test absolute-value assertions route through `crates/cfd-2d/src/scalar.rs`
  and `eunomia::{FloatElement,NumericElement}`. Face geometry and solver
  velocity vectors use `leto::geometry::Vector2`; power-law and hybrid flux
  calculators retain diffusion as `T` instead of narrowing through `f64`.
- **Boundary**: This slice is limited to the FVM scalar and vector-provider
  seam. It deliberately preserves existing `Field2D` and `StructuredGrid2D`
  storage surfaces because ndarray/Leto storage replacement is a separate
  provider migration. It does not remove the cfd-2d manifest's direct
  `num-traits` dependency because direct residues remain outside this slice.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-2d --check` passed.
  `cargo check -p cfd-2d --no-default-features` passed. `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings` passed. Focused `cargo
  nextest run -p cfd-2d --no-default-features fvm --status-level fail` passed
  25/25 tests. `git diff --check` passed for the touched FVM and PM artifact
  files. A direct-provider scan over `crates/cfd-2d/src/solvers/fvm` found no
  `num_traits`, `.to_f64()`, `T::from_*`, `T::zero()`, `T::one()`, scalar
  `.abs()`, `nalgebra`, or `ndarray` residue.
- **Residual risk**: Direct `num-traits` remains in cfd-2d problem setup,
  immersed-boundary tests, streamtube partitioning, momentum, energy,
  turbulence, and f64-only/test scalar surfaces; full cfd-2d direct
  `num-traits` removal remains a larger crate-level Eunomia migration before
  the manifest dependency can be dropped.

---

## Sprint 2026-07-04: cfd-2d Vorticity-Stream Eunomia Scalar Seam
- **Resolved**: `crates/cfd-2d/src/physics/vorticity_stream.rs` no longer uses
  direct `num_traits::FromPrimitive`, `T::from_*`, `T::zero()`, `T::one()`, or
  scalar `.abs()` in the touched APIs and tests. Default tolerances,
  timestep/SOR constants, zero/one identities, Laplacian gradient factors, SOR
  residuals, upwind sign checks, velocity-recovery denominators,
  boundary-vorticity factors, convergence checks, and continuity-test
  absolute-value assertions route through `crates/cfd-2d/src/scalar.rs` and
  `eunomia::{FloatElement,NumericElement}`.
- **Boundary**: This slice is limited to the vorticity-stream scalar-provider
  seam. It deliberately preserves existing nalgebra `RealField`, `Vector2`,
  `Array2D`, and `StructuredGrid2D` surfaces because Leto/nalgebra and storage
  replacement are separate provider migrations. It does not remove the cfd-2d
  manifest's direct `num-traits` dependency because direct residues remain
  outside this slice.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-2d --check` passed.
  `cargo check -p cfd-2d --no-default-features` passed. `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings` passed. Focused `cargo
  nextest run -p cfd-2d --no-default-features vorticity --status-level fail`
  passed 4/4 tests. `git diff --check` passed for the touched
  vorticity-stream and PM artifact files. A direct-provider scan over
  `crates/cfd-2d/src/physics/vorticity_stream.rs` found no `num_traits`,
  `FromPrimitive`, `T::from_*`, `T::zero()`, `T::one()`, or scalar `.abs()`
  residue.
- **Residual risk**: Direct `num-traits` remains in cfd-2d problem setup, FVM,
  immersed-boundary tests, streamtube partitioning, momentum, energy,
  turbulence, and f64-only/test scalar surfaces; full cfd-2d direct
  `num-traits` removal remains a larger crate-level Eunomia migration before
  the manifest dependency can be dropped.

---

## Sprint 2026-07-04: cfd-2d Acoustic-Drift Eunomia Scalar Seam
- **Resolved**: `crates/cfd-2d/src/physics/acoustics/gorkov.rs` and
  `crates/cfd-2d/src/solvers/drift_diffusion_2d.rs` no longer use direct
  `num_traits::{Float,FromPrimitive}`, `T::from_*`, `T::zero()`, `T::one()`,
  `Float::`, or scalar `.abs()` in the touched APIs. ARF material constants,
  compressibility identities, contrast-factor constants, standing-wave sine
  evaluation, drift-diffusion zero/one identities, under-relaxation constants,
  Patankar upwind max terms, pivot-threshold checks, and convergence residual
  magnitudes route through `crates/cfd-2d/src/scalar.rs` and
  `eunomia::{FloatElement,NumericElement}`.
- **Boundary**: This slice is limited to the Gor'kov acoustic-force and
  drift-diffusion scalar-provider seam. It deliberately preserves existing
  nalgebra `RealField`, `Array2D`, and NS-FVM field/grid surfaces because
  Leto/nalgebra and storage replacement are separate provider migrations. It
  does not remove the cfd-2d manifest's direct `num-traits` dependency because
  direct residues remain outside this slice.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-2d --check` passed.
  `cargo check -p cfd-2d --no-default-features` passed. `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings` passed. Focused `cargo
  nextest run -p cfd-2d --no-default-features gorkov --status-level fail`
  passed 4/4 tests. Focused `cargo nextest run -p cfd-2d --no-default-features
  drift_diffusion --status-level fail` passed 1/1 test. `git diff --check`
  passed for the touched acoustic-drift and PM artifact files. A direct-provider
  scan over both files found no `num_traits`, `FromPrimitive`, `T::from_*`,
  `T::zero()`, `T::one()`, `Float::`, or scalar `.abs()` residue.
- **Residual risk**: Direct `num-traits` remains in cfd-2d problem setup, FVM,
  vorticity-stream, immersed-boundary tests, streamtube partitioning, momentum,
  energy, turbulence, and f64-only/test scalar surfaces; full cfd-2d direct
  `num-traits` removal remains a larger crate-level Eunomia migration before
  the manifest dependency can be dropped.

---

## Sprint 2026-07-04: cfd-2d TVD Eunomia Scalar Seam
- **Resolved**: `crates/cfd-2d/src/schemes/tvd/{mod.rs,muscl.rs,quick.rs}`
  no longer use direct `num_traits::{FromPrimitive,ToPrimitive}`,
  `T::from_*`, `T::zero()`, `T::one()`, `.to_f64()`, or scalar `.abs()` in the
  touched APIs. TVD limiter identities, min/max/absolute-value operations,
  MUSCL slope ratios, QUICK interpolation constants, and MUSCL2/MUSCL3 boundary
  constants route through `crates/cfd-2d/src/scalar.rs` and
  `eunomia::{FloatElement,NumericElement}`. `FluxLimiter::apply` remains
  generic in `T` instead of round-tripping the limiter ratio through `f64`.
- **Boundary**: This slice is limited to the TVD/MUSCL/QUICK scalar-provider
  seam. It deliberately preserves existing nalgebra `RealField` and `Grid2D`
  surfaces because Leto/nalgebra and storage replacement are separate provider
  migrations. It does not remove the cfd-2d manifest's direct `num-traits`
  dependency because direct residues remain outside this slice.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-2d --check` passed.
  `cargo check -p cfd-2d --no-default-features` passed. `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings` passed. Focused `cargo
  nextest run -p cfd-2d --no-default-features tvd --status-level fail` passed
  32/32 tests. Focused `cargo nextest run -p cfd-2d --no-default-features
  muscl --status-level fail` passed 14/14 tests. `git diff --check` passed for
  the touched TVD and PM artifact files. A direct-provider scan over
  `crates/cfd-2d/src/schemes/tvd` found no `num_traits`, `FromPrimitive`,
  `ToPrimitive`, `T::from_*`, `.to_f64()`, `T::zero()`, `T::one()`, or scalar
  `.abs()` residue.
- **Residual risk**: Direct `num-traits` remains in cfd-2d problem setup,
  drift-diffusion, FVM, vorticity-stream, acoustics, immersed-boundary tests,
  streamtube partitioning, momentum, energy, turbulence, and f64-only/test
  scalar surfaces; full cfd-2d direct `num-traits` removal remains a larger
  crate-level Eunomia migration before the manifest dependency can be dropped.

---

## Sprint 2026-07-04: cfd-2d LBM Eunomia Scalar Seam
- **Resolved**: `crates/cfd-2d/src/solvers/lbm` no longer uses direct
  `num_traits::{Float,FromPrimitive}`, `T::from_*`, `T::zero()`, `T::one()`,
  `Float::`, or scalar `.abs()` in the touched APIs and module tests. Solver
  defaults, low-Mach validation, grid-index conversion, boundary
  density/velocity reconstruction, passive-scalar boundary handling,
  streaming push-zeroing, macroscopic/lattice assertions, Carreau-Yasuda
  assertions, and Shan-Chen pseudopotential/force calculations route through
  `crates/cfd-2d/src/scalar.rs` and
  `eunomia::{FloatElement,NumericElement,CastFrom}`.
- **Boundary**: This slice is limited to the LBM scalar-provider seam. It
  deliberately preserves existing nalgebra `RealField`, `Vector2`, boxed
  collision-operator dispatch, and flat `Vec<T>` distribution storage because
  Leto/nalgebra, execution/backend, and memory-provider replacements are
  separate provider migrations. It does not remove the cfd-2d manifest's direct
  `num-traits` dependency because direct residues remain outside this slice.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-2d --check` passed.
  `cargo check -p cfd-2d --no-default-features` passed. `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings` passed. Focused `cargo
  nextest run -p cfd-2d --no-default-features lbm --status-level fail` passed
  31/31 tests. `git diff --check` passed for the touched LBM and PM artifact
  files. A direct-provider scan over `crates/cfd-2d/src/solvers/lbm` found no
  `num_traits`, `FromPrimitive`, `T::from_*`, `T::zero()`, `T::one()`,
  `Float::`, or scalar `.abs()` residue.
- **Residual risk**: Direct `num-traits` remains in cfd-2d problem setup,
  drift-diffusion, FVM, vorticity-stream, immersed-boundary tests, streamtube
  partitioning, TVD schemes, momentum, energy, turbulence, and f64-only/test
  scalar surfaces; full cfd-2d direct `num-traits` removal remains a larger
  crate-level Eunomia migration before the manifest dependency can be dropped.

---

## Sprint 2026-07-04: cfd-2d SIMPLE Eunomia Scalar Seam
- **Resolved**: `crates/cfd-2d/src/solvers/simple/{algorithm.rs,momentum.rs,pressure.rs}`
  no longer use direct `num_traits::{FromPrimitive,ToPrimitive}`, `T::from_*`,
  `T::zero()`, `T::one()`, scalar `.abs()`, or direct finite checks in the
  touched APIs and module tests. Patankar default constants, workspace zero
  initialization, stagnant-cell safeguards, pressure-correction identities,
  neighbor-count conversion, finite checks, and module test absolute-value
  assertions route through `crates/cfd-2d/src/scalar.rs` and
  `eunomia::{FloatElement,NumericElement}`.
- **Boundary**: This slice is limited to the legacy SIMPLE scalar-provider
  seam. It deliberately preserves existing nalgebra `RealField`, `DVector`,
  sparse-matrix, and momentum-solver boundaries because Leto/nalgebra and
  storage replacement are separate provider migrations. It does not remove the
  cfd-2d manifest's direct `num-traits` dependency because direct residues
  remain outside this slice.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-2d --check` passed.
  `cargo check -p cfd-2d --no-default-features` passed. `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings` passed. Focused `cargo
  nextest run -p cfd-2d --no-default-features simple --status-level fail`
  passed 19/19 tests. `git diff --check` passed for the touched SIMPLE and PM
  artifact files. A direct-provider scan over
  `crates/cfd-2d/src/solvers/simple` found no `num_traits`, `FromPrimitive`,
  `ToPrimitive`, `T::from_*`, `<T as FromPrimitive>`, `.to_f64()`,
  `T::zero()`, `T::one()`, `unwrap_or(T::one())`, `Float::`, or scalar `.abs()`
  residue.
- **Residual risk**: Direct `num-traits` remains in cfd-2d vorticity-stream,
  drift-diffusion, acoustics, momentum, turbulence, energy,
  streamtube/problem/schemes, FVM, and f64-only/test scalar surfaces; full
  cfd-2d direct `num-traits` removal remains a larger crate-level Eunomia
  migration before the manifest dependency can be dropped.

---

## Sprint 2026-07-04: cfd-2d Venturi Eunomia Scalar Seam
- **Resolved**: `crates/cfd-2d/src/solvers/venturi_flow/{mod.rs,solver.rs}`
  no longer use direct `num_traits::{Float,FromPrimitive,ToPrimitive}`,
  `T::from_*`, `.to_f64()`, `T::zero()`, `T::one()`, `Float::`, or scalar
  `.abs()` in the touched APIs and module tests. ISO/default constants, beta
  clamping, domain bounds, analytical Bernoulli/viscous constants,
  energy-dissipation guards, stretched-grid index conversion, generic
  sine/min/max dispatch, diagnostic f64 conversion, throat-column selection,
  inlet/outlet/throat reductions, and validation error metrics route through
  `crates/cfd-2d/src/scalar.rs` and
  `eunomia::{FloatElement,NumericElement}`.
- **Boundary**: This slice is limited to the Venturi scalar-provider seam. It
  deliberately preserves existing nalgebra `RealField` and NS-FVM storage/
  solver boundaries because Leto/nalgebra and storage replacement are separate
  provider migrations. It does not remove the cfd-2d manifest's direct
  `num-traits` dependency because direct residues remain outside this slice.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-2d --check` passed.
  `cargo check -p cfd-2d --no-default-features` passed. `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings` passed. Focused `cargo
  nextest run -p cfd-2d --no-default-features venturi --status-level fail`
  passed 13/13 tests. `git diff --check` passed for the touched Venturi and PM
  artifact files. A direct-provider scan over
  `crates/cfd-2d/src/solvers/venturi_flow` found no `num_traits`,
  `FromPrimitive`, `ToPrimitive`, `T::from_*`, `<T as FromPrimitive>`,
  `.to_f64()`, `T::zero()`, `T::one()`, `unwrap_or(T::one())`, `Float::`, or
  scalar `.abs()` residue.
- **Residual risk**: Direct `num-traits` remains in cfd-2d vorticity-stream,
  LBM, drift-diffusion, acoustics, momentum, turbulence, energy,
  streamtube/problem/schemes, FVM, and f64-only/test scalar surfaces; full
  cfd-2d direct `num-traits` removal remains a larger crate-level Eunomia
  migration before the manifest dependency can be dropped.

---

## Sprint 2026-07-04: cfd-2d Branching-Flow Eunomia Scalar Seam
- **Resolved**: `crates/cfd-2d/src/solvers/bifurcation_flow.rs` and
  `crates/cfd-2d/src/solvers/n_furcation_flow.rs` no longer use direct
  `num_traits::{Float,FromPrimitive,ToPrimitive}`, `T::from_*`, `T::zero()`,
  or `Float::` in the touched APIs and module tests. Geometry constants,
  branch-index conversion, zero identities, generic sin/cos dispatch, flux
  accumulation, and mass-balance absolute-value normalization route through
  `crates/cfd-2d/src/scalar.rs` and
  `eunomia::{FloatElement,NumericElement}`.
- **Boundary**: This slice is limited to the bifurcation/N-furcation
  scalar-provider seam. It deliberately preserves existing nalgebra
  `RealField` and NS-FVM storage/solver boundaries because Leto/nalgebra and
  storage replacement are separate provider migrations. It does not remove the
  cfd-2d manifest's direct `num-traits` dependency because direct residues
  remain outside this slice.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-2d --check` passed.
  `cargo check -p cfd-2d --no-default-features` passed. `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings` passed. Focused `cargo
  nextest run -p cfd-2d --no-default-features bifurcation --status-level fail`
  passed 9/9 tests. Focused `cargo nextest run -p cfd-2d --no-default-features
  n_furcation --status-level fail` passed 1/1 tests. `git diff --check` passed
  for the touched branching and PM artifact files. A direct-provider scan over
  both branching files found no `num_traits`, `FromPrimitive`, `ToPrimitive`,
  `T::from_*`, `<T as FromPrimitive>`, `.to_f64()`, `T::zero()`, `T::one()`,
  `unwrap_or(T::one())`, or `Float::` residue.
- **Residual risk**: Direct `num-traits` remains in cfd-2d vorticity-stream,
  LBM, drift-diffusion, acoustics, momentum, turbulence,
  energy, streamtube/problem/schemes, FVM, and f64-only/test scalar surfaces;
  full cfd-2d direct `num-traits` removal remains a larger crate-level Eunomia
  migration before the manifest dependency can be dropped.

---

## Sprint 2026-07-04: cfd-2d Cross-Junction Eunomia Scalar Seam
- **Resolved**: `crates/cfd-2d/src/solvers/cross_junction_flow.rs` no longer
  uses direct `num_traits::{Float,FromPrimitive,ToPrimitive}`, `T::from_*`,
  `T::zero()`, or scalar `.abs()` in the touched APIs and module tests.
  Geometry constants, zero identities for port flux accumulation, mass-balance
  absolute-value normalization, and bounding-box test checks route through
  `crates/cfd-2d/src/scalar.rs` and
  `eunomia::{FloatElement,NumericElement}`.
- **Boundary**: This slice is limited to the cross-junction scalar-provider
  seam. It deliberately preserves existing nalgebra `RealField` boundaries
  because Leto/nalgebra replacement is a separate provider migration. It does
  not remove the cfd-2d manifest's direct `num-traits` dependency because
  direct residues remain outside this slice.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-2d --check` passed.
  `cargo check -p cfd-2d --no-default-features` passed. `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings` passed. Focused `cargo
  nextest run -p cfd-2d --no-default-features cross_junction --status-level
  fail` passed 5/5 tests. `git diff --check` passed for the touched
  cross-junction and PM artifact files. A direct-provider scan over
  `crates/cfd-2d/src/solvers/cross_junction_flow.rs` found no `num_traits`,
  `FromPrimitive`, `ToPrimitive`, `T::from_*`, `<T as FromPrimitive>`,
  `.to_f64()`, `T::zero()`, `T::one()`, `unwrap_or(T::one())`, `Float::`, or
  scalar `.abs()` residue.
- **Residual risk**: Direct `num-traits` remains in cfd-2d LBM,
  drift-diffusion, and f64-only/test
  scalar surfaces; full cfd-2d direct `num-traits` removal remains a larger
  crate-level Eunomia migration before the manifest dependency can be dropped.

---

## Sprint 2026-07-04: cfd-2d Poiseuille Eunomia Scalar Seam
- **Resolved**: `crates/cfd-2d/src/solvers/poiseuille/{mod.rs,numerics.rs}`
  no longer use direct `num_traits::{Float,FromPrimitive}`, `T::from_*`,
  `T::zero()`, `T::one()`, or scalar `.abs()` in the touched APIs and module
  tests. Configuration defaults, grid-index conversion, zero/one identities,
  tridiagonal workspaces, harmonic-mean constants, shear-rate magnitudes,
  viscosity residual norms, Thomas pivot thresholds, error diagnostics, and
  test absolute-value checks route through `crates/cfd-2d/src/scalar.rs` and
  `eunomia::{FloatElement,NumericElement}`.
- **Boundary**: This slice is limited to the Poiseuille scalar-provider seam.
  It deliberately preserves existing nalgebra `RealField` boundaries because
  Leto/nalgebra replacement is a separate provider migration. It does not
  remove the cfd-2d manifest's direct `num-traits` dependency because direct
  residues remain outside this slice.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-2d --check` passed.
  `cargo check -p cfd-2d --no-default-features` passed. `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings` passed. Focused `cargo
  nextest run -p cfd-2d --no-default-features poiseuille --status-level fail`
  passed 7/7 tests. `git diff --check` passed for the touched Poiseuille
  files. A direct-provider scan over `crates/cfd-2d/src/solvers/poiseuille`
  found no `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_*`, `<T as
  FromPrimitive>`, `.to_f64()`, `T::zero()`, `T::one()`, `unwrap_or(T::one())`,
  `Float::`, `f64::abs`, or scalar `.abs()` residue.
- **Residual risk**: Direct `num-traits` remains in cfd-2d vorticity-stream,
  LBM solvers, acoustics, momentum, turbulence,
  physics, and tests; full cfd-2d direct `num-traits` removal remains a larger
  crate-level Eunomia migration before the manifest dependency can be dropped.

---

## Sprint 2026-07-04: cfd-2d Serpentine/Scalar Transport Eunomia Scalar Seam
- **Resolved**: `crates/cfd-2d/src/solvers/serpentine_flow/{mod.rs,solver.rs}`
  and `crates/cfd-2d/src/solvers/scalar_transport_2d.rs` no longer use direct
  `num_traits::{Float,FromPrimitive,ToPrimitive}`, `T::from_*`, `T::zero()`,
  `T::one()`, or scalar `.abs()` in the touched APIs and module tests. Scalar
  constants, identities, Peclet and mixing calculations, Fourier-mode index
  conversion, transport relaxation constants, coefficient max/abs operations,
  outlet averaging, validator diagnostics, and test absolute-value checks route
  through `crates/cfd-2d/src/scalar.rs` and
  `eunomia::{FloatElement,NumericElement}`.
- **Boundary**: This slice is limited to the serpentine/scalar-transport
  scalar-provider seam. It deliberately preserves existing nalgebra
  `RealField` and grid/field storage surfaces because Leto/nalgebra and
  ndarray replacement are separate provider migrations. It does not remove the
  cfd-2d manifest's direct `num-traits` dependency because direct residues
  remain outside this slice.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-2d --check` passed.
  `cargo check -p cfd-2d --no-default-features` passed. `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings` passed. Focused `cargo
  nextest run -p cfd-2d --no-default-features serpentine --status-level fail`
  passed 5/5 tests. A `scalar_transport` nextest filter found no matching
  tests, so transport-specific evidence is compile/clippy plus serpentine
  integration coverage. `git diff --check` passed for the touched
  serpentine/transport files. A direct-provider scan over
  `crates/cfd-2d/src/solvers/serpentine_flow` and
  `crates/cfd-2d/src/solvers/scalar_transport_2d.rs` found no `num_traits`,
  `FromPrimitive`, `ToPrimitive`, `T::from_*`, `<T as FromPrimitive>`,
  `.to_f64()`, `T::zero()`, `T::one()`, `unwrap_or(T::one())`, `Float::`, or
  scalar `.abs()` residue.
- **Residual risk**: Direct `num-traits` remains in cfd-2d vorticity-stream,
  LBM solvers, acoustics, momentum,
  turbulence, physics, and tests; full cfd-2d direct `num-traits` removal
  remains a larger crate-level Eunomia migration before the manifest dependency
  can be dropped.

---

## Sprint 2026-07-04: cfd-2d FDM Eunomia Scalar Seam
- **Resolved**: `crates/cfd-2d/src/solvers/fdm/{advection_diffusion.rs,config.rs,diffusion.rs,linear_solver.rs,poisson.rs,mod.rs,tests_advection_diffusion_mms.rs,tests_poisson_mms.rs}`
  no longer use direct `num_traits::FromPrimitive`, `T::from_*`, `T::zero()`,
  `T::one()`, or scalar `.abs()` in the touched FDM APIs and MMS checks.
  Scalar constants, identities, singular-diagonal guards, Gauss-Seidel
  residual magnitudes, finite-difference constants, grid-index conversions, and
  MMS absolute-value checks route through `crates/cfd-2d/src/scalar.rs` and
  `eunomia::{FloatElement,NumericElement}`.
- **Correctness finding**: Wiring the existing FDM MMS source tests into
  `mod.rs` exposed an advection-diffusion sign defect. The stencil now assembles
  the documented steady operator `u dot grad(phi) - alpha laplacian(phi) = S`
  instead of its negation.
- **Boundary**: This slice is limited to the FDM scalar-provider seam and the
  directly exposed FDM MMS correctness defect. It deliberately preserves
  existing nalgebra `RealField`, `DVector`, and sparse-matrix surfaces because
  Leto/nalgebra replacement is a separate provider migration. It does not
  remove the cfd-2d manifest's direct `num-traits` dependency because direct
  residues remain outside this slice.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-2d --check` passed.
  `cargo check -p cfd-2d --no-default-features` passed. `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings` passed. Focused `cargo
  nextest run -p cfd-2d --no-default-features fdm --status-level fail` passed
  2/2 tests. `git diff --check` passed for the touched FDM files. A
  direct-provider scan over `crates/cfd-2d/src/solvers/fdm` found no
  `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_*`, `<T as
  FromPrimitive>`, `.to_f64()`, `T::zero()`, `T::one()`, `unwrap_or(T::one())`,
  `Float::`, or scalar `.abs()` residue.
- **Residual risk**: Direct `num-traits` remains in cfd-2d vorticity-stream,
  LBM solvers, acoustics, momentum,
  turbulence, physics, and tests; full cfd-2d direct `num-traits` removal
  remains a larger crate-level Eunomia migration before the manifest dependency
  can be dropped.

---

## Sprint 2026-07-04: cfd-2d Network Eunomia Scalar Seam
- **Resolved**: `crates/cfd-2d/src/network/{build.rs,channel.rs,coupled.rs,postprocess.rs,projection.rs,reference.rs,solve.rs,types.rs}`
  no longer use direct `num_traits::{Float,FromPrimitive,ToPrimitive}`,
  `T::from_*`, `.to_f64()`, `T::zero()`, or `T::one()` in the touched network
  APIs. Scalar constants, identities, absolute values, min/max operations,
  finite checks, Anderson relaxation constants, projection summaries, channel
  diagnostics, and f64 reporting route through `crates/cfd-2d/src/scalar.rs`
  and `eunomia::{FloatElement,NumericElement}`.
- **Boundary**: This slice is limited to the network scalar-provider seam. It
  deliberately preserves existing nalgebra `RealField` and ndarray storage
  surfaces because Leto/nalgebra and ndarray replacement are separate provider
  migrations. It does not remove the cfd-2d manifest's direct `num-traits`
  dependency because direct residues remain outside this slice.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-2d --check` passed.
  `cargo check -p cfd-2d --no-default-features` passed. `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings` passed. Focused `cargo
  nextest run -p cfd-2d --no-default-features network --status-level fail`
  passed 22/22 tests. `git diff --check` passed for the touched network files.
  A direct-provider scan over `crates/cfd-2d/src/network` found no
  `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_*`, `<T as
  FromPrimitive>`, `.to_f64()`, `T::zero()`, `T::one()`, `unwrap_or(T::one())`,
  or `Float::` residue.
- **Residual risk**: Direct `num-traits` remains in cfd-2d drift-diffusion,
  LBM/FVM, turbulence, physics, and other solver/test surfaces; full
  cfd-2d direct `num-traits`
  removal remains a larger crate-level Eunomia migration before the manifest
  dependency can be dropped.

---

## Sprint 2026-07-04: cfd-2d NS-FVM Eunomia Scalar Seam
- **Resolved**: `crates/cfd-2d/src/solvers/ns_fvm/{field.rs,solver/mod.rs}`,
  `solver/momentum/{u_equation.rs,v_equation.rs}`,
  `solver/pressure/{poisson.rs,rhie_chow.rs}`, and
  `solver/velocity_interpolation.rs` no longer use direct
  `num_traits::{Float,FromPrimitive,ToPrimitive}`, `T::from_*`, `.to_f64()`,
  `T::zero()`, or `T::one()` in the touched NS-FVM APIs. Scalar constants,
  identities, absolute values, min/max/sqrt/finite checks, pressure-correction
  thresholds, turbulence initialization constants, velocity-interpolation
  grid-index conversion, and f64 diagnostics route through
  `crates/cfd-2d/src/scalar.rs` and `eunomia::{FloatElement,NumericElement}`.
- **Boundary**: This slice is limited to the NS-FVM scalar-provider seam. It
  deliberately preserves the existing nalgebra `RealField` and vector surfaces
  because Leto/nalgebra replacement is a separate provider migration. It does
  not remove the cfd-2d manifest's direct `num-traits` dependency because
  direct residues remain outside this slice.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-2d --check` passed.
  `cargo check -p cfd-2d --no-default-features` passed. `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings` passed. Focused `cargo
  nextest run -p cfd-2d --no-default-features ns_fvm --status-level fail`
  passed 2/2 tests. `git diff --check` passed for the touched NS-FVM files. A
  direct-provider scan over `crates/cfd-2d/src/solvers/ns_fvm` found no
  `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_*`, `<T as
  FromPrimitive>`, `.to_f64()`, `T::zero()`, `T::one()`, `unwrap_or(T::one())`,
  or `Float::` residue.
- **Residual risk**: Direct `num-traits` remains in cfd-2d network,
  turbulence, LBM/FDM/FVM, solver, physics, and other test surfaces; full
  cfd-2d direct `num-traits` removal remains a larger crate-level Eunomia
  migration before the manifest dependency can be dropped.

---

## Sprint 2026-07-04: cfd-2d Continuity Eunomia Scalar Seam
- **Resolved**: `crates/cfd-2d/src/solvers/continuity.rs` no longer uses
  direct `num_traits::FromPrimitive`, `T::from_*`, `T::zero()`, or scalar
  method `.abs()` in the shared forward, central, and face continuity-residual
  helpers. Zero identities and the central-difference factor route through
  `crates/cfd-2d/src/scalar.rs`, backed by Eunomia. Absolute-value reductions
  use `eunomia::NumericElement`.
- **Boundary**: This slice is limited to the continuity residual SSOT. It does
  not remove the cfd-2d manifest's direct `num-traits` dependency because
  direct residues remain in other cfd-2d modules. A network scalar cleanup was
  inspected but deferred because `Network2DSolver` stores
  `NavierStokesSolver2D<T>`, whose type definition still owns direct
  `num_traits::{Float,FromPrimitive}` bounds; the NS-FVM scalar-bound migration
  must precede a clean network signature migration.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-2d --check` passed.
  `cargo check -p cfd-2d --no-default-features` passed. `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings` passed. Focused `cargo
  nextest run -p cfd-2d --no-default-features continuity --status-level fail`
  passed 6/6 tests. `git diff --check` passed for the touched file. A
  direct-provider scan over `crates/cfd-2d/src/solvers/continuity.rs` found no
  `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_*`, `<T as
  FromPrimitive>`, `.to_f64()`, `T::zero()`, `T::one()`, or method `.abs()`
  residue.
- **Residual risk**: Direct `num-traits` remains in cfd-2d network,
  turbulence, LBM/FDM/FVM/NS-FVM, solver, physics, and other test surfaces;
  full cfd-2d direct `num-traits` removal remains a larger crate-level
  Eunomia migration before the manifest dependency can be dropped.

---

## Sprint 2026-07-04: cfd-2d SIMPLEC/PIMPLE Eunomia Scalar Seam
- **Resolved**: `crates/cfd-2d/src/simplec_pimple/{config,algorithms,diagnostics,solver,interpolation,simplec,pimple}.rs`
  and `crates/cfd-2d/tests/simplec_pimple_validation.rs` no longer use direct
  `num_traits::{FromPrimitive,ToPrimitive}`, `T::from_*`, `.to_f64()`,
  `T::zero()`, or `T::one()` in the touched SIMPLEC/PIMPLE APIs and focused
  validation helpers. Defaults, validation identities, residual computations,
  adaptive-step constants, Rhie-Chow face caches, pressure extrapolation,
  boundary zeroing, PIMPLE unity relaxation, and Ghia reference interpolation
  route through `crates/cfd-2d/src/scalar.rs`, backed by Eunomia, or local
  Eunomia-backed test helpers.
- **Boundary**: This slice is limited to the SIMPLEC/PIMPLE module and its
  focused validation test. It preserves existing nalgebra `RealField`/`Vector2`
  surfaces because this increment is Eunomia scalar-provider replacement, not
  the later Leto/nalgebra replacement. It does not remove the crate manifest's
  direct `num-traits` dependency because direct residues remain in other
  cfd-2d modules.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-2d --check` passed.
  `cargo check -p cfd-2d --no-default-features` passed. `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings` passed. Focused `cargo
  nextest run -p cfd-2d --no-default-features simplec_pimple --status-level
  fail` passed 4/4 tests. Direct-provider scans over the touched module/test
  found no `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_*`,
  `<T as FromPrimitive>`, `.to_f64()`, `T::zero()`, `T::one()`, or
  `unwrap_or(T::one())` residue.
- **Residual risk**: Direct `num-traits` remains in cfd-2d network,
  turbulence, LBM/FDM/FVM/NS-FVM, solver, physics, and other test surfaces;
  full cfd-2d direct `num-traits` removal remains a larger crate-level
  Eunomia migration before the manifest dependency can be dropped.

---

## Sprint 2026-07-04: cfd-2d PISO Eunomia Scalar Seam
- **Resolved**: `crates/cfd-2d/src/piso_algorithm/{config,convergence,predictor,corrector,solver}.rs`
  no longer use direct `num_traits::{FromPrimitive,ToPrimitive}`,
  `T::from_*`, `.to_f64()`, `T::zero()`, or `T::one()` in the touched PISO
  APIs. Defaults, tolerances, residual monitors, predictor/corrector
  workspaces, relaxation identities, hybrid differencing, Rhie-Chow face-flux
  tiny thresholds, duration thresholds, and logging conversions route through
  `crates/cfd-2d/src/scalar.rs`, backed by Eunomia, and explicit
  `eunomia::{FloatElement,NumericElement}` dispatch.
- **Boundary**: This slice is limited to the PISO module. It preserves the
  existing nalgebra `RealField`/`Vector2` surfaces because this increment is
  Eunomia scalar-provider replacement, not the later Leto/nalgebra
  replacement. It does not remove the crate manifest's direct `num-traits`
  dependency because direct residues remain in other cfd-2d modules.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-2d --check` passed.
  `cargo check -p cfd-2d --no-default-features` passed. `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings` passed. Focused `cargo
  nextest run -p cfd-2d --no-default-features piso --status-level fail`
  passed 6/6 tests. A direct-provider scan over
  `crates/cfd-2d/src/piso_algorithm` found no `num_traits`, `FromPrimitive`,
  `ToPrimitive`, `T::from_*`, `<T as FromPrimitive>`, `.to_f64()`,
  `T::zero()`, `T::one()`, or `unwrap_or(T::one())` residue.
- **Residual risk**: Direct `num-traits` remains in cfd-2d network,
  turbulence, LBM/FDM/FVM/NS-FVM, solver, physics, and test surfaces; full
  cfd-2d direct `num-traits` removal remains a larger crate-level Eunomia
  migration before the manifest dependency can be dropped.

---

## Sprint 2026-07-04: cfd-2d Pressure-Velocity Eunomia Scalar Seam
- **Resolved**: `crates/cfd-2d/src/pressure_velocity/{coefficients,config,solver,rhie_chow,faces,correction,pressure}.rs`
  no longer use direct `T::zero()`/`T::one()` scalar construction in the
  touched pressure-velocity APIs. Coefficient defaults, validation bounds,
  pressure-correction assembly and scatter, solver workspaces, and Rhie-Chow
  coefficient guards route identities through `crates/cfd-2d/src/scalar.rs`,
  backed by Eunomia. Finite checks in the touched solver/Rhie-Chow paths use
  `eunomia::NumericElement`.
- **Boundary**: This slice is limited to the `pressure_velocity` module's
  scalar-provider residue. It preserves the existing nalgebra
  `RealField`/`Vector2`/`DVector` surfaces because this increment is Eunomia
  scalar-provider replacement, not the later Leto/nalgebra replacement. It
  does not remove the crate manifest's direct `num-traits` dependency because
  direct residues remain in other cfd-2d modules.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-2d --check` passed.
  `cargo check -p cfd-2d --no-default-features` passed. `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings` passed. Focused `cargo
  nextest run -p cfd-2d --no-default-features pressure_velocity
  --status-level fail` passed 16/16 tests. A direct-provider scan over
  `crates/cfd-2d/src/pressure_velocity` found no `num_traits`,
  `FromPrimitive`, `ToPrimitive`, `T::from_*`, `.to_f64()`, `T::zero()`, or
  `T::one()` residue.
- **Residual risk**: Direct `num-traits` remains in cfd-2d SIMPLEC/PIMPLE,
  network, turbulence, LBM/FDM/FVM/NS-FVM, solver, physics, and test surfaces;
  full cfd-2d direct `num-traits` removal remains a larger crate-level
  Eunomia migration before the manifest dependency can be dropped.

---

## Sprint 2026-07-04: cfd-2d Time Integration Eunomia Scalar Seam
- **Resolved**: `crates/cfd-2d/src/schemes/time/{mod,explicit,implicit,multistep}.rs`
  and `crates/cfd-2d/src/schemes/time/adaptive/{mod,integrator}.rs` no longer
  use direct `num_traits::{FromPrimitive,ToPrimitive}` bounds, `T::from_f64`,
  `T::zero()`, `T::one()`, or `unwrap_or(T::one())` in the touched
  time-integration APIs. RK, implicit fixed-point, BDF/Adams-Bashforth, CFL,
  Richardson-error, and adaptive-step constants route through
  `crates/cfd-2d/src/scalar.rs`, backed by Eunomia, and abs/min/max/powf
  dispatch is explicit through `eunomia::{NumericElement,FloatElement}`.
- **Boundary**: This slice is limited to the cfd-2d time-integration scheme
  surface. It preserves the existing nalgebra `DVector`/`RealField` API because
  this increment is Eunomia scalar-provider replacement, not the later
  Leto/nalgebra replacement. It does not remove the crate manifest's direct
  `num-traits` dependency because direct residues remain in other cfd-2d
  modules.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-2d --check` passed.
  `cargo check -p cfd-2d --no-default-features` passed. `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings` passed. Focused `cargo
  nextest run -p cfd-2d --no-default-features time --status-level fail`
  passed 29/29 tests. A direct-provider scan over
  `crates/cfd-2d/src/schemes/time` found no `num_traits`, `FromPrimitive`,
  `ToPrimitive`, `T::from_*`, `<T as FromPrimitive>`, `.to_f64()`,
  `T::zero()`, `T::one()`, or `unwrap_or(T::one())` residue.
- **Residual risk**: Direct `num-traits` remains in cfd-2d SIMPLEC/PIMPLE,
  network, turbulence, LBM/FDM/FVM/NS-FVM, solver, physics, and test surfaces;
  full cfd-2d direct `num-traits` removal remains a larger crate-level Eunomia
  migration before the manifest dependency can be dropped.

---

## Sprint 2026-07-04: cfd-2d Discretization Eunomia Scalar Seam
- **Resolved**: `crates/cfd-2d/src/discretization/{convection,extended_stencil}.rs`
  and
  `crates/cfd-2d/src/physics/momentum/coefficient_corrections/quick.rs` no
  longer use direct `num_traits::FromPrimitive` scalar construction or
  `T::zero()`/`T::one()` constructors in the touched first-order upwind,
  central, hybrid, power-law, QUICK, and MUSCL paths. Scalar constants,
  zero/one values, absolute values, and max operations route through
  `crates/cfd-2d/src/scalar.rs`, backed by Eunomia.
- **Boundary**: This slice is limited to the finite-volume discretization
  strategy surface and its QUICK momentum correction caller. It preserves the
  existing `Box<dyn ConvectionScheme<T>>` factory shape because that dynamic
  dispatch already existed and the migration target here is scalar-provider
  replacement, not API redesign. It does not remove the crate manifest's
  direct `num-traits` dependency because direct residues remain in other
  cfd-2d modules.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-2d --check` passed.
  `cargo check -p cfd-2d --no-default-features` passed. `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings` passed. Focused `cargo
  nextest run -p cfd-2d --no-default-features discretization --status-level
  fail` passed 8/8 tests. A direct-provider scan over the touched files found
  no `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_*`,
  `<T as FromPrimitive>`, `.to_f64()`, `T::zero()`, or `T::one()` residue.
- **Residual risk**: Direct `num-traits` remains in cfd-2d network, PISO,
  LBM, turbulence, NS-FVM, scalar transport, solver, and test surfaces; full
  cfd-2d direct `num-traits` removal remains a larger crate-level Eunomia
  migration before the manifest dependency can be dropped.

---

## Sprint 2026-07-04: cfd-2d Spatial/WENO Eunomia Scalar Seam
- **Resolved**: `crates/cfd-2d/src/schemes/{constants,central,grid,upwind,weno_helpers,weno,weno_z}.rs`,
  `crates/cfd-2d/src/schemes/tvd/{mod,quick}.rs`,
  `crates/cfd-2d/src/physics/momentum/coefficients.rs`, and
  `crates/cfd-2d/src/physics/momentum/coefficient_corrections/weno_z.rs`
  no longer use direct `num_traits::{FromPrimitive,ToPrimitive}` scalar
  construction in the touched central, upwind, TVD/QUICK, WENO5/WENO9, and
  WENO-Z momentum correction paths. Scalar constants, zero values, WENO
  `powi` calls, WENO-Z absolute values, deferred-correction relaxation
  factors, and limiter constants route through `crates/cfd-2d/src/scalar.rs`,
  backed by Eunomia.
- **Boundary**: This slice is limited to the spatial/WENO scheme and WENO-Z
  momentum coefficient surface already touched by the cfd-2d Eunomia
  migration. It does not remove the crate manifest's direct `num-traits`
  dependency because direct residues remain in other cfd-2d modules.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-2d --check` passed.
  `cargo check -p cfd-2d --no-default-features` passed. `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings` passed. Focused `cargo
  nextest run -p cfd-2d --no-default-features weno --status-level fail`
  passed 5/5 tests. A direct-provider scan over the touched files found no
  `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_*`,
  `<T as FromPrimitive>`, or `.to_f64()` residue.
- **Residual risk**: Direct `num-traits` remains in cfd-2d discretization,
  PISO, network, LBM, turbulence, NS-FVM, solver, and test surfaces; full
  cfd-2d direct `num-traits` removal remains a larger crate-level Eunomia
  migration before the manifest dependency can be dropped.

---

## Sprint 2026-07-04: cfd-2d Field Container Eunomia Scalar Seam
- **Resolved**: `crates/cfd-2d/src/fields.rs` no longer imports direct
  `num_traits` or uses direct `FromPrimitive`, `Float`, `T::from_*`,
  `T::zero()`, or `T::one()` for field constants, `Field2D::zeros`,
  `SimulationFields` construction/reset, velocity-magnitude reductions, or
  Reynolds-number averaging. Scalar construction and identity values now route
  through `crates/cfd-2d/src/scalar.rs`, backed by Eunomia.
- **Boundary**: This slice is limited to the cfd-2d field-container surface.
  It deliberately keeps `copy_from`, velocity accessors, and zero-copy field
  views under the narrower `RealField + Copy` bound because those methods do
  not construct scalar constants. It does not remove the crate manifest's
  direct `num-traits` dependency because broad direct residues remain outside
  `fields.rs`.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-2d --check` passed.
  `cargo check -p cfd-2d --no-default-features` passed. `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings` passed. Focused `cargo
  nextest run -p cfd-2d --no-default-features piso_algorithm --status-level
  fail` passed 6/6 tests. A direct-provider scan over
  `crates/cfd-2d/src/fields.rs` found no `num_traits`, `FromPrimitive`,
  `ToPrimitive`, `Float::`, `T::from_*`, `T::zero()`, `T::one()`, or
  `.to_f64()` residue.
- **Residual risk**: Direct `num-traits` remains in cfd-2d schemes, network,
  PISO, LBM, turbulence, integration tests, and other modules; full cfd-2d
  direct `num-traits` removal remains a larger crate-level Eunomia migration.

---

## Sprint 2026-07-04: cfd-3d Direct num-traits Ownership Removal
- **Resolved**: `crates/cfd-3d/src/lib.rs` no longer calls
  `num_traits::Float::abs` in root Chebyshev assertions, and
  `crates/cfd-3d/Cargo.toml` no longer declares `num-traits`.
- **Boundary**: This closes direct `num-traits` ownership for `cfd-3d`. It
  does not claim full provider graph removal because `num-traits` still enters
  transitively through upstream/provider paths such as `approx`, `nalgebra`,
  `gaia`, `half`, and `num-complex`.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-3d --check` passed.
  `cargo check -p cfd-3d --no-default-features` passed. Focused `cargo
  nextest run -p cfd-3d --no-default-features chebyshev --status-level fail`
  passed 20/20 tests. `cargo clippy -p cfd-3d --no-default-features --lib --
  -D warnings` passed. A direct-provider scan over
  `crates/cfd-3d/{src,tests,examples}` plus `crates/cfd-3d/Cargo.toml` found
  no direct `num_traits`, `num-traits`, `FromPrimitive`, `ToPrimitive`,
  `num_traits::Float`, `Float::`, `from_f64_or_one`,
  `<T as From<f64>>::from`, `T::zero()`, `T::one()`, or `T::from_` residue.
  `cargo tree -p cfd-3d -e normal -i num-traits` still reports transitive
  provider/nalgebra paths.
- **Residual risk**: Remaining cfd-3d Atlas migration work is nalgebra
  geometry/vector/matrix replacement with Leto/Gaia and upstream transitive
  provider graph cleanup, not direct cfd-3d `num-traits` ownership.

---

## Sprint 2026-07-04: cfd-3d Bifurcation/Venturi Eunomia Scalar Seam
- **Resolved**: `crates/cfd-3d/src/bifurcation/{analysis,geometry,solver,types,validation}.rs`
  and `crates/cfd-3d/src/venturi/{analysis,solver,types,validation}.rs` no
  longer import or bound direct
  `num_traits::{Float,FromPrimitive,ToPrimitive}` and no longer call
  `num_traits::Float::*`, `<T as From<f64>>::from`, `T::zero()`, `T::one()`,
  `T::from_*`, or `T::from_f64_or_one` in the touched module APIs.
  Constants, zero/one values, abs/min/max/sqrt/ln/powf/powi/cos, scalar
  construction, pressure coefficients, Richardson/GCI metrics, Picard
  viscosity-change checks, pressure slicing, flow extraction, wall-shear
  estimates, and validation thresholds now route through `cfd-3d::scalar`,
  backed by Eunomia.
- **Boundary**: This is the Bifurcation/Venturi scalar-provider cleanup. It
  preserves current nalgebra geometry/vector/matrix and FEM boundaries because
  those belong to the later Leto/Gaia storage and geometry migration. It also
  preserves the existing `SafeFromF64` fluid-trait bound where the upstream
  fluid contract still requires it.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-3d --check` passed.
  `cargo check -p cfd-3d --no-default-features` passed. Focused `cargo
  nextest run -p cfd-3d --no-default-features bifurcation venturi
  --status-level fail` passed 35/35 tests; nextest marked one
  mesh-convergence test slow at 16.7s, below the 30s defect threshold in the
  active gate. `cargo clippy -p cfd-3d --no-default-features --lib -- -D
  warnings` passed. A targeted scan over `crates/cfd-3d/src/bifurcation`,
  `crates/cfd-3d/src/venturi`, and `crates/cfd-3d/src/scalar.rs` found no
  direct `num_traits`, `FromPrimitive`, `ToPrimitive`, `num_traits::Float`,
  `Float::`, `from_f64_or_one`, `<T as From<f64>>::from`, `T::zero()`,
  `T::one()`, or `T::from_` residue.
- **Residual risk**: Direct cfd-3d `num-traits` ownership is now closed. Full
  Bifurcation/Venturi provider completion still requires nalgebra
  geometry/vector/matrix and FEM matrix/vector replacement with the Atlas
  geometry/storage stack.

---

## Sprint 2026-07-04: cfd-3d Trifurcation Eunomia Scalar Seam
- **Resolved**: `crates/cfd-3d/src/trifurcation/{geometry,solver,validation}.rs`
  no longer import or bound direct
  `num_traits::{Float,FromPrimitive,ToPrimitive}` and no longer call
  `num_traits::Float::*`, `<T as From<f64>>::from`, `T::zero()`, `T::one()`,
  or `T::from_f64_or_one` in the touched Trifurcation geometry, solver, and
  validation paths. SDF trigonometry, Murray-law powers, solver defaults,
  Picard viscosity-change checks, boundary-flow extraction, pressure drops,
  wall-shear estimates, mass conservation, and validation thresholds now route
  through `cfd-3d::scalar`, backed by Eunomia.
- **Boundary**: This is the Trifurcation scalar-provider cleanup. It preserves
  the current nalgebra `Point3`/`Vector3`/matrix and FEM boundaries because
  those belong to the later Leto/Gaia storage and geometry migration. It also
  preserves the existing `SafeFromF64` fluid-trait bound where the upstream
  fluid contract still requires it.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-3d` completed. `cargo
  check -p cfd-3d --no-default-features` passed. Focused `cargo nextest run
  -p cfd-3d --no-default-features trifurcation --status-level fail` passed
  6/6 tests. `cargo clippy -p cfd-3d --no-default-features --lib -- -D
  warnings` passed. A targeted scan over `crates/cfd-3d/src/trifurcation` and
  `crates/cfd-3d/src/scalar.rs` found no direct `num_traits`,
  `FromPrimitive`, `ToPrimitive`, `num_traits::Float`, `Float::`,
  `from_f64_or_one`, `<T as From<f64>>::from`, `T::zero()`, or `T::one()`
  residue.
- **Residual risk**: Direct cfd-3d `num-traits` ownership is now closed. Full
  Trifurcation provider completion still requires nalgebra geometry/vector and
  FEM matrix/vector replacement with the Atlas geometry/storage stack.

---

## Sprint 2026-07-04: cfd-3d Serpentine Eunomia Scalar Seam
- **Resolved**: `crates/cfd-3d/src/serpentine/{solver,validation}.rs` no
  longer import or bound direct `num_traits::{Float,FromPrimitive,ToPrimitive}`
  and no longer call `num_traits::Float::*`, `T::zero()`, `T::one()`, or
  `T::from_f64_or_one` for Serpentine configuration defaults, Dean-number
  analysis, pressure-drop validation, Picard viscosity-change checks, or
  solution defaults. Constants and elementary scalar operations now route
  through `cfd-3d::scalar`, backed by Eunomia.
- **Boundary**: This is the Serpentine scalar-provider cleanup. It preserves
  the current nalgebra `Vector3`/matrix geometry and FEM boundaries because
  those belong to the later Leto/Gaia storage and vector migration. It also
  preserves the existing `SafeFromF64` fluid-trait bound where the upstream
  fluid contract still requires it.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-3d` completed. `cargo
  check -p cfd-3d --no-default-features` passed. Focused `cargo nextest run
  -p cfd-3d --no-default-features serpentine --status-level fail` passed 7/7
  tests. `cargo clippy -p cfd-3d --no-default-features --lib -- -D warnings`
  passed. A targeted scan over `crates/cfd-3d/src/serpentine` and
  `crates/cfd-3d/src/scalar.rs` found no direct `num_traits`,
  `FromPrimitive`, `ToPrimitive`, `num_traits::Float`, `Float::`,
  `T::zero()`, `T::one()`, or `from_f64_or_one` residue.
- **Residual risk**: Direct cfd-3d `num-traits` ownership is now closed. Full
  Serpentine provider completion still requires nalgebra `Vector3` and FEM
  matrix/vector replacement with the Atlas geometry/storage stack.

---

## Sprint 2026-07-04: cfd-3d IBM Eunomia Scalar Seam
- **Resolved**: `crates/cfd-3d/src/ibm/{forcing,interpolation,solver}.rs` no
  longer import or bound direct `num_traits::{FromPrimitive,ToPrimitive}` and
  no longer call `num_traits::Float::*`, `T::zero()`, `T::one()`,
  `T::from_usize`, or `.to_isize()`. Direct and feedback forcing, the
  Roma/Peskin interpolation kernels, and IBM stencil index arithmetic now
  route scalar constants and elementary operations through `cfd-3d::scalar`,
  backed by Eunomia.
- **Boundary**: This is the IBM scalar-provider cleanup. It preserves the
  current nalgebra `Vector3` geometry/vector boundary because vector
  replacement belongs to the Leto/Gaia pass. The `Box<dyn ForcingMethod<T>>`
  dispatch remains unchanged; this slice removes direct scalar-provider
  ownership rather than redesigning IBM forcing dispatch.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-3d` completed. `cargo
  check -p cfd-3d --no-default-features` passed. Focused `cargo nextest run
  -p cfd-3d --no-default-features ibm --status-level fail` passed 12/12
  tests. `cargo clippy -p cfd-3d --no-default-features --lib -- -D warnings`
  passed. A targeted scan over `crates/cfd-3d/src/ibm` and
  `crates/cfd-3d/src/scalar.rs` found no direct `num_traits`,
  `FromPrimitive`, `ToPrimitive`, `num_traits::Float`, `T::zero()`,
  `T::one()`, `T::from_usize`, `.to_isize()`, or `Float::` residue.
- **Residual risk**: Direct cfd-3d `num-traits` ownership is now closed. Full
  IBM provider completion still requires nalgebra `Vector3` replacement with
  the Atlas geometry/vector stack.

---

## Sprint 2026-07-04: cfd-3d Spectral Eunomia Scalar Seam
- **Resolved**: `crates/cfd-3d/src/spectral/{chebyshev,poisson,solver}.rs`
  and `spectral/chebyshev_tests.rs` no longer import or call direct
  `num_traits::{FromPrimitive,Float}` or `num_traits::Float::*`. Chebyshev
  collocation constants, differentiation-matrix signs, quadrature constants,
  interpolation zero checks, Poisson boundary row constants, and spectral
  configuration tolerances now route through `cfd-3d::scalar`, backed by
  Eunomia `FloatElement`/`NumericElement`.
- **Boundary**: This is the spectral scalar-provider cleanup. It preserves the
  current nalgebra `DMatrix`/`DVector` implementation and `RealField` bounds
  because replacing spectral matrix storage/operations belongs to the Leto
  dense-matrix migration. The spectral module already routes FFT work through
  Apollo; no `rustfft` residue was present in this scope.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-3d` completed. `cargo
  check -p cfd-3d --no-default-features` passed. Focused `cargo nextest run
  -p cfd-3d --no-default-features spectral --status-level fail` passed 41/41
  tests. `cargo clippy -p cfd-3d --no-default-features --lib -- -D warnings`
  passed. A targeted scan over `crates/cfd-3d/src/spectral` and
  `crates/cfd-3d/src/scalar.rs` found no direct `num_traits`,
  `FromPrimitive`, `ToPrimitive`, `num_traits::Float`, `T::zero()`,
  `T::one()`, or `Float::` residue.
- **Residual risk**: Full spectral provider completion still requires Leto
  replacement for nalgebra `DMatrix`/`DVector` and the public `RealField`
  surface. Direct cfd-3d `num-traits` ownership is now closed.

---

## Sprint 2026-07-04: cfd-3d Level-Set Eunomia Scalar Seam
- **Resolved**: `crates/cfd-3d/src/level_set/{weno,advection,solver}.rs` no
  longer import or bound direct `num_traits::{FromPrimitive,Float}`. WENO5-Z
  weights, SSPRK3 coefficients, transport input validation, narrow-band limits,
  and reinitialization/Godunov math now route through `cfd-3d::scalar`, backed
  by Eunomia `FloatElement` and `NumericElement`.
- **Boundary**: This is the level-set module scalar-provider cleanup. Direct
  cfd-3d `num-traits` ownership is now closed by the later root
  lib-test/manifest slice. This entry still preserves the current
  `nalgebra::Vector3` boundary pending the larger Leto/Gaia geometry
  migration.
- **Evidence tier**: static source audit, compile-time integration, focused
  empirical nextest, and lib clippy. `cargo fmt -p cfd-3d --check` passed.
  `cargo check -p cfd-3d --no-default-features` passed. Focused `cargo
  nextest run -p cfd-3d --no-default-features level_set --status-level fail`
  passed 13/13 tests. `cargo clippy -p cfd-3d --no-default-features --lib --
  -D warnings` passed. A targeted scan over `crates/cfd-3d/src/level_set` and
  `crates/cfd-3d/src/scalar.rs` found no direct `num_traits`,
  `FromPrimitive`, `ToPrimitive`, `num_traits::Float`, `T::zero()`,
  `T::one()`, or `Float::` residue.
- **Residual risk**: `cargo clippy -p cfd-3d --no-default-features
  --all-targets -- -D warnings` is still blocked by pre-existing lint debt in
  unrelated cfd-3d test/module code (`poiseuille_test`, `fem_tests`,
  `smagorinsky_test`, `blueprint_integration`, `vof_tests`,
  `robustness_tests`, `bifurcation`, `venturi`, `trifurcation`, and VOF
  modules). Full provider completion still requires the remaining cfd-3d
  direct `num-traits` cleanup plus Leto/Gaia replacement of nalgebra geometry
  and storage surfaces.

---

## Sprint 2026-07-04: cfd-validation Direct num-traits Removal
- **Resolved**: `cfd-validation` no longer declares `num-traits` directly and
  no longer contains direct `num_traits::{FromPrimitive,ToPrimitive,Float}`,
  `T::from_*`, generic `.to_f64()`, or `num_traits::Float` residue in its
  source, tests, or package manifest. The final direct residues in
  `benchmarks::{cavity,step,poiseuille_bifurcation}` and
  `tests/complex_boundary_mms_validation.rs` now route scalar construction,
  diagnostic conversion, absolute values, powers, and max/min dispatch through
  Eunomia `FloatElement` or the crate-local `cfd-validation::scalar` adapter.
- **Boundary**: This is the direct scalar-provider cleanup for the validation
  crate. It intentionally preserves the current `nalgebra`/`nalgebra-sparse`
  matrix/vector/storage boundary and the remaining transitive `num-traits`
  graph through nalgebra, approx, Leto, Eunomia, Gaia, and other providers.
  Those are separate Leto/Gaia/upstream provider migration steps.
- **Evidence tier**: static source/dependency audit, compile-time integration,
  and empirical nextest. `cargo fmt -p cfd-validation --check` passed. `cargo
  check -p cfd-validation` passed. Full `cargo nextest run -p cfd-validation
  --status-level fail` passed 431/431 tests. A focused scan over
  `crates/cfd-validation/src`, `crates/cfd-validation/tests`, and
  `crates/cfd-validation/Cargo.toml` found no direct `num_traits`,
  `num-traits`, `FromPrimitive`, `ToPrimitive`, `T::from_f64`,
  `T::from_usize`, `T::from_u32`, `T::from_i32`, generic `.to_f64()`, or
  `num_traits::Float` residue. `cargo tree -p cfd-validation -e normal -i
  num-traits` still shows transitive paths only, primarily through nalgebra,
  approx, Leto, Eunomia, Gaia, and provider crates.
- **Residual risk**: Full validation-provider completion still requires
  replacing nalgebra/nalgebra-sparse storage/vector surfaces with Leto and mesh
  ownership with Gaia, plus resolving transitive provider dependencies in their
  owning upstream crates.

---

## Sprint 2026-07-04: cfd-2d Pressure-Velocity and Validation 2D Benchmark Eunomia Cleanup
- **Resolved**: The touched cfd-2d pressure-velocity/SIMPLEC solver surface no
  longer exposes direct `num_traits::{FromPrimitive,ToPrimitive}` bounds.
  `simplec_pimple::{config,solver,algorithms,simplec,pimple,interpolation,
  diagnostics}`, `pressure_velocity::{config,solver,pressure,correction,
  faces,rhie_chow}`, and the core `physics::momentum::{solver,solve}` impls
  now use Eunomia `FloatElement`/`NumericElement` or `cfd-2d::scalar` for
  scalar constants and conversions in the touched paths. The dependent
  `cfd-validation::benchmarks::{bifurcation,venturi,trifurcation,serpentine}`
  wrappers no longer name direct `num_traits` bounds and route touched
  constants, f64 diagnostics, square roots, powers, and zero values through
  `cfd-validation::scalar`.
- **Boundary**: This slice covers the cfd-2d pressure-velocity/SIMPLEC
  solver-facing cone and the four simple 2D validation benchmark wrappers. It
  does not claim full cfd-2d momentum cleanup: direct provider residue remains
  in lower momentum coefficient, boundary, interpolation, MUSCL, TVD, and WENO
  helper modules outside the touched core solver/solve impls. It also
  preserves the current nalgebra storage/vector boundary pending Leto
  migration.
- **Evidence tier**: static source audit, compile-time integration, and
  empirical nextest. `cargo fmt -p cfd-2d -p cfd-validation --check` passed.
  `cargo check -p cfd-validation` passed. Focused `cargo nextest run -p
  cfd-2d pressure_velocity simplec_pimple --status-level fail` passed 20/20
  tests. Full `cargo nextest run -p cfd-2d --status-level fail` passed
  567/567 tests with 27 skipped. Full `cargo nextest run -p cfd-validation
  --status-level fail` passed 431/431 tests. Targeted scans found no direct
  `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_f64`,
  `T::from_usize`, or generic `.to_f64()` residue in the touched cfd-2d cone,
  and no direct `num_traits`, `SafeFromF64`, `T::from_*`,
  `from_f64_or_*`, or `try_from_f64` residue in the four validation wrappers.
- **Residual risk**: Full cfd-2d provider completion still requires migrating
  the lower momentum helper modules and replacing nalgebra vector/matrix
  storage with Leto. Full validation benchmark completion still requires
  cavity, backward-facing step, disabled poiseuille-bifurcation, and broader
  Leto/Gaia benchmark storage and mesh ownership work.

---

## Sprint 2026-07-04: cfd-validation 3D Benchmark Eunomia Cleanup
- **Resolved**: `cfd-validation::benchmarks::BenchmarkConfig::default` now
  constructs tolerance and Reynolds-number defaults through
  `cfd-validation::scalar`, backed by Eunomia `FloatElement`.
  `benchmarks::threed::{bifurcation,venturi,serpentine}` no longer carry
  direct `num_traits::{Float,FromPrimitive,ToPrimitive}` or
  `cfd_core::conversion::SafeFromF64` bounds in the validation wrappers.
  Bifurcation and Venturi benchmark constants and validation predicates now
  route scalar construction and absolute values through the same crate-local
  adapter.
- **Boundary**: This slice covers
  `crates/cfd-validation/src/benchmarks/mod.rs`,
  `crates/cfd-validation/src/benchmarks/threed/bifurcation.rs`,
  `crates/cfd-validation/src/benchmarks/threed/serpentine.rs`, and
  `crates/cfd-validation/src/benchmarks/threed/venturi.rs`. It preserves the
  current benchmark geometry, mesh, solver, and nalgebra/Leto storage
  boundaries; this is a scalar-provider cleanup, not the full benchmark
  storage or mesh ownership migration.
- **Evidence tier**: static source audit, compile-time integration, and
  empirical nextest. `cargo fmt -p cfd-validation --check` passed. `cargo
  check -p cfd-validation` passed. Focused `cargo nextest run -p
  cfd-validation bifurcation_flow_3d serpentine_flow_3d venturi_flow_3d`
  passed 3/3 tests. Full `cargo nextest run -p cfd-validation --status-level
  fail` passed 431/431 tests. A targeted scan over the touched benchmark files
  found no direct `num_traits`, `FromPrimitive`, `ToPrimitive`, `SafeFromF64`,
  `from_f64_or_one`, `T::from_f64`, `T::from_i32`, `T::from_usize`, or
  `num_traits::Float::abs` residue.
- **Residual risk**: The wrappers still require `std::convert::From<f64>`
  because the current `cfd-3d::{bifurcation,serpentine,venturi}` solver
  constructors require `T: From<f64>`. Full removal of that conversion
  contract belongs in `cfd-3d`. Remaining `cfd-validation` benchmark work
  includes nalgebra/Leto storage seams, Gaia mesh ownership, and direct
  provider residue outside the touched 3D benchmark wrappers.

---

## Sprint 2026-07-04: cfd-validation Literature and Numerical Eunomia Cleanup
- **Resolved**: `cfd-validation::literature::blood_flow_1d` no longer imports
  direct `num_traits::{FromPrimitive,ToPrimitive}` or uses
  `num_traits::Float` for scalar construction, Poiseuille resistance powers,
  flow/pressure absolute values, mass-error checks, or f64 diagnostics. These
  operations now route through `cfd-validation::scalar`, including the new
  `scalar::powi` adapter over Eunomia `FloatElement::powi`.
  `cfd-validation::numerical::linear_solver` no longer imports direct
  `num_traits::{Float,FromPrimitive}` for generic solver validation bounds or
  threshold construction; tolerances and convergence rates now use the same
  crate-local Eunomia scalar seam.
- **Boundary**: This slice covers
  `crates/cfd-validation/src/literature/blood_flow_1d.rs`,
  `crates/cfd-validation/src/numerical/linear_solver.rs`, and
  `crates/cfd-validation/src/scalar.rs`. It preserves the current
  `nalgebra::RealField` and nalgebra vector/matrix storage boundaries because
  those are the later Leto replacement seam, not the scalar-provider seam.
- **Evidence tier**: static source audit, compile-time integration, and
  empirical nextest. `cargo fmt -p cfd-validation` passed. `cargo fmt -p
  cfd-validation --check` passed. `cargo check -p cfd-validation` passed.
  Focused `cargo nextest run -p cfd-validation blood_flow literature chapman
  patankar` passed 6/6 tests. Full `cargo nextest run -p cfd-validation`
  passed 431/431 tests. A combined targeted scan over
  `crates/cfd-validation/src/literature`,
  `crates/cfd-validation/src/numerical/linear_solver.rs`, and
  `crates/cfd-validation/src/scalar.rs` found no direct `num_traits`,
  `FromPrimitive`, `ToPrimitive`, `num_traits::Float`, `T::from_*`, or generic
  `.to_f64()` residue.
- **Residual risk**: `cfd-validation` still has Atlas-replaced provider
  residue outside this slice, primarily benchmark/test direct providers and
  nalgebra storage/vector seams. Those remain tracked for later Leto/Gaia
  migration; this slice does not claim full crate completion.

---

## Sprint 2026-07-04: cfd-validation Geometry Eunomia Cleanup
- **Resolved**: `cfd-validation::geometry::{serpentine_2d,venturi}` and
  `cfd-validation::geometry::threed::serpentine` no longer use direct
  `T::from_f64`, `T::from_usize`, `num_traits::Zero`, or RealField elementary
  dispatch for the touched geometry scalar constants, sine/cosine, square
  roots, powers, absolute values, and max selection. These operations now
  route through `cfd-validation::scalar`, backed by Eunomia
  `FloatElement`/`NumericElement`.
- **Boundary**: This slice covers
  `crates/cfd-validation/src/geometry/serpentine_2d.rs`,
  `crates/cfd-validation/src/geometry/venturi.rs`, and
  `crates/cfd-validation/src/geometry/threed/serpentine.rs`. It preserves the
  current `nalgebra::RealField` trait boundary and existing geometry trait
  shapes. Leto/Gaia storage and mesh ownership remain separate provider
  migration steps.
- **Evidence tier**: static source audit, compile-time integration, and
  empirical nextest. `cargo fmt -p cfd-validation` passed. `cargo check -p
  cfd-validation` passed. Focused `cargo nextest run -p cfd-validation
  geometry` passed 11/11 tests. Full `cargo nextest run -p cfd-validation`
  passed 431/431 tests. A geometry-wide scan found no `SafeFromF64`,
  `from_f64_or_*`, `try_from_f64`, `FromPrimitive`, `ToPrimitive`,
  `num_traits`, `T::from_f64`, `T::from_usize`, `T::from_u32`,
  `NumericElement::to_f64`, `num_complex`, `ComplexField`,
  `nalgebra::scalar`, `.to_f64().unwrap`, or generic `.to_f64()` residue under
  `crates/cfd-validation/src/geometry`.
- **Residual risk**: `cfd-validation` still has direct Atlas-replaced provider
  residue outside geometry scalar math, including
  `tests/complex_boundary_mms_validation.rs`,
  `src/numerical/linear_solver.rs`, `src/literature/blood_flow_1d.rs`, and
  benchmark modules. Geometry-specific residual risk is the later Leto/Gaia
  storage and mesh ownership migration.

---

## Sprint 2026-07-04: cfd-validation Richardson Eunomia Cleanup
- **Resolved**: `cfd-validation::manufactured::richardson::{core,
  validation,analysis}` no longer imports direct `num_traits` conversion or
  float traits, `nalgebra::ComplexField`, or test-only `nalgebra::scalar`
  helpers. Richardson order estimation, extrapolation, MMS validation grids,
  GCI calculations, report scoring, and confidence diagnostics now route
  scalar construction, powers, logarithms, square roots, absolute values,
  finite checks, and f64 diagnostics through `cfd-validation::scalar`, backed
  by Eunomia `FloatElement`/`NumericElement`.
- **Boundary**: This slice covers
  `crates/cfd-validation/src/manufactured/richardson/core.rs`,
  `crates/cfd-validation/src/manufactured/richardson/validation.rs`, and
  `crates/cfd-validation/src/manufactured/richardson/analysis.rs`. It
  preserves the current `nalgebra::RealField` generic boundary; replacing that
  with Leto remains a separate provider migration step.
- **Evidence tier**: static source audit, compile-time integration, and
  empirical nextest. `cargo fmt -p cfd-validation` passed. `cargo check -p
  cfd-validation` passed. Focused `cargo nextest run -p cfd-validation
  richardson` passed 16/16 tests. Full `cargo nextest run -p cfd-validation`
  passed 431/431 tests. A manufactured-wide scan found no `SafeFromF64`,
  `from_f64_or_*`, `try_from_f64`, `FromPrimitive`, `ToPrimitive`,
  `num_traits`, `T::from_f64`, `T::from_usize`, `T::from_u32`,
  `NumericElement::to_f64`, `num_complex`, `ComplexField`,
  `nalgebra::scalar`, `.to_f64().unwrap`, or generic `.to_f64()` residue under
  `crates/cfd-validation/src/manufactured`.
- **Residual risk**: `cfd-validation` still has direct Atlas-replaced provider
  residue outside manufactured MMS, including numerical linear-solver,
  benchmark modules, `literature/blood_flow_1d.rs`, literature/storage
  boundaries, and geometry storage/mesh seams.

---

## Sprint 2026-07-04: cfd-validation Manufactured Eunomia Cleanup
- **Resolved**: `cfd-validation::manufactured::{advection,
  advection_diffusion,burgers,diffusion,navier_stokes}` no longer imports
  direct `num_traits::FromPrimitive`, cfd-core `SafeFromF64`, or
  `nalgebra::ComplexField` for scalar constants and scalar transcendental
  functions. Manufactured advection, diffusion, Burgers, coupled
  advection-diffusion, polynomial Navier-Stokes, Taylor-Green, and Kovasznay
  scalar math now route through `cfd-validation::scalar`, backed by Eunomia
  `FloatElement`/`NumericElement`.
- **Boundary**: This slice covers
  `crates/cfd-validation/src/manufactured/advection.rs`,
  `crates/cfd-validation/src/manufactured/advection_diffusion.rs`,
  `crates/cfd-validation/src/manufactured/burgers.rs`,
  `crates/cfd-validation/src/manufactured/diffusion.rs`, and
  `crates/cfd-validation/src/manufactured/navier_stokes.rs`. It preserves the
  current `nalgebra::RealField` and `Vector2` storage boundary; replacing that
  with Leto remains a separate provider migration step.
- **Evidence tier**: static source audit, compile-time integration, and
  empirical nextest. `cargo fmt -p cfd-validation --check` passed. `cargo
  check -p cfd-validation` passed. Focused `cargo nextest run -p
  cfd-validation manufactured` passed 50/50 tests. Full `cargo nextest run -p
  cfd-validation` passed 431/431 tests. Focused scans over the five touched
  files found no `SafeFromF64`, `from_f64_or_*`, `try_from_f64`,
  `FromPrimitive`, `ToPrimitive`, `num_traits`, `T::from_f64`,
  `T::from_usize`, `T::from_u32`, `NumericElement::to_f64`, `num_complex`,
  `ComplexField`, `.to_f64().unwrap`, or generic `.to_f64()` residue.
- **Residual risk**: Richardson manufactured residue is closed by the later
  Richardson cleanup slice. Broader `cfd-validation` residue remains in
  numerical linear-solver, benchmark modules, `literature/blood_flow_1d.rs`,
  literature/storage boundaries, and geometry/storage seams.

---

## Sprint 2026-07-04: cfd-validation Chapman-Enskog Eunomia Cleanup
- **Resolved**: `cfd-validation::literature::chapman_enskog` no longer imports
  direct `num_traits::{FromPrimitive,ToPrimitive}` or calls
  `T::from_f64`/generic `.to_f64()` conversion. Chapman-Enskog viscosity,
  thermal-conductivity, expected-accuracy, and report values now route through
  `cfd-validation::scalar`, backed by Eunomia `FloatElement`.
- **Boundary**: This slice covers
  `crates/cfd-validation/src/literature/chapman_enskog.rs`. It preserves the
  current `nalgebra::RealField` bound as the storage/math boundary pending a
  later Leto replacement and does not claim full literature-module provider
  completion.
- **Evidence tier**: static source audit, compile-time integration, and
  empirical nextest. `cargo fmt -p cfd-validation --check` passed. Focused
  `cargo nextest run -p cfd-validation chapman` passed 2/2 tests. Full `cargo
  nextest run -p cfd-validation` passed 431/431 tests. Focused scans over
  `chapman_enskog.rs` found no `SafeFromF64`, `from_f64_or_*`,
  `try_from_f64`, `FromPrimitive`, `ToPrimitive`, `num_traits`,
  `T::from_f64`, `T::from_usize`, `T::from_u32`, `NumericElement::to_f64`,
  `num_complex`, `ComplexField`, `.to_f64().unwrap`, or generic `.to_f64()`
  residue.
- **Residual risk**: `cfd-validation` still has direct Atlas-replaced provider
  residue outside this slice, including numerical linear-solver, benchmark,
  `literature/blood_flow_1d.rs`, the `literature/patankar_1980.rs` storage
  boundary, geometry/storage, and remaining manufactured/Richardson modules.

---

## Sprint 2026-07-04: cfd-validation Time Integration Eunomia Cleanup
- **Resolved**: `cfd-validation::time_integration::{integrators,
  stability_analysis}` no longer imports direct `num_traits::{FromPrimitive,
  ToPrimitive}` or cfd-core `SafeFromF64`. Runge-Kutta constants, CFL
  scenario scalars, von-Neumann wave-number interpolation factors, and f64
  report diagnostics now route through `cfd-validation::scalar`, backed by
  Eunomia `FloatElement`/`NumericElement`. The local von-Neumann spatial
  operators use `eunomia::Complex` rather than `num_complex::Complex`.
- **Boundary**: This slice covers
  `crates/cfd-validation/src/time_integration/integrators.rs` and
  `crates/cfd-validation/src/time_integration/stability_analysis.rs`.
  `results.rs` and `validation.rs` were already on the crate-local scalar
  helper where they construct generic constants. The current
  `nalgebra::DVector`/`DMatrix` storage boundary is intentionally preserved;
  replacing that with Leto remains a separate provider migration step.
- **Evidence tier**: static source audit, compile-time integration, and
  empirical nextest. `cargo fmt -p cfd-validation` passed. `cargo check -p
  cfd-validation` passed. Focused `cargo nextest run -p cfd-validation
  time_integration` passed 12/12 tests. Full `cargo nextest run -p
  cfd-validation` passed 429/429 tests. Focused scans over
  `time_integration` found no `SafeFromF64`, `from_f64_or_*`,
  `try_from_f64`, `FromPrimitive`, `ToPrimitive`, `num_traits`,
  `T::from_f64`, `T::from_usize`, `T::from_u32`, `NumericElement::to_f64`,
  `num_complex`, `ComplexField`, or `.to_f64().unwrap` residue.
- **Residual risk**: `cfd-validation` still has direct Atlas-replaced
  provider residue outside this slice, including numerical linear-solver,
  benchmark, literature, geometry/storage, and remaining
  manufactured/Richardson modules.

---

## Sprint 2026-07-04: cfd-validation Error Metrics Eunomia Cleanup
- **Resolved**: `cfd-validation::error_metrics::{norms,normalized,statistics}`
  no longer imports direct `num_traits::FromPrimitive` or calls
  `T::from_f64`/`T::from_usize`. Error metric denominator construction,
  relative-error tolerance construction, absolute differences, and variance
  square roots now route through `cfd-validation::scalar`, backed by Eunomia
  `FloatElement`/`NumericElement`.
- **Boundary**: This slice covers
  `crates/cfd-validation/src/error_metrics/norms.rs`,
  `crates/cfd-validation/src/error_metrics/normalized.rs`, and
  `crates/cfd-validation/src/error_metrics/statistics.rs`. It intentionally
  preserves the current `nalgebra::RealField`/`Vector3` storage and vector
  norm boundary; Leto replacement for that storage remains a separate
  migration step.
- **Evidence tier**: static source audit, compile-time integration, and
  empirical nextest. `cargo fmt -p cfd-validation` passed. `cargo check -p
  cfd-validation` passed. Focused `cargo nextest run -p cfd-validation
  error_metrics` passed 21/21 tests. Full `cargo nextest run -p
  cfd-validation` passed 429/429 tests. Focused scans over `error_metrics`
  found no `SafeFromF64`, `from_f64_or_*`, `try_from_f64`, `FromPrimitive`,
  `ToPrimitive`, `num_traits`, `T::from_f64`, `T::from_usize`,
  `T::from_u32`, `ComplexField`, or `.to_f64().unwrap` residue.
- **Residual risk**: `cfd-validation` still has direct Atlas-replaced
  provider residue outside this slice, including numerical linear-solver,
  benchmark, literature, geometry/storage, and remaining
  manufactured/Richardson modules.

---

## Sprint 2026-07-04: cfd-validation Eunomia Compile-Blocker Cleanup
- **Resolved**: The active `cfd-validation` compile blocker group is closed.
  `literature::patankar_1980` now uses Eunomia scalar construction for
  reference constants and integer-rational grid reference lookup instead of
  direct `num_traits`/`SafeFromF64` conversions. `manufactured::{
  advanced_physics,multi_physics,reynolds_stress_mms}` now use the
  crate-local Eunomia scalar adapter for constants and transcendental math in
  the touched paths, and the adapter exposes `cosh`/`tanh` for migrated
  hyperbolic calls. Richardson analysis now has the `FloatElement` bound
  required by `RichardsonMmsResult<T>` methods.
- **Boundary**: This slice covers
  `crates/cfd-validation/src/scalar.rs`,
  `crates/cfd-validation/src/literature/patankar_1980.rs`,
  `crates/cfd-validation/src/manufactured/advanced_physics.rs`,
  `crates/cfd-validation/src/manufactured/multi_physics.rs`,
  `crates/cfd-validation/src/manufactured/reynolds_stress_mms.rs`, and the
  Richardson analysis bound in
  `crates/cfd-validation/src/manufactured/richardson/analysis.rs`. It follows
  the earlier conservation/convergence/geometry manufactured cleanup and does
  not claim full `cfd-validation` provider completion.
- **Evidence tier**: static source audit, compile-time integration, and
  empirical nextest. `cargo fmt -p cfd-validation` passed. `cargo check -p
  cfd-validation` passed. `cargo nextest run -p cfd-validation` passed 429/429
  tests. Focused scans over the fully migrated files found no `SafeFromF64`,
  `from_f64_or_*`, `try_from_f64`, `T::from_f64`, `T::from_usize`,
  `ComplexField`, or `<T as FromPrimitive>` residue.
- **Residual risk**: `cfd-validation` still has direct Atlas-replaced provider
  residue outside this slice, including `numerical/linear_solver.rs`,
  benchmark modules, literature validations, and remaining
  manufactured/Richardson modules. The local `scalar::cbrt` helper still maps
  through `powf(1/3)` because Eunomia lacks a dedicated `cbrt` method; promote
  that upstream if additional consumers require a named cube root provider
  surface.

---

## Sprint 2026-07-04: cfd-io Checkpoint Validator Eunomia Constructor Cleanup
- **Resolved**: `cfd-io::checkpoint::CheckpointValidator` no longer uses
  direct `T::from_f64` scalar construction in the mass-conservation validator.
  Grid spacing and the central-difference denominator now route through
  Eunomia `FloatElement`, with the denominator scalar constructed once before
  the interior stencil loop.
- **Boundary**: This slice covers
  `crates/cfd-io/src/checkpoint/validator.rs`. It preserves the existing
  exact-dimension check that rejects zero or non-exactly-convertible grid
  dimensions before constructing generic spacings. It does not address
  provider residue in other CFDrs crates.
- **Evidence tier**: static source audit, compile-time integration, empirical
  nextest, lint, doctest, and rustdoc. `cargo fmt -p cfd-io --check`,
  `cargo check -p cfd-io`, `cargo clippy -p cfd-io --all-targets -- -D
  warnings`, `cargo nextest run -p cfd-io` (3/3), `cargo test --doc -p
  cfd-io`, and `cargo doc -p cfd-io --no-deps` passed. A focused scan over
  `crates/cfd-io/Cargo.toml`, `crates/cfd-io/src`, and `crates/cfd-io/tests`
  found no direct `num_traits`, `num-traits`, `FromPrimitive`, `ToPrimitive`,
  `T::from_f64`, `T::from_usize`, `T::from_u32`, `.to_f64().unwrap`,
  `nalgebra::try_convert`, `ndarray`, or `nalgebra` residue.
- **Residual risk**: The broader workspace still has direct Atlas-replaced
  provider residue in `cfd-2d`, `cfd-3d`, `cfd-core`, `cfd-math`, and
  `cfd-validation`; those remain separate crate-scoped migration slices.

---

## Sprint 2026-07-04: cfd-1d Direct num-traits Removal and Resistance/Vascular Eunomia Cleanup
- **Resolved**: `cfd-1d` no longer declares or directly references
  `num-traits`. The resistance scalar contract now uses Eunomia
  `FloatElement`/`NumericElement` for scalar construction, diagnostics, and
  transcendental/math operations. The touched hydraulic resistance,
  serpentine, slug-flow, Bessel/Womersley, structured-tree, bifurcation,
  network blueprint/sink, solver-analysis, and package-test seams no longer
  import `num_traits`, use `FromPrimitive`/`ToPrimitive`, call
  `T::from_f64`/`T::from_usize`/`T::from_u32`, or bridge through
  `nalgebra::try_convert`.
- **Boundary**: This closes direct `num-traits` ownership for the `cfd-1d`
  crate. It does not remove transitive `num-traits` through `approx`,
  `nalgebra`, `nalgebra-sparse`, `half`/Eunomia/Leto/Hephaestus,
  `num-complex`, or other provider stacks. It also does not replace the
  remaining nalgebra/nalgebra-sparse storage and solve boundaries; those stay
  as Leto-backed dense/sparse migration work.
- **Evidence tier**: static source audit, compile-time integration, empirical
  nextest, and dependency-tree audit. `cargo fmt -p cfd-1d --check` passed.
  `cargo check -p cfd-1d` passed. `cargo nextest run -p cfd-1d` passed
  725/725 tests with 3 skipped. A full scan over `crates/cfd-1d/Cargo.toml`,
  `crates/cfd-1d/src`, and `crates/cfd-1d/tests` found no direct
  `num_traits`, `num-traits`, `FromPrimitive`, `ToPrimitive`, `T::from_f64`,
  `T::from_usize`, `T::from_u32`, `nalgebra::try_convert`, or
  `.to_f64().unwrap` residue.
- **Residual risk**: `cargo clippy -p cfd-1d --all-targets -- -D warnings`
  remains blocked by existing all-target lint debt outside this provider
  cleanup, including `blueprint_solve_trace.rs`, `adversarial_tests.rs`,
  `resistance_model_validation.rs`, `medical_millifluidic_screening.rs`,
  `geometry_integration_demo.rs`, cell-separation tests/modules,
  droplet-regime tests, entrance-model tests, matrix-assembly tests, and
  venturi coefficient tests. Broader Atlas migration work remains for
  Leto/nalgebra-sparse storage replacement and Hephaestus higher-level GPU
  kernel ownership.

---

## Sprint 2026-07-04: cfd-1d Solver-Core Eunomia Boundary Cleanup
- **Resolved**: `cfd-1d` solver core no longer owns direct
  `num_traits::{Float, FromPrimitive, ToPrimitive}` at the shared primary
  network-solver scalar contract. `NetworkSolveScalar` is now expressed through
  nalgebra's current linear-system boundary plus Eunomia `RealField` and the
  cfd-core conversion traits. Anderson acceleration, linear solve
  equilibration, Jacobi preconditioning, convergence checks, SPD detection,
  residual norms, finite-vector validation, and f64 diagnostics route through
  Eunomia `FloatElement`/`NumericElement`.
- **Boundary**: This slice covers
  `crates/cfd-1d/src/solver/core/{mod,anderson_acceleration,linear_system,convergence,solver_detection,geometry}.rs`.
  It does not remove the active nalgebra/nalgebra-sparse matrix/vector storage
  boundary; that requires a later Leto-backed sparse/dense linear-system
  replacement.
- **Evidence tier**: compile-time integration, empirical nextest, and static
  source audit. `cargo fmt -p cfd-1d --check` passed. `cargo check -p cfd-1d`
  passed. `cargo nextest run -p cfd-1d` passed 725/725 tests with 3 skipped.
  A focused scan over the migrated solver-core files found no direct
  `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_f64`,
  `nalgebra::try_convert`, or `<T as Float>` residue.
- **Residual risk**: cfd-1d still has direct provider residue in vascular
  Bessel/Womersley, resistance scalar traits, tests/benches, and the
  nalgebra/nalgebra-sparse storage boundary. The broader CFDrs migration still
  requires replacing those storage and solver backends with Leto/Hephaestus
  where provider functionality is available.

---

## Sprint 2026-07-04: cfd-1d Network Wrapper Eunomia Boundary Cleanup
- **Resolved**: `domain/network/wrapper.rs` no longer owns direct
  `num_traits` scalar construction or nalgebra scalar bridge conversions.
  `EdgeProperties::from`, network characteristic length, Picard resistance
  updates, blood hematocrit propagation, Pries phase-separation bridge values,
  coefficient validation, and parallel conductance estimation now route scalar
  construction, f64 extraction, finite checks, and absolute values through
  `SafeFromF64`/`SafeFromUsize` plus Eunomia `NumericElement`. Adjacent
  `MatrixAssembler` no longer carries a `FromPrimitive` bound, and
  `NetworkProblem`/`NetworkSolver` use the existing `NetworkSolveScalar`
  contract where wrapper methods now require the provider scalar surface.
- **Boundary**: This slice covers
  `crates/cfd-1d/src/domain/network/wrapper.rs`,
  `crates/cfd-1d/src/solver/core/matrix_assembly.rs`,
  `crates/cfd-1d/src/solver/core/problem.rs`, and the solver finite-check
  adjustments in `crates/cfd-1d/src/solver/core/mod.rs`. It does not remove
  the legacy `NetworkSolveScalar` compatibility bounds in `solver/core/mod.rs`.
- **Evidence tier**: compile-time integration, empirical nextest, and static
  source audit. `cargo fmt -p cfd-1d --check` passed. `cargo check -p cfd-1d`
  passed. `cargo nextest run -p cfd-1d` passed 725/725 tests with 3 skipped.
  A focused scan over `wrapper.rs`, `matrix_assembly.rs`, and `problem.rs`
  found no direct `num_traits`, `FromPrimitive`, `ToPrimitive`,
  `T::from_f64`, `T::from_usize`, or `nalgebra::try_convert` residue.
- **Residual risk**: This wrapper cleanup is superseded by the later
  solver-core cleanup for `NetworkSolveScalar`. Additional provider work
  remains in vascular Bessel/Womersley, resistance scalar traits,
  tests/benches, and Leto storage replacement.

---

## Sprint 2026-07-04: cfd-1d Network Blueprint/Sink Eunomia Boundary Cleanup
- **Resolved**: The canonical `network_from_blueprint` construction path no
  longer uses direct `T::from_f64`/`T::from_usize` scalar construction for
  blueprint-derived resistances, areas, physical constants, boundary values,
  serpentine segment counts, flow probes, and blood defaults. These values now
  route through `SafeFromF64`/`SafeFromUsize`, and the flow-invariant policy
  check uses Eunomia `NumericElement::abs` for generic absolute value.
  `NetworkBuilderSink` carries the Atlas conversion bounds required by the
  canonical blueprint entry point.
- **Boundary**: This slice covers
  `crates/cfd-1d/src/domain/network/builder/blueprint_conversion.rs` and
  `crates/cfd-1d/src/domain/network/sink.rs`. It intentionally preserves the
  current direct `FromPrimitive` bound inherited from `Network::new` in
  `domain/network/wrapper.rs`; removing that requires migrating the broader
  wrapper impl and its remaining direct scalar conversions.
- **Evidence tier**: compile-time integration and empirical nextest.
  `cargo fmt -p cfd-1d` passed. `cargo check -p cfd-1d` passed. `cargo
  nextest run -p cfd-1d` passed 725/725 tests with 3 skipped.
- **Residual risk**: `domain/network/wrapper.rs` still contains direct
  `FromPrimitive`, `T::from_f64`, `T::from_usize`, `nalgebra::try_convert`,
  and generic `.abs()` residue. Broader `cfd-1d` provider work also remains in
  solver-core, vascular Bessel/Womersley, tests/benches, and the current
  nalgebra/nalgebra-sparse storage boundaries.

---

## Sprint 2026-07-04: cfd-1d Domain Components Eunomia Boundary Cleanup
- **Resolved**: `cfd-1d` domain components no longer carry direct
  `num_traits` conversion/math bounds. The component trait pressure-drop
  calculation uses Eunomia `NumericElement::abs`; factory defaults and
  component constants use the existing Atlas provider conversion seam; and
  channel, membrane, mixer, pump, valve, and sensor implementations no longer
  import `FromPrimitive` or `Float`.
- **Boundary**: This slice covers
  `crates/cfd-1d/src/domain/components/{mod,channels,factory,membranes,mixers,pumps,sensors,valves}.rs`.
  It preserves the current nalgebra `RealField` and resistance-model
  boundaries for later Leto/Eunomia work.
- **Evidence tier**: compile-time integration, empirical nextest, static
  source audit, and lint regression audit for the touched file. `cargo fmt -p
  cfd-1d --check` passed. `cargo check -p cfd-1d` passed. `cargo nextest run
  -p cfd-1d` passed 725/725 tests with 3 skipped. A focused component scan
  found no direct `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_f64`,
  `T::from_usize`, `Float::`, or generic `.abs()` residue. `cargo clippy -p
  cfd-1d --all-targets -- -D warnings` no longer reports the touched
  `channels.rs` item-order lint.
- **Residual risk**: Full `cfd-1d` all-target clippy remains blocked by
  existing lint debt outside this slice in tests/examples, domain-network,
  cell-separation, vascular, solver-core, and benches. Broader `cfd-1d`
  direct-provider residue remains in solver-core, domain-network,
  vascular Bessel/Womersley, tests, and remaining nalgebra/nalgebra-sparse
  storage boundaries.

---

## Sprint 2026-07-04: cfd-1d Channel/Branching/Analysis Eunomia Boundary Cleanup
- **Resolved**: The next coherent `cfd-1d` provider seam no longer depends on
  direct `num_traits` bounds. Channel flow-regime classification, channel
  flow-resistance constants and powers, channel geometry perimeter math,
  Poiseuille shape factors, branching network solver bounds, and network
  pressure/flow/resistance/performance analysis paths now route scalar
  construction, powers, square roots, absolute values, and scalar-to-f64
  display/oracle conversion through `SafeFromF64`, Eunomia `FloatElement`, and
  Eunomia `NumericElement`.
- **Boundary**: This slice covers `crates/cfd-1d/src/domain/channel`, the
  `domain/junctions/branching` solver/physics/validation cone, and
  `solver/analysis` aggregates/analyzers. It does not remove the direct
  `cfd-1d` manifest `num-traits` dependency because other domains still use
  `FromPrimitive`, `ToPrimitive`, and `Float`.
- **Evidence tier**: compile-time integration, empirical nextest, and static
  source audit. `cargo check -p cfd-1d` passed. `cargo nextest run -p cfd-1d`
  passed 725/725 tests with 3 skipped. Focused source scans found no direct
  `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_f64`, generic
  `.to_f64()`, or `Float::` residue in the touched channel/branching/analyzer
  cone or in `solver/analysis`.
- **Residual risk**: Full `cfd-1d` all-target clippy remains blocked by
  unrelated existing lint debt in examples, tests, and cell-separation/
  resistance modules. The separate domain-components direct `num-traits`
  residue is closed by the later domain-components slice; broader direct
  provider residue remains in solver-core, domain-network,
  vascular Bessel/Womersley, and other package areas outside this slice.

---

## Sprint 2026-07-04: cfd-1d Murray's-Law Eunomia Boundary Cleanup
- **Resolved**: `cfd-1d` Murray's-law vascular geometry no longer depends on
  `num_traits::FromPrimitive`. `MurraysLaw` and `OptimalBifurcation` route
  scalar constants, power functions, absolute value, and inverse cosine
  through Eunomia `FloatElement`/`NumericElement`. This uses the current
  Eunomia `FloatElement::acos` surface, which already has value-semantic
  tests in the provider checkout.
- **Boundary**: This slice covers only
  `crates/cfd-1d/src/physics/vascular/murrays_law/{law,bifurcation}.rs` plus
  the required bound propagation in `physics/vascular/bifurcation.rs`.
  `cfd-1d` still declares `num-traits` because many other modules still use
  `FromPrimitive`, `ToPrimitive`, and `Float`.
- **Evidence tier**: static source audit plus formatting. `cargo fmt -p
  cfd-1d --check` passed, and a focused scan over the Murray's-law files found
  no `num_traits`, `FromPrimitive`, `ToPrimitive`, or unqualified generic
  `acos`/`powf`/`powi`/`abs` residue. The provider inverse-cosine contract was
  re-verified in the Eunomia checkout with `cargo nextest run -p eunomia acos`
  (2/2).
- **Residual risk**: `cargo check -p cfd-1d` remains blocked by unrelated
  dirty-tree errors outside this slice, including missing `SafeFromF64` bounds,
  generic `abs`/`powf` ambiguity, and stale `to_f64().unwrap_or(...)` call
  sites. Full `cfd-1d` nextest evidence is pending behind that cleanup.

---

## Sprint 2026-07-04: cfd-math Direct Num-Traits Boundary Cleanup
- **Resolved**: `cfd-math` no longer declares a direct `num-traits`
  dependency. The remaining package-owned `num_traits` source/test residue was
  replaced with Eunomia `NumericElement` and `FloatElement` in the GPU
  Laplacian operator and AMG integration-test helper.
- **Boundary**: This slice removes direct `cfd-math` ownership of
  `num-traits` only. It preserves current nalgebra and nalgebra-sparse solver
  boundaries, and it does not remove transitive `num-traits` entries through
  `approx`, Leto, Eunomia, Hermes, nalgebra, nalgebra-sparse, simba, or
  num-complex.
- **Evidence tier**: compile-time integration, empirical nextest, and static
  source/dependency audit. `cargo fmt -p cfd-math --check`, `cargo check -p
  cfd-math`, `cargo check -p cfd-math --features gpu`, `cargo clippy -p
  cfd-math --all-targets -- -D warnings`, `cargo clippy -p cfd-math
  --features gpu --all-targets -- -D warnings`, and `cargo nextest run -p
  cfd-math` passed (333/333). A focused source/manifest scan over `cfd-math`
  found no direct `num_traits`, `num-traits`, `FromPrimitive`, or
  `ToPrimitive` residue.
- **Residual risk**: `cargo tree -p cfd-math -e normal -i num-traits` still
  resolves transitive ownership through provider and nalgebra stacks. Full
  graph cleanup remains sequenced behind the remaining Leto/Eunomia/nalgebra
  replacement slices.

---

## Sprint 2026-07-04: cfd-core Direct Num-Traits Boundary Cleanup
- **Resolved**: `cfd-core` no longer declares a direct `num-traits`
  dependency for the touched scalar conversion boundary. `management::
  conversion` and mesh centroid scalar construction now route through Eunomia
  `FloatElement` instead of `num_traits::FromPrimitive`.
- **Boundary**: This slice removes direct `cfd-core` ownership of the
  conversion trait dependency in the touched cone only. It does not remove
  transitive `num-traits` entries pulled through providers, WGPU/Naga, Leto,
  Eunomia, nalgebra, nalgebra-sparse, simba, or num-complex.
- **Evidence tier**: compile-time integration, empirical nextest, and static
  source/dependency audit. `cargo fmt -p cfd-core --check`, `cargo check -p
  cfd-core`, `cargo clippy -p cfd-core --all-targets -- -D warnings`, and
  `cargo nextest run -p cfd-core` passed (229/229). A focused source/manifest
  scan over `crates/cfd-core` found no direct `num_traits`, `num-traits`,
  `FromPrimitive`, or `ToPrimitive` residue.
- **Residual risk**: `cargo tree -p cfd-core -e normal -i num-traits` still
  resolves transitive ownership through provider and nalgebra stacks. Full
  graph cleanup remains sequenced behind the remaining Leto/Eunomia/nalgebra
  replacement slices.

---

## Sprint 2026-07-04: Root Workspace Provider Declaration Cleanup
- **Resolved**: The root workspace no longer declares stale direct
  `wgpu 0.19` or `num-complex` provider dependencies. WGPU ownership for the
  migrated GPU boundary is now expressed through `hephaestus-wgpu`, and direct
  complex vocabulary at migrated CFDrs call sites is owned by Eunomia/Apollo
  provider APIs instead of a workspace-level num-complex entry.
- **Boundary**: This is a manifest/settings cleanup only. It does not remove
  active direct `nalgebra` or `num-traits` package dependencies whose call
  sites still require real crate-level migration.
- **Evidence tier**: static manifest/dependency audit plus compile-time
  metadata/package verification. `cargo metadata --no-deps --format-version 1`
  and `cargo check -p cfd-suite` passed. A focused manifest scan found no
  direct `wgpu`, `dep:wgpu`, `wgpu.workspace`, `num-complex`, or `num_complex`
  declarations under `Cargo.toml` or `crates/*/Cargo.toml`. `cargo tree
  --workspace -e normal -i wgpu@0.19.4` returned no matching package.
- **Residual risk**: `num-complex` still resolves transitively through
  provider/nalgebra stacks; full graph cleanup remains sequenced behind the
  remaining Leto/Eunomia and nalgebra replacement slices.

---

## Sprint 2026-07-04: cfd-core Hephaestus WGPU ABI Migration
- **Resolved**: `cfd-core` no longer owns a direct `wgpu 0.19` dependency.
  Its `gpu` feature depends only on `hephaestus-wgpu`, `GpuContext` acquires a
  `hephaestus_wgpu::WgpuDevice`, and raw transitional WGPU descriptor/buffer/
  pipeline symbols route through the provider-owned `hephaestus_wgpu::wgpu`
  re-export.
- **Boundary**: This is an ABI/provider handoff for the current raw WGPU GPU
  module. The compute algorithms are still CFDrs-owned kernels; deeper
  Hephaestus kernel abstraction replacement remains separate work.
- **Evidence tier**: compile-time integration, empirical nextest, and static
  dependency audit. `cargo fmt -p cfd-core --check`, `cargo check -p
  cfd-core`, `cargo clippy -p cfd-core --all-targets -- -D warnings`, and
  `cargo nextest run -p cfd-core` passed (228/228). `cargo tree -p cfd-core -e
  normal -i wgpu@0.19.4` returned no matching package; `cargo tree -p cfd-core
  -e normal -i wgpu@26.0.1` shows `wgpu -> hephaestus-wgpu -> cfd-core`.
  Focused scans find no `dep:wgpu`, `use wgpu::`,
  `wgpu::Maintain`, `pollster::`, or `global_id(` residue in `crates/cfd-core`.
- **Residual risk**: Full GPU replacement still needs moving CFDrs-owned
  kernels onto higher-level Hephaestus kernel APIs where the provider surface
  exists.

---

## Sprint 2026-07-04: cfd-optim Direct Nalgebra Dev-Dependency Cleanup
- **Resolved**: `cfd-optim` no longer declares a direct `nalgebra`
  dev-dependency. Focused source and manifest scans over the package find no
  direct `nalgebra`, `DMatrix`, `DVector`, `Vector2`, `Vector3`, `RealField`,
  `ndarray`, `num_traits`, `num-complex`, `num_complex`, `rayon`, or `tokio`
  residue. Package clippy cleanup replaced a cloned singleton slice with
  `std::slice::from_ref`, moved test modules after items, and made the
  computed SDT acoustic energy density and acoustic contrast factors part of
  `SdtMetrics`.
- **Boundary**: `cargo tree -p cfd-optim -e normal -i nalgebra` still resolves
  nalgebra transitively through upstream CFDrs/provider crates including
  `cfd-1d`, `cfd-2d`, `cfd-3d`, `cfd-core`, `cfd-math`, `cfd-validation`,
  Gaia, and `nalgebra-sparse`. This slice removes only direct `cfd-optim`
  ownership of that dependency.
- **Evidence tier**: compile-time integration plus empirical nextest and static
  source/dependency audit. `cargo fmt -p cfd-optim --check`, `cargo check -p
  cfd-optim`, `cargo clippy -p cfd-optim --all-targets -- -D warnings`, and
  `cargo nextest run -p cfd-optim` passed with 121/121 tests. Focused scans
  found no direct cfd-optim Atlas residue listed above.
- **Residual risk**: The broader CFDrs Atlas migration still needs upstream
  nalgebra removal from the crates listed in the boundary section before
  `cfd-optim` can be graph-clean transitively.

---

## Sprint 2026-07-04: cfd-io Provider Reference Sync
- **Resolved**: `crates/cfd-io/agents.md` no longer describes the stale
  pre-provider architecture where `cfd-io` directly depended on `cfd-core` and
  `cfd-math`, exposed older VTK module topology, or used obsolete checkpoint
  APIs. The crate-local reference now records the current provider boundary:
  Leto for dense checkpoint/binary arrays, Eunomia for scalar bounds, Consus
  for optional HDF5, RITK for optional VTK, and a local file-format
  `Error`/`Result` type.
- **Boundary**: This is documentation synchronization for an already migrated
  `cfd-io` code slice. Optional `vtk` still pulls its own upstream stack under
  `--all-features`; the default normal `cfd-io` graph remains free of
  nalgebra.
- **Evidence tier**: static documentation/code consistency audit plus
  compile-time and empirical package verification. `cargo fmt -p cfd-io
  --check`, `cargo check -p cfd-io`, `cargo check -p cfd-io --all-features`,
  `cargo clippy -p cfd-io --all-targets -- -D warnings`, and `cargo nextest
  run -p cfd-io` passed. `cargo tree -p cfd-io -e normal -i nalgebra`
  returned no matching package. Focused scans found no cfd-io source or
  manifest `nalgebra`, `DMatrix`, `DVector`, direct `num_traits`,
  `num-traits`, `FromPrimitive`, `ToPrimitive`, `cfd_core::error`, `anyhow`,
  or `bytemuck` residue.
- **Residual risk**: Other CFDrs crates still retain direct nalgebra,
  num-traits, GPU, and concurrency-provider residue. Those remain active under
  the broader Atlas migration goal.

---

## Sprint 2026-07-04: cfd-3d FEM Leto Solution Storage Boundary
- **Resolved**: `cfd-3d::fem::StokesFlowSolution` no longer stores velocity
  and pressure as nalgebra `DVector`. Both fields now use `FemDofVector`, a
  Leto `Array1`-backed FEM degree-of-freedom vector with explicit dense-slice,
  nalgebra-boundary conversion, norm, and checked difference-norm operations.
  FEM solver warm starts, projection previous-state setup, Anderson
  acceleration, and Venturi Picard convergence now use explicit conversion
  only at the current sparse linear-solver boundary.
- **Boundary**: This does not migrate the sparse linear solvers or element
  matrix internals. `cfd-3d::fem::{solver,projection_solver}` still allocate
  nalgebra `DVector` for RHS/linear-system solve boundaries, and
  `cfd-3d::fem::element`/`shape_functions` still retain nalgebra dense matrix
  internals pending the larger Leto matrix-storage slice.
- **Evidence tier**: compile-time integration plus empirical focused FEM tests
  and static residue audit. `cargo fmt -p cfd-3d --check` passed. `cargo
  check -p cfd-3d` passed. `cargo clippy -p cfd-3d --lib -- -D warnings`
  passed. `cargo nextest run -p cfd-3d fem` passed 40/40 tests. Focused scans
  found no `pub velocity: DVector`, `pub pressure: DVector`, stale
  `StokesFlowSolution`/`DVector` storage comments, direct velocity `rows`
  calls, or direct `&velocity - &velocity` convergence arithmetic.
- **Residual risk**: Full cfd-3d FEM provider completion still requires Leto
  matrix/vector storage for element matrices, shape-function gradients, sparse
  RHS/solution vectors, and the cfd-math linear-solver API.

---

## Sprint 2026-07-03: cfd-3d FEM Solver/Projection Eunomia Scalar Boundary
- **Resolved**: `cfd-3d::fem::{solver,projection_solver,solution}` no longer
  use direct `T::zero()`/`T::one()` identity construction in the targeted FEM
  solver, projection, and solution-blending paths. These identities now route
  through the FEM-local `fem::scalar` SSOT backed by Eunomia
  `NumericElement`/`FloatElement`.
- **Boundary**: This deliberately does not migrate nalgebra-backed
  `DVector`/matrix storage in FEM solution, solver, projection, or sparse
  assembly paths. That remains the next Leto-backed FEM storage slice.
- **Evidence tier**: compile-time integration plus empirical focused FEM tests
  and static residue audit. `cargo fmt -p cfd-3d --check` passed. `cargo
  check -p cfd-3d` passed. `cargo clippy -p cfd-3d --lib -- -D warnings`
  passed. `cargo nextest run -p cfd-3d fem` passed 40/40 tests. Focused scans
  found no direct `T::zero()`, `T::one()`, `num_traits`, `FromPrimitive`,
  `ToPrimitive`, `Float::`, or `Float` trait residue in `solver.rs`,
  `projection_solver.rs`, or `solution.rs`.
- **Residual risk**: Full cfd-3d FEM provider completion still requires Leto
  storage for FEM matrices/vectors and solver/projection call-site updates.

---

## Sprint 2026-07-03: cfd-3d FEM Element Eunomia Scalar Boundary
- **Resolved**: `cfd-3d::fem::element` no longer uses direct
  `num_traits::{Float, FromPrimitive}`, `num_traits::Float::abs`,
  exact-half/constant `expect` construction, or old `T::zero` identities in
  its scalar construction and elementary math paths. `ElementMatrices` and
  `FluidElement` now route zero initialization, tetrahedral volume absolute
  values, degeneracy tolerance, half factors, and sixth-volume constants
  through Eunomia `FloatElement`/`NumericElement` plus the FEM-local
  `fem::scalar` SSOT.
- **Boundary**: This deliberately does not migrate `DMatrix`/`DVector`
  storage in the element matrix and derivative paths. That belongs to the
  larger Leto-backed FEM storage slice with solver/projection call-site
  updates.
- **Evidence tier**: compile-time integration plus empirical focused FEM
  element tests and static residue audit. `cargo fmt -p cfd-3d --check`
  passed. `cargo check -p cfd-3d` passed. `cargo clippy -p cfd-3d --lib --
  -D warnings` passed. `cargo nextest run -p cfd-3d element` passed 5/5
  tests. `cargo check -p moirai-executor --tests` passed after an earlier
  transient build against a shifting dirty Moirai tree exposed a stale
  `worker_numa_nodes` initializer error. Focused scans found no direct
  `num_traits`, `FromPrimitive`, `ToPrimitive`, `Float::`, `Float` trait
  bound, old scalar construction, `T::zero()`, `T::one()`, or optional
  `RealField` extrema residue in `element.rs`.
- **Residual risk**: Full cfd-3d FEM provider completion still requires Leto
  storage for FEM matrices/vectors and follow-on Eunomia cleanup in `solver`,
  `projection_solver`, and `solution`.

---

## Sprint 2026-07-03: cfd-3d FEM Problem-Validation Eunomia Scalar Boundary
- **Resolved**: `cfd-3d::fem::problem_validation` no longer uses direct
  `num_traits::Float` dispatch or old `T::zero` identity checks for physical
  invariant validation. Finite checks and positive-value guards now route
  through Eunomia `FloatElement`/`NumericElement` plus the FEM-local
  `fem::scalar` SSOT. `StokesFlowProblem::validate` now carries the explicit
  Eunomia scalar contract required by this validation boundary.
- **Boundary**: This does not migrate FEM element assembly, solver/projection
  assembly, solution scalar identities, or nalgebra matrix/vector storage.
- **Evidence tier**: compile-time integration plus empirical focused FEM
  validation tests and static residue audit. `cargo fmt -p cfd-3d --check`
  passed. `cargo check -p cfd-3d` passed. `cargo clippy -p cfd-3d --lib --
  -D warnings` passed. `cargo nextest run -p cfd-3d validate` passed 10/10
  tests. Focused scans found no direct `num_traits`, `FromPrimitive`,
  `ToPrimitive`, `Float::`, `Float` trait bound, `T::zero()`, `T::one()`,
  direct scalar construction, or optional `RealField` extrema residue in
  `problem_validation.rs`/`problem.rs`.
- **Residual risk**: Full cfd-3d FEM provider completion still requires Leto
  storage for FEM matrices/vectors and follow-on Eunomia cleanup in `element`,
  `solver`, `projection_solver`, and `solution`.

---

## Sprint 2026-07-03: cfd-3d FEM Mesh-Helper Eunomia Scalar Boundary
- **Resolved**: `cfd-3d::fem::{mesh_utils,mid_node_cache}` no longer imports
  direct `num_traits::{Float, FromPrimitive}` for P2 extraction, mid-edge
  midpoint matching, tet corner orientation, or mesh-scale extrema. Midpoint
  constants, zero comparison, extrema initialization, and component-wise
  min/max now route through `fem::scalar` plus Eunomia
  `FloatElement`/`NumericElement`.
- **Boundary**: This does not migrate FEM problem validation, element assembly,
  solver/projection assembly, solution scalar identities, or nalgebra
  matrix/vector storage. Non-scalar domain fallbacks in mesh topology ordering
  and index lookup remain outside this scalar-provider slice.
- **Evidence tier**: compile-time integration plus empirical focused FEM
  mesh-helper tests and static residue audit. `cargo fmt -p cfd-3d --check`
  passed. `cargo check -p cfd-3d` passed. `cargo clippy -p cfd-3d --lib --
  -D warnings` passed. `cargo nextest run -p cfd-3d mesh_utils` passed 5/5
  tests. `cargo nextest run -p cfd-3d mid_node_cache` passed 2/2 tests.
  Focused scans found no direct `num_traits`, `FromPrimitive`, `ToPrimitive`,
  `Float::`, `Float` trait bound, old scalar half-construction, `T::zero()`,
  `T::one()`, or `nalgebra::RealField` optional-extrema residue in the touched
  files.
- **Residual risk**: Full cfd-3d FEM provider completion still requires Leto
  storage for FEM matrices/vectors and follow-on Eunomia cleanup in
  `problem_validation`, `element`, `solver`, `projection_solver`, and
  `solution`.

---

## Sprint 2026-07-03: cfd-3d FEM Stabilization/Boundary Eunomia Scalar Boundary
- **Resolved**: `cfd-3d::fem::{stabilization,boundary_classifier}` no longer
  uses direct `num_traits::{FromPrimitive, Float}`, old `T::zero`, `T::one`,
  `nalgebra::RealField::{min_value,max_value}` optional extrema fallbacks, or
  direct `Scalar::from_f64` construction in the touched scalar paths.
  Stabilization parameters, element-size calculations, and z-axial boundary
  classification now route scalar constants and elementary operations through
  `fem::scalar` plus Eunomia `FloatElement`/`NumericElement`.
- **Boundary**: This does not migrate FEM mesh extraction, mid-node cache
  scalar residue, solver/projection assembly, element matrices, or solution
  storage. Those paths still retain nalgebra matrix/vector and adjacent scalar
  provider work.
- **Evidence tier**: compile-time integration plus empirical focused tests and
  static residue audit. `cargo fmt -p cfd-3d --check` passed. `cargo check -p
  cfd-3d` passed. `cargo clippy -p cfd-3d --lib -- -D warnings` passed.
  `cargo nextest run -p cfd-3d stabilization` passed 11/11 tests. `cargo
  nextest run -p cfd-3d boundary` passed 12/12 tests. Focused scans found no
  `num_traits`, `FromPrimitive`, `ToPrimitive`, `Float::`, `T::from_f64`,
  `unwrap_or`, `unwrap_or_else`, `T::zero()`, `T::one()`,
  `RealField::max_value`, `RealField::min_value`, or direct
  `Scalar::from_f64` residue in the touched files.
- **Residual risk**: Full cfd-3d FEM provider completion still requires Leto
  storage for FEM matrices/vectors and follow-on Eunomia cleanup in
  `mesh_utils` and `mid_node_cache`.

---

## Sprint 2026-07-03: cfd-3d FEM Eunomia Scalar Boundary
- **Resolved**: `cfd-3d::fem::{config,quadrature,stress,shape_functions}` no
  longer constructs analytical constants through direct
  `num_traits::FromPrimitive`, `num_traits::Float`, old `T::from_f64`,
  `T::zero`, `T::one`, or conversion fallbacks. A private FEM-local
  `fem::scalar` module is now the scalar-construction SSOT for that bounded
  FEM cone. `cfd-2d::network::reference` and
  `cfd-validation::literature::blood_flow_1d` now carry the Eunomia
  `RealField` bound required by `cfd-1d::NetworkSolveScalar`, closing the
  dev/test build mismatch exposed by the cfd-3d FEM nextest target.
- **Boundary**: This does not migrate FEM matrix/vector storage away from
  nalgebra. `stabilization`, `mesh_utils`, `mid_node_cache`,
  `boundary_classifier`, solver assembly, projection solver, and solution
  storage still contain nalgebra and/or direct scalar-provider residue. The
  cfd-1d `NetworkSolveScalar` trait still preserves `FromPrimitive`,
  `ToPrimitive`, and `num_traits::Float` as a compatibility bridge for the
  current 1D solver boundary.
- **Evidence tier**: compile-time integration plus empirical focused FEM and
  blood-flow validation tests and static residue audit. `cargo check -p
  cfd-3d` passed. `cargo check -p cfd-validation` passed. `cargo clippy -p
  cfd-3d --lib -- -D warnings` passed. `cargo nextest run -p cfd-3d fem`
  passed 40/40 tests. `cargo nextest run -p cfd-validation blood_flow` passed
  2/2 tests. `cargo fmt -p cfd-2d -p cfd-3d -p cfd-validation --check`
  passed. Focused scans found no `num_traits`, `FromPrimitive`,
  `ToPrimitive`, `Float::`, `T::from_f64`, `unwrap_or`, `unwrap_or_else`,
  `T::zero()`, or `T::one()` residue in the touched FEM scalar cone.
  `cargo clippy -p cfd-validation --lib -- -D warnings` was cancelled after
  prolonged shared build-lock contention; no diagnostics were produced before
  cancellation.
- **Residual risk**: Full cfd-3d FEM provider completion still requires Leto
  matrix/vector storage for element matrices, solver/projection assembly, and
  solution vectors, plus follow-on Eunomia cleanup in the remaining FEM scalar
  modules.

---

## Sprint 2026-07-03: cfd-1d Transient Eunomia Boundary
- **Resolved**: `cfd-1d::solver::core::transient::{composition,droplets}` no
  longer uses direct `num_traits::FromPrimitive` construction or
  `num_traits::Float` dispatch in the migrated scalar surface. Timepoint
  tolerances, blood edge-transport constants, segmented advection exponentials,
  absolute-value checks, hematocrit clamps, count-to-scalar averaging, and
  droplet split defaults now route through Eunomia `FloatElement` and
  `NumericElement`. The segmented substep count now uses a checked Eunomia
  scalar conversion and returns a typed configuration error instead of silently
  defaulting to one substep.
- **Boundary**: The transient subsystem still inherits cfd-1d's broader
  nalgebra `RealField` scalar contract. Segmented bifurcation hematocrit
  splitting still calls the current Pries phase-separation API, whose
  algorithm contract consumes `f64`; the scalar conversion into that boundary
  now routes through Eunomia rather than num-traits.
- **Evidence tier**: compile-time integration plus empirical transient
  regression tests and static residue audit. `cargo fmt -p cfd-1d --check`,
  `cargo check -p cfd-1d`, and `cargo clippy -p cfd-1d --lib -- -D warnings`
  passed. `cargo nextest run -p cfd-1d --test transient_composition_parity`
  passed 20/20 tests, `cargo nextest run -p cfd-1d --test
  transient_droplet_parity` passed 9/9 tests, and `cargo nextest run -p
  cfd-1d --test transient_literature_validation` passed 5/5 tests. Focused
  scans found no `ToPrimitive`, `FromPrimitive`, `num_traits::Float`, direct
  `num_traits`, `T::from_f64`, `T::from_usize`, direct `to_usize`, or silent
  substep-count fallback residue in the touched transient composition/droplet
  cone.
- **Residual risk**: Full cfd-1d Eunomia/Leto completion still requires moving
  the nalgebra `RealField` API boundary and deciding whether the Pries
  phase-separation API should expose a Eunomia-generic contract instead of its
  current f64 algorithm boundary.

---

## Sprint 2026-07-03: cfd-1d Anderson cfd-math Boundary
- **Resolved**: `cfd-1d::solver::core` no longer owns a separate nalgebra
  `DMatrix`/LU Anderson least-squares implementation. The nonlinear Picard
  loop constructs one `cfd_math::nonlinear_solver::AndersonAccelerator` per
  solve, and `SolverWorkspace` no longer stores duplicate Anderson residual
  and iterate histories.
- **Boundary**: The 1D linear system still uses nalgebra `DVector` and
  `nalgebra_sparse::CsrMatrix`; this slice converts only at that boundary so
  the Anderson subproblem is owned by the Leto-backed cfd-math implementation.
  `NetworkSolveScalar` records the temporary requirement that primary solver
  callers satisfy both the legacy nalgebra scalar boundary and Eunomia's scalar
  contract.
- **Evidence tier**: compile-time integration plus empirical focused solver
  tests and static residue audit. `cargo fmt -p cfd-1d -p cfd-math --check`
  passed. `cargo check -p cfd-1d` passed. `cargo nextest run -p cfd-1d
  solver` passed 49/49 focused tests. Static scans found no local `DMatrix`,
  `nalgebra::linalg::LU`, SVD wording, or workspace Anderson residual/iterate
  queues in the touched 1D Anderson/workspace path.
- **Residual risk**: Full cfd-1d Leto migration still requires replacing the
  nalgebra `DVector`/`DMatrix`/`nalgebra_sparse` linear-system boundary and
  removing remaining direct `num_traits` conversions in transient composition
  and other domain modules.

---

## Sprint 2026-07-03: cfd-math Anderson Leto Boundary
- **Resolved**: `cfd_math::nonlinear_solver::AndersonAccelerator` no longer
  owns nalgebra `DVector`/`DMatrix` storage. Its public `compute_next` boundary,
  QR columns, R matrix, history buffers, and normal-equations subproblem now
  use Leto `Array1`/`Array2` with Eunomia scalar operations. The
  `cfd-2d::network::coupled` resistance mixer now passes Leto resistance
  vectors directly into Anderson.
- **Boundary**: This closes the cfd-2d network/cfd-math Anderson `DVector`
  seam. The cfd-3d bifurcation, serpentine, trifurcation, and Venturi Picard
  solvers still store FEM velocities as nalgebra `DVector`; those call sites
  now go through a single explicit `cfd-3d::atlas_anderson` conversion
  boundary until `StokesFlowSolution` storage is migrated.
- **Evidence tier**: compile-time integration plus empirical focused
  regression tests and static residue audit. `cargo fmt -p cfd-math -p cfd-2d
  -p cfd-3d --check` passed. `cargo check -p cfd-math`, `cargo check -p
  cfd-2d`, and `cargo check -p cfd-3d` passed. `cargo clippy -p cfd-math
  --lib -- -D warnings` passed. `cargo nextest run -p cfd-math anderson`
  passed 5/5 focused tests. `cargo nextest run -p cfd-2d network` passed
  22/22 focused tests. Static scans found no `nalgebra`, `DVector`,
  `DMatrix`, `num_traits`, or SVD residue in
  `cfd-math/src/nonlinear_solver/anderson.rs`, and no `DVector` residue in
  `cfd-2d/src/network/coupled.rs`.
- **Residual risk**: Full Leto CPU migration still requires the cfd-3d FEM
  solution velocity storage boundary, remaining cfd-math linear/time-stepping
  `DVector` APIs, cfd-2d validation benchmark/LES-DES helper residue, and
  direct GPU replacement through Hephaestus.

---

## Sprint 2026-07-03: cfd-2d Constants Validation Eunomia Boundary
- **Resolved**: `cfd-2d::physics::turbulence::constants_validation` no longer
  depends on nalgebra `RealField` or num-traits scalar construction/conversion
  at the migrated boundary. `TurbulenceConstantsValidator`,
  `ConstantsValidationResult`, and DNS sensitivity analysis now use Eunomia
  `RealField`/`FloatElement`/`NumericElement` for constants, scalar
  construction, square roots, absolute values, maxima, and display conversion.
- **Boundary**: This closes the local constants-validation scalar seam only.
  It does not claim full cfd-2d provider migration. Remaining turbulence and
  solver provider work includes validation benchmark/LES-DES helpers, the
  cfd-2d network/cfd-math Anderson `DVector` seam, the adjacent
  `cfd-validation` RSM MMS scalar oracle, and direct GPU replacement through
  Hephaestus.
- **Evidence tier**: compile-time integration plus empirical focused
  regression tests and static residue audit. `cargo fmt -p cfd-2d --check`
  passed. `cargo check -p cfd-2d` passed. `cargo nextest run -p cfd-2d
  constants` passed 11/11 focused tests. `cargo nextest run -p cfd-2d rans`
  passed 5/5 focused tests. `cargo nextest run -p cfd-2d macroscopic` passed
  4/4 focused tests. Focused scan over `constants_validation` found no
  `nalgebra::RealField`, `num_traits`, `FromPrimitive`, `ToPrimitive`,
  fallible `T::from_f64`, legacy `T::zero()`, or `T::one()` residue.
- **Residual risk**: `cargo clippy -p cfd-2d --all-targets -- -D warnings`
  still fails on broader all-target lint debt outside this slice, including
  example/test field reassignment, manual multiple-of checks, map identity,
  format-arg, and LBM boundary identity-op cleanups.

---

## Sprint 2026-07-03: cfd-2d Shared RANS Leto/Eunomia Boundary
- **Resolved**: The shared `cfd-2d::physics::turbulence` RANS model
  scalar/vector contract no longer depends on nalgebra vectors or
  num-traits-style scalar construction at the migrated boundary.
  `TurbulenceModel<T>` now binds on Eunomia `RealField`; k-epsilon,
  realizable C_mu, k-omega SST, Spalart-Allmaras, wall boundary conditions,
  wall treatment, RANS validation, and SIMPLE turbulence coupling use Eunomia
  scalar construction/math with Leto `Vector2` velocity updates.
- **Boundary**: This closes the shared RANS/Spalart-Allmaras vector contract.
  Remaining cfd-2d turbulence provider work is in validation scalar helpers
  and direct GPU execution through Hephaestus. The adjacent `cfd-validation`
  RSM MMS scalar oracle is a separate validation-crate migration seam.
- **Evidence tier**: compile-time integration plus empirical value-semantic
  turbulence tests and static residue scans. `cargo fmt -p cfd-2d --check`
  passed. `cargo check -p cfd-2d` passed. `cargo nextest run -p cfd-2d
  k_epsilon` passed 39/39 focused tests. `cargo nextest run -p cfd-2d
  k_omega` passed 19/19 focused tests. `cargo nextest run -p cfd-2d
  spalart_allmaras` passed 21/21 focused tests. `cargo nextest run -p
  cfd-2d wall` passed 51/51 focused tests. `cargo nextest run -p cfd-2d
  rans` passed 5/5 focused tests. Focused scans over the migrated RANS
  boundary report no nalgebra vector residue and no nalgebra/num-traits scalar
  residue except the explicitly transitional validation compatibility bound.
- **Residual risk**: `cargo clippy -p cfd-2d --all-targets -- -D warnings`
  was not completed because unrelated active builds held the shared
  `D:\atlas\target` lock. Full turbulence migration still requires
  Hephaestus GPU replacement and cleanup of the validation scalar-helper
  compatibility seams.

---

## Sprint 2026-07-03: cfd-2d Reynolds-Stress Eunomia Boundary
- **Resolved**: Native `cfd-2d::physics::turbulence::reynolds_stress` scalar
  helper code no longer depends on nalgebra `RealField` or
  `num_traits::FromPrimitive`. RSM model construction, tensor storage,
  production, diffusion, curvature, wall-reflection, pressure-strain, and
  optimized transport math now use Eunomia `RealField` constants,
  construction, and math while preserving the Leto `Array2` storage boundary.
- **Boundary**: This closes the native RSM scalar-helper seam only. The
  shared RANS/Spalart-Allmaras contract is closed by the later Leto/Eunomia
  item; the adjacent `cfd-validation` RSM MMS oracle still has a separate
  scalar-provider migration seam.
- **Evidence tier**: compile-time integration plus empirical value-semantic
  Reynolds-stress tests and static residue scan. `cargo fmt -p cfd-2d`
  completed. `cargo check -p cfd-2d` passed. `cargo nextest run -p cfd-2d
  reynolds_stress` passed 10/10 focused tests. Focused scan over the RSM
  native boundary reports no `FromPrimitive`, direct `nalgebra::RealField`,
  `T::zero()`, `T::one()`, `T::max_value`, legacy generic `.max(T::ZERO)`,
  `DMatrix`, or tuple-indexing residue.
- **Residual risk**: Full cfd-2d turbulence migration still requires the
  adjacent `cfd-validation` RSM MMS scalar oracle, validation scalar-helper
  cleanup, and direct GPU execution through Hephaestus.

---

## Sprint 2026-07-03: cfd-2d Reynolds-Stress Leto Boundary
- **Resolved**: `cfd-2d::physics::turbulence::reynolds_stress` no longer
  exposes nalgebra `DMatrix` at its Reynolds-stress tensor, transport, helper,
  or validation-fixture boundary. `ReynoldsStressTensor`, initialization and
  dissipation setup, production and diffusion helpers, transport velocity and
  scalar helpers, co-located tests, comprehensive Reynolds-stress validation
  fixtures, and the adjacent `cfd-validation` MMS L2-error oracle now use Leto
  `Array2`.
- **Boundary**: This closes the Reynolds-stress matrix storage and validation
  seam only. The broader RANS/Spalart-Allmaras contract is closed by the later
  Leto/Eunomia item; the adjacent `cfd-validation` RSM MMS scalar oracle still
  has a separate scalar-provider migration seam, and raw GPU execution still
  requires the Hephaestus replacement.
- **Evidence tier**: compile-time integration plus empirical value-semantic
  Reynolds-stress tests and static residue scan. `cargo fmt -p cfd-2d -p
  cfd-validation --check` passed. `cargo check -p cfd-2d` passed. `cargo
  nextest run -p cfd-2d reynolds_stress` passed 10/10 focused tests. `cargo
  nextest run -p cfd-validation reynolds_stress` passed 9/9 focused tests.
  Focused scan over the migrated RSM boundary reports no `DMatrix`,
  `DMatrix::`, nalgebra shape methods, `from_element`, or tuple-indexing
  residue.
- **Residual risk**: Full cfd-2d turbulence migration still requires the
  adjacent `cfd-validation` RSM MMS scalar oracle, validation scalar-helper
  cleanup, and direct GPU execution through Hephaestus.

---

## Sprint 2026-07-03: cfd-2d WALE Leto/Eunomia Boundary
- **Resolved**: `cfd-2d::physics::turbulence::les_smagorinsky::wale` no
  longer exposes `Field2D<nalgebra::Vector2<T>>` or direct
  `num_traits::FromPrimitive` at its velocity-field/scalar-provider boundary.
  WALE SGS viscosity and boundary-gradient recovery now consume Leto `Array2`
  velocity component fields and use Eunomia `RealField`/`NumericElement`
  scalar operations.
- **Boundary**: This closes the standalone WALE LES model seam only. The
  broader RANS `TurbulenceModel` contract is closed by the later Leto/Eunomia
  item; the adjacent `cfd-validation` RSM MMS scalar oracle and raw GPU
  execution remain separate seams.
- **Evidence tier**: compile-time integration plus empirical value-semantic
  WALE tests and static residue scan. `cargo fmt -p cfd-2d` completed.
  `cargo check -p cfd-2d` passed. `cargo nextest run -p cfd-2d wale` passed
  4/4 focused tests, including mismatched velocity-component shape rejection.
  Focused scan over the touched WALE file reports no `Field2D`, `Vector2`,
  `nalgebra`, `num_traits`, `FromPrimitive`, old identity methods, fallible
  scalar-construction residue, or direct method `sqrt`.
- **Residual risk**: Full cfd-2d turbulence migration still requires the
  adjacent `cfd-validation` RSM MMS scalar oracle, validation scalar-helper
  cleanup, and direct GPU replacement through Hephaestus.

---

## Sprint 2026-07-03: cfd-2d Shared LES/DES Leto Boundary
- **Resolved**: The shared `cfd-2d::physics::turbulence` LES/DES matrix
  boundary no longer exposes nalgebra `DMatrix`. `LESTurbulenceModel`,
  Smagorinsky model state, strain/viscosity/dynamic helpers, GPU helper
  boundary, DES model fields, DES length-scale utilities, and LES/DES
  validation and benchmark call sites now use Leto `Array2`.
- **Boundary**: This closes the shared LES/DES matrix API and state seam only.
  The RANS/Spalart-Allmaras velocity contract is closed by the later
  Leto/Eunomia item, and raw GPU execution still needs the broader Hephaestus
  replacement. The adjacent `cfd-validation` RSM MMS scalar oracle remains a
  separate turbulence migration slice.
- **Evidence tier**: compile-time integration plus empirical value-semantic
  LES/DES tests and static residue scan. `cargo fmt -p cfd-2d --check`
  passed. `cargo check -p cfd-2d` passed. `cargo nextest run -p cfd-2d
  les_smagorinsky` passed 43/43 focused tests. `cargo nextest run -p cfd-2d
  des` passed 14/14 focused tests. Focused scan over the touched LES/DES files
  reports no `DMatrix`, `DMatrix::`, or nalgebra matrix shape-method residue.
- **Residual risk**: Full cfd-2d turbulence migration still requires the
  adjacent `cfd-validation` RSM MMS scalar oracle, validation scalar-helper
  cleanup, and direct GPU replacement through Hephaestus.

---

## Sprint 2026-07-03: cfd-2d MILES Leto Boundary
- **Resolved**: `cfd-2d::physics::turbulence::les_smagorinsky::miles` no
  longer exposes nalgebra `DMatrix` or direct `num_traits::FromPrimitive` at
  its standalone 2x2 velocity-gradient API boundary. Config defaults now use
  Eunomia scalar constants, shock detection and numerical flux use Leto
  `Array2`, and co-located value tests cover the migrated provider seam.
- **Boundary**: This closes only the standalone MILES LES model seam. The
  shared LES/DES matrix, WALE, and RANS boundaries are now separately
  migrated; the adjacent `cfd-validation` RSM MMS scalar oracle, validation
  scalar-helper cleanup, and raw GPU execution still own provider migration
  work.
- **Evidence tier**: compile-time integration plus empirical value-semantic
  MILES tests and static residue scan. `cargo fmt -p cfd-2d --check` passed.
  `cargo check -p cfd-2d` passed. `cargo nextest run -p cfd-2d miles` passed
  7/7 focused tests. Focused scan over the touched file reports no `DMatrix`,
  `nalgebra`, `num_traits`, `FromPrimitive`, old identity methods, or
  fallible scalar-construction residue.
- **Residual risk**: Full cfd-2d turbulence migration still requires the
  adjacent `cfd-validation` RSM MMS scalar oracle, validation scalar-helper
  cleanup, and the remaining direct GPU path through Hephaestus.

---

## Sprint 2026-07-03: cfd-2d Sigma/Vreman Leto Boundary
- **Resolved**: `cfd-2d::physics::turbulence::les_smagorinsky::{sigma,
  vreman}` no longer expose nalgebra `DMatrix` or direct
  `num_traits::FromPrimitive` at their standalone 2x2 velocity-gradient and
  SGS-stress API boundary. Config defaults now use Eunomia scalar constants,
  viscosity/stress/invariant calculations use Leto `Array2`, and co-located
  value tests cover the migrated provider seam.
- **Boundary**: This closes only the Sigma/Vreman standalone LES model seam.
  The shared LES/DES matrix, WALE, and RANS boundaries are now separately
  migrated; the adjacent `cfd-validation` RSM MMS scalar oracle, validation
  scalar-helper cleanup, and raw GPU execution still own provider migration
  work.
- **Evidence tier**: compile-time integration plus empirical value-semantic
  Sigma/Vreman tests and static residue scan. `cargo fmt -p cfd-2d --check`
  passed. `cargo check -p cfd-2d` passed. `cargo nextest run -p cfd-2d
  sigma` passed 7/7 focused tests. `cargo nextest run -p cfd-2d vreman`
  passed 6/6 focused tests. Focused scans over the two touched files report
  no `DMatrix`, `nalgebra`, `num_traits`, `FromPrimitive`, old identity
  methods, or fallible scalar-construction residue.
- **Residual risk**: Full cfd-2d turbulence migration still requires the
  adjacent `cfd-validation` RSM MMS scalar oracle, validation scalar-helper
  cleanup, and the remaining direct GPU path through Hephaestus.

---

## Sprint 2026-07-03: cfd-3d Poisson Leto Boundary
- **Resolved**: `cfd-3d::spectral::{poisson,solver}` no longer exposes
  nalgebra `DVector` at the public Poisson source/solution boundary.
  `PoissonProblem::source_term`, `SpectralSolution::u`, and
  `PoissonSolver::solve` now use Leto `Array1`, and the validation/property/
  domain/robustness tests plus the spectral Poisson example exercise that
  provider-native seam.
- **Boundary**: This closes the public Poisson storage boundary only. The
  current Poisson dense operator assembly and LU solve still use nalgebra
  `DMatrix`/`DVector` internally; `spectral::chebyshev` remains a separate
  nalgebra vector/matrix migration slice.
- **Evidence tier**: compile-time integration plus empirical analytical
  Poisson tests and static residue scan. `cargo fmt -p cfd-3d --check`
  passed. `cargo check -p cfd-3d` passed. `cargo nextest run -p cfd-3d
  --test poisson_validation` passed 7/7 tests. `cargo nextest run -p cfd-3d
  spectral_poisson` passed 4/4 focused tests. `cargo nextest run -p cfd-3d
  test_poisson_zero_rhs_zero_solution` passed 1/1 focused test. Focused scan
  over the touched Poisson public boundary files shows no `source_term:
  DVector`, `Result<DVector`, or `solve(&DVector` residue.
- **Residual risk**: The later 2026-07-04 Chebyshev and Poisson provider
  slices close the internal spectral dense-matrix/vector ownership called out
  here. Full objective completion still requires migrating the upstream
  `cfd-core::SolverConfig` scalar bound, remaining workspace nalgebra/
  num-traits surfaces, and raw GPU/provider code.

---

## Sprint 2026-07-03: cfd-3d Fourier Leto/Eunomia Boundary
- **Resolved**: `cfd-3d::spectral::fourier` no longer exposes nalgebra
  `DVector` or nalgebra `Complex` at the transform boundary. The Apollo-backed
  Fourier wrapper now accepts and returns Leto `Array1` values and uses Eunomia
  `Complex` for spectral coefficients. The direct Fourier validation suite,
  internal adversarial spectral round-trip tests, and the library Fourier smoke
  test now exercise the provider-native seam.
- **Boundary**: This closes the Fourier transform/spectral-derivative seam.
  The later 2026-07-04 Chebyshev and Poisson provider slices close the
  spectral dense array/matrix ownership that remained here; `spectral::solver`
  still carries a nalgebra scalar alias through upstream `cfd-core::SolverConfig`.
- **Evidence tier**: compile-time integration plus empirical analytical
  Fourier tests and static residue scan. `cargo fmt -p cfd-3d --check` passed.
  `cargo check -p cfd-3d` passed. `cargo nextest run -p cfd-3d --test
  fourier_validation` passed 12/12 tests. Focused scan over the touched Fourier
  boundary files reports no legacy `DVector`/nalgebra-complex residue.
- **Residual risk**: Full objective completion still requires replacing the
  upstream `cfd-core::SolverConfig` nalgebra scalar contract, broader workspace
  nalgebra/num-traits surfaces, and raw wgpu/GPU code with
  Leto/Eunomia/Hephaestus-backed surfaces.

---

## Sprint 2026-07-03: Moirai Text Residue and Eunomia Bound Propagation
- **Resolved**: Direct Rayon/Tokio crate-source references are gone from
  `crates/`. Stale `rayon::par_iter`/Rayon guidance in `cfd-optim`,
  `cfd-math`, `cfd-3d`, and `cfd-validation` now names Moirai or removes the
  obsolete guard wording. `cfd-2d::network::reference` and
  `cfd-validation::literature::blood_flow_1d` now use Eunomia
  `FloatElement`/`NumericElement` helpers for the scalar constants and f64
  extraction paths touched by the cfd-core fluid provider migration.
- **Boundary**: This closes source-level Rayon/Tokio residue, not every
  concurrency architecture concern. The workspace still has broad
  nalgebra/num-traits surfaces and direct GPU/wgpu surfaces that must be
  migrated through Leto/Eunomia and Hephaestus in later coordinated slices.
- **Evidence tier**: static source audit plus compile-time integration and
  empirical focused tests. `cargo fmt -p cfd-math --check`, `cargo fmt -p
  cfd-3d --check`, `cargo fmt -p cfd-optim --check`, and `cargo fmt -p cfd-2d
  -p cfd-validation --check` passed. `cargo check -p cfd-math`, `cargo check
  -p cfd-validation`, and `cargo check -p cfd-optim` passed. `cargo nextest
  run -p cfd-2d reference_trace` passed 3/3 tests, `cargo nextest run -p
  cfd-validation blood_flow` passed 2/2 tests, and `cargo nextest run -p
  cfd-optim` passed 121/121 tests. Source scans report no Rayon/Tokio or
  `par_iter`-style hits under `crates/`, and no old `T::from_f64`,
  `num_traits::Zero::zero`, or optional `to_f64().unwrap*` conversions in the
  two touched scalar-boundary files.
- **Residual risk**: Full objective completion still requires replacing
  nalgebra/ndarray storage and scalar contracts with Leto/Eunomia across the
  solver and validation crates, and replacing direct wgpu/GPU code with
  Hephaestus-backed GPU surfaces.

---

## Sprint 2026-07-03: cfd-core Fluid Identity Dispatch
- **Resolved**: `cfd-core::physics::fluid::{traits,validation,temperature,
  properties,newtonian}` no longer uses `T::zero()`, `T::one()`, or method
  `sqrt()` in the touched scalar paths. Fluid-state Mach-number, validation,
  temperature-dependent fluid, base property, and Newtonian/ideal-gas logic now
  route identities and square root through Eunomia `NumericElement`.
- **Boundary**: This is a scalar API dispatch cleanup, not the full fluid trait
  replacement. The touched files still import `nalgebra::RealField` because
  `FluidTrait<T>`, `FluidState<T>`, and several public fluid structs remain
  consumed by downstream crates under the existing nalgebra scalar contract.
- **Evidence tier**: compile-time integration plus empirical fluid and
  downstream resistance tests with static source audit. `cargo fmt -p
  cfd-core --check` passed. `cargo check -p cfd-core --no-default-features`
  passed. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core
  --no-default-features fluid` passed 45/45 tests with 153 skipped. `cargo
  check -p cfd-1d` passed. `cargo nextest run -p cfd-1d resistance` passed
  173/173 tests with 555 skipped. Focused scan over the touched fluid files
  reports no `T::zero()`, `T::one()`, or method `sqrt()` hits.
- **Residual risk**: Full fluid-module scalar-provider migration still
  requires moving the shared `FluidTrait<T>`/`FluidState<T>` and remaining
  public fluid structs off `nalgebra::RealField` with downstream cfd-1d,
  cfd-2d, cfd-3d, and validation bound updates in the same coordinated slice.

---

## Sprint 2026-07-03: cfd-core Non-Newtonian Eunomia Scalars
- **Resolved**: `cfd-core::physics::fluid::non_newtonian` no longer imports
  direct `num_traits::FromPrimitive`, calls `T::from_f64`, uses old
  `T::zero()`/`T::one()` identity methods, or contains silent
  scalar-conversion fallbacks. Bingham, Casson, Carreau-Yasuda,
  Herschel-Bulkley, and Power-law models now use Eunomia
  `FloatElement`/`NumericElement` for constants, identities, square roots,
  powers, and exponentials.
- **Boundary**: This preserves the existing public fluid model and
  `FluidTrait<T>` API. `nalgebra::RealField` remains in the fluid module
  because the shared fluid trait scalar contract is still a cross-crate
  boundary consumed by cfd-1d, cfd-2d, cfd-3d, and validation crates.
- **Evidence tier**: compile-time integration plus empirical value-semantic
  non-Newtonian and downstream resistance tests with static source audit.
  `cargo fmt -p cfd-core --check` passed. `cargo check -p cfd-core
  --no-default-features` passed. `cargo check -p cfd-core` passed. `cargo
  nextest run -p cfd-core --no-default-features non_newtonian` passed 5/5
  tests with 193 skipped. `cargo check -p cfd-1d` passed. `cargo nextest run
  -p cfd-1d resistance` passed 173/173 tests with 555 skipped. Focused source
  scan reports no direct scalar-provider or fallback-conversion residue under
  `crates/cfd-core/src/physics/fluid/non_newtonian`.
- **Residual risk**: Full fluid-module scalar-provider migration still
  requires moving `FluidTrait<T>`, `FluidState<T>`, and remaining fluid
  structs off `nalgebra::RealField` after downstream generic bounds are
  updated together.

---

## Sprint 2026-07-03: cfd-1d Resistance Scalar Facade
- **Resolved**: `cfd-1d::physics::resistance` no longer carries direct
  per-model `nalgebra::RealField`/`num_traits::FromPrimitive` imports, direct
  `T::from_f64` constants, or silent
  `nalgebra::try_convert(...).unwrap_or(...)` conversion fallbacks. The
  resistance calculator, factory, geometry helpers, and model implementations
  route scalar constants and f64 extraction through `ResistanceScalar`,
  `scalar_from_f64`, and `scalar_to_f64`.
- **Boundary**: This is a centralized-provider slice, not a full Eunomia
  replacement. `physics/resistance/models/traits.rs` still bridges through
  `nalgebra::RealField` and `num_traits::FromPrimitive` because the current
  `cfd_core::physics::fluid::FluidTrait<T>` contract requires
  `T: nalgebra::RealField + Copy`.
- **Evidence tier**: compile-time integration plus empirical resistance tests
  and static source audit. `cargo fmt -p cfd-1d` passed. `cargo check -p
  cfd-1d` passed. `cargo nextest run -p cfd-1d resistance` passed 173/173
  tests with 555 skipped. A focused residue scan over
  `crates/cfd-1d/src/physics/resistance` reports provider hits only in
  `models/traits.rs`.
- **Residual risk**: Full Eunomia replacement for resistance requires migrating
  the shared `cfd_core::physics::fluid::FluidTrait<T>` scalar contract first.

---

## Sprint 2026-07-03: cfd-3d Turbulence Model Eunomia Scalars
- **Resolved**: `cfd-3d::physics::turbulence` model implementations for
  Smagorinsky, Dynamic Smagorinsky, Dynamic Gradient Smagorinsky, AMD, Vreman,
  WALE, Sigma, DES, Spalart-Allmaras, Mixing Length, and shared SGS energy
  helpers now use Eunomia provider contracts for scalar constants and
  transcendental/math dispatch. The Sigma closed-form eigensolver now calls
  `FloatElement::acos`, supplied by Eunomia in this slice.
- **Evidence tier**: compile-time provider integration, value-semantic
  turbulence tests, and static source audit. `cargo check -p cfd-3d` passed.
  `cargo nextest run -p cfd-3d turbulence` passed 46/46. Touched-crate
  rustfmt passed. Focused scans across the edited turbulence model files found
  no direct `nalgebra`, `num_traits`, `FromPrimitive`, `RealField`,
  `T::one()`, `T::zero()`, or `Float::` residue.
- **Residual risk**: The `cfd-3d::physics::turbulence` scalar-provider
  migration is closed. Broader CFDrs provider migration remains incomplete
  outside this module, including other direct nalgebra/ndarray/GPU/execution
  surfaces.

## Sprint 2026-07-03: cfd-3d Turbulence Helper Eunomia Scalars
- **Resolved**: `cfd-3d::physics::turbulence::field_ops` no longer imports
  or bounds helper operations on direct `nalgebra::RealField`,
  `num_traits::FromPrimitive`, or `num_traits::Float`. Derivative, velocity
  gradient, strain, matrix-contract, vorticity, and trace-free helper paths
  now use Leto `Vector3` and Eunomia `NumericElement` identities and square
  root.
- **Resolved**: `cfd-3d::physics::turbulence::filter_ops` no longer imports
  or bounds filter helpers on direct nalgebra or num-traits scalar contracts.
  Filter moments use Eunomia `NumericElement` zero/one values, and stencil
  counts use checked `i32::try_from` followed by Eunomia `CastFrom<i32>`.
- **Evidence tier**: compile-time provider integration, value-semantic
  turbulence tests, and static source audit. `cargo check -p cfd-3d` passed.
  `cargo nextest run -p cfd-3d turbulence` passed 46/46. Touched-crate
  rustfmt passed. Focused scans found no `nalgebra`, `num_traits`,
  `FromPrimitive`, `RealField`, old `T::one()`/`T::zero()`, or `Float::`
  residue in `field_ops.rs` and `filter_ops.rs`. Touched-file `git diff
  --check` reported only LF-to-CRLF normalization warnings.
- **Residual risk**: Closed by the later model-level turbulence scalar-provider
  pass; no direct scalar-provider residue remains under
  `crates/cfd-3d/src/physics/turbulence`.

## Sprint 2026-07-03: Fluid Dynamics VelocityField Leto Storage
- **Resolved**: `cfd-core::physics::fluid_dynamics::VelocityField<T>` now
  stores Leto `geometry::Vector3<T>` components instead of nalgebra vectors.
  `FlowField`, `VelocityField`, `PressureField`, `ScalarField`, RANS trait
  surfaces, and turbulence trait surfaces use Eunomia `NumericElement` where
  they depend on the migrated field representation.
- **Resolved**: Direct 3D spectral DNS/forcing/diagnostics consumers,
  turbulence field/filter helpers, IBM NUFFT field construction boundaries,
  and 3D Taylor-Green, forced-turbulence, and NUFFT-coupling validation
  benchmark consumers now construct Leto-backed `VelocityField` values.
  Validation analytical vectors and IBM marker/probe APIs that remain
  nalgebra-based convert at the boundary instead of leaking nalgebra into the
  shared velocity-field storage contract.
- **Evidence tier**: compile-time provider integration, targeted
  value-semantic/differential validation tests, and static source audit.
  `cargo check -p cfd-core`, `cargo check -p cfd-3d`, and `cargo check -p
  cfd-validation` passed. `cargo nextest run -p cfd-core
  fluid_dynamics::operations` passed 3/3. `cargo nextest run -p cfd-3d
  spectral` passed 41/41. `cargo nextest run -p cfd-3d nufft` passed 2/2.
  `cargo nextest run -p cfd-validation taylor_green` passed 19/19. `cargo
  nextest run -p cfd-validation forced_turbulence` passed 2/2. `cargo
  nextest run -p cfd-validation nufft_coupling` passed 2/2. Touched-crate
  rustfmt passed. Focused residue scans found no direct nalgebra
  `Vector3`/`RealField` residue in the migrated core and spectral files.
  `git diff --check` reported only pre-existing LF-to-CRLF normalization
  warnings.
- **Residual risk**: Broader CFDrs Leto migration remains incomplete. Direct
  nalgebra vectors and dense matrix/vector types remain in unrelated domains,
  including FEM tests/solver surfaces, VOF, cascade, validation analytical
  and error-metric vectors, IBM marker/probe APIs, and larger dense
  solver/operator paths.

## Sprint 2026-07-03: cfd-core Fluid Dynamics Operations Eunomia Scalars
- **Resolved**: `cfd-core::physics::fluid_dynamics::operations` no longer
  imports or bounds flow operations on direct `num_traits::FromPrimitive`.
  Vorticity, divergence, kinetic-energy, and enstrophy operations now use
  Eunomia `NumericElement` for scalar zero/one identities and half/two
  constants.
- **Boundary**: This slice preserves the current nalgebra `Vector3` and
  `VelocityField<T>` storage boundary, and preserves the existing Moirai
  parallel-slice execution. It is a scalar-provider migration, not the full
  Leto vector/storage replacement.
- **Evidence tier**: compile-time provider integration, value-semantic flow
  operation tests, and static source audit. `cargo check -p cfd-core` passed.
  `cargo nextest run -p cfd-core fluid_dynamics::operations` passed 3/3.
  Touched-file rustfmt passed. Focused scans found no direct `num_traits`,
  `FromPrimitive`, `T::from_f64`, or old `unwrap_or_else` scalar-fallback
  patterns in `operations.rs`.
- **Residual risk**: Broader CFDrs Atlas provider migration remains
  incomplete. The local scalar-construction holdout under `fluid_dynamics` is
  closed, but full Leto replacement still requires migrating `VelocityField<T>`,
  `operations.rs`, RANS/turbulence traits, and Rhie-Chow away from nalgebra
  vector/storage contracts.

## Sprint 2026-07-03: cfd-core Fluid Dynamics Service Eunomia Scalars
- **Resolved**: `cfd-core::physics::fluid_dynamics::service` no longer imports
  or bounds pipe-flow scalar formulas on direct `num_traits::{Float,
  FromPrimitive}`. Reynolds/Prandtl helpers, pipe pressure-drop,
  friction-factor selection, Colebrook-White iteration, and Haaland evaluation
  now use Eunomia `FloatElement`/`NumericElement` for constants, powers,
  square roots, logarithms, and absolute convergence checks.
- **Boundary**: This slice preserves `ConstantPropertyFluid<T>`'s existing
  nalgebra scalar/storage boundary. It is a scalar-provider migration for the
  service formulas, not a full Leto replacement of fluid storage.
- **Evidence tier**: compile-time provider integration, value-semantic
  service tests, and static source audit. `cargo check -p cfd-core` passed.
  `cargo nextest run -p cfd-core fluid_dynamics::service` passed 2/2.
  Touched-file rustfmt passed. Focused scans found no direct `num_traits`,
  `FromPrimitive`, `Float::`, `T::from_f64`, or old `unwrap_or_else`
  scalar-fallback patterns in `service.rs`.
- **Residual risk**: Broader CFDrs Atlas provider migration remains
  incomplete. Local fluid-dynamics scalar holdout remains in `operations.rs`;
  full Leto replacement also requires migrating `VelocityField<T>` and vector
  operations away from nalgebra storage.

## Sprint 2026-07-03: cfd-core Flow Regime Eunomia Scalars
- **Resolved**: `cfd-core::physics::fluid_dynamics::flow_regimes` no longer
  imports direct `nalgebra::RealField` or `num_traits::ToPrimitive`. Reynolds,
  Mach, and combined classifier paths now use Eunomia `RealField` plus
  `NumericElement::to_f64`.
- **Resolved**: The old silent `unwrap_or(0.0)` conversion fallback was
  removed from flow classification. The
  `FluidDynamicsService::flow_regime` wrapper now accepts the same Eunomia
  real-scalar contract as the classifier.
- **Boundary**: This slice is limited to flow-regime classification and its
  service wrapper. `FluidDynamicsService` still has broader direct
  `num_traits::{Float, FromPrimitive}` usage in pipe pressure-drop and
  friction-factor formulas, and other fluid-dynamics modules still retain
  nalgebra vector/storage boundaries.
- **Evidence tier**: compile-time provider integration, value-semantic
  threshold tests, and static source audit. `cargo check -p cfd-core` passed.
  `cargo nextest run -p cfd-core flow_regime` passed 3/3. Touched-file
  rustfmt passed. Focused scans found no direct `num_traits`, `ToPrimitive`,
  direct `nalgebra::RealField`, `.to_f64()`, or `unwrap_or(0.0)` residue in
  `flow_regimes.rs`.
- **Residual risk**: Broader CFDrs Atlas provider migration remains
  incomplete. The next local fluid-dynamics scalar holdouts are
  `operations.rs` and the service pipe-flow friction-factor path; full Leto
  replacement also requires migrating `VelocityField<T>` and vector operations
  away from nalgebra storage.

## Sprint 2026-07-03: cfd-3d Spectral Diagnostics Eunomia Conversion
- **Resolved**: `cfd-3d::spectral::diagnostics` no longer imports or bounds
  local scalar conversion on direct `num_traits::{FromPrimitive, ToPrimitive}`.
  The helper that stages velocity components into Leto `Array3<f64>` buffers
  for Apollo FFTs now uses Eunomia `NumericElement::to_f64`.
- **Boundary**: This slice is limited to
  `crates/cfd-3d/src/spectral/diagnostics.rs`. The remaining
  `nalgebra::RealField` and `nalgebra::Vector3` references are inherited from
  `cfd-core::physics::fluid_dynamics::VelocityField<T>`, which still owns the
  current velocity storage boundary.
- **Evidence tier**: compile-time provider integration, existing
  value-semantic spectral diagnostics tests, and static source audit. `cargo
  check -p cfd-3d` passed. `cargo nextest run -p cfd-3d diagnostics` passed
  5/5. Touched-file rustfmt passed. Focused scans found no direct
  `num_traits`, `FromPrimitive`, `ToPrimitive`, `.to_f64()`, or obsolete
  `component_name` residue in `diagnostics.rs`.
- **Residual risk**: Broader CFDrs Atlas provider migration remains
  incomplete. `cfd-core::VelocityField<T>` still stores nalgebra vectors, and
  larger `cfd-3d` spectral Chebyshev/Poisson/solver surfaces still use
  nalgebra dense matrix/vector types.

## Sprint 2026-07-03: cfd-2d MRT/Carreau-Yasuda Eunomia Scalars
- **Resolved**: `RelaxationMatrix<T>`, `MrtCollision<T>`,
  `CarreauYasudaModel<T>`, and `CarreauYasudaBgk<T>` no longer depend on
  direct `nalgebra::RealField`, `num_traits::{Float, FromPrimitive}`, generic
  `T::from_f64`, `T::zero`, `T::one`, or direct `Float::sqrt`/`Float::powf`
  dispatch. Constants, zero/one identities, square roots, real powers, max
  clamps, MRT moment matrices, and local tau construction now route through
  Eunomia-backed scalar helpers and `CastFrom<i32>` for D2Q9 velocity
  components.
- **Boundary**: This closes the local `cfd-2d` LBM scalar-provider residue
  recorded by the previous LBM seam slice. Broader `cfd-2d` and CFDrs
  `nalgebra`/`num_traits` residues remain outside LBM.
- **Evidence tier**: compile-time provider integration, value-semantic LBM and
  Carreau-Yasuda tests, and static source audit. `cargo check -p cfd-2d`
  passed. `cargo nextest run -p cfd-2d lbm` passed 31/31. `cargo nextest run
  -p cfd-2d carreau_yasuda` passed 5/5. Touched-file rustfmt passed. Focused
  scans found no direct `nalgebra::RealField`, `num_traits`, `FromPrimitive`,
  `T::from_f64`, `T::zero`, `T::one`, `from_i32`, or `Float::` residue in the
  migrated MRT/Carreau-Yasuda files.
- **Residual risk**: The broader CFDrs provider migration remains incomplete.
  Direct `nalgebra`, `num_traits`, `ndarray`, raw `wgpu`, `tokio`, and `rayon`
  residues remain outside this bounded LBM scalar-provider closeout.

## Sprint 2026-07-03: cfd-2d LBM Eunomia Scalar Seam
- **Resolved**: `cfd-2d` LBM macroscopic extraction, D2Q9 equilibrium,
  collision trait bounds, and BGK collision constants no longer depend on
  direct `nalgebra::RealField` or `num_traits::FromPrimitive`. The migrated
  seam uses Eunomia `FloatElement` and local scalar helpers for density,
  velocity, pressure, stress, kinetic energy, vorticity, equilibrium weights,
  lattice constants, BGK viscosity, and omega construction.
- **Boundary**: This slice is limited to the authoritative LBM
  macroscopic/equilibrium/BGK seam. MRT and Carreau-Yasuda
  `CollisionOperator` impls now satisfy the migrated trait bound; the later
  MRT/Carreau-Yasuda follow-up closes their broader inherent moment and
  rheology helper scalar-provider residue.
- **Evidence tier**: compile-time provider integration, value-semantic LBM
  tests, and static source audit. `cargo check -p cfd-2d` passed. `cargo
  nextest run -p cfd-2d lbm` passed 31/31. Touched-file rustfmt and
  `git diff --check` passed. Focused scans found no direct `nalgebra`,
  `RealField`, `num_traits`, `FromPrimitive`, or `ToPrimitive` residue in
  `macroscopic.rs`, `lattice.rs`, `collision/traits.rs`, or
  `collision/bgk.rs`.
- **Residual risk**: The broader CFDrs provider migration remains incomplete.
  Direct `nalgebra`, `num_traits`, `ndarray`, raw `wgpu`, `tokio`, and `rayon`
  residues remain outside this bounded LBM slice.

## Sprint 2026-07-03: cfd-3d Wall Functions Eunomia Scalars
- **Resolved**: `cfd-3d::physics::turbulence::wall_functions` no longer
  imports or bounds local scalar math on direct `nalgebra::RealField`,
  `num_traits::FromPrimitive`, or `num_traits::Float`. Spalding wall-law
  constants now use Eunomia `FloatElement::from_f64`; `exp`, `ln`, `sqrt`,
  `abs`, and max clamps use Eunomia `FloatElement`/`NumericElement`.
- **Boundary**: This slice is limited to
  `crates/cfd-3d/src/physics/turbulence/wall_functions.rs`. Other turbulence
  models still retain direct `nalgebra` and `num_traits` residues.
- **Evidence tier**: compile-time provider integration, value-semantic
  wall-law tests, and static source audit. `cargo check -p cfd-3d` passed.
  `cargo nextest run -p cfd-3d wall_functions` passed 3/3. Touched-file
  rustfmt and `git diff --check` passed. Focused scans found no direct
  `nalgebra`, `RealField`, `num_traits`, `FromPrimitive`, or `ToPrimitive`
  residue in `wall_functions.rs`.
- **Residual risk**: Direct `GpuContext` replacement with Hephaestus was
  inspected but not patched in this slice because CFDrs resolves both
  `wgpu@0.19.4` and Hephaestus' `wgpu@26.0.1`; replacing raw `wgpu` device
  fields requires a broader WGPU API-version alignment. The broader CFDrs
  provider migration remains incomplete outside this file.

## Sprint 2026-07-03: cfd-3d Apollo Fourier Eunomia Bounds
- **Resolved**: The Apollo-backed `cfd-3d::spectral::fourier` wrapper no
  longer imports or bounds local scalar conversion on direct
  `num_traits::{FromPrimitive, ToPrimitive}`. Fourier input/output conversion
  now uses Eunomia `NumericElement::to_f64` and `FloatElement::from_f64`,
  while retaining the existing Apollo FFT + Leto array execution path and
  legacy `DVector` compatibility boundary.
- **Boundary**: This slice also propagated the required `FloatElement` bound
  through the downstream `cfd-2d` LBM/Poiseuille/network and `cfd-validation`
  3D blood benchmark call sites that consume migrated `cfd-core` blood or
  cavitation APIs, plus `cfd-3d` DES turbulence blood-model viscosity usage.
- **Evidence tier**: compile-time provider integration, value-semantic Fourier
  integration tests, and static source audit. `cargo check -p cfd-3d` passed.
  `cargo check -p cfd-validation` passed. `cargo nextest run -p cfd-3d --test
  fourier_validation` passed 12/12. Touched-file rustfmt passed. Focused scans
  found no direct `num_traits`, `FromPrimitive`, or `ToPrimitive` residue in
  `crates/cfd-3d/src/spectral/fourier.rs` or
  `crates/cfd-3d/src/physics/turbulence/des.rs`.
- **Residual risk**: This is not a full CFDrs scalar-provider closeout.
  Remaining direct `num_traits`, `nalgebra`, direct `wgpu`, and execution
  provider residues exist outside the touched files and remain later migration
  work.

## Sprint 2026-07-03: cfd-core Fåhræus-Lindqvist Eunomia Scalars
- **Resolved**: `FahraeuasLindqvist<T>` no longer depends on direct
  `num_traits::FromPrimitive`, direct generic `T::from_f64()`, direct
  `T::zero()`, direct `T::one()`, direct generic `powf()`, direct generic
  `exp()`, or direct generic `abs()` for local scalar construction and
  microvascular viscosity formulas. Pries/Secomb exponent formulas, `mu_45`,
  relative-viscosity clamping, and tube hematocrit now route scalar constants
  and math through Eunomia `FloatElement`/`NumericElement`.
- **Boundary**: This slice is limited to
  `crates/cfd-core/src/physics/fluid/blood/fahraeus_lindqvist.rs`. Together
  with the prior Cross/Casson/Carreau slices, local blood-model scalar
  `num_traits` construction is closed. The broader `Fluid<T>` trait still
  carries an inherited `nalgebra::RealField` boundary for later provider work.
- **Evidence tier**: compile-time provider integration, existing
  value-semantic blood tests, and static source audit. `cargo check -p
  cfd-core` passed. `cargo nextest run -p cfd-core fahraeus_lindqvist` passed
  3/3. `cargo nextest run -p cfd-core blood` passed 24/24. Touched-file
  rustfmt and touched-file `git diff --check` passed. Focused
  `fahraeus_lindqvist.rs` residue scan found no direct `num_traits`,
  `FromPrimitive`, generic `T::from_f64`, `T::zero`, `T::one`, generic `powf`,
  generic `exp`, or generic `abs` residue. Broader blood residue scan now
  matches only concrete `f64` helper/test expressions.
- **Residual risk**: `cargo clippy -p cfd-core --all-targets -- -D warnings`
  remains blocked by pre-existing unrelated lints in
  `crates/cfd-core/src/physics/boundary/applicator.rs` and
  `crates/cfd-core/src/physics/fluid_dynamics/rhie_chow.rs`.

## Sprint 2026-07-03: cfd-core Casson/Carreau Blood Eunomia Scalars
- **Resolved**: `CassonBlood<T>`, `CarreauYasudaBlood<T>`, and
  `BloodModel<T>` no longer require direct `num_traits::FromPrimitive` for
  local scalar construction or model dispatch. Casson constants, hematocrit
  scaling, temperature correction, square-root apparent-viscosity formula, and
  validation constants now route through Eunomia `FloatElement`/
  `NumericElement`. Carreau-Yasuda constants, zero/one identities, and real
  powers remain on Eunomia after the dispatch-bound fix.
- **Boundary**: This slice is limited to
  `crates/cfd-core/src/physics/fluid/blood/{casson,carreau_yasuda,mod}.rs`.
  The concrete `temperature_viscosity_factor(f64)` helper remains intentionally
  concrete. Fåhræus-Lindqvist and the broader `Fluid<T>` trait still carry
  residual blood-fluid provider migration work.
- **Evidence tier**: compile-time provider integration, existing
  value-semantic blood tests, and static source audit. `cargo check -p
  cfd-core` passed. `cargo nextest run -p cfd-core casson` passed 12/12.
  `cargo nextest run -p cfd-core carreau_yasuda` passed 4/4. `cargo nextest
  run -p cfd-core blood` passed 24/24. Touched-file rustfmt and touched-file
  `git diff --check` passed. Focused touched-file residue scan found no direct
  `num_traits`, `FromPrimitive`, generic `T::from_f64`, `T::zero`, `T::one`,
  generic `powf`, generic `sqrt`, or generic `exp` residue; the only remaining
  match is the intentionally concrete `f64` temperature helper.
- **Residual risk**: `cargo clippy -p cfd-core --all-targets -- -D warnings`
  remains blocked by pre-existing unrelated lints in
  `crates/cfd-core/src/physics/boundary/applicator.rs` and
  `crates/cfd-core/src/physics/fluid_dynamics/rhie_chow.rs`.

## Sprint 2026-07-03: cfd-core Cross Blood Eunomia Scalars
- **Resolved**: `CrossBlood<T>` no longer depends on direct
  `num_traits::FromPrimitive`, direct `T::from_f64()`, direct `T::zero()`,
  direct `T::one()`, or direct generic `powf()` for its local constructor and
  apparent-viscosity formula. Scalar constants and math now route through
  Eunomia `FloatElement`/`NumericElement`.
- **Boundary**: This slice is limited to
  `crates/cfd-core/src/physics/fluid/blood/cross.rs`. The surrounding
  `Fluid<T>` trait still carries `nalgebra::RealField`, so CrossBlood retains
  that inherited trait bound until the fluid trait is migrated. Other blood
  models still retain direct `num_traits` construction.
- **Evidence tier**: compile-time provider integration, existing
  value-semantic blood tests, and static source audit. `cargo check -p
  cfd-core` passed. `cargo nextest run -p cfd-core cross` passed 1/1 test.
  `cargo nextest run -p cfd-core blood` passed 24/24 tests. Focused
  `cross.rs` residue scan found no `num_traits`, `FromPrimitive`, direct
  `T::from_f64`, direct `T::zero`, direct `T::one`, or direct `powf` residue.
  Touched-file rustfmt and touched-file `git diff --check` passed.

## Sprint 2026-07-03: cfd-core Cavitation Scalar Closeout
- **Resolved**: `CavitationModel<T>` and `ZgbParams<T>` no longer depend on
  `nalgebra::RealField`, direct `num_traits::FromPrimitive`, direct
  `T::from_f64()`, direct `T::zero()`, or direct `T::one()` construction.
  Mass-transfer constants and math now route through Eunomia
  `FloatElement`/`NumericElement`.
- **Resolved**: `heterogeneous_inception_threshold_pa` no longer exposes a
  fake generic `RealField` adapter while converting every input to `f64` and
  converting the result back. The legacy coarse adapter now uses the concrete
  `f64` contract already used by the selective cavitation model.
- **Coverage added**: Value-semantic tests cover Kunz vaporization and
  condensation rates, Schnerr-Sauer zero-density behavior, ZGB degenerate
  bubble-radius rejection, default model construction, and legacy
  heterogeneous-adapter parity with the equivalent selective-cavitation input.
- **Boundary**: This slice is limited to
  `crates/cfd-core/src/physics/cavitation/{models,heterogeneous_nucleation}.rs`.
  The full CFDrs Atlas migration remains incomplete outside cavitation.
- **Evidence tier**: compile-time provider integration, value-semantic focused
  tests, and static source audit. `cargo check -p cfd-core` passed. `cargo
  nextest run -p cfd-core cavitation` passed 40/40 tests. Focused cavitation
  residue scan found no `nalgebra`, `RealField`, `num_traits`,
  `FromPrimitive`, `ToPrimitive`, direct `T::from_f64`, direct `T::zero`,
  direct `T::one`, `to_subset`, or `try_convert` residue. Touched-file rustfmt
  and touched-file `git diff --check` passed.

## Sprint 2026-07-03: cfd-core Nuclei Transport Eunomia Scalars
- **Resolved**: `nuclei_adjusted_vapor_pressure`,
  `NucleiTransportConfig<T>`, and `NucleiTransport<T>` no longer depend on
  `nalgebra::RealField`, direct `T::from_f64()`, direct `T::zero()`, direct
  `T::one()`, or direct generic `exp()`. Nuclei transport scalar constants,
  zero/one identities, and exponential decay now route through Eunomia
  `FloatElement`/`NumericElement`.
- **Boundary**: This slice is limited to
  `crates/cfd-core/src/physics/cavitation/nuclei_transport.rs`. Remaining
  cavitation scalar holdouts were then `models.rs` and
  `heterogeneous_nucleation.rs`; the current cavitation closeout removed both
  holdouts.
- **Evidence tier**: compile-time provider integration, existing
  value-semantic focused tests, and static source audit. `cargo check -p
  cfd-core` passed. `cargo nextest run -p cfd-core cavitation` passed 35/35
  tests. Focused `nuclei_transport.rs` residue scan found no `nalgebra`,
  `RealField`, `num_traits`, `FromPrimitive`, `ToPrimitive`, direct
  `T::from_f64`, direct `T::zero`, direct `T::one`, or direct generic `exp`
  residue. Touched-file rustfmt and touched-file `git diff --check` passed.

## Sprint 2026-07-03: cfd-core Venturi Cavitation Eunomia Scalars
- **Resolved**: `VenturiCavitation<T>` no longer depends on
  `nalgebra::RealField`, `num_traits::FromPrimitive`, direct `T::zero()`,
  direct `T::one()`, or direct `T::from_f64()` construction. Venturi scalar
  constants and math now route through Eunomia `FloatElement`/`NumericElement`.
- **Boundary**: This slice is limited to
  `crates/cfd-core/src/physics/cavitation/venturi.rs`. Remaining cavitation
  scalar holdouts were `models.rs`, `nuclei_transport.rs`, and
  `heterogeneous_nucleation.rs`; the later nuclei-transport slice removed
  `nuclei_transport.rs` from the active residual list.
- **Evidence tier**: compile-time provider integration, existing
  value-semantic focused tests, and static source audit. `cargo check -p
  cfd-core` passed. `cargo nextest run -p cfd-core cavitation` passed 35/35
  tests. Focused migrated-production scan found no `nalgebra`, `RealField`,
  `num_traits`, `FromPrimitive`, `ToPrimitive`, direct `T::from_f64`, direct
  `T::zero`, direct `T::one`, direct `powf`, direct `powi`, or direct `tan`
  residue in `venturi.rs`; the only `powi` hit is an existing `f64` test
  oracle expression. Touched-file rustfmt and touched-file `git diff --check`
  passed.

## Sprint 2026-07-03: cfd-core Cavitation Eunomia Scalar Cone
- **Resolved**: The touched cavitation Rayleigh-Plesset, biological damage,
  regime-analysis, cavitation-number, and material-damage surfaces no longer
  depend on `nalgebra::RealField`, direct `num_traits::FromPrimitive`, direct
  `T::zero()`, direct `T::one()`, or direct `T::from_f64()` construction.
  Scalar constants, powers, square roots, exponentials, finite checks, and
  scalar min/max now route through Eunomia `FloatElement`/`NumericElement`.
- **Coverage added**: Closed-form value tests now cover cavitation-number
  definition, zero-velocity large-index behavior, pressure-recovery scaling,
  Hammitt erosion pressure-ratio power, Rayleigh collapse impact pressure, and
  pit-depth hardness normalization.
- **Boundary**: This slice is limited to
  `crates/cfd-core/src/physics/cavitation/{rayleigh_plesset,bio_damage,number,damage}.rs`
  and `crates/cfd-core/src/physics/cavitation/regimes/`. Remaining cavitation
  scalar holdouts at that point were `models.rs`, `venturi.rs`,
  `nuclei_transport.rs`, and `heterogeneous_nucleation.rs`; the later Venturi
  slice removed `venturi.rs` from the active residual list.
- **Evidence tier**: compile-time provider integration, value-semantic focused
  tests, and static source audit. `cargo check -p cfd-core` passed. `cargo
  nextest run -p cfd-core cavitation` passed 35/35 tests. Focused residue scan
  across the migrated cavitation files found no `nalgebra`, `RealField`,
  `num_traits`, `FromPrimitive`, `ToPrimitive`, direct `T::from_f64`, direct
  `T::from_u64`, direct `T::zero`, direct `T::one`, direct `powf`, or direct
  `powi` residue. Touched-file rustfmt and touched-file `git diff --check`
  passed. Full `cargo clippy -p cfd-core --all-targets -- -D warnings`
  remains blocked by unrelated existing lints in
  `physics/boundary/applicator.rs` and `physics/fluid_dynamics/rhie_chow.rs`.

## Sprint 2026-07-03: cfd-core Hemolysis Eunomia Scalars
- **Resolved**: `HemolysisCalculator<T>` and `PlateletActivation<T>` no longer
  depend on `nalgebra::RealField`, `num_traits::FromPrimitive`, direct
  `T::zero()`, or direct `T::one()`. Hemolysis scalar constants and zero/one
  identities now route through Eunomia `FloatElement`/`NumericElement`.
- **Boundary**: This slice is limited to `cfd-core::physics::hemolysis`.
  Broader `cfd-core` compute, boundary, cavitation, fluid, mesh, and
  fluid-dynamics modules still retain nalgebra contracts for later Leto,
  Gaia, and Eunomia provider slices.
- **Evidence tier**: compile-time provider integration, focused
  value-semantic tests, and static source audit. `cargo check -p cfd-core`
  passed. `cargo nextest run -p cfd-core hemolysis` passed 9/9 tests. Focused
  hemolysis residue scan found no `nalgebra`, `RealField`, `num_traits`,
  `FromPrimitive`, `ToPrimitive`, `T::zero`, or `T::one` residue. Touched-file
  rustfmt and touched-file `git diff --check` passed.

## Sprint 2026-07-03: cfd-2d Grid/FVM Leto-Eunomia Boundary
- **Resolved**: `Grid2D::cell_center`, `StructuredGrid2D`,
  `UnstructuredGrid2D`, and `AdaptiveGrid2D::effective_resolution` no longer
  depend on nalgebra vector/scalar contracts or direct num-traits scalar
  fallbacks. Grid centers use `leto::geometry::Vector2`, structured grid
  spacing and cell-center construction use Eunomia-backed helpers, and
  `FvmSolver<T>` no longer declares `nalgebra::RealField`.
- **Consumer updates**: Momentum rotating-wall boundary code, FDM/LBM
  consumers, and validation benchmark/cross-fidelity code that read grid
  centers now use Leto vector indexing instead of nalgebra `.x`/`.y` fields.
- **Boundary**: This slice closes the grid/FVM inherited `RealField` blocker
  but does not claim full `cfd-2d` nalgebra removal. Momentum assembly,
  pressure-velocity, SIMPLE/PIMPLE, and other solver surfaces still retain
  nalgebra/num-traits contracts for later Leto/Gaia/Eunomia provider slices.
  `tests_poisson_mms.rs` contains an MMS test function, but `cargo nextest
  list -p cfd-2d --lib` does not list it because the file is not wired into the
  library module tree; the zero-test filter is not used as evidence.
- **Evidence tier**: compile-time provider/consumer integration, focused
  value-semantic tests, and static source audit. `cargo check -p cfd-2d` and
  `cargo check -p cfd-validation` passed. `cargo nextest run -p cfd-2d --lib
  grid` passed 21/21 tests. `cargo nextest run -p cfd-2d --lib fvm` passed
  25/25 tests. Focused grid/FVM residue scans found no `RealField`,
  `FromPrimitive`, `ToPrimitive`, `num_traits`, `nalgebra::Vector2`, direct
  `T::from_f64`, direct `T::from_usize`, or stale scalar fallback residue.
  Touched-file rustfmt and touched-file `git diff --check` passed.

## Sprint 2026-07-03: cfd-2d FVM Leto Vector Boundary
- **Resolved**: `cfd-2d::solvers::fvm::Face` no longer stores nalgebra
  `Vector2`, and `FvmSolver::solve` now accepts `Field2D<leto::geometry::Vector2<T>>`.
  FVM face normal normalization and face flux dot products route through Leto's
  fixed-vector provider surface.
- **Provider extension**: Leto geometry now exposes `Vector2<T>` as the
  fixed-size 2-D vector alias and implements generic fixed-vector
  `norm_squared`, `norm`, `try_normalize`, and `normalize` methods.
- **Boundary**: This slice did not migrate `StructuredGrid2D<T>`; the later
  grid/FVM provider slice closes that inherited `RealField` blocker. The
  focused FVM scan shows no nalgebra `Vector2` import or direct num-traits
  scalar residue in the migrated FVM files.
- **Evidence tier**: compile-time provider/consumer integration plus
  value-semantic focused tests and static source audit. `cargo check -p leto`
  passed. `cargo check -p cfd-2d` passed. `cargo nextest run -p leto
  fixed_vector_norm_and_normalization_are_value_semantic` passed 1/1 test.
  `cargo nextest run -p cfd-2d --lib fvm` passed 25/25 focused tests.

## Sprint 2026-07-03: cfd-2d FVM Eunomia Scalar Constants
- **Resolved**: `cfd-2d::solvers::fvm` no longer uses direct
  `num_traits::FromPrimitive`, direct `T::from_f64`, direct `T::from_usize`,
  or the previous silent `unwrap_or(1e-3)` diffusion fallback. `FvmConfig`
  default grid spacing, timestep, convergence tolerance, CFL number,
  relaxation factor, and diffusion coefficient now use Eunomia
  `FloatElement::from_f64`. `FvmSolver` carries `FloatElement` on its owner
  bounds because it stores `FvmConfig<T>` and uses Eunomia-backed grid-index
  conversion for face centers. Flux schemes use Eunomia scalar constants, and
  mixed Eunomia/nalgebra scalar method overlaps explicitly dispatch `abs` via
  `NumericElement` and `powf` via `FloatElement`.
- **Coverage added**: `default_fvm_config_preserves_physical_constants`
  asserts the default grid dimensions and scalar constants directly.
  `test_power_law_rejects_nonfinite_diffusion` and
  `test_hybrid_rejects_nonfinite_diffusion` assert the new typed error
  contract for invalid diffusion coefficients.
- **Boundary**: This slice only removes FVM direct num-traits scalar
  construction/fallbacks. `RealField` remains for the existing nalgebra
  vector, ordering, absolute value, max, and power contracts; `Face<T>` still
  stores nalgebra `Vector2`.
- **Evidence tier**: compile-time dependency-chain verification plus
  value-semantic focused tests. Focused scan of
  `crates/cfd-2d/src/solvers/fvm` shows no
  `num_traits`, `FromPrimitive`, `ToPrimitive`, direct `T::from_f64`, direct
  `T::from_usize`, `.to_subset()`, or stale diffusion fallback residue, and
  touched-file rustfmt is clean. `cargo check -p cfd-2d` and `cargo check -p
  cfd-3d` passed, and `cargo nextest run -p cfd-2d --lib fvm` passed 25/25
  focused tests.

## Sprint 2026-07-03: cfd-2d CFL Eunomia Scalars
- **Resolved**: `cfd-2d::stability::CFLCalculator` no longer imports
  `nalgebra::RealField` or `num_traits::FromPrimitive`. CFL scalar constants,
  absolute values, and literal construction now use Eunomia
  `FloatElement`/`NumericElement`.
- **Boundary**: This slice only covers the standalone CFL calculator.
  Additional `cfd-2d` stability/time-scheme, turbulence, validation, and solver
  surfaces still contain nalgebra/num-traits usage.
- **Evidence tier**: compile-time integration plus empirical focused tests and
  static source audit. Focused CFL scan shows no `RealField`, `num_traits`,
  `FromPrimitive`, `ToPrimitive`, `T::from_f64`, `T::zero`, or `T::one`
  residue. `cargo check -p cfd-2d` passed, and `cargo nextest run -p cfd-2d
  --lib cfl` passed 9/9 tests. A broader non-library nextest filter hit
  pre-existing `simplec_pimple_validation` generic-bound debt unrelated to this
  file. No runtime performance claim is made.

## Sprint 2026-07-03: cfd-core Material Eunomia Traits
- **Resolved**: `SolidProperties`, `InterfaceProperties`, `ElasticSolid`,
  `WettingProperties`, and `FluidSolidInterface` no longer depend on
  `nalgebra::RealField`; scalar math and constant construction now use Eunomia
  `FloatElement`/`NumericElement`.
- **Boundary**: `MaterialDatabase` still requires `RealField` because it stores
  `Box<dyn Fluid<T>>`; fluid, hemolysis, fluid-dynamics, boundary, geometry,
  mesh, and solver nalgebra surfaces remain open.
- **Evidence tier**: compile-time integration plus empirical focused tests and
  static source audit. Focused material scan shows no `nalgebra::RealField`,
  `RealField`, `num_traits`, `FromPrimitive`, `ToPrimitive`, or `Float` residue
  in the migrated solid/interface files. `cargo check -p cfd-core` passed, and
  `cargo nextest run -p cfd-core material` passed 4/4 tests. No runtime
  performance claim is made.

## Sprint 2026-07-03: cfd-core Velocity Leto Vector
- **Resolved**: `cfd-core::physics::values::Velocity` now stores
  `leto::geometry::Vector3<T>` instead of `nalgebra::Vector3<T>`, and its
  generic contract is Eunomia `FloatElement`/`NumericElement` instead of
  `nalgebra::RealField`. `PhysicalParameters::gravity` now uses the same Leto
  vector type and no longer requires `RealField` for its own methods.
- **Provider extension**: Leto geometry now derives Serde for `Point2`,
  `Point3`, `Vector3`, `UnitVector3`, and `Isometry3`, preserving CFDrs'
  serialized value-object boundary while using provider-owned vector storage.
- **Boundary**: `ProblemAggregate` and `SimulationAggregate` still carry
  `RealField` because their `Domain<T>` and fluid-property contracts have not
  been migrated in this slice. Material, hemolysis, and fluid-dynamics
  nalgebra-bound surfaces remain open.
- **Evidence tier**: compile-time integration plus static source audit.
  Touched-file rustfmt passed. Focused residue scan shows no nalgebra
  `Vector3`/`RealField` in `Velocity` or `PhysicalParameters`. `cargo check -p
  cfd-core` passed after compiling the modified local Leto provider. `cargo
  nextest run -p cfd-core --lib` passed 183/183 tests. Downstream `cargo check
  -p cfd-2d`, `cargo check -p cfd-3d`, and `cargo check -p cfd-validation`
  passed. No runtime performance claim is made.

## Sprint 2026-07-03: cfd-core Physics Value Eunomia Scalars
- **Resolved**: Scalar-only physics value wrappers no longer use
  `nalgebra::RealField` or `nalgebra::ComplexField`. `Temperature`,
  `Pressure`, `ReynoldsNumber`, and `DimensionlessNumber` now depend on
  Eunomia `FloatElement`/`NumericElement` for scalar construction, zero,
  absolute value, and square root. Immediate aggregate owners now declare the
  same bound where they store these value wrappers.
- **Boundary**: `Velocity` remains on `nalgebra::Vector3`, so its `RealField`
  bound is intentionally preserved for the Leto/Gaia vector replacement slice.
  Additional `RealField` use in material, hemolysis, and broader aggregate
  contracts remains open.
- **Evidence tier**: compile-time integration plus static source audit.
  Touched-file rustfmt passed. A focused scan of the four scalar value-wrapper
  files returns no `nalgebra::RealField`, `RealField`, or `ComplexField`
  matches. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core
  --lib` passed 183/183 tests. No runtime performance claim is made.
- **Residual risk**: Full nalgebra removal from `cfd-core` still requires
  replacing vector/matrix contracts and material/hemolysis scalar bounds with
  Atlas-owned providers.

## Sprint 2026-07-03: cfd-math WGPU Boundary Reduction
- **Resolved**: `cfd-math` no longer declares a direct optional `wgpu`
  dependency or imports WGPU symbols to collect GPU dispatch metrics. The root
  package `gpu` feature no longer activates `dep:wgpu` directly. GPU
  synchronization and timestamp-query capability now flow through
  `cfd_core::compute::gpu::GpuContext`.
- **Boundary**: This slice does not replace the raw WGPU kernels. It removes
  downstream WGPU ownership from the math crate so `cfd-core` remains the single
  current raw-WGPU boundary being migrated toward Hephaestus.
- **Evidence tier**: compile-time integration and static manifest/source audit
  plus empirical focused tests. A focused source/manifest scan shows no direct
  `wgpu`, `dep:wgpu`, `TIMESTAMP_QUERY`, or `Maintain::Wait` residue in
  `cfd-math`; `cargo check -p cfd-math --features gpu` passed; `cargo nextest
  run -p cfd-math --features gpu --lib` passed 298/298 tests; `cargo check -p
  cfd-suite --features gpu` passed; inverse WGPU trees show `wgpu@0.19.4` and
  `wgpu@26.0.1` reach `cfd-math` only through `cfd-core`. No runtime GPU
  performance claim is made.
- **Residual risk**: `cfd-core::compute::gpu` still owns raw WGPU buffers,
  command encoders, pipelines, and WGSL kernels. Full GPU replacement still
  requires moving those surfaces onto Hephaestus abstractions.

## Sprint 2026-07-03: cfd-core Hephaestus GPU Probe
- **Resolved**: `cfd-core::compute::ComputeBackend::detect_gpu_support` no
  longer performs raw `wgpu::Instance` adapter probing itself. The `gpu`
  feature now activates `hephaestus-wgpu`, and availability is determined by
  `hephaestus_wgpu::WgpuDevice::try_default`, making Hephaestus the active
  device-acquisition provider for the compute backend.
- **Dependency audit**: `cargo tree --workspace -i hephaestus-wgpu` shows
  `hephaestus-wgpu v0.11.0 (D:\atlas\repos\hephaestus\crates\hephaestus-wgpu)
  -> cfd-core`. `cargo tree --workspace -i rustfft` reports no matching
  package. `cargo tree --workspace -i ndarray` still resolves only through
  `numpy -> cfd-python`.
- **Evidence tier**: compile-time integration plus empirical focused tests and
  static dependency audit. `cargo check -p cfd-core --features gpu` passed.
  `cargo nextest run -p cfd-core --features gpu --lib` passed 183/183 tests.
- **Residual risk**: CFDrs still owns raw WGPU buffers, command encoders,
  pipelines, and WGSL kernels in `cfd-core::compute::gpu`. Full GPU replacement
  still requires moving those kernel surfaces onto Hephaestus `WgpuDevice`/
  `WgpuBuffer`/kernel abstractions or completing the recorded WGPU version
  normalization.

## Sprint 2026-07-03: cfd-1d Non-Python ndarray Path Removal
- **Resolved**: `cfd-1d` no longer declares the unused `sprs` dependency that
  pulled `ndarray v0.17.2` into the active 1D/3D dependency graph. The root
  workspace no longer declares unused `ndarray` as a shared dependency.
  cfd-1d tests that construct `ConstantPropertyFluid::water_20c()` now declare
  the required `eunomia::FloatElement` bound explicitly.
- **Dependency audit**: `cargo tree -p cfd-1d -i ndarray` and `cargo tree -p
  cfd-3d -i ndarray` report no matching package. `cargo tree --workspace -i
  ndarray` shows the remaining workspace path as `ndarray v0.16.1 -> numpy
  v0.22.1 -> cfd-python`.
- **Evidence tier**: manifest/lock static audit plus compile-time integration
  and empirical focused tests. Touched-file rustfmt passed. `cargo update -p
  sprs` removed `sprs`, `ndarray v0.17.2`, and stale transitive packages.
  `cargo check -p cfd-1d` passed. `cargo nextest run -p cfd-1d` passed 725/725
  tests with 3 skipped.
- **Residual risk**: cfd-1d still uses `nalgebra`/`nalgebra-sparse`; replacing
  those surfaces with Leto sparse/vector types remains open. Workspace
  `ndarray` remains through the Python `numpy` boundary, which needs a separate
  binding API decision.

## Sprint 2026-07-03: cfd-3d Apollo ndarray Removal
- **Resolved**: CFDrs no longer resolves `apollo-fft`/`apollo-nufft` from the
  older ndarray-backed Apollo Git revision. The workspace now patches Apollo,
  Leto, Eunomia, Mnemosyne, Moirai, Hermes, and Melinoe to the side-by-side
  Atlas checkouts required by the provider migration. cfd-3d FEM scalar
  constants and spectral complex scalar multiplication now use Eunomia
  contracts, and cfd-validation propagates the required `FloatElement` bounds
  for the cfd-3d nextest dev-dependency graph.
- **Dependency audit**: `cargo tree -p apollo-fft | rg ndarray` and `cargo
  tree -p apollo-nufft | rg ndarray` return no matches. A source/manifest/lock
  scan of `D:\atlas\repos\apollo` returns no `ndarray` hits. The remaining
  active `ndarray` in `cargo tree -p cfd-3d -i ndarray` is via
  `sprs -> cfd-1d`, not Apollo.
- **Evidence tier**: compile-time integration plus empirical test execution and
  static dependency/source audit. Touched-file rustfmt passed. `cargo check -p
  cfd-3d` passed. `cargo nextest run -p cfd-3d` passed 394/394 tests; one
  existing mesh-convergence test was reported slow at 16.394s.
- **Residual risk**: Leto sparse/operator replacement and full Hephaestus GPU
  kernel replacement remain separate provider-migration slices. The non-Apollo
  `sprs -> cfd-1d` ndarray path remained open at this Apollo slice boundary,
  and Apollo still contained rustfft/num-complex validation and benchmark
  residue outside this ndarray-removal slice.

## Sprint 2026-07-03: cfd-math Eunomia Geometric Multigrid
- **Resolved**: `cfd-math::linear_solver::preconditioners::multigrid::gmg`
  no longer imports direct `num_traits::FromPrimitive` or
  `cfd_core::conversion::SafeFromF64`, and no longer uses direct
  `T::from_f64(...)`, `T::from_usize(...)`, or `T::from_f64_or_one(...)`
  scalar construction. AMG's stale direct `FromPrimitive` bound/import was
  removed after the multigrid provider contracts no longer required it.
- **Boundary**: This slice preserves existing nalgebra dense/sparse/vector
  multigrid surfaces. Leto replacement remains a later provider migration
  slice.
- **Evidence tier**: compile-time integration plus empirical focused GMG/AMG
  tests and static source audit. Touched-file rustfmt passed. `cargo check -p
  cfd-math` passed. `cargo nextest run -p cfd-math geometric_multigrid
  poisson_matrix restrict_residual fas_solve amg` passed 11/11 tests.
  Multigrid-wide static scan found no direct `num_traits`, `FromPrimitive`,
  `ToPrimitive`, direct `T::from_f64`, direct `T::from_usize`, conversion
  fallback, `from_f64_or`, `SafeFromF64`, stale `rayon`/`tokio`, `rustfft`,
  `ndarray`, or `num_complex` hits in `multigrid`.
- **Residual risk**: The broader `cfd-math` and workspace still contain
  nalgebra surfaces and non-multigrid provider migration work. Raw GPU paths
  still need Hephaestus migration, and broader CPU array/vector surfaces still
  need Leto replacement.

## Sprint 2026-07-03: cfd-math Eunomia Multigrid Interpolation
- **Resolved**: `cfd-math::linear_solver::preconditioners::multigrid::interpolation`
  has been changed to remove direct `num_traits::{FromPrimitive, ToPrimitive}`
  imports and direct scalar conversion/fallback paths. Interpolation scalar
  constants, index-distance conversion, quality row-sum extraction, constant
  preservation error, sparsity ratio, and absolute-value dispatch now route
  through Eunomia helpers.
- **Boundary**: This slice preserves the existing nalgebra `DVector` and
  sparse/vector surfaces. GMG transfer, AMG residual bounds, and full Leto
  migration remain separate provider slices.
- **Evidence tier**: compile-time integration plus empirical focused
  interpolation tests and static source audit. Touched-file rustfmt passed.
  `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math
  interpolation` passed 14/14 tests. Static scan found no direct `num_traits`,
  `FromPrimitive`, `ToPrimitive`, direct `T::from_f64`, direct
  `T::from_usize`, conversion fallback, `from_f64_or`, `SafeFromF64`, stale
  `rayon`, direct `as f64`, or `.to_f64()` fallback hits in
  `multigrid/interpolation.rs`.
- **Residual risk**: `multigrid/gmg` and `amg.rs` still retain direct provider
  residue; nalgebra sparse/vector surfaces remain open for Leto migration.

## Sprint 2026-07-03: cfd-math Eunomia Multigrid Coarsening
- **Resolved**: `cfd-math::linear_solver::preconditioners::multigrid::coarsening`
  no longer imports direct `num_traits::{FromPrimitive, ToPrimitive}` or uses
  direct `T::from_f64(...)`, `T::from_usize(...)`, `T::max_value(...)`, or
  `to_f64().unwrap_or(...)` scalar conversion/fallback paths. Coarsening
  constants, count-to-scalar conversion, max-distance sentinels, quality
  extraction, and strength-matrix absolute-value dispatch now route through
  Eunomia.
- **Boundary**: This slice preserves existing nalgebra sparse/vector surfaces.
  It does not migrate multigrid interpolation, GMG transfer, restriction, or
  AMG's remaining direct `FromPrimitive` bound.
- **Evidence tier**: compile-time integration plus empirical focused
  value-semantic coarsening tests and static source audit. Touched-file rustfmt
  passed. `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math
  coarsening` passed 10/10 tests. Static scan found no direct `num_traits`,
  `FromPrimitive`, `ToPrimitive`, direct `T::from_f64`, direct
  `T::from_usize`, conversion fallback, `from_f64_or`, `SafeFromF64`, or stale
  `rayon` hits in `multigrid/coarsening`.
- **Residual risk**: `multigrid/interpolation.rs`, `multigrid/gmg`, and
  `amg.rs` still retain direct provider residue. The broader cfd-math
  multigrid API still uses nalgebra sparse/vector types pending Leto migration.

## Sprint 2026-07-03: cfd-math Eunomia Multigrid Smoothers
- **Resolved**: `cfd-math::linear_solver::preconditioners::multigrid::smoothers`
  no longer uses direct `T::from_f64(...).unwrap_or_else(...)` conversion
  fallbacks or nalgebra absolute-value dispatch for smoother diagonal
  thresholds, Chebyshev eigenvalue defaults, Chebyshev recurrence constants, or
  smoother update thresholds. The touched AMG owner paths now construct
  coarsening thresholds, smoother relaxation values, and complexity filters
  through Eunomia.
- **Boundary**: This slice preserves the existing nalgebra `DVector` and
  `SparseMatrix` surfaces and keeps AMG's `FromPrimitive` bound because deeper
  coarsening/interpolation contracts still require it.
- **Evidence tier**: compile-time integration plus empirical focused
  value-semantic smoother tests and static source audit. Touched-file rustfmt
  passed. `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math
  test_gauss_seidel_smoother test_jacobi_smoother
  test_symmetric_gauss_seidel test_sor_smoother test_chebyshev_smoother`
  passed 5/5 tests. Static scan found no direct scalar-conversion fallback,
  stale `SafeFromF64`, `from_f64_or`, direct `T::from_usize`, or stale `rayon`
  hits in the touched smoother/AMG files, and no direct `num_traits` provider
  residue in `smoothers.rs`.
- **Residual risk**: `amg.rs` still imports and bounds direct
  `num_traits::FromPrimitive` for deeper coarsening/interpolation routines.
  Other multigrid modules, the raw GPU operator path, and nalgebra sparse/vector
  surfaces remain open provider-migration work.

## Sprint 2026-07-03: cfd-math Eunomia Convergence Monitor
- **Resolved**: `cfd-math::linear_solver::traits` no longer imports
  `cfd_core::conversion::SafeFromF64`, uses direct `T::from_f64(...)`, or uses
  fallback `T::from_f64_or(...)` conversion for convergence-factor exponent,
  CG theoretical-bound factor, or validation safety multiplier construction.
  Convergence-factor exponentiation now dispatches through Eunomia
  `FloatElement::powf`.
- **Boundary**: This slice preserves the existing nalgebra `DVector` and
  `LinearOperator` trait surfaces. Leto vector migration remains a later
  provider slice.
- **Evidence tier**: compile-time integration plus empirical focused
  value-semantic convergence tests and static source audit. Touched-file
  rustfmt passed. `cargo check -p cfd-math` passed. `cargo nextest run -p
  cfd-math convergence_factor cg_theoretical_bound validate_convergence`
  passed 4/4 tests, including the existing AMG convergence-factor match.
  Static scan found no direct `num_traits`, `FromPrimitive`, `ToPrimitive`,
  `Float::`, direct `T::from_f64`, `T::from_usize`, `SafeFromF64`,
  `from_f64_or`, conversion fallback, or stale `rayon` hits in
  `linear_solver/traits.rs`.
- **Residual risk**: Other `linear_solver` modules still contain direct
  provider residue, including the raw `operators/gpu.rs` path, multigrid, and
  tests. The core solver trait surfaces still use nalgebra vector types
  pending later Leto migration.

## Sprint 2026-07-03: cfd-math Eunomia Linear Operators
- **Resolved**: `cfd-math::linear_solver::operators::{poisson,momentum}` no
  longer use direct `num_traits::FromPrimitive` provider bounds or direct
  `T::from_f64(...).unwrap_or_else(num_traits::Zero::zero)` fallbacks for
  finite-difference constants. The touched operators now construct constants
  through Eunomia `FloatElement`.
- **Boundary**: This slice preserves the existing nalgebra `DVector` operator
  API. Leto vector migration and Hephaestus GPU-operator replacement remain
  later provider slices.
- **Evidence tier**: compile-time integration plus empirical focused
  value-semantic operator tests and static source audit. Touched-file rustfmt
  passed. `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math
  laplacian_center_impulse poisson_center_impulse` passed 2/2 tests. `cargo
  nextest run -p cfd-math momentum_1d_applies momentum_2d_applies
  energy_2d_applies` passed 3/3 tests. Static scan found no direct
  `num_traits`, `FromPrimitive`, `ToPrimitive`, `Float::`, direct
  `T::from_f64`, `T::from_usize`, conversion fallback, `from_f64_or`, or
  stale `rayon` hits in
  `linear_solver/operators/{poisson.rs,momentum.rs}`.
- **Residual risk**: Other `linear_solver` modules still contain direct
  provider residue, including the raw `operators/gpu.rs` path, multigrid,
  traits, and tests. The touched CPU operators still use nalgebra `DVector`
  surfaces pending later Leto migration.

## Sprint 2026-07-03: cfd-math Eunomia Schwarz/Cholesky Preconditioners
- **Resolved**: `cfd-math::linear_solver::preconditioners::schwarz` no longer
  imports or requires direct `num_traits::FromPrimitive`, and
  `cfd-math::linear_solver::preconditioners::cholesky` no longer uses direct
  `T::from_f64(...).unwrap_or_else(...)` for symmetry tolerance construction.
  Cholesky symmetry residual absolute-value dispatch is now explicit through
  Eunomia `NumericElement`.
- **Boundary**: This slice preserves existing nalgebra
  `CsrMatrix`/`DVector` preconditioner APIs. Leto sparse/vector migration
  remains a later provider slice.
- **Evidence tier**: compile-time integration plus empirical focused
  preconditioner tests and static source audit. Touched-file rustfmt passed.
  `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math cholesky`
  passed 4/4 tests. `cargo nextest run -p cfd-math schwarz` passed 1/1 test.
  Static scan found no direct `num_traits`, `FromPrimitive`, `ToPrimitive`,
  `Float::`, direct `T::from_f64`, `T::from_usize`, conversion fallback,
  `from_f64_or`, or stale `rayon` hits in
  `linear_solver/preconditioners/{schwarz.rs,cholesky.rs}`.
- **Residual risk**: Other `linear_solver` modules still contain direct
  `num-traits` usage, including operators, multigrid, traits, and tests.
  Schwarz still carries nalgebra sparse surfaces pending later Leto migration;
  Cholesky moved to Leto CSR in Sprint 1.96.166.

## Sprint 2026-07-03: cfd-math Eunomia SSOR Preconditioner
- **Resolved**: `cfd-math::linear_solver::preconditioners::ssor` no longer
  imports `num_traits::FromPrimitive` or uses direct `T::from_f64` conversion
  fallback for default relaxation construction.
- **Boundary**: This slice preserves the existing nalgebra
  `CsrMatrix`/`DVector` preconditioner API. Leto sparse/vector migration
  remains a later provider slice.
- **Evidence tier**: compile-time integration plus empirical focused SSOR
  tests and static source audit. Touched-file rustfmt passed. `cargo check -p
  cfd-math` passed. `cargo nextest run -p cfd-math ssor` passed 4/4 tests.
  Static scan found no direct `num_traits`, `FromPrimitive`, `ToPrimitive`,
  `Float::`, direct `T::from_f64`, `T::from_usize`, conversion fallback,
  `from_f64_or`, or stale `rayon` hits in
  `linear_solver/preconditioners/ssor.rs`.
- **Residual risk**: Other `linear_solver` modules still contain direct
  `num-traits` usage, including Schwarz preconditioners, operators,
  multigrid, traits, and tests. SSOR still carries nalgebra sparse/vector
  surfaces pending later Leto migration.

## Sprint 2026-07-03: cfd-math Eunomia Basic Preconditioners
- **Resolved**: `cfd-math::linear_solver::preconditioners::basic` no longer
  imports `num_traits::FromPrimitive` or uses direct `T::from_f64`,
  conversion-fallback, or direct nalgebra absolute-value dispatch for Jacobi
  diagonal tolerances, SOR omega bounds, or 1D-Poisson omega construction.
- **Boundary**: This slice preserves the existing nalgebra
  `CsrMatrix`/`DVector` preconditioner API. Leto sparse/vector migration
  remains a later provider slice.
- **Evidence tier**: compile-time integration plus empirical focused
  preconditioner tests and static source audit. Touched-file rustfmt passed.
  `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math jacobi`
  passed 5/5 tests. `cargo nextest run -p cfd-math sor` passed 6/6 tests.
  Static scan found no direct `num_traits`, `FromPrimitive`, `ToPrimitive`,
  `Float::`, direct `T::from_f64`, `T::from_usize`, conversion fallback,
  `from_f64_or`, or stale `rayon` hits in
  `linear_solver/preconditioners/basic.rs`.
- **Residual risk**: Other `linear_solver` modules still contain direct
  `num-traits` usage, including SSOR/Schwarz preconditioners, operators,
  multigrid, traits, and tests. Basic preconditioners still carry nalgebra
  sparse/vector surfaces pending later Leto migration.

## Sprint 2026-07-03: cfd-math Eunomia GMRES Chain
- **Resolved**: `cfd-math::linear_solver` GMRES, `LinearSolverChain`,
  `DirectSparseSolver`, `BlockDiagonalPreconditioner`, and
  `SimplePreconditioner` no longer import direct
  `num_traits::{Float, FromPrimitive, ToPrimitive}` or use direct
  `Float::abs`, `T::from_f64`, `T::from_usize`, or conversion fallback paths
  in the touched provider boundary.
- **Boundary**: This slice preserves nalgebra `DVector`/`DMatrix`/`RealField`
  APIs and the existing rsparse f64-backed direct sparse LU backend. The
  rsparse backend still converts through f64 by design until the direct solver
  is replaced by an Atlas-owned Leto/solver backend.
- **Evidence tier**: compile-time integration plus empirical focused
  linear-solver tests and static source audit. Touched-file rustfmt passed.
  `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math gmres`
  passed 21/21 tests. `cargo nextest run -p cfd-math direct_solver` passed
  3/3 tests. `cargo nextest run -p cfd-math block_preconditioner` passed 2/2
  tests. Static scan found no direct `num_traits`, `FromPrimitive`,
  `ToPrimitive`, `Float::`, direct `T::from_f64`, `T::from_usize`, conversion
  fallback, or stale `rayon` hits in `gmres`, `chain.rs`, `direct_solver.rs`,
  or `block_preconditioner.rs`.
- **Residual risk**: Other `linear_solver` modules still contain direct
  `num-traits` usage, including operators, basic/SSOR/Schwarz preconditioners,
  multigrid, and tests. The touched direct solver remains f64-backed through
  rsparse and must be replaced by an Atlas-owned backend in a later Leto/Apollo
  solver migration.

## Sprint 2026-07-03: cfd-math Eunomia Linear-Solver Config
- **Resolved**: `IterativeSolverConfig::default` no longer imports
  `num_traits::FromPrimitive` or uses direct `T::from_f64(...).expect(...)`
  construction for the default tolerance. CG and BiCGSTAB default
  constructors plus their linear-solver trait impl bounds now use Eunomia
  `FloatElement`; GMRES default construction declares the same provider bound.
  The config SpMV doc now names Moirai instead of Rayon.
- **Boundary**: This slice preserves the current linear-solver nalgebra
  `RealField`/`DVector` surface and does not migrate the rest of
  `linear_solver` away from direct `num-traits`.
- **Evidence tier**: compile-time integration plus empirical focused
  linear-solver tests and static source audit. Touched-file rustfmt passed.
  `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math
  default_solver` passed 2/2 tests. `cargo nextest run -p cfd-math
  test_gmres_configurable_trait` passed 1/1 test. Static scan found no direct
  `num_traits`, `FromPrimitive`, direct `T::from_f64`, or stale `rayon` hits
  in `linear_solver/config.rs`, `bicgstab/mod.rs`, or
  `conjugate_gradient/mod.rs`.
- **Residual risk**: `cfd-math::linear_solver` still contains direct
  `num-traits` usage in GMRES internals, direct solvers, operators,
  preconditioners, multigrid modules, and tests. Linear solvers still carry
  nalgebra vectors/matrices pending later Leto migration.

## Sprint 2026-07-03: cfd-math Eunomia SIMD Scalars
- **Resolved**: `cfd-math::simd` no longer imports
  `num_traits::FromPrimitive` or uses direct `T::from_f64`,
  conversion-fallback, or direct nalgebra square-root dispatch for the CFD SIMD
  central-difference constants and field norm.
- **Boundary**: This slice preserves the existing Moirai-backed parallel slice
  execution and nalgebra `RealField` API. Leto vector/scalar API cleanup
  remains a later provider slice.
- **Evidence tier**: compile-time integration plus empirical SIMD tests and
  static source audit. Touched-file rustfmt passed. `cargo check -p cfd-math`
  passed. `cargo nextest run -p cfd-math simd` passed 26/26 tests. Static scan
  found no direct `num_traits`, `FromPrimitive`, `ToPrimitive`, `Signed`,
  `Float::`, direct `T::from_f64`, conversion fallback, `T::epsilon`, or
  `rayon` hits in `crates/cfd-math/src/simd`.
- **Residual risk**: `cfd-math` still contains direct `num-traits` in
  linear-solver and test modules. SIMD still carries nalgebra `RealField` and
  nalgebra vector wrappers pending later Leto/eunomia scalar-bound cleanup.

## Sprint 2026-07-03: cfd-math Eunomia Sparse Scalars
- **Resolved**: `cfd-math::sparse` no longer imports
  `num_traits::{Float, FromPrimitive, Signed}` or uses direct `Float::...`,
  `T::from_f64`, conversion-fallback, `T::epsilon`, or direct signed absolute
  value paths for sparse stencil constants, Frobenius norms,
  condition-estimate singular-diagonal handling, or diagonal dominance checks.
- **Boundary**: This slice preserves the existing `nalgebra_sparse::CsrMatrix`
  sparse API and `DVector` SpMV surface. Leto sparse/vector migration remains
  a later provider slice.
- **Evidence tier**: compile-time integration plus empirical sparse tests and
  static source audit. Touched-file rustfmt passed. `cargo check -p cfd-math`
  passed. `cargo nextest run -p cfd-math sparse` passed 15/15 tests. Static
  scan found no direct `num_traits`, `FromPrimitive`, `ToPrimitive`, `Signed`,
  `Float::`, direct `T::from_f64`, conversion fallback, `T::epsilon`, or
  `rayon` hits in `crates/cfd-math/src/sparse`.
- **Residual risk**: `cfd-math` still contains direct `num-traits` in
  linear-solver, simd, and test modules. Sparse still carries
  `nalgebra_sparse::CsrMatrix` and nalgebra `DVector` pending later Leto
  sparse/vector migration.

## Sprint 2026-07-03: cfd-math Eunomia Nonlinear Solver Scalars
- **Resolved**: `cfd-math::nonlinear_solver` no longer imports
  `num_traits::{Float, FromPrimitive}` or uses direct `Float::...`,
  `T::from_f64`, conversion-fallback, or `T::epsilon` paths for
  Anderson/JFNK default constants, QR norms, diagonal checks,
  finite-difference perturbation safeguards, EW forcing clamps, Givens
  rotations, or back-substitution checks.
- **Boundary update**: The later Leto vector migration is now complete for
  `cfd-math::nonlinear_solver`; Anderson/JFNK no longer expose nalgebra
  `DVector`/`DMatrix` in this cone.
- **Evidence tier**: compile-time integration plus empirical nonlinear-solver
  tests and static source audit. Touched-file rustfmt passed. `cargo check -p
  cfd-math` passed. `cargo nextest run -p cfd-math nonlinear_solver` passed
  9/9 tests. Static scan found no direct `num_traits`, `FromPrimitive`,
  `Float::`, direct `T::from_f64`, conversion fallback, or `T::epsilon` hits
  in `crates/cfd-math/src/nonlinear_solver`.
- **Residual risk**: `cfd-math` still contains direct provider residue in
  linear-solver, sparse, simd, and test modules.

## Sprint 2026-07-03: cfd-math Eunomia Pressure-Velocity Scalars
- **Resolved**: `cfd-math::pressure_velocity` no longer imports
  nalgebra `RealField` or `num_traits::FromPrimitive`, and no longer uses
  direct `T::from_f64`/conversion `expect` construction or old
  `T::zero()`/`T::one()` identities for SIMPLE default tolerance, relaxation
  constants, or validation bounds.
- **Boundary**: This slice preserves the existing `SIMPLEConfig` and
  `SolveResult` APIs while changing their scalar contract to Eunomia
  `RealField`.
- **Evidence tier**: compile-time integration plus empirical pressure-velocity
  tests and static source audit. `cargo fmt -p cfd-math --check` passed.
  `cargo check -p cfd-math --no-default-features --lib` passed. `cargo
  nextest run -p cfd-math --no-default-features pressure_velocity
  --status-level fail` passed 3/3 tests with 341 skipped. `cargo clippy -p
  cfd-math --no-default-features --all-targets -- -D warnings` passed. Static
  scan found no nalgebra, `DVector`, `DMatrix`, old scalar identities,
  `num_traits`, `num_complex`, ndarray, rayon, tokio, rustfft,
  `FromPrimitive`, `ToPrimitive`, or `From<f64>` hits in
  `crates/cfd-math/src/pressure_velocity`.
- **Residual risk**: `cfd-math` still contains direct old-provider
  scalar/storage surfaces in linear-solver, sparse, SIMD, and test modules.
  The high-order and nonlinear-solver cones are closed by later July 4
  provider slices.

## Sprint 2026-07-03: cfd-math Eunomia Iterator Scalars
- **Resolved**: `cfd-math::iterators` no longer imports
  `num_traits::FromPrimitive` or uses direct `T::from_f64`, `T::from_usize`,
  or conversion-fallback scalar construction for stencil coefficients,
  iterator mean/variance counts, or standard-deviation square-root dispatch.
- **Boundary**: This slice preserves the iterator extension trait API and
  nalgebra `RealField` scalar surface. It also removes the prior
  second-derivative zero placeholder by returning real `[1, -2, 1]`
  coefficients for every declared 3-point stencil pattern.
- **Evidence tier**: compile-time integration plus empirical iterator tests
  and static source audit. Touched-file `rustfmt --check` passed. `cargo check
  -p cfd-math` passed. `cargo nextest run -p cfd-math iterators` passed 7/7
  tests. Static scan found no `num_traits`, `FromPrimitive`, `ToPrimitive`,
  `T::from_f64`, `T::from_usize`, `T::from_f32`, or conversion-fallback hits
  in `crates/cfd-math/src/iterators`.
- **Residual risk**: `cfd-math` still contains direct old-provider surfaces in
  linear-solver, sparse, SIMD, and test modules. The iterator cone itself no
  longer carries nalgebra `RealField` or old scalar identities; later July 4
  slices closed the high-order, interpolation, and nonlinear-solver cones.

## Sprint 2026-07-03: cfd-math Eunomia WENO Scalars
- **Resolved**: `cfd-math::high_order::weno` no longer imports
  `num_traits::FromPrimitive` or uses direct `T::from_f64` or
  conversion-fallback scalar construction for WENO5/WENO7 epsilon defaults,
  linear weights, ENO reconstruction coefficients, or smoothness-indicator
  constants.
- **Boundary**: This slice preserves the existing WENO5/WENO7 reconstruction
  formulas and nalgebra `RealField` scalar surface. Leto/eunomia
  scalar-bound cleanup remains separate.
- **Evidence tier**: compile-time integration plus empirical WENO tests and
  static source audit. Touched-file `rustfmt --check` passed. `cargo check -p
  cfd-math` passed. `cargo nextest run -p cfd-math weno` passed 6/6 tests.
  Static scan found no `num_traits`, `FromPrimitive`, `ToPrimitive`,
  `T::from_f64`, `T::from_usize`, `T::from_f32`, or conversion-fallback hits
  in `crates/cfd-math/src/high_order`.
- **Residual risk**: `cfd-math` still contains direct `num-traits` in
  linear-solver, sparse, simd, and test modules. The WENO scalar-bound residue
  from this older slice is closed by Sprint 2026-07-04 WENO Eunomia Scalars.

## Sprint 2026-07-03: cfd-math Eunomia Interpolation Scalars
- **Resolved**: `cfd-math::interpolation` no longer imports nalgebra
  `RealField` or `num_traits::FromPrimitive`, and no longer uses old
  `T::zero()`/`T::one()` identities, direct `T::from_f64`, or
  conversion-fallback scalar construction in the interpolation trait, linear,
  Lagrange, or cubic-spline implementations.
- **Boundary**: This slice preserves the existing interpolation trait APIs and
  current algorithms while changing their scalar contract to Eunomia
  `RealField`. Lagrange now rejects duplicate/non-increasing nodes before
  basis denominator division.
- **Evidence tier**: compile-time integration plus empirical interpolation
  tests and static source audit. `cargo fmt -p cfd-math --check` passed.
  `cargo check -p cfd-math --no-default-features --lib` passed. `cargo
  nextest run -p cfd-math --no-default-features interpolation --status-level
  fail` passed 15/15 tests with 328 skipped. `cargo clippy -p cfd-math
  --no-default-features --all-targets -- -D warnings` passed. Static scan
  found no nalgebra, `DVector`, `DMatrix`, old scalar identities,
  `num_traits`, `num_complex`, ndarray, rayon, tokio, rustfft,
  `FromPrimitive`, `ToPrimitive`, or `From<f64>` hits in
  `crates/cfd-math/src/interpolation`.
- **Residual risk**: `cfd-math` still contains direct old-provider
  scalar/storage surfaces in linear-solver, sparse, SIMD, and test modules.
  Later July 4 slices closed high-order and nonlinear-solver provider residue.

## Sprint 2026-07-03: cfd-math Eunomia Differentiation Scalars
- **Resolved**: `cfd-math::differentiation` finite-difference and gradient
  operator files no longer import `num_traits::FromPrimitive` or use direct
  `T::from_f64`, `T::from_f32`, `From<f64>`, or conversion-fallback scalar
  construction for stencil constants, gradient/divergence/curl constants, SIMD
  helper scalar staging, or 2D Laplacian constants.
- **Boundary**: Superseded by the 2026-07-04 differentiation follow-up, which
  replaced the preserved nalgebra `DVector`/`Vector3` surfaces with Leto
  `Array1`/`Vector3` while preserving the derivative formulas.
- **Evidence tier**: compile-time integration plus empirical differentiation
  tests and static source audit. Touched-file `rustfmt --check` passed.
  `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math
  differentiation` passed 12/12 tests. Static scan found no `num_traits`,
  `FromPrimitive`, `ToPrimitive`, `T::from_f64`, `T::from_usize`, `From<f64>`,
  `T::from_f32`, or conversion-fallback hits in
  `crates/cfd-math/src/differentiation/{finite_difference.rs,gradient.rs}`.
- **Residual risk**: Superseded for differentiation by the 2026-07-04 Leto
  follow-up. `cfd-math` still contains provider residue in linear-solver,
  sparse, and SIMD modules. Later July 4 slices closed high-order and
  nonlinear-solver provider residue.

## Sprint 2026-07-03: cfd-math Eunomia Integration Scalars
- **Resolved**: `cfd-math::integration` no longer imports
  `num_traits::FromPrimitive` or nalgebra `RealField`, and no longer uses
  direct `T::from_f64`, `T::from_usize`, `From<f64>`, `T::zero()`,
  `T::one()`, or conversion-fallback scalar construction for quadrature
  constants, interval counts, adaptive tolerance/error dispatch, and
  tetrahedral quadrature constants.
- **Boundary**: This slice preserves the existing `Quadrature` trait APIs and
  current quadrature algorithms while changing their scalar contract to
  Eunomia `RealField`.
- **Evidence tier**: compile-time integration plus empirical
  integration-named tests and static source audit. `cargo fmt -p cfd-math
  --check` passed. `cargo check -p cfd-math --no-default-features --lib`
  passed. `cargo nextest run -p cfd-math --no-default-features integration
  --status-level fail` passed 12/12 tests with 330 skipped. `cargo clippy -p
  cfd-math --no-default-features --all-targets -- -D warnings` passed. Static
  scan found no nalgebra, `DVector`, `DMatrix`, old scalar identities,
  `num_traits`, `num_complex`, ndarray, rayon, tokio, or rustfft hits in
  `crates/cfd-math/src/integration`.
- **Residual risk**: `cfd-math` still contains direct old-provider surfaces in
  linear-solver, sparse, SIMD, and test modules. Later July 4 slices closed the
  high-order, interpolation, and nonlinear-solver provider cones; remaining
  nalgebra storage replacement is outside this integration cone.

## Sprint 2026-07-03: cfd-2d Eunomia Scheme Complex
- **Resolved**: `cfd-2d::schemes::SpatialDiscretization::amplification_factor`
  no longer returns or constructs `num_complex::Complex<f64>`, and the crate
  no longer directly depends on `num-complex`.
- **Boundary**: This slice keeps existing finite-volume/finite-difference
  formulas and retains direct `num-traits` where current 2D APIs still require
  `FromPrimitive`/`ToPrimitive`. The new `FloatElement` bounds are explicit on
  paths that already construct provider-owned `cfd-core` grids or solver
  configs.
- **Evidence tier**: compile-time integration plus empirical regression tests
  and static source audit. `cargo fmt --package cfd-2d --package cfd-3d
  --package cfd-validation --check` passed. `cargo check -p cfd-2d`, `cargo
  check -p cfd-3d`, and `cargo check -p cfd-validation` passed. `cargo nextest
  run -p cfd-2d` passed 563/563 tests with 27 skipped. Static scan found no
  `num_complex`, `num-complex`, or `NumComplex` hits under `crates/cfd-2d`
  source or manifest files.
- **Residual risk**: Full cfd-2d provider migration still needs direct
  `num-traits` and nalgebra-owned matrix/vector paths removed. This slice
  clears the direct cfd-2d `num-complex` boundary and the dependency-chain
  compile blockers.

## Sprint 2026-07-03: cfd-3d Apollo/Leto Adapter Compile Repair
- **Resolved**: The private cfd-3d Atlas array adapter no longer imports
  nonexistent Apollo Leto-native FFT symbols or feature-gated
  `leto::MnemosyneStorage`, and cfd-3d/cfd-validation provider bounds no longer
  block `cargo nextest run -p cfd-2d`.
- **Boundary**: Algorithm storage stays Leto-owned. The adapter converts at the
  private Apollo boundary because the reachable Apollo FFT/NUFFT APIs still
  expose ndarray arrays.
- **Evidence tier**: compile-time integration plus empirical cfd-2d regression
  tests. `cargo check -p cfd-3d` passed. `cargo check -p cfd-validation`
  passed. `cargo nextest run -p cfd-2d` passed 563/563 tests with 27 skipped.
- **Residual risk**: Full cfd-3d ndarray elimination still depends on Apollo
  exposing reachable Leto-native FFT/NUFFT signatures or moving the conversion
  boundary upstream into Apollo.

## Sprint 2026-07-03: cfd-math/cfd-validation Eunomia Stability Complex
- **Resolved**: `cfd-math::time_stepping::stability` no longer exposes
  `num_complex::Complex<f64>` in the von Neumann spatial-operator callback API,
  and `cfd-validation::time_integration::stability_analysis` no longer imports
  or constructs `num_complex::Complex` for advection, diffusion, or Burgers
  stability operators.
- **Boundary**: This slice preserves the explicit RK stability polynomial. Since
  the analyzer already rejects non-explicit Butcher tableaux, `(I - zA)x = 1`
  now evaluates by forward substitution on the strictly lower-triangular `A`
  instead of allocating a dense `num_complex` matrix and inverting it through
  nalgebra.
- **Evidence tier**: compile-time integration plus empirical value-semantic
  stability tests and static source audit. Touched-file rustfmt passed.
  `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math stability`
  passed 5/5 tests. Static scan found no `num_complex`, `num-complex`, or
  `NumComplex` hits in `crates/cfd-math` or `crates/cfd-validation`.
- **Residual risk**: `num-complex` still resolves transitively through
  `nalgebra`/`simba` and provider crates. The direct cfd-2d scheme boundary and
  downstream cfd-validation compile blocker were cleared in the later cfd-2d
  Eunomia slice. Full removal still requires continuing the nalgebra/Leto
  provider migration.

## Sprint 2026-07-03: cfd-math Eunomia Stability Scalars
- **Resolved**: `cfd-math::time_stepping::stability` no longer imports
  `num_traits::ToPrimitive` or uses direct
  `T::from_f64(...).unwrap_or_else(num_traits::Zero::zero)` fallback
  conversions for analyzer defaults, CFL thresholds, RK order checks, or von
  Neumann amplification outputs.
- **Boundary**: This slice preserves the existing nalgebra `DMatrix`/`DVector`
  Butcher-tableau API. Leto replacement for that matrix/vector surface remains
  a separate migration item.
- **Evidence tier**: compile-time integration plus empirical value-semantic
  stability tests and static source audit. Touched-file `rustfmt --check`
  passed. `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math
  stability` passed 5/5 tests. `git diff --check` passed. Static scan found no
  `num_traits`, `ToPrimitive`, `FromPrimitive`, `T::from_f64`, or conversion
  fallback hits under `crates/cfd-math/src/time_stepping/stability`.
- **Residual risk**: Package-level `cargo fmt --package cfd-math --check` is
  still blocked by unrelated existing formatting drift outside the touched
  stability files. Broader cfd-math still contains direct `num-traits`,
  nalgebra matrix/vector surfaces, and GPU provider gaps.

## Sprint 2026-07-03: cfd-math Eunomia Runge-Kutta Scalars
- **Resolved**: `cfd-math::time_stepping::runge_kutta` no longer uses direct
  `T::from_f64(...).unwrap_or_else(num_traits::Zero::zero)` fallback
  conversions for RK3, RK4, or Carpenter-Kennedy low-storage RK4 constants.
- **Correctness fix**: `LowStorageRK4` now applies the Carpenter-Kennedy 2N
  residual recurrence `r_i = a_i r_{i-1} + dt f(t_i, u_i)`,
  `u_{i+1} = u_i + b_i r_i`. The prior implementation combined the solution
  accumulator with the stage state directly and failed the zero-RHS invariant.
- **Boundary**: This slice preserves the existing `DVector`-based
  `TimeStepper` API. Leto replacement for vector storage remains a separate
  migration item.
- **Evidence tier**: compile-time integration plus empirical value-semantic
  Runge-Kutta tests and static source audit. Touched-file `rustfmt --check`
  passed. `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math
  runge_kutta` passed 5/5 tests. Static scan found no `num_traits`,
  `FromPrimitive`, `ToPrimitive`, `T::from_f64`, or conversion-fallback hits in
  `crates/cfd-math/src/time_stepping/runge_kutta.rs`.
- **Residual risk**: `cfd-math` still contains direct `num-traits` in other
  time-stepping, high-order, and linear-solver modules. The
  `TimeStepper` API still uses nalgebra `DVector` pending the later Leto vector
  migration.

## Sprint 2026-07-03: cfd-math Eunomia Adaptive Scalars
- **Resolved**: `cfd-math::time_stepping::adaptive` no longer uses direct
  `T::from_f64(...).unwrap_or_else(num_traits::Zero::zero)` fallback
  conversions for adaptive-stepper defaults, PI controller gains, rejection
  scaling, clamp bounds, or Dormand-Prince tableau coefficients.
- **Boundary**: This slice preserves the existing `DVector`-based
  `TimeStepper`/`EmbeddedMethod` API. Leto replacement for vector storage
  remains a separate migration item.
- **Evidence tier**: compile-time integration plus empirical adaptive
  time-stepping tests and static source audit. Touched-file `rustfmt --check`
  passed. `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math
  adaptive` passed 3/3 tests. Static scan found no `num_traits`,
  `FromPrimitive`, `ToPrimitive`, `T::from_f64`, or conversion-fallback hits in
  `crates/cfd-math/src/time_stepping/adaptive.rs`.
- **Residual risk**: `cfd-math` still contains direct `num-traits` in other
  time-stepping, high-order, and linear-solver modules. The
  time-stepper state API still uses nalgebra `DVector` pending the later Leto
  vector migration.

## Sprint 2026-07-03: cfd-math Eunomia IMEX Scalars
- **Resolved**: `cfd-math::time_stepping::imex` no longer uses direct
  `T::from_f64(...).unwrap_or_else(num_traits::Zero::zero)` fallback
  conversions for Newton tolerance, ARS343 gamma/delta construction, tableau
  coefficients, or explicit/implicit solution weights.
- **Boundary**: This slice preserves the existing `DVector`/`DMatrix`-based
  `IMEXTimeStepper` API. Leto replacement for matrix/vector storage remains a
  separate migration item.
- **Evidence tier**: compile-time integration plus empirical IMEX tests and
  static source audit. Touched-file `rustfmt --check` passed. `cargo check -p
  cfd-math` passed. `cargo nextest run -p cfd-math imex` passed 5/5 tests.
  Static scan found no `num_traits`, `FromPrimitive`, `ToPrimitive`,
  `T::from_f64`, `T::from_usize`, or conversion-fallback hits in
  `crates/cfd-math/src/time_stepping/imex.rs`.
- **Residual risk**: `cfd-math` still contains direct `num-traits` in
  linear-solver, sparse, simd, and test modules.
  IMEX still uses nalgebra `DVector` and `DMatrix` pending the later
  Leto matrix/vector migration.

## Sprint 2026-07-03: cfd-math Eunomia Exponential Scalars
- **Resolved**: `cfd-math::time_stepping::exponential` no longer imports
  `num_traits::FromPrimitive` or uses direct
  `T::from_f64(...).unwrap_or_else(num_traits::Zero::zero)` fallback
  conversions for ERK4 stage weights, the RK combination denominator,
  phi-function small-argument thresholds, or scaling/squaring factorial
  conversions.
- **Boundary**: This slice preserves the existing `DVector`/`DMatrix`-based
  exponential integrator API. Leto replacement for matrix/vector storage
  remains a separate migration item.
- **Evidence tier**: compile-time integration plus empirical exponential
  time-stepping tests and static source audit. Touched-file `rustfmt --check`
  passed. `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math
  time_stepping::exponential::` passed 2/2 tests. Static scan found no
  `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_f64`, `T::from_usize`,
  or conversion-fallback hits in
  `crates/cfd-math/src/time_stepping/exponential.rs`.
- **Residual risk**: `cfd-math` still contains direct `num-traits` in
  linear-solver, sparse, simd, and test modules.
  Exponential time stepping still uses nalgebra `DVector` and
  `DMatrix` pending the later Leto matrix/vector migration.

## Sprint 2026-07-03: cfd-math Eunomia RKC Scalars
- **Resolved**: `cfd-math::time_stepping::rk_chebyshev` no longer imports
  `num_traits::FromPrimitive` or uses direct `T::from_f64`/`T::from_usize`
  conversions for RKC defaults, Chebyshev recurrence constants, adaptive
  step-halving, error normalization, clamp bounds, or error exponent dispatch.
- **Boundary**: This slice preserves the existing `DVector`-based
  `RhsFunction`/`RungeKuttaChebyshev` API. Leto replacement for vector storage
  remains a separate migration item.
- **Evidence tier**: compile-time integration plus empirical RKC tests and
  static source audit. Touched-file `rustfmt --check` passed. `cargo check -p
  cfd-math` passed. `cargo nextest run -p cfd-math rk_chebyshev` passed 4/4
  tests. Static scan found no `num_traits`, `FromPrimitive`, `ToPrimitive`,
  `T::from_f64`, `T::from_usize`, or conversion-fallback hits in
  `crates/cfd-math/src/time_stepping/rk_chebyshev.rs`.
- **Residual risk**: `cfd-math` still contains direct `num-traits` in
  linear-solver, sparse, simd, and test modules.
  The time-stepper state API still uses nalgebra `DVector` pending the
  later Leto vector migration.

## Sprint 2026-07-03: cfd-core Eunomia Rhie-Chow Interpolation
- **Resolved**: `cfd-core::physics::fluid_dynamics::rhie_chow` no longer imports `num_traits::FromPrimitive` or uses direct `T::from_f64(...).unwrap_or_else(...)` fallback conversions for the default Rhie-Chow relaxation factor or face-interpolation `2` constant.
- **Boundary**: This slice preserves the existing Rhie-Chow interpolation formulas and public type. The scalar conversion capability is now expressed through Eunomia `FloatElement`; broader fluid-dynamics finite-difference operations remain a separate direct `num-traits` migration slice.
- **Evidence tier**: compile-time integration plus empirical value-semantic Rhie-Chow tests and static source audit. `rustfmt --check --edition 2021 crates\cfd-core\src\physics\fluid_dynamics\rhie_chow.rs` passed. `cargo check -p cfd-core --no-default-features` passed. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core --no-default-features rhie_chow` passed 2/2 tests. `cargo nextest run -p cfd-core --no-default-features` passed 176/176 tests. Static scan found no `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_`, or conversion-fallback hits in `physics/fluid_dynamics/rhie_chow.rs`.
- **Residual risk**: Full cfd-core Eunomia migration still requires remaining direct numeric conversion modules outside this Rhie-Chow slice, including `physics/fluid_dynamics/operations.rs`, `physics/cavitation/number.rs`, and the shared `management::conversion` facade used by downstream crates. Full Leto CPU replacement for nalgebra-owned fluid-dynamics vectors and matrices remains a separate migration.

---

## Sprint 2026-07-03: cfd-core Eunomia Boundary Geometry
- **Resolved**: `cfd-core::physics::boundary::geometry` no longer uses direct `T::from_f64(...).unwrap_or_else(...)` fallback conversions for sphere and cylinder measure constants.
- **Boundary**: This slice preserves the existing public geometry types and formulas. Only `measure()` gains the Eunomia `FloatElement` bound because it constructs scalar constants; `contains_point()` and `dimension()` remain available under the narrower existing `RealField + Copy` bound.
- **Evidence tier**: compile-time integration plus empirical value-semantic geometry tests and static source audit. `rustfmt --check --edition 2021 crates\cfd-core\src\physics\boundary\geometry.rs` passed. `cargo check -p cfd-core --no-default-features` passed. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core --no-default-features boundary::geometry` passed 4/4 tests. `cargo nextest run -p cfd-core --no-default-features` passed 174/174 tests. Static scan found no `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_`, or conversion-fallback hits in `physics/boundary/geometry.rs`.
- **Residual risk**: Full cfd-core Eunomia migration still requires remaining direct numeric conversion modules outside this boundary geometry slice, including `physics/fluid_dynamics/rhie_chow.rs`, `physics/fluid_dynamics/operations.rs`, `physics/cavitation/number.rs`, and the shared `management::conversion` facade used by downstream crates. Full Leto CPU replacement for nalgebra-owned boundary/geometry surfaces remains a separate migration.

---

## Sprint 2026-07-03: cfd-core Eunomia Boundary Ghost Cells
- **Resolved**: `cfd-core::physics::boundary::ghost_cells` no longer imports `num_traits::{FromPrimitive, ToPrimitive}` or uses direct `T::from_f64(...).unwrap_or_else(...)` conversions for ghost-cell coefficients or Robin singularity tolerances.
- **Boundary**: This slice preserves the existing Dirichlet, Neumann, and Robin ghost-cell formulas while making scalar constants provider-owned. It also makes the documented degenerate Robin coefficient case explicit by returning `BoundaryErrorKind::RobinSingularity` before dividing by zero when both `alpha` and `beta` are zero.
- **Evidence tier**: compile-time integration plus empirical value-semantic ghost-cell tests and static source audit. `cargo check -p cfd-core --no-default-features` passed. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core --no-default-features ghost_cells` passed 5/5 tests. `cargo nextest run -p cfd-core --no-default-features` passed 170/170 tests. Touched-file `rustfmt --check` passed. Static scan found no `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_`, or conversion-fallback hits in `physics/boundary/ghost_cells.rs`.
- **Residual risk**: Full cfd-core Eunomia migration still requires remaining direct numeric conversion modules outside this boundary ghost-cell slice, including boundary geometry helpers, mesh operations, cavitation, fluid dynamics, and blood/non-Newtonian fluid modules. Full Leto CPU replacement for nalgebra-owned boundary/geometry surfaces remains a separate migration.

---

## Sprint 2026-07-03: cfd-core Eunomia Staggered Grid
- **Resolved**: `cfd-core::geometry::staggered` no longer imports `num_traits::FromPrimitive` or uses direct `T::from_usize`/`T::from_f64` conversions for staggered-grid dimensions, coordinate indices, or half-cell offsets.
- **Boundary**: This slice preserves existing uniform and stretched staggered-grid formulas. Integer dimensions and indices now assert exact representability before conversion so the migration does not silently round grid topology through the scalar provider.
- **Evidence tier**: compile-time integration plus empirical value-semantic staggered-grid tests and static source audit. `cargo check -p cfd-core --no-default-features` passed. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core --no-default-features staggered` passed 5/5 tests. `cargo nextest run -p cfd-core --no-default-features` passed 167/167 tests. Touched-file `rustfmt --check` and `git diff --check` passed. Static scan found no `num_traits`, `FromPrimitive`, `T::from_`, `from_usize`, or conversion-fallback hits in `geometry/staggered.rs`.
- **Residual risk**: Full cfd-core Eunomia migration still requires remaining direct numeric conversion modules outside this geometry slice, including boundary geometry/ghost-cell helpers, mesh operations, cavitation, fluid dynamics, and blood/non-Newtonian fluid modules. Full Leto CPU replacement for nalgebra-owned geometry surfaces remains a separate migration.

---

## Sprint 2026-07-03: cfd-core Eunomia Boundary Time Functions
- **Resolved**: `cfd-core::physics::boundary::{time_dependent,applicator}` no longer use direct `num_traits::FromPrimitive` conversions or silent `T::one()` fallbacks for `2π`, ghost-cell reflection constants, or time-function exponential/trigonometric dispatch.
- **Boundary**: This slice preserves existing boundary condition APIs and formulas while making the scalar math provider explicit. The `FloatElement` bound now propagates through boundary applicator, concrete applicator, specification, and manager surfaces that evaluate time-dependent boundary conditions or ghost-cell formulas.
- **Evidence tier**: compile-time integration plus empirical value-semantic boundary tests and static source audit. `cargo check -p cfd-core --no-default-features` passed. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core --no-default-features boundary` passed 30/30 tests. `cargo nextest run -p cfd-core --no-default-features` passed 167/167 tests. Touched-file `rustfmt --check` passed. Static scan found no `num_traits`, `FromPrimitive`, `T::from_f64`, conversion-fallback, bare `.exp()`, or bare `.sin()` hits in the touched boundary files.
- **Residual risk**: Full cfd-core Eunomia migration still requires remaining direct numeric conversion modules outside this boundary slice, including boundary ghost-cell geometry helpers, geometry operations/staggered-grid, fluid dynamics, cavitation, and blood/non-Newtonian fluid modules. Full Leto CPU replacement for nalgebra-owned boundary/geometry surfaces remains a separate migration.

---

## Sprint 2026-07-03: cfd-core Eunomia Temperature-Dependent Fluids
- **Resolved**: `cfd-core::physics::fluid::temperature` no longer imports `num_traits::FromPrimitive`, no longer carries stale conversion bounds on polynomial fluid models, and no longer uses direct scalar conversion fallbacks for the Sutherland exponent.
- **Boundary**: This slice preserves the existing temperature-dependent fluid APIs and formulas. Arrhenius, Andrade, and Sutherland math dispatch now uses Eunomia `FloatElement` explicitly so the touched module no longer depends on nalgebra method resolution for `exp`/`powf`.
- **Evidence tier**: compile-time integration plus empirical value-semantic temperature-fluid tests and static source audit. `cargo check -p cfd-core --no-default-features` passed. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core --no-default-features fluid::temperature` passed 6/6 tests. `cargo nextest run -p cfd-core --no-default-features` passed 162/162 tests. Touched-file `rustfmt --check` passed. Static scan found no `num_traits`, `FromPrimitive`, `T::from_f64`, conversion-fallback, bare `.exp()`, or bare `.powf(...)` hits in `physics/fluid/temperature.rs`.
- **Residual risk**: Full cfd-core Eunomia migration still requires remaining direct numeric conversion modules outside temperature-dependent fluids, including geometry operations/staggered-grid, fluid dynamics, cavitation, boundary, and blood/non-Newtonian fluid modules. Full Leto CPU replacement for nalgebra-owned physics surfaces remains a separate migration.

---

## Sprint 2026-07-03: cfd-core Eunomia Hemolysis Constants
- **Resolved**: `cfd-core::physics::hemolysis::{calculator,trauma}` no longer imports `num_traits::FromPrimitive` or uses direct scalar conversion fallbacks for NIH/MIH/exposure-time constants or platelet activation defaults.
- **Boundary**: This slice preserves the existing hemolysis calculator APIs and clinical-index formulas. It also fixes `PlateletActivation::activation_probability` to compute `1 - exp(-k * excess_stress * exposure_time)`; the prior implementation evaluated `1 - exp(k * excess_stress * exposure_time)`, producing negative probabilities above threshold.
- **Evidence tier**: compile-time integration plus empirical value-semantic hemolysis tests and static source audit. `cargo check -p cfd-core --no-default-features` passed. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core --no-default-features hemolysis` passed 9/9 tests. `cargo nextest run -p cfd-core --no-default-features` passed 160/160 tests. Touched-file `rustfmt --check` passed. Static scan found no `num_traits`, `FromPrimitive`, `T::from_f64`, or conversion-fallback hits in `physics/hemolysis/{calculator,trauma}.rs`.
- **Residual risk**: Full cfd-core Eunomia migration still requires remaining direct numeric conversion modules outside hemolysis, including geometry operations/staggered-grid, fluid dynamics, cavitation, boundary, and blood/non-Newtonian fluid modules. Full Leto CPU replacement for nalgebra-owned physics surfaces remains a separate migration.

---

## Sprint 2026-07-03: cfd-core Eunomia Mesh Quality Thresholds
- **Resolved**: `cfd-core::geometry::mesh::quality` no longer imports `num_traits::FromPrimitive` or uses direct scalar conversion fallbacks for aspect-ratio, skewness, or orthogonality quality thresholds.
- **Boundary**: This slice preserves the existing quality-level strict comparison contract and recommendation text. `nalgebra::RealField` remains because mesh quality statistics still use the current geometry scalar boundary.
- **Evidence tier**: compile-time integration plus empirical value-semantic mesh quality tests and static source audit. `cargo check -p cfd-core --no-default-features` passed. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core --no-default-features mesh_quality` passed 3/3 tests. `cargo nextest run -p cfd-core --no-default-features` passed 158/158 tests. Touched-file `rustfmt --check` passed. Static scan found no `num_traits`, `FromPrimitive`, `T::from_f64`, or conversion-fallback hits in `geometry/mesh/quality.rs`.
- **Residual risk**: Full cfd-core Eunomia migration still requires remaining direct numeric conversion modules outside mesh quality, including geometry operations/staggered-grid and cavitation/hemolysis modules. Full Leto CPU replacement for nalgebra-owned geometry surfaces remains a separate migration.

---

## Sprint 2026-07-03: cfd-core Eunomia CPU Backend
- **Resolved**: `cfd-core::compute::cpu` no longer imports `num_traits::FromPrimitive`, no longer uses the local `safe_f64_to_t` fallback conversion helper, and no longer requires scalar-conversion bounds for `CpuBuffer` storage operations.
- **Boundary**: This slice preserves the existing CPU advection discretization and `KernelParams` API. The advection kernel still accepts `f64` domain parameters through the existing trait contract, but conversion into the generic scalar type is now owned by Eunomia `FloatElement` rather than direct `num-traits`.
- **Evidence tier**: compile-time integration plus empirical value-semantic backend tests and static source audit. `cargo check -p cfd-core --no-default-features` passed. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core --no-default-features test_cpu_advection_kernel_linear_exactness` passed 1/1 test. `cargo nextest run -p cfd-core --no-default-features` passed 155/155 tests. Touched-file `rustfmt --check` passed. Static scan found no `num_traits`, `FromPrimitive`, `safe_f64_to_t`, or conversion-fallback hits in `compute/cpu.rs`.
- **Residual risk**: Full cfd-core Eunomia migration still requires remaining direct numeric conversion modules outside `compute/cpu`. Full Leto CPU replacement for nalgebra-owned vector/matrix surfaces remains a separate migration.

---

## Sprint 2026-07-03: cfd-core Eunomia Time Integrators
- **Resolved**: `cfd-core::compute::time::integrators` no longer uses direct `num_traits::FromPrimitive` for RK2/RK4 scalar constants, Crank-Nicolson half-step constants, or implicit integrator default tolerances.
- **Boundary**: This slice preserves the existing time integrator APIs and formulas. `ForwardEuler` remains conversion-free. RK2/RK4/Crank-Nicolson now require Eunomia `FloatElement` only where scalar constants are needed. `BackwardEuler` and `CrankNicolson` defaults now use exact Eunomia constants instead of silently substituting `T::one()` if conversion failed.
- **Evidence tier**: compile-time integration plus empirical value-semantic tests and static source audit. `cargo check -p cfd-core --no-default-features` passed. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core --no-default-features compute::time::integrators` passed 6/6 tests. `cargo nextest run -p cfd-core --no-default-features` passed 155/155 tests. Touched-file `rustfmt --check` passed. Static scan found no `num_traits`, `FromPrimitive`, `T::from_f64`, or conversion-fallback hits in `compute/time/integrators.rs`.
- **Residual risk**: Full cfd-core Eunomia migration still requires remaining direct numeric conversion modules outside `compute/time`. Full nalgebra replacement remains a separate Leto-backed state/vector migration.

---

## Sprint 2026-07-03: cfd-core Eunomia Time-Step Controllers
- **Resolved**: `cfd-core::compute::time::controllers` no longer uses direct `num_traits::FromPrimitive` or `num_traits::Float` for controller defaults, power functions, min/max clamping, or integration-order conversion.
- **Boundary**: This slice preserves the physical controller formulas and default values, but changes `calculate_dt` on both controllers to return `Result<T>` so zero, overflowing, or unsupported integration orders are surfaced as typed configuration errors instead of silently falling back to order one. There were no in-repo call sites to update.
- **Evidence tier**: compile-time integration plus empirical value-semantic tests and static source audit. `cargo check -p cfd-core --no-default-features` passed. `cargo nextest run -p cfd-core --no-default-features compute::time::controllers` passed 4/4 tests. `cargo nextest run -p cfd-core --no-default-features` passed 149/149 tests. Touched-file `rustfmt --check` passed. Static scan found no `num_traits`, `FromPrimitive`, `T::from_f64`, or conversion-fallback hits in `compute/time/controllers.rs`.
- **Residual risk**: Full cfd-core Eunomia migration still requires remaining direct numeric conversion modules. Time integrators still use direct `num-traits` and remain a separate slice because they return existing `Result` values for conversion failures.

---

## Sprint 2026-07-03: cfd-core Eunomia Solver Config Defaults
- **Resolved**: `cfd-core::compute::solver::config` no longer uses direct `num_traits::FromPrimitive` conversions or silent zero fallbacks for default convergence tolerances, relative tolerances, time step, CFL number, or linear solver tolerance.
- **Boundary**: This slice preserves the existing default values and builder API while routing builder construction through the canonical `Default` implementation. It adds value-semantic tests for the changed defaults. `cfd-core` still retains direct `num-traits` in modules outside this cone.
- **Evidence tier**: compile-time integration plus empirical value-semantic tests and static source audit. `cargo check -p cfd-core --no-default-features` passed. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core --no-default-features solver_config` passed 2/2 tests. `cargo nextest run -p cfd-core --no-default-features` passed 146/146 tests. Touched-file `rustfmt --check` passed. Static scan found no `num_traits`, `FromPrimitive`, `T::from_f64`, or zero-conversion fallback hits in `compute/solver/config.rs`.
- **Residual risk**: Full cfd-core Eunomia migration still requires remaining direct numeric conversion modules. Solver runtime topology still uses existing nalgebra scalar boundaries pending later Leto-backed replacement.

---

## Sprint 2026-07-03: cfd-core Eunomia Abstraction Defaults
- **Resolved**: `cfd-core::abstractions::{state,problem}` no longer use direct `num_traits::FromPrimitive` conversions for default time-step, reference-pressure, or reference-temperature constants.
- **Boundary**: This slice preserves `FieldState::new`, `ProblemParameters::default`, and `ProblemBuilder::new` behavior while narrowing field-state accessor/mutator methods back to `RealField + Copy`. `cfd-core` still retains direct `num-traits` in other modules outside this cone, including solver config, boundary helpers, fluid dynamics, cavitation, hemolysis, and blood/non-Newtonian models.
- **Evidence tier**: compile-time integration plus empirical abstraction tests and static source audit. `cargo check -p cfd-core --no-default-features` passed. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core --no-default-features abstractions` passed 2/2 tests. `cargo nextest run -p cfd-core --no-default-features` passed 144/144 tests. Touched-file `rustfmt --check` passed. Static scan found no `num_traits`, `FromPrimitive`, `T::from_f64`, or zero/one conversion-fallback hits in `abstractions/{state,problem}.rs`.
- **Residual risk**: Full cfd-core Eunomia migration still requires the remaining direct numeric conversion modules. Full nalgebra replacement remains a later Leto-backed state/vector migration.

---

## Sprint 2026-07-03: cfd-core Eunomia Fluid Validation Thresholds
- **Resolved**: `cfd-core::physics::fluid::validation` no longer uses direct `num_traits::FromPrimitive` conversions or silent zero/one fallbacks for property bounds, Reynolds/Prandtl limits, temperature limits, or pressure limits.
- **Boundary**: This slice preserves the existing validation thresholds and error behavior. `cfd-core` still retains direct `num-traits` in other modules outside this cone, including solver config, fluid dynamics, cavitation, hemolysis, and blood models.
- **Evidence tier**: compile-time integration plus empirical validation tests and static source audit. `cargo check -p cfd-core --no-default-features` passed. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core --no-default-features physics::fluid::validation` passed 2/2 tests. `cargo nextest run -p cfd-core --no-default-features` passed 144/144 tests. Touched-file `rustfmt --check` and `git diff --check` passed. Static scan found no `num_traits`, `FromPrimitive`, `T::from_f64`, or conversion-fallback hits in `validation.rs`.
- **Residual risk**: Full cfd-core Eunomia migration still requires the remaining direct numeric conversion modules. Full nalgebra replacement remains a later Leto-backed state/vector migration.

---

## Sprint 2026-07-03: cfd-core Eunomia Material/Fluid Constants
- **Resolved**: `cfd-core` material constructors (`ElasticSolid::{steel, aluminum}`, `FluidSolidInterface::water_air`, `MaterialDatabase::with_common_materials`) and immediate fluid constructor dependencies (`ConstantPropertyFluid::{water_20c, air_20c}`, fluid database constructors, `IdealGas`) no longer use direct `num_traits::FromPrimitive` scalar conversions or silent zero/one conversion fallbacks.
- **Boundary**: `nalgebra::RealField` remains because these modules still implement existing cfd-core physics traits over nalgebra scalar bounds. `cfd-core` still retains direct `num-traits` for other modules outside this cone, including blood/non-Newtonian, hemolysis, fluid dynamics, and management conversion helpers.
- **Evidence tier**: compile-time integration plus empirical regression tests and static source audit. `cargo check -p cfd-core --no-default-features` passed. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core --no-default-features` passed 144/144 tests. Touched-file `rustfmt --check` and touched-file `git diff --check` passed. Static scan over the touched cone found no `num_traits`, `FromPrimitive`, `T::from_f64`, or conversion-fallback hits.
- **Residual risk**: Full cfd-core Eunomia migration still needs the remaining direct numeric conversion modules. Full nalgebra removal still requires later Leto-backed vector/state replacements rather than scalar-constant rewrites.

---

## Sprint 2026-07-03: cfd-core Eunomia Physics Value Boundary
- **Resolved**: `cfd-core` physics value objects (`Velocity`, `Temperature`, `Pressure`, `ReynoldsNumber`, `DimensionlessNumber`) and their dependent management aggregates no longer use direct `num_traits::FromPrimitive` bounds or `T::from_f64(...).unwrap_or_else(T::zero/one)` scalar constant conversions.
- **Boundary**: `nalgebra::RealField` remains in this slice because the touched value objects still store `nalgebra::Vector3` values and use nalgebra vector/scalar methods. The cfd-core manifest still retains `num-traits` because `management::conversion` and other modules outside the touched cone still use direct `num-traits`.
- **Evidence tier**: compile-time integration plus empirical value-semantic regression tests and static source audit. `cargo check -p cfd-core --no-default-features` passed. `cargo nextest run -p cfd-core --no-default-features` passed 144/144 tests. Touched-file `rustfmt --check` passed. Static scan over the touched cone found no `num_traits`, `FromPrimitive`, `T::from_f64`, or silent fallback conversion hits.
- **Residual risk**: Full cfd-core Eunomia migration still needs `management::conversion`, remaining numeric helper modules, and eventual Leto replacement for nalgebra-owned vector/matrix state. Package-level `cargo fmt --check --package cfd-core` is still blocked by pre-existing unrelated formatting drift in `error.rs`, `physics/boundary/ghost_cells.rs`, and `physics/fluid_dynamics/operations.rs`.

---

## Sprint 2026-07-03: cfd-python Eunomia Blood Binding Boundary
- **Resolved**: `cfd-python` no longer directly depends on `num-traits`; blood-model PyO3 wrappers use Eunomia `NumericElement` for exposed field conversions and pass `f64` shear-rate inputs directly into the Rust-owned blood models.
- **Boundary**: NumPy remains the external Python ABI surface and cfd-python still resolves provider-transitive numeric crates through its internal CFD dependencies. The cfd-python source/manifest boundary no longer owns `num-traits`, `FromPrimitive`, or `ToPrimitive`. Getter methods now expose the actual model fields instead of silently falling back to hardcoded defaults.
- **Evidence tier**: compile-time integration plus static source/dependency audit. `cargo check -p cfd-python` passed; `cargo nextest run -p cfd-python --no-tests pass` confirmed 0 cfd-python test binaries; `cargo fmt -p cfd-python --check` passed; static scan found no `num_traits`, `num-traits`, `FromPrimitive`, or `ToPrimitive` hits under `crates/cfd-python`. Dependency-chain clippy blockers found during the cfd-python clippy gate were mechanically fixed across cfd-1d, cfd-2d, cfd-3d, cfd-validation, and cfd-python; `cargo clippy -p cfd-python --all-targets -- -D warnings` now passes.
- **Residual risk**: Other CFDrs crates still use direct `num-traits`/`num-complex`; those remain for later Eunomia/provider migration slices. Package-wide cfd-2d formatting still has pre-existing unrelated drift outside the touched files.

---

## Sprint 2026-07-03: cfd-io Eunomia Scalar Boundary
- **Resolved**: `cfd-io` no longer directly depends on `num-traits`; checkpoint validation, binary vector/matrix helpers, and CSV helper types now use Eunomia `RealField`.
- **Boundary**: The workspace still has `num-traits` for other crates, and cfd-io still resolves it transitively through upstream provider crates (`leto`/`eunomia`/`half`/`ndarray`). The cfd-io source and manifest boundary is now Eunomia-owned. Checkpoint mass-conservation validation now rejects zero or non-exactly-convertible mesh dimensions instead of silently falling back to unit spacing.
- **Evidence tier**: compile-time integration plus empirical value-semantic regression tests. `cargo check -p cfd-io` passed; `cargo check -p cfd-io --all-features` passed; `cargo nextest run -p cfd-io --no-fail-fast` passed 3/3 checkpoint roundtrip/property tests; `cargo clippy -p cfd-io --all-targets --all-features -- -D warnings` passed; `cargo fmt -p cfd-io --check` passed; static scan found no `num_traits`, `num-traits`, `FromPrimitive`, or `ToPrimitive` hits under `crates/cfd-io`.
- **Residual risk**: Other CFDrs crates still use direct `num-traits` and `num-complex`; cfd-io still has provider-transitive `num-traits` through Leto/Eunomia. Those remain for later Eunomia/provider migration slices.

---

## Sprint 2026-07-03: cfd-io Leto Checkpoint/Binary Boundary
- **Resolved**: `cfd-io` no longer directly depends on or normally resolves `nalgebra`; checkpoint fields and binary vector/matrix helpers now use Leto `Array1`/`Array2`, and CSV scalar bounds no longer use `nalgebra::RealField`.
- **Boundary**: `cfd-io` now owns a local file-format `Error`/`Result` type so the I/O crate does not pull `cfd-core` and its remaining nalgebra graph just to report I/O failures. This is a breaking API alignment for callers that referenced `cfd_core::error::Error` through cfd-io return types.
- **Evidence tier**: compile-time integration plus empirical value-semantic regression tests. `cargo check -p cfd-io` passed; `cargo nextest run -p cfd-io --no-fail-fast` passed 3/3 checkpoint roundtrip/property tests; `cargo check -p cfd-io --all-features` passed; `cargo fmt -p cfd-io --check` passed; `cargo tree -p cfd-io -e normal -i nalgebra` reported no matching package; static scan found no `nalgebra`, `DMatrix`, `DVector`, or `RealField` hits under `crates/cfd-io`.
- **Residual risk**: Optional `vtk` feature dependencies may still resolve nalgebra through RITK/Gaia/Burn; the default cfd-io graph no longer does.

---

## Sprint 2026-07-03: cfd-python Leto PyO3 Boundary
- **Resolved**: `cfd-python` no longer directly depends on `ndarray` or `nalgebra`; 2D NumPy-return helper paths construct dense Rust-owned data with Leto `Array2` and copy into NumPy only at the PyO3 boundary.
- **Boundary**: NumPy remains the required Python-visible array representation. Leto currently depends transitively on `ndarray` through its compatibility feature in the shared workspace dependency, but cfd-python no longer owns ndarray/nalgebra arrays directly.
- **Evidence tier**: compile-time integration plus static source/dependency audit: `cargo check -p cfd-python` passed; `cargo nextest run -p cfd-python --no-tests pass` confirmed there are 0 cfd-python test binaries; `cargo fmt -p cfd-python --check` passed; static scan found no `ndarray`, `nalgebra`, or `DMatrix` source/manifest hits in `crates/cfd-python`; lockfile now lists `leto` instead of direct `nalgebra`/`ndarray` under `cfd-python`.
- **Residual risk**: cfd-python still uses `num-traits` directly in blood-model bindings; that remains for the later Eunomia numeric-trait migration.

---

## Sprint 2026-07-03: Atlas GPU/Concurrency Provider Migration
- **Resolved**: CFDrs no longer resolves `pollster`; `cfd-core` GPU context creation, GPU capability probing, and Poisson residual readback, plus `cfd-math` GPU operator sync dispatch, now synchronize through `moirai::block_on`.
- **Boundary**: Hephaestus cannot directly replace `cfd-core::compute::gpu::GpuContext` in the same slice because CFDrs depends on `wgpu 0.19` and Hephaestus' `hephaestus-wgpu` provider depends on `wgpu 26.0`. The raw `wgpu::Device`, `wgpu::Queue`, `wgpu::Buffer`, and pipeline types are not ABI/API-compatible across those major versions.
- **Evidence tier**: compile-time integration plus empirical regression tests: `cargo check -p cfd-core --features gpu` passed; `cargo nextest run -p cfd-core --features gpu --lib` passed 151/151 tests; `cargo check -p cfd-math --features gpu` passed; `cargo nextest run -p cfd-math --features gpu --lib` passed 279/279 tests; static scan found no `pollster` usage in manifests, lockfile, or Rust sources.
- **Residual risk**: Full GPU-provider migration requires a coordinated CFDrs WGPU upgrade or a higher-level rewrite of GPU kernels onto Hephaestus buffer/kernel abstractions before raw context ownership can move to `hephaestus-wgpu`.

---

## Sprint 2026-07-02: Atlas Array Provider Migration
- **Resolved**: `cfd-3d` spectral DNS, forcing, diagnostics, Fourier wrapper, and IBM NUFFT code now own dense 1D/3D arrays through Leto instead of direct `ndarray` algorithm storage.
- **Boundary**: The reachable Apollo FFT/NUFFT APIs still accept and return `ndarray`, so `crates/cfd-3d/src/atlas_array.rs` contains the only conversion boundary for this slice.
- **Evidence tier**: compile-time integration plus empirical value-semantic regression tests: `cargo check -p cfd-3d --no-default-features` passed; `cargo nextest run -p cfd-3d --no-default-features --lib` passed 206/206 tests.
- **Residual risk**: Full ndarray removal from `cfd-3d` depends on publishing or consuming Apollo Leto-native FFT/NUFFT APIs; Apollo ndarray conversion remains boundary-only in this slice.

---

# Gap Audit: cfd-schematics

## Phase 1: Foundation (Current)
- Auditing `cfd-schematics`...
- Goal: Identify redundancy, shared component duplication, and structural integrity violations.
- Checking for placeholders, temporary workarounds, and approximations.

## Findings
### 1. Dual-Path Topology Generation (SSOT & DRY Violation)
- `interface::presets` (e.g., `bifurcation.rs`, `serpentine.rs`) directly instantiate `NetworkBlueprint`, manually inserting nodes and channels.
- `topology::presets` combined with `topology::factory::BlueprintTopologyFactory` does exactly the same thing but via declarative `BlueprintTopologySpec`.
- **Gap**: `interface::presets` must be completely removed or refactored into thin wrappers that call `topology::presets` and `BlueprintTopologyFactory`.

### 2. Leaking Physics Constants (SOC Violation)
- `BLOOD_MU` (3.5e-3) and `shah_london_resistance` / `hp_resistance` formulas are hardcoded into schematic generation (`cfd-schematics/src/interface/presets/*` and `factory.rs`).
- **Gap**: Schematic generation (topology/geometry mappings) should not calculate fluidic resistances. That is the domain of `cfd-1d`. The schematic layer should only provide geometry (`length_m`, `width_m`, `height_m`) and let downstream solvers compute resistances based on their own runtime fluid models.

---

# Finding 2026-07-10: cfd-1d double-trifurcation Picard non-convergence (test regression)

- **Symptom**: `cfd-suite` integration test `tests/cross_fidelity_blueprint.rs::cross_fidelity_blueprint_complex_branching` fails. `solve_reference_trace` → cfd-1d `Network2DSolver` returns `MaxIterationsExceeded: Maximum iterations (10000) exceeded` on the `double_trifurcation_cif_venturi_rect` network (blood ρ=1060, µ=3.5e-3, Q=1e-7 m³/s). The sibling `cross_fidelity_blueprint_bifurcation` passes.
- **Isolation**: Reproduces with `-p cfd-suite --no-default-features` (gpu OFF), so it is independent of the concurrent cfd-core `gpu/turbulence_compute` directory-module refactor that currently breaks the `gpu` feature build (`E0583 turbulence_compute` — peer WIP, not this finding).
- **Locus**: `crates/cfd-1d/src/solver/core/mod.rs` Picard fixed-point loop (line ~358) with Anderson acceleration + Newton fallback; `has_converged_dual` never satisfied within 10000 iters for this stiff nonlinear network. Solver correctly returns a typed error (no silent bad result).
- **Hypothesis (unverified)**: migration regression from the Leto-CSR assembly push (`d58d1fe3`/`1d768895`) subtly altering the assembled matrix, or a genuine conditioning/stiffness limit the current Picard+Anderson path cannot crack (cf. recent peer work `0d101352` "enhance Anderson QR collapse detection" — actively-evolving convergence path).
- **DoR to resolve**: (1) differential-test the assembled cfd-1d matrix + rhs for this network against the pre-migration commit (parent of `d58d1fe3`) to classify regression vs. genuine stiffness; (2) capture the Picard residual trajectory (diverge / plateau / slow-decay) via a bounded-iteration diagnostic; (3) if genuine stiffness, verify the Newton fallback engages on Picard stagnation. Blocked from immediate fix: convergence path is under active concurrent peer edit (`0d101352`); coordinate before touching `solver/core`.
- **Evidence tier**: empirical (reproduced deterministically, gpu-independent). Not test-gaming: the test asserts real mass-conservation physics; the fix must be in the solver/assembly, never a weakened tolerance or raised iteration cap.
