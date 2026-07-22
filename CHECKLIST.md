# CFDrs Work Checklist

Target version: `0.3.0` (pre-1.0 breaking provider-boundary release).

- [x] CFD-HYPERION-OPTICAL-1 [patch] [arch]: replace the raw 405-nm
      Beer-Lambert report expression with direct Hyperion ownership.
  - [x] Align the Aequitas, Proteus, and Hyperion Git source identities.
  - [x] Retain empirical coefficient and hematocrit policy locally, but move
        coefficient/path validation and transmission evaluation to Hyperion.
  - [x] Delete the raw production exponential and add analytical consumer
        regressions for zero path, hematocrit policy, and finite attenuation.
  - [x] Pass the cfd-optim compile, configured Nextest (132/132),
        warning-denied all-target Clippy, doctest (2 passed, 3 existing
        ignored), warning-denied Rustdoc, dependency identity, and production-
        residue gates. The linked verification lane mounted the canonical,
        gitignored contract fixtures so the full suite ran without exclusions.

- [x] CFD-IRIS-COLOR-1 [major] [arch]: make Iris the normalized color-law
  owner for schematic analysis overlays.
  - [x] Add Iris as a direct dependency and remove the local map enum and
    blue-red, Viridis, and grayscale formulas without a compatibility layer.
  - [x] Accept borrowed or owned node/edge maps through `Cow`, keep map storage
    private, validate finite values once, and precompute scalar ranges.
  - [x] Migrate every cfd-1d, cfd-2d, and root example call site to the direct
    Iris contract and fallible builders.
  - [x] Verify exact colors, byte quantization, borrowed identity, range
    normalization, constant fields, missing IDs, and non-finite rejection;
    focused Nextest passes 176/176, warning-denied all-target/all-feature
    Clippy passes, all affected examples compile, 16 doctests pass, and
    warning-denied Rustdoc is clean.
  - [x] Execute the real Venturi 2D solver and inspect its pressure-overlay PNG;
    computed high/low pressure colors, units, axes, and legend are coherent.
  - [x] Attempt major SemVer classification. The tool's isolated dependency
    graph cannot build cfd-core because CFDrs and Proteus/Hephaestus resolve
    different pre-existing Aequitas and Leto Git revisions; the intentional
    removal and builder changes remain classified `[major]`.

- [x] CFD-LAPLACIAN-PROVIDER-1 [major] [arch]: replace the complete cfd-math
  two-dimensional CPU/GPU solver pair with Leto/Hephaestus execution.
  - [x] Add the Leto typed stencil and CPU implementation upstream.
  - [x] Derive Hephaestus parameters from the shared boundary and polarity.
  - [x] Delete the cfd-math CPU formula and cfd-core test-only oracle.
  - [x] Add full-grid anisotropic Neumann quadratic regressions.
  - [x] Pass focused formatting, all-feature and CPU-only checks,
    warning-denied Clippy, configured Nextest (622/622), doctests (6 passed,
    2 existing ignored), warning-denied Rustdoc, and the updated example check.
    `cargo-semver-checks` was attempted but its nightly Rustdoc subprocess
    remained blocked by a long-lived shared-target Leto IDE check; the public
    constructor break is classified `[major]` in the ADR and changelog.

- [x] CFD-SPARSE-DIRECT-OWNERSHIP-1 [patch]: Reconcile the stale
  rsparse-removal work without changing solver semantics.
  - [x] Restore the rsparse workspace/package dependency and exact sparse-LU
    implementation.
  - [x] Record `LETO-SPARSE-DIRECT-1` as the required upstream replacement
    boundary.
  - [x] Run focused formatting, warning-denied Clippy, and direct/consumer
    Nextest gates; synchronize PM evidence. Package formatting and cfd-math
    all-target/all-feature Clippy pass; direct-solver Nextest passes 4/4; the
    cfd-2d independent direct/GMRES regression passes 1/1. Workspace-wide
    formatting remains blocked before source inspection by Windows error 206.

- [x] CFD-EXAMPLE-CLIPPY-1 [patch]: Replace false root validation reports with
  canonical cfd-1d/cfd-2d computations, remove three unreferenced reports that
  contained preset results, and restore warning-denied root example builds.
  Evidence: four retained examples execute with calculated outputs; root
  all-target Clippy passes with warnings denied; focused cfd-1d nextest and
  docs pass. Next increment: CFD-3D-BIFURCATION-BOUNDARIES-1 owns labeled SDF
  terminal facets before a 3D bifurcation example can return.

- [x] `cfd-core` [major]: Delete the remaining public WGPU adapter/feature
  capability surface from `GpuContext`; acquire and query through Hephaestus
  traits, preserve the seven-binding downlevel limit contract, and serialize
  only process-global GPU tests through Nextest. Acceptance: empty raw API
  consumer scan; value-semantic limit regression; complete cfd-core GPU,
  cfd-math GPU, cfd-2d GPU, and root integration nextest; warning-denied
  diagnostics; doctest/rustdoc; and pre-1.0 SemVer classification.

- [x] `cfd-core` [major]: Make the Hephaestus buffer handle crate-private and
  delete `GpuBuffer::buffer`, which exposed a raw WGPU buffer despite having no
  external consumer. Evidence: GPU-enabled cfd-core nextest 244/244 and source
  audit shows the handle is used only by the typed Laplacian dispatch.
- [x] `cfd-core`/`cfd-2d` [major]: Remove public raw WGPU device, queue, and
  limit fields from `GpuContext`; replace `GpuPoissonSolver::new(device,
  queue, ...)` with `from_context(&GpuContext, ...)`. Evidence: GPU-enabled
  accelerated-solver nextest 2/2, cfd-core nextest 244/244, and exact consumer
  scans contain no old constructor or context device/queue access.
- [x] `cfd-core` [patch]: Route `GpuContext::synchronize` through the
  Hephaestus `ComputeDevice` contract instead of calling the raw WGPU polling
  API. Evidence: touched-file `rustfmt --check`, GPU-enabled core nextest
  (244/244), and warning-denied GPU all-target Clippy pass; the GPU compute
  tree contains no direct `device.poll` call.
- [x] `cfd-schematics` [patch]: Replace the unnamed serpentine-Venturi tuple
  with a typed geometry record and repair all warning-denied test/example
  diagnostics. Exact bit equality preserves the direct-value contract for
  geometry bounds and bilateral phase direction. Evidence: all-target,
  all-feature Clippy passes with warnings denied; focused nextest and doctests
  pass.
- [x] DEP-657-01 [patch]: Pin Moirai to merged `main`
  `5ead788c70c728d971237d7afa0b915ea7cf87e3` and regenerate the locked Atlas
  provider graph. Evidence: locked metadata and all-feature `cfd-schematics`
  check pass with Moirai 0.4 and Themis 0.10. Merged as PR #292.
- [x] DEP-657-02 [patch]: Preserve the peer's merged Leto `main` source pin
  `6aedde0c7835238867d6f3cd17b030f7e69cb6f2` and update the accompanying
  Moirai pin to `5ead788c70c728d971237d7afa0b915ea7cf87e3`. Evidence: locked
  metadata plus all-feature `cfd-schematics` check and warning-denied Clippy
  pass after merging current CFDrs `main`.
- [x] `cfd-1d`/`cfd-3d` [patch]: Delete the two stale scalar-file entries from
  `xtask/legacy_surface.allowlist`; exact source scans are clean and the legacy
  migration audit reports no cleanup candidates while retaining clean status.
- [x] `cfd-math`/`cfd-1d`/`cfd-2d` [major]: Delete the dead
  `spmv_parallel`, `use_parallel_spmv`, and `with_parallel_spmv` compatibility
  surface after Leto became the single CSR SpMV owner; migrate every live
  solver initializer to the canonical operation and delete the duplicate
  benchmark/example. Evidence: three-package compile passes; focused SpMV
  nextest passes 7/7; production static audit has no removed symbols or stale
  Rayon wording. Full warning-denied SpMV gates were not rerun in the
  compute-dispatch slice; `cfd-core` compute-dispatch clippy is now clean.
- [x] `cfd-core` [patch]: Remove the residual compute-dispatch CPU downgrade
  path. `ComputeDispatcher` now executes the requested backend, reports
  unsupported kernel/backend combinations as `UnsupportedOperation`, and
  propagates GPU provider failures instead of recomputing on CPU. Evidence:
  focused dispatcher nextest passes 4/4; `cfd-core` all-target clippy passes
  with warnings denied; exact migration-target scan finds only audit-tool
  tokens plus the root Moirai replacement comment; legacy migration audit
  reports allowlist status clean.
- [x] `cfd-core` [arch]: Delete the unconsumed raw-WGPU
  `GpuPipelineManager`, `GpuContext::create_compute_pipeline_with_layout`, and
  obsolete `GpuKernel<T>` trait after all live operations moved to Hephaestus.
  Evidence: static consumer and raw-pipeline audits are empty; cfd-core
  nextest passes 243/243; GPU/no-default and downstream cfd-2d all-target checks,
  warning-denied clippy, doctests 3/3, docs, and migration audit pass.
- [x] `cfd-core`/`cfd-2d` [arch]: Consolidate Smagorinsky viscosity, DES grid
  scale, and rectangular wall distance behind one Hephaestus turbulence facade;
  delete raw-WGPU/fake-generic kernels, duplicate buffer caches, unused DES
  velocity inputs, hardcoded speedup estimates, the nondeterministic wall-clock
  test, and the precision-changing f64 GPU bridge. Evidence: focused provider
  nextest passes 4/4; full `cfd-core` passes 243/243; full `cfd-2d` passes
  570/570 with 27 existing skips; warning-denied dual-package clippy, legacy
  benchmark compilation, doctests, and migration audit pass.
- [x] `cfd-core` [arch]: Replace the cosmetic generic/raw-WGPU pressure type
  with real Hephaestus `f32` weighted-Jacobi and residual operations plus
  validated `PressureConfig`; correct edge/corner Neumann application and
  consolidate shared 3D dispatch construction. Evidence: focused pressure
  nextest passes 6/6; full `cfd-core` nextest passes 247/247; GPU and
  no-default checks pass; all-target clippy passes with warnings denied;
  doctests pass 3/3; docs, migration allowlist, and static audits are clean.
- [x] `cfd-core` [arch]: Replace the cosmetic generic/raw-WGPU velocity type
  with real Hephaestus `f32` correction and divergence-source operations plus
  validated `VelocityConfig`; delete the unsupported compute-trait and raw
  pipeline surfaces. Evidence: focused velocity nextest passes 5/5; full
  `cfd-core` nextest passes 242/242; GPU and no-default checks pass; all-target
  clippy passes with warnings denied; doctests pass 3/3; docs are
  warning-clean; migration allowlist and raw-provider/fake-generic audits are
  clean.
- [x] `cfd-core` [arch]: Replace the cosmetic generic/raw-WGPU diffusion type
  with a real Hephaestus `f32` kernel and validated `DiffusionConfig`; delete
  the unsupported compute-trait body and unbound raw-pipeline surface.
  Completion requires exact constant/quadratic/boundary/partial-workgroup
  tests, typed invalid input and stability tests, clean provider/fake-generic
  residue audits, package format/check/clippy/nextest/doctest/doc gates, and
  synchronized artifacts. Evidence: focused diffusion nextest passes 4/4;
  full `cfd-core` nextest passes 238/238; GPU and no-default checks pass;
  all-target clippy passes with warnings denied; doctests pass 3/3; docs are
  warning-clean; migration allowlist and raw-provider/fake-generic audits are
  clean.
- [x] `cfd-core` [arch]: Replace the cosmetic generic/raw-WGPU advection type
  with a real Hephaestus `f32` kernel and validated `AdvectionConfig`; delete
  the unsupported compute-trait body and separate raw-pipeline test. Completion
  requires exact directional/boundary/partial-workgroup tests, typed invalid
  input and CFL tests, clean provider/fake-generic residue audits, package
  format/check/clippy/nextest/doctest/doc gates, and synchronized artifacts.
  Evidence: focused advection nextest passes 6/6; full `cfd-core` nextest passes
  234/234; GPU and no-default checks pass; all-target clippy passes with
  warnings denied; doctests pass 3/3; docs are warning-clean; migration
  allowlist and raw-provider/fake-generic audits are clean.
- [x] `cfd-core`/`cfd-math` [arch]: Route the 2D GPU Laplacian through
  Hephaestus typed multi-storage dispatch, remove raw WGPU orchestration and
  silent CPU fallback, and narrow the downstream operator to the WGSL kernel's
  real `f32` contract. Completion requires exact/differential boundary tests,
  typed invalid-contract tests, static residue audit, package checks, clippy,
  nextest, doctests, and warning-clean docs. Evidence: GPU/no-default checks
  pass; focused Laplacian nextest passes 10/10; full `cfd-core` nextest passes
  231/231; full `cfd-math` nextest passes 362/362; all-target GPU clippy passes
  with warnings denied; doctests pass 6/6 with 3 intentionally ignored; docs
  are warning-clean; raw pipeline/staging/fallback and fake-generic audits are
  clean. The pass also resolves one pre-existing Anderson clippy diagnostic.
  `cargo-semver-checks` cannot select a registry baseline because `cfd-core`
  and `cfd-math` are not published; the breaking API delta is recorded in the
  ADR and CHANGELOG migration section.
- [x] Extend the provider boundary with Aequitas `Length<f32>` so callers
  cannot pass dimensionless grid spacing. CFDrs validates finite positive
  metre values before dispatch, while Hephaestus performs the single
  quantity-to-POD conversion. Evidence: GPU-feature checks and warning-denied
  all-target Clippy pass; focused Laplacian Nextest passes 13/13; doctests pass
  6/6 with 2 intentionally ignored.
- [x] `cfd-core` [arch]: Delete duplicated raw-WGPU field addition and scalar
  multiplication kernels and route the fallible `GpuFieldOps` arithmetic
  facade through Hephaestus typed elementwise operations. Exact tests cover
  partial workgroups through length 257, empty inputs, and typed length errors;
  static audit finds no arithmetic WGSL/raw WGPU/fallback residue. Evidence in
  `D:/atlas/repos/CFDrs`: no-default and GPU checks pass, all-target GPU clippy
  passes, full cfd-core nextest passes 230/230 with no skips, doctests pass 5/5,
  and package docs are warning-clean.
- [x] Remove direct `num-traits` from the `cfd-1d` and `cfd-3d` scalar
  seams. `Cargo.toml` no longer declares `num-traits` in the workspace
  dependency catalog, `crates/cfd-1d/Cargo.toml` and
  `crates/cfd-3d/Cargo.toml` no longer depend on it directly, and
  `Cfd1dScalar`/`Cfd3dScalar` now provide `zero()`/`one()` through Eunomia
  `NumericElement` constants instead of `num_traits::{Zero,One}` bounds.
  Evidence in `D:/atlas/repos/CFDrs`: touched-file rustfmt passed,
  `rustup run nightly cargo check -p cfd-1d -p cfd-3d` passed, targeted
  residue scan found no `num_traits`/direct `num-traits` hits in the touched
  manifests and scalar cones, and `rustup run nightly cargo nextest run -p
  cfd-1d -p cfd-3d --status-level fail` passed 1122/1122 with one existing
  slow 3D mesh-convergence validation. Residual: package-wide
  `cargo fmt -p cfd-1d -p cfd-3d --check` and all-targets clippy remain
  blocked by pre-existing unrelated formatting/lint drift outside this scalar
  dependency slice.
- [x] Move `cfd-1d::solver::core::SolverWorkspace` reusable vector storage
  off nalgebra `DVector`. `SolverWorkspace::{rhs,last_solution,
  linear_solution}` now store `leto::Array1<T>`, `MatrixAssembler::assemble`
  returns a Leto RHS array, and `validate_linear_system`/residual checks
  consume the Leto RHS directly. The one remaining nalgebra vector crossing in
  this cone is centralized in private `vector_bridge`, which feeds the current
  nalgebra-returning linear-system and Anderson Picard solution boundary.
  Evidence in `D:/atlas/repos/CFDrs`: `rustup run nightly cargo check -p
  cfd-1d --tests`, `rustup run nightly cargo clippy -p cfd-1d --lib -- -D
  warnings`, focused solver-core nextest over matrix assembly, Dirichlet
  enforcement, linearization, and core solver cases (10/10), and `rustup run
  nightly cargo doc -p cfd-1d --no-deps` completed with the same 11
  pre-existing rustdoc link warnings. Residual: `linear_system.rs` still
  exposes nalgebra `DVector` solutions, dense `DMatrix` fallback, and the
  nalgebra-sparse matrix input boundary; `matrix_assembly.rs` still returns
  `nalgebra_sparse::CsrMatrix`; `solver_detection.rs` and
  `anderson_acceleration.rs` still operate on nalgebra solution vectors at the
  current linear-solver boundary.
- [x] Move `cfd-1d::solver::core::NetworkState` public storage off nalgebra
  `DVector`. `NetworkState::{pressures,flow_rates}` now expose Leto
  `Array1<T>` storage, `new`/`from_network` construct Leto arrays directly,
  and focused unit coverage asserts shape, values, clone preservation, and
  time mutation semantics over the Leto-backed state. Evidence in
  `D:/atlas/repos/CFDrs`: `rustup run nightly cargo check -p cfd-1d --tests`,
  `rustup run nightly cargo clippy -p cfd-1d --lib -- -D warnings`, focused
  `rustup run nightly cargo nextest run -p cfd-1d -E
  'test(new_allocates_leto_state_vectors_with_zero_values) |
  test(cloned_state_preserves_leto_vector_values)'` (2/2), and `rustup run
  nightly cargo doc -p cfd-1d --no-deps` completed with the same 11
  pre-existing rustdoc link warnings. The upstream Atlas blocker was closed by
  adding `FloatElement::acos` in `D:/atlas/repos/eunomia`; targeted residue
  scan found no `DVector` in `crates/cfd-1d/src/solver/core/state.rs`.
  Residual: the broader 1D linear-system, sparse matrix assembly,
  solver-detection, and Anderson solution-vector boundary still own nalgebra
  `DVector`/`DMatrix`/`nalgebra_sparse` surfaces.
- [x] Move `cfd-1d::solver::core::ConvergenceChecker` off the nalgebra
  `DVector` public boundary. `ConvergenceChecker::{check,has_converged,
  has_converged_dual}` now accepts Leto `Array1<T>` values, computes L2 norms
  through direct value semantics, and returns typed `InvalidInput` on
  mismatched vector lengths. The existing 1D network solver workspace remains
  nalgebra-backed for this slice and converts once at the convergence-check
  call boundary. Evidence in `D:/atlas/repos/CFDrs`: `rustup run nightly cargo
  check -p cfd-1d --tests`, `rustup run nightly cargo clippy -p cfd-1d --lib
  --test solver_core_tests -- -D warnings`, and `rustup run nightly cargo
  nextest run -p cfd-1d convergence` (20/20). `rustup run nightly cargo doc
  -p cfd-1d --no-deps` generated docs but still emitted 11 pre-existing
  broken/private intra-doc link warnings outside this convergence slice.
  Residual: the broader 1D linear-system, sparse matrix assembly,
  solver-detection, and Anderson solution-vector boundary still own nalgebra
  `DVector`/`DMatrix`/`nalgebra_sparse` surfaces.
- [x] Remove direct root `num-traits` and `simba` dependencies. `Cargo.toml`
  no longer declares `num-traits` or `simba` in `[workspace.dependencies]` or
  the root package `[dependencies]`; current source and member manifests use
  Eunomia/Atlas scalar contracts instead of direct `num_traits` or `simba`
  imports. Evidence in `D:/atlas/repos/CFDrs`: targeted source/member
  manifest residue scan, `rustup run nightly cargo metadata --no-deps
  --format-version 1`. Residual: lockfile `num-traits` and `simba` entries
  remain as transitive upstream dependencies; `rustup run nightly cargo check
  -p cfd-suite --no-default-features` was queued behind concurrent shared
  Cargo cache/build locks and then terminated without compiler diagnostics;
  nextest did not run.
- [x] Remove direct root/`cfd-3d` `crossbeam` dependency. `Cargo.toml` no
  longer declares `crossbeam` in `[workspace.dependencies]`, and
  `crates/cfd-3d/Cargo.toml` no longer depends on `crossbeam.workspace`.
  Current source and member manifests have no direct Crossbeam imports or
  workspace dependency users. Evidence in `D:/atlas/repos/CFDrs`: targeted
  source/member manifest residue scan, `rustup run nightly cargo metadata
  --no-deps --format-version 1`, `rustup run nightly cargo check -p cfd-3d
  --no-default-features`, `rustup run nightly cargo nextest run -p cfd-3d
  --no-default-features` (394/394), and `git diff --check -- Cargo.toml
  crates/cfd-3d/Cargo.toml CHANGELOG.md CHECKLIST.md backlog.md
  gap_audit.md`. Residual: transitive Crossbeam dependencies, if any, remain
  owned by upstream crates; the focused nextest run marked existing
  mesh-convergence validation slow at 19.7s.
- [x] Close the public sparse/linear-solver Atlas boundary residual.
  `crates/cfd-math/src/sparse` now exposes `leto_ops::CsrMatrix` through the
  public `SparseMatrix<T>` alias and direct builder/operation paths, and
  `LinearOperator`, `Preconditioner`, `LinearSolver`, direct solver,
  solver-chain, and migrated fixtures no longer expose nalgebra `DVector` or
  `nalgebra_sparse::CsrMatrix` in the requested cone. `cfd-validation`
  numerical solver validation now stores computed/analytical vectors as
  `leto::Array1<T>` and constructs sparse test matrices through the Leto CSR
  boundary. Evidence in `D:/atlas/repos/CFDrs`: `rustup run nightly cargo
  check -p cfd-math`, `rustup run nightly cargo check -p cfd-math --tests`,
  `rustup run nightly cargo clippy -p cfd-math --all-targets -- -D
  warnings`, `rustup run nightly cargo doc -p cfd-math --no-deps`, `rustup
  run nightly cargo nextest run -p cfd-math` (361/361), `rustup run nightly
  cargo check -p cfd-validation`, `rustup run nightly cargo clippy -p
  cfd-validation --all-targets -- -D warnings`, `rustup run nightly cargo doc
  -p cfd-validation --no-deps`, and a targeted residue scan with no
  `nalgebra_sparse::CsrMatrix`, public nalgebra-sparse re-export, `DVector`,
  `row_offsets()`, `try_from_csr_data`, `CooMatrix`, or `nalgebra::` matches
  under the migrated sparse/linear-solver/validation files. Residual:
  `rustup run nightly cargo nextest run -p cfd-validation` still fails in
  existing venturi cross-fidelity convergence tests
  `option2_selected_45um_geometry_routes_to_fallback_and_converges` and
  `microventuri_35um_case_produces_converged_informative_2d_result`; those are
  outside the sparse/linear-solver boundary.
- [x] Move cfd-math IncompleteCholesky preconditioner to Leto CSR construction.
  `crates/cfd-math/src/linear_solver/preconditioners/cholesky.rs` now owns and
  factorizes `leto_ops::CsrMatrix`, replaces nalgebra `get_entry` with Leto
  CSR row-value lookup, constructs IC(0) factors with `CsrMatrix::from_parts`,
  and performs forward/backward substitution through Leto CSR row views.
  Source and integration preconditioner edge tests now pass Leto CSR into
  `IncompleteCholesky::new`. Evidence in `D:/atlas/repos/CFDrs`: `rustup run
  nightly cargo fmt -p cfd-math --check`, `rustup run nightly cargo check -p
  cfd-math --lib`, `rustup run nightly cargo check -p cfd-math
  --all-targets`, `rustup run nightly cargo clippy -p cfd-math --lib --tests
  -- -D warnings`, `rustup run nightly cargo clippy -p cfd-math --all-targets
  -- -D warnings`, `rustup run nightly cargo nextest run -p cfd-math
  cholesky --status-level fail` (5/5), `rustup run nightly cargo nextest run
  -p cfd-math preconditioner --status-level fail` (76/76), and targeted
  Cholesky residue scans. Residual: Schwarz, direct solver, remaining
  transitional solver fixtures, and the shared `crate::sparse::SparseMatrix`
  solver matrix boundary still use nalgebra-sparse.
- [x] Move cfd-math ILU preconditioner family to Leto CSR construction.
  `crates/cfd-math/src/linear_solver/preconditioners/ilu` now stores and
  factorizes `leto_ops::CsrMatrix` for ILU(0), ILU(k), and triangular solve
  application. ILU construction sites in source tests, integration edge tests,
  `LinearSolverChain`, and Schwarz local solves now pass Leto CSR matrices or
  convert once at the remaining shared solver/Schwarz boundary. The ILU module
  has no direct `nalgebra_sparse`, `row_offsets()`, `try_from_csr_data`,
  `get_entry()`, or `SparseEntry` residue. Evidence in `D:/atlas/repos/CFDrs`:
  `rustup run nightly cargo fmt -p cfd-math --check`, `rustup run nightly
  cargo check -p cfd-math --lib`, `rustup run nightly cargo check -p
  cfd-math --all-targets`, `rustup run nightly cargo clippy -p cfd-math
  --lib --tests -- -D warnings`, `rustup run nightly cargo clippy -p
  cfd-math --all-targets -- -D warnings`, `rustup run nightly cargo nextest
  run -p cfd-math ilu --status-level fail` (21/21), `rustup run nightly cargo
  nextest run -p cfd-math preconditioner --status-level fail` (76/76),
  `rustup run nightly cargo nextest run -p cfd-math "linear_solver::tests"
  --status-level fail` (53/53), and targeted ILU residue scans. Residual:
  Schwarz, direct solver, remaining transitional solver fixtures, and the
  shared `crate::sparse::SparseMatrix` solver matrix boundary still use
  nalgebra-sparse.
- [x] Move cfd-math SSOR preconditioner to Leto CSR construction.
  `crates/cfd-math/src/linear_solver/preconditioners/ssor.rs` now owns
  `leto_ops::CsrMatrix`, uses Leto CSR row views for forward/backward sweeps,
  and has no direct `nalgebra_sparse`, `row_offsets()`, `try_from_csr_data`,
  or `get_entry()` residue. Source preconditioner edge tests convert their
  still-transitional solver CSR fixtures once into Leto CSR before
  constructing SSOR, and the SSOR vector-length regression now asserts the
  typed `InvalidConfiguration` messages. Evidence in `D:/atlas/repos/CFDrs`:
  `rustup run nightly cargo fmt -p cfd-math --check`, `rustup run nightly
  cargo check -p cfd-math --lib`, `rustup run nightly cargo check -p
  cfd-math --all-targets`, `rustup run nightly cargo clippy -p cfd-math
  --lib --tests -- -D warnings`, `rustup run nightly cargo clippy -p
  cfd-math --all-targets -- -D warnings`, `rustup run nightly cargo nextest
  run -p cfd-math ssor --status-level fail` (5/5), `rustup run nightly cargo
  nextest run -p cfd-math preconditioner --status-level fail` (76/76), and
  targeted residue scans. Residual: Schwarz, direct solver,
  integration-test fixtures, and the shared `crate::sparse::SparseMatrix`
  solver matrix boundary still use nalgebra-sparse.
- [x] Move cfd-math basic preconditioners to Leto CSR construction.
  `crates/cfd-math/src/linear_solver/preconditioners/basic.rs` now takes
  `leto_ops::CsrMatrix` for Jacobi and SOR construction and has no direct
  `nalgebra_sparse`, `SparseMatrixExt`, `row_offsets()`, `try_from_csr_data`,
  or `get_entry()` residue. Jacobi diagonal extraction uses Leto CSR
  `diagonal()`, and SOR row sweeps use Leto CSR row views while
  `Preconditioner::apply_to` remains on `leto::Array1`. Source linear-solver
  tests and `core_solver_tests.rs` convert their still-transitional solver CSR
  fixtures once into Leto CSR before constructing Jacobi/SOR. Evidence in
  `D:/atlas/repos/CFDrs`: `rustup run nightly cargo fmt -p cfd-math --check`,
  `rustup run nightly cargo check -p cfd-math --lib`, `rustup run nightly
  cargo check -p cfd-math --test core_solver_tests`, `rustup run nightly
  cargo check -p cfd-math --all-targets`, `rustup run nightly cargo clippy -p
  cfd-math --lib --tests -- -D warnings`, `rustup run nightly cargo clippy -p
  cfd-math --test core_solver_tests -- -D warnings`, `rustup run nightly
  cargo clippy -p cfd-math --all-targets -- -D warnings`, `rustup run
  nightly cargo nextest run -p cfd-math "linear_solver::tests" --status-level
  fail` (53/53), `rustup run nightly cargo nextest run -p cfd-math --test
  core_solver_tests --status-level fail` (4/4), `rustup run nightly cargo
  nextest run -p cfd-math preconditioner --status-level fail` (76/76), and
  targeted residue scans. Residual: Schwarz, direct solver, and the shared
  `crate::sparse::SparseMatrix` solver matrix boundary still use
  nalgebra-sparse.
- [x] Move cfd-math AMG/coarsening sparse boundary to Leto CSR.
  `crates/cfd-math/src/linear_solver/preconditioners/multigrid` now owns a
  local `leto_ops::CsrMatrix` sparse alias plus construction/value helpers,
  and AMG setup/coarsening/interpolation/smoothers/cycles route sparse
  products, transpose, SpMV, and row access through Leto CSR APIs. `crates/
  cfd-math/benches/{coarsening_bench,algebraic_distance_bench}.rs` and
  `tests/amg_coarsening_tests.rs` now construct Leto CSR matrices directly.
  `tests/amg_integration_test.rs` now constructs the AMG preconditioner from
  Leto CSR while the still-transitional Krylov solver matrix is assembled
  through the existing sparse builder, not direct nalgebra-sparse APIs.
  Evidence in `D:/atlas/repos/CFDrs`: `rustup run nightly cargo fmt -p
  cfd-math --check`, `rustup run nightly cargo check -p cfd-math --lib`,
  `rustup run nightly cargo check -p cfd-math --test amg_coarsening_tests`,
  `rustup run nightly cargo check -p cfd-math --test amg_integration_test`,
  `rustup run nightly cargo check -p cfd-math --bench coarsening_bench`,
  `rustup run nightly cargo check -p cfd-math --bench
  algebraic_distance_bench`, focused clippy for cfd-math lib, the two AMG
  tests, and both migrated benches with `-D warnings`, `rustup run nightly
  cargo nextest run -p cfd-math --test amg_integration_test --status-level
  fail` (5/5), `rustup run nightly cargo nextest run -p cfd-math amg
  --status-level fail` (6/6), `rustup run nightly cargo nextest run -p
  cfd-math multigrid::coarsening --status-level fail` (10/10), and a targeted
  AMG/coarsening sparse residue scan showing no `nalgebra_sparse`,
  `CooMatrix`, `try_from_csr_data`, `get_entry`, `crate::sparse::spmv`,
  `spmv_array`, or `row_offsets()` residue in the migrated boundary.
  Residual: cfd-math still exposes nalgebra-sparse through the broader
  `crate::sparse::SparseMatrix` solver/direct/preconditioner boundary; the
  `LinearSolverChain` AMG tier still converts that boundary once into Leto CSR
  until the whole solver matrix surface moves.
- [x] Migrate cfd-math SpMV/CG benchmarks to direct Leto CSR.
  `crates/cfd-math/benches/{spmv_bench,cg_bench,math_benchmarks}.rs` now
  construct `leto_ops::CsrMatrix` values directly with
  `CsrMatrix::from_parts` instead of `nalgebra_sparse::CsrMatrix` and
  `try_from_csr_data`. `spmv_bench.rs` now measures the direct Leto CSR
  `LinearOperator::apply` path rather than the legacy `cfd_math::sparse::spmv`
  nalgebra boundary. Evidence in `D:/atlas/repos/CFDrs`: `rustup run nightly
  cargo fmt -p cfd-math --check`, `rustup run nightly cargo check -p
  cfd-math --bench spmv_bench`, `rustup run nightly cargo check -p cfd-math
  --bench cg_bench`, `rustup run nightly cargo check -p cfd-math --bench
  math_benchmarks`, focused clippy for each migrated bench with `-D warnings`,
  `rustup run nightly cargo check -p cfd-math --all-targets`, `rustup run
  nightly cargo clippy -p cfd-math --all-targets -- -D warnings`, `rustup run
  nightly cargo nextest run -p cfd-math sparse --status-level fail` (19/19),
  and a targeted scan showing no `nalgebra_sparse`, `CooMatrix`, `DVector`,
  `DMatrix`, `cfd_math::sparse::spmv`, or `try_from_csr_data` residue in the
  migrated benchmark files. Residual after Sprint 1.96.162 moved the AMG
  coarsening/algebraic-distance boundary to Leto CSR: the broader
  solver/direct/preconditioner sparse matrix surface still uses the
  transitional nalgebra-sparse `SparseMatrix<T>` boundary.
- [x] Add cfd-math Leto CSR linear-operator path. `crates/cfd-math/src/
  sparse/operations.rs` now implements `LinearOperator<T>` directly for
  `leto_ops::CsrMatrix<T>` and routes the legacy `nalgebra_sparse::CsrMatrix`
  SpMV path through the same `try_leto_spmv` helper after conversion. `crates/
  cfd-math/tests/simple_gmres_tests.rs` now constructs Leto CSR matrices
  directly and verifies GMRES residuals through the Atlas-native operator path,
  with no nalgebra sparse/vector residue in that integration test. Evidence in
  `D:/atlas/repos/CFDrs`: `rustup run nightly cargo fmt -p cfd-math --check`,
  `rustup run nightly cargo check -p cfd-math --test simple_gmres_tests`,
  `rustup run nightly cargo nextest run -p cfd-math --test
  simple_gmres_tests --status-level fail` (3/3), `rustup run nightly cargo
  clippy -p cfd-math --test simple_gmres_tests -- -D warnings`, `rustup run
  nightly cargo check -p cfd-math --lib`, `rustup run nightly cargo check -p
  cfd-math --all-targets`, `rustup run nightly cargo clippy -p cfd-math
  --all-targets -- -D warnings`, `rustup run nightly cargo nextest run -p
  cfd-math sparse --status-level fail` (19/19), `rustup run nightly cargo
  nextest run -p cfd-math gmres --status-level fail` (21/21), and a targeted
  `simple_gmres_tests.rs` scan showing no `nalgebra_sparse`, `CooMatrix`,
  `nalgebra::`, `DVector`, `DMatrix`, or `cfd_math::sparse` residue. Residual:
  cfd-math still exposes the public nalgebra-sparse CSR boundary in many
  production preconditioner/sparse APIs and tests; this slice opens the
  direct Leto CSR solver path and migrates one integration test to prove it.
- [x] Close remaining cfd-math Leto storage-slice residue. `crates/cfd-math/
  src/nonlinear_solver/linalg.rs` no longer exposes mutable `StorageMut` slice
  helpers, and `crates/cfd-math/src/nonlinear_solver/anderson.rs` now performs
  dense pivot swaps through direct `Array1`/`Array2` indexing.
  `crates/cfd-math/src/linear_solver/preconditioners/multigrid/
  interpolation.rs` now evaluates interpolation constant-preservation error
  through direct `Array1` indexing, and `smoothers.rs` test assertions use
  direct indexed comparisons instead of `storage().as_slice()`. Evidence in
  `D:/atlas/repos/CFDrs`: `rustup run nightly cargo fmt -p cfd-math
  --check`, `rustup run nightly cargo check -p cfd-math --lib`, `rustup run
  nightly cargo nextest run -p cfd-math nonlinear_solver multigrid
  --status-level fail` (46/46), `rustup run nightly cargo clippy -p cfd-math
  --lib --tests -- -D warnings`, `rustup run nightly cargo check -p cfd-math
  --all-targets`, `rustup run nightly cargo clippy -p cfd-math --all-targets
  -- -D warnings`, and a cfd-math `src`/`tests` scan showing no
  `leto::Storage`, `StorageMut`, `.storage().as_slice()`, `as_slice_mut()`,
  `vector_slice_mut`, or `matrix_slice_mut` residue. Residual cfd-math
  provider work is now outside the Leto storage-slice cleanup and includes
  `nalgebra_sparse::CsrMatrix`, nalgebra dense/test oracles, and remaining
  scalar/provider boundaries.
- [x] Remove cfd-math sparse/basic preconditioner Leto storage-slice
  dependencies. `crates/cfd-math/src/sparse/operations.rs` now stages SpMV
  input/output and row/column scaling buffers through direct `Array1` indexing
  before delegating to Leto CSR operations, removing the remaining
  `leto::Storage`, `.storage().as_slice()`, and output `as_slice_mut()`
  assumptions from the sparse operation boundary. `crates/cfd-math/src/
  linear_solver/preconditioners/basic.rs` now indexes the Leto diagonal
  directly when building Jacobi inverse diagonals. Evidence in
  `D:/atlas/repos/CFDrs`: `rustup run nightly cargo fmt -p cfd-math
  --check`, `rustup run nightly cargo check -p cfd-math --lib`, `rustup run
  nightly cargo nextest run -p cfd-math sparse preconditioner --status-level
  fail` (95/95), `rustup run nightly cargo clippy -p cfd-math --lib --tests
  -- -D warnings`, `rustup run nightly cargo check -p cfd-math
  --all-targets`, `rustup run nightly cargo clippy -p cfd-math --all-targets
  -- -D warnings`, and targeted scans showing no `Storage`,
  `.storage().as_slice()`, `as_slice_mut()`, or obsolete SpMV contiguity
  diagnostic residue in the migrated sparse/basic files. Residual
  storage-slice owners are `crates/cfd-math/src/nonlinear_solver/linalg.rs`
  mutable `StorageMut` helpers and multigrid interpolation/smoother internals.
- [x] Remove cfd-math GPU operator Leto storage-slice dependency.
  `crates/cfd-math/src/linear_solver/operators/gpu.rs` now validates Leto
  `Array1` input/output lengths at the operator boundary, stages GPU upload
  and readback buffers through direct `Array1` indexing, and writes results
  back by index instead of importing `leto::Storage`, borrowing
  `.storage().as_slice()`, or requiring output `as_slice_mut()` contiguity.
  The existing execution path remains the Hephaestus-backed `cfd-core`
  `GpuContext`/`GpuBuffer`/`Laplacian2DKernel` path. Evidence in
  `D:/atlas/repos/CFDrs`: `rustup run nightly cargo fmt -p cfd-math
  --check`, `rustup run nightly cargo check -p cfd-math --features gpu`,
  `rustup run nightly cargo nextest run -p cfd-math --features gpu
  linear_solver::operators --status-level fail` (5/5), `rustup run nightly
  cargo clippy -p cfd-math --features gpu --lib -- -D warnings`,
  `rustup run nightly cargo check -p cfd-math --features gpu --all-targets`,
  `rustup run nightly cargo clippy -p cfd-math --features gpu --all-targets
  -- -D warnings`, and a targeted scan showing no `Storage`,
  `.storage().as_slice()`, `as_slice_mut()`, or obsolete contiguity diagnostic
  residue in `operators/gpu.rs`. Residual storage-slice owners now include
  `crates/cfd-math/src/sparse/operations.rs` and multigrid
  smoother/interpolation internals; broader GPU provider work remains in raw
  WGPU context/kernel ownership outside this operator.
- [x] Remove cfd-math finite-difference operator Leto storage-slice dependencies.
  `crates/cfd-math/src/linear_solver/operators/{poisson,momentum}.rs` now
  read and write `Array1` values through direct indexing for the 2D Laplacian,
  3D Poisson, 1D/2D momentum, and 2D energy operators. This removes
  `leto::Storage`, `.storage().as_slice()`, and output `as_slice_mut()`
  assumptions from those CPU finite-difference linear operators. Evidence in
  `D:/atlas/repos/CFDrs`: `rustup run nightly cargo fmt -p cfd-math --check`,
  `rustup run nightly cargo check -p cfd-math --lib`, `rustup run nightly
  cargo nextest run -p cfd-math "linear_solver::operators" --status-level
  fail` (5/5), `rustup run nightly cargo clippy -p cfd-math --lib --tests --
  -D warnings`, package `rustup run nightly cargo check -p cfd-math
  --all-targets`, package `rustup run nightly cargo clippy -p cfd-math
  --all-targets -- -D warnings`, and a targeted scan showing no `Storage`,
  `.storage().as_slice()`, or `as_slice_mut()` residue in the migrated
  operator files. Residual storage-slice owners now include
  `crates/cfd-math/src/sparse/operations.rs`,
  `crates/cfd-math/src/linear_solver/operators/gpu.rs`, and multigrid
  smoother/interpolation internals.
- [x] Remove cfd-math nonlinear linalg immutable Leto storage-slice dependency.
  `crates/cfd-math/src/nonlinear_solver/linalg.rs` now evaluates `dot`,
  `add`, `sub`, `add_scaled`, `add_scaled_in_place`, and `scale` through
  direct `Array1` indexing instead of an immutable `vector_slice` helper backed
  by `leto::Storage`. `anderson.rs` now indexes the Anderson `gamma` vector
  directly, so the immutable storage slice helper was deleted. Evidence in
  `D:/atlas/repos/CFDrs`: `rustup run nightly cargo fmt -p cfd-math --check`,
  `rustup run nightly cargo check -p cfd-math --lib`, `rustup run nightly
  cargo nextest run -p cfd-math nonlinear_solver --status-level fail` (9/9),
  `rustup run nightly cargo clippy -p cfd-math --lib --tests -- -D warnings`,
  package `rustup run nightly cargo check -p cfd-math --all-targets`, package
  `rustup run nightly cargo clippy -p cfd-math --all-targets -- -D warnings`,
  and a targeted scan showing no immutable `.storage().as_slice()`,
  `leto::Storage`, or `vector_slice` residue in nonlinear linalg/Anderson.
  Residual mutable dense-workspace helpers still use `StorageMut`; remaining
  immutable storage-slice owners include sparse operations, linear operators,
  GPU operator, and multigrid internals.
- [x] Remove cfd-math production SIMD vector Leto storage-slice dependency.
  `crates/cfd-math/src/simd/vector.rs` now implements `simd_mul`,
  `simd_dot`, `simd_norm`, and `par_map` through direct `Array1` indexing
  while still dispatching work through Moirai `Adaptive` index mapping and
  reduction. This removes the `leto::Storage` import and
  `.storage().as_slice()` calls from the production SIMD vector extension.
  Evidence in `D:/atlas/repos/CFDrs`: `rustup run nightly cargo fmt -p
  cfd-math --check`, `rustup run nightly cargo check -p cfd-math --lib`,
  `rustup run nightly cargo nextest run -p cfd-math simd::vector
  --status-level fail` (1/1), `rustup run nightly cargo clippy -p cfd-math
  --lib --tests -- -D warnings`, package `rustup run nightly cargo check -p
  cfd-math --all-targets`, package `rustup run nightly cargo clippy -p
  cfd-math --all-targets -- -D warnings`, `rustup run nightly cargo nextest
  run -p cfd-math simd --status-level fail` (26/26), and a targeted scan
  showing no `Storage` or `.storage().as_slice()` residue in
  `src/simd/vector.rs`. Residual source-level storage-slice owners now include
  `crates/cfd-math/src/nonlinear_solver/linalg.rs`,
  `crates/cfd-math/src/sparse/operations.rs`,
  `crates/cfd-math/src/linear_solver/operators/{poisson,momentum,gpu}.rs`,
  and multigrid smoother/interpolation internals.
- [x] Remove cfd-math SIMD integration test Leto storage-slice bridge.
  `crates/cfd-math/tests/simd_tests.rs` now reads the Leto `spmv` result by
  direct `Array1` indexing and feeds the existing SIMD slice API from those
  values, removing the `leto::Storage` import and `.storage().as_slice()`
  conversion from the remaining cfd-math integration test residue scan.
  Evidence in `D:/atlas/repos/CFDrs`: `rustup run nightly cargo fmt -p
  cfd-math --check`, `rustup run nightly cargo check -p cfd-math --test
  simd_tests`, `rustup run nightly cargo nextest run -p cfd-math --test
  simd_tests --status-level fail` (12/12), `rustup run nightly cargo clippy
  -p cfd-math --test simd_tests -- -D warnings`, package `rustup run nightly
  cargo check -p cfd-math --all-targets`, package `rustup run nightly cargo
  clippy -p cfd-math --all-targets -- -D warnings`, `rustup run nightly cargo
  nextest run -p cfd-math simd --status-level fail` (26/26), and targeted
  scans showing no `DVector`, nalgebra vector import, local preconditioner
  bridge, `Storage`, or storage-slice conversion residue in
  `crates/cfd-math/tests`. Residual: cfd-math still has production/test
  provider boundaries outside the integration-test vector bridge layer,
  including `nalgebra_sparse::CsrMatrix`, dense nalgebra test oracles, and
  source-level Leto storage-slice internals.
- [x] Migrate cfd-math AMG integration vectors to Leto arrays.
  `crates/cfd-math/tests/amg_integration_test.rs` now constructs exact
  solutions, RHS vectors, solver outputs, AMG cycle outputs, and two-grid
  preconditioner work buffers as `leto::Array1` values directly instead of
  round-tripping through nalgebra `DVector` and Leto storage-slice conversions.
  The two-grid convergence check still uses nalgebra `DMatrix` and
  `SymmetricEigen` as the dense energy-norm oracle; that remaining dense
  eigenvalue provider gap is explicit. Evidence in `D:/atlas/repos/CFDrs`:
  `rustup run nightly cargo fmt -p cfd-math --check`, `rustup run nightly
  cargo check -p cfd-math --test amg_integration_test`, `rustup run nightly
  cargo nextest run -p cfd-math --test amg_integration_test --status-level
  fail` (5/5), `rustup run nightly cargo clippy -p cfd-math --test
  amg_integration_test -- -D warnings`, package `rustup run nightly cargo
  check -p cfd-math --all-targets`, package `rustup run nightly cargo clippy
  -p cfd-math --all-targets -- -D warnings`, `rustup run nightly cargo
  nextest run -p cfd-math amg --status-level fail` (6/6), and a targeted
  residue scan showing no `DVector`, nalgebra vector import, `Storage`,
  storage-slice conversion, or local preconditioner bridge residue in the AMG
  integration test. Residual: `amg_integration_test.rs` still uses
  `nalgebra_sparse::CsrMatrix` for sparse storage and nalgebra dense
  `DMatrix`/`SymmetricEigen` for the energy-norm oracle; the integration-test
  storage-slice holdout was closed by Sprint 1.96.153.
- [x] Migrate cfd-math preconditioner edge-case integration tests to Leto arrays.
  `crates/cfd-math/tests/preconditioner_edge_cases.rs` now constructs
  `leto::Array1` RHS/solution buffers directly for ILU(0), ILU(k),
  repeated-application, extreme-value, and sparsity-preservation coverage
  instead of round-tripping through a local nalgebra `DVector` bridge helper.
  Preconditioner application now stays at the Leto boundary. Evidence in
  `D:/atlas/repos/CFDrs`: `rustup run nightly cargo fmt -p cfd-math --check`,
  `rustup run nightly cargo check -p cfd-math --test preconditioner_edge_cases`,
  `rustup run nightly cargo nextest run -p cfd-math --test
  preconditioner_edge_cases --status-level fail` (6/6), `rustup run nightly
  cargo clippy -p cfd-math --test preconditioner_edge_cases -- -D warnings`,
  package `rustup run nightly cargo check -p cfd-math --all-targets`, package
  `rustup run nightly cargo clippy -p cfd-math --all-targets -- -D warnings`,
  broader `rustup run nightly cargo nextest run -p cfd-math preconditioner
  --status-level fail` (76/76), and a targeted residue scan showing no
  `DVector`, nalgebra vector import, `Storage`, storage-slice conversion, or
  local preconditioner bridge residue in the migrated test. Residual:
  `preconditioner_edge_cases.rs` still uses the current
  `nalgebra_sparse::CsrMatrix` matrix storage boundary; integration-test
  vector bridge holdouts were closed by Sprints 1.96.152 and 1.96.153.
- [x] Migrate cfd-math linear-solver source test module tree to Leto arrays.
  `crates/cfd-math/src/linear_solver/tests/{mod,edge_case_tests,
  adversarial_solver_tests,extended_edge_case_tests}.rs` now constructs and
  verifies `leto::Array1` values directly instead of using local nalgebra
  `DVector` bridge macros. Residual checks use the Leto SpMV/helper path, and
  the touched solver/sparse cone no longer calls old scalar helpers
  `T::zero()`, `T::one()`, `default_epsilon()`, or `to_subset()`. Evidence in
  `D:/atlas/repos/CFDrs`: `rustup run nightly cargo fmt -p cfd-math --check`,
  `rustup run nightly cargo check -p cfd-math --lib`, `rustup run nightly
  cargo check -p cfd-math --all-targets`, `rustup run nightly cargo nextest
  run -p cfd-math linear_solver::tests --status-level fail` (53/53), `rustup
  run nightly cargo nextest run -p cfd-math linear_solver --status-level fail`
  (176/176), `rustup run nightly cargo clippy -p cfd-math --lib --tests -- -D
  warnings`, `rustup run nightly cargo clippy -p cfd-math --all-targets -- -D
  warnings`, plus current targeted scans showing no source-test
  `DVector`/`Storage` bridge residue and no old scalar helper residue in the
  searched solver/sparse cone. Residual: integration-test vector bridge
  holdouts were closed by Sprints 1.96.152 and 1.96.153; production
  sparse/scalar boundaries still include `nalgebra_sparse::CsrMatrix` and
  transitional nalgebra scalar/test utilities.
- [x] Migrate cfd-math core solver validation tests to Leto arrays.
  `crates/cfd-math/tests/core_solver_tests.rs` now constructs
  `leto::Array1` RHS/solution buffers directly for BiCGSTAB, GMRES,
  preconditioner integration, and condition-number robustness coverage instead
  of round-tripping through a local nalgebra `DVector` bridge macro.
  Preconditioner application and residual checks now stay at the Leto array
  boundary and use the existing Leto SpMV path with explicit residual
  thresholds. Evidence in `D:/atlas/repos/CFDrs`: `rustup run nightly cargo
  fmt -p cfd-math --check`, `rustup run nightly cargo check -p cfd-math --test
  core_solver_tests`, `rustup run nightly cargo nextest run -p cfd-math --test
  core_solver_tests --status-level fail` (4/4), `rustup run nightly cargo
  clippy -p cfd-math --test core_solver_tests -- -D warnings`, package
  `rustup run nightly cargo check -p cfd-math --all-targets`, package
  `rustup run nightly cargo clippy -p cfd-math --all-targets -- -D warnings`,
  and a targeted residue scan showing no `DVector`, nalgebra vector,
  `Storage`, storage-slice conversion, local solve/preconditioner bridge,
  matrix-vector `&a * &x`, or vector `.norm()` residue in the core solver test
  module. Residual: `core_solver_tests.rs` still uses the current
  `nalgebra_sparse::CsrMatrix`/`CooMatrix` matrix storage boundary, and broader
  cfd-math integration/adversarial/preconditioner diagnostics still contain
  nalgebra `DVector` bridges.
- [x] Migrate cfd-math simple GMRES integration tests to Leto arrays.
  `crates/cfd-math/tests/simple_gmres_tests.rs` now constructs
  `leto::Array1` RHS/solution buffers directly for basic, restarted, and
  preconditioned GMRES integration tests instead of round-tripping through a
  local nalgebra `DVector` bridge macro. Residual checks now use the Leto SpMV
  path and assert value-semantic residual thresholds for all three tests.
  Evidence in `D:/atlas/repos/CFDrs`: `rustup run nightly cargo fmt -p
  cfd-math --check`, `rustup run nightly cargo check -p cfd-math --test
  simple_gmres_tests`, `rustup run nightly cargo nextest run -p cfd-math
  --test simple_gmres_tests --status-level fail` (3/3), `rustup run nightly
  cargo clippy -p cfd-math --test simple_gmres_tests -- -D warnings`, package
  `rustup run nightly cargo check -p cfd-math --all-targets`, package
  `rustup run nightly cargo clippy -p cfd-math --all-targets -- -D warnings`,
  touched-file `git diff --check`, and a targeted residue scan with no
  `DVector`, nalgebra vector, `Storage`, storage-slice conversion,
  matrix-vector `&a * &x`, or vector `.norm()` residue in the simple GMRES
  test module. Residual: `simple_gmres_tests.rs` still uses the current
  `nalgebra_sparse::CsrMatrix`/`CooMatrix` matrix storage boundary, and broader
  cfd-math integration/adversarial/core/preconditioner tests still contain
  nalgebra `DVector` bridge diagnostics.
- [x] Migrate cfd-math matrix-free solver tests to Leto arrays.
  `crates/cfd-math/src/linear_solver/matrix_free/tests.rs` now constructs
  `leto::Array1` RHS/solution buffers directly for CG and GMRES matrix-free
  identity/scaled-operator tests instead of round-tripping through nalgebra
  `DVector`. The operator-size mismatch test now asserts the exact typed
  `InvalidConfiguration` message at the Leto boundary. Evidence in
  `D:/atlas/repos/CFDrs`: `rustup run nightly cargo fmt -p cfd-math --check`,
  `rustup run nightly cargo check -p cfd-math --all-targets`, `rustup run
  nightly cargo nextest run -p cfd-math matrix_free --status-level fail`
  (4/4), `rustup run nightly cargo clippy -p cfd-math --all-targets -- -D
  warnings`, touched-file `git diff --check`, and a targeted residue scan with
  no `DVector`, nalgebra, `Storage`, storage-slice conversion, or local
  solve-bridge macro residue in the matrix-free test module. Residual:
  cfd-math still has nalgebra `DVector` bridges in broader linear-solver
  integration/adversarial/core/preconditioner test diagnostics, plus the
  shared `nalgebra_sparse::CsrMatrix` and transitional `nalgebra::RealField`
  provider boundaries.
- [x] Migrate cfd-math BiCGSTAB workspaces to Leto arrays.
  `crates/cfd-math/src/linear_solver/bicgstab/mod.rs` now accepts
  `leto::Array1` RHS/solution buffers for preconditioned and unpreconditioned
  BiCGSTAB, owns residual/search/operator/preconditioned workspaces as Leto
  arrays, and applies operators/preconditioners without the legacy nalgebra
  vector bridge. Shared Leto vector helpers now live in
  `crates/cfd-math/src/linear_solver/array_ops.rs` and are reused by CG and
  BiCGSTAB. `LinearSolverChain` now keeps its final BiCGSTAB fallback on the
  Leto RHS/solution path, and the obsolete bridge helpers were removed from
  `linear_solver/traits.rs`. Evidence in `D:/atlas/repos/CFDrs`: `rustup run
  nightly cargo fmt -p cfd-math --check`, `rustup run nightly cargo check -p
  cfd-math --all-targets`, `rustup run nightly cargo nextest run -p cfd-math
  bicgstab --status-level fail` (24/24), `rustup run nightly cargo nextest run
  -p cfd-math linear_solver --status-level fail` (176/176), `rustup run
  nightly cargo nextest run -p cfd-math --test amg_integration_test
  --status-level fail` (5/5), `rustup run nightly cargo clippy -p cfd-math
  --all-targets -- -D warnings`, `git diff --check` for the touched files, and
  targeted scans showing no legacy bridge helper residue plus no `DVector`/
  nalgebra vector-operation residue in the migrated BiCGSTAB/chain/traits
  source. Residual: cfd-math still retains the shared
  `nalgebra_sparse::CsrMatrix` storage/provider boundary, transitional
  `nalgebra::RealField` scalar bounds, and nalgebra `DVector` in remaining
  matrix-free/preconditioner/integration test diagnostics.
- [x] Migrate cfd-math Conjugate Gradient workspaces to Leto arrays.
  `crates/cfd-math/src/linear_solver/conjugate_gradient/mod.rs` now stores the
  reusable CG workspace as `leto::Array1` buffers, runs preconditioned and
  unpreconditioned CG directly on Leto RHS/solution arrays, applies operators
  and preconditioners through the Leto trait boundary, and no longer calls the
  legacy nalgebra-vector bridge helpers. CG benchmarks in
  `crates/cfd-math/benches/{cg_bench,math_benchmarks}.rs` now construct Leto
  vectors at the measured CG API. Co-located tests assert value semantics for
  solved systems and exact dimension/max-iteration errors. Evidence in
  `D:/atlas/repos/CFDrs`: `cargo fmt -p cfd-math --check`, `cargo check -p
  cfd-math --all-targets`, `cargo clippy -p cfd-math --all-targets -- -D
  warnings`, focused `cargo nextest run -p cfd-math conjugate --status-level
  fail` (13/13), and broader `cargo nextest run -p cfd-math linear_solver
  --status-level fail` (176/176) passed. A targeted scan found no `DVector`,
  legacy bridge helper calls, `num_traits`, Leto storage-slice conversion,
  nalgebra vector math calls, or nalgebra `copy_from` residue in the migrated
  CG code and CG benchmark call sites. Residual after the BiCGSTAB follow-up:
  the shared linear-solver trait family still carries the transitional
  `nalgebra::RealField` scalar bound and sparse storage remains on
  `nalgebra_sparse::CsrMatrix`.
- [x] Migrate cfd-math Schwarz preconditioner local apply paths to Leto arrays.
  `crates/cfd-math/src/linear_solver/preconditioners/schwarz.rs` now exposes
  additive and multiplicative Schwarz application as `leto::Array1` methods,
  extracts local RHS buffers directly into Leto arrays, solves local ILU
  problems without nalgebra `DVector` bridges, and applies the default
  `Preconditioner::apply_to` path without converting through nalgebra vectors.
  Residual/output length mismatches now return typed configuration errors, and
  the regression tests assert exact identity-preconditioned values plus exact
  mismatch messages. Evidence in `D:/atlas/repos/CFDrs`: `cargo fmt -p
  cfd-math --check`, `cargo check -p cfd-math --all-targets`, `cargo clippy
  -p cfd-math --all-targets -- -D warnings`, focused `cargo nextest run -p
  cfd-math schwarz --status-level fail` (3/3), and broader `cargo nextest run
  -p cfd-math preconditioner --status-level fail` (76/76) passed. A targeted
  scan found no `DVector`, `Storage`, `num_traits`, `FromPrimitive`,
  local-RHS conversion bridge, `DMatrix`, or `ndarray` residue in
  `schwarz.rs`. Residual: Schwarz still stores/constructs local matrices
  through the shared `nalgebra_sparse::CsrMatrix` boundary and remains under
  the global `Preconditioner<T>` nalgebra scalar bound until the sparse
  preconditioner stack moves fully to Leto/Eunomia storage and scalar traits.
- [x] Migrate cfd-math ILU triangular solve workspaces to Leto arrays.
  `crates/cfd-math/src/linear_solver/preconditioners/ilu/{types,triangular}.rs`
  now runs ILU forward/backward substitution over `leto::Array1` residual,
  intermediate, and solution buffers without constructing nalgebra `DVector`
  workspaces. `IncompleteLU::apply_to` validates residual/output lengths with
  typed configuration errors, and the U-solve diagonal identity routes through
  Eunomia `NumericElement`. `ilu/tests.rs` now covers exact typed
  residual/output length mismatch errors at the Leto boundary. Evidence in
  `D:/atlas/repos/CFDrs`: `cargo fmt -p cfd-math --check`, `cargo check -p
  cfd-math --all-targets`, `cargo clippy -p cfd-math --all-targets -- -D
  warnings`, focused `cargo nextest run -p cfd-math ilu --status-level fail`
  (21/21), and broader `cargo nextest run -p cfd-math preconditioner
  --status-level fail` (74/74) passed. A targeted scan found no `DVector`,
  `DMatrix`, `Storage`, old scalar identities, fallback wording,
  `component_mul`, `rows_mut`, row-view residue, `num_traits`, or
  `num_complex` in `ilu/types.rs` or `ilu/triangular.rs`. Residual:
  IncompleteLU still stores the shared `nalgebra_sparse::CsrMatrix` LU factor
  boundary and participates in the global `Preconditioner<T>` nalgebra scalar
  bound until the wider sparse/preconditioner stack moves fully to Leto and
  Eunomia.
- [x] Migrate cfd-math deflation preconditioner vector state to Leto arrays.
  `crates/cfd-math/src/linear_solver/preconditioners/deflation.rs` now stores
  deflation eigenvectors as `leto::Array1`, accepts added eigenpairs at the
  Leto boundary, computes projection coefficients without nalgebra `DVector`
  workspaces, validates residual/output/eigenvector lengths, and rejects zero
  eigenvalues before division. `preconditioners/edge_case_tests.rs` now covers
  Leto-native deflation projection, typed length-mismatch errors, and
  zero-eigenvalue rejection. Evidence in `D:/atlas/repos/CFDrs`: `cargo fmt
  -p cfd-math --check`, `cargo check -p cfd-math --all-targets`, `cargo
  clippy -p cfd-math --all-targets -- -D warnings`, focused
  `cargo nextest run -p cfd-math deflation --status-level fail` (3/3), and
  broader `cargo nextest run -p cfd-math preconditioner --status-level fail`
  (73/73) passed. A targeted scan found no `DVector`, `DMatrix`, `Storage`,
  old scalar identities, fallback wording, `component_mul`, `rows_mut`,
  row-view residue, `num_traits`, or `num_complex` in `deflation.rs`.
  Residual: `DeflationPreconditioner` still wraps the base preconditioner as
  the existing `Box<dyn Preconditioner<T>>`, and the global
  `Preconditioner<T>` trait remains constrained by nalgebra `RealField` until
  the broader linear-solver trait family moves fully to Eunomia.
- [x] Migrate cfd-math basic preconditioner vector internals to Leto arrays.
  `crates/cfd-math/src/linear_solver/preconditioners/basic.rs` now stores the
  Jacobi inverse diagonal as `leto::Array1`, applies Identity/Jacobi/SOR at
  the Leto residual/output boundary, validates matrix-sized Jacobi/SOR inputs
  and all outputs, and routes scalar identities and diagonal tolerance through
  Eunomia. `linear_solver/tests/mod.rs` now has a value-semantic regression
  that matches the exact typed length-mismatch errors for Identity, Jacobi,
  and SOR. Evidence in `D:/atlas/repos/CFDrs`: `cargo fmt -p cfd-math
  --check`, `cargo check -p cfd-math --all-targets`, `cargo clippy -p
  cfd-math --all-targets -- -D warnings`, focused nextest for
  `basic_preconditioners_reject_mismatched_leto_vector_lengths` (1/1), and
  `cargo nextest run -p cfd-math preconditioner --status-level fail` (70/70)
  passed. A targeted scan found no `DVector`, `DMatrix`, old scalar
  identities, fallback wording, `component_mul`, `rows_mut`, or row-view
  residue in `basic.rs`. Residual: basic matrix-backed preconditioners still
  store the shared `nalgebra_sparse::CsrMatrix` boundary and remain
  constrained by the global `Preconditioner<T>` nalgebra scalar bound until
  the broader preconditioner trait family moves fully to Eunomia/Leto sparse
  storage.
- [x] Migrate cfd-math IncompleteCholesky preconditioner solve workspaces to
  Leto arrays. `crates/cfd-math/src/linear_solver/preconditioners/cholesky.rs`
  now performs forward/backward triangular substitution over `leto::Array1`
  buffers and applies into the caller-provided Leto output array without
  allocating nalgebra `DVector` residual, intermediate, or solution
  workspaces. Cholesky residual/output shape mismatches now return typed
  configuration errors before substitution. The IC(0) factorization square
  root dispatch is routed through Eunomia `NumericElement`. Evidence in
  `D:/atlas/repos/CFDrs`: `cargo fmt -p cfd-math --check`, `cargo check -p
  cfd-math --all-targets`, `cargo clippy -p cfd-math --all-targets -- -D
  warnings`, and `cargo nextest run -p cfd-math cholesky --status-level fail`
  (5/5) passed. A targeted scan found no `DVector`, `DMatrix`, `Storage`, old
  scalar identities, silent fallback wording, ambiguous `.sqrt()`, or residual
  `DVector` workspace construction in `cholesky.rs`. Follow-up:
  Sprint 1.96.166 moved IncompleteCholesky factor storage to Leto CSR.
  Residual sparse-provider work is now Schwarz, direct solver, and the shared
  solver matrix boundary.
- [x] Migrate cfd-math SSOR preconditioner sweep workspaces to Leto arrays.
  `crates/cfd-math/src/linear_solver/preconditioners/ssor.rs` now runs
  forward/backward SOR sweeps directly over `leto::Array1` buffers and applies
  into the caller-provided Leto output array without allocating nalgebra
  `DVector` residual or solution workspaces. SSOR residual/output shape
  mismatches now return typed configuration errors before sweeping. Evidence
  in `D:/atlas/repos/CFDrs`: `cargo fmt -p cfd-math --check`, `cargo check
  -p cfd-math --all-targets`, `cargo clippy -p cfd-math --all-targets -- -D
  warnings`, and `cargo nextest run -p cfd-math ssor --status-level fail`
  (5/5) passed. A targeted scan found no `DVector`, `DMatrix`, `Storage`,
  nalgebra row-view operations, `component_mul`, old scalar identities, silent
  fallback wording, or clone fallback residue in `ssor.rs`. Residual: SSOR
  still stores the current shared `nalgebra_sparse::CsrMatrix` matrix boundary
  and participates in the global `Preconditioner<T>` nalgebra scalar bound
  until the broader preconditioner trait family moves fully to Eunomia.
- [x] Migrate cfd-math block/SIMPLE preconditioner vector internals to Leto
  arrays. `crates/cfd-math/src/linear_solver/block_preconditioner.rs` now
  stores Jacobi diagonal inverses and SIMPLE Schur diagonal state as
  `leto::Array1`, exposes direct `apply` methods over Leto arrays, and applies
  the `Preconditioner::apply_to` trait without constructing nalgebra
  `DVector` bridges. Size mismatches now return typed configuration errors
  instead of cloning the input vector. Evidence in `D:/atlas/repos/CFDrs`:
  `cargo fmt -p cfd-math --check`, `cargo check -p cfd-math --all-targets`,
  `cargo clippy -p cfd-math --all-targets -- -D warnings`, and `cargo nextest
  run -p cfd-math block_preconditioner --status-level fail` (4/4) passed.
  A targeted source scan found no `DVector`, `DMatrix`, nalgebra import,
  nalgebra row-view operations, `component_mul`, old scalar identities, silent
  fallback wording, or clone fallback in `block_preconditioner.rs`. Residual:
  the `Preconditioner` trait itself still carries the broader
  `nalgebra::RealField` scalar bound, and sparse matrix storage still aliases
  `nalgebra_sparse::CsrMatrix`.
- [x] Migrate cfd-math GMRES workspace and solver-chain GMRES tiers to Leto
  arrays. `crates/cfd-math/src/linear_solver/gmres/{solver,arnoldi}.rs` now
  keeps Arnoldi basis work vectors, preconditioned work vectors, and residual
  checks in `leto::Array1` instead of nalgebra `DVector`, and GMRES
  `solve_preconditioned`/`solve_unpreconditioned` now accept Leto RHS and
  solution arrays directly. `LinearSolverChain` now routes its GMRES+AMG,
  GMRES+block, unpreconditioned GMRES, and GMRES+ILU tiers through the same
  Leto arrays instead of converting to nalgebra vectors first. Evidence in
  `D:/atlas/repos/CFDrs`: `cargo check -p cfd-math --lib`, `cargo check -p
  cfd-math --tests`, `cargo check -p cfd-math --all-targets`, `cargo fmt -p
  cfd-math --check`, `cargo clippy -p cfd-math --all-targets -- -D warnings`,
  `cargo nextest run -p cfd-math gmres --status-level fail` (21/21), `cargo
  check -p cfd-2d --no-default-features --lib`, and `cargo nextest run -p
  cfd-2d --no-default-features momentum --status-level fail` (53/53) passed.
  Targeted scans found no `DVector` or legacy operator/preconditioner bridge
  calls under `linear_solver/gmres`. Residual after the CG/BiCGSTAB follow-up:
  cfd-math still has nalgebra scalar bounds and `nalgebra_sparse` storage
  pending later Leto/Eunomia slices.
- [x] Migrate cfd-2d momentum vectors and scalar seam off direct nalgebra.
  `crates/cfd-2d/src/physics/momentum/**` now assembles, stores, solves, and
  applies momentum RHS/solution vectors as `leto::Array1` instead of nalgebra
  `DVector`. Momentum boundary helpers mutate Leto RHS buffers directly,
  momentum GMRES/direct solve paths call the Leto-native solver APIs, and the
  obsolete `crates/cfd-2d/src/linear_solver_bridge.rs` module was removed.
  `Cfd2dScalar` no longer requires `nalgebra::RealField`, and `cfd-2d` no
  longer declares direct `nalgebra`/`nalgebra-sparse` dependencies. Evidence in
  `D:/atlas/repos/CFDrs`: `cargo fmt -p cfd-2d --check`, `cargo check -p
  cfd-2d --no-default-features --lib`, `cargo clippy -p cfd-2d
  --no-default-features --all-targets -- -D warnings`, `cargo nextest run -p
  cfd-2d --no-default-features momentum --status-level fail` (53/53), and
  targeted cfd-2d nalgebra/DVector/DMatrix/bridge source/manifest scans passed.
  Residual: cfd-2d still resolves nalgebra/nalgebra-sparse transitively through
  upstream crates including cfd-1d, cfd-core, cfd-math, cfd-schematics, and
  Gaia while those owners continue migrating to Leto/Eunomia/Gaia provider
  boundaries.
- [x] Migrate cfd-2d pressure-velocity correction vectors to Leto arrays.
  `crates/cfd-2d/src/pressure_velocity/{pressure,correction,faces}.rs` now
  cache, assemble, solve, and scatter pressure-correction RHS/solution vectors
  as `leto::Array1` instead of nalgebra `DVector`. The pressure-velocity solve
  dispatch now calls the Leto-native iterative solver trait and direct sparse
  fallback directly, and RHS diagnostics use `leto_ops::norm_l2`. Evidence in
  `D:/atlas/repos/CFDrs`: `cargo fmt -p cfd-2d --check`, `cargo check -p
  cfd-2d --no-default-features --lib`, `cargo clippy -p cfd-2d
  --no-default-features --all-targets -- -D warnings`, `cargo nextest run -p
  cfd-2d --no-default-features pressure_velocity --status-level fail` (16/16),
  and a targeted `pressure_velocity` `DVector`/nalgebra/bridge residue scan
  passed. Residual: direct cfd-2d source/manifest nalgebra ownership is now
  removed; nalgebra still resolves transitively through upstream crates while
  those owners continue migrating to Leto/Eunomia/Gaia.
- [x] Migrate cfd-2d SIMPLE pressure-correction vectors to Leto arrays.
  `crates/cfd-2d/src/solvers/simple/{algorithm,pressure}.rs` now stores and
  solves pressure-correction RHS and `p_prime` as `leto::Array1` instead of
  nalgebra `DVector`. The SIMPLE pressure solve now calls the Leto-native
  `IterativeLinearSolver::solve` boundary directly instead of routing through
  the local nalgebra conversion bridge. Evidence in `D:/atlas/repos/CFDrs`:
  `cargo fmt -p cfd-2d --check`, `cargo check -p cfd-2d --no-default-features
  --lib`, `cargo clippy -p cfd-2d --no-default-features --all-targets -- -D
  warnings`, `cargo nextest run -p cfd-2d --no-default-features simple
  --status-level fail` (19/19), and a targeted `solvers/simple`
  `DVector`/nalgebra/bridge residue scan passed. Residual: cfd-2d still
  has no direct cfd-2d source/manifest nalgebra ownership; nalgebra still
  resolves transitively through upstream crates while those owners continue
  migrating to Leto/Eunomia/Gaia.
- [x] Migrate cfd-2d FDM RHS and Gauss-Seidel vectors to Leto arrays.
  `crates/cfd-2d/src/solvers/fdm/{linear_solver,poisson,advection_diffusion}.rs`
  now use `leto::Array1` for RHS and solution vectors instead of nalgebra
  `DVector`. Poisson and advection-diffusion stencil assembly mutate Leto
  arrays directly, and `solve_gauss_seidel` returns a Leto solution array.
  Evidence in `D:/atlas/repos/CFDrs`: `cargo fmt -p cfd-2d --check`,
  `cargo check -p cfd-2d --no-default-features --lib`, `cargo clippy -p
  cfd-2d --no-default-features --all-targets -- -D warnings`, `cargo nextest
  run -p cfd-2d --no-default-features fdm --status-level fail` (2/2), and a
  targeted `solvers/fdm` `DVector`/nalgebra residue scan passed. Residual:
  direct cfd-2d source/manifest nalgebra ownership is now removed; nalgebra
  still resolves transitively through upstream crates while those owners
  continue migrating to Leto/Eunomia/Gaia.
- [x] Migrate cfd-2d time-integration vectors to Leto arrays.
  `crates/cfd-2d/src/schemes/time/**` now exposes and tests time-step state as
  `leto::Array1` through the `StateVector<T>` domain alias instead of nalgebra
  `DVector`. Explicit, implicit, multistep, adaptive-controller, adaptive
  integrator, and time tests now use provider-native Leto vector arithmetic,
  indexing, shape, and `leto_ops::norm_l2` convergence checks. Evidence in
  `D:/atlas/repos/CFDrs`: `cargo fmt -p cfd-2d --check`, `cargo check -p
  cfd-2d --no-default-features --lib`, `cargo clippy -p cfd-2d
  --no-default-features --all-targets -- -D warnings`, `cargo nextest run -p
  cfd-2d --no-default-features time --status-level fail` (29/29), and
  targeted `schemes/time` `DVector`/nalgebra residue scans passed. Residual:
  direct cfd-2d source/manifest nalgebra ownership is now removed; nalgebra
  still resolves transitively through upstream crates while those owners
  continue migrating to Leto/Eunomia/Gaia.
- [x] Migrate compact cfd-2d `DMatrix` residue to Leto arrays.
  `crates/cfd-2d/src/physics/immersed_boundary.rs`,
  `crates/cfd-2d/src/schemes/grid.rs`, the dependent scheme callers/tests, and
  `crates/cfd-2d/examples/blood_venturi.rs` now use `leto::Array2` shape and
  indexing instead of nalgebra `DMatrix`. Evidence in `D:/atlas/repos/CFDrs`:
  `cargo fmt -p cfd-2d --check`, `cargo check -p cfd-2d
  --no-default-features --lib`, `cargo check -p cfd-2d --no-default-features
  --example blood_venturi`, `cargo clippy -p cfd-2d --no-default-features --lib
  -- -D warnings`, and `cargo nextest run -p cfd-2d --no-default-features
  immersed_boundary schemes upwind muscl --status-level fail` (60/60) passed.
  Targeted scans found no cfd-2d source/example/test/bench `DMatrix` residue and
  no nalgebra-style `Grid2D.data` tuple access. Residual: cfd-2d still retains
  nalgebra `DVector`/sparse linear-system boundaries and other non-DMatrix
  nalgebra use while those provider seams continue migrating to Leto/Eunomia.
- [x] Migrate cfd-1d vascular Bessel/Womersley complex math to Eunomia.
  `crates/cfd-1d/src/physics/vascular/bessel.rs` and the Womersley profile
  now use `eunomia::Complex` instead of nalgebra `Complex`/`ComplexField`;
  Bessel convergence uses Eunomia's native complex `norm()`; and a targeted
  vascular residue scan found no remaining nalgebra complex imports or
  `ComplexField::modulus()` calls in the Womersley path. Evidence in
  `D:/atlas/repos/CFDrs`: `cargo check -p cfd-1d --no-default-features
  --lib`, `cargo fmt -p cfd-1d --check`, `cargo clippy -p cfd-1d
  --no-default-features --lib -- -D warnings`, and `cargo nextest run -p
  cfd-1d --no-default-features bessel womersley --status-level fail`
  (26/26) passed. Residual: `cfd-1d` still depends on nalgebra for its
  network sparse/dense linear-system boundary and scalar seam while that
  storage layer continues to migrate to Leto/Eunomia.
- [x] Migrate `cfd-math::linear_solver::LinearOperator::apply` to Leto
  vectors. The public operator trait now accepts `leto::Array1<T>` input and
  output buffers, and `apply_transpose` follows the same boundary. Sparse CSR,
  identity/scaled, Poisson, momentum, energy, and GPU operator adapters now
  implement the Leto vector API; CG and BiCGSTAB use an internal legacy bridge
  for their current nalgebra workspaces; and the Laplacian
  benchmark now calls the Leto operator boundary. Evidence in
  `D:/atlas/repos/CFDrs`: `cargo check -p cfd-math --no-default-features
  --all-targets`, `cargo fmt -p cfd-math --check`, `cargo clippy -p cfd-math
  --no-default-features --all-targets -- -D warnings`, `cargo nextest run -p
  cfd-math --no-default-features conjugate_gradient bicgstab gmres operator
  poisson momentum matrix_free core_solver simple_gmres --status-level fail`
  (80/80), and a targeted DVector operator-signature residue scan passed.
  Residual after the CG/BiCGSTAB follow-up: nalgebra sparse storage and some
  preconditioner internals still retain local `DVector`/`CsrMatrix` conversion
  bridges; cfd-validation result/error storage still uses `DVector`; broad
  cfd-validation nextest still has the existing venturi cross-fidelity
  convergence failures.
- [x] Migrate `cfd-math::linear_solver::Preconditioner::apply_to` to Leto
  vectors. The public preconditioner trait now accepts `leto::Array1<T>`
  residual and output buffers; cfd-math concrete preconditioners, AMG,
  Schwarz local solves, deflation, and preconditioner tests call the Leto
  boundary; and cfd-1d network Jacobi preconditioning now implements the same
  Leto public trait. Evidence in `D:/atlas/repos/CFDrs`: `cargo check -p
  cfd-math --no-default-features --all-targets`, `cargo check -p cfd-1d
  --no-default-features --lib`, `cargo fmt -p cfd-math -p cfd-1d --check`,
  `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings`, `cargo clippy -p cfd-1d --no-default-features --lib -- -D
  warnings`, `cargo nextest run -p cfd-math --no-default-features
  conjugate_gradient bicgstab gmres preconditioner ilu matrix_free
  core_solver simple_gmres --status-level fail` (131/131), and a targeted
  `apply_to` DVector signature residue scan passed. Residual:
  `LinearOperator::apply`, preconditioner internals backed by nalgebra sparse
  matrices, and iterative solver workspaces still retain nalgebra
  `DVector`/`CsrMatrix` conversion bridges; cfd-validation result/error
  storage still uses `DVector`; broad cfd-validation nextest still has the
  existing venturi cross-fidelity convergence failures.
- [x] Migrate `cfd-math::linear_solver::IterativeLinearSolver::solve` to Leto
  vectors and update downstream crate call sites. The public iterative solver
  trait now accepts `leto::Array1<T>` RHS/result buffers; CG, BiCGSTAB, and
  GMRES bridge only inside their current nalgebra workspaces; cfd-math tests
  call the Leto boundary; cfd-1d network solves, cfd-2d momentum/pressure
  solves, and cfd-3d FEM projection solve through local Leto bridges. Evidence
  in `D:/atlas/repos/CFDrs`: `cargo fmt -p cfd-math --check`, `cargo check
  -p cfd-math --no-default-features --all-targets`, `cargo clippy -p
  cfd-math --no-default-features --all-targets -- -D warnings`, `cargo
  nextest run -p cfd-math --no-default-features conjugate_gradient bicgstab
  gmres matrix_free core_solver simple_gmres --status-level fail` (61/61),
  `cargo fmt -p cfd-math -p cfd-1d -p cfd-2d -p cfd-3d --check`, `cargo check
  -p cfd-1d --no-default-features --lib`, `cargo check -p cfd-2d
  --no-default-features --lib`, `cargo check -p cfd-3d --no-default-features
  --lib`, `cargo clippy -p cfd-1d --no-default-features --lib -- -D
  warnings`, `cargo clippy -p cfd-2d --no-default-features --all-targets --
  -D warnings`, `cargo clippy -p cfd-3d --no-default-features --lib -- -D
  warnings`, `cargo clippy -p cfd-validation --no-default-features
  --all-targets -- -D warnings`, targeted `solve` DVector call-site residue
  scan, and `git diff --check` passed. Residual: `LinearOperator::apply`,
  `Preconditioner::apply_to`, preconditioners, and internal solver workspaces
  still expose nalgebra `DVector`/`CsrMatrix`; cfd-validation numerical result
  storage still uses `DVector`; broad cfd-validation nextest still has the
  existing venturi cross-fidelity convergence failures.
- [x] Migrate `cfd-math::linear_solver::LinearSolver::solve_system` to Leto
  vectors. The public trait and the CG, BiCGSTAB, and GMRES implementations
  now accept `leto::Array1<T>` RHS/initial-guess vectors and return
  `leto::Array1<T>` results; the iterative algorithms bridge internally to
  the current nalgebra workspaces until the remaining operator and
  preconditioner traits move; and cfd-validation numerical solver validation
  calls the new Leto public API while retaining nalgebra only for existing
  validation error metrics. Evidence in `D:/atlas/repos/CFDrs`: `cargo fmt -p
  cfd-math -p cfd-validation --check`, `cargo check -p cfd-math
  --no-default-features --lib`, `cargo check -p cfd-validation
  --no-default-features --lib`, `cargo clippy -p cfd-math
  --no-default-features --all-targets -- -D warnings`, `cargo clippy -p
  cfd-validation --no-default-features --lib -- -D warnings`, `cargo clippy
  -p cfd-validation --no-default-features --all-targets -- -D warnings`,
  `cargo clippy -p cfd-2d --no-default-features --all-targets -- -D
  warnings`, `cargo nextest run -p cfd-math --no-default-features
  conjugate_gradient bicgstab gmres --status-level fail` (58/58), targeted
  `solve_system` signature/residue scans, and `git diff --check` passed.
  Residual: `IterativeLinearSolver`, `LinearOperator`, `Preconditioner`,
  preconditioners, and internal iterative workspaces still expose nalgebra
  `DVector`/`CsrMatrix`; cfd-validation result/error metric storage still
  uses `DVector`; and broad cfd-validation nextest still has the existing
  venturi cross-fidelity convergence failures recorded below.
- [x] Migrate cfd-validation SpMV benchmark callers and validation scalar
  bounds to the Leto/Eunomia provider seams. `profile_matrix_operations`,
  `profile_algorithm_performance`, and `benchmark_spmv` now call
  `cfd_math::sparse::spmv` with `leto::Array1` instead of nalgebra
  `DVector`; `LinearSolverValidator` uses the crate-local
  `ValidationScalar` provider contract for sparse linear-operator calls; and
  1D blood-flow literature validations use the same seam for cfd-1d network
  solver dispatch. Evidence in `D:/atlas/repos/CFDrs`: `cargo fmt -p
  cfd-validation --check`, `cargo check -p cfd-validation
  --no-default-features --lib`, `cargo clippy -p cfd-validation
  --no-default-features --lib -- -D warnings`, `cargo clippy -p
  cfd-validation --no-default-features --all-targets -- -D warnings`,
  `cargo clippy -p cfd-2d --no-default-features --all-targets -- -D
  warnings`, targeted SpMV DVector-residue scan over cfd-validation, `cargo
  nextest run -p cfd-validation --no-default-features benchmark
  --status-level fail` (40/40), and `git diff --check` passed. The broad
  `cargo nextest run -p cfd-validation --no-default-features --status-level
  fail` built successfully but failed in two existing venturi cross-fidelity
  convergence tests unrelated to this SpMV/scalar-bound slice:
  `microventuri_35um_case_produces_converged_informative_2d_result` and
  `option2_selected_45um_geometry_routes_to_fallback_and_converges`.
  Residual: cfd-validation still has the venturi 2D fallback convergence
  baseline failures, and cfd-math public solver/preconditioner traits still
  expose nalgebra `DVector`/`CsrMatrix` boundaries.
- [x] Migrate the `LinearSolverChain` public vector boundary and downstream
  cfd-2d/cfd-3d direct assembly consumers to Leto arrays. The chain now
  accepts and returns `leto::Array1<T>` for `solve` and `solve_with_guess`;
  cfd-3d FEM converts at a private FEM Leto bridge before chain dispatch and
  before `SparseMatrixBuilder::build_with_rhs`; cfd-2d momentum/pressure
  direct fallbacks share a private direct-solver bridge; cfd-2d and cfd-3d
  crate scalar seams now carry the Leto real-scalar provider bound; and
  cfd-1d row-equilibration/scalar bounds were tightened so 2D network coupling
  compiles through the migrated sparse provider APIs. Evidence in
  `D:/atlas/repos/CFDrs`: `cargo fmt -p cfd-math -p cfd-1d -p cfd-2d -p
  cfd-3d --check`, `cargo check -p cfd-math --no-default-features --lib`,
  `cargo check -p cfd-1d --no-default-features --lib`, `cargo check -p
  cfd-2d --no-default-features --lib`, `cargo check -p cfd-3d
  --no-default-features --lib`, `cargo nextest run -p cfd-math
  --no-default-features chain direct_solver core_solver simple_gmres
  --status-level fail` (4/4), `cargo clippy -p cfd-math
  --no-default-features --all-targets -- -D warnings`, `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings`, and `cargo clippy -p
  cfd-3d --no-default-features --lib -- -D warnings`. Residual:
  cfd-validation still blocks `cargo clippy -p cfd-2d --no-default-features
  --all-targets -- -D warnings` because it calls public Leto SpMV with
  nalgebra vectors and its generic validation/1D-reference paths lack the
  propagated Leto scalar bounds.
- [x] Migrate `cfd-math::linear_solver::DirectSparseSolver` to Leto RHS/result
  vectors. The direct sparse solver public `solve` API now accepts
  `leto::Array1<T>` and returns `leto::Array1<T>`; its rsparse RHS conversion,
  solution construction, finite checks, and Leto dense fallback all operate on
  Leto arrays; the obsolete DVector dense-fallback wrapper was removed; and
  `LinearSolverChain` performs the only remaining conversion at its still
  nalgebra-facing boundary. Evidence in `D:/atlas/repos/CFDrs`: `cargo fmt -p
  cfd-math --check`, `cargo check -p cfd-math --no-default-features --lib`,
  `cargo nextest run -p cfd-math --no-default-features direct_solver chain
  core_solver simple_gmres --status-level fail` (4/4), `cargo clippy -p
  cfd-math --no-default-features --all-targets -- -D warnings`, `cargo doc -p
  cfd-math --no-default-features --no-deps`, `cargo test --doc -p cfd-math
  --no-default-features` (3 passed, 3 ignored), targeted direct-solver
  DVector-signature residue scan, and `git diff --check`. Residual:
  `LinearSolverChain`, iterative solver traits, preconditioners, and sparse
  storage still expose nalgebra `DVector`/`CsrMatrix` boundaries.
- [x] Migrate `cfd-math::sparse::SparseMatrixBuilder::build_with_rhs` to Leto
  RHS vectors. The Dirichlet column-elimination assembly API now accepts
  `leto::Array1<T>` and mutates that provider-owned RHS directly; `build()`
  no longer fabricates a dummy nalgebra vector; sparse tests include a
  value-semantic Dirichlet fixture proving the Leto RHS receives the eliminated
  column contribution; and the direct/block-preconditioner tests use Leto for
  matrix assembly while keeping nalgebra `DVector` only at the still-current
  solver boundary. Evidence in `D:/atlas/repos/CFDrs`: `cargo fmt -p
  cfd-math --check`, `cargo check -p cfd-math --no-default-features --lib`,
  `cargo nextest run -p cfd-math --no-default-features sparse direct_solver
  block_preconditioner --status-level fail` (25/25), `cargo clippy -p
  cfd-math --no-default-features --all-targets -- -D warnings`, targeted
  `build_with_rhs` residue scan showing `Array1`, and `git diff --check`.
  Residual: the public sparse storage alias and the linear-solver/
  preconditioner traits still expose `nalgebra_sparse::CsrMatrix` and
  nalgebra `DVector`.
- [x] Migrate public `cfd-math::sparse` SpMV wrappers to Leto vectors.
  The then-public `spmv`, `spmv_parallel`, and `try_spmv` accepted
  `leto::Array1` input/output vectors; the compatibility entry point was
  subsequently removed by the 2026-07-10 provider-ownership closure above.
  The remaining nalgebra `DVector` SpMV bridge is private to the
  `LinearOperator for CsrMatrix` implementation because that trait boundary is
  still nalgebra-based. Sparse tests, GMRES/AMG integration tests, interpolation
  quality checks, and the SpMV benchmark now call the public SpMV API with Leto
  arrays. Evidence in `D:/atlas/repos/CFDrs`: `cargo fmt -p cfd-math
  --check`, `cargo check -p cfd-math --no-default-features --lib`, `cargo
  nextest run -p cfd-math --no-default-features sparse spmv interpolation amg
  simple_gmres core_solver --status-level fail` (40/40), `cargo clippy -p
  cfd-math --no-default-features --all-targets -- -D warnings`, public SpMV
  signature scan showing `Array1` and only the private `try_spmv_dvector`
  retains `DVector`, and `git diff --check`. Residual: the public
  `LinearOperator`/solver/preconditioner trait family still exposes nalgebra
  `DVector`, and sparse storage still aliases `nalgebra_sparse::CsrMatrix`.
- [x] Migrate `cfd-math::sparse::SparseMatrixExt` diagonal and scaling
  vectors to Leto arrays. `SparseMatrixExt::diagonal` now returns
  `leto::Array1`, and `set_diagonal`, `scale_rows`, and `scale_columns` now
  accept `leto::Array1` instead of nalgebra `DVector`. The Jacobi
  preconditioner consumes the provider-backed Leto diagonal through
  `Storage::as_slice`, while retaining nalgebra `DVector` only for the current
  public `Preconditioner` trait boundary. Evidence in `D:/atlas/repos/CFDrs`:
  `cargo fmt -p cfd-math --check`, `cargo check -p cfd-math
  --no-default-features --lib`, `cargo nextest run -p cfd-math
  --no-default-features sparse basic --status-level fail` (21/21), `cargo
  clippy -p cfd-math --no-default-features --all-targets -- -D warnings`,
  targeted residue scan showing no `DVector` signatures for
  `SparseMatrixExt` diagonal/scaling methods, and `git diff --check` for the
  touched sparse/basic files. Residual: the linear-solver/preconditioner trait
  boundary still exposes nalgebra `DVector`, and sparse storage still aliases
  `nalgebra_sparse::CsrMatrix`.
- [x] Migrate `cfd-math::linear_solver::preconditioners::multigrid`
  smoother and cycle vector paths to Leto. `MultigridLevel`,
  `AMGHierarchy`, and `MultigridSmoother` now use `leto::Array1` through the
  `MultigridVector` alias; Jacobi, Gauss-Seidel, symmetric Gauss-Seidel,
  SSOR, and Chebyshev smoothers compute residuals through the new
  Leto-array `sparse::spmv_array` bridge; V/W/F cycles and coarsest solves
  now consume/return Leto vectors; and AMG V-cycle recursion bridges to
  nalgebra `DVector` only at the current public `Preconditioner::apply_to`
  boundary. Evidence in `D:/atlas/repos/CFDrs`: `cargo fmt -p cfd-math
  --check`, `cargo check -p cfd-math --no-default-features --lib`,
  `cargo nextest run -p cfd-math --no-default-features multigrid::cycles
  smoothers --status-level fail` (10/10), `cargo clippy -p cfd-math
  --no-default-features --all-targets -- -D warnings`, clean provider-residue
  scan over `smoothers.rs`/`cycles.rs`, and `git diff --check` for the
  migrated sparse/multigrid files. Residual: public sparse and linear-solver
  APIs still expose `nalgebra_sparse::CsrMatrix`, nalgebra `DVector`, and
  nalgebra `RealField`; AMG keeps a bridge at `Preconditioner::apply_to` until
  that trait boundary moves to Leto/Eunomia.
- [x] Migrate `cfd-math::linear_solver::preconditioners::multigrid::gmg`
  to Leto/Eunomia. `GeometricMultigrid`, `NonlinearOperator`, FAS/linear solve
  vectors, Poisson hierarchy matrices, transfer operators, Jacobi relaxation,
  residual computation, and GMG tests now use `leto::Array2`/`Array1` and
  Eunomia `RealField`/`NumericElement` instead of nalgebra `DMatrix`/`DVector`
  and `num_traits` conversions. Local helpers provide provider-backed matvec,
  residual, vector add/subtract, and L2 norm operations without preserving a
  nalgebra compatibility layer. Evidence in `D:/atlas/repos/CFDrs`: `cargo fmt
  -p cfd-math --check`, `cargo check -p cfd-math --no-default-features
  --lib`, `cargo nextest run -p cfd-math --no-default-features gmg
  --status-level fail` (5/5), `cargo clippy -p cfd-math
  --no-default-features --all-targets -- -D warnings`, clean GMG
  provider-residue scan, and `git diff --check`. Residual: the public
  linear-solver/preconditioner boundary still exposes nalgebra
  `DVector`/`RealField`.
- [x] Migrate `cfd-math::linear_solver::gmres` internal dense workspace to
  Leto/Eunomia. GMRES Krylov basis, Hessenberg matrix, Givens rotation
  coefficients, least-squares RHS, and triangular-solve result now use
  `leto::Array2`/`Array1` instead of nalgebra `DMatrix`/internal `DVector`
  storage. Givens math uses Eunomia `RealField`/`NumericElement`, and the
  tiered `LinearSolverChain` carries the Eunomia real-field bound needed by
  the migrated GMRES constructor. Evidence in `D:/atlas/repos/CFDrs`:
  `cargo fmt -p cfd-math --check`, `cargo check -p cfd-math
  --no-default-features --lib`, `cargo nextest run -p cfd-math
  --no-default-features gmres --status-level fail` (21/21), `cargo clippy
  -p cfd-math --no-default-features --all-targets -- -D warnings`, clean
  GMRES provider-residue scan for `DMatrix`/`ndarray`/`num_traits`/
  `num_complex`/`rayon`/`tokio`/`rustfft`/`wgpu`, and `git diff --check`.
  Residual: GMRES still crosses nalgebra `DVector` and `nalgebra::RealField`
  at the current public `LinearOperator`/`Preconditioner`/`LinearSolver`
  trait boundary.
- [x] Migrate `cfd-math::linear_solver::preconditioners::multigrid::restriction`
  dense transfer utilities to Leto arrays. Restriction construction,
  validation, vector restriction, and Galerkin projection now expose/consume
  `leto::Array2`/`Array1`; `restrict_matrix` delegates the dense
  `R * A * P` products to `leto_ops::MatrixProduct` instead of nalgebra
  chained multiplication. The module's tests now assert the exact `P^T v` and
  `P^T A P` fixture values. Evidence in `D:/atlas/repos/CFDrs`: `cargo fmt
  -p cfd-math --check`, `cargo check -p cfd-math --no-default-features
  --lib`, `cargo nextest run -p cfd-math --no-default-features restriction
  --status-level fail` (7/7), `cargo clippy -p cfd-math
  --no-default-features --all-targets -- -D warnings`, and a clean provider
  residue scan over `restriction.rs`. Residual: broader sparse/linear-solver
  APIs still expose nalgebra `DVector` and `nalgebra_sparse::CsrMatrix`.
- [x] Consolidate `cfd-math::linear_solver` legacy CSR-to-dense solves behind
  a Leto provider bridge. `linear_solver::dense_bridge` now builds row-major
  `leto::Array2`/`Array1` values from the current
  `nalgebra_sparse::CsrMatrix`/nalgebra `DVector` API and solves through
  `leto_ops::solve`. `DirectSparseSolver` dense fallback and multigrid cycle
  coarsest small-system solves both consume that bridge; the old local
  Gaussian-elimination helper in multigrid cycles was removed.
  `LinearSolverChain` carries the Leto real-scalar bound required by the
  direct-solver path. Evidence in `D:/atlas/repos/CFDrs`: `cargo fmt -p
  cfd-math --check`, `cargo check -p cfd-math --no-default-features --lib`,
  `cargo nextest run -p cfd-math --no-default-features multigrid::cycles
  direct_solver --status-level fail` (9/9), and `cargo clippy -p cfd-math
  --no-default-features --all-targets -- -D warnings` pass. Residual: the
  public sparse/linear-solver API still exposes
  `nalgebra_sparse::CsrMatrix`/`DVector`, and the primary sparse direct solve
  path still uses `rsparse` until a Leto-owned sparse direct-solver provider
  replaces it.
- [x] Centralize `cfd-math::sparse` CSR bridging and route sparse builder
  construction through Leto validation. Added a single `sparse::bridge` module
  for nalgebra/Leto CSR conversion and DVector-to-Leto view creation;
  `SparseMatrixBuilder::{build,build_with_rhs,build_parallel}` and
  `ParallelAssembly::block_diagonal` now construct/validate CSR storage through
  `leto_ops::CsrMatrix` before converting back at the legacy
  `nalgebra_sparse::CsrMatrix` boundary. Sparse assembly, patterns, and
  Schwarz subdomain extraction now carry the Leto scalar bound required by that
  provider path. Evidence in `D:/atlas/repos/CFDrs`: `cargo check -p
  cfd-math --no-default-features --lib`, `cargo nextest run -p cfd-math
  --no-default-features sparse --status-level fail` (18/18), and `cargo
  clippy -p cfd-math --no-default-features --all-targets -- -D warnings`
  pass. Residual: `cfd-math::sparse::SparseMatrix` still re-exports
  `nalgebra_sparse::CsrMatrix`; replacing that public type remains a larger
  linear-solver API migration.
- [x] Move remaining `cfd-math::sparse::SparseMatrixExt` value operations to
  the Leto CSR provider. `diagonal`, scalar scaling, row scaling, column
  scaling, Frobenius norm, diagonal dominance, and condition-estimate now
  bridge the current `nalgebra_sparse::CsrMatrix` boundary into
  `leto_ops::CsrMatrix` instead of carrying CFDrs-local CSR traversal loops.
  `JacobiPreconditioner::new` now requires the Leto scalar contract because its
  diagonal extraction is provider-backed. Evidence in `D:/atlas/repos/CFDrs`:
  `cargo fmt -p cfd-math --check`, `cargo check -p cfd-math
  --no-default-features --lib`, `cargo nextest run -p cfd-math
  --no-default-features sparse --status-level fail` (18/18), and `cargo
  clippy -p cfd-math --no-default-features --all-targets -- -D warnings`
  pass. Residual: the sparse/linear-solver public API still exposes
  `nalgebra_sparse::CsrMatrix` and nalgebra `DVector`; the next larger slice is
  to replace that storage boundary with `leto_ops::CsrMatrix` and
  `leto::Array1`.
- [x] Consume the Leto-ops CSR×CSR provider in `cfd-math` sparse/AMG.
  `try_sparse_sparse_mul` now bridges validated `nalgebra_sparse::CsrMatrix`
  storage into `leto_ops::spgemm`, and AMG recompute/setup Galerkin products
  call the fallible Leto-backed path instead of local or nalgebra chained
  sparse multiplication. `try_sparse_transpose` now delegates AMG restriction
  construction (`R = P^T`) to `leto_ops::CsrMatrix::transpose`, removing the
  `nalgebra_sparse::transpose_as_csc` dependency from AMG setup. `try_spmv`
  now delegates the current `DVector`/CSR boundary to `leto_ops::spmv_into`,
  and the legacy entry points delegated to that provider call. The redundant
  parallel-named wrapper was subsequently deleted. Evidence in `D:/atlas/repos/CFDrs`: `cargo fmt -p
  cfd-math --check`, `cargo check -p cfd-math --no-default-features --lib`,
  `cargo nextest run -p cfd-math --no-default-features sparse --status-level
  fail` (17/17), `cargo nextest run -p cfd-math --no-default-features
  interpolation --status-level fail` (15/15), `cargo nextest run -p cfd-math
  --no-default-features amg --status-level fail` (6/6), and `cargo clippy -p cfd-math
  --no-default-features --all-targets -- -D warnings` pass. Residual:
  sparse/linear-solver storage and vector APIs still expose
  `nalgebra_sparse::CsrMatrix`/`DVector`; next slice should migrate those
  boundaries to `leto_ops::CsrMatrix`/`leto::Array1`.
- [x] Close the upstream Leto-ops CSR×CSR provider gap blocking the
  `cfd-math::sparse`/AMG migration. `leto_ops::spgemm` now provides the real
  CSR×CSR product needed for AMG Galerkin operators, and `CsrRow::nnz` covers
  sparse-pattern row cardinality checks without CFDrs-local access helpers.
  Evidence in `D:/atlas/repos/leto`: `cargo fmt -p leto-ops --check`, `cargo
  check -p leto-ops`, `cargo nextest run -p leto-ops --test ops_tests sparse
  --status-level fail` (14/14), `cargo clippy -p leto-ops --all-targets -- -D
  warnings`, and `cargo doc -p leto-ops --no-deps` pass. Residual: consume this
  provider surface in `crates/cfd-math/src/sparse` and then migrate the
  linear-solver/AMG DVector contracts.
- [x] Complete the `cfd-math::simd` Leto/Eunomia provider follow-up.
  `SimdVectorOps` is now implemented for Leto `Array1` instead of nalgebra
  `DVector`, `sparse_matvec` routes through `leto_ops::CsrMatrix`/`spmv`, and
  SIMD generic bounds now use Eunomia scalar traits. The external SIMD
  integration test builds its residual matrix with `leto_ops::CsrMatrix`
  instead of `nalgebra_sparse`. Evidence: `cargo fmt -p cfd-math --check`,
  `cargo check -p cfd-math --no-default-features --lib`, `cargo nextest run
  -p cfd-math --no-default-features simd --status-level fail` (26/26, 318
  skipped), `cargo clippy -p cfd-math --no-default-features --all-targets --
  -D warnings`, and a clean SIMD provider-residue scan. Residual: sparse and
  linear-solver surfaces still carry broader provider work.
- [x] Complete the `cfd-math::nonlinear_solver` Leto vector API follow-up.
  `JfnkConfig`, `JfnkSolver`, `JvpOperator`, GMRES work vectors, and JFNK
  tests now use Leto `Array1`/`Array2` and Eunomia `RealField`/`FloatElement`
  instead of nalgebra `DVector`/`DMatrix`/`RealField`. Anderson and JFNK now
  share one `nonlinear_solver::linalg` helper module for Leto vector math,
  bounded history operations, dense small-system helpers, norms, and scaled
  updates. Evidence: `cargo fmt -p cfd-math --check`, `cargo check -p
  cfd-math --no-default-features --lib`, `cargo nextest run -p cfd-math
  --no-default-features nonlinear --status-level fail` (9/9, 335 skipped),
  `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings`, and a clean nonlinear-solver provider-residue scan. Residual:
  sparse and linear-solver surfaces still carry broader provider work.
- [x] Complete the `cfd-math::high_order::dg` Leto dense-array
  follow-up. DG basis, solution, flux, operator, limiter, solver, and
  time-integration surfaces now use Leto `Array1`/`Array2` instead of nalgebra
  `DVector`/`DMatrix`; DG projection, derivative, RHS, and implicit Newton
  correction solves now route through Leto dense solve helpers and surface
  solver errors instead of silently falling back. DG examples and DG-related
  benchmarks compile against Leto arrays. Evidence: `cargo fmt -p cfd-math
  --check`, `cargo check -p cfd-math --no-default-features --lib`, `cargo
  check -p cfd-math --no-default-features --bench dg_benchmarks`, `cargo
  check -p cfd-math --no-default-features --bench flux_alloc_bench`, `cargo
  nextest run -p cfd-math --no-default-features dg --status-level fail`
  (62/62, 282 skipped), `cargo clippy -p cfd-math --no-default-features
  --all-targets -- -D warnings`, `cargo test --doc -p cfd-math
  --no-default-features` (3 passed, 3 ignored), and a clean
  high-order/DG-bench provider-residue scan. Residual: sparse and
  linear-solver surfaces still carry broader provider work.
- [x] Complete the `cfd-math::high_order::spectral` Leto dense-array
  follow-up. `SpectralElement`, `SpectralMesh1D`,
  `SpectralDiffOp`, `SpectralInterp`, `SpectralQuadrature`, `SpectralFilter`,
  and spectral time-integration helpers now expose and consume Leto `Array1`
  and `Array2` values instead of nalgebra `DVector`/`DMatrix`. Spectral
  assembly now uses Leto arrays for local matrices, optional local RHS values,
  and dense CSR materialization. Shared local helpers own derivative
  matrix-vector multiplication, dot products, vector construction, and
  stiffness assembly. `SpectralInterp::l2_projection` now returns
  `Result<Array1<f64>>` and maps Leto solve failures to typed solver errors
  instead of silently falling back to interpolation. Evidence: `cargo fmt -p
  cfd-math --check`, `cargo check -p cfd-math --no-default-features --lib`,
  `cargo nextest run -p cfd-math --no-default-features spectral
  --status-level fail` (13/13, 331 skipped), `cargo clippy -p cfd-math
  --no-default-features --all-targets -- -D warnings`, and a clean
  provider-residue scan over `crates/cfd-math/src/high_order/spectral`.
  Residual: sparse and linear-solver surfaces still carry broader provider
  work.
- [x] Complete the `cfd-math::high_order::weno` Eunomia scalar-provider
  follow-up. `WENO5`, `WENO7`, and the `WenoReconstruction` constructors now
  bind on `eunomia::RealField`/`FloatElement`; WENO constant conversion and
  nonlinear-weight squaring stay on the existing Eunomia helper path, and the
  final nalgebra scalar import was removed from the WENO cone. Evidence:
  `cargo fmt -p cfd-math --check`, `cargo check -p cfd-math
  --no-default-features --lib`, `cargo nextest run -p cfd-math
  --no-default-features weno --status-level fail` (6/6, 338 skipped),
  `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings`, and a clean WENO provider-residue scan. Residual: sparse and
  linear-solver surfaces still carry broader nalgebra/Leto provider work
  outside this WENO slice.
- [x] Complete the `cfd-2d`/`cfd-3d`/`cfd-validation` Eunomia-Leto
  boundary unblock slice. `cfd-2d` now exposes `Cfd2dScalar` as its
  crate-level scalar contract, routes the NS-FVM boundary/fluid path through
  that seam, and keeps direct nalgebra scalar residue isolated to
  `crates/cfd-2d/src/scalar.rs`. cfd-2d moving-wall tests and NS-FVM inlet
  setup now construct Leto `Vector3` boundary values. The LES Smagorinsky GPU
  update path now refreshes Yoshizawa SGS kinetic-energy/dissipation
  diagnostics after GPU viscosity computation, matching the CPU path.
  `cfd-3d` now has `Cfd3dScalar`, and FEM problem/stress/validation plus DES
  cfd-core fluid/boundary integrations use that seam; boundary-condition
  vector construction in bifurcation/cascade/serpentine/trifurcation/venturi
  now emits Leto vectors. `cfd-validation` now has `ValidationScalar` for
  cross-fidelity cfd-1d/2d/3d calls, and Casson/Womersley analytical wrappers
  plus 2D/3D solver-backed benchmarks use it. Verified with `rustup run
  nightly cargo fmt -p cfd-2d --check`, `cargo check -p cfd-2d
  --no-default-features --features gpu --lib`, `cargo check -p cfd-2d
  --no-default-features --lib`, `cargo clippy -p cfd-2d
  --no-default-features --features gpu --lib -- -D warnings`, `cargo nextest
  run -p cfd-2d --no-default-features --features gpu --lib --status-level
  fail` (518/518, 1 skipped), `cargo check -p cfd-3d --no-default-features
  --lib`, `cargo clippy -p cfd-3d --no-default-features --lib -- -D
  warnings`, `cargo check -p cfd-validation --no-default-features --lib`, and
  `cargo clippy -p cfd-validation --no-default-features --lib -- -D
  warnings`. Follow-up all-targets cleanup moved the remaining cfd-2d
  examples and integration tests onto `Cfd2dScalar`/Leto boundary values and
  resolved the crate's all-target clippy blockers without changing solver
  semantics. Additional evidence: `cargo clippy -p cfd-2d
  --no-default-features --features gpu --all-targets -- -D warnings` and
  `cargo nextest run -p cfd-2d --no-default-features --features gpu
  --status-level fail` (572/572, 27 skipped) pass. cfd-validation geometry
  follow-up: `Geometry2D`/`Geometry3D`, all 2D/3D geometry implementations,
  and the directly dependent bifurcation/serpentine/trifurcation/venturi
  benchmark wrappers now use Leto point/vector types and Eunomia scalar
  contracts instead of nalgebra `RealField`/point exports. Test boundary
  velocity values in the validation suite that target migrated cfd-core
  boundary APIs now construct Leto vectors. Evidence: `cargo fmt -p
  cfd-validation --check`, `cargo check -p cfd-validation
  --no-default-features --lib`, `cargo clippy -p cfd-validation
  --no-default-features --all-targets -- -D warnings`, and `cargo nextest run
  -p cfd-validation --no-default-features geometry --status-level fail`
  (11/11, 420 skipped) pass; the migrated geometry/benchmark scan finds no
  `nalgebra::RealField`, nalgebra point re-export, or `T::zero`/`T::one`
  residue. Analytical follow-up: `src/analytical/**` and `src/solutions/mod.rs`
  now return Leto `Vector3` values and bind on Eunomia scalar contracts instead
  of nalgebra analytical trait surfaces; generic identities route through
  `cfd_validation::scalar`. Evidence: `cargo check -p cfd-validation
  --no-default-features --lib`, focused `cargo nextest run -p cfd-validation
  --no-default-features analytical --status-level fail` (20/20, 411 skipped),
  `cargo clippy -p cfd-validation --no-default-features --all-targets -- -D
  warnings`, and a clean analytical/solutions residue scan.
  Error-metrics follow-up: `src/error_metrics/**` now uses Eunomia scalar
  contracts and Leto `Vector3` magnitudes instead of nalgebra `RealField`/
  `Vector3` surfaces; generic reductions route through the crate-local
  Eunomia scalar helpers. Evidence: `cargo fmt -p cfd-validation --check`,
  `cargo check -p cfd-validation --no-default-features --lib`, focused
  `cargo nextest run -p cfd-validation --no-default-features error_metrics
  --status-level fail` (21/21, 410 skipped), `cargo clippy -p
  cfd-validation --no-default-features --all-targets -- -D warnings`, and a
  clean `src/error_metrics` nalgebra/identity residue scan.
  Convergence follow-up: `src/convergence/**` now uses Eunomia scalar
  contracts and crate-local scalar helpers instead of nalgebra `RealField`,
  `T::zero()`, and `T::one()`; Richardson MMS result holders that embed the
  migrated `ConvergenceStudy<T>` carry the required Eunomia field bound.
  Evidence: `cargo check -p cfd-validation --no-default-features --lib`,
  focused `cargo nextest run -p cfd-validation --no-default-features
  convergence --status-level fail` (27/27, 404 skipped), all-target clippy,
  and a clean `src/convergence` nalgebra/identity residue scan.
  Edge-case testing follow-up: `src/edge_case_testing/**` now binds on
  Eunomia `RealField` instead of nalgebra `RealField`. Evidence: `cargo check
  -p cfd-validation --no-default-features --lib`, focused `cargo nextest run
  -p cfd-validation --no-default-features edge_case --status-level fail`
  (15/15, 416 skipped), all-target clippy, and a clean edge-case/convergence/
  error-metrics residue scan.
  Manufactured follow-up: `src/manufactured/**` now binds on Eunomia scalar
  contracts and crate-local scalar helpers instead of nalgebra `RealField`,
  `T::zero()`, and `T::one()`; the remaining cfd-2d Poisson bridge solver
  requirement is localized as `cfd_2d::Cfd2dScalar`. Evidence: `cargo check
  -p cfd-validation --no-default-features --lib`, focused `cargo nextest run
  -p cfd-validation --no-default-features manufactured --status-level fail`
  (50/50, 381 skipped), all-target clippy, and a clean manufactured residue
  scan.
  Conservation follow-up: `src/conservation/**` now uses Leto `Array1`/
  `Array2` flow fields and Eunomia scalar contracts instead of nalgebra
  `DMatrix`/`DVector`/`RealField`; conservation tolerances, reports, mass,
  momentum, energy, angular-momentum, vorticity, and GCL checks route through
  crate-local scalar helpers. Evidence: `cargo check -p cfd-validation
  --no-default-features --lib`, focused `cargo nextest run -p cfd-validation
  --no-default-features conservation --status-level fail` (18/18, 413
  skipped), all-target clippy, and a clean conservation residue scan.
  Time-integration follow-up: `cfd-math::time_stepping::stability` now accepts
  Leto `Array2`/`Array1` RK Butcher tableaus and binds on Eunomia `RealField`,
  while `cfd-validation/src/time_integration/**` uses Leto `Array1` ODE states
  for Euler/RK2/RK4 validation, stability reporting, and edge-case tests.
  Evidence: `cargo check -p cfd-math --no-default-features --lib`, `cargo
  check -p cfd-validation --no-default-features --lib`, focused `cargo nextest
  run -p cfd-math --no-default-features stability --status-level fail` (5/5,
  328 skipped), focused `cargo nextest run -p cfd-validation
  --no-default-features time_integration --status-level fail` (12/12, 419
  skipped), all-target clippy for both crates, and a clean migrated-cone
  residue scan.
  cfd-math time-stepping follow-up: `time_stepping::{traits,runge_kutta,
  adaptive,rk_chebyshev,exponential,imex}` now use Leto-backed `TimeState<T>`
  and shared `TimeMatrix<T>` where dense matrices are required; scalar
  contracts route through Eunomia and dense matrix operations route through
  Leto-ops (`matexp` and `solve`) instead of nalgebra `DVector`/`DMatrix`/LU
  bridges. `rk4_bench` and `imex_bench` target the migrated Leto-backed API.
  Evidence: `cargo fmt -p cfd-math --check`, `cargo check -p cfd-math
  --no-default-features --lib`, focused nextest filters for `runge_kutta`
  (5/5), `adaptive` (3/3), `chebyshev` (5/5), `exponential` (6/6), and
  `imex` (5/5), `cargo clippy -p cfd-math --no-default-features
  --all-targets -- -D warnings`, and residue scans showing no nalgebra,
  `DMatrix`, `DVector`, `num_traits`, `num_complex`, or ndarray hits in
  `crates/cfd-math/src/time_stepping/**` or the RK4/IMEX benches. Residual:
  broader `cfd-math` still has nalgebra-owned sparse and linear-solver
  surfaces outside the time-stepping/differentiation/SIMD cones. Differentiation
  follow-up: `cfd-math::differentiation` now returns Leto `Array1` derivative
  vectors, uses Leto `Vector3` for gradient/divergence/curl surfaces, routes
  scalar constants through Eunomia, and replaces the type-suffixed
  `first_derivative_simd_f32` API with `FiniteDifference<f32>::
  first_derivative_simd`. Evidence: `cargo fmt -p cfd-math --check`, `cargo
  check -p cfd-math --no-default-features --lib`, focused `cargo nextest run
  -p cfd-math --no-default-features differentiation --status-level fail`
  (12/12, 323 skipped), `cargo clippy -p cfd-math --no-default-features
  --all-targets -- -D warnings`, and a clean differentiation source scan for
  nalgebra, `DVector`, old scalar identities, `num_traits`, `num_complex`,
  ndarray, and the removed type-suffixed SIMD name.
- [x] Complete the `cfd-1d` Eunomia/Leto scalar-boundary unblock slice.
  Introduced `Cfd1dScalar` as the crate-local scalar contract carrying the
  remaining nalgebra linear-system requirement plus the Eunomia scalar provider
  contract required by migrated `cfd-core` fluid, boundary, and domain APIs.
  Routed `cfd-1d` domain, network, resistance, vascular, solver, transient,
  and analysis scalar bounds through that contract, switched the 1D domain
  `contains_1d` signature to Leto `Point1<T>`, and qualified ambiguous
  Eunomia scalar operations where nalgebra and Eunomia overlap. Verified with
  `rustup run nightly cargo fmt -p cfd-1d --check`, `cargo check -p cfd-1d
  --no-default-features --lib`, `cargo nextest run -p cfd-1d
  --no-default-features --status-level fail` (725/725, 3 skipped), `cargo
  clippy -p cfd-1d --no-default-features --lib -- -D warnings`, and a scan
  showing nalgebra `RealField` remains only in `crates/cfd-1d/src/scalar.rs`.
  Attempted downstream `cargo check -p cfd-2d --no-default-features
  --features gpu --lib` now reaches `cfd-2d` and is blocked by `cfd-2d`'s own
  nalgebra scalar bounds around migrated `cfd-core` boundary/fluid contracts.
  `cfd-1d --all-targets` clippy remains blocked by existing unrelated
  tests/examples/cell-separation lint debt.
- [x] Complete the `cfd-core::compute::gpu::poisson_solver` Hephaestus kernel
  orchestration slice. `GpuPoissonSolver` now wraps the existing WGPU
  device/queue in `hephaestus_wgpu::WgpuDevice`, compiles Jacobi, red-black,
  and residual WGSL entry points through `WgslMultiStorageKernel`, and routes
  Poisson field/source/residual buffers through Hephaestus `ComputeDevice`
  upload/allocation/download APIs. The solver now stores and validates the
  constructor `nx`, `ny`, `dx`, and `dy` instead of inferring a square grid from
  `phi.len()`. Verified with `rustup run nightly cargo fmt -p cfd-core
  --check`, `cargo check -p cfd-core --features gpu`, `cargo check -p cfd-core
  --no-default-features`, `cargo clippy -p cfd-core --features gpu
  --all-targets -- -D warnings`, full `cargo nextest run -p cfd-core
  --features gpu --status-level fail` (231/231), and a clean direct scan of
  `poisson_solver.rs` for raw buffer allocation, raw pipeline/layout fields,
  manual map/poll readback, and `futures`/`mpsc` mapping channels. Attempted
  downstream `cargo check -p cfd-2d --no-default-features --features gpu --lib`
  now reaches `cfd-2d` and is still blocked before the accelerated Poisson
  consumer by `cfd-2d` Eunomia/nalgebra trait-bound migration errors.
- [x] Complete the `cfd-core::compute::gpu` Hephaestus/Eunomia provider slice.
  `GpuBuffer<T>` now allocates, uploads, downloads, and writes through the
  Hephaestus `WgpuDevice` `ComputeDevice` contract, while raw `wgpu::Buffer`
  access is limited to the existing bind-group construction points. GPU buffer,
  pipeline, and kernel scalar bounds now use `eunomia::RealField` instead of
  nalgebra. Verified with `rustup run nightly cargo fmt -p cfd-core --check`,
  `cargo check -p cfd-core --features gpu`, `cargo check -p cfd-core
  --no-default-features`, `cargo clippy -p cfd-core --features gpu
  --all-targets -- -D warnings`, full `cargo nextest run -p cfd-core
  --no-default-features --status-level fail` (201/201), full `cargo nextest
  run -p cfd-core --features gpu --status-level fail` (231/231), and a clean
  scan over `crates/cfd-core/src/compute/{gpu,traits.rs}` for nalgebra scalar
  residue, old numeric identity calls, and `num_traits`. Residual: CFDrs still
  owns WGSL CFD kernels and raw bind-group/pipeline assembly; the next GPU
  increment should move reusable CFD kernel/pipeline orchestration into
  Hephaestus instead of adding downstream compatibility wrappers.
- [x] Complete the `cfd-core::geometry::{mesh,staggered}` Leto/Eunomia storage
  slice. `Mesh<T>` now stores `leto::geometry::Point3<T>`, mesh operations use
  `leto::geometry::Vector3<T>` and `leto::FixedMatrix<T, 3, 3>` for transforms,
  and mesh/staggered scalar contracts use `eunomia::RealField`. Leto now
  provides the missing `FixedMatrix<T, 3, 3> * Vector3<T>` provider operation
  with value-semantic coverage. Verified with `rustup run nightly cargo fmt -p
  cfd-core --check`, `cargo check -p cfd-core --no-default-features`, full
  `cargo nextest run -p cfd-core --no-default-features --status-level fail`
  (201/201), `cargo clippy -p cfd-core --no-default-features --all-targets --
  -D warnings`, provider `cargo fmt/check/clippy -p leto`, full provider
  `cargo nextest run -p leto --status-level fail` (171/171), and a clean scan
  over `crates/cfd-core/src/geometry/{mesh,staggered.rs}` for nalgebra and old
  scalar helper residue. Residual: this closes mesh geometry storage, not the
  full Gaia topology replacement; the next mesh increment should map or replace
  CFDrs `Mesh<T>`/`Element` topology with Gaia-owned mesh primitives.
- [x] Complete the `cfd-core` geometry/boundary/fluid-management provider
  contract slice. `geometry::shapes::Domain`, problem/aggregate management,
  solver factory/traits/configuration helpers, material/fluid-dynamics
  services, and the boundary condition/geometry/applicator family now use
  `eunomia::RealField` for scalar contracts. Geometry and boundary point/vector
  storage now use `leto::geometry::{Point1, Point2, Point3, Vector3}` where
  those contracts were in the migrated cone. Leto now provides the missing
  `Point1<T>` primitive, conditional `Eq` derives for fixed geometry values,
  and serde `std`/`alloc` feature wiring. Verified with `rustup run nightly
  cargo fmt -p cfd-core --check`, `cargo check -p cfd-core
  --no-default-features`, `cargo nextest run -p cfd-core
  --no-default-features --status-level fail` (201/201), `cargo clippy -p
  cfd-core --no-default-features --all-targets -- -D warnings`, provider
  `cargo check -p leto`, `cargo nextest run -p leto --status-level fail`
  (170/170), `cargo clippy -p leto --all-targets -- -D warnings`, and a clean
  scan over the migrated cone for nalgebra scalar imports, old `T::zero()` /
  `T::one()` calls, nalgebra geometry imports, `default_epsilon`, and nalgebra
  geometry helper calls.
- [x] Complete the `cfd-core::compute::{traits,cpu,dispatch}` Eunomia
  execution-boundary slice and the `abstractions::problem` gravity-vector
  cleanup. `ComputeKernel`, `ComputeBuffer`, `CpuBuffer`, `CpuAdvectionKernel`,
  and `ComputeDispatcher` now use `eunomia::RealField`; CPU zero values use
  Eunomia identities instead of nalgebra `T::zero()`. `ProblemParameters::
  gravity` and `ProblemBuilder::gravity` now use `leto::geometry::Vector3<T>`.
  Verified with `cargo fmt -p cfd-core --check`, `cargo check -p cfd-core
  --no-default-features`, `cargo clippy -p cfd-core --no-default-features
  --all-targets -- -D warnings`, focused `cargo nextest run -p cfd-core
  --no-default-features problem_builder cpu compute dispatcher backend
  --status-level fail` (21/21), and clean scans of the touched files for
  nalgebra vector ownership, `DVector`, `num_traits`, ndarray, rayon, tokio,
  rustfft, and `T::zero()`. Residual: `Problem<T>` still imports
  `nalgebra::RealField` because `Domain<T>` and `FluidTrait<T>` still own that
  scalar contract.
- [x] Complete the `cfd-core::abstractions::state` Leto/Eunomia provider slice.
  `FieldData::Scalar` now stores `leto::Array1<T>` instead of nalgebra
  `DVector<T>`, vector fields use `leto::geometry::Vector3<T>`, and
  `SimulationState`/`FieldState` bind on `eunomia::RealField`. The slice
  required a Leto provider gap fix: owned Leto arrays now serialize and
  deserialize with validated layout/storage reconstruction. Verified with
  `cargo fmt -p leto --check`, focused `cargo nextest run -p leto
  owned_array_round_trips_shape_and_values_through_serde --status-level fail`
  (1/1), `cargo clippy -p leto --all-targets -- -D warnings`, `cargo fmt -p
  cfd-core --check`, `cargo check -p cfd-core --no-default-features`, `cargo
  clippy -p cfd-core --no-default-features --all-targets -- -D warnings`,
  focused `cargo nextest run -p cfd-core --no-default-features field_state
  field_data test_core_state --status-level fail` (3/3), and a clean scan of
  `crates/cfd-core/src/abstractions/state.rs` for nalgebra, `DVector`,
  `num_traits`, ndarray, rayon, tokio, wgpu, and cuda.
- [x] Complete the `cfd-core::compute::time` Leto/Eunomia provider slice.
  `TimeIntegrator::State` for `ForwardEuler`, `RungeKutta2`, `RungeKutta4`,
  `BackwardEuler`, and `CrankNicolson` is now `leto::Array1<T>` instead of
  nalgebra `DVector<T>`. The time-step controllers now use Eunomia
  `RealField`/`FloatElement`/`NumericElement` constants and power functions
  instead of nalgebra/num-traits conversion and epsilon APIs. Verified with
  `cargo fmt -p cfd-core --check`, `cargo check -p cfd-core
  --no-default-features`, `cargo check -p cfd-core --no-default-features
  --tests`, `cargo clippy -p cfd-core --no-default-features --lib -- -D
  warnings`, `cargo clippy -p cfd-core --all-targets -- -D warnings`, focused
  `cargo nextest run -p cfd-core --no-default-features time forward_euler
  runge_kutta backward_euler crank adaptive_controller variable_controller
  --status-level fail` (17/17), and a clean provider scan over
  `crates/cfd-core/src/compute/time` for `nalgebra`, `DVector`, `ndarray`,
  `num_traits`, `rayon`, `tokio`, `wgpu`, and `cuda`. The follow-up feature
  hygiene is also closed: `turbulence_benchmark` is now marked with
  `required-features = ["gpu"]`, so `cargo clippy -p cfd-core
  --no-default-features --all-targets -- -D warnings` passes.
- [x] Complete the `cfd-core::compute::solver::config` Eunomia scalar slice.
  `SolverConfig` and its convergence, numerical, linear, network, and builder
  configuration types now bind on `eunomia::RealField`; default relaxation uses
  the Eunomia numeric constant instead of nalgebra's `T::one()`. Downstream
  config holders in cfd-2d FDM/SIMPLE/pressure-velocity/vorticity, cfd-3d FEM,
  and cfd-validation MMS Richardson now carry the Eunomia config bound, and
  `cfd-3d::spectral::solver` no longer imports `NalgebraRealField`. Verified
  with `cargo fmt -p cfd-core -p cfd-2d -p cfd-3d -p cfd-validation --check`,
  `cargo check -p cfd-core --no-default-features`, `cargo check -p cfd-2d
  --no-default-features --lib`, `cargo check -p cfd-3d --no-default-features
  --lib`, `cargo check -p cfd-3d --no-default-features --test property_tests
  --test domain_solver_validation --test poisson_validation`, `cargo check -p
  cfd-validation --no-default-features --lib`, scoped clippy for cfd-core,
  cfd-2d, cfd-3d, and cfd-validation, focused `cargo nextest run -p cfd-core
  --no-default-features solver_config --status-level fail` (2/2), focused
  `cargo nextest run -p cfd-2d --no-default-features fdm pressure_velocity
  simple vorticity --status-level fail` (41/41), focused `cargo nextest run
  -p cfd-3d --no-default-features spectral poisson fem --status-level fail`
  (89/89), and scans showing `spectral::solver` has no nalgebra/ndarray/num/
  rustfft/runtime-provider residue. Residual solver work remains in
  `cfd_core::compute::solver::{traits,convergence,direct,iterative,monitor}`
  because those traits still depend on the nalgebra-bound
  `abstractions::Problem<T>` contract.
- [x] Complete the `cfd-3d::spectral::poisson` Leto/Leto-ops provider slice.
  `PoissonSolver` now builds Chebyshev Laplacian terms as `leto::Array2`,
  uses `leto_ops::MatrixProduct::kron` for tensor-product assembly, applies
  boundary rows in Leto storage, solves through `leto_ops::MatrixSolve`, and
  returns `leto::Array1`. `PoissonProblem` no longer requires a nalgebra
  scalar bound. The only remaining direct nalgebra hit in the spectral Poisson
  cone is the `NalgebraRealField` alias in `spectral::solver`, inherited from
  `cfd_core::compute::solver::SolverConfig`. Verified with `cargo fmt -p
  cfd-3d --check`, `cargo check -p cfd-3d --no-default-features --lib`,
  `cargo check -p cfd-3d --no-default-features --test robustness_tests --test
  poisson_validation`, `cargo clippy -p cfd-3d --no-default-features --lib
  --test robustness_tests --test poisson_validation -- -D warnings`, and
  focused `cargo nextest run -p cfd-3d --no-default-features poisson
  --status-level fail` (12/12). Next spectral increment: migrate
  `cfd-core::compute::solver::SolverConfig` off nalgebra scalar bounds so
  `spectral::solver` can drop `NalgebraRealField`.
- [x] Complete the `cfd-3d::spectral::chebyshev` Leto/Eunomia provider slice.
  `ChebyshevPolynomial` now stores its differentiation operators as
  `leto::Array2`, accepts/returns `leto::Array1` for first and second
  derivative application, validates vector lengths with
  `cfd_core::Error::DimensionMismatch`, and uses `eunomia::RealField` for its
  scalar contract. `spectral::basis` now uses Eunomia `RealField`; the spectral
  solver documents the temporary dual Eunomia/nalgebra scalar bound required by
  the current Poisson LU boundary. Co-located Chebyshev tests, `lib.rs`
  Chebyshev tests, and robustness Chebyshev tests no longer construct
  nalgebra `DVector`s. Verified with `cargo fmt -p cfd-3d --check`, `cargo
  check -p cfd-3d --no-default-features --lib`, `cargo check -p cfd-3d
  --no-default-features --test robustness_tests --test poisson_validation`,
  `cargo clippy -p cfd-3d --no-default-features --lib --test
  robustness_tests --test poisson_validation -- -D warnings`, focused `cargo
  nextest run -p cfd-3d --no-default-features chebyshev --status-level fail`
  (20/20), and focused `cargo nextest run -p cfd-3d --no-default-features
  poisson --status-level fail` (12/12). The follow-on Poisson slice now closes
  the former nalgebra `DMatrix`/`DVector` LU/Kronecker boundary; the remaining
  `NalgebraRealField` alias in `spectral::solver` is inherited from upstream
  `cfd_core::compute::solver::SolverConfig`.
- [x] Complete the `cfd-3d::level_set` Atlas provider slice. `level_set::{
  solver,advection,weno}` now uses `leto::geometry::Vector3` for velocity
  storage/transport and a level-set-local `LevelSetScalar` seam backed by
  `eunomia::RealField`/`FloatElement`; the dead cfd-3d crate-level
  `CfdScalar` nalgebra-bound trait was removed. `tests/level_set_tests.rs` and
  the level-set call sites in `tests/robustness_tests.rs` now construct Leto
  vectors. Verified with `cargo fmt -p cfd-3d --check`, `cargo check -p
  cfd-3d --no-default-features --lib`, `cargo check -p cfd-3d
  --no-default-features --test level_set_tests`, `cargo check -p cfd-3d
  --no-default-features --test robustness_tests`, `cargo clippy -p cfd-3d
  --no-default-features --lib --test level_set_tests --test robustness_tests
  -- -D warnings`, focused `cargo nextest run -p cfd-3d --no-default-features
  level_set --status-level fail` (13/13), and a clean forbidden-provider scan
  over `crates/cfd-3d/src/level_set` plus `tests/level_set_tests.rs` for
  `nalgebra`, `DMatrix`, `DVector`, `CfdScalar`, `crate::scalar`,
  `num_traits`, `ndarray`, `rayon`, `tokio`, `rustfft`, `wgpu`, and `cuda`.
  Residual cfd-3d provider work remains in FEM, spectral Chebyshev/Poisson,
  IBM, turbulence, validation, and unrelated robustness-test nalgebra surfaces.
- [x] Complete the `cfd-3d::vof` Atlas provider slice. `vof::cavitation_solver`
  now exposes `CavitationField = leto::Array2<f64>` for pressure, density,
  volume-fraction, inception, damage, bubble-radius, nuclei, and
  sonoluminescence fields, with explicit row-major Leto field offsets and
  VOF-alpha offset mapping. `tests/cavitation_solver_validation.rs` now builds
  cavitation dense fields with Leto `Array2`, and `vof::scalar::VofScalar` now
  uses `eunomia::RealField` instead of `nalgebra::RealField`. Verified with
  `cargo fmt -p cfd-3d --check`, `cargo check -p cfd-3d
  --no-default-features --lib`, `cargo check -p cfd-3d --no-default-features
  --test cavitation_solver_validation`, `cargo clippy -p cfd-3d
  --no-default-features --lib --test cavitation_solver_validation -- -D
  warnings`, focused `cargo nextest run -p cfd-3d --no-default-features
  cavitation --status-level fail` (23/23), focused `cargo nextest run -p
  cfd-3d --no-default-features vof --status-level fail` (42/42), and a clean
  forbidden-provider scan over `crates/cfd-3d/src/vof` plus the VOF/cavitation
  tests for `nalgebra`, `DMatrix`, `DVector`, `num_traits`, `ndarray`, `rayon`,
  `tokio`, `rustfft`, `wgpu`, and `cuda`. Residual cfd-3d provider work is now
  outside VOF: FEM, spectral Chebyshev/Poisson, level-set, IBM, turbulence, and
  validation surfaces still require Leto/Gaia/Eunomia/Hephaestus/Moirai slices
  before the cfd-3d manifest can drop broader legacy providers.
- [x] Migrate the `cfd-3d::vof` cavitation velocity-provider surface from
  nalgebra `Vector3` to Leto. `vof::{cavitation_solver,bubble_dynamics}` and
  `tests/cavitation_solver_validation.rs` now use `leto::geometry::Vector3`
  for cavitation velocity input and bubble update calls; the cavitation solver
  copies directly into the Leto-owned `VofSolver` velocity buffer without a
  nalgebra-to-Leto boundary conversion. Verified with `cargo fmt -p cfd-3d
  --check`, `cargo check -p cfd-3d --no-default-features --lib`, `cargo check
  -p cfd-3d --no-default-features --test cavitation_solver_validation`, `cargo
  clippy -p cfd-3d --no-default-features --lib --test
  cavitation_solver_validation -- -D warnings`, focused `cargo nextest run -p
  cfd-3d --no-default-features cavitation --status-level fail` (23/23), and a
  static scan showing no direct nalgebra `Vector3` residue in the cavitation
  VOF cone. The remaining cavitation dense-field work from this item is closed
  by the later Leto `CavitationField` slice above.
- [x] Migrate the non-cavitation `cfd-3d::vof` vector-provider surface from
  nalgebra `Vector3` to Leto. `vof::{solver,reconstruction,initialization,
  plic_geometry,advection}` and `tests/vof_tests.rs` now use
  `leto::geometry::Vector3`. Verified with
  `cargo fmt -p cfd-3d --check`, `cargo check -p cfd-3d --no-default-features
  --test vof_tests`, `cargo check -p cfd-3d --no-default-features --test
  robustness_tests`, `cargo clippy -p cfd-3d --no-default-features --lib
  --test vof_tests -- -D warnings`, focused `cargo nextest run -p cfd-3d
  --no-default-features vof --status-level fail` (42/42), `git diff --check`
  over the touched files, and a static scan showing no direct nalgebra
  `Vector3` residue in the migrated VOF files. The cavitation dense-storage
  residue from this item is closed by the later Leto `CavitationField` slice
  above.
- [x] Close direct nalgebra `Vector2` ownership in `cfd-validation` source and
  tests. `manufactured::navier_stokes`, `conservation::{momentum,
  angular_momentum,mod}`, the vorticity-stream benchmark, and MMS/conservation
  test consumers now use `leto::geometry::Vector2` with indexed component
  access. Verified with `cargo fmt -p cfd-validation --check`, `cargo check -p
  cfd-validation --no-default-features --tests`, `cargo clippy -p
  cfd-validation --no-default-features --lib -- -D warnings`, focused `cargo
  nextest run -p cfd-validation --no-default-features manufactured mms
  conservation taylor momentum angular --status-level fail` (104/104), and a
  direct scan over `crates/cfd-validation/{src,tests}` showing no direct
  nalgebra `Vector2` imports/usages. Residual validation nalgebra remains in
  dense `DMatrix`, 3D `Vector3`, and geometry point/vector surfaces that need
  separate Leto/Gaia provider slices.
- [x] Close direct nalgebra `Vector2` ownership in `cfd-2d` source, tests,
  examples, and benches. `physics::vorticity_stream`, `piso_algorithm::
  corrector`, `solvers::lbm::solver`, `physics::immersed_boundary`,
  `examples/blood_venturi.rs`, and `benches/solver_benchmarks.rs` now use
  `leto::geometry::Vector2` with indexed component access. The downstream
  `cfd-validation::benchmarks::vorticity_stream` consumer now uses the same Leto
  vector contract. Verified with `cargo fmt -p cfd-2d -p cfd-validation
  --check`, `cargo check -p cfd-2d --no-default-features --examples --benches`,
  `cargo check -p cfd-validation --no-default-features`, focused `cargo clippy
  -p cfd-2d --no-default-features --example blood_venturi --bench
  solver_benchmarks -- -D warnings`, `cargo clippy -p cfd-validation
  --no-default-features --lib -- -D warnings`, focused `cargo nextest run -p
  cfd-2d --no-default-features vorticity corrector lbm immersed --status-level
  fail` (44/44), and a direct scan over cfd-2d source/tests/examples/benches
  plus the touched cfd-validation benchmark showing no direct nalgebra
  `Vector2` imports/usages. Residual cfd-2d nalgebra remains in scalar
  `RealField`, boundary `Vector3`, dense `DVector`/`DMatrix`, and
  nalgebra-sparse storage surfaces.
- [x] Migrate the `cfd-2d` pressure-velocity/SIMPLEC/PIMPLE vector workspace
  family from nalgebra `Vector2` to Leto. `fields.rs`,
  `pressure_velocity::{solver,pressure,rhie_chow}`, `simplec_pimple::{
  algorithms,diagnostics,interpolation,pimple,simplec,solver}`,
  `physics::momentum::interpolation` reference tests, and
  `tests/simplec_pimple_validation.rs` now use `leto::geometry::Vector2` for
  2D velocity workspaces, caches, correction buffers, face velocities, and
  validation setup. Component access now uses Leto indexing instead of nalgebra
  `.x`/`.y` fields. Verified with `cargo fmt -p cfd-2d --check`, `cargo check
  -p cfd-2d --no-default-features`, `cargo clippy -p cfd-2d
  --no-default-features --lib -- -D warnings`, focused `cargo nextest run -p
  cfd-2d --no-default-features pressure_velocity simplec pimple rhie
  --status-level fail` (29/29), and a direct scan over the touched family
  showing no `nalgebra::Vector2` or `use nalgebra::{RealField, Vector2}`
  residue. Remaining nalgebra in this family is scalar `RealField`, `Vector3`
  boundary-condition input, and `DVector` linear-solver storage.
- [x] Move the `cfd-2d` Ghia pressure-correction helper test off direct
  nalgebra vector construction. `tests/ghia_cavity_simplec_validation.rs` now
  uses `leto::geometry::Vector2` for the local divergence-free velocity helper
  instead of `nalgebra::Vector2`. Verified with `cargo fmt -p cfd-2d --check`,
  `cargo check -p cfd-2d --no-default-features --test
  ghia_cavity_simplec_validation`, focused `cargo nextest run -p cfd-2d
  --no-default-features test_pressure_correction_basic --status-level fail`
  (1/1), and a file scan showing no direct `nalgebra`/`nalgebra-sparse`
  residue in `tests/ghia_cavity_simplec_validation.rs`. Remaining cfd-2d test
  nalgebra ownership is concentrated in `tests/simplec_pimple_validation.rs`
  and the production SIMPLEC/PIMPLE/pressure-velocity vector workspaces.
- [x] Migrate the `cfd-2d::physics::turbulence::validation` scalar contract
  from direct nalgebra bounds to Eunomia. `validation::{mod,rans,les_des,
  benchmarks}` no longer imports or bounds `nalgebra::RealField`; the
  validation scalar contract now uses `eunomia::RealField`/`FloatElement`, and
  the validation arrays/vectors stay on Leto `Array2`/`Vector2`. Verified with
  `cargo fmt -p cfd-2d --check`, `cargo check -p cfd-2d
  --no-default-features`, `cargo clippy -p cfd-2d --no-default-features --lib
  -- -D warnings`, focused `cargo nextest run -p cfd-2d
  --no-default-features turbulence validation --status-level fail` (207/207),
  and a direct scan over `crates/cfd-2d/src/physics/turbulence/validation`
  showing no `nalgebra`, `nalgebra-sparse`, `nalgebra_sparse`, `DVector`, or
  `DMatrix` residue. This does not remove cfd-2d's direct
  `nalgebra`/`nalgebra-sparse` manifest entries because other cfd-2d modules
  and tests still own those provider surfaces.
- [x] Remove the direct `num-traits` dependency from `cfd-2d`. The final
  residue pass removed direct `num_traits` imports/bounds from
  `physics::turbulence::validation::{mod,les_des,benchmarks}`,
  `physics::immersed_boundary` test code, and
  `tests/ghia_cavity_simplec_validation.rs`; `crates/cfd-2d/Cargo.toml` no
  longer declares `num-traits`. Turbulence validation construction now uses
  `eunomia::{FloatElement,RealField}` bounds, existing Leto `Array2` validation
  arrays, and concrete `f64` validation/test arithmetic where the API is
  already concrete. Verified with `cargo fmt -p cfd-2d --check`, `cargo check
  -p cfd-2d --no-default-features --tests`, `cargo clippy -p cfd-2d
  --no-default-features --lib -- -D warnings`, focused `cargo nextest run -p
  cfd-2d --no-default-features turbulence validation immersed ghia
  --status-level fail` (213/213), `cargo tree -p cfd-2d --depth 1` showing no
  direct `num-traits`, `git diff --check`, and a direct scan over
  `crates/cfd-2d/{src,tests,Cargo.toml}` showing no `num_traits`,
  `FromPrimitive`, `ToPrimitive`, or `num-traits` residue. `num-traits` remains
  in the workspace through other crates/transitive dependencies; cfd-2d still
  has direct `nalgebra`/`nalgebra-sparse` provider work. This entry supersedes
  older cfd-2d slice-local residual notes that said the cfd-2d manifest still
  declared `num-traits`.
- [x] Migrate the `cfd-2d::physics::momentum` setup/boundary scalar-provider
  seam to Eunomia. `physics::momentum::setup` no longer imports or bounds
  direct `num_traits::{FromPrimitive,ToPrimitive}`. `physics::momentum::boundary`
  and `boundary::directional` now route no-slip/slip/symmetry/periodic/outflow
  zero/one row entries, zero-gradient comparisons, quadratic wall-extrapolation
  constants, and corner consistency absolute-value checks through
  `crates/cfd-2d/src/scalar.rs` and `eunomia::FloatElement` instead of direct
  `num_traits::FromPrimitive`, `T::from_f64(...).unwrap_or_else`, `T::zero()`,
  `T::one()`, or scalar `.abs()` usage. Verified with `cargo fmt -p cfd-2d
  --check`, `cargo check -p cfd-2d --no-default-features`, `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings`, focused `cargo nextest
  run -p cfd-2d --no-default-features boundary --status-level fail` (42/42),
  `git diff --check`, and a direct-provider scan over the touched files showing
  no direct scalar-provider residue. The cfd-2d manifest still declares
  `num-traits` before the final crate-level cfd-2d dependency removal pass.
- [x] Migrate the `cfd-2d::physics::momentum::interpolation` scalar-provider
  seam to Eunomia. The coefficient-aware Rhie-Chow interpolation path now
  routes the tiny diagonal floor, zero initialization for face arrays, and the
  harmonic face coefficient factor through `crates/cfd-2d/src/scalar.rs` and
  `eunomia::FloatElement` instead of direct `num_traits::FromPrimitive`,
  `T::from_f64(...).unwrap_or_else`, `T::zero()`, or `T::one()` usage.
  Verified with `cargo fmt -p cfd-2d --check`, `cargo check -p cfd-2d
  --no-default-features`, `cargo clippy -p cfd-2d --no-default-features --lib
  -- -D warnings`, focused `cargo nextest run -p cfd-2d --no-default-features
  coefficient_aware_interpolation_matches_exact_reference --status-level fail`
  (1/1), `git diff --check`, and a direct-provider scan over
  `crates/cfd-2d/src/physics/momentum/interpolation.rs` showing no direct
  `num_traits` or scalar-constructor residue. The cfd-2d manifest still
  declares `num-traits` because direct residues remain in Ghia validation tests,
  immersed-boundary test code, turbulence validation, and f64-only/test scalar
  surfaces.
- [x] Migrate the `cfd-2d` problem/streamtube scalar-provider seam to Atlas
  providers. `problem.rs` now stores incompressible problem and solution
  velocity fields with `leto::geometry::Vector2` and routes initial pressure,
  velocity-magnitude maxima, and pressure maxima through
  `crates/cfd-2d/src/scalar.rs`/Eunomia instead of local nalgebra vector
  storage or direct `T::zero()` folds. `physics::streamtube::partitioning`
  now uses `eunomia::{FloatElement,NumericElement}` and the crate-local scalar
  adapter instead of `num_traits::{Float,FromPrimitive}` for constants,
  absolute values, and square roots. Added value-semantic problem tests for
  provider-owned initialization and maxima. Verified with `cargo fmt -p
  cfd-2d --check`, `cargo check -p cfd-2d --no-default-features`, `cargo
  clippy -p cfd-2d --no-default-features --lib -- -D warnings`, focused
  `cargo nextest run -p cfd-2d --no-default-features problem streamtube
  separating --status-level fail` (4/4), `git diff --check`, and a direct
  provider scan over both touched files showing no direct residue. The
  remaining `problem.rs` `nalgebra::RealField` bound is inherited from
  `cfd-core` boundary/fluid types and remains an upstream provider blocker.
  The cfd-2d manifest still declares `num-traits` because direct residues
  remain in immersed-boundary tests, momentum setup/interpolation/boundary,
  turbulence validation, and f64-only/test scalar surfaces.
- [x] Migrate the `cfd-2d::solvers::fvm` scalar and vector-provider surface to
  Atlas providers. `solvers::fvm::{config,flux,geometry,solver}` now route FVM
  defaults, face generation constants, residual magnitudes, max terms, flux
  Peclet calculations, nonfinite diffusion validation, and module test
  absolute-value assertions through `crates/cfd-2d/src/scalar.rs` and
  `eunomia::{FloatElement,NumericElement}` instead of direct `num_traits`,
  `T::from_*`, `.to_f64()`, `T::zero()`, `T::one()`, or scalar `.abs()`
  usage. FVM face geometry and solver velocity vectors now use
  `leto::geometry::Vector2` instead of `nalgebra::Vector2`. Power-law and
  hybrid flux calculators keep diffusion in `T` instead of round-tripping
  through `f64`. Verified with `cargo fmt -p cfd-2d --check`, `cargo check -p
  cfd-2d --no-default-features`, `cargo clippy -p cfd-2d
  --no-default-features --lib -- -D warnings`, focused `cargo nextest run -p
  cfd-2d --no-default-features fvm --status-level fail` (25/25), `git diff
  --check`, and a direct-provider scan over `crates/cfd-2d/src/solvers/fvm`
  showing no direct residue. The cfd-2d manifest still declares `num-traits`
  because direct residues remain outside this slice in problem setup,
  immersed-boundary tests, streamtube partitioning, momentum, energy,
  turbulence, and f64-only/test scalar surfaces.
- [x] Migrate the `cfd-2d::physics::vorticity_stream` scalar-provider surface
  to Eunomia. `VorticityStreamConfig` and `VorticityStreamSolver` now route
  default tolerances, timestep/SOR constants, zero/one identities, Laplacian
  gradient factors, SOR residuals, upwind sign checks, velocity-recovery
  denominators, boundary-vorticity factors, convergence checks, and the
  continuity test absolute-value assertion through `crates/cfd-2d/src/scalar.rs`
  and `eunomia::{FloatElement,NumericElement}` instead of direct
  `num_traits::FromPrimitive`, `T::from_*`, `T::zero()`, `T::one()`, or scalar
  `.abs()` usage. Verified with `cargo fmt -p cfd-2d --check`, `cargo check -p
  cfd-2d --no-default-features`, `cargo clippy -p cfd-2d
  --no-default-features --lib -- -D warnings`, focused `cargo nextest run -p
  cfd-2d --no-default-features vorticity --status-level fail` (4/4), `git
  diff --check`, and a direct-provider scan over
  `crates/cfd-2d/src/physics/vorticity_stream.rs` showing no direct residue.
  The cfd-2d manifest still declares `num-traits` because direct residues
  remain outside this slice in problem setup, FVM, immersed-boundary tests,
  streamtube partitioning, momentum, energy, turbulence, and f64-only/test
  scalar surfaces.
- [x] Migrate the `cfd-2d` acoustic-drift scalar-provider surface to Eunomia.
  `physics::acoustics::gorkov` and `solvers::drift_diffusion_2d` now route
  ARF material constants, compressibility identities, contrast-factor
  constants, standing-wave sine evaluation, drift-diffusion zero/one
  identities, under-relaxation constants, Patankar upwind max terms,
  pivot-threshold checks, and convergence residual magnitudes through
  `crates/cfd-2d/src/scalar.rs` and
  `eunomia::{FloatElement,NumericElement}` instead of direct
  `num_traits::{Float,FromPrimitive}`, `T::from_*`, `T::zero()`, `T::one()`,
  `Float::`, or scalar `.abs()` usage. Verified with `cargo fmt -p cfd-2d
  --check`, `cargo check -p cfd-2d --no-default-features`, `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings`, focused `cargo nextest
  run -p cfd-2d --no-default-features gorkov --status-level fail` (4/4),
  focused `cargo nextest run -p cfd-2d --no-default-features
  drift_diffusion --status-level fail` (1/1), `git diff --check`, and a
  direct-provider scan over both files showing no direct residue. The cfd-2d
  manifest still declares `num-traits` because direct residues remain outside
  this slice in problem setup, FVM, vorticity-stream, immersed-boundary tests,
  streamtube partitioning, momentum, energy, turbulence, and f64-only/test
  scalar surfaces.
- [x] Migrate the `cfd-2d::schemes::tvd` scalar-provider family to Eunomia.
  `schemes::tvd::{mod,muscl,quick}` now route TVD limiter identities,
  min/max/absolute-value operations, MUSCL slope ratios, QUICK interpolation
  constants, and MUSCL2/MUSCL3 boundary constants through
  `crates/cfd-2d/src/scalar.rs` and
  `eunomia::{FloatElement,NumericElement}` instead of direct
  `num_traits::{FromPrimitive,ToPrimitive}`, `T::from_*`, `T::zero()`,
  `T::one()`, `.to_f64()`, or scalar `.abs()` usage. The limiter path no
  longer round-trips `r` through `f64`; `FluxLimiter::apply` remains generic
  in `T`. Verified with `cargo fmt -p cfd-2d --check`, `cargo check -p
  cfd-2d --no-default-features`, `cargo clippy -p cfd-2d
  --no-default-features --lib -- -D warnings`, focused `cargo nextest run -p
  cfd-2d --no-default-features tvd --status-level fail` (32/32), focused
  `cargo nextest run -p cfd-2d --no-default-features muscl --status-level
  fail` (14/14), `git diff --check`, and a direct-provider scan over
  `crates/cfd-2d/src/schemes/tvd` showing no direct residue. The cfd-2d
  manifest still declares `num-traits` because direct residues remain outside
  this slice in problem setup, drift-diffusion, FVM, vorticity-stream,
  acoustics, immersed-boundary tests, streamtube partitioning, momentum,
  energy, turbulence, and f64-only/test scalar surfaces.
- [x] Migrate the `cfd-2d::solvers::lbm` scalar-provider surface to Eunomia.
  LBM solver configuration, boundary reconstruction, streaming push-zeroing,
  macroscopic/lattice helpers, Carreau-Yasuda assertions, and Shan-Chen
  pseudopotential force calculations now route scalar construction, zero/one
  identities, integer lattice-direction casts, `sqrt`, `exp`, max, and
  absolute-value checks through `crates/cfd-2d/src/scalar.rs` and
  `eunomia::{FloatElement,NumericElement,CastFrom}` instead of direct
  `num_traits::{Float,FromPrimitive}`, `T::from_*`, `T::zero()`, `T::one()`,
  `Float::`, or scalar `.abs()` usage. Verified with `cargo fmt -p cfd-2d
  --check`, `cargo check -p cfd-2d --no-default-features`, `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings`, focused `cargo nextest
  run -p cfd-2d --no-default-features lbm --status-level fail` (31/31),
  `git diff --check`, and a direct-provider scan over
  `crates/cfd-2d/src/solvers/lbm` showing no direct residue. The cfd-2d
  manifest still declares `num-traits` because direct residues remain outside
  this slice in problem setup, drift-diffusion, FVM, vorticity-stream,
  immersed-boundary tests, streamtube partitioning, TVD schemes, momentum,
  energy, turbulence, and f64-only/test scalar surfaces.
- [x] Migrate the `cfd-2d::solvers::simple` scalar-provider surface to Eunomia.
  `solvers::simple::{algorithm,momentum,pressure}` now route Patankar default
  constants, workspace zero initialization, stagnant-cell safeguards,
  pressure-correction identities, neighbor-count conversion, finite checks,
  and module test absolute-value assertions through
  `crates/cfd-2d/src/scalar.rs` and `eunomia::{FloatElement,NumericElement}`
  instead of direct `num_traits::{FromPrimitive,ToPrimitive}`, `T::from_*`,
  `T::zero()`, `T::one()`, scalar `.abs()`, or direct finite checks. Verified
  with `cargo fmt -p cfd-2d --check`, `cargo check -p cfd-2d
  --no-default-features`, `cargo clippy -p cfd-2d --no-default-features --lib
  -- -D warnings`, focused `cargo nextest run -p cfd-2d
  --no-default-features simple --status-level fail` (19/19), `git diff
  --check`, and a direct-provider scan over `crates/cfd-2d/src/solvers/simple`
  showing no direct residue. The cfd-2d manifest still declares `num-traits`
  because direct residues remain outside this slice in vorticity-stream,
  drift-diffusion, acoustics, momentum, turbulence, energy,
  streamtube/problem/schemes, FVM, and f64-only/test scalar surfaces.
- [x] Migrate the `cfd-2d::solvers::venturi_flow` scalar-provider surface to
  Eunomia. `venturi_flow::{mod,solver}` now routes ISO/default constants,
  beta clamping, domain bounds, analytical Bernoulli/viscous constants,
  energy-dissipation guards, stretched-grid index conversion, generic
  sine/min/max dispatch, diagnostic f64 conversion, throat-column selection,
  inlet/outlet/throat reductions, and validation error metrics through
  `crates/cfd-2d/src/scalar.rs` and `eunomia::{FloatElement,NumericElement}`
  instead of direct `num_traits::{Float,FromPrimitive,ToPrimitive}`,
  `T::from_*`, `T::zero()`, `T::one()`, `Float::`, `.to_f64()`, or scalar
  `.abs()` usage. Verified with `cargo fmt -p cfd-2d --check`, `cargo check
  -p cfd-2d --no-default-features`, `cargo clippy -p cfd-2d
  --no-default-features --lib -- -D warnings`, focused `cargo nextest run -p
  cfd-2d --no-default-features venturi --status-level fail` (13/13), `git
  diff --check`, and a direct-provider scan over `crates/cfd-2d/src/solvers/
  venturi_flow` showing no direct residue. The cfd-2d manifest still declares
  `num-traits` because direct residues remain outside this slice in
  vorticity-stream, LBM, drift-diffusion, acoustics, momentum,
  turbulence, energy, streamtube/problem/schemes, FVM, and f64-only/test scalar
  surfaces.
- [x] Migrate the `cfd-2d` branching-flow scalar-provider family to Eunomia.
  `solvers::{bifurcation_flow,n_furcation_flow}` now route geometry constants,
  branch-index conversion, zero identities, generic sin/cos dispatch, flux
  accumulation, and mass-balance absolute-value normalization through
  `crates/cfd-2d/src/scalar.rs` and `eunomia::{FloatElement,NumericElement}`
  instead of direct `num_traits::{Float,FromPrimitive,ToPrimitive}`,
  `T::from_*`, `T::zero()`, or `Float::` usage. Verified with `cargo fmt -p
  cfd-2d --check`, `cargo check -p cfd-2d --no-default-features`, `cargo
  clippy -p cfd-2d --no-default-features --lib -- -D warnings`, focused
  `cargo nextest run -p cfd-2d --no-default-features bifurcation
  --status-level fail` (9/9), focused `cargo nextest run -p cfd-2d
  --no-default-features n_furcation --status-level fail` (1/1), `git diff
  --check`, and a direct-provider scan over both branching files showing no
  direct residue. The cfd-2d manifest still declares `num-traits` because
  direct residues remain outside this slice in vorticity-stream, LBM,
  drift-diffusion, acoustics, momentum, turbulence, energy,
  streamtube/problem/schemes, FVM, and f64-only/test scalar surfaces.
- [x] Migrate the `cfd-2d::solvers::cross_junction_flow` scalar-provider
  surface to Eunomia. `CrossJunctionGeometry`, `CrossJunctionSolver2D`, flux
  accumulation, mass-balance normalization, and module test absolute-value
  checks now route through `crates/cfd-2d/src/scalar.rs` and
  `eunomia::{FloatElement,NumericElement}` instead of direct
  `num_traits::{Float,FromPrimitive,ToPrimitive}`, `T::from_*`, `T::zero()`,
  or scalar `.abs()` usage. Verified with `cargo fmt -p cfd-2d --check`,
  `cargo check -p cfd-2d --no-default-features`, `cargo clippy -p cfd-2d
  --no-default-features --lib -- -D warnings`, focused `cargo nextest run -p
  cfd-2d --no-default-features cross_junction --status-level fail` (5/5),
  `git diff --check`, and a direct-provider scan over
  `crates/cfd-2d/src/solvers/cross_junction_flow.rs` showing no direct
  residue. The cfd-2d manifest still declares `num-traits` because direct
  residues remain outside this slice in LBM, drift-diffusion,
  and f64-only/test scalar surfaces.
- [x] Migrate the `cfd-2d::solvers::poiseuille` scalar-provider surface to
  Eunomia. `poiseuille::{mod,numerics}` now routes configuration defaults,
  grid-index conversion, zero/one identities, tridiagonal workspaces,
  harmonic-mean constants, shear-rate magnitudes, viscosity residual norms,
  Thomas pivot thresholds, error diagnostics, and module test absolute-value
  checks through `crates/cfd-2d/src/scalar.rs` and
  `eunomia::{FloatElement,NumericElement}` instead of direct
  `num_traits::{Float,FromPrimitive}`, `T::from_*`, `T::zero()`, `T::one()`,
  or scalar `.abs()` usage. Verified with `cargo fmt -p cfd-2d --check`,
  `cargo check -p cfd-2d --no-default-features`, `cargo clippy -p cfd-2d
  --no-default-features --lib -- -D warnings`, focused `cargo nextest run -p
  cfd-2d --no-default-features poiseuille --status-level fail` (7/7), `git
  diff --check`, and a direct-provider scan over
  `crates/cfd-2d/src/solvers/poiseuille` showing no direct residue. The
  cfd-2d manifest still declares `num-traits` because direct residues remain
  outside this slice in vorticity-stream, LBM solvers, acoustics, momentum,
  turbulence, physics, and tests.
- [x] Migrate the `cfd-2d::solvers::serpentine_flow` and
  `scalar_transport_2d` scalar-provider surface to Eunomia.
  `serpentine_flow::{mod,solver}` and `scalar_transport_2d` now route scalar
  construction, zero/one identities, Peclet/mixing constants, Fourier-mode
  index conversion, transport relaxation constants, coefficient max/abs
  operations, outlet averaging, validator diagnostics, and module test
  absolute-value checks through `crates/cfd-2d/src/scalar.rs` and
  `eunomia::{FloatElement,NumericElement}` instead of direct
  `num_traits::{Float,FromPrimitive,ToPrimitive}`, `T::from_*`, `T::zero()`,
  `T::one()`, or scalar `.abs()` usage. Verified with `cargo fmt -p cfd-2d
  --check`, `cargo check -p cfd-2d --no-default-features`, `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings`, focused `cargo nextest
  run -p cfd-2d --no-default-features serpentine --status-level fail` (5/5),
  `git diff --check`, and a direct-provider scan over
  `crates/cfd-2d/src/solvers/serpentine_flow` and
  `crates/cfd-2d/src/solvers/scalar_transport_2d.rs` showing no direct
  residue. A `scalar_transport` nextest filter has no matching tests, so
  scalar transport is covered here by compile/clippy and the serpentine tests
  that exercise it. The cfd-2d manifest still declares `num-traits` because
  direct residues remain outside this slice in vorticity-stream, LBM solvers,
  acoustics, momentum, turbulence,
  physics, and tests.
- [x] Migrate the `cfd-2d::solvers::fdm` scalar-provider surface to Eunomia.
  `fdm::{advection_diffusion,config,diffusion,linear_solver,poisson}` now
  route scalar construction, zero/one identities, singular-diagonal guards,
  Gauss-Seidel residual magnitudes, finite-difference constants, and MMS test
  absolute-value checks through `crates/cfd-2d/src/scalar.rs` and
  `eunomia::{FloatElement,NumericElement}` instead of direct
  `num_traits::FromPrimitive`, `T::from_*`, `T::zero()`, `T::one()`, or scalar
  `.abs()` usage. The existing FDM MMS source tests are now wired under
  `#[cfg(test)]`, and the advection-diffusion stencil sign now matches the
  documented steady operator `u dot grad(phi) - alpha laplacian(phi) = S`.
  Verified with `cargo fmt -p cfd-2d --check`, `cargo check -p cfd-2d
  --no-default-features`, `cargo clippy -p cfd-2d --no-default-features --lib
  -- -D warnings`, focused `cargo nextest run -p cfd-2d --no-default-features
  fdm --status-level fail` (2/2), `git diff --check`, and a direct-provider
  scan over `crates/cfd-2d/src/solvers/fdm` showing no direct residue. The
  cfd-2d manifest still declares `num-traits` because direct residues remain
  outside this slice in vorticity-stream, LBM solvers, acoustics, momentum,
  turbulence, physics, and tests.
- [x] Migrate the `cfd-2d::network` scalar-provider surface to Eunomia.
  `network::{build,channel,coupled,postprocess,projection,reference,solve,
  types}` now route scalar construction, zero/one identities, abs/min/max,
  finite checks, Anderson relaxation constants, projection summaries,
  channel diagnostics, and f64 reporting through `crates/cfd-2d/src/scalar.rs`
  and `eunomia::{FloatElement,NumericElement}` instead of direct
  `num_traits::{Float,FromPrimitive,ToPrimitive}`, `T::from_*`, `.to_f64()`,
  `T::zero()`, or `T::one()` usage. Verified with `cargo fmt -p cfd-2d
  --check`, `cargo check -p cfd-2d --no-default-features`, `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings`, focused `cargo nextest
  run -p cfd-2d --no-default-features network --status-level fail` (22/22),
  `git diff --check`, and a direct-provider scan over
  `crates/cfd-2d/src/network` showing no direct residue. The cfd-2d manifest
  still declares `num-traits` because direct residues remain outside this
  slice in drift-diffusion, LBM/FVM, turbulence, physics, and other
  solver/test surfaces.
- [x] Migrate the `cfd-2d::solvers::ns_fvm` scalar-provider surface to
  Eunomia. `field`, `solver::{mod,velocity_interpolation}`,
  `solver::momentum::{u_equation,v_equation}`, and
  `solver::pressure::{poisson,rhie_chow}` now route scalar construction,
  zero/one identities, abs/min/max/sqrt/finite checks, interpolation
  grid-index conversion, pressure-correction thresholds, turbulence seeds, and
  f64 diagnostics through `crates/cfd-2d/src/scalar.rs` and
  `eunomia::{FloatElement,NumericElement}` instead of direct
  `num_traits::{Float,FromPrimitive,ToPrimitive}`, `T::from_*`, `.to_f64()`,
  `T::zero()`, or `T::one()` usage. Verified with `cargo fmt -p cfd-2d
  --check`, `cargo check -p cfd-2d --no-default-features`, `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings`, focused `cargo nextest
  run -p cfd-2d --no-default-features ns_fvm --status-level fail` (2/2),
  `git diff --check`, and a direct-provider scan over
  `crates/cfd-2d/src/solvers/ns_fvm` showing no direct residue. The cfd-2d
  manifest still declares `num-traits` because direct residues remain outside
  this slice in network, turbulence, LBM/FDM/FVM, solver, physics, and other
  test surfaces.
- [x] Migrate the `cfd-2d::solvers::continuity` scalar-provider surface to
  Eunomia. The shared continuity-residual SSOT now routes zero identities,
  central-difference constants, and absolute-value reductions through
  `crates/cfd-2d/src/scalar.rs` and `eunomia::{FloatElement,NumericElement}`
  instead of direct `num_traits::FromPrimitive`, `T::from_*`, `T::zero()`, or
  scalar `.abs()` usage. Verified with `cargo fmt -p cfd-2d --check`, `cargo
  check -p cfd-2d --no-default-features`, `cargo clippy -p cfd-2d
  --no-default-features --lib -- -D warnings`, focused `cargo nextest run -p
  cfd-2d --no-default-features continuity --status-level fail` (6/6), and a
  direct-provider scan over `crates/cfd-2d/src/solvers/continuity.rs` showing
  no `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_*`, `<T as
  FromPrimitive>`, `.to_f64()`, `T::zero()`, `T::one()`, or method `.abs()`
  residue. `crates/cfd-2d/Cargo.toml` still declares `num-traits` because
  direct residues remain outside this slice in network, turbulence,
  LBM/FDM/FVM/NS-FVM, solver, physics, and other test surfaces.
- [x] Migrate the `cfd-2d::simplec_pimple` scalar-provider surface to
  Eunomia. `simplec_pimple::{config,algorithms,diagnostics,solver,
  interpolation,simplec,pimple}` and
  `tests/simplec_pimple_validation.rs` now route defaults, validation
  bounds, residual identities, adaptive-step constants, pressure-correction
  caches, Rhie-Chow face arrays, pressure extrapolation averages, boundary
  zeroing, PIMPLE unity relaxation, and validation reference interpolation
  through the crate-local Eunomia scalar adapter or local Eunomia-backed test
  helpers instead of direct `num_traits::{FromPrimitive,ToPrimitive}`,
  `T::from_*`, `.to_f64()`, `T::zero()`, or `T::one()` usage. Verified with
  `cargo fmt -p cfd-2d --check`, `cargo check -p cfd-2d
  --no-default-features`, `cargo clippy -p cfd-2d --no-default-features --lib
  -- -D warnings`, focused `cargo nextest run -p cfd-2d
  --no-default-features simplec_pimple --status-level fail` (4/4), and
  direct-provider scans over the touched module/test showing no `num_traits`,
  `FromPrimitive`, `ToPrimitive`, `T::from_*`, `<T as FromPrimitive>`,
  `.to_f64()`, `T::zero()`, `T::one()`, or `unwrap_or(T::one())` residue.
  `crates/cfd-2d/Cargo.toml` still declares `num-traits` because direct
  residues remain outside this slice in network, turbulence, LBM/FDM/FVM/
  NS-FVM, solver, physics, and other test surfaces.
- [x] Migrate the `cfd-2d::piso_algorithm` scalar-provider surface to Eunomia.
  `piso_algorithm::{config,convergence,predictor,corrector,solver}` now route
  PISO defaults, tolerances, residual fallbacks, pressure-correction
  workspaces, predictor workspaces, relaxation identities, diffusion constants,
  hybrid differencing max/abs operations, Rhie-Chow tiny thresholds, duration
  thresholds, and logging conversions through the crate-local Eunomia scalar
  adapter and `eunomia::{FloatElement,NumericElement}` instead of direct
  `num_traits::{FromPrimitive,ToPrimitive}`, `T::from_*`, `.to_f64()`,
  `T::zero()`, or `T::one()` usage. Verified with `cargo fmt -p cfd-2d
  --check`, `cargo check -p cfd-2d --no-default-features`, `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings`, focused `cargo nextest
  run -p cfd-2d --no-default-features piso --status-level fail` (6/6), and a
  direct-provider scan over `crates/cfd-2d/src/piso_algorithm` showing no
  `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_*`,
  `<T as FromPrimitive>`, `.to_f64()`, `T::zero()`, `T::one()`, or
  `unwrap_or(T::one())` residue. `crates/cfd-2d/Cargo.toml` still declares
  `num-traits` because broad direct residues remain outside this slice in
  network, turbulence, LBM/FDM/FVM/NS-FVM, solver, physics, and test surfaces.
- [x] Migrate the `cfd-2d::pressure_velocity` scalar-provider surface to
  Eunomia. `pressure_velocity::{coefficients,config,solver,rhie_chow,faces,
  correction,pressure}` now route zero/one identities, pressure-correction
  workspaces, coefficient defaults, validation bounds, Rhie-Chow coefficient
  guards, and pressure-correction scatter values through the crate-local
  Eunomia scalar adapter instead of direct `T::zero()`/`T::one()`
  construction. Finite checks in the touched pressure-velocity paths use
  `eunomia::NumericElement` where the code previously relied on scalar
  methods. Verified with `cargo fmt -p cfd-2d --check`, `cargo check -p
  cfd-2d --no-default-features`, `cargo clippy -p cfd-2d
  --no-default-features --lib -- -D warnings`, focused `cargo nextest run -p
  cfd-2d --no-default-features pressure_velocity --status-level fail`
  (16/16), and a direct-provider scan over `crates/cfd-2d/src/pressure_velocity`
  showing no `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_*`,
  `.to_f64()`, `T::zero()`, or `T::one()` residue. `crates/cfd-2d/Cargo.toml`
  still declares `num-traits` because broad direct residues remain outside
  this slice in network, turbulence, LBM/FDM/FVM/NS-FVM, solver, physics, and
  test surfaces.
- [x] Migrate the `cfd-2d::schemes::time` scalar-provider surface to Eunomia.
  `schemes::time::{explicit,implicit,multistep,adaptive}` now route RK,
  implicit fixed-point, BDF/Adams-Bashforth, CFL, Richardson-error, and
  adaptive-step scalar constants through the crate-local Eunomia scalar
  adapter instead of direct `num_traits::{FromPrimitive,ToPrimitive}` bounds,
  `T::from_f64`, `T::zero()`, `T::one()`, or fallback `unwrap_or(T::one())`
  construction. Ambiguous abs/min/max/powf operations use explicit
  `eunomia::{NumericElement,FloatElement}` dispatch. Verified with `cargo
  fmt -p cfd-2d --check`, `cargo check -p cfd-2d --no-default-features`,
  `cargo clippy -p cfd-2d --no-default-features --lib -- -D warnings`,
  focused `cargo nextest run -p cfd-2d --no-default-features time
  --status-level fail` (29/29), and a direct-provider scan over
  `crates/cfd-2d/src/schemes/time` showing no `num_traits`, `FromPrimitive`,
  `ToPrimitive`, `T::from_*`, `<T as FromPrimitive>`, `.to_f64()`,
  `T::zero()`, `T::one()`, or `unwrap_or(T::one())` residue.
  `crates/cfd-2d/Cargo.toml` still declares `num-traits` because broad direct
  residues remain outside this slice in SIMPLEC/PIMPLE, network, turbulence,
  LBM/FDM/FVM/NS-FVM, solver, physics, and test surfaces.
- [x] Migrate the `cfd-2d::discretization` scalar-provider surface to Eunomia.
  `discretization::{convection,extended_stencil}` and
  `physics::momentum::coefficient_corrections::quick` now route convection
  coefficients, QUICK/MUSCL extended-stencil constants, limiter zero/one
  values, absolute values, and max operations through the crate-local Eunomia
  scalar adapter instead of direct `num_traits::FromPrimitive` construction or
  `T::zero()`/`T::one()` scalar constructors. Verified with `cargo fmt -p
  cfd-2d --check`, `cargo check -p cfd-2d --no-default-features`, `cargo
  clippy -p cfd-2d --no-default-features --lib -- -D warnings`, focused
  `cargo nextest run -p cfd-2d --no-default-features discretization
  --status-level fail` (8/8), and a direct-provider scan over the touched
  files showing no `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_*`,
  `<T as FromPrimitive>`, `.to_f64()`, `T::zero()`, or `T::one()` residue.
  `crates/cfd-2d/Cargo.toml` still declares `num-traits` because direct
  residues remain outside this slice in network, PISO, LBM, turbulence,
  NS-FVM, scalar transport, solver, and test surfaces.
- [x] Migrate the `cfd-2d` spatial/WENO scalar-provider surface to Eunomia.
  `schemes::{constants,central,grid,upwind,weno_helpers,weno,weno_z}`,
  `schemes::tvd::{mod,quick}`, and
  `physics::momentum::{coefficients,coefficient_corrections::weno_z}` now
  route scalar constants, zero values, WENO5/WENO9 smoothness/weight powers,
  WENO-Z momentum corrections, central/upwind coefficients, and TVD/QUICK
  limiter constants through the crate-local Eunomia scalar adapter instead of
  direct `num_traits::{FromPrimitive,ToPrimitive}` construction. Verified with
  `cargo fmt -p cfd-2d --check`, `cargo check -p cfd-2d
  --no-default-features`, `cargo clippy -p cfd-2d --no-default-features --lib
  -- -D warnings`, focused `cargo nextest run -p cfd-2d
  --no-default-features weno --status-level fail` (5/5), and a direct-provider
  scan over the touched files showing no `num_traits`, `FromPrimitive`,
  `ToPrimitive`, `T::from_*`, `<T as FromPrimitive>`, or `.to_f64()`
  residue. `crates/cfd-2d/Cargo.toml` still declares `num-traits` because
  direct residues remain outside this slice in discretization, PISO, network,
  LBM, turbulence, NS-FVM, solver, and test surfaces.
- [x] Migrate the `cfd-2d::fields` scalar-provider surface to Eunomia.
  Field constants, `Field2D::zeros`, `SimulationFields::new`,
  `with_fluid`, `reset`, velocity-magnitude reductions, and Reynolds-number
  averaging now route scalar construction, zero/one values, square root, max,
  and grid-size conversion through `crates/cfd-2d/src/scalar.rs` backed by
  `eunomia::{FloatElement,NumericElement}`. The `copy_from` and velocity
  accessor methods remain on a narrower `RealField + Copy` impl so non-scalar
  construction paths are not over-constrained. Verified with `cargo fmt -p
  cfd-2d --check`, `cargo check -p cfd-2d --no-default-features`, `cargo
  clippy -p cfd-2d --no-default-features --lib -- -D warnings`, focused
  `cargo nextest run -p cfd-2d --no-default-features piso_algorithm
  --status-level fail` (6/6), and a direct-provider scan over
  `crates/cfd-2d/src/fields.rs` showing no `num_traits`, `FromPrimitive`,
  `ToPrimitive`, `Float::`, `T::from_*`, `T::zero()`, `T::one()`, or
  `.to_f64()` residue. `crates/cfd-2d/Cargo.toml` still declares
  `num-traits` because broad direct residues remain in schemes, network,
  PISO, LBM, turbulence, and other cfd-2d modules.
- [x] Remove direct `num-traits` ownership from `cfd-3d`. The last root
  `lib.rs` Chebyshev assertions now use inherent `f64::abs`, and
  `crates/cfd-3d/Cargo.toml` no longer declares `num-traits`. Verified with
  `cargo fmt -p cfd-3d --check`, `cargo check -p cfd-3d
  --no-default-features`, focused `cargo nextest run -p cfd-3d
  --no-default-features chebyshev --status-level fail` (20/20), `cargo
  clippy -p cfd-3d --no-default-features --lib -- -D warnings`, and a full
  direct-provider scan over `crates/cfd-3d/{src,tests,examples}` plus
  `crates/cfd-3d/Cargo.toml` showing no direct `num_traits`, `num-traits`,
  `FromPrimitive`, `ToPrimitive`, `num_traits::Float`, `Float::`,
  `from_f64_or_one`, `<T as From<f64>>::from`, `T::zero()`, `T::one()`, or
  `T::from_` residue. `cargo tree -p cfd-3d -e normal -i num-traits` still
  shows transitive provider/nalgebra paths through `approx`, `nalgebra`,
  `gaia`, `half`, `num-complex`, and related upstream dependencies.
- [x] Migrate the `cfd-3d::{bifurcation,venturi}` scalar-provider seams to
  Eunomia. `bifurcation::{analysis,geometry,solver,types,validation}` and
  `venturi::{analysis,solver,types,validation}` no longer import or bound
  direct `num_traits::{Float,FromPrimitive,ToPrimitive}` or call
  `num_traits::Float::*`, `<T as From<f64>>::from`, `T::zero()`,
  `T::one()`, `T::from_*`, or `T::from_f64_or_one` in the touched module
  APIs. Constants, zero/one values, abs/min/max/sqrt/ln/powf/powi/cos,
  usize/f64 construction, pressure coefficients, Picard viscosity-change
  checks, Richardson/GCI metrics, flow extraction, pressure slicing, wall
  shear, and validation thresholds now route through `cfd-3d::scalar` backed
  by Eunomia. Verified with `cargo fmt -p cfd-3d --check`, `cargo check -p
  cfd-3d --no-default-features`, focused `cargo nextest run -p cfd-3d
  --no-default-features bifurcation venturi --status-level fail` (35/35; one
  16.7s nextest-slow test, below the 30s defect threshold), `cargo clippy -p
  cfd-3d --no-default-features --lib -- -D warnings`, and a targeted scan
  over `crates/cfd-3d/src/{bifurcation,venturi}` plus
  `crates/cfd-3d/src/scalar.rs` showing no direct `num_traits`,
  `FromPrimitive`, `ToPrimitive`, `num_traits::Float`, `Float::`,
  `from_f64_or_one`, `<T as From<f64>>::from`, `T::zero()`, `T::one()`, or
  `T::from_` residue. Direct cfd-3d `num-traits` ownership is now closed;
  remaining provider work is transitive nalgebra/geometry/vector/matrix
  replacement.
- [x] Migrate the `cfd-3d::trifurcation` scalar-provider seam to Eunomia.
  `trifurcation::{geometry,solver,validation}` no longer import or bound
  direct `num_traits::{Float,FromPrimitive,ToPrimitive}` or call
  `num_traits::Float::*`, `<T as From<f64>>::from`, `T::zero()`,
  `T::one()`, or `T::from_f64_or_one` for geometry SDF bounds, Murray-law
  checks, solver configuration, Picard viscosity-change checks, boundary-flow
  extraction, pressure drops, wall-shear estimates, or validation thresholds.
  Constants, zero/one values, abs/min/max/sqrt, sin/cos, integer powers, and
  f64 conversion now route through `cfd-3d::scalar` backed by Eunomia.
  Verified with `cargo fmt -p cfd-3d`, `cargo check -p cfd-3d
  --no-default-features`, focused `cargo nextest run -p cfd-3d
  --no-default-features trifurcation --status-level fail` (6/6), `cargo
  clippy -p cfd-3d --no-default-features --lib -- -D warnings`, and a
  targeted scan over `crates/cfd-3d/src/trifurcation` plus
  `crates/cfd-3d/src/scalar.rs` showing no direct `num_traits`,
  `FromPrimitive`, `ToPrimitive`, `num_traits::Float`, `Float::`,
  `from_f64_or_one`, `<T as From<f64>>::from`, `T::zero()`, or `T::one()`
  residue. Direct cfd-3d `num-traits` ownership is now closed; remaining
  provider work is transitive nalgebra/geometry/vector replacement.
- [x] Migrate the `cfd-3d::serpentine` scalar-provider seam to Eunomia.
  `serpentine::{solver,validation}` no longer import or bound direct
  `num_traits::{Float,FromPrimitive,ToPrimitive}` or call
  `num_traits::Float::*`, `T::zero()`, `T::one()`, or
  `T::from_f64_or_one` for Serpentine configuration, Dean-number analysis,
  pressure-drop validation, Picard viscosity-change checks, or solution
  defaults. Constants, zero/one values, abs/max/sqrt, usize-to-scalar
  construction, and f64 diagnostics now route through `cfd-3d::scalar` backed
  by Eunomia. Verified with `cargo fmt -p cfd-3d`, `cargo check -p cfd-3d
  --no-default-features`, focused `cargo nextest run -p cfd-3d
  --no-default-features serpentine --status-level fail` (7/7), `cargo clippy
  -p cfd-3d --no-default-features --lib -- -D warnings`, and a targeted scan
  over `crates/cfd-3d/src/serpentine` plus `crates/cfd-3d/src/scalar.rs`
  showing no direct `num_traits`, `FromPrimitive`, `ToPrimitive`,
  `num_traits::Float`, `Float::`, `T::zero()`, `T::one()`, or
  `from_f64_or_one` residue. Direct cfd-3d `num-traits` ownership is now
  closed; remaining provider work is transitive nalgebra geometry/vector
  replacement.
- [x] Migrate the `cfd-3d::ibm` scalar-provider seam to Eunomia. IBM
  forcing, interpolation kernels, and grid-index conversion no longer import
  or bound direct `num_traits::{FromPrimitive,ToPrimitive}` or call
  `num_traits::Float::*`, `T::zero()`, `T::one()`, `T::from_usize`, or
  `.to_isize()`. Constants, zero/one values, abs/max/sqrt/cos/floor, and
  index-to-scalar construction now route through `cfd-3d::scalar`; floored
  grid-index conversion now returns a typed `InvalidConfiguration` error for
  non-finite or out-of-range coordinates instead of silently falling back to
  zero. Verified with `cargo fmt -p cfd-3d`, `cargo check -p cfd-3d
  --no-default-features`, focused `cargo nextest run -p cfd-3d
  --no-default-features ibm --status-level fail` (12/12), `cargo clippy -p
  cfd-3d --no-default-features --lib -- -D warnings`, and a targeted scan over
  `crates/cfd-3d/src/ibm` plus `crates/cfd-3d/src/scalar.rs` showing no direct
  `num_traits`, `FromPrimitive`, `ToPrimitive`, `num_traits::Float`,
  `T::zero()`, `T::one()`, `T::from_usize`, `.to_isize()`, or `Float::`
  residue. Direct cfd-3d `num-traits` ownership is now closed; remaining
  provider work is transitive nalgebra geometry/vector replacement.
- [x] Migrate the `cfd-3d::spectral` direct scalar-provider residue to
  Eunomia. `spectral::{chebyshev,poisson,solver}` and the co-located
  Chebyshev tests no longer import `num_traits`, `FromPrimitive`,
  `ToPrimitive`, or `num_traits::Float`; scalar constants, signs, zero/one
  values, and generic absolute values now route through `cfd-3d::scalar`
  backed by Eunomia, while f64 analytical-test cos/max calls use inherent f64
  methods. Verified with `cargo fmt -p cfd-3d`, `cargo check -p cfd-3d
  --no-default-features`, focused `cargo nextest run -p cfd-3d
  --no-default-features spectral --status-level fail` (41/41), `cargo clippy
  -p cfd-3d --no-default-features --lib -- -D warnings`, and a targeted scan
  over `crates/cfd-3d/src/spectral` plus `crates/cfd-3d/src/scalar.rs`
  showing no direct `num_traits`, `FromPrimitive`, `ToPrimitive`,
  `num_traits::Float`, `T::zero()`, `T::one()`, or `Float::` residue.
  Residual spectral provider work is the nalgebra `DMatrix`/`DVector` and
  `RealField` storage/API surface, which remains for the Leto matrix pass.
- [x] Migrate the `cfd-3d::level_set` scalar-provider seam to Eunomia. The
  WENO5-Z reconstruction, SSPRK3 advection, narrow-band update, and
  reinitialization helpers now use `crates/cfd-3d/src/scalar.rs` backed by
  `eunomia::{FloatElement,NumericElement}` instead of direct
  `num_traits::{FromPrimitive,Float}` bounds, `T::zero()`/`T::one()`, or
  `num_traits::Float::*` dispatch. Verified with `cargo fmt -p cfd-3d
  --check`, `cargo check -p cfd-3d --no-default-features`, focused `cargo
  nextest run -p cfd-3d --no-default-features level_set --status-level fail`
  (13/13), `cargo clippy -p cfd-3d --no-default-features --lib -- -D
  warnings`, and a targeted scan over `crates/cfd-3d/src/level_set` plus
  `crates/cfd-3d/src/scalar.rs` showing no direct `num_traits`,
  `FromPrimitive`, `ToPrimitive`, `num_traits::Float`, `T::zero()`,
  `T::one()`, or `Float::` residue. Direct cfd-3d `num-traits` ownership is
  now closed. `cargo clippy -p cfd-3d
  --no-default-features --all-targets -- -D warnings` remains blocked by
  pre-existing lint debt in unrelated cfd-3d tests and modules, so it is not
  evidence against this slice.
- [x] Remove direct `num-traits` ownership from `cfd-validation`. The remaining
  validation residues in `benchmarks::{cavity,step,poiseuille_bifurcation}`
  and `tests/complex_boundary_mms_validation.rs` now use Eunomia
  `FloatElement` or the crate-local `cfd-validation::scalar` adapter for
  scalar construction, f64 diagnostics, absolute values, powers, and max/min
  dispatch instead of `num_traits::{FromPrimitive,ToPrimitive,Float}` or
  `T::from_*` fallback conversions. `crates/cfd-validation/Cargo.toml` no
  longer declares `num-traits`. Verified with `cargo fmt -p cfd-validation
  --check`, `cargo check -p cfd-validation`, full `cargo nextest run -p
  cfd-validation --status-level fail` (431/431), a focused source/test/manifest
  scan showing no direct `num_traits`, `num-traits`, `FromPrimitive`,
  `ToPrimitive`, `T::from_*`, generic `.to_f64()`, or `num_traits::Float`
  residue under `crates/cfd-validation`, and `cargo tree -p cfd-validation -e
  normal -i num-traits` showing only transitive provider/nalgebra paths.
  Residual validation provider work is now nalgebra/nalgebra-sparse storage,
  Leto/Gaia mesh/vector ownership, and transitive upstream providers rather
  than direct `num-traits` ownership.
- [x] Migrate the cfd-2d pressure-velocity/SIMPLEC scalar-provider boundary
  and the dependent `cfd-validation` 2D blood benchmark wrappers to Eunomia.
  `cfd-2d::{simplec_pimple,pressure_velocity}` plus the core
  `physics::momentum::{solver,solve}` impls no longer expose direct
  `num_traits::{FromPrimitive,ToPrimitive}` bounds in the touched solver
  surface; construction constants, Rhie-Chow defaults, pressure-solver
  tolerances, and pressure extrapolation averaging now route through Eunomia
  `FloatElement`/`NumericElement` or the cfd-2d scalar adapter. The
  validation wrappers `benchmarks::{bifurcation,venturi,trifurcation,
  serpentine}` now route touched constants, f64 diagnostics, square roots,
  powers, and zero values through `cfd-validation::scalar` and no longer name
  direct `num_traits` bounds. Verified with `cargo fmt -p cfd-2d -p
  cfd-validation --check`, `cargo check -p cfd-validation`, focused `cargo
  nextest run -p cfd-2d pressure_velocity simplec_pimple --status-level fail`
  (20/20), full `cargo nextest run -p cfd-2d --status-level fail` (567/567,
  27 skipped), full `cargo nextest run -p cfd-validation --status-level fail`
  (431/431), and targeted scans showing no direct `num_traits`,
  `FromPrimitive`, `ToPrimitive`, `T::from_f64`, `T::from_usize`, or generic
  `.to_f64()` residue in the touched cfd-2d solver cone and no direct
  `num_traits`/`SafeFromF64`/`T::from_*` residue in the four validation
  wrappers.
- [x] Migrate `cfd-validation::benchmarks` 3D scalar-provider residue to the
  Eunomia scalar seam. `BenchmarkConfig::default` now constructs tolerance and
  Reynolds-number defaults through `cfd-validation::scalar`, while
  `benchmarks::threed::{bifurcation,venturi,serpentine}` no longer carry
  direct `num_traits::{Float,FromPrimitive,ToPrimitive}` or
  `cfd_core::conversion::SafeFromF64` bounds in the validation wrappers.
  Bifurcation and Venturi benchmark constants and validation predicates now
  route scalar construction and absolute values through the crate-local
  Eunomia adapter. Verified with `cargo fmt -p cfd-validation --check`,
  `cargo check -p cfd-validation`, focused `cargo nextest run -p
  cfd-validation bifurcation_flow_3d serpentine_flow_3d venturi_flow_3d`
  (3/3), full `cargo nextest run -p cfd-validation --status-level fail`
  (431/431), and a targeted scan over the touched benchmark files showing no
  direct `num_traits`, `SafeFromF64`, `from_f64_or_one`, `T::from_*`, or
  `num_traits::Float::abs` residue. The wrapper bounds still include
  `std::convert::From<f64>` because the current `cfd-3d` solver constructors
  require that upstream contract.
- [x] Migrate `cfd-validation::literature` blood-flow validation and
  `cfd-validation::numerical::linear_solver` scalar thresholds to the Eunomia
  scalar seam. `literature::blood_flow_1d` no longer imports direct
  `num_traits::{FromPrimitive,ToPrimitive}` or uses `num_traits::Float` for
  Poiseuille resistance powers, flow/pressure absolute values, scalar
  construction, or f64 diagnostics; these now route through
  `cfd-validation::scalar`, including the new `scalar::powi` adapter backed by
  `FloatElement::powi`. `numerical::linear_solver` now uses `FloatElement`
  bounds and crate-local scalar constants for tolerances and convergence-rate
  metadata instead of direct `num_traits::{Float,FromPrimitive}`. Verified with
  `cargo fmt -p cfd-validation`, `cargo fmt -p cfd-validation --check`,
  `cargo check -p cfd-validation`, focused `cargo nextest run -p
  cfd-validation blood_flow literature chapman patankar` (6/6), full `cargo
  nextest run -p cfd-validation` (431/431), and a combined targeted scan over
  `crates/cfd-validation/src/literature`,
  `crates/cfd-validation/src/numerical/linear_solver.rs`, and
  `crates/cfd-validation/src/scalar.rs` showing no direct `num_traits`,
  `FromPrimitive`, `ToPrimitive`, `num_traits::Float`, `T::from_*`, or generic
  `.to_f64()` residue. Remaining `cfd-validation` provider work is primarily
  benchmarks, tests, nalgebra storage/vector seams, and Leto/Gaia ownership.
- [x] Migrate `cfd-validation::geometry` scalar construction and elementary
  math to the Eunomia scalar seam. `geometry::{serpentine_2d,venturi}` and
  `geometry::threed::serpentine` no longer call direct `T::from_f64`,
  `T::from_usize`, `num_traits::Zero`, or nalgebra/RealField elementary
  methods for the touched scalar constants, sine/cosine, square roots, powers,
  absolute values, and max selection. These now route through
  `cfd-validation::scalar` with explicit `FloatElement` bounds. Verified with
  `cargo fmt -p cfd-validation`, `cargo check -p cfd-validation`, focused
  `cargo nextest run -p cfd-validation geometry` (11/11), full `cargo
  nextest run -p cfd-validation` (431/431), and a geometry-wide residue scan
  showing no targeted direct `SafeFromF64`, `FromPrimitive`, `ToPrimitive`,
  `num_traits`, `ComplexField`, `num_complex`, `T::from_*`,
  `from_f64_or_*`, `try_from_f64`, `nalgebra::scalar`, or generic
  `.to_f64()` residue under `crates/cfd-validation/src/geometry`. The
  remaining geometry work is Leto/Gaia storage and mesh ownership, not direct
  scalar-provider conversion.
- [x] Migrate `cfd-validation::manufactured::richardson` to the Eunomia
  scalar seam. Richardson core extrapolation, MMS validation, and validation
  analysis no longer import direct `num_traits::{FromPrimitive,ToPrimitive}`,
  `num_traits::Float`, `nalgebra::ComplexField`, or test-only
  `nalgebra::scalar` helpers for scalar construction, absolute values,
  powers, logarithms, square roots, finite checks, and f64 diagnostics. These
  operations now route through `cfd-validation::scalar`, backed by Eunomia
  `FloatElement`/`NumericElement`, while the existing `nalgebra::RealField`
  generic storage boundary is preserved for later Leto replacement. Verified
  with `cargo fmt -p cfd-validation`, `cargo check -p cfd-validation`,
  focused `cargo nextest run -p cfd-validation richardson` (16/16), full
  `cargo nextest run -p cfd-validation` (431/431), and a manufactured-wide
  scan showing no direct `SafeFromF64`, `FromPrimitive`, `ToPrimitive`,
  `num_traits`, `ComplexField`, `num_complex`, `T::from_*`,
  `from_f64_or_*`, `try_from_f64`, `nalgebra::scalar`, or generic
  `.to_f64()` residue under `crates/cfd-validation/src/manufactured`.
  Residual direct provider work in `cfd-validation` remains in numerical
  linear-solver, benchmark modules, `literature/blood_flow_1d.rs`,
  literature/storage boundaries, and geometry storage/mesh seams.
- [x] Migrate the `cfd-validation::manufactured` scalar/trig MMS cluster to
  Eunomia. `advection`, `advection_diffusion`, `burgers`, `diffusion`, and
  `navier_stokes` no longer import direct `num_traits::FromPrimitive`,
  cfd-core `SafeFromF64`, or `nalgebra::ComplexField` for scalar
  construction and transcendental math. Constants, square roots, exponentials,
  sine, and cosine now route through `cfd-validation::scalar`, backed by
  Eunomia `FloatElement`/`NumericElement`, while the current `nalgebra`
  `RealField`/`Vector2` storage boundary is preserved for later Leto
  migration. Verified with `cargo fmt -p cfd-validation --check`, `cargo
  check -p cfd-validation`, focused `cargo nextest run -p cfd-validation
  manufactured` (50/50), full `cargo nextest run -p cfd-validation`
  (431/431), and a focused scan over the five touched files showing no direct
  `SafeFromF64`, `FromPrimitive`, `ToPrimitive`, `num_traits`,
  `ComplexField`, `T::from_*`, `from_f64_or_*`, or generic `.to_f64()`
  residue. The Richardson manufactured modules are now closed by the later
  Richardson slice; residual direct provider work in `cfd-validation` remains
  in numerical linear-solver, benchmark modules, `literature/blood_flow_1d.rs`,
  literature/storage boundaries, and broader geometry/storage seams.
- [x] Migrate `cfd-validation::literature::chapman_enskog` to the Eunomia
  scalar seam. Chapman-Enskog transport coefficient construction and report
  formatting no longer import direct `num_traits::{FromPrimitive,ToPrimitive}`
  or call `T::from_f64`/generic `.to_f64()`. Generic scalar conversion now
  routes through the crate-local Eunomia-backed scalar helper with an explicit
  `FloatElement` bound. Added value-semantic tests for the Air 300 K viscosity
  and thermal-conductivity constants plus report scalar values. Verified with
  `cargo fmt -p cfd-validation --check`, focused `cargo nextest run -p
  cfd-validation chapman` (2/2), full `cargo nextest run -p cfd-validation`
  (431/431), and a focused file scan showing no direct `SafeFromF64`,
  `FromPrimitive`, `ToPrimitive`, `num_traits`, `num_complex`,
  `T::from_*`, or generic `.to_f64()` residue. Residual direct provider work
  remains in `cfd-validation` numerical linear-solver, benchmark,
  `literature/blood_flow_1d.rs`, literature storage boundary
  `literature/patankar_1980.rs`, geometry/storage, and remaining
  manufactured/Richardson modules.
- [x] Migrate the `cfd-validation::time_integration` conversion seam to
  Eunomia. `integrators` and `stability_analysis` no longer import direct
  `num_traits`, `FromPrimitive`, `ToPrimitive`, or `SafeFromF64`; RK constants,
  CFL/von-Neumann scalar setup, and display diagnostics now route through the
  crate-local Eunomia scalar helper. `stability_analysis` also uses
  `eunomia::Complex` in the local von-Neumann spatial operators instead of
  `num_complex::Complex`. The current `nalgebra::DVector`/`DMatrix` storage
  boundary is preserved for later Leto migration. Verified with `cargo fmt -p
  cfd-validation`, `cargo check -p cfd-validation`, focused `cargo nextest run
  -p cfd-validation time_integration` (12/12), full `cargo nextest run -p
  cfd-validation` (429/429), and a focused time-integration scan showing no
  direct `SafeFromF64`, `FromPrimitive`, `ToPrimitive`, `num_traits`,
  `num_complex`, `T::from_*`, or `.to_f64().unwrap` residue. Residual direct
  provider work in `cfd-validation` remains in numerical linear-solver,
  benchmark, literature, geometry/storage, and remaining
  manufactured/Richardson modules.
- [x] Migrate the `cfd-validation::error_metrics` scalar metrics group to
  Eunomia. `norms`, `normalized`, and `statistics` no longer import
  `num_traits::FromPrimitive` or call `T::from_f64`/`T::from_usize`; length and
  tolerance construction plus absolute-value and square-root operations now
  route through the crate-local Eunomia scalar helper with explicit
  `FloatElement` bounds. Verified with `cargo fmt -p cfd-validation`, `cargo
  check -p cfd-validation`, focused `cargo nextest run -p cfd-validation
  error_metrics` (21/21), full `cargo nextest run -p cfd-validation`
  (429/429), and a focused error-metrics scan showing no direct
  `num_traits`/`FromPrimitive`/`ToPrimitive`/`T::from_*` residue. Residual
  direct provider work in `cfd-validation` remains in numerical linear-solver,
  benchmark, literature, geometry/storage, and remaining
  manufactured/Richardson modules.
- [x] Close the active `cfd-validation` Eunomia compile blocker group. The
  Patankar SIMPLE literature validation path, advanced-physics MMS, coupled
  multi-physics MMS, Reynolds-stress MMS, and Richardson analysis bound seam
  now use the crate-local Eunomia scalar helper for generic scalar
  construction, transcendental functions, square roots, min/max, and f64
  diagnostics where this slice touched them. `cfd-validation::scalar` now also
  exposes `cosh`/`tanh` so migrated manufactured solutions do not call
  `nalgebra::ComplexField` directly for hyperbolic functions. Verified with
  `cargo fmt -p cfd-validation`, `cargo check -p cfd-validation`, `cargo
  nextest run -p cfd-validation` (429/429), and focused scans showing no
  `SafeFromF64`, `from_f64_or_*`, `try_from_f64`, `T::from_f64`,
  `T::from_usize`, `ComplexField`, or `<T as FromPrimitive>` residue in
  `literature/patankar_1980.rs`, `manufactured/advanced_physics.rs`,
  `manufactured/multi_physics.rs`, and
  `manufactured/reynolds_stress_mms.rs`. Residual direct provider work in
  `cfd-validation` remains in numerical linear-solver, benchmark, literature,
  and remaining manufactured/Richardson modules.
- [x] Remove direct `num-traits` ownership from `cfd-1d` and finish the
  current resistance/vascular Eunomia scalar seam. `ResistanceScalar` now
  requires Eunomia `FloatElement`; resistance geometry/models, serpentine,
  slug-flow, Womersley, Bessel, structured-tree, and bifurcation paths use
  explicit `FloatElement`/`NumericElement` dispatch for scalar construction,
  powers, logarithms, trigonometry, square roots, and absolute values where
  nalgebra is still the storage/linear-algebra boundary. The remaining
  network blueprint/sink bounds, solver-analysis average denominators, and
  package tests were also moved off direct `num_traits`; `crates/cfd-1d/Cargo.toml`
  no longer declares `num-traits`. Verified with `cargo fmt -p cfd-1d
  --check`, `cargo check -p cfd-1d`, `cargo nextest run -p cfd-1d` (725/725,
  3 skipped), and a full cfd-1d source/test/manifest scan showing no direct
  `num_traits`, `num-traits`, `FromPrimitive`, `ToPrimitive`,
  `T::from_f64`, `T::from_usize`, `T::from_u32`, `nalgebra::try_convert`, or
  `.to_f64().unwrap` residue. `cargo tree -p cfd-1d -e normal -i num-traits`
  still shows transitive `num-traits` through `approx`, `nalgebra`,
  `nalgebra-sparse`, `half`/Eunomia/Leto/Hephaestus, `num-complex`, and other
  provider stacks. `cargo clippy -p cfd-1d --all-targets -- -D warnings`
  remains blocked by existing all-target lint debt outside this provider
  cleanup, primarily tests/examples, cell-separation, droplet-regime,
  entrance-model tests, matrix-assembly tests, and venturi coefficient tests.
- [x] Migrate the `cfd-1d` solver-core scalar contract to Eunomia. The shared
  `NetworkSolveScalar` bound no longer names direct `num_traits` traits;
  `anderson_acceleration`, `linear_system`, `convergence`, `solver_detection`,
  and `geometry` route scalar construction, finite checks, absolute values,
  square roots, f64 diagnostics, and solver-method detection through
  `FloatElement`, `NumericElement`, and cfd-core conversion traits. Verified
  with `cargo fmt -p cfd-1d --check`, `cargo check -p cfd-1d`, `cargo nextest
  run -p cfd-1d` (725/725, 3 skipped), and a focused solver-core scan showing
  no direct `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_f64`,
  `nalgebra::try_convert`, or `<T as Float>` residue in the migrated files.
  Residual direct provider work remains in cfd-1d vascular Bessel/Womersley,
  resistance model traits, package tests/benches, and nalgebra/nalgebra-sparse
  storage boundaries.
- [x] Migrate the `cfd-1d` network wrapper scalar/provider seam to Eunomia.
  `EdgeProperties::from`, network characteristic length, Picard resistance
  refresh, blood hematocrit propagation, bifurcation phase-separation bridges,
  coefficient validation, and parallel edge conductance now use
  `SafeFromF64`/`SafeFromUsize` and Eunomia `NumericElement` instead of direct
  `FromPrimitive`, `T::from_f64`, `T::from_usize`, `nalgebra::try_convert`,
  or generic `.abs()` dispatch. Adjacent solver bounds were tightened:
  `NetworkProblem`/`NetworkSolver` use `NetworkSolveScalar`, and
  `MatrixAssembler` no longer requires `FromPrimitive`. Verified with `cargo
  fmt -p cfd-1d --check`, `cargo check -p cfd-1d`, `cargo nextest run -p
  cfd-1d` (725/725, 3 skipped), and a focused scan over
  `wrapper.rs`, `matrix_assembly.rs`, and `problem.rs` showing no direct
  `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_f64`,
  `T::from_usize`, or `nalgebra::try_convert` residue. Residual provider work
  remains in `solver/core/mod.rs`'s `NetworkSolveScalar` compatibility
  contract and broader solver-core/vascular/tests storage seams.
- [x] Migrate the `cfd-1d` network blueprint/sink scalar construction seam to
  the Atlas conversion provider. `network_from_blueprint` now constructs
  blueprint resistances, areas, boundary conditions, serpentine lengths,
  blood defaults, roughness constants, and flow probes through
  `SafeFromF64`/`SafeFromUsize`; the generic quadratic-coefficient policy
  check uses Eunomia `NumericElement::abs`; and `NetworkBuilderSink` carries
  the Atlas conversion bounds required by the canonical blueprint entry point.
  Verified with `cargo fmt -p cfd-1d`, `cargo check -p cfd-1d`, and
  `cargo nextest run -p cfd-1d` (725/725, 3 skipped). Residual direct
  `FromPrimitive` remains in this boundary only because `Network::new` still
  lives in the broader `wrapper.rs` impl requiring it; the next increment is
  the `domain/network/wrapper.rs` provider cleanup.
- [x] Migrate the `cfd-1d` domain-components provider seam from direct
  `num_traits` bounds to Eunomia. Component trait pressure-drop math now uses
  `NumericElement::abs`; channel, membrane, mixer, pump, valve, factory, and
  sensor bounds use `SafeFromF64`/`SafeFromUsize` or `NumericElement` instead
  of `FromPrimitive`/`Float`; and component constants route through the
  provider conversion helpers. Verified with `cargo fmt -p cfd-1d --check`,
  `cargo check -p cfd-1d`, `cargo nextest run -p cfd-1d` (725/725, 3 skipped),
  and a focused component-cone scan showing no direct `num_traits`,
  `FromPrimitive`, `ToPrimitive`, `T::from_f64`, `T::from_usize`, `Float::`,
  or generic `.abs()` residue. `cargo clippy -p cfd-1d --all-targets -- -D
  warnings` still fails on existing all-target lint debt outside this slice
  (tests/examples, domain-network, cell-separation, vascular, solver-core).
- [x] Migrate the next `cfd-1d` channel, branching, and network-analysis
  provider seam from direct `num_traits` bounds to Eunomia. Channel flow
  regime classification, flow resistance, geometry perimeter math,
  Poiseuille shape factors, branching network solver bounds, and
  solver-analysis pressure/flow/resistance/performance aggregates now use
  `SafeFromF64`, `FloatElement`, and `NumericElement` instead of
  `FromPrimitive`, `ToPrimitive`, `num_traits::Float`, `T::from_f64`, or
  generic `to_f64().unwrap_or(...)` patterns. Verified with `cargo fmt -p
  cfd-1d`, `cargo check -p cfd-1d`, `cargo nextest run -p cfd-1d` (725/725,
  3 skipped), a focused channel/branching/analyzer scan showing no direct
  `num_traits` residue, and a focused `solver/analysis` scan showing no direct
  `num_traits` residue. `cargo clippy -p cfd-1d --all-targets -- -D warnings`
  remains blocked by existing unrelated all-target lint debt in examples,
  tests, and cell-separation/resistance modules.
- [x] Migrate the `cfd-1d` Murray's-law vascular geometry provider boundary
  from `num_traits::FromPrimitive` to Eunomia. `MurraysLaw` and
  `OptimalBifurcation` now use `FloatElement`/`NumericElement` for scalar
  construction, power operations, absolute value, and inverse cosine, relying
  on Eunomia's existing `FloatElement::acos` contract. The vascular
  bifurcation network builder propagates the required Eunomia bound where it
  calls the Murray's-law API. Verified with `cargo fmt -p cfd-1d --check` and
  a focused source scan showing no `num_traits`, `FromPrimitive`,
  `ToPrimitive`, or unqualified generic `acos`/`powf`/`powi`/`abs` residue in
  `crates/cfd-1d/src/physics/vascular/murrays_law/{law,bifurcation}.rs`.
  The upstream inverse-cosine provider contract was verified with `cargo
  nextest run -p eunomia acos` in `D:\atlas\repos\eunomia` (2/2).
  `cargo check -p cfd-1d` remains blocked by existing unrelated dirty-tree
  errors in other `cfd-1d` modules (`SafeFromF64` bounds, numeric method
  ambiguity, and stale `to_f64().unwrap_or(...)` call sites).
- [x] Remove direct `cfd-math` `num-traits` ownership from the remaining
  package residue. The GPU Laplacian operator now uses Eunomia
  `NumericElement` for scalar extraction, AMG integration-test constants use
  Eunomia `FloatElement`, and the `cfd-math` manifest no longer declares
  `num-traits`. The touched GPU path also propagates Hephaestus GPU
  synchronization errors directly and no longer wraps synchronous GPU work in
  a stale async/Moirai bridge. Verified with `cargo fmt -p cfd-math --check`,
  `cargo check -p cfd-math`, `cargo check -p cfd-math --features gpu`,
  `cargo clippy -p cfd-math --all-targets -- -D warnings`, `cargo clippy -p
  cfd-math --features gpu --all-targets -- -D warnings`, `cargo nextest run
  -p cfd-math` (333/333), and a focused cfd-math residue scan showing no
  direct `num_traits`, `num-traits`, `FromPrimitive`, or `ToPrimitive` hits.
  `num-traits` remains transitive through provider/nalgebra stacks.
- [x] Remove direct `cfd-core` `num-traits` ownership from the touched scalar
  conversion boundary. `management::conversion` and mesh centroid operations
  now construct scalars through Eunomia `FloatElement`, and the `cfd-core`
  manifest no longer declares `num-traits`. Verified with `cargo fmt -p
  cfd-core --check`, `cargo check -p cfd-core`, `cargo clippy -p cfd-core
  --all-targets -- -D warnings`, `cargo nextest run -p cfd-core` (229/229),
  and a focused `cfd-core` scan showing no direct `num_traits`,
  `num-traits`, `FromPrimitive`, or `ToPrimitive` residue. `num-traits`
  remains transitive through upstream provider/nalgebra stacks.
- [x] Remove stale root workspace provider declarations for direct `wgpu 0.19`
  and `num-complex`. The root workspace no longer advertises unowned direct
  WGPU or num-complex providers; GPU routing remains through
  `hephaestus-wgpu`, and complex values at migrated CFDrs call sites remain
  owned by Eunomia/Apollo provider APIs. Verified with `cargo metadata
  --no-deps --format-version 1`, `cargo check -p cfd-suite`, a manifest scan
  showing no direct `wgpu`, `dep:wgpu`, `wgpu.workspace`, `num-complex`, or
  `num_complex` declarations under `Cargo.toml`/`crates/*/Cargo.toml`, and
  `cargo tree --workspace -e normal -i wgpu@0.19.4` returning no matching
  package. `num-complex` still resolves transitively through provider/nalgebra
  stacks and remains a broader Leto/Eunomia replacement residual.
- [x] Migrate the `cfd-core` GPU crate-level WGPU ABI to Hephaestus. The
  `gpu` feature now depends only on `hephaestus-wgpu`; raw descriptor, buffer,
  kernel, pipeline, Poisson, turbulence, and GPU performance-test WGPU symbols
  resolve through `hephaestus_wgpu::wgpu`; and `GpuContext` acquires devices
  through `hephaestus_wgpu::WgpuDevice` instead of owning adapter selection
  directly. Verified with `cargo fmt -p cfd-core --check`, `cargo check -p
  cfd-core`, `cargo clippy -p cfd-core --all-targets -- -D warnings`,
  `cargo nextest run -p cfd-core` (228/228), `cargo tree -p cfd-core -e
  normal -i wgpu@0.19.4` returning no matching package, and `cargo tree -p
  cfd-core -e normal -i wgpu@26.0.1` showing `wgpu -> hephaestus-wgpu ->
  cfd-core`.
- [x] Remove the unused direct `cfd-optim` nalgebra dev-dependency and clean
  the package gate without widening scope. `cfd-optim` source, tests, and
  manifest now have no direct `nalgebra`, `ndarray`, `num_traits`, `rayon`, or
  `tokio` textual residue; deterministic ranking uses `std::slice::from_ref`,
  clippy-sensitive test modules sit after items, and computed SDT acoustic
  energy/contrast values propagate into `SdtMetrics` instead of remaining as
  underscore-prefixed unused fields. Verified with `cargo fmt -p cfd-optim
  --check`, `cargo check -p cfd-optim`, `cargo clippy -p cfd-optim
  --all-targets -- -D warnings`, `cargo nextest run -p cfd-optim` (121/121),
  focused source/provider-residue scans, and `cargo tree -p cfd-optim -e
  normal -i nalgebra` showing residual nalgebra only through upstream
  transitive crates.
- [x] Sync `crates/cfd-io/agents.md` with the completed Leto/Eunomia I/O
  provider boundary. The crate-local reference now records that `cfd-io` has no
  direct internal solver/core/math dependencies, uses Leto for dense
  checkpoint/binary payloads, Eunomia for scalar bounds, Consus for optional
  HDF5, and RITK for optional VTK. Verified with `cargo fmt -p cfd-io
  --check`, `cargo check -p cfd-io`, `cargo check -p cfd-io --all-features`,
  `cargo clippy -p cfd-io --all-targets -- -D warnings`, `cargo nextest run
  -p cfd-io` (3/3), `cargo tree -p cfd-io -e normal -i nalgebra` (no package
  match), and focused cfd-io provider-residue scans.
- [x] Close the residual `cfd-io` checkpoint-validator scalar constructor
  residue. `check_mass_conservation` now routes domain spacing and the central
  difference factor through Eunomia `FloatElement` instead of direct
  `T::from_f64` calls, while preserving the existing exact-grid-dimension
  validation. Verified with `cargo fmt -p cfd-io --check`, `cargo check -p
  cfd-io`, `cargo clippy -p cfd-io --all-targets -- -D warnings`, `cargo
  nextest run -p cfd-io` (3/3), `cargo test --doc -p cfd-io`, `cargo doc -p
  cfd-io --no-deps`, and a cfd-io source/test/manifest scan showing no direct
  `num_traits`, `num-traits`, `FromPrimitive`, `ToPrimitive`, `T::from_f64`,
  `T::from_usize`, `T::from_u32`, `nalgebra::try_convert`, `ndarray`, or
  `nalgebra` residue.
- [x] Replace `cfd-3d::fem::StokesFlowSolution` velocity/pressure storage from
  nalgebra `DVector` to a Leto-backed `FemDofVector`. Solver and projection
  paths now convert explicitly at the current sparse linear-solver `DVector`
  boundary, Anderson acceleration consumes/returns the Leto-backed FEM DOF
  vector, and the Venturi Picard convergence check uses a checked
  provider-owned difference norm instead of direct `DVector` subtraction.
  Verified with `cargo fmt -p cfd-3d --check`, `cargo check -p cfd-3d`,
  `cargo clippy -p cfd-3d --lib -- -D warnings`, `cargo nextest run -p
  cfd-3d fem` (40/40), and focused scans showing no `DVector` velocity or
  pressure fields remain on `StokesFlowSolution`. Residual FEM provider work
  remains in element `DMatrix`/`DVector` internals and solver/projection sparse
  linear-system `DVector` boundaries.
- [x] Replace the `cfd-3d::fem` solver/projection/solution scalar-identity
  residue with the FEM-local Eunomia scalar SSOT. `FemSolver`,
  `ProjectionSolver`, and `StokesFlowSolution::blend` now route direct
  additive/multiplicative identities through `fem::scalar` while preserving the
  current nalgebra `DVector`/matrix storage boundary for the later Leto FEM
  storage slice. Verified with `cargo fmt -p cfd-3d --check`, `cargo check -p
  cfd-3d`, `cargo clippy -p cfd-3d --lib -- -D warnings`, `cargo nextest run
  -p cfd-3d fem` (40/40), and a focused residue scan over `solver.rs`,
  `projection_solver.rs`, and `solution.rs`. Residual FEM provider work now
  remains in nalgebra-backed FEM matrix/vector storage and the broader
  solver/projection Leto storage boundary.
- [x] Replace the `cfd-3d::fem::element` scalar construction/math residue with
  the FEM-local Eunomia scalar SSOT while preserving the current nalgebra
  matrix-storage boundary for a later Leto slice. `ElementMatrices` and
  `FluidElement` now use Eunomia `FloatElement`/`NumericElement` plus
  `fem::scalar` for zero initialization, volume absolute value, degeneracy
  tolerance, half factors, and sixth-volume constants. Removed direct
  `num_traits::{Float, FromPrimitive}`, `num_traits::Float::abs`,
  exact-half/constant `expect` construction, and old `T::zero` identities from
  the touched element cone. Verified with `cargo fmt -p cfd-3d --check`,
  `cargo check -p cfd-3d`, `cargo clippy -p cfd-3d --lib -- -D warnings`,
  `cargo nextest run -p cfd-3d element` (5/5), `cargo check -p
  moirai-executor --tests` after a transient upstream build race, and a
  focused residue scan over `element.rs`. Residual FEM provider work remains
  in solver/projection assembly, solution scalar identities, and nalgebra-backed
  FEM matrix/vector storage.
- [x] Replace the `cfd-3d::fem` problem-validation scalar residue with the
  FEM-local Eunomia scalar SSOT. `problem_validation` now uses Eunomia
  `FloatElement`/`NumericElement` for finite checks and positive-value
  guards, and `StokesFlowProblem::validate` carries the explicit Eunomia
  scalar bound required by that validation path. Removed direct
  `num_traits::Float` dispatch and old `T::zero` identity checks from the
  touched validation cone. Verified with `cargo fmt -p cfd-3d --check`,
  `cargo check -p cfd-3d`, `cargo clippy -p cfd-3d --lib -- -D warnings`,
  `cargo nextest run -p cfd-3d validate` (10/10), and a focused residue scan
  over `problem_validation.rs`/`problem.rs`. Residual FEM provider work remains
  in element assembly, solver/projection assembly, solution scalar identities,
  and nalgebra-backed FEM matrix/vector storage.
- [x] Replace the `cfd-3d::fem` mesh-helper scalar residue in P2 extraction
  and mid-node lookup with the FEM-local Eunomia scalar SSOT. `mesh_utils` and
  `mid_node_cache` now use `fem::scalar` plus Eunomia
  `FloatElement`/`NumericElement` for midpoint constants, zero comparison,
  extrema initialization, and component-wise min/max. Removed direct
  `num_traits::{Float, FromPrimitive}`, exact-half `expect` construction, and
  old scalar identity/extrema calls from the touched mesh-helper cone.
  Verified with `cargo fmt -p cfd-3d --check`, `cargo check -p cfd-3d`,
  `cargo clippy -p cfd-3d --lib -- -D warnings`, `cargo nextest run -p
  cfd-3d mesh_utils` (5/5), `cargo nextest run -p cfd-3d mid_node_cache`
  (2/2), and a focused residue scan over the touched helper files. Residual
  FEM provider work remains in problem validation, element assembly,
  solver/projection assembly, solution scalar identities, and nalgebra-backed
  FEM matrix/vector storage.
- [x] Replace the next `cfd-3d::fem` scalar residue in stabilization and
  axial boundary classification with the FEM-local Eunomia scalar SSOT.
  `StabilizationParameters`, element-size helpers, and
  `AxialBoundaryClassifier` now use `fem::scalar` plus Eunomia
  `FloatElement`/`NumericElement` for constants, zero/one identities,
  extrema, absolute values, square roots, and scalar min/max. Removed direct
  `num_traits::{FromPrimitive, Float}`, `T::zero`, `T::one`, and
  `nalgebra::RealField::{min_value,max_value}` fallback residue from the two
  target files. Verified with `cargo fmt -p cfd-3d --check`, `cargo check -p
  cfd-3d`, `cargo clippy -p cfd-3d --lib -- -D warnings`, `cargo nextest run
  -p cfd-3d stabilization` (11/11), `cargo nextest run -p cfd-3d boundary`
  (12/12), and a focused residue scan over the touched files. Residual FEM
  provider work remains in `mesh_utils`, `mid_node_cache`, solver/projection
  assembly, and nalgebra-backed FEM matrix/vector storage.
- [x] Replace the `cfd-3d::fem` constant-heavy scalar construction boundary in
  configuration, quadrature, stress/strain, and P2 shape functions with a
  FEM-local Eunomia scalar SSOT. Added `fem::scalar` and removed direct
  `num_traits::{FromPrimitive, Float}` construction/dispatch, old
  `T::from_f64`, `T::zero`, `T::one`, and conversion-fallback residue from
  `config.rs`, `quadrature.rs`, `stress.rs`, and `shape_functions.rs`.
  Propagated the `cfd-1d::NetworkSolveScalar` Eunomia `RealField` bound into
  `cfd-2d::network::reference` and `cfd-validation::literature::blood_flow_1d`
  so cfd-3d dev/test builds compile against the current provider contract.
  Verified with `cargo check -p cfd-3d`, `cargo check -p cfd-validation`,
  `cargo clippy -p cfd-3d --lib -- -D warnings`, `cargo nextest run -p
  cfd-3d fem` (40/40), `cargo nextest run -p cfd-validation blood_flow`
  (2/2), `cargo fmt -p cfd-2d -p cfd-3d -p cfd-validation --check`, and a
  focused residue scan over the touched FEM scalar cone. `cargo clippy -p
  cfd-validation --lib -- -D warnings` was cancelled after prolonged shared
  build-lock contention; `cargo check -p cfd-validation` and the focused
  blood-flow tests passed. Residual provider work remains in FEM nalgebra
  matrix/vector storage, stabilization/mesh/boundary scalar residue, and the
  upstream cfd-1d `NetworkSolveScalar` dependence on num-traits.
- [x] Replace the `cfd-1d` transient composition/droplet scalar conversion
  surface from direct `num_traits` construction and `Float` dispatch to
  Eunomia. `SimulationTimeConfig`, `BloodEdgeTransportConfig`,
  `MixtureComposition`, segmented blood transport, and transient droplet split
  policy now use `FloatElement`/`NumericElement` for constants, absolute
  values, exponentials, clamping, and count-to-scalar conversion. Verified with
  `cargo fmt -p cfd-1d --check`, `cargo check -p cfd-1d`, `cargo clippy -p
  cfd-1d --lib -- -D warnings`, `cargo nextest run -p cfd-1d --test
  transient_composition_parity` (20/20), `cargo nextest run -p cfd-1d --test
  transient_droplet_parity` (9/9), `cargo nextest run -p cfd-1d --test
  transient_literature_validation` (5/5), and focused scans showing no
  `ToPrimitive`, `FromPrimitive`, `num_traits::Float`, direct `num_traits`,
  `T::from_f64`, `T::from_usize`, direct `to_usize`, or silent substep-count
  fallback residue in the touched transient composition/droplet cone. Residual
  provider work remains in the inherited nalgebra `RealField` scalar boundary
  and the Pries phase-separation algorithm's f64 API contract.
- [x] Replace the `cfd-1d` nonlinear solver's local Anderson least-squares
  implementation with the Leto-backed `cfd-math` accelerator. The 1D Picard
  loop now constructs one `cfd_math::nonlinear_solver::AndersonAccelerator`
  per nonlinear solve, converts only at the current nalgebra `DVector`
  linear-system boundary, and no longer stores duplicate Anderson residual or
  iterate queues in `SolverWorkspace`. `NetworkSolveScalar` records the
  temporary scalar bridge between the nalgebra linear-system boundary and the
  Eunomia scalar surface required by the canonical accelerator. Verified with
  `cargo fmt -p cfd-1d -p cfd-math --check`, `cargo check -p cfd-1d`, `cargo
  nextest run -p cfd-1d solver` (49/49 passed), and focused scans showing no
  local `DMatrix`/LU/SVD Anderson solve or workspace Anderson residual queues.
  Residual provider work remains in broader cfd-1d nalgebra vector/sparse
  linear-system storage and transient composition direct `num_traits`
  conversions.
- [x] Replace the cfd-math Anderson acceleration vector/matrix boundary from
  nalgebra to Leto. `AndersonAccelerator::compute_next`, QR state, history
  buffers, and the normal-equations subproblem now use Leto `Array1`/`Array2`
  plus Eunomia scalar math, and `cfd-2d::network::coupled` now passes Leto
  resistance vectors into the mixer. The 3D bifurcation, serpentine,
  trifurcation, and Venturi Picard solvers route their current FEM `DVector`
  velocity fields through one explicit `cfd-3d::atlas_anderson` conversion
  boundary until FEM solution storage is migrated. Verified with `cargo fmt
  -p cfd-math -p cfd-2d -p cfd-3d --check`, `cargo check -p cfd-math`,
  `cargo check -p cfd-2d`, `cargo check -p cfd-3d`, `cargo clippy -p
  cfd-math --lib -- -D warnings`, `cargo nextest run -p cfd-math anderson`
  (5/5 passed), `cargo nextest run -p cfd-2d network` (22/22 passed), and
  focused residue scans showing no nalgebra/num-traits residue in
  `cfd-math::nonlinear_solver::anderson` and no `DVector` residue in
  `cfd-2d::network::coupled`. Residual provider work remains in cfd-3d FEM
  solution velocity storage, other cfd-math linear/time-stepping `DVector`
  APIs, validation benchmark/LES-DES helpers, and direct GPU replacement
  through Hephaestus.
- [x] Replace `cfd-2d::physics::turbulence::constants_validation` scalar
  helpers from nalgebra/num-traits to Eunomia. `TurbulenceConstantsValidator`,
  `ConstantsValidationResult`, and DNS sensitivity analysis now use Eunomia
  `RealField`/`FloatElement`/`NumericElement` for construction, constants,
  square roots, absolute values, maxima, and display conversion. Direct
  `nalgebra::RealField`, `num_traits::{FromPrimitive, ToPrimitive}`,
  fallible `T::from_f64(...).expect(...)`, and legacy `T::zero()` residue was
  removed from the constants-validation module. Verified with `cargo fmt -p
  cfd-2d --check`, `cargo check -p cfd-2d`, `cargo nextest run -p cfd-2d
  constants` (11/11 passed), `cargo nextest run -p cfd-2d rans` (5/5
  passed), `cargo nextest run -p cfd-2d macroscopic` (4/4 passed), and a
  focused residue scan over `constants_validation`. `cargo clippy -p cfd-2d
  --all-targets -- -D warnings` still fails on broader all-target lint debt
  outside this slice, including example/test field reassignment, manual
  multiple-of checks, map identity, format-arg, and LBM boundary identity-op
  cleanups. Residual provider work remains in validation benchmark/LES-DES
  helpers, the cfd-2d network/cfd-math Anderson `DVector` seam, the adjacent
  `cfd-validation` RSM MMS scalar oracle, and direct GPU replacement through
  Hephaestus.
- [x] Replace the shared `cfd-2d::physics::turbulence` RANS model
  scalar/vector contract from nalgebra/num-traits to Leto/Eunomia.
  `TurbulenceModel<T>` now binds on Eunomia `RealField`; k-epsilon,
  realizable C_mu, k-omega SST, Spalart-Allmaras, wall boundary conditions,
  wall treatment, RANS validation, and SIMPLE turbulence coupling use Eunomia
  scalar construction/math with Leto `Vector2` velocity updates. Verified with
  `cargo fmt -p cfd-2d --check`, `cargo check -p cfd-2d`, `cargo nextest run
  -p cfd-2d k_epsilon` (39/39 passed), `cargo nextest run -p cfd-2d
  k_omega` (19/19 passed), `cargo nextest run -p cfd-2d spalart_allmaras`
  (21/21 passed), `cargo nextest run -p cfd-2d wall` (51/51 passed), `cargo
  nextest run -p cfd-2d rans` (5/5 passed), and focused residue scans.
  `cargo clippy -p cfd-2d --all-targets -- -D warnings` was not completed
  because unrelated active builds held the shared `D:\atlas\target` lock.
  Residual turbulence provider work remains in validation scalar helpers, the
  adjacent `cfd-validation` RSM MMS scalar oracle, and direct GPU replacement
  through Hephaestus.
- [x] Replace native `cfd-2d::physics::turbulence::reynolds_stress` scalar
  helper bounds from nalgebra `RealField` and `num_traits::FromPrimitive` to
  Eunomia `RealField`. RSM model construction, tensor storage, production,
  diffusion, curvature, wall-reflection, pressure-strain, and optimized
  transport math now use Eunomia constants/construction/math (`T::ZERO`,
  `T::ONE`, `T::MAX_VALUE`, `T::from_f64`, `sqrt`, `exp`, `clamp`, and
  `max_scalar`) while preserving the existing Leto `Array2` storage boundary.
  Verified with `cargo fmt -p cfd-2d`, `cargo check -p cfd-2d`, `cargo
  nextest run -p cfd-2d reynolds_stress` (10/10 passed), and a focused
  residue scan showing no `FromPrimitive`, direct `nalgebra::RealField`,
  `T::zero()`, `T::one()`, `T::max_value`, legacy generic `.max(T::ZERO)`,
  `DMatrix`, or tuple-indexing residue in the RSM native boundary. The shared
  RANS/Spalart-Allmaras contract is closed by the later Leto/Eunomia item;
  residual turbulence provider work is tracked in validation scalar helpers
  and Hephaestus GPU replacement.
- [x] Replace `cfd-2d::physics::turbulence::reynolds_stress` matrix storage,
  transport, and validation boundaries from nalgebra `DMatrix` to Leto
  `Array2`. `ReynoldsStressTensor`, initialization/dissipation setup,
  production and diffusion helpers, transport velocity/scalar helpers,
  co-located tests, the comprehensive Reynolds-stress validation fixtures, and
  the adjacent `cfd-validation` MMS L2-error oracle now use provider-owned
  arrays. Verified with `cargo fmt -p cfd-2d -p cfd-validation --check`,
  `cargo check -p cfd-2d`, `cargo nextest run -p cfd-2d reynolds_stress`
  (10/10 passed), `cargo nextest run -p cfd-validation reynolds_stress` (9/9
  passed), and a focused residue scan showing no `DMatrix`,
  `DMatrix::`, nalgebra shape methods, `from_element`, or tuple-indexing
  residue in the migrated RSM boundary. Remaining cfd-2d turbulence provider
  work includes validation scalar helpers, the adjacent `cfd-validation` RSM
  MMS scalar oracle, and direct GPU replacement through Hephaestus.
- [x] Replace `cfd-2d::physics::turbulence::les_smagorinsky::wale`
  velocity-field and scalar-provider boundary from `Field2D<nalgebra::Vector2>`
  plus direct `num_traits::FromPrimitive` constants to Leto `Array2`
  component fields and Eunomia `RealField`/`NumericElement`. WALE SGS
  viscosity and gradient recovery now consume provider-owned velocity
  component arrays, validate mismatched shapes/out-of-bounds indices/invalid
  spacings through typed `cfd_core::error::Result`, and preserve the existing
  second-order one-sided boundary-gradient tests. Verified with `cargo fmt -p
  cfd-2d`, `cargo check -p cfd-2d`, `cargo nextest run -p cfd-2d wale`
  (4/4 passed), and a focused residue scan showing no `Field2D`, `Vector2`,
  `nalgebra`, `num_traits`, `FromPrimitive`, old identity methods, fallible
  scalar-construction residue, or direct method `sqrt` in the touched file.
  Remaining cfd-2d turbulence provider work includes validation scalar
  helpers, the adjacent `cfd-validation` RSM MMS scalar oracle, and direct GPU
  replacement through Hephaestus.
- [x] Replace the shared `cfd-2d::physics::turbulence` LES/DES matrix
  boundary from nalgebra `DMatrix` to Leto `Array2`. `LESTurbulenceModel`,
  Smagorinsky model state, strain/viscosity/dynamic helpers, the GPU helper
  boundary, DES model fields, DES length-scale utilities, and LES/DES
  validation and benchmark call sites now use provider-owned matrix storage.
  Verified with `cargo fmt -p cfd-2d --check`, `cargo check -p cfd-2d`,
  `cargo nextest run -p cfd-2d les_smagorinsky` (43/43 passed), `cargo
  nextest run -p cfd-2d des` (14/14 passed), and a focused residue scan
  showing no `DMatrix`/`DMatrix::`/nalgebra shape-method residue in the
  touched LES/DES files. Remaining cfd-2d turbulence provider work includes
  validation scalar helpers, the adjacent `cfd-validation` RSM MMS scalar
  oracle, and direct GPU replacement through Hephaestus.
- [x] Replace `cfd-2d::physics::turbulence::les_smagorinsky::miles`
  standalone 2x2 velocity-gradient API from nalgebra `DMatrix` and direct
  `num_traits::FromPrimitive` constants to Leto `Array2` and Eunomia
  `RealField`/`NumericElement`/`FloatElement`. MILES config defaults, shock
  detection, numerical flux, applicability validation, and co-located value
  tests now use provider-owned scalar and matrix surfaces. Verified with
  `cargo fmt -p cfd-2d --check`, `cargo check -p cfd-2d`, `cargo nextest run
  -p cfd-2d miles` (7/7 passed), and a focused residue scan showing no
  `DMatrix`, `nalgebra`, `num_traits`, `FromPrimitive`, old identity methods,
  or fallible scalar-construction residue in the touched file. Remaining
  cfd-2d turbulence provider work includes validation scalar helpers, the
  adjacent `cfd-validation` RSM MMS scalar oracle, and direct GPU replacement
  through Hephaestus.
- [x] Replace `cfd-2d::physics::turbulence::les_smagorinsky::{sigma,
  vreman}` standalone 2x2 gradient/stress APIs from nalgebra `DMatrix` and
  direct `num_traits::FromPrimitive` constants to Leto `Array2` and Eunomia
  `RealField`/`NumericElement`. Sigma and Vreman configs now use
  provider-owned constants, viscosity/stress/invariant methods accept and
  return `Array2`, and co-located value tests exercise the provider-native
  seam. Verified with `cargo fmt -p cfd-2d --check`, `cargo check -p
  cfd-2d`, `cargo nextest run -p cfd-2d sigma` (7/7 passed), `cargo nextest
  run -p cfd-2d vreman` (6/6 passed), and a focused residue scan showing no
  `DMatrix`, `nalgebra`, `num_traits`, `FromPrimitive`, old identity methods,
  or fallible scalar-construction residue in the two touched files. Remaining
  cfd-2d turbulence provider work includes validation scalar helpers, the
  adjacent `cfd-validation` RSM MMS scalar oracle, and direct GPU replacement
  through Hephaestus.
- [x] Replace the `cfd-3d::spectral::{poisson,solver}` public Poisson
  storage boundary from nalgebra `DVector` to Leto `Array1`. `PoissonProblem`
  and `SpectralSolution` now own Leto arrays, `PoissonSolver::solve` accepts
  and returns Leto arrays with an explicit source-size check, and the direct
  Poisson validation tests, property tests, domain validation tests,
  robustness zero-RHS test, and spectral Poisson example use provider arrays
  at the public seam. Verified with `cargo fmt -p cfd-3d --check`, `cargo
  check -p cfd-3d`, `cargo nextest run -p cfd-3d --test poisson_validation`
  (7/7 passed), `cargo nextest run -p cfd-3d spectral_poisson` (4/4 passed),
  and `cargo nextest run -p cfd-3d test_poisson_zero_rhs_zero_solution` (1/1
  passed). The residual dense Poisson assembly/LU internals still use
  nalgebra `DMatrix`/`DVector`; Chebyshev remains a separate nalgebra
  migration slice.
- [x] Replace the `cfd-3d::spectral::fourier` public transform seam from
  nalgebra `DVector`/`Complex` to Leto `Array1` and Eunomia `Complex` while
  preserving Apollo-backed execution and CFDrs normalization. `FourierTransform`
  and `SpectralDerivative` now accept and return `leto::Array1` values; the
  direct Fourier validation tests, internal adversarial spectral tests, and
  library Fourier smoke test use Leto arrays and Eunomia complex magnitudes.
  Verified with `cargo fmt -p cfd-3d --check`, `cargo check -p cfd-3d`, and
  `cargo nextest run -p cfd-3d --test fourier_validation` (12/12 passed).
  Focused seam scan reports no `DVector`/nalgebra-complex residue in
  `crates/cfd-3d/src/spectral/fourier.rs`,
  `crates/cfd-3d/tests/fourier_validation.rs`,
  `crates/cfd-3d/src/tests/adversarial_tests.rs`, or the updated library
  smoke test.
- [x] Clear the remaining direct Rayon/Tokio crate-source residue and repair
  downstream Eunomia scalar bounds exposed by the provider migration. Stale
  Rayon references in `cfd-optim` progress docs, `cfd-math` docs/agent notes,
  `cfd-3d` FEM/cascade comments, and `cfd-validation` scaling comments now
  point to Moirai or remove the obsolete Rayon guard wording; a full
  `crates/` scan now reports no `rayon`, `tokio`, `par_iter`,
  `par_iter_mut`, `into_par_iter`, `par_chunks`, or `par_chunks_mut` hits.
  `cfd-2d::network::reference` now constructs/extracts reference-solve
  scalars through Eunomia `FloatElement`/`NumericElement`, and
  `cfd-validation::literature::blood_flow_1d` uses Eunomia scalar helpers for
  Carreau-Yasuda blood-flow constants and diagnostic f64 extraction. Verified
  with `cargo fmt -p cfd-math --check`, `cargo fmt -p cfd-3d --check`,
  `cargo fmt -p cfd-optim --check`, `cargo fmt -p cfd-2d -p
  cfd-validation --check`, `cargo check -p cfd-math`, `cargo check -p
  cfd-validation`, `cargo check -p cfd-optim`, `cargo nextest run -p cfd-2d
  reference_trace` (3/3 passed, 590 skipped), `cargo nextest run -p
  cfd-validation blood_flow` (2/2 passed, 427 skipped), and `cargo nextest
  run -p cfd-optim` (121/121 passed). Focused scalar residue scan over the two
  touched scalar-boundary files reports no `T::from_f64`,
  `num_traits::Zero::zero`, or optional `to_f64().unwrap*` conversions.
- [x] Replace `cfd-core::physics::fluid::{traits,validation,temperature,
  properties,newtonian}` nalgebra-style scalar identity dispatch with Eunomia
  `NumericElement`. `FluidState::mach_number`, fluid validation thresholds,
  polynomial/Arrhenius/Andrade/Sutherland temperature paths,
  `FluidProperties` derived-property checks, and Newtonian/ideal-gas
  validation now use provider-owned `ZERO`/`ONE` identities and
  `NumericElement::sqrt` instead of `T::zero()`, `T::one()`, or method
  `sqrt()`. Verified with `cargo fmt -p cfd-core --check`, `cargo check -p
  cfd-core --no-default-features`, `cargo check -p cfd-core`, `cargo nextest
  run -p cfd-core --no-default-features fluid` (45/45 passed, 153 skipped),
  `cargo check -p cfd-1d`, and `cargo nextest run -p cfd-1d resistance`
  (173/173 passed, 555 skipped). Focused scan over the touched fluid files now
  has no `T::zero()`, `T::one()`, or method `sqrt()` hits.
- [x] Replace `cfd-core::physics::fluid::non_newtonian` direct
  `num_traits::FromPrimitive` scalar construction and fallback conversions
  with Eunomia `FloatElement`/`NumericElement`. Bingham, Casson,
  Carreau-Yasuda, Herschel-Bulkley, and Power-law models now construct
  constants through `FloatElement::from_f64`, use provider-owned zero/one
  identities, and dispatch `sqrt`/`powf`/`exp` through Eunomia instead of
  direct nalgebra/num-traits method paths. Verified with `cargo fmt -p
  cfd-core --check`, `cargo check -p cfd-core --no-default-features`, `cargo
  check -p cfd-core`, `cargo nextest run -p cfd-core --no-default-features
  non_newtonian` (5/5 passed, 193 skipped), `cargo check -p cfd-1d`, and
  `cargo nextest run -p cfd-1d resistance` (173/173 passed, 555 skipped).
  Focused source scan reports no `num_traits`, `FromPrimitive`,
  `T::from_f64`, `T::zero()`, `T::one()`, `num_traits::Zero`,
  `nalgebra::try_convert`, or `unwrap_or` residue in the non-Newtonian module.
- [x] Replace `cfd-1d::physics::resistance` model-local scalar-provider
  imports and conversion fallbacks with the centralized resistance scalar
  facade. Hagen-Poiseuille, Darcy-Weisbach, rectangular, entrance, membrane,
  junction-loss, slug-flow, serpentine, Venturi, calculator, factory, and
  geometry paths now use `ResistanceScalar`, `scalar_from_f64`, and
  `scalar_to_f64` from `models::traits` instead of direct per-model
  `nalgebra::RealField`, `num_traits::FromPrimitive`, `T::from_f64`, or
  `nalgebra::try_convert(...).unwrap_or(...)` patterns. Verified with
  `cargo fmt -p cfd-1d`, `cargo check -p cfd-1d`, and `cargo nextest run -p
  cfd-1d resistance` (173/173 passed, 555 skipped). Focused residue scan now
  reports provider hits only in `physics/resistance/models/traits.rs`;
  complete Eunomia replacement for this module remains blocked by the current
  `cfd_core::physics::fluid::FluidTrait<T>: nalgebra::RealField + Copy`
  contract.
- [x] Replace the `cfd-3d::physics::turbulence` LES/RANS model scalar-provider
  seam for Smagorinsky, Dynamic Smagorinsky, Dynamic Gradient Smagorinsky,
  AMD, Vreman, WALE, Sigma, DES, Spalart-Allmaras, Mixing Length, and shared
  SGS energy helpers. These models now use Eunomia `FloatElement`/
  `NumericElement` constants and math dispatch instead of direct
  `nalgebra::RealField`, `num_traits::FromPrimitive`, or
  `num_traits::Float`; Sigma specifically routes inverse cosine through
  Eunomia `FloatElement::acos`. Verified with `cargo fmt -p cfd-3d`,
  `cargo check -p cfd-3d`, and `cargo nextest run -p cfd-3d turbulence`
  (46/46 passed). Focused scans across the edited turbulence files show no
  direct `nalgebra`, `num_traits`, `FromPrimitive`, `RealField`, old
  `T::one()`/`T::zero()`, or `Float::` residue across
  `crates/cfd-3d/src/physics/turbulence`.
- [x] Replace `cfd-3d::physics::turbulence::{field_ops,filter_ops}` direct
  nalgebra/num-traits scalar contracts with Eunomia provider contracts.
  Derivative, gradient, strain, vorticity, matrix, filter-moment, and resolved
  stress helper operations now use `NumericElement` identities, square root,
  and Leto `Vector3`; box-filter stencil counts convert through checked
  `i32::try_from` plus Eunomia `CastFrom<i32>`. Verified with `cargo fmt -p
  cfd-3d`, `cargo check -p cfd-3d`, `cargo nextest run -p cfd-3d
  turbulence` (46/46 passed), focused residue scans showing no `nalgebra`,
  `num_traits`, `FromPrimitive`, or `RealField` in the two helper files, and
  touched-file `git diff --check` with only LF-to-CRLF normalization warnings.
- [x] Replace shared `cfd-core::physics::fluid_dynamics::VelocityField`
  storage from nalgebra `Vector3` to Leto `geometry::Vector3`, and move
  `FlowField`/RANS/turbulence field bounds to Eunomia `NumericElement`.
  Updated direct spectral DNS/forcing/diagnostics, turbulence field/filter
  helpers, IBM NUFFT construction boundaries, and 3D validation benchmark
  consumers to build Leto-backed velocity components. Analytical and
  marker/probe APIs that still expose nalgebra vectors now convert only at
  the boundary. Verified with rustfmt on touched crates, `cargo check -p
  cfd-core`, `cargo check -p cfd-3d`, `cargo check -p cfd-validation`,
  `cargo nextest run -p cfd-core fluid_dynamics::operations` (3/3 passed),
  `cargo nextest run -p cfd-3d spectral` (41/41 passed), `cargo nextest run
  -p cfd-3d nufft` (2/2 passed), `cargo nextest run -p cfd-validation
  taylor_green` (19/19 passed), `cargo nextest run -p cfd-validation
  forced_turbulence` (2/2 passed), `cargo nextest run -p cfd-validation
  nufft_coupling` (2/2 passed), and focused residue scans showing no direct
  nalgebra `Vector3`/`RealField` residue in the migrated core/spectral files.
- [x] Replace `cfd-core::physics::fluid_dynamics::operations` direct
  `num_traits::FromPrimitive` scalar construction with Eunomia
  `NumericElement`. Vorticity, divergence, kinetic-energy, and enstrophy
  operations now use provider-owned zero/one identities, two, and half
  factors while preserving the current nalgebra `Vector3` storage boundary and
  Moirai parallel-slice execution. Added value-semantic tests for constant
  velocity zero divergence/vorticity, affine central-difference divergence,
  and kinetic-energy half-norm semantics. Verified with rustfmt on the touched
  file, `cargo check -p cfd-core`, `cargo nextest run -p cfd-core
  fluid_dynamics::operations` (3/3 passed), and focused residue scans showing
  no direct `num_traits`, `FromPrimitive`, `T::from_f64`, or old
  `unwrap_or_else` scalar fallback patterns in `operations.rs`.
- [x] Replace `cfd-core::physics::fluid_dynamics::service` pipe-flow scalar
  construction and math dispatch from direct `num_traits::{Float,
  FromPrimitive}` to Eunomia `FloatElement`/`NumericElement`. Reynolds and
  Prandtl helpers no longer require `num_traits::Float`; pipe pressure-drop,
  laminar friction factor, Blasius, Colebrook-White, and Haaland constants,
  powers, square roots, logarithms, and absolute convergence checks now route
  through provider-owned Eunomia APIs. Added a value-semantic laminar
  friction-factor test and preserved the existing Colebrook-White validation.
  Verified with rustfmt on the touched file, `cargo check -p cfd-core`,
  `cargo nextest run -p cfd-core fluid_dynamics::service` (2/2 passed), and
  focused residue scans showing no direct `num_traits`, `FromPrimitive`,
  `Float::`, `T::from_f64`, or old `unwrap_or_else` scalar fallback patterns
  in `service.rs`.
- [x] Replace `cfd-core::physics::fluid_dynamics::flow_regimes` direct
  `nalgebra::RealField`/`num_traits::ToPrimitive` scalar conversion with
  Eunomia `RealField`/`NumericElement`. Reynolds, Mach, and combined
  classification now convert through `NumericElement::to_f64` without silent
  `unwrap_or(0.0)` fallback behavior, and the
  `FluidDynamicsService::flow_regime` wrapper now accepts the migrated Eunomia
  real-scalar contract. Added value-semantic tests for Reynolds thresholds,
  Mach thresholds, and hypersonic Mach priority. Verified with rustfmt on the
  touched files, `cargo check -p cfd-core`, `cargo nextest run -p cfd-core
  flow_regime` (3/3 passed), and focused residue scans showing no direct
  `num_traits`, `ToPrimitive`, direct `nalgebra::RealField`, `.to_f64()`, or
  `unwrap_or(0.0)` residue in `flow_regimes.rs`.
- [x] Replace `cfd-3d::spectral::diagnostics` direct
  `num_traits::{FromPrimitive, ToPrimitive}` scalar conversion with Eunomia
  `NumericElement` on the Apollo/Leto diagnostics path. Kinetic-energy and
  enstrophy spectra now build Leto `Array3<f64>` inputs for Apollo FFTs via
  `NumericElement::to_f64`, removing the obsolete fallible conversion branch
  and preserving the existing value-semantic diagnostics tests. Verified with
  rustfmt on the touched file, `cargo check -p cfd-3d`, `cargo nextest run -p
  cfd-3d diagnostics` (5/5 passed), and focused residue scans showing no
  direct `num_traits`, `FromPrimitive`, `ToPrimitive`, `.to_f64()`, or
  obsolete `component_name` residue in `diagnostics.rs`. The remaining
  `nalgebra::RealField`/`Vector3` use is the current `VelocityField<T>` storage
  boundary owned by `cfd-core`, not local diagnostics conversion logic.
- [x] Close the remaining `cfd-2d` LBM MRT and Carreau-Yasuda scalar-provider
  holdouts. `RelaxationMatrix<T>`, `MrtCollision<T>`, `CarreauYasudaModel<T>`,
  and `CarreauYasudaBgk<T>` now use Eunomia `FloatElement`/`NumericElement`
  plus local scalar helpers instead of direct `nalgebra::RealField`,
  `num_traits::{Float, FromPrimitive}`, `T::from_f64`, `T::zero`, `T::one`,
  or direct `Float::sqrt`/`Float::powf` dispatch. Lattice velocity integer
  conversion routes through Eunomia `CastFrom<i32>`. Verified with rustfmt on
  touched files, `cargo check -p cfd-2d`, `cargo nextest run -p cfd-2d lbm`
  (31/31 passed), `cargo nextest run -p cfd-2d carreau_yasuda` (5/5 passed),
  and focused residue scans across the migrated MRT/Carreau-Yasuda files.
- [x] Replace `cfd-2d` LBM macroscopic extraction, D2Q9 equilibrium,
  collision trait seam, and BGK collision scalar bounds from direct
  `nalgebra::RealField`/`num_traits::FromPrimitive` to Eunomia
  `FloatElement` plus local scalar helpers. LBM density, velocity, pressure,
  stress, kinetic-energy, vorticity, equilibrium weights/constants, and BGK
  viscosity/omega construction now route provider-owned scalar constants
  through `crate::scalar`. MRT and Carreau-Yasuda `CollisionOperator` impls
  were bridged with `FloatElement` only where required to satisfy the migrated
  trait seam; the later MRT/Carreau-Yasuda follow-up closes that residual
  scalar-provider work. Verified with rustfmt on touched files, `cargo check
  -p cfd-2d`, `cargo nextest run -p cfd-2d lbm` (31/31 passed), and focused
  residue scans showing no direct `nalgebra`/`RealField`/`num_traits` imports
  or bounds in `macroscopic.rs`, `lattice.rs`, `collision/traits.rs`, or
  `collision/bgk.rs`. Nextest compile also required a narrow provider-bound
  repair in `simplec_pimple_validation.rs` and a value-equivalent concrete
  `f64::abs` disambiguation in `poiseuille/mod.rs`.
- [x] Replace `cfd-3d::physics::turbulence::wall_functions` direct
  `nalgebra::RealField` and `num_traits::{FromPrimitive, Float}` scalar
  construction/math dispatch with Eunomia `FloatElement`/`NumericElement`.
  Spalding wall-law constants, exponentials, logarithm, square root,
  absolute-value convergence tests, and min/max clamps now use provider-owned
  scalar APIs while preserving the existing formulas and value-semantic tests.
  Verified with rustfmt on the touched file, `cargo check -p cfd-3d`, `cargo
  nextest run -p cfd-3d wall_functions` (3/3 passed), focused residue scans
  showing no direct `nalgebra`/`RealField`/`num_traits` imports or bounds in
  `wall_functions.rs`, and touched-file `git diff --check`. Direct
  `GpuContext` replacement with Hephaestus was inspected but remains a
  separate migration item because CFDrs currently resolves both `wgpu@0.19.4`
  and Hephaestus' `wgpu@26.0.1`, making raw-device field replacement a
  broader WGPU API-version alignment.
- [x] Replace the Apollo-backed `cfd-3d::spectral::fourier` wrapper's local
  scalar conversion bounds from direct `num_traits::{FromPrimitive,
  ToPrimitive}` to Eunomia `FloatElement`/`NumericElement`. This removes
  direct num-traits conversion residue from the Fourier wrapper while keeping
  the existing Apollo FFT + Leto array path and legacy `DVector` public
  boundary. Propagated the required Eunomia `FloatElement` bound through
  downstream provider-migrated blood/cavitation consumers in
  `cfd-2d` LBM/Poiseuille/network code, `cfd-3d` DES turbulence, and
  `cfd-validation` 3D blood benchmarks so the focused 3D gate builds against
  the migrated `cfd-core` contracts. Verified with rustfmt on touched files,
  `cargo check -p cfd-3d`, `cargo check -p cfd-validation`, `cargo nextest
  run -p cfd-3d --test fourier_validation` (12/12 passed), and residue scans
  showing no direct `num_traits`/`FromPrimitive`/`ToPrimitive` in
  `cfd-3d/src/spectral/fourier.rs` or `cfd-3d/src/physics/turbulence/des.rs`.
- [x] Replace `cfd-core::physics::fluid::blood::FahraeuasLindqvist`
  local scalar construction and math dispatch from direct
  `num_traits::FromPrimitive` to Eunomia `FloatElement`/`NumericElement`.
  Plasma-viscosity defaults, significance threshold, Pries/Secomb exponent
  formulas, `μ_45` fit, final relative-viscosity clamp, and tube hematocrit
  correlation now use provider-owned scalar constants, zero/one identities,
  `powf`, `exp`, `abs`, and max dispatch. Verified with touched-file rustfmt,
  `cargo check -p cfd-core`, `cargo nextest run -p cfd-core
  fahraeus_lindqvist` (3/3 passed), `cargo nextest run -p cfd-core blood`
  (24/24 passed), focused file residue scan, broader blood residue scan, and
  touched-file `git diff --check`. `cargo clippy -p cfd-core --all-targets --
  -D warnings` remains blocked only by pre-existing unrelated lints in
  `physics/boundary/applicator.rs` and `physics/fluid_dynamics/rhie_chow.rs`.
  Residual blood-fluid provider work is now the broader fluid trait's inherited
  `RealField` boundary; local blood-model scalar `num_traits` construction is
  closed.
- [x] Replace `cfd-core::physics::fluid::blood::{CassonBlood,
  CarreauYasudaBlood, BloodModel}` direct scalar construction and dispatch
  bounds from `num_traits::FromPrimitive` to Eunomia
  `FloatElement`/`NumericElement`. Casson normal/custom/hematocrit
  constructors, temperature viscosity correction, apparent viscosity square
  roots, validation constants, Carreau-Yasuda constants/powers, and the shared
  `BloodModel` selector now route generic scalar math through Eunomia. Verified
  with touched-file rustfmt, `cargo check -p cfd-core`, `cargo nextest run -p
  cfd-core casson` (12/12 passed), `cargo nextest run -p cfd-core
  carreau_yasuda` (4/4 passed), `cargo nextest run -p cfd-core blood` (24/24
  passed), focused touched-file residue scan, and touched-file `git diff
  --check`. `cargo clippy -p cfd-core --all-targets -- -D warnings` remains
  blocked only by pre-existing unrelated lints in
  `physics/boundary/applicator.rs` and `physics/fluid_dynamics/rhie_chow.rs`.
  Residual blood-fluid scalar holdout remains in Fåhræus-Lindqvist and the
  broader fluid trait's inherited `RealField` boundary.
- [x] Replace `cfd-core::physics::fluid::blood::CrossBlood` local scalar
  construction from direct `num_traits::FromPrimitive` to Eunomia
  `FloatElement`/`NumericElement`. `normal_blood` now builds literature
  constants through provider-owned scalar construction, and the Cross apparent
  viscosity denominator uses Eunomia `powf` plus zero/one identities. The
  surrounding `Fluid<T>` trait still imposes the broader `RealField` boundary,
  so that trait migration remains separate. Verified with touched-file rustfmt,
  `cargo check -p cfd-core`, `cargo nextest run -p cfd-core cross` (1/1
  passed), `cargo nextest run -p cfd-core blood` (24/24 passed), focused
  `cross.rs` residue scan, and touched-file `git diff --check`.
- [x] Close the remaining `cfd-core::physics::cavitation` scalar holdouts.
  `CavitationModel<T>` and `ZgbParams<T>` now use Eunomia
  `FloatElement`/`NumericElement`, and the legacy
  `heterogeneous_inception_threshold_pa` adapter is now honestly `f64`
  instead of a fake generic that converted `T` to `f64` and back. Added
  value-semantic tests for Kunz vaporization/condensation, Schnerr-Sauer zero
  bubble density, ZGB degenerate radius rejection, default model construction,
  and heterogeneous legacy-adapter parity with the equivalent selective input.
  Verified with touched-file rustfmt, `cargo check -p cfd-core`, `cargo
  nextest run -p cfd-core cavitation` (40/40 passed), focused cavitation
  residue scan, and touched-file `git diff --check`. No
  `nalgebra::RealField`, direct `num_traits`, direct `T::from_f64`,
  `T::zero`, `T::one`, `to_subset`, or `try_convert` residue remains under
  `crates/cfd-core/src/physics/cavitation/`.
- [x] Replace `cfd-core::physics::cavitation::nuclei_transport` scalar
  contracts from `nalgebra::RealField` to Eunomia
  `FloatElement`/`NumericElement`. The affine nuclei-adjusted vapor-pressure
  closure, default transport config, dissolution sink, generation source, net
  reaction, diffusion accessor, and 1D exponential dissolution update now use
  provider-owned constants, zero/one identities, and exponential evaluation.
  Verified with touched-file rustfmt, `cargo check -p cfd-core`, `cargo
  nextest run -p cfd-core cavitation` (35/35 passed), focused
  `nuclei_transport.rs` residue scan, and touched-file `git diff --check`.
  Residual cavitation `RealField` usage was then limited to `models.rs` and
  `heterogeneous_nucleation.rs`; the current cavitation closeout removed both
  holdouts.
- [x] Replace `cfd-core::physics::cavitation::VenturiCavitation` from
  `nalgebra::RealField` plus direct `num_traits::FromPrimitive` to Eunomia
  `FloatElement`/`NumericElement`. Throat/outlet continuity, Bernoulli throat
  pressure, cavitation number, pressure recovery, loss coefficient, choked
  flow, Nurick cavity length, closure position, and conical cavity volume now
  use provider-owned constants, integer powers, tangent, absolute value, and
  scalar identities. Verified with touched-file rustfmt, `cargo check -p
  cfd-core`, `cargo nextest run -p cfd-core cavitation` (35/35 passed),
  migrated-cone residue scan, and touched-file `git diff --check`. Residual
  cavitation `RealField` usage was then limited to `models.rs`,
  `nuclei_transport.rs`, and `heterogeneous_nucleation.rs`; Sprint 1.96.88
  subsequently removed the `nuclei_transport.rs` holdout.
- [x] Replace the current touched `cfd-core::physics::cavitation` scalar cone
  from `nalgebra::RealField`/direct `num_traits::FromPrimitive` to Eunomia
  `FloatElement`/`NumericElement`. `RayleighPlesset`,
  `SonoluminescenceEstimate`, `CellularMembraneMechanics`,
  `CellularInjuryProfile`, `CavitationRegimeClassifier`,
  `CavitationRegimeAnalysis`, `CavitationNumber`, and `CavitationDamage` now
  use Eunomia-owned constants, finite checks, sqrt/pow/exp, and scalar
  min/max where needed. Added value-semantic tests for cavitation number,
  pressure recovery, Hammitt erosion, Rayleigh collapse impact pressure, and
  pit-depth scaling. Verified with a focused cavitation residue scan,
  touched-file rustfmt, `cargo check -p cfd-core`, `cargo nextest run -p
  cfd-core cavitation` (35/35 passed), and touched-file `git diff --check`.
  `cargo clippy -p cfd-core --all-targets -- -D warnings` is blocked by
  pre-existing unrelated lints in
  `crates/cfd-core/src/physics/boundary/applicator.rs` and
  `crates/cfd-core/src/physics/fluid_dynamics/rhie_chow.rs`. Residual
  cavitation `RealField` usage remained at that point in `models.rs`,
  `venturi.rs`, `nuclei_transport.rs`, and `heterogeneous_nucleation.rs`;
  Sprint 1.96.87 subsequently removed the `venturi.rs` holdout.
- [x] Replace `cfd-core::physics::hemolysis` calculator and platelet activation
  scalar contracts from `nalgebra::RealField` plus direct `num_traits`
  construction to Eunomia `FloatElement`/`NumericElement`. NIH, MIH,
  hemoglobin release, exposure-time estimation, and platelet activation now
  use provider-owned constants, zero/one values, and exponential evaluation.
  Verified with a focused hemolysis residue scan, touched-file rustfmt,
  `cargo check -p cfd-core`, `cargo nextest run -p cfd-core hemolysis` (9/9
  passed), and touched-file `git diff --check`. Residual `cfd-core`
  `RealField` usage remains in compute, boundary, fluid, cavitation, mesh, and
  fluid-dynamics surfaces for later provider slices.
- [x] Replace the `cfd-2d` grid/FVM scalar and cell-center boundary from
  `nalgebra::RealField`/`nalgebra::Vector2` plus direct num-traits conversion
  to Atlas providers. `Grid2D::cell_center`, `StructuredGrid2D`,
  `UnstructuredGrid2D`, and `AdaptiveGrid2D` now use
  `leto::geometry::Vector2` and Eunomia scalar construction; `FvmSolver` no
  longer inherits `RealField` from the structured-grid contract. Validation and
  FDM/LBM consumers that read grid centers now use Leto vector indexing.
  Verified with `cargo check -p cfd-2d`, `cargo check -p cfd-validation`,
  `cargo nextest run -p cfd-2d --lib grid` (21/21 passed), `cargo nextest run
  -p cfd-2d --lib fvm` (25/25 passed), focused grid/FVM residue scans, touched
  rustfmt, and touched-file `git diff --check`. A direct nextest filter for
  `test_poisson_mms_dirichlet_quadratic_interior_accuracy` currently lists zero
  tests because `tests_poisson_mms.rs` is not included in the library module
  tree; that is not counted as verification for this slice.
- [x] Replace `cfd-2d::solvers::fvm` face geometry and velocity-field vector
  storage from nalgebra `Vector2` to `leto::geometry::Vector2`, adding direct
  `leto` consumption to `cfd-2d`. Leto geometry now exposes the `Vector2`
  provider alias and generic fixed-vector norm/normalization methods used by
  FVM face construction. FVM flux schemes remain on Eunomia scalar contracts
  and no longer need nalgebra vector math. The later grid/FVM provider slice
  removed the `StructuredGrid2D<T>` `RealField` blocker. Verified with
  `cargo check -p leto`, `cargo check -p cfd-2d`,
  `cargo nextest run -p leto
  fixed_vector_norm_and_normalization_are_value_semantic` (1/1 passed),
  `cargo nextest run -p cfd-2d --lib fvm` (25/25 passed), focused FVM residue
  scans, and touched-file `git diff --check`.
- [x] Replace `cfd-2d::solvers::fvm` scalar constants from
  `num_traits::FromPrimitive` to Eunomia `FloatElement`/`NumericElement`.
  `FvmConfig` now stores/builds default grid spacing, timestep, tolerance, CFL,
  relaxation, and diffusion constants through provider-owned scalar
  construction, with value-semantic default coverage. `FvmSolver` now declares
  the same `FloatElement` owner bound where it stores `FvmConfig<T>` and uses
  Eunomia-backed exact-index conversion for face centers. FVM flux schemes now
  use Eunomia scalar constants and reject non-finite diffusion coefficients with
  value-semantic error coverage. `RealField` remains only for the existing
  nalgebra vector/solver/flux math boundary. Ambiguous scalar math introduced
  by the Eunomia/RealField overlap now explicitly dispatches `abs` through
  `NumericElement` and `powf` through `FloatElement`. Focused FVM residue scan
  and touched-file rustfmt are clean. Verified with `cargo check -p cfd-2d`,
  `cargo check -p cfd-3d` (covering the migrated spectral dependency chain),
  and
  `cargo nextest run -p cfd-2d --lib fvm` (25/25 passed).
- [x] Replace `cfd-2d::stability::CFLCalculator` from
  `nalgebra::RealField`/`num_traits::FromPrimitive` to Eunomia
  `FloatElement`/`NumericElement`. CFL advection, diffusion, Peclet, explicit
  Euler, QUICK, and max-stable-dt calculations now use Eunomia constants,
  absolute values, and scalar construction. Verified with focused CFL residue
  scan, touched-file rustfmt, `cargo check -p cfd-2d`, and
  `cargo nextest run -p cfd-2d --lib cfl` (9/9 passed). A broader
  `cargo nextest run -p cfd-2d cfl` compiled unrelated integration tests and
  hit the existing `simplec_pimple_validation` generic-bound debt; it is not
  caused by the CFL migration.
- [x] Replace `cfd-core::physics::material` solid/interface trait and value
  contracts from nalgebra `RealField` to Eunomia `FloatElement`/
  `NumericElement`. `SolidProperties`, `InterfaceProperties`, `ElasticSolid`,
  `WettingProperties`, and `FluidSolidInterface` no longer import or require
  `RealField`; `MaterialDatabase` retains `RealField` only because it still
  stores `Box<dyn Fluid<T>>`. Added value-semantic shear-modulus, material
  constructor, contact-angle, and adhesion-energy tests. Verified with a
  focused material residue scan, touched-file rustfmt, `cargo check -p
  cfd-core`, and `cargo nextest run -p cfd-core material` (4/4 passed).
  Residual material-fluid coupling remains until the fluid trait/database
  migrate.
- [x] Replace `cfd-core::physics::values::Velocity` and
  `PhysicalParameters::gravity` from nalgebra `Vector3` storage to
  `leto::geometry::Vector3`, add the direct `leto` dependency to `cfd-core`,
  and remove the `Velocity`/`PhysicalParameters` `RealField` bound where the
  Leto/Eunomia contracts are sufficient. Extended the Leto provider geometry
  surface with Serde derives so serialized CFDrs value objects keep their
  existing boundary contract. Verified with focused provider-residue scan,
  touched-file rustfmt, `cargo check -p cfd-core`, `cargo nextest run -p
  cfd-core --lib` (183/183 passed), `cargo check -p cfd-2d`, `cargo check -p
  cfd-3d`, and `cargo check -p cfd-validation`. Residual `RealField` remains
  in `ProblemAggregate` through `Domain<T>`/fluid contracts and in
  material/hemolysis/fluid-dynamics surfaces.
- [x] Replace `cfd-core::physics::values` scalar wrappers for temperature,
  pressure, Reynolds number, and generic dimensionless numbers from
  `nalgebra::RealField`/`nalgebra::ComplexField` to Eunomia
  `FloatElement`/`NumericElement`, propagating the required bounds through the
  immediate aggregate structs that store those wrappers. Verified with focused
  scalar-provider residue scan, touched-file rustfmt, `cargo check -p
  cfd-core`, and `cargo nextest run -p cfd-core --lib` (183/183 passed).
  Residual `nalgebra::RealField` remains in `Velocity` because it still owns
  `nalgebra::Vector3`, and in material/hemolysis/aggregate paths that still
  depend on vector or nalgebra-owned contracts.
- [x] Remove `cfd-math`'s direct WGPU dependency by adding
  `GpuContext::synchronize` and `GpuContext::supports_timestamp_queries` to
  `cfd-core`, routing `cfd-math::linear_solver::operators::gpu` metrics through
  those methods, and dropping `dep:wgpu` from the `cfd-math/gpu` feature plus
  the root package `gpu` feature. Verified with a focused WGPU residue scan,
  `rustfmt --edition 2021 --check` on touched Rust files, `cargo check -p
  cfd-math --features gpu`, `cargo nextest run -p cfd-math --features gpu
  --lib` (298/298 passed), `cargo check -p cfd-suite --features gpu`, inverse
  WGPU trees showing both WGPU versions flow through `cfd-core`, and `git diff
  --check`. Residual raw WGPU buffer, pipeline, shader, and command-encoder
  ownership remains in `cfd-core::compute::gpu` for the next Hephaestus
  migration slice.
- [x] Route `cfd-core::compute::ComputeBackend::detect_gpu_support` through the
  Atlas Hephaestus provider by adding optional `hephaestus-wgpu` ownership to
  the `gpu` feature and using `hephaestus_wgpu::WgpuDevice::try_default` for
  device acquisition. Verified with `cargo check -p cfd-core --features gpu`,
  `cargo nextest run -p cfd-core --features gpu --lib` (183/183 passed), and
  `cargo tree --workspace -i hephaestus-wgpu` showing the active provider path
  `hephaestus-wgpu -> cfd-core`. Residual raw WGPU buffer/pipeline/shader
  ownership remains in `cfd-core::compute::gpu` and `cfd-math` operator
  adapters for the next Hephaestus migration slice.
- [x] Remove the stale `sprs` dependency from `cfd-1d` and the unused root
  workspace `ndarray` dependency declaration, eliminating the non-Python active
  `ndarray` graph path (`sprs -> cfd-1d`). Add the explicit cfd-1d Eunomia
  dependency needed by tests that name `FloatElement` and propagate the bound
  into generic water-property test helpers. Verified with touched-file rustfmt,
  `cargo update -p sprs` removing `sprs`, `ndarray v0.17.2`, and related stale
  transitive packages, `cargo tree -p cfd-1d -i ndarray` and `cargo tree -p
  cfd-3d -i ndarray` reporting no matching package, `cargo check -p cfd-1d`,
  and `cargo nextest run -p cfd-1d` (725/725 passed, 3 skipped). Residual
  workspace `ndarray v0.16.1` is only via `numpy -> cfd-python`.
- [x] Replace CFDrs' Apollo consumption path with the side-by-side Atlas Apollo
  provider checkout so active `apollo-fft`/`apollo-nufft` no longer resolve the
  older ndarray-backed Git revision; update the cfd-3d FEM and spectral call
  sites for Eunomia scalar/complex contracts and propagate required
  `FloatElement` bounds through cfd-validation dev-dependency code. Verified
  with touched-file rustfmt, `cargo check -p cfd-3d`, `cargo nextest run -p
  cfd-3d` (394/394 passed; one existing mesh-convergence test reported slow at
  16.394s), `cargo tree -p apollo-fft | rg ndarray` and `cargo tree -p
  apollo-nufft | rg ndarray` returning no matches, and an Apollo source/
  manifest/lock scan returning no `ndarray` hits. Residual `ndarray` in the
  `cfd-3d` graph is via `sprs -> cfd-1d`, not Apollo; Leto sparse/operator
  replacement remains a later provider slice.
- [x] Replace `cfd-math` geometric multigrid scalar construction and transfer
  weights from direct `num_traits::FromPrimitive`,
  `cfd_core::conversion::SafeFromF64`, `T::from_f64(...)`,
  `T::from_usize(...)`, and `T::from_f64_or_one(...)` usage to Eunomia
  `FloatElement`/`NumericElement`; remove the now-stale AMG `FromPrimitive`
  bound/import and add value-semantic Poisson stencil and full-weighting
  restriction tests. Verified with touched-file rustfmt, a multigrid-wide
  residue scan showing no direct `num_traits`, `FromPrimitive`, `ToPrimitive`,
  direct `T::from_f64`, direct `T::from_usize`, fallback, `from_f64_or`,
  `SafeFromF64`, stale `rayon`/`tokio`, `rustfft`, `ndarray`, or
  `num_complex` hits, `cargo check -p cfd-math`, and `cargo nextest run -p
  cfd-math geometric_multigrid poisson_matrix restrict_residual fas_solve amg`
  (11/11 passed).
- [x] Replace `cfd-math` multigrid interpolation scalar conversions and
  quality metrics from direct `num_traits::{FromPrimitive, ToPrimitive}`,
  direct `T::from_f64(...)`, direct `T::from_usize(...)`, direct index
  `as f64` conversion, and `to_f64().unwrap_or(...)` fallback usage to
  Eunomia `FloatElement`/`NumericElement`; remove touched test debug output
  and direct index-distance casts. Verified with touched-file rustfmt, a
  focused static residue scan showing no direct `num_traits`, `FromPrimitive`,
  `ToPrimitive`, direct `T::from_f64`, direct `T::from_usize`, fallback,
  `from_f64_or`, `SafeFromF64`, stale `rayon`, direct `as f64`, or
  `.to_f64()` fallback hits in `multigrid/interpolation.rs`, `cargo check -p
  cfd-math`, and `cargo nextest run -p cfd-math interpolation` (14/14 passed).
- [x] Replace `cfd-math` multigrid coarsening algorithm and quality-analysis
  scalar conversions from direct `num_traits::{FromPrimitive, ToPrimitive}`,
  `T::from_f64(...)`, `T::from_usize(...)`, `T::max_value(...)`, and
  `to_f64().unwrap_or(...)` usage to Eunomia `FloatElement`/`NumericElement`;
  route strength-matrix absolute-value dispatch through Eunomia and add a
  value-semantic strength-matrix grid-connectivity test. Verified with
  touched-file rustfmt, a focused coarsening residue scan showing no direct
  `num_traits`, `FromPrimitive`, `ToPrimitive`, direct `T::from_f64`,
  direct `T::from_usize`, fallback, `from_f64_or`, `SafeFromF64`, or stale
  `rayon` hits in `multigrid/coarsening`, `cargo check -p cfd-math`, and
  `cargo nextest run -p cfd-math coarsening` (10/10 passed).
- [x] Replace `cfd-math` multigrid smoother scalar thresholds, Chebyshev
  eigenvalue/default constants, and immediate AMG smoother-owner constants
  with Eunomia `FloatElement`/`NumericElement`, removing direct
  `T::from_f64(...).unwrap_or_else(...)` fallback usage from
  `linear_solver::preconditioners::multigrid::smoothers` and the touched AMG
  owner paths. Added value-semantic smoother/eigenvalue tests for
  Gauss-Seidel, Jacobi, symmetric Gauss-Seidel, SOR, and Chebyshev. Verified
  with touched-file rustfmt, focused smoother/AMG residue scans, `cargo check
  -p cfd-math`, and `cargo nextest run -p cfd-math
  test_gauss_seidel_smoother test_jacobi_smoother
  test_symmetric_gauss_seidel test_sor_smoother test_chebyshev_smoother`
  (5/5 passed).
- [x] Replace `cfd-math` linear-solver convergence monitor scalar
  conversions from `cfd_core::conversion::SafeFromF64`,
  `T::from_f64(...)`, and `T::from_f64_or(...)` fallback usage to Eunomia
  `FloatElement`; disambiguate convergence-factor `powf` through Eunomia and
  add value-semantic coverage for geometric convergence factor, CG theoretical
  bound, and theoretical-bound validation rejection. Verified with touched-file
  rustfmt, a focused `linear_solver/traits.rs` residue scan showing no direct
  `num_traits`, `FromPrimitive`, `ToPrimitive`, `Float::`, direct
  `T::from_f64`, `T::from_usize`, `SafeFromF64`, `from_f64_or`, fallback, or
  stale `rayon` hits, `cargo check -p cfd-math`, and `cargo nextest run -p
  cfd-math convergence_factor cg_theoretical_bound validate_convergence` (4/4
  passed, including the existing AMG convergence-factor match).
- [x] Replace `cfd-math` Poisson/Laplacian and momentum/energy linear
  operator scalar constants and provider bounds from direct
  `num_traits::FromPrimitive`/`T::from_f64(...).unwrap_or_else(...)` usage to
  Eunomia `FloatElement`; add value-semantic finite-difference operator tests
  for 2D Laplacian, 3D Poisson, 1D momentum, 2D momentum, and 2D energy center
  stencils. Verified with touched-file rustfmt, a focused
  `operators/{poisson.rs,momentum.rs}` residue scan showing no direct
  `num_traits`, `FromPrimitive`, `ToPrimitive`, `Float::`, direct
  `T::from_f64`, `T::from_usize`, fallback, `from_f64_or`, or stale `rayon`
  hits, `cargo check -p cfd-math`, `cargo nextest run -p cfd-math
  laplacian_center_impulse poisson_center_impulse` (2/2 passed), and `cargo
  nextest run -p cfd-math momentum_1d_applies momentum_2d_applies
  energy_2d_applies` (3/3 passed).
- [x] Replace `cfd-math` Schwarz/IncompleteCholesky preconditioner provider
  residue by removing Schwarz's stale direct `num_traits::FromPrimitive`
  bound/import and routing Cholesky symmetry tolerance plus symmetry residual
  absolute-value dispatch through Eunomia `FloatElement`/`NumericElement`;
  preserve existing nalgebra CSR/vector preconditioner surfaces for later Leto
  migration. Verified with touched-file rustfmt, a focused
  `schwarz.rs`/`cholesky.rs` residue scan showing no direct `num_traits`,
  `FromPrimitive`, `ToPrimitive`, `Float::`, direct `T::from_f64`,
  `T::from_usize`, fallback, `from_f64_or`, or stale `rayon` hits, `cargo
  check -p cfd-math`, `cargo nextest run -p cfd-math cholesky` (4/4 passed),
  and `cargo nextest run -p cfd-math schwarz` (1/1 passed).
- [x] Replace `cfd-math` SSOR preconditioner default relaxation construction
  and provider bounds from direct `num_traits::FromPrimitive` and
  `T::from_f64(...).unwrap_or_else(...)` usage to Eunomia `FloatElement`;
  preserve the existing nalgebra CSR/vector preconditioner surface for later
  Leto migration. Verified with touched-file rustfmt, a focused `ssor.rs`
  residue scan showing no direct `num_traits`, `FromPrimitive`, `ToPrimitive`,
  `Float::`, direct `T::from_f64`, `T::from_usize`, fallback,
  `from_f64_or`, or stale `rayon` hits, `cargo check -p cfd-math`, and
  `cargo nextest run -p cfd-math ssor` (4/4 passed).
- [x] Replace `cfd-math` basic Jacobi/SOR preconditioner scalar constants,
  diagonal tolerances, omega construction, and absolute-value dispatch from
  direct `num_traits::FromPrimitive`/`T::from_f64(...).unwrap_or_else(...)`
  usage to Eunomia `FloatElement`/`NumericElement`; preserve the existing
  nalgebra CSR/vector preconditioner surface for later Leto migration.
  Verified with touched-file rustfmt, a focused `basic.rs` residue scan showing
  no direct `num_traits`, `FromPrimitive`, `ToPrimitive`, `Float::`, direct
  `T::from_f64`, `T::from_usize`, fallback, `from_f64_or`, or stale `rayon`
  hits, `cargo check -p cfd-math`, `cargo nextest run -p cfd-math jacobi`
  (5/5 passed), and `cargo nextest run -p cfd-math sor` (6/6 passed).
- [x] Replace `cfd-math` GMRES, linear-solver chain, direct sparse solver, and
  block/SIMPLE preconditioner provider bounds and scalar safeguards from direct
  `num_traits::{Float, FromPrimitive, ToPrimitive}` usage to Eunomia
  `FloatElement`/`NumericElement`; preserve the current nalgebra vector/matrix
  and rsparse f64-backed direct-LU surfaces for later Leto/backend migration.
  Verified with touched-file rustfmt, focused source scan showing no direct
  `num_traits`, `FromPrimitive`, `ToPrimitive`, `Float::`, direct
  `T::from_f64`, `T::from_usize`, fallback, or stale `rayon` hits in the
  touched slice, `cargo check -p cfd-math`, `cargo nextest run -p cfd-math
  gmres` (21/21 passed), `cargo nextest run -p cfd-math direct_solver` (3/3
  passed), and `cargo nextest run -p cfd-math block_preconditioner` (2/2
  passed).
- [x] Replace `cfd-math` iterative linear-solver config default tolerance
  construction from direct `num_traits::FromPrimitive`/
  `T::from_f64(...).expect(...)` usage to Eunomia `FloatElement`, update
  stale parallel SpMV wording from Rayon to Moirai, and propagate
  `FloatElement` default-construction bounds through CG, BiCGSTAB, and GMRES.
  Verified with touched-file rustfmt, `cargo check -p cfd-math`, `cargo
  nextest run -p cfd-math default_solver` (2/2 passed), `cargo nextest run -p
  cfd-math test_gmres_configurable_trait` (1/1 passed), and a focused source
  scan showing no direct `num_traits`, `FromPrimitive`, direct `T::from_f64`,
  or stale `rayon` hits in `linear_solver/config.rs`, `bicgstab/mod.rs`, or
  `conjugate_gradient/mod.rs`.
- [x] Replace `cfd-math` CFD SIMD central-difference constants and field-norm
  square-root dispatch from direct `num_traits::FromPrimitive`/
  `T::from_f64(...).unwrap_or_else(...)` usage to Eunomia
  `FloatElement`/`NumericElement`; preserve the existing Moirai-backed
  parallel slice execution and nalgebra `RealField` API for later Leto scalar
  cleanup. Verified with touched-file rustfmt, `cargo check -p cfd-math`,
  `cargo nextest run -p cfd-math simd` (26/26 passed), and a SIMD source scan
  showing no direct `num_traits`, `FromPrimitive`, `ToPrimitive`, `Signed`,
  `Float::`, direct `T::from_f64`, conversion fallback, `T::epsilon`, or
  `rayon` hits.
- [x] Replace `cfd-math` sparse stencil constants, Frobenius norm dispatch,
  diagonal dominance checks, and condition-estimate absolute-value/threshold
  handling from direct `num_traits::{Float, FromPrimitive, Signed}` usage to
  Eunomia `FloatElement`/`NumericElement`; update sparse parallel SpMV docs to
  name the existing Moirai slice adapter instead of Rayon. Verified with
  touched-file rustfmt, `cargo check -p cfd-math`, `cargo nextest run -p
  cfd-math sparse` (15/15 passed), and a sparse source scan showing no direct
  `num_traits`, `FromPrimitive`, `ToPrimitive`, `Signed`, `Float::`, direct
  `T::from_f64`, conversion fallback, `T::epsilon`, or `rayon` hits.
- [x] Replace `cfd-math` Anderson/JFNK nonlinear solver defaults,
  finite-difference perturbation safeguards, QR norms/diagonal checks, EW
  forcing clamp math, Givens rotations, and back-substitution checks from
  direct `num_traits::{Float, FromPrimitive}` usage to Eunomia
  `FloatElement`/`NumericElement`; preserve the current nalgebra
  `DVector`/`DMatrix` solver API for the later Leto migration. Verified with
  touched-file rustfmt, `cargo check -p cfd-math`, `cargo nextest run -p
  cfd-math nonlinear_solver` (9/9 passed), and a nonlinear-solver source scan
  showing no direct `num_traits`, `FromPrimitive`, `Float::`, direct
  `T::from_f64`, conversion fallback, or `T::epsilon` hits.
- [x] Replace `cfd-math` SIMPLE pressure-velocity defaults from direct
  nalgebra `RealField`, `num_traits::FromPrimitive`,
  `T::from_f64(...).expect(...)`, and old `T::zero()`/`T::one()` identities
  to Eunomia `RealField`/`FloatElement`; `SIMPLEConfig::new` no longer
  requires any scalar-conversion trait. Verified with `cargo fmt -p cfd-math
  --check`, `cargo check -p cfd-math --no-default-features --lib`, `cargo
  nextest run -p cfd-math --no-default-features pressure_velocity
  --status-level fail` (3/3 passed, 341 skipped), `cargo clippy -p cfd-math
  --no-default-features --all-targets -- -D warnings`, and a
  pressure-velocity source scan showing no nalgebra, `DVector`, `DMatrix`, old
  scalar identities, `num_traits`, `num_complex`, ndarray, rayon, tokio,
  rustfft, `FromPrimitive`, `ToPrimitive`, or `From<f64>` hits.
- [x] Replace `cfd-math` iterator stencil coefficients and iterator
  statistics count/math conversions from direct
  `num_traits::FromPrimitive`/`T::from_f64`/`T::from_usize` fallback
  construction to Eunomia `FloatElement`/`NumericElement`; removed the
  previous zero-vector placeholder branch for second-derivative coefficients
  by making every declared 3-point stencil return the real `[1, -2, 1]`
  second-derivative coefficients. Verified with touched-file
  `rustfmt --check`, `cargo check -p cfd-math`, `cargo nextest run -p
  cfd-math iterators` (7/7 passed), and an iterator source scan showing no
  `num_traits`, `FromPrimitive`, `ToPrimitive`, direct `T::from_f64`/
  `T::from_usize`/`T::from_f32`, or conversion-fallback hits.
- [x] Replace `cfd-math` WENO5/WENO7 epsilon defaults, linear weights, ENO
  reconstruction coefficients, and smoothness-indicator constants from direct
  `num_traits::FromPrimitive`/`T::from_f64(...).unwrap_or_else(...)` fallback
  construction to Eunomia `FloatElement`; denominator squaring now uses a
  local multiplication helper to avoid `FloatElement`/nalgebra `powi`
  ambiguity. Verified with touched-file `rustfmt --check`, `cargo check -p
  cfd-math`, `cargo nextest run -p cfd-math weno` (6/6 passed), and a
  `crates/cfd-math/src/high_order` scan showing no `num_traits`,
  `FromPrimitive`, `ToPrimitive`, direct `T::from_f64`/`T::from_usize`/
  `T::from_f32`, or conversion-fallback hits.
- [x] Replace `cfd-2d` scheme amplification complex values from direct
  `num_complex::Complex<f64>` to `eunomia::Complex<f64>`, remove the direct
  `num-complex` manifest dependency, and propagate the required
  `eunomia::FloatElement` bounds through 2D network and solver paths that
  construct `cfd-core` staggered grids or solver configs; verified with
  `cargo fmt --package cfd-2d --package cfd-3d --package cfd-validation
  --check`, `cargo check -p cfd-2d`, `cargo check -p cfd-3d`, `cargo check -p
  cfd-validation`, `cargo nextest run -p cfd-2d` (563/563 passed, 27 skipped),
  and a `crates/cfd-2d` source/manifest scan showing no `num_complex`,
  `num-complex`, or `NumComplex` hits.
- [x] Replace `cfd-math`/`cfd-validation` time-stepping stability analysis
  complex arithmetic from direct `num_complex::Complex<f64>` to
  `eunomia::Complex<f64>`, remove their direct `num-complex` manifest
  dependencies, and replace the explicit RK complex matrix inverse with
  forward substitution over `(I - zA)x = 1`; verified with touched-file
  rustfmt, `cargo check -p cfd-math`, `cargo nextest run -p cfd-math stability`
  (5/5 passed), and a source/manifest scan showing no `num_complex`,
  `num-complex`, or `NumComplex` hits in `crates/cfd-math` or
  `crates/cfd-validation`. The downstream validation compile blocker was later
  cleared by the cfd-2d/cfd-validation `FloatElement` propagation slice.
- [x] Replace `cfd-math` stability analyzer scalar constants and diagnostic
  numeric conversions from direct `num_traits::ToPrimitive` and
  `T::from_f64(...).unwrap_or_else(num_traits::Zero::zero)` fallbacks to
  Eunomia `FloatElement`/`NumericElement` helpers; verified with touched-file
  `rustfmt --check`, `cargo check -p cfd-math`, `cargo nextest run -p cfd-math
  stability` (5/5 passed), `git diff --check`, and a touched stability source
  scan showing no `num_traits`, `ToPrimitive`, `FromPrimitive`, `T::from_f64`,
  or conversion-fallback hits. Package-level `cargo fmt --package cfd-math
  --check` remains blocked by unrelated existing formatting drift outside the
  touched stability files.
- [x] Replace `cfd-math` Runge-Kutta scalar constants from direct
  `T::from_f64(...).unwrap_or_else(num_traits::Zero::zero)` fallbacks to
  Eunomia `FloatElement`, and fix `LowStorageRK4` to use the
  Carpenter-Kennedy 2N residual recurrence instead of mutating the solution
  accumulator incorrectly; verified with touched-file `rustfmt --check`,
  `cargo check -p cfd-math`, `cargo nextest run -p cfd-math runge_kutta`
  (5/5 passed), and a touched source scan showing no `num_traits`,
  `FromPrimitive`, `ToPrimitive`, `T::from_f64`, or conversion-fallback hits.
- [x] Replace `cfd-math` adaptive time-stepper scalar constants and
  Dormand-Prince tableau coefficients from direct
  `T::from_f64(...).unwrap_or_else(num_traits::Zero::zero)` fallbacks to
  Eunomia `FloatElement`; verified with touched-file `rustfmt --check`,
  `cargo check -p cfd-math`, `cargo nextest run -p cfd-math adaptive`
  (3/3 passed), and a touched source scan showing no `num_traits`,
  `FromPrimitive`, `ToPrimitive`, `T::from_f64`, or conversion-fallback hits.
- [x] Replace `cfd-math` Runge-Kutta-Chebyshev scalar constants, stage-count
  conversions, adaptive error-control constants, and transcendental/math
  dispatch from direct `num_traits::FromPrimitive`/nalgebra method calls to
  Eunomia `FloatElement`/`NumericElement`; verified with touched-file
  `rustfmt --check`, `cargo check -p cfd-math`, `cargo nextest run -p cfd-math
  rk_chebyshev` (4/4 passed), and a touched source scan showing no
  `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_f64`,
  `T::from_usize`, or conversion-fallback hits.
- [x] Replace `cfd-math` IMEX Newton tolerance, ARS343 gamma/delta constants,
  tableau coefficients, and solution weights from direct
  `T::from_f64(...).unwrap_or_else(num_traits::Zero::zero)` fallbacks and
  nalgebra scalar `sqrt` dispatch to Eunomia `FloatElement`/`NumericElement`;
  verified with touched-file `rustfmt --check`, `cargo check -p cfd-math`,
  `cargo nextest run -p cfd-math imex` (5/5 passed), and a touched source scan
  showing no `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_f64`,
  `T::from_usize`, or conversion-fallback hits.
- [x] Replace `cfd-math` exponential time-differencing scalar constants, ERK4
  stage weights, phi-function small-argument thresholds, and series factorial
  conversions from direct `num_traits::FromPrimitive`/
  `T::from_f64(...).unwrap_or_else(num_traits::Zero::zero)` fallbacks to
  Eunomia `FloatElement`; verified with touched-file `rustfmt --check`,
  `cargo check -p cfd-math`, `cargo nextest run -p cfd-math
  time_stepping::exponential::` (2/2 passed), and a touched source scan
  showing no `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_f64`,
  `T::from_usize`, or conversion-fallback hits.
- [x] Replace `cfd-math` integration quadrature scalar constants,
  interval-count conversions, adaptive tolerance/error dispatch, and
  tetrahedral quadrature constants from direct `num_traits::FromPrimitive`,
  `From<f64>`, `T::from_f64`, `T::from_usize`, nalgebra `RealField`, and old
  `T::zero()`/`T::one()` identities to Eunomia `RealField`/
  `FloatElement`/`NumericElement`; verified with `cargo fmt -p cfd-math
  --check`, `cargo check -p cfd-math --no-default-features --lib`, `cargo
  nextest run -p cfd-math --no-default-features integration --status-level
  fail` (12/12 passed, 330 skipped), `cargo clippy -p cfd-math
  --no-default-features --all-targets -- -D warnings`, and an integration
  source scan showing no nalgebra, `DVector`, `DMatrix`, old scalar identities,
  `num_traits`, `num_complex`, ndarray, rayon, tokio, or rustfft hits.
- [x] Replace `cfd-math` finite-difference and gradient scalar constants,
  SIMD helper scalar staging, and 2D Laplacian constants from direct
  `num_traits::FromPrimitive`, `From<f64>`, `T::from_f64`, and `T::from_f32`
  construction to Eunomia `FloatElement`/`NumericElement`; follow-up Leto
  migration replaced `DVector` differentiation results with Leto `Array1`,
  gradient vector surfaces with Leto `Vector3`, and the type-suffixed
  `first_derivative_simd_f32` API with `FiniteDifference<f32>::
  first_derivative_simd`. Verified with touched-file `rustfmt --check`,
  cfd-math no-default library check, focused no-default differentiation
  nextest (12/12), no-default all-target clippy, and a differentiation source
  scan showing no nalgebra, `DVector`, old scalar identities, `num_traits`,
  `num_complex`, ndarray, or removed type-suffixed SIMD name.
- [x] Replace `cfd-math` interpolation trait, linear, Lagrange, and
  cubic-spline scalar contracts from nalgebra `RealField`, old
  `T::zero()`/`T::one()` identities, direct `num_traits::FromPrimitive`,
  `T::from_f64`, and conversion-fallback construction to Eunomia
  `RealField`/`FloatElement`; Lagrange now rejects duplicate/non-increasing
  nodes before basis denominator division. Verified with `cargo fmt -p
  cfd-math --check`, `cargo check -p cfd-math --no-default-features --lib`,
  `cargo nextest run -p cfd-math --no-default-features interpolation
  --status-level fail` (15/15 passed, 328 skipped), `cargo clippy -p cfd-math
  --no-default-features --all-targets -- -D warnings`, and an interpolation
  source scan showing no nalgebra, `DVector`, `DMatrix`, old scalar identities,
  `num_traits`, `num_complex`, ndarray, rayon, tokio, rustfft,
  `FromPrimitive`, `ToPrimitive`, or `From<f64>` hits.
- [x] Replace `cfd-core` Rhie-Chow interpolation scalar constants from direct `num_traits::FromPrimitive` fallback conversions to Eunomia `FloatElement`, and add value-semantic pressure-correction tests for both velocity components; verified with `rustfmt --check --edition 2021 crates\cfd-core\src\physics\fluid_dynamics\rhie_chow.rs`, `cargo check -p cfd-core --no-default-features`, `cargo check -p cfd-core`, `cargo nextest run -p cfd-core --no-default-features rhie_chow` (2/2 passed), `cargo nextest run -p cfd-core --no-default-features` (176/176 passed), and a touched-file source scan showing no `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_`, or conversion-fallback hits in `physics/fluid_dynamics/rhie_chow.rs`.
- [x] Replace `cfd-core` boundary geometry measure constants from direct `T::from_f64(...).unwrap_or_else(...)` fallback conversions to Eunomia `FloatElement`, while keeping `contains_point`/`dimension` on their narrower existing scalar bounds; verified with `rustfmt --check --edition 2021 crates\cfd-core\src\physics\boundary\geometry.rs`, `cargo check -p cfd-core --no-default-features`, `cargo check -p cfd-core`, `cargo nextest run -p cfd-core --no-default-features boundary::geometry` (4/4 passed), `cargo nextest run -p cfd-core --no-default-features` (174/174 passed), and a touched-file source scan showing no `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_`, or conversion-fallback hits in `physics/boundary/geometry.rs`.
- [x] Replace `cfd-core` boundary ghost-cell scalar constants and Robin singularity reporting from direct `num_traits::{FromPrimitive, ToPrimitive}` to Eunomia `FloatElement`/`NumericElement`, and reject degenerate Robin coefficients before division; verified with `cargo check -p cfd-core --no-default-features`, `cargo check -p cfd-core`, `cargo nextest run -p cfd-core --no-default-features ghost_cells` (5/5 passed), `cargo nextest run -p cfd-core --no-default-features` (170/170 passed), touched-file `rustfmt --check`, and a source scan showing no `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_`, or conversion-fallback hits in `physics/boundary/ghost_cells.rs`.
- [x] Replace `cfd-core` staggered-grid scalar coordinate conversions from direct `num_traits::FromPrimitive` to Eunomia `FloatElement`, with exact-representability assertions for integer grid indices; verified with `cargo check -p cfd-core --no-default-features`, `cargo check -p cfd-core`, `cargo nextest run -p cfd-core --no-default-features staggered` (5/5 passed), `cargo nextest run -p cfd-core --no-default-features` (167/167 passed), touched-file `rustfmt --check`, `git diff --check`, and a source scan showing no `num_traits`, `FromPrimitive`, `T::from_`, `from_usize`, or conversion-fallback hits in `geometry/staggered.rs`.
- [x] Replace `cfd-core` boundary time-function and ghost-cell scalar constants/math dispatch from direct `num_traits::FromPrimitive` fallback conversions and nalgebra method dispatch to Eunomia `FloatElement`, propagating the boundary applicator/specification/manager bounds and adding value-semantic boundary tests; verified with `cargo check -p cfd-core --no-default-features`, `cargo check -p cfd-core`, `cargo nextest run -p cfd-core --no-default-features boundary` (30/30 passed), `cargo nextest run -p cfd-core --no-default-features` (167/167 passed), touched-file `rustfmt --check`, and a source scan showing no `num_traits`, `FromPrimitive`, `T::from_f64`, conversion-fallback, bare `.exp()`, or bare `.sin()` hits in the touched boundary files.
- [x] Replace `cfd-core` temperature-dependent fluid scalar constants and math calls from direct `num_traits::FromPrimitive`/nalgebra method dispatch to Eunomia `FloatElement`, remove stale conversion bounds from polynomial models, and add value-semantic polynomial/Andrade/Sutherland tests; verified with `cargo check -p cfd-core --no-default-features`, `cargo check -p cfd-core`, `cargo nextest run -p cfd-core --no-default-features fluid::temperature` (6/6 passed), `cargo nextest run -p cfd-core --no-default-features` (162/162 passed), touched-file `rustfmt --check`, and a source scan showing no `num_traits`, `FromPrimitive`, `T::from_f64`, conversion-fallback, bare `.exp()`, or bare `.powf(...)` hits in `physics/fluid/temperature.rs`.
- [x] Replace `cfd-core` hemolysis calculator and platelet activation scalar constants from direct `num_traits::FromPrimitive` fallback conversions to Eunomia `FloatElement`, and fix platelet activation to use the decaying exponential probability contract; verified with `cargo check -p cfd-core --no-default-features`, `cargo check -p cfd-core`, `cargo nextest run -p cfd-core --no-default-features hemolysis` (9/9 passed), `cargo nextest run -p cfd-core --no-default-features` (160/160 passed), touched-file `rustfmt --check`, and a source scan showing no `num_traits`, `FromPrimitive`, `T::from_f64`, or conversion-fallback hits in `physics/hemolysis/{calculator,trauma}.rs`.
- [x] Replace `cfd-core` mesh quality threshold constants from direct `num_traits::FromPrimitive` conversions to Eunomia `FloatElement`, and add value-semantic quality-level and recommendation tests; verified with `cargo check -p cfd-core --no-default-features`, `cargo check -p cfd-core`, `cargo nextest run -p cfd-core --no-default-features mesh_quality` (3/3 passed), `cargo nextest run -p cfd-core --no-default-features` (158/158 passed), touched-file `rustfmt --check`, and a source scan showing no `num_traits`, `FromPrimitive`, `T::from_f64`, or conversion-fallback hits in `geometry/mesh/quality.rs`.
- [x] Replace `cfd-core` CPU backend scalar parameter conversion from direct `num_traits::FromPrimitive` and the local `safe_f64_to_t` fallback helper to Eunomia `FloatElement`, while narrowing `CpuBuffer` impl bounds to conversion-free `RealField + Copy`; verified with `cargo check -p cfd-core --no-default-features`, `cargo check -p cfd-core`, `cargo nextest run -p cfd-core --no-default-features test_cpu_advection_kernel_linear_exactness` (1/1 passed), `cargo nextest run -p cfd-core --no-default-features` (155/155 passed), touched-file `rustfmt --check`, and a source scan showing no `num_traits`, `FromPrimitive`, `safe_f64_to_t`, or conversion fallback hits in `compute/cpu.rs`.
- [x] Replace `cfd-core` time integrator scalar constants from direct `num_traits::FromPrimitive` conversions to Eunomia `FloatElement`, remove silent implicit-solver default fallback tolerances, and add value-semantic explicit/implicit integrator tests; verified with `cargo check -p cfd-core --no-default-features`, `cargo check -p cfd-core`, `cargo nextest run -p cfd-core --no-default-features compute::time::integrators` (6/6 passed), `cargo nextest run -p cfd-core --no-default-features` (155/155 passed), touched-file `rustfmt --check`, and a source scan showing no `num_traits`, `FromPrimitive`, `T::from_f64`, or conversion-fallback hits in `compute/time/integrators.rs`.
- [x] Replace `cfd-core` time-step controller default constants and runtime math helpers from direct `num_traits::{FromPrimitive, Float}` usage to Eunomia `FloatElement`, and change invalid integration-order handling from silent fallback to typed `Result` errors; verified with `cargo check -p cfd-core --no-default-features`, `cargo nextest run -p cfd-core --no-default-features compute::time::controllers` (4/4 passed), `cargo nextest run -p cfd-core --no-default-features` (149/149 passed), touched-file `rustfmt --check`, and a source scan showing no `num_traits`, `FromPrimitive`, `T::from_f64`, or conversion-fallback hits in `compute/time/controllers.rs`.
- [x] Replace `cfd-core` solver configuration default scalar constants from direct `num_traits::FromPrimitive` conversions to Eunomia `FloatElement`, consolidate builder default construction through `SolverConfig::default`, and add value-semantic default/builder parity tests; verified with `cargo check -p cfd-core --no-default-features`, `cargo check -p cfd-core`, `cargo nextest run -p cfd-core --no-default-features solver_config` (2/2 passed), `cargo nextest run -p cfd-core --no-default-features` (146/146 passed), touched-file `rustfmt --check`, and a source scan showing no `num_traits`, `FromPrimitive`, `T::from_f64`, or zero-conversion fallback hits in `compute/solver/config.rs`.
- [x] Replace `cfd-core` abstraction default scalar constants from direct `num_traits::FromPrimitive` conversions to Eunomia `FloatElement`, and narrow field-state methods so only constructors/defaults carry scalar-conversion bounds; verified with `cargo check -p cfd-core --no-default-features`, `cargo check -p cfd-core`, `cargo nextest run -p cfd-core --no-default-features abstractions` (2/2 passed), `cargo nextest run -p cfd-core --no-default-features` (144/144 passed), touched-file `rustfmt --check`, and a source scan showing no `num_traits`, `FromPrimitive`, `T::from_f64`, or zero/one conversion-fallback hits in `abstractions/{state,problem}.rs`.
- [x] Replace `cfd-core` fluid property validation thresholds from direct `num_traits::FromPrimitive` conversions to Eunomia `FloatElement`, preserving existing Reynolds/Prandtl/temperature/pressure/property bound validation behavior; verified with `cargo check -p cfd-core --no-default-features`, `cargo check -p cfd-core`, `cargo nextest run -p cfd-core --no-default-features physics::fluid::validation` (2/2 passed), `cargo nextest run -p cfd-core --no-default-features` (144/144 passed), touched-file `rustfmt --check`, touched-file `git diff --check`, and a source scan showing no `num_traits`, `FromPrimitive`, `T::from_f64`, or conversion-fallback hits in `validation.rs`.
- [x] Replace `cfd-core` material and fluid constant constructors from direct `num_traits::FromPrimitive` conversions to Eunomia `FloatElement`, including common solids, water-air interface constants, constant-property fluids, ideal air, and the material database dependency cone; verified with `cargo check -p cfd-core --no-default-features`, `cargo check -p cfd-core`, `cargo nextest run -p cfd-core --no-default-features` (144/144 passed), touched-file `rustfmt --check`, touched-file `git diff --check`, and source scans over the touched cone showing no `num_traits`, `FromPrimitive`, `T::from_f64`, or conversion-fallback hits.
- [x] Replace `cfd-core` physics value objects and dependent management aggregates from direct `num_traits::FromPrimitive` scalar constants to Eunomia `FloatElement`, remove silent conversion fallbacks in default physical parameters, and make `PhysicalParameters::with_reynolds` surface invalid Reynolds values; verified with `cargo check -p cfd-core --no-default-features`, `cargo nextest run -p cfd-core --no-default-features` (144/144 passed), touched-file `rustfmt --check`, and source scans over the touched cone showing no `num_traits`, `FromPrimitive`, `T::from_f64`, or fallback conversion hits. Package-level `cargo fmt --check --package cfd-core` remains blocked by pre-existing unrelated formatting drift in `error.rs`, `ghost_cells.rs`, and `fluid_dynamics/operations.rs`.
- [x] Replace `cfd-python` blood-model PyO3 numeric conversions with Eunomia `NumericElement`, remove the crate's direct `num-traits` dependency, pass f64 shear rates directly into Rust-owned blood models, and remove silent getter fallback defaults; verified with `cargo check -p cfd-python`, `cargo nextest run -p cfd-python --no-tests pass` (0 test binaries), `cargo fmt -p cfd-python --check`, and a cfd-python source/manifest scan showing no `num_traits`, `num-traits`, `FromPrimitive`, or `ToPrimitive` hits. Dependency-chain clippy blockers in cfd-1d, cfd-2d, cfd-3d, cfd-validation, and cfd-python were mechanically fixed; full `cargo clippy -p cfd-python --all-targets -- -D warnings` now passes.
- [x] Replace `cfd-io` checkpoint, binary, and CSV scalar bounds/conversions with Eunomia `RealField`, remove the crate's direct `num-traits` dependency, and replace silent mass-conservation conversion fallbacks with explicit finite/exact mesh-dimension conversion rejection; verified with `cargo check -p cfd-io`, `cargo check -p cfd-io --all-features`, `cargo nextest run -p cfd-io --no-fail-fast` (3/3 passed), `cargo clippy -p cfd-io --all-targets --all-features -- -D warnings`, `cargo fmt -p cfd-io --check`, and cfd-io source/manifest scans showing no `num_traits`, `num-traits`, `FromPrimitive`, or `ToPrimitive` hits.
- [x] Replace `cfd-io` checkpoint and binary dense matrix/vector payloads with Leto arrays, remove the crate's direct and normal-transitive `nalgebra` dependency path by dropping unused `cfd-math` and `cfd-core` error coupling, and preserve checkpoint row-major roundtrip/checksum behavior; verified with `cargo check -p cfd-io`, `cargo nextest run -p cfd-io --no-fail-fast` (3/3 passed), `cargo check -p cfd-io --all-features`, `cargo fmt -p cfd-io --check`, `cargo tree -p cfd-io -e normal -i nalgebra` (no package match), source scans, and `git diff --check`.
- [x] Move `cfd-python` 2D PyO3 NumPy-return helper paths off direct `ndarray`/`nalgebra` ownership by constructing Leto `Array2` values in Rust and copying to NumPy only at the Python boundary; verified with `cargo check -p cfd-python`, `cargo nextest run -p cfd-python --no-tests pass` (0 test binaries), `cargo fmt -p cfd-python --check`, touched-file diff whitespace check, and no `ndarray`/`nalgebra`/`DMatrix` source hits in `crates/cfd-python`.
- [x] Route `cfd-core` and `cfd-math` GPU synchronous wait/probe boundaries through Moirai instead of `pollster`, remove `pollster` from the workspace dependency graph, and record the Hephaestus blocker: CFDrs still uses `wgpu 0.19` while Hephaestus' WGPU provider is on `wgpu 26.0`; verified with `cargo check -p cfd-core --features gpu`, `cargo nextest run -p cfd-core --features gpu --lib` (151/151 passed), `cargo check -p cfd-math --features gpu`, `cargo nextest run -p cfd-math --features gpu --lib` (279/279 passed), no `pollster` hits in manifests/lockfile/Rust sources, and touched-file rustfmt/diff checks.
- [x] Migrate the `cfd-3d` spectral DNS, forcing, diagnostics, Fourier wrapper, and IBM NUFFT dense-array surface from direct `ndarray` ownership to Leto arrays, with a private Apollo boundary adapter for the currently published Apollo ndarray API; verified with `cargo check -p cfd-3d --no-default-features` and `cargo nextest run -p cfd-3d --no-default-features --lib` (206/206 passed).
- [x] Enable Moirai parallel, Mnemosyne memory integration, and Mellinoe branding through the CFDrs Moirai dependency while preserving the opt-in `cfd-core/mnemosyne` process allocator to avoid `cfd-validation` allocator conflicts; verified with `cargo check --workspace`.
- [x] Audit `cfd-schematics` crate for redundancy and consolidate shared components.
- [x] Enhance deep vertical hierarchical file tree in `cfd-schematics` with DIP, SRP, SSOT, SOC.
- [ ] Optimize performance and memory efficiency representation in `cfd-schematics`.
- [x] Ensure no placeholders, stubs, approximations, or simplifications exist in `cfd-schematics` (implemented native math handling for `n_furcation`).
- [x] Sync README, backlog, and checklist artifacts.
- [x] Integrate Tyche Latin-hypercube sampling into `cfd-optim`, including
  deterministic seed differentiation, invalid-count rejection, value-semantic
  consumer tests, provider-residue scans, and package gates.

## Proteus temperature-response increment

- [x] Advance Aequitas and Proteus to their merged temperature-response
  revisions.
- [x] Replace the local polynomial-fluid density equation with the generic
  Proteus `TemperatureLaw`, using a linear density strategy and zero-sized
  constant heat-capacity and conductivity strategies.
- [x] Propagate material-domain failures through `cfd_core::Error::InvalidInput`.
- [x] Document the ownership decision and verify the independent closed-form
  density oracle plus negative-response rejection.
