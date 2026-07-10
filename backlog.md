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
# CFDrs Backlog

## Structural Improvements
- [x] `cfd-core` [arch]: Consolidate GPU add and scalar multiply on the
  Hephaestus elementwise SSOT. Delete the raw WGPU pipelines, local WGSL,
  staging/polling code, and silent CPU fallbacks; retain one fallible
  `GpuFieldOps` domain surface with exact-value and typed-error verification.
  Completed with clean GPU/no-default compilation, all-target clippy, 230/230
  nextest, 5/5 doctests, warning-clean docs, and zero arithmetic provider
  residue.
- [x] `cfd-1d`/`cfd-3d` [patch]: Remove direct `num-traits` scalar
  dependency ownership from the crate scalar seams. The workspace dependency
  catalog and both crate manifests no longer declare `num-traits`, and
  `Cfd1dScalar`/`Cfd3dScalar` source identity values from Eunomia
  `NumericElement` constants. Evidence: touched-file rustfmt, dual-package
  check, targeted direct-dependency/source residue scan, and dual-package
  nextest 1122/1122 with one existing slow 3D mesh-convergence validation.
  Residual: package-wide fmt/clippy are still blocked by unrelated baseline
  formatting and lint debt outside the touched scalar/manifests cone.
- [x] `cfd-1d` [patch]: Move solver-core reusable workspace vector storage
  to Leto arrays. `SolverWorkspace::{rhs,last_solution,linear_solution}` now
  stores `leto::Array1<T>` values, matrix assembly returns a Leto RHS array,
  validation/residual helpers read Leto RHS storage directly, and the private
  `vector_bridge` module owns the remaining nalgebra<->Leto conversions at
  the current linear-solver/Anderson solution boundary. Evidence: cfd-1d
  test-target check, cfd-1d lib clippy, focused solver-core nextest (10/10),
  and cfd-1d doc generation completed with the same unrelated 11 rustdoc link
  warnings. Residual: `linear_system.rs` still exposes nalgebra solution and
  dense fallback storage plus the nalgebra-sparse matrix input boundary;
  `matrix_assembly.rs` still returns nalgebra-sparse CSR; solver detection and
  Anderson acceleration still operate on nalgebra solution vectors.
- [x] `cfd-1d` [patch]: Move the solver-core `NetworkState` public storage
  boundary to Leto arrays. `NetworkState::{pressures,flow_rates}` now exposes
  `leto::Array1<T>` instead of nalgebra `DVector<T>`, and state construction
  populates Leto arrays directly from network pressure/flow slices. Focused
  unit tests assert Leto array shapes, zero/value semantics, clone
  preservation, and time mutation. Evidence: cfd-1d test-target check,
  cfd-1d lib clippy, focused state nextest (2/2), and cfd-1d doc generation
  completed with the same unrelated 11 rustdoc link warnings; targeted state
  residue scan found no `DVector`. The upstream Atlas blocker was closed by
  adding `FloatElement::acos` in Eunomia. Residual: the broader cfd-1d
  linear-system matrix/solution, sparse-assembly, solver-detection, and
  Anderson solution-vector boundaries still use nalgebra storage.
- [x] `cfd-1d` [patch]: Move the solver-core convergence checker public vector
  boundary to Leto arrays. `ConvergenceChecker::{check,has_converged,
  has_converged_dual}` now accepts `leto::Array1<T>` instead of nalgebra
  `DVector<T>`, computes L2 norms from values directly, and returns typed
  `InvalidInput` for mismatched vector lengths. The current network solver
  still owns nalgebra workspaces and converts once at the checker call site.
  Evidence: cfd-1d test-target check, focused cfd-1d clippy for lib plus
  `solver_core_tests`, focused convergence nextest 20/20, and a focused
  residue scan showing no `DVector` in `convergence.rs` or the migrated checker
  tests. Residual: cfd-1d doc generation still reports unrelated rustdoc link
  warnings, and the broader linear-system matrix/solution, sparse-assembly,
  solver-detection, and Anderson solution-vector boundaries still use nalgebra
  storage.
- [x] `cfd-suite` [patch]: Remove direct `num-traits` and `simba`
  dependencies from the root package and workspace dependency catalog. Source
  and member manifests now rely on Eunomia/Atlas scalar contracts without
  direct legacy numeric dependencies from CFDrs. Evidence: targeted
  source/member manifest residue scan and root package metadata. Residual:
  `num-traits` and `simba` remain in `Cargo.lock` only through transitive
  dependencies owned by upstream crates; root package check/nextest were
  blocked by concurrent shared Cargo cache/build locks and did not run.
- [x] `cfd-suite`/`cfd-3d` [patch]: Remove direct `crossbeam` ownership from
  the workspace dependency catalog and `cfd-3d` manifest. Current source and
  member manifests have no direct Crossbeam imports or workspace dependency
  users, leaving CFDrs concurrency ownership on Moirai-facing surfaces.
  Evidence: targeted source/member manifest residue scan, root package
  metadata, `cfd-3d` no-default-features check, `cfd-3d`
  no-default-features nextest (394/394), and diff whitespace check. Residual:
  this does not claim the full transitive graph is Crossbeam-free; upstream
  crates own any remaining transitive Crossbeam dependencies. The focused
  nextest run marked existing mesh-convergence validation slow at 19.7s.
- [x] `cfd-math`/`cfd-validation` [patch]: Close the public sparse and
  linear-solver Atlas boundary residual. `cfd-math::sparse::SparseMatrix<T>`
  now aliases `leto_ops::CsrMatrix<T>`; sparse builders, sparse operations,
  `LinearOperator`, preconditioners, direct solver, and solver-chain callers
  consume Leto CSR and Leto `Array1` vectors in the requested public cone.
  `cfd-validation::numerical` now keeps validation solution storage on
  `leto::Array1<T>` and uses the Leto CSR sparse boundary for test cases and
  SpMV benchmarking. Evidence: cfd-math check, cfd-math test-target check,
  cfd-math all-target clippy, cfd-math doc, cfd-math nextest 361/361,
  cfd-validation check, cfd-validation all-target clippy, cfd-validation doc,
  and a targeted sparse/linear-solver/validation residue scan with no exposed
  nalgebra sparse/vector matches. Residual: package cfd-validation nextest is
  still blocked by the existing venturi cross-fidelity convergence failures in
  `option2_selected_45um_geometry_routes_to_fallback_and_converges` and
  `microventuri_35um_case_produces_converged_informative_2d_result`, outside
  this sparse/linear-solver boundary.
- [x] `cfd-math` [patch]: Move IncompleteCholesky to Leto CSR construction.
  `src/linear_solver/preconditioners/cholesky.rs` now stores
  `leto_ops::CsrMatrix`, performs symmetry and factorization lookups through
  Leto CSR row views, constructs IC(0) factors with `CsrMatrix::from_parts`,
  and applies forward/backward substitutions through Leto row views. Source
  and integration preconditioner edge tests now pass Leto CSR into
  `IncompleteCholesky::new`. Evidence: cfd-math fmt, lib check, all-target
  check, lib/tests clippy, all-target clippy, cholesky-filter nextest (5/5),
  preconditioner-filter nextest (76/76), and targeted Cholesky residue scans.
  Residual: Schwarz, direct solver, remaining transitional solver fixtures,
  and the shared solver sparse matrix boundary still use nalgebra-sparse.
- [x] `cfd-math` [patch]: Move ILU preconditioner family to Leto CSR
  construction. `src/linear_solver/preconditioners/ilu` now factorizes and
  applies ILU(0)/ILU(k) using `leto_ops::CsrMatrix`, with triangular solves
  walking Leto CSR row views. Source tests, integration edge tests,
  `LinearSolverChain`, and Schwarz local-solver construction now pass Leto CSR
  matrices or convert once at the remaining shared solver/Schwarz boundary.
  Evidence: cfd-math fmt, lib check, all-target check, lib/tests clippy,
  all-target clippy, ilu-filter nextest (21/21), preconditioner-filter nextest
  (76/76), linear_solver::tests nextest (53/53), and targeted ILU residue
  scans. Residual: Schwarz, direct solver, remaining transitional solver
  fixtures, and the shared solver sparse matrix boundary still use
  nalgebra-sparse.
- [x] `cfd-math` [patch]: Move SSOR preconditioner to Leto CSR construction.
  `src/linear_solver/preconditioners/ssor.rs` now stores
  `leto_ops::CsrMatrix` and performs forward/backward sweeps through Leto CSR
  row views. Source preconditioner edge tests convert transitional solver CSR
  fixtures once before constructing SSOR and assert exact typed errors for
  SSOR vector-length mismatches. Evidence: cfd-math fmt, lib check,
  all-target check, lib/tests clippy, all-target clippy, ssor-filter nextest
  (5/5), preconditioner-filter nextest (76/76), and targeted SSOR residue
  scans. Residual: Schwarz, direct solver, integration-test
  fixtures, and the shared solver sparse matrix boundary still use
  nalgebra-sparse.
- [x] `cfd-math` [patch]: Move basic preconditioners to Leto CSR construction.
  `src/linear_solver/preconditioners/basic.rs` now constructs Jacobi and SOR
  preconditioners from `leto_ops::CsrMatrix`, using Leto CSR diagonal/row
  access and Leto `Array1` application with explicit length validation. The
  source linear-solver tests and `tests/core_solver_tests.rs` convert their
  still-transitional solver CSR fixtures once before constructing Jacobi/SOR.
  Evidence: cfd-math fmt, lib check, core solver test check, all-target check,
  lib/tests clippy, core solver test clippy, all-target clippy,
  linear_solver::tests nextest (53/53), core_solver_tests nextest (4/4),
  preconditioner nextest (76/76), and targeted basic-preconditioner residue
  scans. Residual: Schwarz, direct solver, and the shared solver sparse matrix
  boundary still use nalgebra-sparse.
- [x] `cfd-math` [patch]: Move AMG/coarsening sparse boundary to Leto CSR.
  The multigrid bounded context now uses a local `leto_ops::CsrMatrix` alias
  and Leto CSR helpers for AMG setup, coarsening, interpolation, smoothers,
  cycles, sparse products, transpose, and SpMV. `coarsening_bench.rs`,
  `algebraic_distance_bench.rs`, `amg_coarsening_tests.rs`, and AMG
  integration preconditioner construction now use Leto CSR rather than direct
  nalgebra-sparse construction. Evidence: cfd-math fmt, lib check, focused
  AMG/coarsening test and bench checks, focused clippy for lib/touched
  tests/benches, AMG integration nextest (5/5), AMG-filter nextest (6/6),
  multigrid::coarsening nextest (10/10), and a targeted AMG/coarsening sparse
  residue scan with no nalgebra sparse/Coo/direct CSR/get_entry/old SpMV
  residue in the migrated boundary. Residual: the broader
  `crate::sparse::SparseMatrix` solver/direct/preconditioner surface remains
  nalgebra-sparse; `LinearSolverChain` converts that boundary once before AMG.
- [x] `cfd-math` [patch]: Migrate SpMV/CG benchmarks to direct Leto CSR.
  `benches/spmv_bench.rs`, `benches/cg_bench.rs`, and the CG section of
  `benches/math_benchmarks.rs` now construct `leto_ops::CsrMatrix` directly
  with `from_parts`. The SpMV benchmark now dispatches through the direct Leto
  CSR `LinearOperator::apply` path instead of the legacy nalgebra sparse
  helper. Evidence: cfd-math fmt, focused bench checks, focused bench clippy
  with `-D warnings`, all-target check, all-target clippy, sparse-filter
  nextest (19/19), and a targeted migrated-benchmark residue scan with no
  nalgebra sparse/vector or old SpMV helper matches. Residual:
  after the AMG/coarsening boundary moved to Leto CSR, the remaining sparse
  provider gap is the broader solver/direct/preconditioner nalgebra-sparse
  surface rather than those two benchmarks.
- [x] `cfd-math` [patch]: Add a direct Leto CSR linear-operator path.
  `src/sparse/operations.rs` now implements `LinearOperator<T>` for
  `leto_ops::CsrMatrix<T>` and shares one `try_leto_spmv` execution helper
  between Atlas-native CSR inputs and the remaining legacy
  `nalgebra_sparse::CsrMatrix` conversion path. `tests/simple_gmres_tests.rs`
  now builds Leto CSR matrices directly and verifies GMRES residuals through
  that operator path instead of constructing nalgebra COO/CSR matrices or
  using nalgebra vector residuals. Evidence: cfd-math fmt, simple_gmres test
  check, simple_gmres nextest (3/3), simple_gmres clippy, cfd-math lib check,
  all-target check, all-target clippy, sparse-filter nextest (19/19),
  gmres-filter nextest (21/21), and a targeted simple GMRES residue scan with
  no nalgebra sparse/vector or old sparse helper residue. Residual:
  cfd-math still exposes the broader nalgebra-sparse public matrix boundary in
  sparse builders, preconditioners, AMG, direct solver, tests, and benches.
- [x] `cfd-math` [patch]: Close remaining Leto storage-slice residue. Nonlinear
  dense pivoting now swaps `Array1`/`Array2` entries through direct indexing,
  AMG interpolation quality validation indexes interpolated Leto arrays
  directly, and multigrid smoother tests assert expected vectors by indexed
  value comparison instead of `storage().as_slice()`. Evidence: cfd-math fmt,
  lib check, focused nonlinear_solver/multigrid nextest (46/46), lib/tests
  clippy, all-target check, all-target clippy, and a cfd-math `src`/`tests`
  residue scan with no `leto::Storage`, `StorageMut`,
  `.storage().as_slice()`, `as_slice_mut()`, `vector_slice_mut`, or
  `matrix_slice_mut` matches. Residual cfd-math provider work is now
  nalgebra/nalgebra-sparse and other Atlas boundaries, not Leto
  storage-slice access.
- [x] `cfd-math` [patch]: Remove sparse/basic preconditioner Leto
  storage-slice dependencies. `src/sparse/operations.rs` now stages SpMV
  input/output and row/column scaling through direct `Array1` indexing before
  delegating to Leto CSR operations, with no `leto::Storage`,
  `.storage().as_slice()`, or output `as_slice_mut()` assumptions. The Jacobi
  preconditioner in `src/linear_solver/preconditioners/basic.rs` now indexes
  the Leto diagonal directly. Evidence: cfd-math fmt, lib check, focused
  sparse/preconditioner nextest (95/95), lib/tests clippy, all-target check,
  all-target clippy, and targeted residue scans over the migrated files.
  Residual storage-slice owners are nonlinear mutable dense-workspace helpers
  and multigrid interpolation/smoother internals.
- [x] `cfd-math` [patch]: Remove GPU operator Leto storage-slice dependency.
  `src/linear_solver/operators/gpu.rs` now validates Leto `Array1`
  input/output lengths, stages GPU upload/readback buffers through direct
  `Array1` indexing, and writes results back by index instead of importing
  `leto::Storage`, borrowing `.storage().as_slice()`, or requiring output
  `as_slice_mut()` contiguity. The operator still routes execution through the
  existing Hephaestus-backed `cfd-core` GPU context/buffer/kernel path.
  Evidence: cfd-math fmt, GPU-feature check, focused GPU-feature
  linear_solver::operators nextest (5/5), GPU-feature lib clippy,
  GPU-feature all-target check, GPU-feature all-target clippy, and a targeted
  scan with no storage-slice residue in `operators/gpu.rs`. Residual
  storage-slice owners are sparse operations and multigrid internals; broader
  GPU provider work remains in raw WGPU context/kernel ownership outside this
  operator.
- [x] `cfd-math` [patch]: Remove finite-difference operator Leto storage-slice
  dependencies. `src/linear_solver/operators/{poisson,momentum}.rs` now read
  and write `Array1` values through direct indexing for 2D Laplacian, 3D
  Poisson, 1D/2D momentum, and 2D energy operators, removing `leto::Storage`,
  `.storage().as_slice()`, and output `as_slice_mut()` assumptions from those
  CPU operator paths. Evidence: cfd-math fmt, lib check, focused
  linear_solver::operators nextest (5/5), lib/tests clippy, all-target check,
  all-target clippy, and a targeted scan with no storage-slice residue in the
  migrated operator files. Residual storage-slice owners are sparse
  operations, GPU operator, and multigrid internals.
- [x] `cfd-math` [patch]: Remove nonlinear linalg immutable Leto storage-slice
  dependency. `src/nonlinear_solver/linalg.rs` now evaluates vector dot,
  add/sub, scaled add, in-place scaled add, and scale through direct
  `Array1` indexing; `anderson.rs` now indexes `gamma` directly, allowing the
  immutable `vector_slice` helper and `leto::Storage` import to be removed.
  Evidence: cfd-math fmt, lib check, nonlinear_solver nextest (9/9),
  lib/tests clippy, all-target check, all-target clippy, and a targeted scan
  with no immutable storage-slice/vector_slice residue in nonlinear
  linalg/Anderson. Residual mutable dense-workspace helpers still use
  `StorageMut`; remaining immutable storage-slice owners are sparse
  operations, linear operators, GPU operator, and multigrid internals.
- [x] `cfd-math` [patch]: Remove production SIMD vector Leto storage-slice
  dependency. `src/simd/vector.rs` now evaluates `simd_mul`, `simd_dot`,
  `simd_norm`, and `par_map` through direct `Array1` indexing while preserving
  Moirai `Adaptive` map/reduce dispatch. Evidence: cfd-math fmt, lib check,
  focused simd::vector nextest (1/1), lib/tests clippy, all-target check,
  all-target clippy, SIMD-filter nextest (26/26), and a targeted scan with no
  `Storage` or `.storage().as_slice()` residue in `src/simd/vector.rs`.
  Residual source-level storage-slice owners are now nonlinear linalg, sparse
  operations, linear operators, GPU operator, and multigrid internals.
- [x] `cfd-math` [patch]: Remove SIMD integration test Leto storage-slice
  bridge. `tests/simd_tests.rs` now indexes the Leto `spmv` result directly
  and passes collected values to the existing SIMD slice API without importing
  `leto::Storage` or calling `.storage().as_slice()`. Evidence: cfd-math fmt,
  simd_tests check, simd_tests nextest (12/12), simd_tests clippy, cfd-math
  all-target check, cfd-math all-target clippy, SIMD-filter nextest (26/26),
  and targeted scans with no integration-test `DVector`/nalgebra vector/local
  preconditioner bridge/`Storage`/storage-slice residue under
  `crates/cfd-math/tests`. Residual: remaining provider work has moved to
  broader production/test boundaries, including `nalgebra_sparse::CsrMatrix`,
  dense nalgebra test oracles, and source-level Leto storage-slice internals.
- [x] `cfd-math` [patch]: Move AMG integration vector paths to Leto arrays.
  AMG exact solutions, RHS values, solver outputs, cycle outputs, and two-grid
  preconditioner work buffers now allocate `leto::Array1` directly and removed
  nalgebra `DVector` plus Leto storage-slice conversions from
  `tests/amg_integration_test.rs`. Evidence: cfd-math fmt,
  amg_integration_test check, amg_integration_test nextest (5/5),
  amg_integration_test clippy, cfd-math all-target check, cfd-math all-target
  clippy, AMG-filter nextest (6/6), and a targeted scan with no
  `DVector`/nalgebra vector/`Storage`/storage-slice/preconditioner bridge
  residue in the migrated test. Residual: this test still uses
  `nalgebra_sparse::CsrMatrix` for sparse storage and nalgebra
  `DMatrix`/`SymmetricEigen` for the dense energy-norm oracle; the integration
  storage-slice residue was closed by Sprint 1.96.153.
- [x] `cfd-math` [patch]: Move preconditioner edge-case integration tests to
  Leto arrays. ILU(0), ILU(k), repeated-application, extreme-value, and
  sparsity-preservation tests now allocate `leto::Array1` RHS/solution buffers
  directly and removed the local nalgebra `DVector` preconditioner bridge.
  Evidence: cfd-math fmt, preconditioner_edge_cases check,
  preconditioner_edge_cases nextest (6/6), preconditioner_edge_cases clippy,
  cfd-math all-target check, cfd-math all-target clippy, broader
  preconditioner nextest (76/76), and a targeted residue scan with no
  `DVector`/nalgebra vector/`Storage`/storage-slice bridge residue in the
  migrated test. Residual: this test still uses the shared
  `nalgebra_sparse::CsrMatrix` matrix boundary; integration-test vector bridge
  residue was closed by Sprints 1.96.152 and 1.96.153.
- [x] `cfd-math` [patch]: Move linear-solver source test module tree to Leto
  arrays. The `src/linear_solver/tests` module tree now allocates
  `leto::Array1` RHS/solution/work buffers directly and removed local
  nalgebra `DVector` solve/preconditioner bridge macros. Solver residual
  checks now use the Leto SpMV/helper path, and the touched solver/sparse cone
  routes scalar constants/tolerances through Eunomia instead of old scalar
  helpers. Evidence: cfd-math fmt, lib check, all-target check,
  linear_solver::tests nextest (53/53), linear_solver nextest (176/176),
  lib/tests clippy, all-target clippy, and current targeted residue scans with
  no source-test `DVector`/`Storage` bridge residue and no old scalar helper
  residue in the searched solver/sparse cone. Residual: integration-test
  vector bridge holdouts were closed by Sprints 1.96.152 and 1.96.153, and
  production sparse/scalar boundaries still include `nalgebra_sparse` plus
  transitional nalgebra utilities.
- [x] `cfd-math` [patch]: Move core solver validation tests to Leto arrays.
  BiCGSTAB, GMRES, preconditioner integration, and condition-number robustness
  tests now allocate `leto::Array1` RHS/solution buffers directly and removed
  the local nalgebra `DVector` solve/preconditioner bridge. Residual checks now
  run through the Leto SpMV path and assert value-semantic thresholds.
  Evidence: cfd-math fmt, core_solver_tests check, core_solver_tests nextest
  (4/4), core_solver_tests clippy, cfd-math all-target check, cfd-math
  all-target clippy, and a targeted core solver test scan with no
  `DVector`/nalgebra vector/`Storage`/storage-slice conversion residue.
  Residual: this test still uses the shared `nalgebra_sparse::CsrMatrix`/
  `CooMatrix` matrix boundary, and broader cfd-math test diagnostics still
  contain nalgebra `DVector` bridges.
- [x] `cfd-math` [patch]: Move simple GMRES integration tests to Leto arrays.
  Basic, restarted, and preconditioned GMRES integration tests now allocate
  `leto::Array1` RHS/solution buffers directly and removed the local
  nalgebra `DVector` solve bridge macro. Residual checks now run through the
  Leto SpMV path and assert value-semantic thresholds. Evidence: cfd-math fmt,
  simple_gmres test check, simple_gmres nextest (3/3), simple_gmres clippy,
  cfd-math all-target check, cfd-math all-target clippy, touched-file
  `git diff --check`, and a targeted simple GMRES test scan with no
  `DVector`/nalgebra vector/`Storage`/storage-slice conversion residue.
  Residual: this test still uses the shared `nalgebra_sparse::CsrMatrix`/
  `CooMatrix` matrix boundary, and broader cfd-math test diagnostics still
  contain nalgebra `DVector` bridges.
- [x] `cfd-math` [patch]: Move matrix-free solver tests to Leto arrays.
  Matrix-free CG/GMRES identity and scaled-operator tests now allocate
  `leto::Array1` RHS/solution buffers directly and no longer use a local
  nalgebra `DVector` bridge macro. The operator-size mismatch test asserts the
  exact typed Leto-boundary error. Evidence: cfd-math fmt, all-target check,
  all-target clippy, focused matrix-free nextest (4/4), touched-file
  `git diff --check`, and a targeted matrix-free test residue scan with no
  `DVector`/nalgebra/`Storage`/storage-slice conversion residue. Residual:
  broader cfd-math linear-solver integration/adversarial/core/preconditioner
  test diagnostics still contain nalgebra `DVector` bridges, and production
  sparse/scalar provider boundaries still include `nalgebra_sparse::CsrMatrix`
  and `nalgebra::RealField`.
- [x] `cfd-math` [patch]: Move BiCGSTAB workspaces and final chain fallback to
  Leto arrays. BiCGSTAB direct solve methods, `IterativeLinearSolver`
  dispatch, residual/search/operator/preconditioned workspaces, and the
  `LinearSolverChain` last-resort fallback now operate on `leto::Array1`
  without nalgebra `DVector` conversion. CG and BiCGSTAB share one Leto vector
  helper module, and obsolete legacy bridge helpers were removed from
  `linear_solver/traits.rs`. Evidence: cfd-math fmt, all-target check,
  all-target clippy, focused BiCGSTAB nextest (24/24), broader linear-solver
  nextest (176/176), AMG integration nextest (5/5), touched-file
  `git diff --check`, and targeted scans with no legacy bridge helper residue
  or migrated BiCGSTAB/chain/traits `DVector` ownership. Residual: cfd-math
  still carries `nalgebra_sparse::CsrMatrix`, transitional
  `nalgebra::RealField` scalar bounds, and nalgebra `DVector` in remaining
  matrix-free/preconditioner/integration test diagnostics.
- [x] `cfd-math` [patch]: Move Conjugate Gradient workspaces to Leto arrays.
  CG reusable work buffers, preconditioned/unpreconditioned solve methods, and
  `IterativeLinearSolver` dispatch now operate directly on `leto::Array1`
  without nalgebra `DVector` or legacy bridge helpers. CG benchmark call sites
  now measure the Leto API. Tests assert solved-system values plus exact
  dimension and max-iteration errors. Evidence: cfd-math fmt, all-target
  check, all-target clippy, focused conjugate nextest (13/13), broader
  linear-solver nextest (176/176), and a targeted CG/benchmark residue scan
  with no `DVector`/legacy-helper/vector-conversion residue. Residual after
  the BiCGSTAB follow-up: the linear-solver trait family still carries the
  transitional `nalgebra::RealField` scalar bound and sparse storage remains
  on `nalgebra_sparse::CsrMatrix`.
- [x] `cfd-math` [patch]: Move Schwarz preconditioner local apply paths to Leto
  arrays. Additive and multiplicative Schwarz now accept and return
  `leto::Array1`; local RHS extraction uses Leto buffers directly; local ILU
  solves no longer require nalgebra `DVector` conversion bridges; and
  `Preconditioner::apply_to` copies from a Leto result. Residual/output length
  mismatches return typed configuration errors with value-semantic coverage.
  Evidence: cfd-math fmt, all-target check, all-target clippy, focused
  Schwarz nextest (3/3), broader preconditioner nextest (76/76), and a
  targeted `schwarz.rs` residue scan with no
  `DVector`/`Storage`/`num_traits`/local-RHS conversion residue. Residual:
  Schwarz still uses the shared `nalgebra_sparse::CsrMatrix` sparse boundary
  and the global `Preconditioner<T>` nalgebra scalar bound.
- [x] `cfd-math` [patch]: Move ILU triangular solve workspaces to Leto arrays.
  ILU forward/backward substitution and `IncompleteLU::apply_to` now use
  `leto::Array1` residual/intermediate/solution buffers without nalgebra
  `DVector` workspaces; residual/output length mismatches return typed
  configuration errors; and U-solve diagonal identity dispatches through
  Eunomia. Evidence: cfd-math fmt, all-target check, all-target clippy,
  focused ILU nextest (21/21), broader preconditioner nextest (74/74), and a
  targeted `ilu/{types,triangular}.rs` residue scan with no
  `DVector`/`DMatrix`/old-scalar-identity/vector-conversion residue.
  Residual: IncompleteLU still stores the shared
  `nalgebra_sparse::CsrMatrix` LU factor boundary and the global
  preconditioner trait still carries the nalgebra `RealField` bound.
- [x] `cfd-math` [patch]: Move deflation preconditioner vector state to Leto
  arrays. Deflation eigenvectors are now stored as `leto::Array1`; added
  eigenpairs use the Leto boundary; projection coefficients are computed
  without nalgebra `DVector` workspaces; residual/output/eigenvector length
  mismatches return typed configuration errors; and zero eigenvalues are
  rejected before division. Evidence: cfd-math fmt, all-target check,
  all-target clippy, focused deflation nextest (3/3), broader preconditioner
  nextest (73/73), and a targeted `deflation.rs` residue scan with no
  `DVector`/`DMatrix`/old-scalar-identity/vector-conversion residue.
  Residual: Deflation still wraps the base preconditioner behind the existing
  `Box<dyn Preconditioner<T>>`, and the global preconditioner trait still
  carries the nalgebra `RealField` bound.
- [x] `cfd-math` [patch]: Move basic preconditioner vector internals to Leto
  arrays. Identity/Jacobi/SOR now apply at the `leto::Array1` residual/output
  boundary; Jacobi stores its inverse diagonal as `leto::Array1`; Jacobi/SOR
  matrix-sized residuals and all outputs are validated with typed
  configuration errors; and scalar constants/tolerances route through
  Eunomia. Evidence: cfd-math fmt, all-target check, all-target clippy,
  focused mismatch-regression nextest (1/1), broader preconditioner nextest
  (70/70), and a targeted source scan with no
  `DVector`/`DMatrix`/old-scalar-identity/basic vector residue in
  `basic.rs`. Residual: basic matrix-backed preconditioners still use the
  shared `nalgebra_sparse::CsrMatrix` boundary and the global
  `Preconditioner<T>` nalgebra scalar bound.
- [x] `cfd-math` [patch]: Move IncompleteCholesky preconditioner solve
  workspaces to Leto arrays. Forward/backward triangular substitution and
  `apply_to` now operate directly on `leto::Array1` buffers and no longer
  allocate nalgebra `DVector` residual/intermediate/solution workspaces;
  mismatched residual/output lengths return typed configuration errors. The
  IC(0) factorization square root now dispatches through Eunomia
  `NumericElement`. Evidence: cfd-math fmt, all-target check, all-target
  clippy, focused Cholesky nextest (5/5), and a targeted source scan with no
  `DVector`/`DMatrix`/workspace-conversion residue in `cholesky.rs`.
  Follow-up: Sprint 1.96.166 moved IncompleteCholesky factor storage to Leto
  CSR. Residual sparse-provider work is now Schwarz, direct solver, and the
  shared solver matrix boundary.
- [x] `cfd-math` [patch]: Move SSOR preconditioner sweep workspaces to Leto
  arrays. `forward_sweep`, `backward_sweep`, and `apply_to` now operate
  directly on `leto::Array1` buffers and no longer allocate nalgebra
  `DVector` residual/solution workspaces; mismatched residual/output lengths
  return typed configuration errors. Evidence: cfd-math fmt, all-target
  check, all-target clippy, focused SSOR nextest (5/5), and a targeted source
  scan with no `DVector`/`DMatrix`/fallback-vector clone residue in `ssor.rs`.
  Residual: SSOR still uses the shared `nalgebra_sparse::CsrMatrix` storage
  boundary and the global `Preconditioner<T>` nalgebra scalar bound.
- [x] `cfd-math` [patch]: Move block/SIMPLE preconditioner vector internals to
  Leto arrays. Diagonal inverse storage, Schur diagonal storage, direct
  `apply` methods, and `Preconditioner::apply_to` implementations in
  `linear_solver/block_preconditioner.rs` now use `leto::Array1` without
  constructing nalgebra `DVector` bridges; dimension mismatches are explicit
  typed errors. Evidence: cfd-math fmt, all-target check, all-target clippy,
  focused block-preconditioner nextest (4/4), and a targeted source scan with
  no `DVector`/`DMatrix`/nalgebra import or fallback-vector clone residue.
  Residual: the global `Preconditioner` trait still retains a
  `nalgebra::RealField` scalar bound and sparse storage still aliases
  `nalgebra_sparse::CsrMatrix`.
- [x] `cfd-math` [patch]: Move GMRES workspace and solver-chain GMRES tiers to
  Leto arrays. GMRES Arnoldi basis extraction, matrix-vector work,
  preconditioner work, residual checks, and public GMRES
  `solve_preconditioned`/`solve_unpreconditioned` calls now use
  `leto::Array1` directly; `LinearSolverChain` keeps GMRES+AMG,
  GMRES+block, unpreconditioned GMRES, and GMRES+ILU tiers on Leto arrays.
  Evidence: cfd-math lib/tests/all-targets check, cfd-math fmt, cfd-math
  all-target clippy, focused GMRES nextest (21/21), cfd-2d no-default lib
  check, focused cfd-2d momentum nextest (53/53), and targeted GMRES
  DVector/legacy-bridge residue scan. Residual after the CG/BiCGSTAB
  follow-up: cfd-math still retains nalgebra scalar bounds and
  nalgebra-sparse storage.
- [x] `cfd-2d` [patch]: Move momentum vectors and the scalar seam off direct
  nalgebra. `physics::momentum` now uses `leto::Array1` RHS/solution buffers,
  boundary helpers mutate Leto RHS arrays, momentum solves call Leto-native
  iterative/direct solver APIs, the obsolete local nalgebra bridge module was
  removed, and `cfd-2d` no longer directly depends on `nalgebra` or
  `nalgebra-sparse`. Evidence: cfd-2d fmt, no-default lib check, no-default
  all-target clippy, focused `momentum` nextest (53/53), direct
  cfd-2d source/manifest nalgebra residue scan, and cargo-tree proof that
  remaining nalgebra is transitive through upstream owners. Residual:
  nalgebra/nalgebra-sparse still resolve transitively through cfd-1d, cfd-core,
  cfd-math, cfd-schematics, and Gaia.
- [x] `cfd-2d` [patch]: Move pressure-velocity correction vectors to Leto
  arrays. `pressure_velocity` now caches RHS and correction solution vectors as
  `leto::Array1`, dispatches directly through the Leto-native iterative solver
  trait and direct sparse fallback, and uses `leto_ops::norm_l2` for RHS
  diagnostics. Evidence: cfd-2d fmt, no-default lib check, no-default
  all-target clippy, focused `pressure_velocity` nextest (16/16), and targeted
  `pressure_velocity` DVector/nalgebra/bridge residue scan. Residual: direct
  cfd-2d source/manifest nalgebra ownership is now removed; nalgebra remains
  transitive through upstream owners.
- [x] `cfd-2d` [patch]: Move SIMPLE pressure-correction vectors to Leto
  arrays. `solvers::simple` now stores RHS and `p_prime` as `leto::Array1`
  and calls the Leto-native `IterativeLinearSolver::solve` boundary directly
  instead of the nalgebra conversion bridge. Evidence: cfd-2d fmt,
  no-default lib check, no-default all-target clippy, focused `simple`
  nextest (19/19), and targeted `solvers/simple` DVector/nalgebra/bridge
  residue scan. Residual: direct cfd-2d source/manifest nalgebra ownership is
  now removed; nalgebra remains transitive through upstream owners.
- [x] `cfd-2d` [patch]: Move FDM RHS and Gauss-Seidel vectors to Leto arrays.
  `solvers::fdm` Poisson, advection-diffusion, and shared Gauss-Seidel solver
  paths now use `leto::Array1` RHS/result vectors instead of nalgebra
  `DVector`. Evidence: cfd-2d fmt, no-default lib check, no-default
  all-target clippy, focused `fdm` nextest (2/2), and targeted `solvers/fdm`
  DVector/nalgebra residue scan. Residual: direct cfd-2d source/manifest
  nalgebra ownership is now removed; nalgebra remains transitive through
  upstream owners.
- [x] `cfd-2d` [patch]: Move time-integration state vectors to Leto arrays.
  `schemes::time` explicit, implicit, multistep, adaptive-controller,
  adaptive-integrator, and test paths now use `leto::Array1` through the
  `StateVector<T>` domain alias instead of nalgebra `DVector`. Fixed-point
  convergence now uses `leto_ops::norm_l2`. Evidence: cfd-2d fmt, no-default
  lib check, no-default all-target clippy, focused `time` nextest (29/29), and
  targeted `schemes/time` DVector/nalgebra residue scans. Residual: direct
  cfd-2d source/manifest nalgebra ownership is now removed; nalgebra remains
  transitive through upstream owners.
- [x] `cfd-2d` [patch]: Move compact `DMatrix` residue to Leto arrays.
  Immersed-boundary force/velocity matrices, `schemes::Grid2D` storage,
  dependent scheme callers/tests, and the `blood_venturi` example now use
  `leto::Array2` with provider-native shape/indexing. Evidence: cfd-2d fmt,
  no-default lib check/clippy, no-default `blood_venturi` example check,
  focused immersed-boundary/schemes/upwind/MUSCL nextest (60/60), targeted
  DMatrix residue scan, and targeted nalgebra-style `Grid2D.data` access scan.
  Residual: direct cfd-2d source/manifest nalgebra ownership is now removed;
  nalgebra remains transitive through upstream owners.
- [x] `cfd-1d` [patch]: Move vascular Bessel/Womersley complex math to
  Eunomia. The Bessel recurrence and Womersley profile now use
  `eunomia::Complex`; the convergence test uses Eunomia's complex `norm()`;
  and the vascular path no longer imports nalgebra `Complex`/`ComplexField`.
  Evidence: cfd-1d lib check/clippy, focused Bessel/Womersley nextest (26/26),
  fmt check, and targeted nalgebra-complex residue scan. Residual: cfd-1d
  still retains nalgebra for network linear-system storage and its transitional
  scalar seam.
- [x] `cfd-math` [patch]: Move the public `LinearOperator::apply` vector
  boundary to Leto. `LinearOperator::apply` and `apply_transpose` now accept
  `leto::Array1` buffers; sparse CSR, identity/scaled, Poisson, momentum,
  energy, and GPU operator adapters implement that boundary; and iterative
  solvers bridge only at their current nalgebra workspaces. Evidence:
  cfd-math all-target check/clippy, focused cfd-math solver/operator nextest
  (80/80), fmt check, and targeted DVector operator-signature residue scan.
  Residual provider work remains in solver workspaces, nalgebra sparse storage,
  and preconditioner internals that still carry local `DVector`/`CsrMatrix`
  bridges.
- [x] `cfd-math`/`cfd-1d` [patch]: Move the public
  `Preconditioner::apply_to` residual/result boundary to Leto. The trait now
  accepts `leto::Array1` buffers; concrete cfd-math preconditioners and their
  tests call that boundary; AMG and Schwarz keep their current matrix
  internals behind Leto vector inputs; and cfd-1d network `DiagJacobi`
  implements the Leto preconditioner contract. Evidence: cfd-math
  all-target check/clippy, cfd-1d lib check/clippy, focused cfd-math
  solver/preconditioner nextest (131/131), fmt check, and targeted
  `apply_to` DVector-signature residue scan. Residual provider work remains
  in `LinearOperator::apply`, nalgebra sparse preconditioner internals, and
  the iterative solver workspaces that still bridge through `DVector`/`CsrMatrix`.
- [x] `cfd-math`/`cfd-1d`/`cfd-2d`/`cfd-3d` [patch]: Move the public
  `IterativeLinearSolver::solve` RHS/result boundary to Leto. The trait now
  accepts Leto arrays, solver implementations bridge internally to the current
  nalgebra workspaces, cfd-math tests construct Leto boundary arrays, cfd-1d
  network solves convert at their solver boundary, cfd-2d momentum/pressure
  solves share `linear_solver_bridge::iterative_solve`, and cfd-3d FEM
  projection uses its FEM Leto bridge. Evidence: cfd-math fmt/check/clippy,
  focused cfd-math solver nextest (61/61), cfd-1d/cfd-2d/cfd-3d focused
  checks and clippy, cfd-validation all-target clippy, targeted DVector
  call-site residue scan, and `git diff --check`. Residual provider work
  remains in `LinearOperator::apply`, `Preconditioner::apply_to`, concrete
  preconditioner internals, and nalgebra sparse storage.
- [x] `cfd-math`/`cfd-validation` [patch]: Move the public iterative solver
  `solve_system` vector boundary to Leto. `LinearSolver::solve_system` and the
  CG, BiCGSTAB, and GMRES implementations now accept Leto RHS/initial-guess
  arrays and return Leto result arrays. Current nalgebra work vectors are
  confined behind solver-local conversion helpers until `LinearOperator`,
  `IterativeLinearSolver`, and `Preconditioner` migrate. cfd-validation
  numerical solver validation now calls the Leto public API and bridges back
  only for its existing DVector-based error metrics. Evidence: cfd-math and
  cfd-validation fmt/check, cfd-math all-target clippy, cfd-validation lib and
  all-target clippy, cfd-2d no-default all-target clippy, focused cfd-math
  CG/BiCGSTAB/GMRES nextest (58/58), targeted `solve_system` signature and
  DVector-residue scans, and `git diff --check`. Residual provider work
  remains in the nalgebra-based iterative operator/preconditioner trait family
  and cfd-validation result/error storage.
- [x] `cfd-validation` [patch]: Move SpMV benchmark callers and validation
  scalar bounds to Leto/Eunomia. Matrix-operation profiling, algorithm
  profiling, and generic SpMV benchmarking now use `leto::Array1` for the
  public `cfd_math::sparse::spmv` API; `LinearSolverValidator` and 1D
  blood-flow literature validation now use the crate-local `ValidationScalar`
  seam needed by migrated cfd-math/cfd-1d provider contracts. Evidence:
  cfd-validation fmt/check, cfd-validation lib and all-target clippy, the
  cfd-2d no-default all-target clippy gate that previously failed in
  cfd-validation, targeted SpMV DVector-residue scan, focused benchmark
  nextest (40/40), and `git diff --check`. Residual: broad cfd-validation
  nextest still fails in two venturi cross-fidelity convergence tests
  unrelated to this provider-boundary slice; cfd-math solver/preconditioner
  traits still expose nalgebra vector/storage boundaries.
- [x] `cfd-math`/`cfd-2d`/`cfd-3d` [patch]: Move the solver-chain and
  FEM/direct-fallback vector boundaries to Leto. `LinearSolverChain::solve`
  and `solve_with_guess` now expose `leto::Array1<T>`; cfd-3d FEM assembly
  uses one private Leto bridge for `DVector` work vectors crossing
  `SparseMatrixBuilder::build_with_rhs`; cfd-2d momentum and pressure direct
  sparse fallbacks use one private bridge into `DirectSparseSolver`; and the
  cfd-2d/cfd-3d scalar seams carry the Leto real-scalar provider contract.
  Evidence: focused fmt/checks for cfd-math/cfd-1d/cfd-2d/cfd-3d, focused
  cfd-math solver-chain nextest (4/4), cfd-math all-target clippy, and
  cfd-2d/cfd-3d no-default lib clippy. Residual provider work remains in the
  nalgebra-based iterative solver traits and cfd-validation, which still
  needs Leto SpMV vectors and propagated Leto scalar bounds before cfd-2d
  all-target clippy can pass.
- [x] `cfd-math` [patch]: Move direct sparse-solver vector API to Leto.
  `DirectSparseSolver::solve` now accepts and returns `leto::Array1<T>`;
  rsparse RHS conversion, sparse solution construction, finite checks, and
  dense fallback all use Leto vectors. The stale DVector dense-fallback wrapper
  was removed, and `LinearSolverChain` performs the only conversion at its
  current nalgebra-facing tier boundary. Evidence: cfd-math fmt/check,
  focused direct-solver/chain/core-solver/simple-GMRES nextest (4/4),
  no-default all-target clippy, cfd-math docs/doctests, targeted
  direct-solver DVector-signature scan, and `git diff --check`. Residual
  cfd-math provider work remains in `LinearSolverChain`, iterative solver
  traits, preconditioners, and sparse storage, which still expose nalgebra
  `DVector`/`CsrMatrix` boundaries.
- [x] `cfd-math` [patch]: Move sparse builder Dirichlet RHS assembly to Leto.
  `SparseMatrixBuilder::build_with_rhs` now accepts `leto::Array1<T>` and
  mutates that provider-owned RHS during column elimination. `build()` no
  longer carries a dummy nalgebra vector, sparse tests assert the exact
  eliminated RHS and matrix rows for a Dirichlet fixture, and direct/block
  preconditioner tests use Leto RHS values for assembly while retaining
  nalgebra only at the still-unmigrated solver boundary. Evidence: cfd-math
  fmt/check, focused sparse/direct-solver/block-preconditioner nextest
  (25/25), no-default all-target clippy, targeted `build_with_rhs` residue
  scan, and `git diff --check`. Residual cfd-math provider work remains in
  the public sparse storage alias and linear-solver/preconditioner traits:
  `nalgebra_sparse::CsrMatrix` and nalgebra `DVector` still need replacement
  by `leto_ops::CsrMatrix` and `leto::Array1`.
- [x] `cfd-math` [patch]: Move public sparse SpMV wrappers to Leto vectors.
  `spmv`, `spmv_parallel`, and `try_spmv` now accept `leto::Array1` input and
  output vectors. The only remaining nalgebra `DVector` SpMV bridge is private
  to `LinearOperator for CsrMatrix`, matching the still-unmigrated
  linear-solver trait boundary. Sparse tests, GMRES/AMG integration tests,
  interpolation quality checks, and the SpMV benchmark now use Leto arrays at
  the public sparse SpMV call sites. Evidence: cfd-math fmt/check, focused
  sparse/spmv/interpolation/AMG/solver nextest (40/40), no-default all-target
  clippy, public SpMV signature scan showing Leto `Array1`, and
  `git diff --check`. Residual cfd-math provider work remains in the public
  `LinearOperator`/solver/preconditioner trait family and the
  `nalgebra_sparse::CsrMatrix` storage alias.
- [x] `cfd-math` [patch]: Move `SparseMatrixExt` diagonal/scaling vector
  surfaces to Leto. `diagonal` now returns `leto::Array1`, and
  `set_diagonal`, `scale_rows`, and `scale_columns` now accept Leto arrays
  instead of nalgebra `DVector`. `JacobiPreconditioner::new` consumes the
  Leto diagonal through `Storage::as_slice`; its stored inverse diagonal
  remains nalgebra-owned only because `Preconditioner::apply_to` is still the
  public `DVector` boundary. Evidence: cfd-math fmt/check, focused sparse
  plus basic-preconditioner nextest (21/21), no-default all-target clippy,
  clean targeted `SparseMatrixExt` DVector-signature residue scan, and
  `git diff --check`. Residual cfd-math provider work remains in the public
  linear-solver/preconditioner trait boundary and the
  `nalgebra_sparse::CsrMatrix` storage alias.
- [x] `cfd-math` [patch]: Move multigrid smoother/cycle vector paths to
  Leto. `MultigridLevel`, `AMGHierarchy`, and `MultigridSmoother` now carry
  `leto::Array1` through the shared `MultigridVector` alias; Jacobi,
  Gauss-Seidel, symmetric Gauss-Seidel, SSOR, and Chebyshev smoothers compute
  residuals with a Leto-array `sparse::spmv_array` bridge; V/W/F cycles and
  coarsest solves now consume/return Leto vectors; and AMG V-cycle recursion
  converts to/from nalgebra `DVector` only at the current
  `Preconditioner::apply_to` boundary. Evidence: cfd-math fmt/check, focused
  multigrid cycle+smoother nextest (10/10), no-default all-target clippy,
  clean smoother/cycle provider-residue scan, and `git diff --check`.
  Residual cfd-math provider work remains in the public sparse/linear-solver
  API boundary: `nalgebra_sparse::CsrMatrix`, nalgebra `DVector`, and
  nalgebra `RealField` still need replacement by Leto/Eunomia contracts.
- [x] `cfd-math` [patch]: Move geometric multigrid (GMG) to Leto/Eunomia.
  `GeometricMultigrid`, `NonlinearOperator`, FAS/linear solves, Poisson
  hierarchy matrices, transfer operators, Jacobi relaxation, residuals, and
  tests now use `leto::Array2`/`Array1` and Eunomia scalar traits instead of
  nalgebra `DMatrix`/`DVector` and `num_traits` conversions. Evidence:
  cfd-math fmt/check, focused GMG nextest (5/5), no-default all-target clippy,
  clean GMG provider-residue scan, and `git diff --check`. Residual cfd-math
  provider work remains in the public linear-solver/preconditioner trait
  boundary.
- [x] `cfd-math` [patch]: Move GMRES internal dense workspace to
  Leto/Eunomia. `gmres::arnoldi` now stores the Krylov basis and Hessenberg
  matrix in `leto::Array2`, `gmres::givens` uses `leto::Array1`/`Array2` for
  rotation coefficients and least-squares state, and the solver workspace no
  longer allocates nalgebra dense matrices for the restarted Krylov subspace.
  `LinearSolverChain` carries the Eunomia real-field bound required by the
  migrated GMRES constructor. Evidence: cfd-math fmt/check, focused GMRES
  nextest (21/21), no-default all-target clippy, clean GMRES provider-residue
  scan for replaced providers, and `git diff --check`. Residual cfd-math
  provider work remains in the public linear-solver trait boundary:
  `LinearOperator`, `Preconditioner`, and `LinearSolver` still consume/return
  nalgebra `DVector` and carry nalgebra `RealField`.
- [x] `cfd-math` [patch]: Move the standalone AMG restriction dense-transfer
  utilities to Leto. `create_*_restriction`, `validate_restriction_operator`,
  `restrict_vector`, and `restrict_matrix` now operate on `leto::Array2`/
  `Array1`; Galerkin projection uses `leto_ops::MatrixProduct` for the dense
  `R * A * P` provider call. Evidence: cfd-math fmt/check, focused
  restriction nextest (7/7), no-default all-target clippy, and a clean
  restriction provider-residue scan. Residual cfd-math provider work remains
  in the broader sparse/linear-solver storage boundary.
- [x] `cfd-math` [patch]: Consolidate legacy CSR-to-dense solves behind a
  Leto dense bridge. `linear_solver::dense_bridge` now owns conversion from
  the current `nalgebra_sparse::CsrMatrix`/nalgebra `DVector` boundary into
  `leto::Array2`/`Array1` and solves with `leto_ops::solve`.
  `DirectSparseSolver` and multigrid cycle coarsest solves both consume that
  bridge instead of local nalgebra/manual dense LU paths. `LinearSolverChain`
  carries the Leto real-scalar bound required by the direct-solver provider
  path. Evidence: cfd-math fmt/check, focused direct-solver and multigrid
  cycle nextest (9/9), and no-default all-target clippy pass. Residual
  cfd-math provider work remains in the public sparse/linear-solver storage
  contracts and in the primary `rsparse` sparse-LU path.
- [x] `cfd-math` [patch]: Centralize sparse Leto/nalgebra CSR bridging and
  route builder construction through Leto provider validation.
  `SparseMatrixBuilder::{build,build_with_rhs,build_parallel}` and
  `ParallelAssembly::block_diagonal` now construct `leto_ops::CsrMatrix`
  first, then cross the remaining legacy boundary in one bridge module.
  Sparse assembly, patterns, and Schwarz subdomain extraction carry the Leto
  scalar bound required by the provider-backed construction path. Evidence:
  cfd-math check, focused sparse nextest (18/18), and no-default all-target
  clippy pass. Residual cfd-math provider work remains in the public
  sparse/linear-solver storage boundary: `SparseMatrix` still aliases
  `nalgebra_sparse::CsrMatrix` and vector APIs still use nalgebra `DVector`.
- [x] `cfd-math` [patch]: Consume Leto-owned CSR utility operations in
  `SparseMatrixExt`. Diagonal extraction, scalar/value scaling, row scaling,
  column scaling, Frobenius norm, diagonal dominance, and the
  condition-estimate heuristic now convert the current
  `nalgebra_sparse::CsrMatrix` bridge into `leto_ops::CsrMatrix` and call the
  provider surface instead of reimplementing CSR loops downstream. Evidence:
  cfd-math fmt/check, focused sparse nextest (18/18), and no-default
  all-target clippy pass. Residual cfd-math provider work remains in the
  sparse/linear-solver storage boundary: `nalgebra_sparse::CsrMatrix` and
  nalgebra `DVector` still need replacement by `leto_ops::CsrMatrix` and
  `leto::Array1`.
- [x] `cfd-math` [patch]: Consume Leto-ops SpGEMM in the sparse/AMG
  Galerkin path. `try_sparse_sparse_mul` now converts the current
  `nalgebra_sparse::CsrMatrix` boundary into `leto_ops::CsrMatrix`, delegates
  CSR×CSR products to `leto_ops::spgemm`, and maps provider or CSR validation
  failures into typed CFDrs errors. `try_sparse_transpose` now delegates CSR
  transpose to `leto_ops::CsrMatrix::transpose`, and AMG hierarchy setup uses
  it for restriction construction instead of `nalgebra_sparse::transpose_as_csc`.
  `try_spmv` now converts the current `DVector`/CSR boundary into Leto views
  and delegates SpMV to `leto_ops::spmv_into`; `spmv_parallel` no longer owns a
  separate CFDrs CSR traversal. Evidence: cfd-math fmt/check, focused sparse
  nextest (17/17), focused interpolation nextest (15/15), focused AMG nextest
  (6/6), and no-default all-target clippy pass. Residual cfd-math provider work
  remains in the sparse/linear-solver storage boundary:
  `nalgebra_sparse::CsrMatrix` and nalgebra `DVector` still need replacement by
  `leto_ops::CsrMatrix` and `leto::Array1`.
- [x] `leto-ops` provider [patch]: Add the CSR×CSR product needed by the
  `cfd-math` sparse/AMG migration. `leto_ops::spgemm` now owns sparse
  matrix-matrix products with sorted CSR output and exact-zero cancellation
  removal, and `CsrRow::nnz` covers sparse row-cardinality consumers. Evidence:
  upstream Leto-ops fmt/check/focused sparse nextest/clippy/doc all pass.
  Next CFDrs increment: consume `spgemm` while migrating `crates/cfd-math/src/sparse`
  and the AMG Galerkin path off `nalgebra_sparse`.
- [x] `cfd-math` [patch]: Move the SIMD helper cone to Leto/Eunomia.
  `SimdVectorOps` now targets `leto::Array1`, sparse SIMD matvec uses
  `leto_ops::CsrMatrix`/`spmv`, and generic SIMD/field/vectorization bounds
  use Eunomia scalar traits instead of nalgebra `RealField`. The SIMD
  integration test now builds its residual matrix through Leto-ops CSR parts
  instead of `nalgebra_sparse`. Evidence: cfd-math fmt/check, focused SIMD
  nextest (26/26, 318 skipped), no-default all-target clippy, and a clean
  SIMD provider-residue scan. Residual cfd-math provider work remains in
  sparse and linear-solver storage/solver surfaces.
- [x] `cfd-math` [patch]: Move the nonlinear-solver vector cone to Leto.
  JFNK configuration, solver state, JvP, restarted GMRES basis/residual work
  vectors, and tests now use `leto::Array1`/`Array2` and Eunomia scalar traits
  instead of nalgebra `DVector`/`DMatrix`/`RealField`. Anderson and JFNK share
  one local `nonlinear_solver::linalg` helper module for Leto vector norms,
  dot products, scaled updates, and small dense solve support. Evidence:
  cfd-math fmt/check, focused nonlinear nextest (9/9, 335 skipped),
  no-default all-target clippy, and a clean `src/nonlinear_solver`
  provider-residue scan. Residual cfd-math provider work remains in sparse
  and linear-solver storage/solver surfaces.
- [x] `cfd-math` [patch]: Move the high-order DG dense-array cone to Leto.
  `DGBasis`, `DGSolution`, numerical fluxes, `DGOperator`, DG limiters,
  `DGSolver`, and DG time integrators now use `leto::Array1`/`Array2` instead
  of nalgebra `DVector`/`DMatrix`. DG mass-matrix projection/derivative/RHS
  solves and implicit Newton correction solves now delegate to Leto dense
  solve helpers and return typed solver errors rather than silent fallbacks.
  The high-order DG docs and DG allocation/performance benches now compile
  against the migrated API. Evidence: cfd-math fmt/check, DG benchmark checks,
  focused DG nextest (62/62, 282 skipped), no-default all-target clippy,
  no-default doctests (3 passed, 3 ignored), and a clean provider-residue scan
  over `crates/cfd-math/src/high_order` plus the DG benches. Residual cfd-math
  provider work remains in sparse and linear-solver storage/solver surfaces.
- [x] `cfd-math` [patch]: Move the high-order spectral dense-array cone to
  Leto. `SpectralElement`,
  `SpectralMesh1D`, `SpectralDiffOp`, `SpectralInterp`,
  `SpectralQuadrature`, `SpectralFilter`, and spectral time-integration
  helpers now use `leto::Array1`/`Array2` instead of nalgebra
  `DVector`/`DMatrix`; spectral assembly now accepts and returns Leto dense
  arrays for local matrices, optional local RHS values, and debug dense CSR
  materialization. Local derivative, stiffness, dot, and matrix-vector paths
  share Leto helpers. `SpectralInterp::l2_projection` now surfaces Leto solve
  failures as typed solver errors instead of silently falling back to
  interpolation. Evidence: cfd-math fmt/check, spectral nextest (13/13, 331
  skipped), no-default all-target clippy, and a clean provider-residue scan
  over `crates/cfd-math/src/high_order/spectral`. Residual cfd-math provider
  work remains in sparse and linear-solver storage/solver surfaces.
- [x] `cfd-math` [patch]: Move the high-order WENO scalar-provider cone to
  Eunomia. `high_order::weno::{WENO5,WENO7,WenoReconstruction}` now use
  `eunomia::RealField`/`FloatElement` instead of nalgebra scalar traits or
  direct `num_traits` conversion bounds, while preserving the existing WENO5
  and WENO7 reconstruction formulas. Evidence: cfd-math fmt/check,
  WENO-focused nextest (6/6, 338 skipped), no-default all-target clippy, and a
  clean `crates/cfd-math/src/high_order/weno` provider-residue scan. Residual
  cfd-math provider work remains in sparse and linear-solver storage/solver
  surfaces.
- [x] `cfd-2d`/`cfd-3d`/`cfd-validation` [patch]: Move the next
  cross-crate boundary cone to crate-level Eunomia/Leto scalar seams.
  `Cfd2dScalar`, `Cfd3dScalar`, and `ValidationScalar` now centralize the
  migrated cfd-core boundary/fluid contracts instead of scattering
  nalgebra-only `RealField` bounds. cfd-2d NS-FVM and moving-wall boundary
  tests, cfd-3d FEM/DES and branch/venturi boundary setup, and
  cfd-validation Casson/Womersley/cross-fidelity benchmark surfaces now consume
  those seams and Leto boundary vectors. The LES Smagorinsky GPU path now
  recomputes SGS energy/dissipation diagnostics after GPU viscosity updates.
  Evidence: cfd-2d fmt/check/clippy library gates, cfd-2d GPU-feature lib
  nextest (518/518, 1 skipped), cfd-3d no-default library check/clippy, and
  cfd-validation no-default library check/clippy. Follow-up all-targets
  cleanup moved `examples/blood_venturi.rs` and
  `tests/simplec_pimple_validation.rs` off nalgebra boundary vectors/old
  scalar bounds and cleared the associated example/test lint blockers.
  Additional evidence: `cargo clippy -p cfd-2d --no-default-features
  --features gpu --all-targets -- -D warnings` and `cargo nextest run -p
  cfd-2d --no-default-features --features gpu --status-level fail` (572/572,
  27 skipped) pass. cfd-validation follow-up moved the geometry module and
  directly dependent 2D benchmark wrappers from nalgebra scalar/point contracts
  to Leto points/vectors plus Eunomia `RealField`; validation test boundary
  velocities targeting migrated cfd-core APIs now construct Leto vectors.
  Additional evidence: `cargo clippy -p cfd-validation --no-default-features
  --all-targets -- -D warnings` and `cargo nextest run -p cfd-validation
  --no-default-features geometry --status-level fail` (11/11, 420 skipped)
  pass, with a clean migrated geometry/benchmark residue scan. Analytical
  follow-up moved `src/analytical/**` and `src/solutions/mod.rs` from nalgebra
  `Vector3`/`RealField` trait surfaces to Leto vectors and Eunomia scalar
  contracts. Additional evidence: focused `cargo nextest run -p cfd-validation
  --no-default-features analytical --status-level fail` (20/20, 411 skipped)
  and a clean analytical/solutions residue scan. Error-metrics follow-up moved
  `src/error_metrics/**` to Eunomia scalar contracts and Leto `Vector3`
  magnitudes, with generic reductions routed through the crate-local Eunomia
  scalar helpers. Additional evidence: focused `cargo nextest run -p
  cfd-validation --no-default-features error_metrics --status-level fail`
  (21/21, 410 skipped), all-target clippy, and a clean error-metrics residue
  scan for nalgebra imports and old `T::zero()`/`T::one()` identities.
  Convergence follow-up moved `src/convergence/**` to Eunomia scalar
  contracts and crate-local scalar helpers; Richardson MMS result holders that
  embed `ConvergenceStudy<T>` now carry the required Eunomia field bound.
  Additional evidence: focused `cargo nextest run -p cfd-validation
  --no-default-features convergence --status-level fail` (27/27, 404 skipped),
  all-target clippy, and a clean convergence residue scan. Edge-case testing
  follow-up moved `src/edge_case_testing/**` from nalgebra scalar bounds to
  Eunomia `RealField`. Additional evidence: focused `cargo nextest run -p
  cfd-validation --no-default-features edge_case --status-level fail` (15/15,
  416 skipped), all-target clippy, and a clean edge-case/convergence/error-
  metrics residue scan.
  Manufactured follow-up moved `src/manufactured/**` from nalgebra scalar
  bounds and `T::zero()`/`T::one()` identities to Eunomia scalar contracts and
  crate-local scalar helpers; the remaining cfd-2d Poisson bridge requirement
  is localized as `cfd_2d::Cfd2dScalar`. Additional evidence: focused `cargo
  nextest run -p cfd-validation --no-default-features manufactured
  --status-level fail` (50/50, 381 skipped), all-target clippy, and a clean
  manufactured residue scan.
  Conservation follow-up moved `src/conservation/**` from nalgebra
  `DMatrix`/`DVector`/`RealField` to Leto `Array1`/`Array2` and Eunomia scalar
  contracts. Additional evidence: focused `cargo nextest run -p
  cfd-validation --no-default-features conservation --status-level fail`
  (18/18, 413 skipped), all-target clippy, and a clean conservation residue
  scan.
  Time-integration follow-up moved `cfd-math::time_stepping::stability` RK
  tableaus to Leto `Array2`/`Array1` plus Eunomia `RealField`, and moved
  `cfd-validation/src/time_integration/**` Euler/RK2/RK4 validation state to
  Leto `Array1`. Additional evidence: focused `cargo nextest run -p cfd-math
  --no-default-features stability --status-level fail` (5/5, 328 skipped),
  focused `cargo nextest run -p cfd-validation --no-default-features
  time_integration --status-level fail` (12/12, 419 skipped), all-target
  clippy for both crates, and a clean migrated-cone residue scan.
  cfd-math time-stepping follow-up moved `time_stepping::{traits,
  runge_kutta,adaptive,rk_chebyshev,exponential,imex}` to Leto-backed
  `TimeState<T>` and shared `TimeMatrix<T>` where dense matrices are required.
  Scalar contracts route through Eunomia, dense matrix exponential delegates to
  `leto_ops::matexp`, and IMEX Newton corrections now use `leto_ops::solve`
  instead of a nalgebra dense LU bridge. `rk4_bench`/`imex_bench` compile
  against the migrated API. Additional evidence: cfd-math fmt/check, focused
  nextest filters for `runge_kutta` (5/5), `adaptive` (3/3), `chebyshev`
  (5/5), `exponential` (6/6), and `imex` (5/5), no-default all-target clippy,
  and source/bench residue scans showing no nalgebra, `DMatrix`, `DVector`,
  `num_traits`, `num_complex`, or ndarray hits in the time-stepping cone.
  Differentiation follow-up moved `cfd-math::differentiation` from nalgebra
  `DVector`/`Vector3` result surfaces to Leto `Array1`/`Vector3`, kept scalar
  constants on Eunomia, and replaced the type-suffixed SIMD helper with
  `FiniteDifference<f32>::first_derivative_simd`. Additional evidence:
  cfd-math fmt/check, focused differentiation nextest (12/12), no-default
  all-target clippy, and a clean differentiation old-provider residue scan.
  Remaining cfd-math provider work is outside time stepping and
  differentiation.
- [x] `cfd-1d` [patch]: Move the 1D scalar boundary to a Eunomia-compatible
  provider contract. `Cfd1dScalar` is now the crate-local scalar seam carrying
  the remaining nalgebra linear-system constraint and the Eunomia
  `RealField`/numeric contract required by migrated cfd-core APIs. Domain,
  network, resistance, vascular, solver, transient, and analysis bounds now
  route through that seam; `NetworkDomain::contains_1d` uses Leto `Point1<T>`.
  Evidence: cfd-1d fmt, no-default library check, no-default nextest (725/725,
  3 skipped), no-default library clippy, and a residue scan showing nalgebra
  `RealField` only in `src/scalar.rs`. Residual work: downstream cfd-2d GPU
  consumer verification now reaches cfd-2d and is blocked by cfd-2d's own
  nalgebra scalar bounds around migrated cfd-core boundary/fluid contracts;
  cfd-1d all-target clippy also has unrelated pre-existing lint debt in
  tests/examples/cell-separation code.
- [x] `cfd-core` [patch]: Move GPU Poisson kernel dispatch and buffer transfer
  to Hephaestus. `GpuPoissonSolver` now owns Hephaestus
  `WgslMultiStorageKernel` instances for Jacobi, red-black, and residual WGSL
  entry points and uses `WgpuDevice` `ComputeDevice` buffers for uploads,
  work buffers, residual allocation, and downloads. The configured grid
  geometry is stored and validated instead of inferred from vector length.
  Evidence: cfd-core fmt, GPU-feature check, no-default check, GPU-feature
  all-target clippy, full GPU-feature nextest (231/231), and a clean
  `poisson_solver.rs` scan for raw WGPU allocation/pipeline/readback residue.
  Downstream cfd-2d GPU consumer verification now reaches cfd-2d and is still
  blocked by cfd-2d Eunomia/nalgebra trait-bound errors in the current
  workspace graph.
- [x] `cfd-core` [patch]: Move the GPU buffer/scalar provider cone to
  Hephaestus/Eunomia. `GpuBuffer<T>` now uses Hephaestus `WgpuDevice` through
  the `ComputeDevice` allocation, upload, download, and write contracts, and
  exposes raw WGPU only for existing bind-group wiring. GPU buffer, pipeline,
  and kernel scalar bounds now use `eunomia::RealField`. Evidence: cfd-core
  fmt, GPU-feature check, no-default check, GPU-feature all-target clippy,
  full no-default nextest (201/201), full GPU-feature nextest (231/231), and a
  clean `compute/{gpu,traits.rs}` residue scan for nalgebra scalar contracts,
  old identity calls, and `num_traits`. Residual work: CFDrs still owns WGSL
  CFD kernel/pipeline assembly; move reusable CFD GPU orchestration into
  Hephaestus in a follow-up rather than layering compatibility shims.
- [x] `cfd-core` [patch]: Move the mesh/staggered geometry storage cone to
  Leto/Eunomia. `Mesh<T>` now stores `leto::geometry::Point3<T>`, mesh
  transforms use Leto `Vector3<T>` and `FixedMatrix<T, 3, 3>`, mesh service/
  refinement/quality/statistics/connectivity bounds use `eunomia::RealField`,
  and `StaggeredGrid2D` no longer imports nalgebra. Upstream Leto provider gap
  closed: `FixedMatrix<T, 3, 3>` now multiplies `Vector3<T>` directly. Evidence:
  cfd-core fmt/check/no-default all-target clippy, full cfd-core no-default
  nextest (201/201), provider Leto fmt/check/clippy/full nextest (171/171), and
  clean mesh/staggered residue scans. Residual work: replace CFDrs-owned mesh
  topology with Gaia primitives in a follow-up instead of adding a compatibility
  adapter.
- [x] `cfd-core` [patch]: Move the geometry/domain, boundary, fluid-service,
  material, management, and solver factory scalar contracts in the migrated
  cone to Eunomia and Leto geometry. `Domain`, `BoundaryCondition`,
  `BoundaryGeometry`, `WallType`, inlet/outlet conditions, boundary
  applicators/managers/specs, fluid-dynamics service helpers, material database,
  management broadcast/conversion/factory, and solver trait/config companion
  modules now use `eunomia::RealField`; boundary/domain point and vector storage
  uses `leto::geometry` in the migrated cone. Upstream Leto provider gaps closed:
  `Point1<T>`, conditional `Eq` derives for fixed geometry values, and serde
  `std`/`alloc` feature propagation. Evidence: cfd-core fmt/check/no-default
  all-target clippy, full cfd-core no-default nextest (201/201), provider
  Leto check/clippy/nextest (170/170), and clean migrated-cone scans. Residual
  work remains outside this cone in mesh/GPU and other still-unmigrated provider
  surfaces.
- [x] `cfd-core` [patch]: Move the compute execution boundary to Eunomia and
  problem gravity to Leto geometry. `compute::{traits,cpu,dispatch}` now bind
  compute kernels, buffers, CPU backend storage, CPU advection, and runtime
  dispatcher methods on `eunomia::RealField`; CPU zero initialization and
  comparisons use Eunomia constants. `ProblemParameters::gravity` now stores
  `leto::geometry::Vector3<T>`. Evidence: cfd-core fmt/check/no-default
  all-target clippy, focused cfd-core compute/problem nextest (21/21), and
  touched-file scans. Residual work: `Problem<T>` still imports nalgebra's
  scalar trait because `Domain<T>` and `FluidTrait<T>` still require it; solver
  trait and convergence surfaces remain separate scalar-bound slices.
- [x] `cfd-core` [patch]: Move `abstractions::state` off nalgebra-owned field
  storage. `FieldData::Scalar` now uses `leto::Array1<T>`, vector fields use
  `leto::geometry::Vector3<T>`, and state scalar contracts use
  `eunomia::RealField`. The required upstream provider gap was closed in Leto
  by adding validated serde support for owned arrays/storage. Evidence:
  focused Leto serde nextest, Leto clippy, cfd-core fmt/check/no-default
  all-target clippy, focused cfd-core state nextest (3/3), and a clean
  `abstractions/state.rs` scan for nalgebra `DVector`, num-traits, ndarray,
  rayon, tokio, wgpu, and cuda. Residual work: `abstractions::problem`,
  `compute::{traits,cpu,gpu,solver}`, mesh, fluid, boundary, and material cones
  still own nalgebra/raw provider residue.
- [x] `cfd-core` [patch]: Move `compute::time` off nalgebra/num-traits.
  `TimeIntegrator` implementations now use `leto::Array1<T>` state and local
  slice-based scaled update, vector sum, and convergence-norm helpers over
  Eunomia scalar traits. Adaptive/variable time-step controllers use Eunomia
  constants and `FloatElement::powf`. Evidence: scoped fmt/check, cfd-core
  no-default library clippy, cfd-core default all-target clippy, focused
  cfd-core time nextest (17/17), no-default test build, and a clean
  `compute/time` scan for nalgebra `DVector`, ndarray, num-traits, rayon,
  tokio, wgpu, and cuda. The GPU benchmark feature gate is now corrected with
  `required-features = ["gpu"]`, so all-target no-default cfd-core clippy
  passes. Residual work: cfd-core still owns nalgebra in solver/problem, mesh,
  fluid, boundary, and GPU provider cones.
- [x] `cfd-core` [patch]: Move the solver configuration scalar boundary to
  Eunomia. `SolverConfig` and nested convergence/numerical/linear/network/
  builder config types in `cfd_core::compute::solver::config` now use
  `eunomia::RealField`; downstream holders in cfd-2d FDM, SIMPLE,
  pressure-velocity, vorticity-stream, cfd-3d FEM, and cfd-validation MMS
  Richardson were updated to consume that bound. `cfd-3d::spectral::solver`
  dropped the inherited `NalgebraRealField` alias. Evidence: scoped fmt,
  cfd-core/cfd-2d/cfd-3d/cfd-validation check and clippy, focused cfd-core
  `solver_config` nextest (2/2), focused cfd-2d FDM/pressure-velocity/SIMPLE/
  vorticity nextest (41/41), focused cfd-3d spectral/Poisson/FEM nextest
  (89/89), and static scans. Residual work: `SolverConfiguration` and
  `cfd_core::compute::solver::{traits,convergence,direct,iterative,monitor}`
  still use nalgebra because the solver trait contract depends on
  `abstractions::Problem<T>`.
- [x] `cfd-3d` [patch]: Complete the spectral Poisson Leto/Leto-ops provider
  slice. `PoissonSolver` now owns its global Laplacian and right-hand side as
  Leto arrays, assembles tensor-product operators with
  `leto_ops::MatrixProduct::kron`, applies boundary conditions in Leto storage,
  and solves through `leto_ops::MatrixSolve` instead of nalgebra
  `DMatrix`/`DVector` LU. `PoissonProblem` no longer requires nalgebra scalar
  bounds. Evidence: `cargo fmt -p cfd-3d --check`, `cargo check -p cfd-3d
  --no-default-features --lib`, `cargo check -p cfd-3d --no-default-features
  --test robustness_tests --test poisson_validation`, `cargo clippy -p
  cfd-3d --no-default-features --lib --test robustness_tests --test
  poisson_validation -- -D warnings`, and focused `cargo nextest run -p
  cfd-3d --no-default-features poisson --status-level fail` (12/12 passed).
  Residual work: `spectral::solver` still imports `NalgebraRealField` only
  because `cfd_core::compute::solver::SolverConfig` remains nalgebra-bound.
- [x] `cfd-3d` [patch]: Complete the spectral Chebyshev Leto/Eunomia provider
  slice. `ChebyshevPolynomial` owns `leto::Array2` differentiation matrices and
  applies first/second derivatives over `leto::Array1` inputs with typed
  dimension errors instead of nalgebra `DVector` multiplication. The
  Chebyshev scalar contract and `spectral::basis` now use `eunomia::RealField`.
  Direct Chebyshev consumers in co-located tests, `src/lib.rs`, and
  `tests/robustness_tests.rs` use Leto arrays. Evidence: `cargo fmt -p
  cfd-3d --check`, `cargo check -p cfd-3d --no-default-features --lib`,
  `cargo check -p cfd-3d --no-default-features --test robustness_tests --test
  poisson_validation`, `cargo clippy -p cfd-3d --no-default-features --lib
  --test robustness_tests --test poisson_validation -- -D warnings`, focused
  `cargo nextest run -p cfd-3d --no-default-features chebyshev --status-level
  fail` (20/20 passed), and focused `cargo nextest run -p cfd-3d
  --no-default-features poisson --status-level fail` (12/12 passed). The
  follow-on Poisson slice now closes the former nalgebra `DMatrix`/`DVector`
  global LU/Kronecker solve; the remaining `spectral::solver`
  `NalgebraRealField` bound is inherited from upstream
  `cfd_core::compute::solver::SolverConfig`.
- [x] `cfd-3d` [patch]: Complete the level-set Atlas provider slice.
  `level_set::{solver,advection,weno}` now uses `leto::geometry::Vector3` for
  velocity storage and transport and a level-set-local `LevelSetScalar` seam
  backed by `eunomia::RealField`/`FloatElement`. The dead crate-level
  `CfdScalar` trait that still pulled `nalgebra::RealField` into
  `crates/cfd-3d/src/scalar.rs` was removed. `tests/level_set_tests.rs` and the
  level-set sections of `tests/robustness_tests.rs` now construct Leto vectors.
  Evidence: `cargo fmt -p cfd-3d --check`, `cargo check -p cfd-3d
  --no-default-features --lib`, `cargo check -p cfd-3d --no-default-features
  --test level_set_tests`, `cargo check -p cfd-3d --no-default-features --test
  robustness_tests`, `cargo clippy -p cfd-3d --no-default-features --lib
  --test level_set_tests --test robustness_tests -- -D warnings`, focused
  `cargo nextest run -p cfd-3d --no-default-features level_set --status-level
  fail` (13/13 passed), and a clean forbidden-provider scan over
  `crates/cfd-3d/src/level_set` plus `crates/cfd-3d/tests/level_set_tests.rs`
  for `nalgebra`, `DMatrix`, `DVector`, `CfdScalar`, `crate::scalar`,
  `num_traits`, `num-traits`, `ndarray`, `rayon`, `tokio`, `rustfft`, `wgpu`,
  and `cuda`. Residual work: cfd-3d still owns legacy providers in FEM,
  spectral Chebyshev/Poisson, IBM, turbulence, validation, and non-level-set
  robustness-test surfaces.
- [x] `cfd-3d` [patch]: Complete the VOF Atlas provider slice.
  `vof::cavitation_solver` now uses `CavitationField = leto::Array2<f64>` for
  pressure, density, volume-fraction, inception, damage, bubble-radius, nuclei,
  and sonoluminescence fields; row-major Leto field offsets are explicit and
  separated from the VOF flat alpha offset.
  `tests/cavitation_solver_validation.rs` now constructs dense cavitation fields with Leto `Array2`, and
  `vof::scalar::VofScalar` now uses `eunomia::RealField` instead of
  `nalgebra::RealField`. Evidence: `cargo fmt -p cfd-3d --check`, `cargo check
  -p cfd-3d --no-default-features --lib`, `cargo check -p cfd-3d
  --no-default-features --test cavitation_solver_validation`, `cargo clippy -p
  cfd-3d --no-default-features --lib --test cavitation_solver_validation -- -D
  warnings`, focused `cargo nextest run -p cfd-3d --no-default-features
  cavitation --status-level fail` (23/23 passed), focused `cargo nextest run
  -p cfd-3d --no-default-features vof --status-level fail` (42/42 passed), and
  a clean forbidden-provider scan over `crates/cfd-3d/src/vof` plus the
  VOF/cavitation tests for `nalgebra`, `DMatrix`, `DVector`, `num_traits`,
  `ndarray`, `rayon`, `tokio`, `rustfft`, `wgpu`, and `cuda`. Residual work:
  cfd-3d legacy-provider ownership is now outside VOF, concentrated in FEM,
  spectral Chebyshev/Poisson, level-set, IBM, turbulence, and validation
  surfaces that still need Leto/Gaia/Eunomia/Hephaestus/Moirai provider slices.
- [x] `cfd-3d` [patch]: Migrate the VOF cavitation velocity-provider surface
  from nalgebra `Vector3` to Leto. `vof::{cavitation_solver,bubble_dynamics}`
  and `tests/cavitation_solver_validation.rs` now use
  `leto::geometry::Vector3` for cavitation velocity fields and bubble dynamics
  update calls. The cavitation solver now copies the Leto velocity slice
  directly into the Leto-owned `VofSolver` buffer instead of converting from
  nalgebra at the boundary. Evidence: `cargo fmt -p cfd-3d --check`, `cargo
  check -p cfd-3d --no-default-features --lib`, `cargo check -p cfd-3d
  --no-default-features --test cavitation_solver_validation`, `cargo clippy -p
  cfd-3d --no-default-features --lib --test cavitation_solver_validation -- -D
  warnings`, focused `cargo nextest run -p cfd-3d --no-default-features
  cavitation --status-level fail` (23/23 passed), and a clean direct nalgebra
  `Vector3` scan over the cavitation VOF cone. The dense cavitation field
  residue from this item is closed by the later Leto `CavitationField` slice;
  broader cfd-3d FEM,
  spectral Chebyshev/Poisson, level-set, and validation tests still need
  Leto/Gaia/Eunomia provider slices before cfd-3d can drop nalgebra.
- [x] `cfd-3d` [patch]: Migrate the non-cavitation VOF vector-provider surface
  from nalgebra `Vector3` to Leto. `vof::{solver,reconstruction,
  initialization,plic_geometry,advection}` and `tests/vof_tests.rs` now use
  `leto::geometry::Vector3`; robustness VOF call sites pass Leto vectors.
  Evidence: `cargo fmt -p cfd-3d --check`,
  `cargo check -p cfd-3d --no-default-features --test vof_tests`, `cargo check
  -p cfd-3d --no-default-features --test robustness_tests`, `cargo clippy -p
  cfd-3d --no-default-features --lib --test vof_tests -- -D warnings`,
  focused `cargo nextest run -p cfd-3d --no-default-features vof
  --status-level fail` (42/42 passed), and a clean direct nalgebra `Vector3`
  scan over the migrated VOF files. The cavitation dense-storage residue from
  this item is closed by the later Leto `CavitationField` slice; broader cfd-3d
  FEM,
  spectral Chebyshev/Poisson, level-set, and validation tests still need
  Leto/Gaia/Eunomia provider slices before cfd-3d can drop nalgebra.
- [x] `cfd-validation` [patch]: Close direct nalgebra `Vector2` ownership in
  validation source and tests. `manufactured::navier_stokes`,
  `conservation::{momentum,angular_momentum,mod}`,
  `benchmarks::vorticity_stream`, and MMS/conservation test consumers now use
  `leto::geometry::Vector2` and Leto indexed component access. Evidence: `cargo
  fmt -p cfd-validation --check`, `cargo check -p cfd-validation
  --no-default-features --tests`, `cargo clippy -p cfd-validation
  --no-default-features --lib -- -D warnings`, focused `cargo nextest run -p
  cfd-validation --no-default-features manufactured mms conservation taylor
  momentum angular --status-level fail` (104/104 passed), and a clean direct
  `nalgebra::Vector2` scan over `crates/cfd-validation/{src,tests}`. Residual
  work: cfd-validation still owns nalgebra dense `DMatrix`, 3D `Vector3`, and
  geometry point/vector surfaces pending Leto/Gaia provider migrations; broad
  `--tests` clippy remains blocked by unrelated baseline lints outside this
  provider slice.
- [x] `cfd-2d` [patch]: Close direct nalgebra `Vector2` ownership across
  source, tests, examples, and benches. `physics::vorticity_stream`,
  `piso_algorithm::corrector`, `solvers::lbm::solver`,
  `physics::immersed_boundary`, `examples/blood_venturi.rs`, and
  `benches/solver_benchmarks.rs` now use `leto::geometry::Vector2`; component
  access uses Leto `[0]`/`[1]` indexing instead of nalgebra `.x`/`.y` fields.
  The downstream `cfd-validation::benchmarks::vorticity_stream` consumer now
  accepts the Leto vector field returned by `VorticityStreamSolver`. Evidence:
  `cargo fmt -p cfd-2d -p cfd-validation --check`, `cargo check -p cfd-2d
  --no-default-features --examples --benches`, `cargo check -p cfd-validation
  --no-default-features`, focused `cargo clippy -p cfd-2d
  --no-default-features --example blood_venturi --bench solver_benchmarks -- -D
  warnings`, `cargo clippy -p cfd-validation --no-default-features --lib -- -D
  warnings`, focused `cargo nextest run -p cfd-2d --no-default-features
  vorticity corrector lbm immersed --status-level fail` (44/44 passed), and a
  clean direct `nalgebra::Vector2` scan over cfd-2d source/tests/examples/benches
  plus the touched cfd-validation benchmark. Residual work: direct cfd-2d
  nalgebra ownership remains for scalar `RealField`, boundary `Vector3`,
  `DVector`/`DMatrix`, and nalgebra-sparse storage/provider surfaces.
- [x] `cfd-2d` [patch]: Migrate the pressure-velocity/SIMPLEC/PIMPLE vector
  workspace family from nalgebra `Vector2` to Leto. `fields.rs`,
  `pressure_velocity::{solver,pressure,rhie_chow}`, `simplec_pimple::{
  algorithms,diagnostics,interpolation,pimple,simplec,solver}`,
  `physics::momentum::interpolation` reference tests, and
  `tests/simplec_pimple_validation.rs` now use `leto::geometry::Vector2` for
  velocity workspaces, Rhie-Chow caches, correction buffers, face velocities,
  and validation setup. Component access uses Leto `[0]`/`[1]` indexing instead
  of nalgebra `.x`/`.y` fields. Evidence: `cargo fmt -p cfd-2d --check`,
  `cargo check -p cfd-2d --no-default-features`, `cargo clippy -p cfd-2d
  --no-default-features --lib -- -D warnings`, focused `cargo nextest run -p
  cfd-2d --no-default-features pressure_velocity simplec pimple rhie
  --status-level fail` (29/29 passed), `git diff --check`, and a clean direct
  scan over the touched family for `nalgebra::Vector2` imports/usages. Residual
  work: this family still has nalgebra scalar `RealField`, `Vector3` boundary
  inputs, and `DVector` linear-solver storage pending upstream cfd-core/cfd-math
  scalar/vector and Leto linear-solver migrations.
- [x] `cfd-2d` [patch]: Move the Ghia pressure-correction helper test off
  direct nalgebra vector construction. `crates/cfd-2d/tests/
  ghia_cavity_simplec_validation.rs` now uses `leto::geometry::Vector2` for the
  local divergence-free velocity helper instead of `nalgebra::Vector2`.
  Evidence: `cargo fmt -p cfd-2d --check`, `cargo check -p cfd-2d
  --no-default-features --test ghia_cavity_simplec_validation`, focused `cargo
  nextest run -p cfd-2d --no-default-features test_pressure_correction_basic
  --status-level fail` (1/1 passed), and a clean direct nalgebra-provider scan
  over the test file. Residual work: `tests/simplec_pimple_validation.rs` and
  the SIMPLEC/PIMPLE plus pressure-velocity production vector workspaces still
  own nalgebra vector/scalar surfaces.
- [x] `cfd-2d` [patch]: Migrate the turbulence validation scalar contract from
  direct nalgebra bounds to Eunomia. `crates/cfd-2d/src/physics/turbulence/
  validation/{mod.rs,rans.rs,les_des.rs,benchmarks.rs}` no longer imports or
  bounds `nalgebra::RealField`; validation scalar contracts use
  `eunomia::RealField`/`FloatElement`, and validation arrays/vectors stay on
  Leto `Array2`/`Vector2`. Evidence: `cargo fmt -p cfd-2d --check`, `cargo
  check -p cfd-2d --no-default-features`, `cargo clippy -p cfd-2d
  --no-default-features --lib -- -D warnings`, focused `cargo nextest run -p
  cfd-2d --no-default-features turbulence validation --status-level fail`
  (207/207 passed), and a clean direct nalgebra-provider scan over the
  validation directory. Residual work: cfd-2d still directly declares
  `nalgebra`/`nalgebra-sparse` because other modules and tests still own
  `RealField`, `DVector`, `DMatrix`, `Vector2`, `Vector3`, and sparse matrix
  provider surfaces.
- [x] `cfd-2d` [patch]: Remove direct package `num-traits` ownership. The
  closing pass migrated `physics::turbulence::validation::{mod,les_des,benchmarks}`
  from `num_traits::{FromPrimitive,ToPrimitive}` bounds to
  `eunomia::{FloatElement,RealField}`, removed the immersed-boundary test's
  unnecessary `ToPrimitive` conversion over an already-`f64` field, removed the
  concrete `f64` Ghia validation `num_traits` imports/bounds, and deleted
  `num-traits.workspace = true` from `crates/cfd-2d/Cargo.toml`. Evidence:
  `cargo fmt -p cfd-2d --check`, `cargo check -p cfd-2d
  --no-default-features --tests`, `cargo clippy -p cfd-2d
  --no-default-features --lib -- -D warnings`, focused `cargo nextest run -p
  cfd-2d --no-default-features turbulence validation immersed ghia
  --status-level fail` (213/213 passed), `cargo tree -p cfd-2d --depth 1`
  showing no direct `num-traits`, `git diff --check`, and a clean direct scan
  over `crates/cfd-2d/{src,tests,Cargo.toml}`. This supersedes older
  cfd-2d slice-local residual notes that said the manifest still owned
  `num-traits`; remaining cfd-2d provider work is now nalgebra/nalgebra-sparse,
  Leto storage/vector migration, and downstream transitive provider cleanup.
- [x] `cfd-2d` [patch]: Migrate the momentum setup/boundary scalar-provider
  seam to Eunomia. `physics::momentum::setup` no longer imports or bounds
  direct `num_traits::{FromPrimitive,ToPrimitive}`. `physics::momentum::boundary`
  and `boundary::directional` now route zero/one row entries, zero-gradient
  comparisons, quadratic wall-extrapolation constants, and corner consistency
  absolute-value checks through `crates/cfd-2d/src/scalar.rs` and
  `eunomia::FloatElement` instead of direct `num_traits::FromPrimitive`,
  `T::from_f64(...).unwrap_or_else`, `T::zero()`, `T::one()`, or scalar
  `.abs()` usage. Evidence: `cargo fmt -p cfd-2d --check`, `cargo check -p
  cfd-2d --no-default-features`, `cargo clippy -p cfd-2d
  --no-default-features --lib -- -D warnings`, focused `cargo nextest run -p
  cfd-2d --no-default-features boundary --status-level fail` (42/42 passed),
  `git diff --check`, and a clean direct-provider scan over the touched files.
  Residual work: direct `num-traits` remains in Ghia validation tests,
  immersed-boundary test code, turbulence validation, and f64-only/test scalar
  surfaces, so the cfd-2d manifest still owns `num-traits`.
- [x] `cfd-2d` [patch]: Migrate the momentum interpolation scalar-provider
  seam to Eunomia. `physics::momentum::interpolation` no longer imports or
  bounds direct `num_traits::FromPrimitive` and no longer calls
  `T::from_f64(...).unwrap_or_else`, `T::zero()`, or `T::one()` in the touched
  Rhie-Chow path. The tiny diagonal floor, zero initialization for face arrays,
  and harmonic face coefficient factor route through
  `crates/cfd-2d/src/scalar.rs` and `eunomia::FloatElement`. Evidence: `cargo
  fmt -p cfd-2d --check`, `cargo check -p cfd-2d --no-default-features`,
  `cargo clippy -p cfd-2d --no-default-features --lib -- -D warnings`,
  focused `cargo nextest run -p cfd-2d --no-default-features
  coefficient_aware_interpolation_matches_exact_reference --status-level fail`
  (1/1 passed), `git diff --check`, and a clean direct-provider scan over the
  touched file. Residual work: direct `num-traits` remains in Ghia validation
  tests, immersed-boundary test code, turbulence validation, and f64-only/test
  scalar surfaces, so the cfd-2d manifest still owns `num-traits`.
- [x] `cfd-2d` [patch]: Migrate the problem/streamtube scalar-provider seam to
  Atlas providers. `problem.rs` now stores incompressible problem and solution
  velocity fields with `leto::geometry::Vector2` and routes initial pressure,
  velocity-magnitude maxima, and pressure maxima through
  `crates/cfd-2d/src/scalar.rs`/Eunomia instead of local nalgebra vector
  storage or direct `T::zero()` folds. `physics::streamtube::partitioning` no
  longer imports or bounds direct `num_traits::{Float,FromPrimitive}` and now
  routes constants, absolute values, and square roots through
  `eunomia::{FloatElement,NumericElement}` plus the crate-local scalar adapter.
  Evidence: `cargo fmt -p cfd-2d --check`, `cargo check -p cfd-2d
  --no-default-features`, `cargo clippy -p cfd-2d --no-default-features --lib
  -- -D warnings`, focused `cargo nextest run -p cfd-2d
  --no-default-features problem streamtube separating --status-level fail`
  (4/4 passed), `git diff --check`, and a clean direct-provider scan over both
  touched files. Residual work: `problem.rs` still inherits
  `nalgebra::RealField` from `cfd-core` boundary/fluid types, and direct
  `num-traits` remains in immersed-boundary tests, momentum
  setup/interpolation/boundary, turbulence validation, and f64-only/test scalar
  surfaces, so the cfd-2d manifest still owns `num-traits`.
- [x] `cfd-2d` [patch]: Migrate the FVM scalar and vector-provider surface to
  Atlas providers. `solvers::fvm::{config,flux,geometry,solver}` no longer
  imports or bounds direct `num_traits`, calls `T::from_*`, `.to_f64()`,
  `T::zero()`, `T::one()`, or scalar `.abs()` in the touched APIs and module
  tests, and no longer imports `nalgebra` for FVM face/solver vectors. FVM
  defaults, face generation constants, residual magnitudes, max terms, flux
  Peclet calculations, nonfinite diffusion validation, and module test
  absolute-value assertions route through `crates/cfd-2d/src/scalar.rs` and
  `eunomia::{FloatElement,NumericElement}`. `geometry` and `solver` now use
  `leto::geometry::Vector2`; power-law and hybrid flux calculators preserve
  the diffusion coefficient in `T` rather than narrowing through `f64`.
  Evidence: `cargo fmt -p cfd-2d --check`, `cargo check -p cfd-2d
  --no-default-features`, `cargo clippy -p cfd-2d --no-default-features --lib
  -- -D warnings`, focused `cargo nextest run -p cfd-2d
  --no-default-features fvm --status-level fail` (25/25 passed), `git diff
  --check`, and a clean direct-provider scan over
  `crates/cfd-2d/src/solvers/fvm`. Residual work: direct `num-traits` still
  remains outside this slice in problem setup, immersed-boundary tests,
  streamtube partitioning, momentum, energy, turbulence, and f64-only/test
  scalar surfaces, so the cfd-2d manifest still owns `num-traits`.
- [x] `cfd-2d` [patch]: Migrate the vorticity-stream scalar-provider surface
  to Eunomia. `physics::vorticity_stream` no longer imports or bounds direct
  `num_traits::FromPrimitive` or calls `T::from_*`, `T::zero()`, `T::one()`,
  or scalar `.abs()` in the touched APIs and tests. Default tolerances,
  timestep/SOR constants, zero/one identities, Laplacian gradient factors, SOR
  residuals, upwind sign checks, velocity-recovery denominators,
  boundary-vorticity factors, convergence checks, and continuity-test
  absolute-value assertions route through `crates/cfd-2d/src/scalar.rs` and
  `eunomia::{FloatElement,NumericElement}`. Evidence: `cargo fmt -p cfd-2d
  --check`, `cargo check -p cfd-2d --no-default-features`, `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings`, focused `cargo nextest
  run -p cfd-2d --no-default-features vorticity --status-level fail` (4/4
  passed), `git diff --check`, and a clean direct-provider scan over
  `crates/cfd-2d/src/physics/vorticity_stream.rs`. Residual work: direct
  `num-traits` still remains outside this slice in problem setup, FVM,
  immersed-boundary tests, streamtube partitioning, momentum, energy,
  turbulence, and f64-only/test scalar surfaces, so the cfd-2d manifest still
  owns `num-traits`.
- [x] `cfd-2d` [patch]: Migrate the acoustic-drift scalar-provider surface to
  Eunomia. `physics::acoustics::gorkov` and `solvers::drift_diffusion_2d` no
  longer import or bound direct `num_traits::{Float,FromPrimitive}` or call
  `T::from_*`, `T::zero()`, `T::one()`, `Float::`, or scalar `.abs()` in the
  touched APIs. ARF material constants, compressibility identities,
  contrast-factor constants, standing-wave sine evaluation, drift-diffusion
  zero/one identities, under-relaxation constants, Patankar upwind max terms,
  pivot-threshold checks, and convergence residual magnitudes route through
  `crates/cfd-2d/src/scalar.rs` and
  `eunomia::{FloatElement,NumericElement}`. Evidence: `cargo fmt -p cfd-2d
  --check`, `cargo check -p cfd-2d --no-default-features`, `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings`, focused `cargo nextest
  run -p cfd-2d --no-default-features gorkov --status-level fail` (4/4
  passed), focused `cargo nextest run -p cfd-2d --no-default-features
  drift_diffusion --status-level fail` (1/1 passed), `git diff --check`, and a
  clean direct-provider scan over both files. Residual work: direct
  `num-traits` still remains outside this slice in problem setup, FVM,
  vorticity-stream, immersed-boundary tests, streamtube partitioning, momentum,
  energy, turbulence, and f64-only/test scalar surfaces, so the cfd-2d manifest
  still owns `num-traits`.
- [x] `cfd-2d` [patch]: Migrate the TVD scheme scalar-provider family to
  Eunomia. `schemes::tvd::{mod,muscl,quick}` no longer imports or bounds
  direct `num_traits::{FromPrimitive,ToPrimitive}` or calls `T::from_*`,
  `T::zero()`, `T::one()`, `.to_f64()`, or scalar `.abs()` in the touched APIs.
  TVD limiter identities, min/max/absolute-value operations, MUSCL slope
  ratios, QUICK interpolation constants, and MUSCL2/MUSCL3 boundary constants
  route through `crates/cfd-2d/src/scalar.rs` and
  `eunomia::{FloatElement,NumericElement}`. `FluxLimiter::apply` stays generic
  in `T` instead of round-tripping the limiter ratio through `f64`. Evidence:
  `cargo fmt -p cfd-2d --check`, `cargo check -p cfd-2d
  --no-default-features`, `cargo clippy -p cfd-2d --no-default-features --lib
  -- -D warnings`, focused `cargo nextest run -p cfd-2d
  --no-default-features tvd --status-level fail` (32/32 passed), focused
  `cargo nextest run -p cfd-2d --no-default-features muscl --status-level
  fail` (14/14 passed), `git diff --check`, and a clean direct-provider scan
  over `crates/cfd-2d/src/schemes/tvd`. Residual work: direct `num-traits`
  still remains outside this slice in problem setup, drift-diffusion, FVM,
  vorticity-stream, acoustics, immersed-boundary tests, streamtube
  partitioning, momentum, energy, turbulence, and f64-only/test scalar
  surfaces, so the cfd-2d manifest still owns `num-traits`.
- [x] `cfd-2d` [patch]: Migrate the LBM scalar-provider surface to Eunomia.
  `solvers::lbm` no longer imports or bounds direct
  `num_traits::{Float,FromPrimitive}` or calls `T::from_*`, `T::zero()`,
  `T::one()`, `Float::`, or scalar `.abs()` in the touched APIs and module
  tests. Solver defaults, low-Mach validation, grid-index conversion,
  boundary density/velocity reconstruction, passive-scalar boundary
  zeroing, streaming push-zeroing, macroscopic/lattice assertions,
  Carreau-Yasuda assertions, and Shan-Chen pseudopotential/force constants
  route through `crates/cfd-2d/src/scalar.rs` and
  `eunomia::{FloatElement,NumericElement,CastFrom}`. Evidence: `cargo fmt -p
  cfd-2d --check`, `cargo check -p cfd-2d --no-default-features`, `cargo
  clippy -p cfd-2d --no-default-features --lib -- -D warnings`, focused
  `cargo nextest run -p cfd-2d --no-default-features lbm --status-level fail`
  (31/31 passed), `git diff --check`, and a clean direct-provider scan over
  `crates/cfd-2d/src/solvers/lbm`. Residual work: direct `num-traits` still
  remains outside this slice in problem setup, drift-diffusion, FVM,
  vorticity-stream, immersed-boundary tests, streamtube partitioning, TVD
  schemes, momentum, energy, turbulence, and f64-only/test scalar surfaces, so
  the cfd-2d manifest still owns `num-traits`.
- [x] `cfd-2d` [patch]: Migrate the SIMPLE scalar-provider surface to Eunomia.
  `solvers::simple::{algorithm,momentum,pressure}` no longer imports or bounds
  direct `num_traits::{FromPrimitive,ToPrimitive}` or calls `T::from_*`,
  `T::zero()`, `T::one()`, scalar `.abs()`, or direct finite checks in the
  touched APIs and module tests. Patankar default constants, workspace zero
  initialization, stagnant-cell safeguards, pressure-correction identities,
  neighbor-count conversion, finite checks, and module test absolute-value
  assertions route through `crates/cfd-2d/src/scalar.rs` and
  `eunomia::{FloatElement,NumericElement}`. Evidence: `cargo fmt -p cfd-2d
  --check`, `cargo check -p cfd-2d --no-default-features`, `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings`, focused `cargo nextest
  run -p cfd-2d --no-default-features simple --status-level fail` (19/19
  passed), `git diff --check`, and a clean direct-provider scan over
  `crates/cfd-2d/src/solvers/simple`. Residual work: direct `num-traits` still
  remains outside this slice in vorticity-stream, drift-diffusion,
  acoustics, momentum, turbulence, energy, streamtube/problem/schemes, FVM, and
  f64-only/test scalar surfaces, so the cfd-2d manifest still owns
  `num-traits`.
- [x] `cfd-2d` [patch]: Migrate the Venturi scalar-provider surface to
  Eunomia. `solvers::venturi_flow::{mod,solver}` no longer imports or bounds
  direct `num_traits::{Float,FromPrimitive,ToPrimitive}` or calls
  `T::from_*`, `.to_f64()`, `T::zero()`, `T::one()`, `Float::`, or scalar
  `.abs()` in the touched APIs and module tests. ISO/default constants, beta
  clamping, domain bounds, analytical Bernoulli/viscous constants,
  energy-dissipation guards, stretched-grid index conversion, generic
  sine/min/max dispatch, diagnostic f64 conversion, throat-column selection,
  inlet/outlet/throat reductions, and validation error metrics route through
  `crates/cfd-2d/src/scalar.rs` and `eunomia::{FloatElement,NumericElement}`.
  Evidence: `cargo fmt -p cfd-2d --check`, `cargo check -p cfd-2d
  --no-default-features`, `cargo clippy -p cfd-2d --no-default-features --lib
  -- -D warnings`, focused `cargo nextest run -p cfd-2d --no-default-features
  venturi --status-level fail` (13/13 passed), `git diff --check`, and a
  clean direct-provider scan over `crates/cfd-2d/src/solvers/venturi_flow`.
  Residual work: direct `num-traits` still remains outside this slice in
  vorticity-stream, LBM, drift-diffusion, acoustics, momentum,
  turbulence, energy, streamtube/problem/schemes, FVM, and f64-only/test scalar
  surfaces, so the cfd-2d manifest still owns `num-traits`.
- [x] `cfd-2d` [patch]: Migrate the branching-flow scalar-provider family to
  Eunomia. `solvers::{bifurcation_flow,n_furcation_flow}` no longer imports or
  bounds direct `num_traits::{Float,FromPrimitive,ToPrimitive}` or calls
  `T::from_*`, `T::zero()`, or `Float::` in the touched APIs and module tests.
  Geometry constants, branch-index conversion, zero identities, generic
  sin/cos dispatch, flux accumulation, and mass-balance absolute-value
  normalization route through `crates/cfd-2d/src/scalar.rs` and
  `eunomia::{FloatElement,NumericElement}`. Evidence: `cargo fmt -p cfd-2d
  --check`, `cargo check -p cfd-2d --no-default-features`, `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings`, focused `cargo nextest
  run -p cfd-2d --no-default-features bifurcation --status-level fail` (9/9
  passed), focused `cargo nextest run -p cfd-2d --no-default-features
  n_furcation --status-level fail` (1/1 passed), `git diff --check`, and a
  clean direct-provider scan over both branching files. Residual work: direct
  `num-traits` still remains outside this slice in vorticity-stream, LBM,
  drift-diffusion, acoustics, momentum, turbulence, energy,
  streamtube/problem/schemes, FVM, and f64-only/test scalar surfaces, so the
  cfd-2d manifest still owns `num-traits`.
- [x] `cfd-2d` [patch]: Migrate the cross-junction scalar-provider surface to
  Eunomia. `solvers::cross_junction_flow` no longer imports or bounds direct
  `num_traits::{Float,FromPrimitive,ToPrimitive}` or calls `T::from_*`,
  `T::zero()`, or scalar `.abs()` in the touched APIs and module tests.
  Geometry constants, flux accumulation identities, mass-balance absolute-value
  normalization, and bounding-box assertions route through
  `crates/cfd-2d/src/scalar.rs` and `eunomia::{FloatElement,NumericElement}`.
  Evidence: `cargo fmt -p cfd-2d --check`, `cargo check -p cfd-2d
  --no-default-features`, `cargo clippy -p cfd-2d --no-default-features --lib
  -- -D warnings`, focused `cargo nextest run -p cfd-2d
  --no-default-features cross_junction --status-level fail` (5/5 passed),
  `git diff --check`, and a clean direct-provider scan over
  `crates/cfd-2d/src/solvers/cross_junction_flow.rs`. Residual work: direct
  `num-traits` still remains outside this slice in LBM,
  drift-diffusion, and f64-only/test scalar surfaces, so
  the cfd-2d manifest still owns `num-traits`.
- [x] `cfd-2d` [patch]: Migrate the Poiseuille scalar-provider surface to
  Eunomia. `solvers::poiseuille::{mod,numerics}` no longer imports or bounds
  direct `num_traits::{Float,FromPrimitive}` or calls `T::from_*`,
  `T::zero()`, `T::one()`, or scalar `.abs()` in the touched APIs and module
  tests. Configuration defaults, grid-index conversion, identities,
  tridiagonal workspaces, harmonic-mean constants, shear-rate magnitudes,
  viscosity residual norms, Thomas pivot thresholds, and error diagnostics
  route through `crates/cfd-2d/src/scalar.rs` and
  `eunomia::{FloatElement,NumericElement}`. Evidence: `cargo fmt -p cfd-2d
  --check`, `cargo check -p cfd-2d --no-default-features`, `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings`, focused `cargo nextest
  run -p cfd-2d --no-default-features poiseuille --status-level fail` (7/7
  passed), `git diff --check`, and a clean direct-provider scan over
  `crates/cfd-2d/src/solvers/poiseuille`. Residual work: direct `num-traits`
  still remains outside this slice in vorticity-stream, LBM solvers, acoustics,
  momentum, turbulence, physics, and tests, so
  the cfd-2d manifest still owns `num-traits`.
- [x] `cfd-2d` [patch]: Migrate serpentine and scalar transport scalar-provider
  surfaces to Eunomia. `solvers::serpentine_flow::{mod,solver}` and
  `solvers::scalar_transport_2d` no longer import or bound direct
  `num_traits::{Float,FromPrimitive,ToPrimitive}` or call `T::from_*`,
  `T::zero()`, `T::one()`, or scalar `.abs()` in the touched APIs and module
  tests. Constants, identities, Peclet/mixing calculations, Fourier-mode index
  conversion, transport relaxation constants, coefficient max/abs operations,
  outlet averaging, and validator diagnostics route through
  `crates/cfd-2d/src/scalar.rs` and
  `eunomia::{FloatElement,NumericElement}`. Evidence: `cargo fmt -p cfd-2d
  --check`, `cargo check -p cfd-2d --no-default-features`, `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings`, focused `cargo nextest
  run -p cfd-2d --no-default-features serpentine --status-level fail` (5/5
  passed), `git diff --check`, and a clean direct-provider scan over
  `crates/cfd-2d/src/solvers/serpentine_flow` plus
  `crates/cfd-2d/src/solvers/scalar_transport_2d.rs`. The `scalar_transport`
  nextest filter has no matching tests, so transport-specific evidence is
  compile/clippy plus serpentine integration coverage. Residual work: direct
  `num-traits` still remains outside this slice in vorticity-stream,
  LBM solvers, acoustics, momentum,
  turbulence, physics, and tests, so the cfd-2d manifest still owns
  `num-traits`.
- [x] `cfd-2d` [patch]: Migrate the FDM scalar-provider surface to Eunomia.
  `solvers::fdm::{advection_diffusion,config,diffusion,linear_solver,poisson}`
  no longer imports or bounds direct `num_traits::FromPrimitive` or calls
  `T::from_*`, `T::zero()`, `T::one()`, or scalar `.abs()` in the touched FDM
  APIs and MMS checks. Constants, identities, singular-diagonal guards,
  Gauss-Seidel residual magnitudes, and finite-difference source/test
  conversions route through `crates/cfd-2d/src/scalar.rs` and
  `eunomia::{FloatElement,NumericElement}`. The existing FDM MMS source tests
  are now wired into the test graph, and the advection-diffusion stencil sign
  now matches the documented steady operator. Evidence: `cargo fmt -p cfd-2d
  --check`, `cargo check -p cfd-2d --no-default-features`, `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings`, focused `cargo nextest
  run -p cfd-2d --no-default-features fdm --status-level fail` (2/2 passed),
  `git diff --check`, and a clean direct-provider scan over
  `crates/cfd-2d/src/solvers/fdm`. Residual work: direct `num-traits` still
  remains outside this slice in vorticity-stream, LBM solvers, acoustics,
  momentum, turbulence, physics, and
  tests, so the cfd-2d manifest still owns `num-traits`.
- [x] `cfd-2d` [patch]: Migrate the network scalar-provider surface to
  Eunomia. `network::{build,channel,coupled,postprocess,projection,reference,
  solve,types}` no longer imports or bounds direct
  `num_traits::{Float,FromPrimitive,ToPrimitive}` or calls `T::from_*`,
  `.to_f64()`, `T::zero()`, or `T::one()` in the touched APIs. Constants,
  identities, abs/min/max, finite checks, Anderson relaxation constants,
  projection summaries, channel diagnostics, and f64 reporting route through
  `crates/cfd-2d/src/scalar.rs` and
  `eunomia::{FloatElement,NumericElement}`. Evidence: `cargo fmt -p cfd-2d
  --check`, `cargo check -p cfd-2d --no-default-features`, `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings`, focused `cargo nextest
  run -p cfd-2d --no-default-features network --status-level fail` (22/22
  passed), `git diff --check`, and a clean direct-provider scan over
  `crates/cfd-2d/src/network`. Residual work: direct `num-traits` still
  remains outside this slice in drift-diffusion, LBM/FVM, turbulence,
  physics, and other solver/test surfaces, so the cfd-2d manifest still owns
  `num-traits`.
- [x] `cfd-2d` [patch]: Migrate the NS-FVM scalar-provider surface to
  Eunomia. `solvers::ns_fvm::{field,solver}` and the momentum, pressure, and
  velocity-interpolation submodules no longer import or bound direct
  `num_traits::{Float,FromPrimitive,ToPrimitive}` or call `T::from_*`,
  `.to_f64()`, `T::zero()`, or `T::one()` in the touched APIs. Constants,
  identities, abs/min/max/sqrt/finite checks, interpolation grid-index
  conversion, pressure-correction thresholds, turbulence seeds, and diagnostic
  formatting route through `crates/cfd-2d/src/scalar.rs` and
  `eunomia::{FloatElement,NumericElement}`. Evidence: `cargo fmt -p cfd-2d
  --check`, `cargo check -p cfd-2d --no-default-features`, `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings`, focused `cargo nextest
  run -p cfd-2d --no-default-features ns_fvm --status-level fail` (2/2
  passed), `git diff --check`, and a clean direct-provider scan over
  `crates/cfd-2d/src/solvers/ns_fvm`. Residual work: direct `num-traits` still
  remains outside this slice in broader cfd-2d network, turbulence,
  LBM/FDM/FVM, solver, physics, and other test surfaces, so the cfd-2d
  manifest still owns `num-traits`.
- [x] `cfd-2d` [patch]: Migrate the continuity residual SSOT to Eunomia.
  `crates/cfd-2d/src/solvers/continuity.rs` no longer uses direct
  `num_traits::FromPrimitive`, `T::from_*`, `T::zero()`, or method `.abs()` in
  the shared forward/central/face continuity residual helpers. Zero values and
  constants route through `crates/cfd-2d/src/scalar.rs`; absolute-value
  reductions use `eunomia::NumericElement`. Evidence: `cargo fmt -p cfd-2d
  --check`, `cargo check -p cfd-2d --no-default-features`, `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings`, focused `cargo nextest
  run -p cfd-2d --no-default-features continuity --status-level fail` (6/6
  passed), `git diff --check`, and a clean direct-provider scan over the
  touched module. Residual work: direct `num-traits` still remains across
  broader cfd-2d network, turbulence, LBM/FDM/FVM/NS-FVM, solver, physics, and
  other test surfaces, so the cfd-2d manifest still owns `num-traits`.
- [x] `cfd-2d` [patch]: Migrate the SIMPLEC/PIMPLE scalar-provider surface to
  Eunomia. `simplec_pimple::{config,algorithms,diagnostics,solver,
  interpolation,simplec,pimple}` and the focused validation test no longer use
  direct `num_traits::{FromPrimitive,ToPrimitive}`, `T::from_*`, `.to_f64()`,
  `T::zero()`, or `T::one()` in the touched defaults, validation bounds,
  residuals, adaptive-step constants, Rhie-Chow face caches, pressure
  extrapolation, boundary zeroing, PIMPLE unity relaxation, and reference-data
  interpolation paths. Production constants route through
  `crates/cfd-2d/src/scalar.rs`; validation helpers route through Eunomia.
  Evidence: `cargo fmt -p cfd-2d --check`, `cargo check -p cfd-2d
  --no-default-features`, `cargo clippy -p cfd-2d --no-default-features --lib
  -- -D warnings`, focused `cargo nextest run -p cfd-2d
  --no-default-features simplec_pimple --status-level fail` (4/4 passed), and
  clean direct-provider scans over the touched module/test. Residual work:
  direct `num-traits` still remains across broader cfd-2d network,
  turbulence, LBM/FDM/FVM/NS-FVM, solver, physics, and other test surfaces, so
  the cfd-2d manifest still owns `num-traits`.
- [x] `cfd-2d` [patch]: Migrate the PISO scalar-provider surface to Eunomia.
  `piso_algorithm::{config,convergence,predictor,corrector,solver}` no longer
  uses direct `num_traits::{FromPrimitive,ToPrimitive}`, `T::from_*`,
  `.to_f64()`, `T::zero()`, or `T::one()` in the touched defaults, tolerances,
  residual monitors, predictor/corrector workspaces, hybrid differencing,
  Rhie-Chow face-flux update, duration threshold, and logging paths. Constants
  route through the crate-local Eunomia scalar adapter, and abs/max/sqrt/f64
  diagnostics use `eunomia::{NumericElement,FloatElement}`. Evidence: `cargo
  fmt -p cfd-2d --check`, `cargo check -p cfd-2d --no-default-features`,
  `cargo clippy -p cfd-2d --no-default-features --lib -- -D warnings`,
  focused `cargo nextest run -p cfd-2d --no-default-features piso
  --status-level fail` (6/6 passed), and a clean direct-provider scan over
  `crates/cfd-2d/src/piso_algorithm`. Residual work: direct `num-traits` still
  remains across broader cfd-2d network, turbulence, LBM/FDM/FVM/NS-FVM,
  solver, physics, and test surfaces, so the cfd-2d manifest still owns
  `num-traits`.
- [x] `cfd-2d` [patch]: Migrate the pressure-velocity scalar-provider surface
  to Eunomia. `pressure_velocity::{coefficients,config,solver,rhie_chow,faces,
  correction,pressure}` no longer uses direct `T::zero()`/`T::one()` scalar
  construction in the touched coefficient defaults, validation bounds,
  pressure-correction assembly, scatter, solver workspace, and Rhie-Chow
  coefficient paths. Identities route through the crate-local Eunomia scalar
  adapter, and finite checks use `eunomia::NumericElement` in the touched
  paths. Evidence: `cargo fmt -p cfd-2d --check`, `cargo check -p cfd-2d
  --no-default-features`, `cargo clippy -p cfd-2d --no-default-features --lib
  -- -D warnings`, focused `cargo nextest run -p cfd-2d
  --no-default-features pressure_velocity --status-level fail` (16/16
  passed), and a clean direct-provider scan over
  `crates/cfd-2d/src/pressure_velocity`. Residual work: direct `num-traits`
  still remains across broader cfd-2d network, turbulence, LBM/FDM/FVM/
  NS-FVM, solver, physics, and test surfaces, so the cfd-2d manifest still
  owns `num-traits`.
- [x] `cfd-2d` [patch]: Migrate the time-integration scalar-provider surface
  to Eunomia. `schemes::time::{explicit,implicit,multistep,adaptive}` no
  longer uses direct `num_traits::{FromPrimitive,ToPrimitive}` bounds,
  `T::from_f64`, `T::zero()`, `T::one()`, or `unwrap_or(T::one())` in the
  touched RK, implicit fixed-point, BDF/Adams-Bashforth, CFL, Richardson-error,
  and adaptive-step paths. Constants route through the crate-local Eunomia
  scalar adapter, and abs/min/max/powf use explicit
  `eunomia::{NumericElement,FloatElement}` dispatch where needed. Evidence:
  `cargo fmt -p cfd-2d --check`, `cargo check -p cfd-2d
  --no-default-features`, `cargo clippy -p cfd-2d --no-default-features --lib
  -- -D warnings`, focused `cargo nextest run -p cfd-2d
  --no-default-features time --status-level fail` (29/29 passed), and a clean
  direct-provider scan over `crates/cfd-2d/src/schemes/time`. Residual work:
  direct `num-traits` still remains across broader cfd-2d SIMPLEC/PIMPLE,
  network, turbulence, LBM/FDM/FVM/NS-FVM, solver, physics, and test surfaces,
  so the cfd-2d manifest still owns `num-traits`.
- [x] `cfd-2d` [patch]: Migrate the finite-volume discretization
  scalar-provider surface to Eunomia. `discretization::{convection,
  extended_stencil}` and
  `physics::momentum::coefficient_corrections::quick` no longer use direct
  `num_traits::FromPrimitive` scalar construction or `T::zero()`/`T::one()`
  constructors for first-order upwind, central, hybrid, power-law, QUICK, and
  MUSCL extended-stencil paths. Constants, zero/one values, absolute values,
  and max operations route through the crate-local Eunomia scalar adapter.
  Evidence: `cargo fmt -p cfd-2d --check`, `cargo check -p cfd-2d
  --no-default-features`, `cargo clippy -p cfd-2d --no-default-features --lib
  -- -D warnings`, focused `cargo nextest run -p cfd-2d --no-default-features
  discretization --status-level fail` (8/8 passed), and a clean direct-provider
  scan over the touched files. Residual work: direct `num-traits` still
  remains across broader cfd-2d network, PISO, LBM, turbulence, NS-FVM, scalar
  transport, solver, and test surfaces, so the cfd-2d manifest still owns
  `num-traits`.
- [x] `cfd-2d` [patch]: Migrate the spatial/WENO scalar-provider surface to
  Eunomia. `schemes::{constants,central,grid,upwind,weno_helpers,weno,weno_z}`,
  `schemes::tvd::{mod,quick}`, and
  `physics::momentum::{coefficients,coefficient_corrections::weno_z}` no
  longer use direct `num_traits::{FromPrimitive,ToPrimitive}` scalar
  construction for the touched central/upwind/TVD/QUICK/WENO5/WENO9 and
  WENO-Z momentum correction paths. Constants, zero values, smoothness powers,
  weights, deferred-correction relaxation factors, and limiter constants route
  through the crate-local Eunomia scalar adapter. Evidence: `cargo fmt -p
  cfd-2d --check`, `cargo check -p cfd-2d --no-default-features`, `cargo
  clippy -p cfd-2d --no-default-features --lib -- -D warnings`, focused
  `cargo nextest run -p cfd-2d --no-default-features weno --status-level
  fail` (5/5 passed), and a clean direct-provider scan over the touched files.
  Residual work: direct `num-traits` still remains across broader cfd-2d
  discretization, network, PISO, LBM, turbulence, NS-FVM, solver, and test
  surfaces, so the cfd-2d manifest still owns `num-traits`.
- [x] `cfd-2d` [patch]: Migrate the field-container scalar-provider surface
  to Eunomia. `fields::{constants,Field2D,SimulationFields}` no longer imports
  direct `num_traits` or uses `FromPrimitive`, `Float`, `T::from_*`,
  `T::zero()`, or `T::one()` for field defaults, reset values, velocity
  magnitude, or Reynolds-number averaging; those now route through the
  crate-local Eunomia scalar adapter. Evidence: `cargo fmt -p cfd-2d
  --check`, `cargo check -p cfd-2d --no-default-features`, `cargo clippy -p
  cfd-2d --no-default-features --lib -- -D warnings`, focused `cargo nextest
  run -p cfd-2d --no-default-features piso_algorithm --status-level fail`
  (6/6 passed), and a clean direct-provider scan over
  `crates/cfd-2d/src/fields.rs`. Residual work: direct `num-traits` still
  remains across broader cfd-2d schemes, network, PISO, LBM, turbulence, and
  test surfaces, so the cfd-2d manifest still owns `num-traits`.
- [x] `cfd-3d` [patch]: Remove direct `num-traits` ownership from the 3D
  package. The final root `lib.rs` Chebyshev assertions now use inherent
  `f64::abs`, and the manifest no longer declares `num-traits`. Evidence:
  `cargo fmt -p cfd-3d --check`, `cargo check -p cfd-3d
  --no-default-features`, focused `cargo nextest run -p cfd-3d
  --no-default-features chebyshev --status-level fail` (20/20 passed),
  `cargo clippy -p cfd-3d --no-default-features --lib -- -D warnings`, and a
  direct-provider scan across `crates/cfd-3d/{src,tests,examples}` plus
  `crates/cfd-3d/Cargo.toml` with no `num_traits`, `num-traits`,
  `FromPrimitive`, `ToPrimitive`, `num_traits::Float`, `Float::`,
  `from_f64_or_one`, `<T as From<f64>>::from`, `T::zero()`, `T::one()`, or
  `T::from_` residue. Residual work: `cargo tree -p cfd-3d -e normal -i
  num-traits` still shows transitive provider/nalgebra paths through
  `approx`, `nalgebra`, `gaia`, `half`, `num-complex`, and related upstream
  dependencies.
- [x] `cfd-3d` [patch]: Migrate the Bifurcation and Venturi scalar-provider
  seams to Eunomia. `bifurcation::{analysis,geometry,solver,types,
  validation}` and `venturi::{analysis,solver,types,validation}` no longer
  use direct `num_traits::{Float,FromPrimitive,ToPrimitive}` or
  `num_traits::Float::*`; scalar constants, zero/one values, elementary math,
  f64/usize construction, pressure coefficients, Richardson/GCI metrics,
  Picard viscosity-change checks, flow extraction, pressure slicing, and
  validation thresholds route through `cfd-3d::scalar`. Evidence: `cargo fmt
  -p cfd-3d --check`, `cargo check -p cfd-3d --no-default-features`,
  focused `cargo nextest run -p cfd-3d --no-default-features bifurcation
  venturi --status-level fail` (35/35 passed; one 16.7s nextest-slow test),
  `cargo clippy -p cfd-3d --no-default-features --lib -- -D warnings`, and a
  targeted Bifurcation/Venturi residue scan. Residual work: direct cfd-3d
  `num-traits` ownership is now closed; nalgebra geometry/vector/matrix
  replacement remains a later Leto/Gaia pass.
- [x] `cfd-3d` [patch]: Migrate the Trifurcation scalar-provider seam to
  Eunomia. `trifurcation::{geometry,solver,validation}` no longer use direct
  `num_traits::{Float,FromPrimitive,ToPrimitive}` or `num_traits::Float::*`;
  scalar constants, zero/one values, abs/min/max/sqrt, sin/cos, integer
  powers, f64 conversion, solver defaults, pressure-drop metrics, and
  validation thresholds route through `cfd-3d::scalar`. Evidence: `cargo fmt
  -p cfd-3d`, `cargo check -p cfd-3d --no-default-features`, focused `cargo
  nextest run -p cfd-3d --no-default-features trifurcation --status-level
  fail` (6/6 passed), `cargo clippy -p cfd-3d --no-default-features --lib --
  -D warnings`, and a targeted Trifurcation residue scan. Residual work:
  direct cfd-3d `num-traits` ownership is now closed; nalgebra
  geometry/vector replacement remains a later Leto/Gaia pass.
- [x] `cfd-3d` [patch]: Migrate the Serpentine scalar-provider seam to
  Eunomia. `serpentine::{solver,validation}` no longer use direct
  `num_traits::{Float,FromPrimitive,ToPrimitive}` or `num_traits::Float::*`;
  scalar constants, zero/one values, abs/max/sqrt, integer-to-scalar
  construction, and f64 diagnostics route through `cfd-3d::scalar`. Evidence:
  `cargo fmt -p cfd-3d`, `cargo check -p cfd-3d --no-default-features`,
  focused `cargo nextest run -p cfd-3d --no-default-features serpentine
  --status-level fail` (7/7 passed), `cargo clippy -p cfd-3d
  --no-default-features --lib -- -D warnings`, and a targeted Serpentine
  residue scan. Residual work: direct cfd-3d `num-traits` ownership is now
  closed; nalgebra geometry/vector
  replacement remains a later Leto/Gaia pass.
- [x] `cfd-3d` [patch]: Migrate the IBM scalar-provider seam to Eunomia.
  `ibm::{forcing,interpolation,solver}` no longer use direct
  `num_traits::{FromPrimitive,ToPrimitive}` or `num_traits::Float::*`; scalar
  constants, zero/one values, abs/max/sqrt/cos/floor, and index scalar
  construction route through `cfd-3d::scalar`. Grid-index conversion now
  returns typed configuration errors for non-finite/out-of-range scaled
  coordinates instead of silently using zero. Evidence: `cargo fmt -p cfd-3d`,
  `cargo check -p cfd-3d --no-default-features`, focused `cargo nextest run
  -p cfd-3d --no-default-features ibm --status-level fail` (12/12 passed),
  `cargo clippy -p cfd-3d --no-default-features --lib -- -D warnings`, and a
  targeted IBM residue scan. Residual work: direct cfd-3d `num-traits`
  ownership is now closed;
  nalgebra geometry/vector replacement remains a later Leto/Gaia pass.
- [x] `cfd-3d` [patch]: Migrate the spectral direct scalar-provider residue to
  Eunomia. `spectral::{chebyshev,poisson,solver}` and the co-located
  Chebyshev tests no longer use direct `num_traits::{FromPrimitive,Float}` or
  `num_traits::Float::*` calls; scalar constants, signs, zero/one values, and
  generic absolute values route through `cfd-3d::scalar`. Evidence: `cargo fmt
  -p cfd-3d`, `cargo check -p cfd-3d --no-default-features`, focused `cargo
  nextest run -p cfd-3d --no-default-features spectral --status-level fail`
  (41/41 passed), `cargo clippy -p cfd-3d --no-default-features --lib -- -D
  warnings`, and a targeted spectral residue scan. Residual work: spectral
  still owns nalgebra `DMatrix`/`DVector` and `RealField` surfaces pending the
  Leto matrix/storage migration; direct cfd-3d `num-traits` ownership is now
  closed.
- [x] `cfd-3d` [patch]: Migrate the level-set scalar-provider seam to
  Eunomia. `level_set::{weno,advection,solver}` now route scalar constants,
  finite checks, min/max/abs/sqrt, and integer powers through the crate-local
  `cfd-3d::scalar` adapter backed by `eunomia::{FloatElement,NumericElement}`
  instead of direct `num_traits::{FromPrimitive,Float}` calls. Evidence:
  `cargo fmt -p cfd-3d --check`, `cargo check -p cfd-3d
  --no-default-features`, focused `cargo nextest run -p cfd-3d
  --no-default-features level_set --status-level fail` (13/13 passed),
  `cargo clippy -p cfd-3d --no-default-features --lib -- -D warnings`, and a
  targeted level-set residue scan. Residual work: direct cfd-3d
  `num-traits` ownership is now closed; nalgebra/vector replacement remains a
  later Leto/Gaia pass.
- [x] `cfd-validation` [patch]: Remove direct `num-traits` ownership from the
  validation package. The last direct residues in
  `benchmarks::{cavity,step,poiseuille_bifurcation}` and
  `tests/complex_boundary_mms_validation.rs` now route scalar construction,
  diagnostic conversion, absolute values, powers, and max/min dispatch through
  Eunomia `FloatElement` or `cfd-validation::scalar`; the package manifest no
  longer declares `num-traits`. Evidence: `cargo fmt -p cfd-validation
  --check`, `cargo check -p cfd-validation`, full `cargo nextest run -p
  cfd-validation --status-level fail` (431/431 passed), focused direct-provider
  residue scan, and `cargo tree -p cfd-validation -e normal -i num-traits`
  showing only transitive provider/nalgebra paths. Residual validation
  provider work remains in nalgebra/nalgebra-sparse storage and Leto/Gaia
  mesh/vector ownership, not direct `num-traits` ownership.
- [x] `cfd-2d`/`cfd-validation` [patch]: Migrate the pressure-velocity
  SIMPLEC/PIMPLE scalar-provider boundary and dependent validation 2D blood
  benchmark wrappers to Eunomia. The touched cfd-2d solver cone no longer
  exposes direct `num_traits::{FromPrimitive,ToPrimitive}` bounds, and
  `cfd-validation::benchmarks::{bifurcation,venturi,trifurcation,serpentine}`
  now use `cfd-validation::scalar` for touched scalar constants, diagnostics,
  powers, square roots, and zero values instead of direct `num_traits`
  contracts. Evidence: `cargo fmt -p cfd-2d -p cfd-validation --check`,
  `cargo check -p cfd-validation`, focused `cargo nextest run -p cfd-2d
  pressure_velocity simplec_pimple --status-level fail` (20/20 passed), full
  `cargo nextest run -p cfd-2d --status-level fail` (567/567 passed, 27
  skipped), full `cargo nextest run -p cfd-validation --status-level fail`
  (431/431 passed), and targeted residue scans. Residual work remains in
  heavier cfd-2d momentum helpers outside the core solver/solve impls and in
  validation benchmark bodies such as cavity, step, and disabled
  poiseuille-bifurcation.
- [x] `cfd-validation` [patch]: Migrate the 3D benchmark scalar-provider
  wrapper seam to Eunomia. `BenchmarkConfig::default` and
  `benchmarks::threed::{bifurcation,venturi,serpentine}` now route touched
  scalar constants and absolute-value validation predicates through
  `cfd-validation::scalar`; the wrappers no longer name direct
  `num_traits::{Float,FromPrimitive,ToPrimitive}` or
  `cfd_core::conversion::SafeFromF64`. Evidence: `cargo fmt -p
  cfd-validation --check`, `cargo check -p cfd-validation`, focused `cargo
  nextest run -p cfd-validation bifurcation_flow_3d serpentine_flow_3d
  venturi_flow_3d` (3/3 passed), full `cargo nextest run -p cfd-validation
  --status-level fail` (431/431 passed), and a targeted benchmark residue
  scan. Residual work: the validation wrappers still require
  `std::convert::From<f64>` because the `cfd-3d` bifurcation, serpentine, and
  Venturi solver constructors require it; broader benchmark storage/vector
  migration to Leto/Gaia remains open.
- [x] `cfd-validation` [patch]: Migrate literature blood-flow validation and
  numerical linear-solver scalar thresholds to Eunomia. `blood_flow_1d` now
  routes scalar construction, f64 diagnostics, absolute values, and
  Poiseuille resistance `powi` through `cfd-validation::scalar` rather than
  direct `num_traits::{FromPrimitive,ToPrimitive}` / `num_traits::Float`.
  `numerical::linear_solver` now uses `FloatElement` bounds and crate-local
  scalar constants for tolerances and convergence-rate metadata. Evidence:
  `cargo fmt -p cfd-validation`, `cargo fmt -p cfd-validation --check`,
  `cargo check -p cfd-validation`, focused `cargo nextest run -p
  cfd-validation blood_flow literature chapman patankar` (6/6 passed), full
  `cargo nextest run -p cfd-validation` (431/431 passed), and a combined
  targeted residue scan over the literature module, `numerical/linear_solver.rs`,
  and `scalar.rs`. Residual `cfd-validation` provider migration remains in
  benchmark/test code plus nalgebra storage/vector seams pending Leto/Gaia.
- [x] `cfd-validation` [patch]: Migrate validation geometry scalar math to
  Eunomia. `geometry::{serpentine_2d,venturi}` and
  `geometry::threed::serpentine` no longer use direct `T::from_f64`,
  `T::from_usize`, `num_traits::Zero`, or RealField elementary dispatch in
  their touched scalar constants, trigonometric calls, square roots, powers,
  absolute values, and max selection; these route through the crate-local
  Eunomia scalar adapter with `FloatElement` bounds. Evidence: `cargo fmt -p
  cfd-validation`, `cargo check -p cfd-validation`, focused `cargo nextest
  run -p cfd-validation geometry` (11/11 passed), full `cargo nextest run -p
  cfd-validation` (431/431 passed), and a geometry-wide targeted provider
  residue scan. Remaining geometry work is Leto/Gaia storage and mesh
  ownership, not scalar-provider construction.
- [x] `cfd-validation` [patch]: Migrate the manufactured Richardson
  extrapolation stack to Eunomia. `manufactured::richardson::{core,
  validation,analysis}` no longer imports direct `num_traits` conversion/float
  traits, `nalgebra::ComplexField`, or test-only `nalgebra::scalar` helpers;
  scalar constants, powers, logarithms, square roots, absolute values, finite
  checks, and f64 diagnostics now route through the crate-local Eunomia scalar
  adapter. Evidence: `cargo fmt -p cfd-validation`, `cargo check -p
  cfd-validation`, focused `cargo nextest run -p cfd-validation richardson`
  (16/16 passed), full `cargo nextest run -p cfd-validation` (431/431
  passed), and a manufactured-wide residue scan with no targeted direct
  provider hits. The `nalgebra::RealField` generic boundary remains for later
  Leto replacement.
- [x] `cfd-validation` [patch]: Migrate the manufactured scalar/trig MMS
  cluster to Eunomia. `manufactured::{advection,advection_diffusion,burgers,
  diffusion,navier_stokes}` no longer imports direct
  `num_traits::FromPrimitive`, cfd-core `SafeFromF64`, or
  `nalgebra::ComplexField` for scalar construction and transcendental math;
  constants, `sin`/`cos`, `exp`, and `sqrt` route through the crate-local
  Eunomia scalar adapter. Evidence: `cargo fmt -p cfd-validation --check`,
  `cargo check -p cfd-validation`, focused `cargo nextest run -p
  cfd-validation manufactured` (50/50 passed), full `cargo nextest run -p
  cfd-validation` (431/431 passed), and focused touched-file residue scans.
  The `nalgebra::RealField`/`Vector2` storage boundary remains for later Leto
  replacement; Richardson manufactured direct-provider work is closed by the
  later Richardson slice.
- [x] `cfd-validation` [patch]: Migrate the Chapman-Enskog literature scalar
  seam to Eunomia. `literature::chapman_enskog` no longer imports direct
  `num_traits::{FromPrimitive,ToPrimitive}` or uses `T::from_f64`/generic
  `.to_f64()` conversion; transport coefficient and report scalar values now
  route through the crate-local Eunomia adapter with explicit `FloatElement`
  bounds. Evidence: `cargo fmt -p cfd-validation --check`, focused `cargo
  nextest run -p cfd-validation chapman` (2/2 passed), full `cargo nextest run
  -p cfd-validation` (431/431 passed), and a focused Chapman-Enskog residue
  scan. Remaining literature/provider work is scoped to `blood_flow_1d.rs`,
  the Patankar storage boundary, and broader validation numerical/storage
  seams.
- [x] `cfd-validation` [patch]: Migrate the `time_integration` conversion seam
  to Eunomia. `integrators` and `stability_analysis` no longer use direct
  `num_traits::{FromPrimitive,ToPrimitive}` or cfd-core `SafeFromF64`
  conversion helpers; Runge-Kutta constants, CFL/von-Neumann scalar inputs,
  wave-number interpolation factors, and stability-report f64 diagnostics now
  route through the crate-local Eunomia scalar adapter. The local
  von-Neumann spatial operators use `eunomia::Complex` instead of
  `num_complex::Complex`. Evidence: `cargo fmt -p cfd-validation`, `cargo
  check -p cfd-validation`, focused `cargo nextest run -p cfd-validation
  time_integration` (12/12 passed), full `cargo nextest run -p
  cfd-validation` (429/429 passed), and focused time-integration residue
  scans. The current nalgebra `DVector`/`DMatrix` storage boundary remains
  tracked for later Leto replacement.
- [x] `cfd-validation` [patch]: Migrate the `error_metrics` scalar metrics
  group to Eunomia. `norms`, `normalized`, and `statistics` no longer depend
  on direct `num_traits::FromPrimitive` or `T::from_f64`/`T::from_usize`;
  scalar construction, absolute values, and square roots now use the
  crate-local Eunomia adapter and explicit `FloatElement` bounds. Evidence:
  `cargo fmt -p cfd-validation`, `cargo check -p cfd-validation`, focused
  `cargo nextest run -p cfd-validation error_metrics` (21/21 passed), full
  `cargo nextest run -p cfd-validation` (429/429 passed), and focused
  error-metrics residue scans.
- [x] `cfd-validation` [patch]: Close the current Eunomia scalar
  compile-blocker group. `PatankarLidDrivenCavity` no longer uses direct
  `num_traits`/`SafeFromF64` conversion for reference constants or grid
  lookup, the advanced and multi-physics manufactured-solution paths now route
  scalar constants and transcendental calls through the crate-local Eunomia
  adapter, Reynolds-stress MMS no longer depends on direct `FromPrimitive` or
  `nalgebra::ComplexField`, and Richardson analysis now carries the
  `FloatElement` bound required by its result API. Evidence: `cargo fmt -p
  cfd-validation`, `cargo check -p cfd-validation`, `cargo nextest run -p
  cfd-validation` (429/429 passed), and focused migrated-file residue scans.
  Residual direct provider work remains in other `cfd-validation` modules,
  including `numerical/linear_solver`, benchmarks, literature validations, and
  remaining manufactured/Richardson files.
- [x] `cfd-1d` [patch]: Remove direct package `num-traits` ownership and
  complete the active resistance/vascular Eunomia scalar seam. Resistance
  model traits now route scalar construction and diagnostics through
  Eunomia/cfd-core provider traits; hydraulic resistance geometry/models,
  serpentine/slug-flow, Bessel/Womersley, structured-tree, and bifurcation
  math use explicit `FloatElement`/`NumericElement` dispatch where nalgebra
  remains the current storage boundary. The remaining direct cfd-1d
  `num_traits` residues in network blueprint/sink bounds, solver-analysis
  averages, and package tests were removed, and `crates/cfd-1d/Cargo.toml`
  no longer declares `num-traits`. Evidence: `cargo fmt -p cfd-1d --check`,
  `cargo check -p cfd-1d`, `cargo nextest run -p cfd-1d` (725/725 passed, 3
  skipped), full cfd-1d direct-provider residue scan, and
  `cargo tree -p cfd-1d -e normal -i num-traits` showing only transitive
  provider/nalgebra graph paths. Broad all-target clippy still fails on
  existing unrelated lint debt outside this cleanup.
- [x] `cfd-1d` [patch]: Migrate the solver-core scalar contract and linear
  solve helpers to Eunomia. `NetworkSolveScalar` no longer carries direct
  `FromPrimitive`, `ToPrimitive`, or `num_traits::Float` compatibility bounds;
  solver-core Anderson acceleration, convergence checks, linear-system
  equilibration/Jacobi preconditioning, SPD detection, residual norms, and f64
  diagnostics now route through `FloatElement`, `NumericElement`, and
  cfd-core conversion traits. Evidence: `cargo fmt -p cfd-1d --check`,
  `cargo check -p cfd-1d`, `cargo nextest run -p cfd-1d` (725/725 passed, 3
  skipped), plus a focused solver-core residue scan. Residual direct provider
  work remains in vascular Bessel/Womersley, resistance scalar traits,
  tests/benches, and the current nalgebra/nalgebra-sparse storage boundary.
- [x] `cfd-1d` [patch]: Migrate the network wrapper scalar/provider seam to
  Eunomia. `EdgeProperties::from`, characteristic-length calculation, Picard
  resistance refresh, blood hematocrit propagation, bifurcation
  phase-separation bridges, coefficient validation, and parallel edge
  conductance now route scalar construction, f64 bridge conversion, finite
  checks, and absolute values through `SafeFromF64`/`SafeFromUsize` and
  Eunomia `NumericElement`. Adjacent solver bounds were aligned so
  `NetworkProblem`/`NetworkSolver` use `NetworkSolveScalar`, and
  `MatrixAssembler` no longer requires `FromPrimitive`. Evidence: `cargo fmt
  -p cfd-1d --check`, `cargo check -p cfd-1d`, `cargo nextest run -p cfd-1d`
  (725/725 passed, 3 skipped), plus a focused wrapper/matrix/problem residue
  scan. Residual work remains in `solver/core/mod.rs`'s direct
  `NetworkSolveScalar` num-traits compatibility contract and broader
  solver-core/vascular/tests/storage seams.
- [x] `cfd-1d` [patch]: Migrate the network blueprint/sink scalar construction
  seam to Atlas provider conversions. The canonical `network_from_blueprint`
  path now uses `SafeFromF64`/`SafeFromUsize` for blueprint scalar constants,
  areas, boundary values, serpentine segment counts, blood defaults, and
  initial flow probes, and uses Eunomia `NumericElement::abs` for the generic
  quadratic-coefficient policy check. `NetworkBuilderSink` propagates the
  required Atlas conversion bounds. Evidence: `cargo fmt -p cfd-1d`, `cargo
  check -p cfd-1d`, and `cargo nextest run -p cfd-1d` (725/725 passed, 3
  skipped). Residual direct `FromPrimitive` in this path is inherited from
  the current `Network::new`/`wrapper.rs` impl and is tracked as the next
  network cleanup seam.
- [x] `cfd-1d` [patch]: Migrate the domain-components provider seam from
  direct `num_traits` bounds to Eunomia. `Component<T>` pressure-drop math,
  component factory defaults, channels, membranes, mixers, pumps, valves, and
  sensors now use `SafeFromF64`, `SafeFromUsize`, and `NumericElement` instead
  of local `FromPrimitive`/`Float` requirements. Evidence: `cargo fmt -p
  cfd-1d --check`, `cargo check -p cfd-1d`, `cargo nextest run -p cfd-1d`
  (725/725 passed, 3 skipped), a focused component residue scan, and a clippy
  rerun confirming the touched `channels.rs` lint was removed. Broad
  all-target clippy still fails on unrelated examples/tests and other
  cfd-1d modules.
- [x] `cfd-1d` [patch]: Migrate the next channel, branching, and
  solver-analysis provider seam from direct `num_traits` bounds to Eunomia.
  `FlowRegime`, channel flow resistance/geometry/shape factors, branching
  network solver bounds, and network-analysis pressure/flow/resistance/
  performance paths now use `SafeFromF64`, `FloatElement`, and
  `NumericElement`. Evidence: `cargo fmt -p cfd-1d`, `cargo check -p cfd-1d`,
  `cargo nextest run -p cfd-1d` (725/725 passed, 3 skipped), plus focused
  residue scans over the touched channel/branching/analyzer and
  `solver/analysis` cones. `cargo clippy -p cfd-1d --all-targets -- -D
  warnings` is still blocked by unrelated pre-existing all-target lint debt.
- [x] `cfd-1d` [patch]: Migrate the Murray's-law vascular geometry provider
  boundary from `num_traits::FromPrimitive` to Eunomia. `MurraysLaw` and
  `OptimalBifurcation` now use `FloatElement`/`NumericElement` for scalar
  construction, power operations, absolute value, and inverse cosine. The
  vascular bifurcation network builder carries the Eunomia bound required by
  its Murray's-law calls. Evidence: `cargo fmt -p cfd-1d --check`, a focused
  residue scan over the Murray's-law files, and upstream `cargo nextest run -p
  eunomia acos` (2/2). `cargo check -p cfd-1d` is still blocked by unrelated
  dirty-tree errors outside this slice, so nextest coverage for this package
  remains pending behind that cleanup.
- [x] `cfd-math` [patch]: Remove direct `num-traits` ownership from the
  remaining package residue. GPU scalar extraction now uses Eunomia
  `NumericElement`, AMG integration-test constants use Eunomia `FloatElement`,
  and the package manifest no longer declares `num-traits`. The touched GPU
  path also propagates Hephaestus synchronization failures instead of ignoring
  the result and no longer drives synchronous GPU work through a stale async
  bridge. Evidence: `cargo fmt -p cfd-math --check`, `cargo check -p
  cfd-math`, `cargo check -p cfd-math --features gpu`, both default and
  GPU-feature clippy gates with warnings denied, `cargo nextest run -p
  cfd-math` (333/333), and a focused residue scan. Remaining `num-traits`
  graph entries are transitive through provider/nalgebra stacks.
- [x] `cfd-core` [patch]: Remove direct `num-traits` ownership from the
  touched scalar conversion boundary. `management::conversion` and mesh
  centroid operations now use Eunomia `FloatElement` for scalar construction,
  and `cfd-core` no longer declares `num-traits` directly. Evidence: `cargo
  fmt -p cfd-core --check`, `cargo check -p cfd-core`, `cargo clippy -p
  cfd-core --all-targets -- -D warnings`, `cargo nextest run -p cfd-core`
  (229/229), and a focused residue scan. Remaining `num-traits` graph entries
  are transitive through provider/nalgebra stacks until the later replacement
  slices land.
- [x] `workspace` [patch]: Remove stale direct workspace provider declarations
  for `wgpu 0.19` and `num-complex`. The root manifest no longer exposes an
  unused direct WGPU provider after the `cfd-core` Hephaestus handoff, and it
  no longer advertises a direct num-complex provider now that migrated complex
  call sites use Eunomia/Apollo vocabulary. Evidence: `cargo metadata
  --no-deps --format-version 1`, `cargo check -p cfd-suite`, focused manifest
  scans, and `cargo tree --workspace -e normal -i wgpu@0.19.4` returning no
  matching package. Residual `num-complex` is transitive through provider and
  nalgebra stacks until the remaining Leto/Eunomia migration slices land.
- [x] `cfd-core`/`hephaestus-wgpu` [minor]: Remove the direct `cfd-core`
  `wgpu 0.19` dependency and route the crate-level GPU ABI through
  Hephaestus. `GpuContext` now acquires the provider-owned `WgpuDevice`, GPU
  modules import WGPU descriptor/buffer/pipeline types through
  `hephaestus_wgpu::wgpu`, and the GPU feature depends only on
  `hephaestus-wgpu`. Evidence: `cargo fmt -p cfd-core --check`, `cargo check
  -p cfd-core`, `cargo nextest run -p cfd-core` (228/228), and dependency-tree
  audits proving `wgpu@0.19.4` is absent while `wgpu@26.0.1` reaches
  `cfd-core` only through Hephaestus. `cargo clippy -p cfd-core --all-targets
  -- -D warnings` is also clean after resolving the package-gate lints exposed
  by this migration.
- [x] `cfd-optim` [patch]: Remove the unused direct nalgebra dev-dependency
  and keep the optimizer package boundary free of direct Atlas-replaced
  providers. The cleanup also resolves package all-target clippy debt by using
  `std::slice::from_ref`, moving test modules after items, and promoting
  computed SDT acoustic energy/contrast values into `SdtMetrics` instead of
  retaining unused underscore-prefixed fields. Verified with the focused
  cfd-optim package gate, source/provider-residue scans, and a dependency-tree
  audit showing remaining nalgebra only through upstream transitive crates.
- [x] `cfd-io` [patch]: Sync the crate-local `agents.md` architecture
  reference with the completed Atlas provider boundary: no direct internal
  solver/core/math dependencies, Leto dense checkpoint/binary arrays, Eunomia
  scalar bounds, Consus-backed optional HDF5, and optional RITK VTK re-export.
  Verified with the focused cfd-io package gate, default normal nalgebra tree
  audit, and provider-residue scans.
- [x] `cfd-io` [patch]: Close the residual checkpoint-validator scalar
  constructor residue. `CheckpointValidator::check_mass_conservation` now
  constructs spacing and central-difference constants through Eunomia
  `FloatElement`, leaving the crate source/test/manifest scan clean for direct
  `num-traits`, direct nalgebra/ndarray, and direct `T::from_*` constructor
  residue. Evidence: `cargo fmt -p cfd-io --check`, `cargo check -p cfd-io`,
  `cargo clippy -p cfd-io --all-targets -- -D warnings`, `cargo nextest run
  -p cfd-io` (3/3), `cargo test --doc -p cfd-io`, `cargo doc -p cfd-io
  --no-deps`, and focused cfd-io provider-residue scan.
- [x] `cfd-3d` [patch]: Replace FEM solution velocity/pressure storage from
  nalgebra `DVector` to a Leto-backed `FemDofVector`. The sparse solver and
  projection kernels now cross an explicit `DVector` conversion boundary only
  where the current linear-solver API still requires it; Anderson acceleration
  and Venturi Picard convergence operate on the provider-owned FEM DOF vector.
  Verified with `cargo fmt -p cfd-3d --check`, `cargo check -p cfd-3d`,
  `cargo clippy -p cfd-3d --lib -- -D warnings`, `cargo nextest run -p
  cfd-3d fem` (40/40), and focused scans. Residual FEM provider work remains
  in element `DMatrix`/`DVector` internals and solver/projection sparse
  linear-system `DVector` boundaries.
- [x] `cfd-3d` [patch]: Replace FEM solver/projection/solution scalar identity
  residue with the FEM-local Eunomia scalar SSOT. `FemSolver`,
  `ProjectionSolver`, and `StokesFlowSolution::blend` now use `fem::scalar`
  for direct zero/one construction while preserving the current nalgebra
  `DVector`/matrix storage boundary for the planned Leto FEM storage slice.
  Verified with `cargo fmt -p cfd-3d --check`, `cargo check -p cfd-3d`,
  `cargo clippy -p cfd-3d --lib -- -D warnings`, `cargo nextest run -p
  cfd-3d fem` (40/40), and a focused residue scan. Residual FEM provider work
  remains in nalgebra-backed FEM matrix/vector storage and solver/projection
  Leto storage migration.
- [x] `cfd-3d` [patch]: Replace FEM element scalar construction/math residue
  with the FEM-local Eunomia scalar SSOT while leaving nalgebra matrix storage
  for the planned Leto storage slice. `ElementMatrices` and `FluidElement` now
  use Eunomia `FloatElement`/`NumericElement` plus `fem::scalar` for zero
  initialization, volume absolute value, degeneracy tolerance, half factors,
  and sixth-volume constants. Verified with `cargo fmt -p cfd-3d --check`,
  `cargo check -p cfd-3d`, `cargo clippy -p cfd-3d --lib -- -D warnings`,
  `cargo nextest run -p cfd-3d element` (5/5), `cargo check -p
  moirai-executor --tests`, and a focused residue scan. Residual provider work
  remains in solver/projection assembly, solution scalar identities, and
  nalgebra-backed FEM matrix/vector storage.
- [x] `cfd-3d` [patch]: Replace FEM problem-validation scalar residue with the
  FEM-local Eunomia scalar SSOT. Physical invariant validation now uses
  Eunomia `FloatElement`/`NumericElement` for finite checks and positive-value
  guards, with `StokesFlowProblem::validate` carrying the explicit Eunomia
  scalar contract. Verified with `cargo fmt -p cfd-3d --check`, `cargo check
  -p cfd-3d`, `cargo clippy -p cfd-3d --lib -- -D warnings`, `cargo nextest
  run -p cfd-3d validate` (10/10), and a focused residue scan. Residual
  provider work remains in element assembly, solver/projection assembly,
  solution scalar identities, and nalgebra-backed FEM matrix/vector storage.
- [x] `cfd-3d` [patch]: Replace FEM mesh-helper scalar residue in P2
  extraction and mid-node lookup with the FEM-local Eunomia scalar SSOT.
  `mesh_utils` and `mid_node_cache` now use `fem::scalar` plus Eunomia
  `FloatElement`/`NumericElement` for midpoint constants, zero comparison,
  extrema initialization, and component min/max. Verified with `cargo fmt -p
  cfd-3d --check`, `cargo check -p cfd-3d`, `cargo clippy -p cfd-3d --lib --
  -D warnings`, `cargo nextest run -p cfd-3d mesh_utils` (5/5), `cargo
  nextest run -p cfd-3d mid_node_cache` (2/2), and a focused residue scan.
  Residual provider work remains in FEM problem validation, element assembly,
  solver/projection assembly, solution scalar identities, and nalgebra-backed
  FEM matrix/vector storage.
- [x] `cfd-3d` [patch]: Replace the next FEM scalar residue in stabilization
  and axial boundary classification with the FEM-local Eunomia scalar SSOT.
  `StabilizationParameters`, directional element-size helpers, and
  `AxialBoundaryClassifier` now use `fem::scalar` plus Eunomia
  `FloatElement`/`NumericElement` for constants, identities, extrema, absolute
  values, square roots, and scalar min/max. Verified with `cargo fmt -p
  cfd-3d --check`, `cargo check -p cfd-3d`, `cargo clippy -p cfd-3d --lib --
  -D warnings`, `cargo nextest run -p cfd-3d stabilization` (11/11), `cargo
  nextest run -p cfd-3d boundary` (12/12), and a focused residue scan.
  Residual provider work remains in `mesh_utils`, `mid_node_cache`,
  solver/projection assembly, and FEM matrix/vector storage.
- [x] `cfd-3d`/`cfd-2d`/`cfd-validation` [patch]: Replace the
  `cfd-3d::fem` constant-heavy scalar construction boundary with a FEM-local
  Eunomia scalar SSOT. `FemConfig`, Keast quadrature, stress/strain helpers,
  and P2 shape functions now construct constants and identities through
  `fem::scalar`/Eunomia instead of direct `num_traits::{FromPrimitive, Float}`,
  old `T::from_f64`, `T::zero`, `T::one`, or conversion fallbacks. Propagated
  the cfd-1d `NetworkSolveScalar` Eunomia `RealField` bound into the cfd-2d
  network reference and cfd-validation 1D blood-flow validation boundaries so
  cfd-3d dev/test builds compile against the current provider contract.
  Verified with `cargo check -p cfd-3d`, `cargo check -p cfd-validation`,
  `cargo clippy -p cfd-3d --lib -- -D warnings`, `cargo nextest run -p
  cfd-3d fem` (40/40), `cargo nextest run -p cfd-validation blood_flow`
  (2/2), formatting checks, and focused residue scans. `cargo clippy -p
  cfd-validation --lib -- -D warnings` was cancelled after prolonged shared
  build-lock contention. Residual provider work remains in FEM nalgebra
  matrix/vector storage, stabilization/mesh/boundary scalar residue, and the
  upstream cfd-1d `NetworkSolveScalar` num-traits compatibility contract.
- [x] `cfd-1d` [patch]: Replace transient composition/droplet direct
  `num_traits` scalar construction and `Float` dispatch with Eunomia
  `FloatElement`/`NumericElement`. Time-control tolerances, edge-transport
  Courant constants, segmented blood advection math, mixture comparisons, and
  droplet split-policy constants now route through Eunomia. Verified with
  `cargo fmt -p cfd-1d --check`, `cargo check -p cfd-1d`, `cargo clippy -p
  cfd-1d --lib -- -D warnings`, and transient composition/droplet/literature
  nextest binaries (20/20, 9/9, 5/5). Focused scans show no `ToPrimitive`,
  `FromPrimitive`, `num_traits::Float`, direct `num_traits`, `T::from_f64`,
  `T::from_usize`, direct `to_usize`, or silent substep-count fallback residue
  in the touched transient composition/droplet cone. Residual provider work
  remains in the inherited nalgebra `RealField` scalar boundary and the Pries
  phase-separation algorithm's f64 API contract.
- [x] `cfd-1d` [patch]: Replace the nonlinear solver's local Anderson
  least-squares implementation with the canonical Leto-backed `cfd-math`
  accelerator. The Picard loop now owns one `AndersonAccelerator` per solve,
  `SolverWorkspace` no longer stores duplicate Anderson residual/iterate
  queues, and `NetworkSolveScalar` documents the temporary bridge between the
  nalgebra linear-system boundary and Eunomia/Leto Anderson. Verified with
  `cargo fmt -p cfd-1d -p cfd-math --check`, `cargo check -p cfd-1d`, `cargo
  nextest run -p cfd-1d solver` (49/49), and focused residue scans. Residual
  provider work remains in cfd-1d nalgebra vector/sparse linear-system storage
  and transient composition direct `num_traits` conversions.
- [x] `cfd-math`/`cfd-2d`/`cfd-3d` [patch]: Replace the Anderson
  acceleration vector/matrix boundary from nalgebra to Leto. `cfd-math`
  Anderson QR state, history buffers, normal-equations solve, and
  `compute_next` now use Leto `Array1`/`Array2`; `cfd-2d::network::coupled`
  now mixes Leto resistance vectors; and four cfd-3d Picard solvers use a
  single `atlas_anderson` conversion boundary while FEM velocity storage still
  owns `DVector`. Verified with package formatting/checks, cfd-math lib
  clippy, `cargo nextest run -p cfd-math anderson` (5/5), `cargo nextest run
  -p cfd-2d network` (22/22), and focused residue scans. Residual provider
  work remains in cfd-3d FEM solution velocity storage, other cfd-math
  linear/time-stepping `DVector` APIs, validation benchmark/LES-DES helpers,
  and direct GPU replacement through Hephaestus.
- [x] `cfd-2d` [patch]: Replace turbulence constants-validation scalar
  helpers from nalgebra/num-traits to Eunomia. `TurbulenceConstantsValidator`,
  `ConstantsValidationResult`, and DNS sensitivity analysis now use Eunomia
  `RealField`/`FloatElement`/`NumericElement` for construction, constants,
  square roots, absolute values, maxima, and display conversion. Verified with
  `cargo fmt -p cfd-2d --check`, `cargo check -p cfd-2d`, `cargo nextest run
  -p cfd-2d constants` (11/11), `cargo nextest run -p cfd-2d rans` (5/5),
  `cargo nextest run -p cfd-2d macroscopic` (4/4), and a focused
  constants-validation residue scan. `cargo clippy -p cfd-2d --all-targets
  -- -D warnings` still fails on broader all-target lint debt outside this
  slice. Residual provider work remains in validation benchmark/LES-DES
  helpers, the cfd-2d network/cfd-math Anderson `DVector` seam, the adjacent
  `cfd-validation` RSM MMS scalar oracle, and direct GPU replacement through
  Hephaestus.
- [x] `cfd-2d` [patch]: Replace the shared RANS turbulence scalar/vector
  contract from nalgebra/num-traits to Leto/Eunomia. `TurbulenceModel<T>` now
  binds on Eunomia `RealField`; k-epsilon, realizable C_mu, k-omega SST,
  Spalart-Allmaras, wall boundary conditions, wall treatment, RANS validation,
  and SIMPLE turbulence coupling use Eunomia scalar construction/math with
  Leto `Vector2` velocity updates. Verified with `cargo fmt -p cfd-2d
  --check`, `cargo check -p cfd-2d`, `cargo nextest run -p cfd-2d k_epsilon`
  (39/39), `cargo nextest run -p cfd-2d k_omega` (19/19), `cargo nextest run
  -p cfd-2d spalart_allmaras` (21/21), `cargo nextest run -p cfd-2d wall`
  (51/51), `cargo nextest run -p cfd-2d rans` (5/5), and focused residue
  scans. `cargo clippy -p cfd-2d --all-targets -- -D warnings` is still
  pending because unrelated active builds held the shared `D:\atlas\target`
  lock. Residual turbulence provider work remains in validation scalar
  helpers, the adjacent `cfd-validation` RSM MMS scalar oracle, and direct GPU
  replacement through Hephaestus.
- [x] `cfd-2d` [patch]: Replace native Reynolds-stress scalar helper bounds
  from nalgebra `RealField` and `num_traits::FromPrimitive` to Eunomia
  `RealField`. RSM model construction, tensor storage, production, diffusion,
  curvature, wall-reflection, pressure-strain, and optimized transport math now
  use Eunomia scalar construction/constants/math while preserving the Leto
  `Array2` storage boundary. Verified with `cargo fmt -p cfd-2d`, `cargo
  check -p cfd-2d`, `cargo nextest run -p cfd-2d reynolds_stress` (10/10),
  and a focused residue scan. The shared RANS/Spalart-Allmaras contract is
  closed by the later Leto/Eunomia item; residual turbulence provider work is
  tracked in validation scalar helpers and Hephaestus GPU replacement.
- [x] `cfd-2d`/`cfd-validation` [patch]: Replace the Reynolds-stress matrix
  storage, transport, and MMS validation boundary from nalgebra `DMatrix` to
  Leto `Array2`. `ReynoldsStressTensor`, initialization/dissipation setup,
  production and diffusion helpers, transport velocity/scalar helpers,
  co-located tests, comprehensive Reynolds-stress validation fixtures, and the
  `cfd-validation` MMS L2-error oracle now use provider-owned arrays.
  Verified with `cargo fmt -p cfd-2d -p cfd-validation --check`, `cargo check
  -p cfd-2d`, `cargo nextest run -p cfd-2d reynolds_stress` (10/10), `cargo
  nextest run -p cfd-validation reynolds_stress` (9/9), and a focused residue
  scan. Residual cfd-2d turbulence provider work remains in validation scalar
  helpers, the adjacent `cfd-validation` RSM MMS scalar oracle, and direct GPU
  replacement through Hephaestus.
- [x] `cfd-2d` [patch]: Replace the WALE LES velocity-field and scalar-provider
  boundary from `Field2D<nalgebra::Vector2>` plus direct
  `num_traits::FromPrimitive` constants to Leto `Array2` component fields and
  Eunomia `RealField`/`NumericElement`. WALE SGS viscosity and boundary
  gradient recovery now consume provider-owned velocity component arrays and
  reject invalid component shapes, indices, and spacings through typed
  `cfd_core::error::Result`. Verified with `cargo fmt -p cfd-2d`, `cargo check
  -p cfd-2d`, `cargo nextest run -p cfd-2d wale` (4/4), and a focused residue
  scan. Residual cfd-2d turbulence provider work remains in validation scalar
  helpers, the adjacent `cfd-validation` RSM MMS scalar oracle, and direct GPU
  replacement through Hephaestus.
- [x] `cfd-2d` [patch]: Replace the shared LES/DES turbulence matrix boundary
  from nalgebra `DMatrix` to Leto `Array2`. The `LESTurbulenceModel` trait,
  Smagorinsky model state, strain/viscosity/dynamic helpers, GPU helper
  boundary, DES fields, DES length-scale utilities, and LES/DES validation and
  benchmark call sites now use provider-owned matrix storage. Verified with
  `cargo fmt -p cfd-2d --check`, `cargo check -p cfd-2d`, `cargo nextest run
  -p cfd-2d les_smagorinsky` (43/43), `cargo nextest run -p cfd-2d des`
  (14/14), and a focused residue scan. Residual turbulence provider work
  remains in validation scalar helpers, the adjacent `cfd-validation` RSM MMS
  scalar oracle, and direct GPU replacement through Hephaestus.
- [x] `cfd-2d` [patch]: Replace the standalone MILES LES 2x2 velocity-gradient
  API from nalgebra `DMatrix` plus direct `num_traits::FromPrimitive`
  constants to Leto `Array2` and Eunomia `RealField`/`NumericElement`/
  `FloatElement`. Config defaults, shock detection, numerical flux,
  applicability validation, and co-located value tests now use provider-owned
  scalar and matrix surfaces. Verified with `cargo fmt -p cfd-2d --check`,
  `cargo check -p cfd-2d`, `cargo nextest run -p cfd-2d miles` (7/7), and a
  focused residue scan. Residual turbulence provider work remains in
  validation scalar helpers and direct GPU replacement through Hephaestus.
- [x] `cfd-2d` [patch]: Replace the standalone Sigma/Vreman LES 2x2
  gradient and SGS-stress APIs from nalgebra `DMatrix` plus direct
  `num_traits::FromPrimitive` constants to Leto `Array2` and Eunomia
  `RealField`/`NumericElement`. Config defaults, viscosity, stress, invariant,
  and co-located value tests now use provider-owned scalar and matrix
  surfaces. Verified with `cargo fmt -p cfd-2d --check`, `cargo check -p
  cfd-2d`, `cargo nextest run -p cfd-2d sigma` (7/7), `cargo nextest run -p
  cfd-2d vreman` (6/6), and focused residue scans. Residual turbulence
  provider work remains in validation scalar helpers, the adjacent
  `cfd-validation` RSM MMS scalar oracle, and direct GPU replacement through
  Hephaestus.
- [x] `cfd-3d` [patch]: Replace the spectral Poisson/SpectralSolver public
  storage boundary from nalgebra `DVector` to Leto `Array1`. `PoissonProblem`,
  `SpectralSolution`, `PoissonSolver::solve`, the Poisson validation/property/
  domain/robustness tests, and the spectral Poisson example now use provider
  arrays at the public seam. Verified with `cargo fmt -p cfd-3d --check`,
  `cargo check -p cfd-3d`, `cargo nextest run -p cfd-3d --test
  poisson_validation` (7/7), `cargo nextest run -p cfd-3d spectral_poisson`
  (4/4), and `cargo nextest run -p cfd-3d
  test_poisson_zero_rhs_zero_solution` (1/1). Residual spectral work remains
  in the internal Poisson dense assembly/LU `DMatrix`/`DVector` bridge and the
  separate Chebyshev nalgebra vector/matrix surface.
- [x] `cfd-3d` [patch]: Replace the spectral Fourier public boundary from
  nalgebra `DVector`/`Complex` to Leto `Array1` and Eunomia `Complex`.
  Apollo remains the transform executor through the existing Atlas FFT helper,
  while Fourier validation and internal spectral smoke tests now use provider
  arrays/scalars directly instead of preserving the legacy nalgebra seam.
  Verified with `cargo fmt -p cfd-3d --check`, `cargo check -p cfd-3d`,
  `cargo nextest run -p cfd-3d --test fourier_validation` (12/12), and a
  focused scan showing no `DVector`/nalgebra-complex residue in the touched
  Fourier boundary. Residual spectral work remains in Chebyshev/Poisson solver
  paths that still expose nalgebra vectors and matrices.
- [x] `cfd-2d`/`cfd-validation`/`cfd-math`/`cfd-3d`/`cfd-optim`
  [patch]: Remove remaining direct Rayon/Tokio crate-source references and
  propagate the Eunomia scalar-provider contract through the exposed cfd-2d
  reference-trace and cfd-validation blood-flow validation boundaries. Stale
  Rayon documentation now names Moirai, `cfd-2d::network::reference` uses
  Eunomia `FloatElement`/`NumericElement` for reference-fluid constants,
  scale factors, resistance summaries, and f64 extraction, and
  `cfd-validation::literature::blood_flow_1d` uses Eunomia helpers for
  Carreau-Yasuda blood-flow constants and diagnostic conversions. Verified
  with package rustfmt/checks, `cargo nextest run -p cfd-2d reference_trace`
  (3/3), `cargo nextest run -p cfd-validation blood_flow` (2/2), `cargo
  nextest run -p cfd-optim` (121/121), and residue scans showing no
  Rayon/Tokio crate-source hits plus no old scalar construction/extraction
  residue in the two touched scalar-boundary files. Residual provider work
  remains the broader nalgebra/num-traits network, solver, FEM, and validation
  surfaces plus direct GPU/wgpu replacement with Hephaestus.
- [x] `cfd-core` [patch]: Replace nalgebra-style identity and square-root
  dispatch in `physics::fluid::{traits,validation,temperature,properties,
  newtonian}` with Eunomia `NumericElement`. Fluid state Mach-number,
  validation, temperature-dependent fluid, base fluid-property, and
  Newtonian/ideal-gas checks now use provider-owned `ZERO`/`ONE` identities
  and `NumericElement::sqrt` instead of `T::zero()`, `T::one()`, or method
  `sqrt()`. Verified with `cargo fmt -p cfd-core --check`, `cargo check -p
  cfd-core --no-default-features`, `cargo check -p cfd-core`, `cargo nextest
  run -p cfd-core --no-default-features fluid` (45/45 passed, 153 skipped),
  `cargo check -p cfd-1d`, `cargo nextest run -p cfd-1d resistance`
  (173/173 passed, 555 skipped), and a focused scan showing no old identity or
  method-sqrt hits in the touched fluid files. Residual `nalgebra::RealField`
  imports remain the broader fluid trait/storage contract for a later
  cross-crate migration.
- [x] `cfd-core` [patch]: Replace
  `physics::fluid::non_newtonian` direct `num_traits::FromPrimitive`
  construction and silent scalar fallback conversions with Eunomia
  `FloatElement`/`NumericElement`. Bingham, Casson, Carreau-Yasuda,
  Herschel-Bulkley, and Power-law models now use provider-owned scalar
  constants, zero/one identities, and `sqrt`/`powf`/`exp` dispatch while
  preserving the existing fluid trait API. Verified with `cargo fmt -p
  cfd-core --check`, `cargo check -p cfd-core --no-default-features`, `cargo
  check -p cfd-core`, `cargo nextest run -p cfd-core --no-default-features
  non_newtonian` (5/5 passed, 193 skipped), `cargo check -p cfd-1d`, `cargo
  nextest run -p cfd-1d resistance` (173/173 passed, 555 skipped), and a
  focused residue scan showing no direct `num_traits`, `FromPrimitive`,
  `T::from_f64`, old identity methods, or conversion-fallback hits under
  `physics/fluid/non_newtonian`. Residual `nalgebra::RealField` remains the
  shared fluid trait scalar boundary for a later cross-crate migration.
- [x] `cfd-1d` [patch]: Centralize `physics::resistance` scalar-provider
  construction and conversion through `ResistanceScalar`, `scalar_from_f64`,
  and `scalar_to_f64` in `models::traits`. Resistance models, calculator,
  factory, and geometry code no longer carry direct per-file
  `nalgebra::RealField`, `num_traits::FromPrimitive`, `T::from_f64`, or
  silent `nalgebra::try_convert(...).unwrap_or(...)` conversion fallbacks.
  Verified with `cargo fmt -p cfd-1d`, `cargo check -p cfd-1d`, `cargo
  nextest run -p cfd-1d resistance` (173/173 passed, 555 skipped), and a
  focused residue scan showing only the centralized
  `physics/resistance/models/traits.rs` bridge. Residual full Eunomia
  replacement is sequenced behind migration of
  `cfd_core::physics::fluid::FluidTrait<T>` away from its inherited
  `nalgebra::RealField` contract.
- [x] `cfd-3d` [patch]: Replace the turbulence model scalar-provider seam for
  Smagorinsky, Dynamic Smagorinsky, Dynamic Gradient Smagorinsky, AMD, Vreman,
  WALE, Sigma, DES, Spalart-Allmaras, Mixing Length, and shared SGS energy
  helpers with Eunomia `FloatElement`/`NumericElement`. Sigma now depends on
  provider-owned `FloatElement::acos` instead of direct `num_traits::Float`
  inverse cosine. Verified with `cargo check -p cfd-3d`, `cargo nextest run
  -p cfd-3d turbulence` (46/46 passed), rustfmt, and focused residue scans
  across `crates/cfd-3d/src/physics/turbulence`; no direct `RealField`,
  `FromPrimitive`, `num_traits::Float`, old `T::one()`/`T::zero()`, or
  `T::from_f64` residue remains in that module.
- [x] `cfd-3d` [patch]: Replace the turbulence helper scalar seam in
  `physics::turbulence::{field_ops,filter_ops}` from direct
  nalgebra/num-traits contracts to Eunomia `NumericElement`/`CastFrom<i32>`
  plus Leto `Vector3`. Gradient, strain, vorticity, filter-moment, and
  resolved-stress helpers no longer import or bound on `nalgebra::RealField`,
  `num_traits::FromPrimitive`, or `num_traits::Float`. Verified with `cargo
  check -p cfd-3d`, `cargo nextest run -p cfd-3d turbulence` (46/46 passed),
  rustfmt, focused residue scans, and touched-file diff checks. The later
  model-level pass closed the remaining turbulence scalar-provider residue.
- [x] `cfd-core`/`cfd-3d`/`cfd-validation` [patch]: Replace shared
  fluid-dynamics velocity storage from nalgebra `Vector3` to Leto
  `geometry::Vector3`, and route the affected field abstractions through
  Eunomia `NumericElement`. Direct spectral DNS/forcing/diagnostics,
  turbulence field/filter helpers, IBM NUFFT construction boundaries, and 3D
  validation benchmark consumers now build Leto-backed `VelocityField`
  components. Analytical validation vectors and IBM marker/probe APIs that
  still use nalgebra convert at the boundary. Verified with `cargo check -p
  cfd-core`, `cargo check -p cfd-3d`, `cargo check -p cfd-validation`,
  focused `cargo nextest run` filters for fluid-dynamics operations, spectral,
  NUFFT, Taylor-Green, forced-turbulence, and NUFFT-coupling coverage,
  rustfmt, residue scans, and `git diff --check`. Residual Leto migration
  remains in unrelated nalgebra domains such as FEM, VOF, cascade, analytical
  validation vectors, and dense matrix/vector solver surfaces.
- [x] `cfd-core` [patch]: Replace
  `physics::fluid_dynamics::operations` direct `num_traits::FromPrimitive`
  scalar construction with Eunomia `NumericElement`. Vorticity, divergence,
  kinetic-energy, and enstrophy operations now use provider-owned identities
  and half/two factors while preserving current nalgebra `Vector3` storage and
  Moirai parallel slice execution. Added value-semantic constant-flow,
  affine-divergence, and kinetic-energy tests. Verified with `cargo check -p
  cfd-core`, `cargo nextest run -p cfd-core fluid_dynamics::operations` (3/3
  passed), rustfmt, focused residue scans, and touched-file diff checks.
  Residual fluid-dynamics provider work is now the nalgebra vector/storage
  boundary in `fields.rs`, `operations.rs`, RANS/turbulence traits, and
  Rhie-Chow.
- [x] `cfd-core` [patch]: Replace `physics::fluid_dynamics::service`
  pipe-flow direct `num_traits::{Float, FromPrimitive}` scalar construction
  and math dispatch with Eunomia `FloatElement`/`NumericElement`. Reynolds and
  Prandtl helpers no longer require num-traits; pipe pressure-drop and
  friction-factor formulas now use Eunomia constants, powers, square roots,
  logarithms, and absolute convergence checks. Added value-semantic laminar
  friction-factor coverage and preserved Colebrook-White validation. Verified
  with `cargo check -p cfd-core`, `cargo nextest run -p cfd-core
  fluid_dynamics::service` (2/2 passed), rustfmt, focused residue scans, and
  touched-file diff checks. Residual fluid-dynamics provider work remains in
  `operations.rs` and the nalgebra storage/vector boundary.
- [x] `cfd-core` [patch]: Replace `physics::fluid_dynamics::flow_regimes`
  direct `nalgebra::RealField`/`num_traits::ToPrimitive` conversion with
  Eunomia `RealField`/`NumericElement`. Reynolds, Mach, and combined flow
  classification now use `NumericElement::to_f64`; the old silent
  `unwrap_or(0.0)` conversion fallback is removed; and
  `FluidDynamicsService::flow_regime` now exposes the migrated Eunomia scalar
  bound. Added value-semantic threshold tests and verified with `cargo check
  -p cfd-core`, `cargo nextest run -p cfd-core flow_regime` (3/3 passed),
  rustfmt, focused residue scans, and touched-file diff checks. Residual
  fluid-dynamics provider work remains in `operations.rs` and the
  service-level pipe-flow friction-factor formulas.
- [x] `cfd-3d` [patch]: Replace the Apollo/Leto spectral diagnostics helper's
  direct `num_traits::{FromPrimitive, ToPrimitive}` scalar conversion with
  Eunomia `NumericElement`. Kinetic-energy and enstrophy spectra now convert
  `VelocityField` components into Leto arrays for Apollo FFTs through
  `NumericElement::to_f64`, and the obsolete fallible conversion branch was
  removed. Verified with `cargo check -p cfd-3d`, `cargo nextest run -p
  cfd-3d diagnostics` (5/5 passed), rustfmt, focused residue scans, and
  touched-file diff checks. Residual work remains in `cfd-core`
  `VelocityField<T>`'s nalgebra `Vector3`/`RealField` storage boundary and the
  larger `cfd-3d` spectral `DMatrix`/`DVector` surfaces.
- [x] `cfd-2d` [patch]: Close the remaining LBM MRT and Carreau-Yasuda
  scalar-provider holdouts. `RelaxationMatrix`, `MrtCollision`, the local
  Carreau-Yasuda rheology model, and its BGK collision operator now use
  Eunomia `FloatElement`/`NumericElement`, local scalar helpers, and
  `CastFrom<i32>` for lattice velocities instead of direct `nalgebra`,
  `num_traits`, or legacy generic scalar constructors. Verified with `cargo
  check -p cfd-2d`, `cargo nextest run -p cfd-2d lbm` (31/31 passed), `cargo
  nextest run -p cfd-2d carreau_yasuda` (5/5 passed), rustfmt, and focused
  residue scans. Residual provider migration remains broad outside the LBM
  scalar-provider seam.
- [x] `cfd-2d` [patch]: Replace the LBM macroscopic/D2Q9/BGK scalar seam
  from direct `nalgebra::RealField` and `num_traits::FromPrimitive` to
  Eunomia `FloatElement` and local scalar helpers. Density, velocity,
  pressure, stress, kinetic-energy, vorticity, equilibrium constants, and BGK
  viscosity/omega construction now use provider-owned scalar APIs. MRT and
  Carreau-Yasuda `CollisionOperator` impls were bridged to the migrated trait
  seam; the later MRT/Carreau-Yasuda follow-up closes their broader
  rheology/moment scalar-provider residue. Verified with `cargo check -p cfd-2d`, `cargo nextest run -p
  cfd-2d lbm` (31/31 passed), rustfmt, focused residue scans, and touched-file
  diff checks. Residual provider migration remains broad across direct
  `nalgebra`, `num_traits`, `ndarray`, raw `wgpu`, and execution-provider
  surfaces.
- [x] `cfd-3d` [patch]: Replace near-wall Spalding wall-law scalar dispatch
  from direct `nalgebra::RealField` and `num_traits::{FromPrimitive, Float}`
  to Eunomia `FloatElement`/`NumericElement`. Wall-law constants,
  `exp`/`ln`/`sqrt`, convergence `abs`, and nonnegative clamps now use
  provider-owned scalar APIs. Verified with `cargo check -p cfd-3d`, `cargo
  nextest run -p cfd-3d wall_functions` (3/3 passed), rustfmt,
  touched-file diff checks, and focused residue scans. Residual provider work
  remains across other turbulence models, nalgebra geometry, direct WGPU, and
  execution-provider surfaces. Direct `GpuContext` Hephaestus replacement is
  sequenced after resolving the current `wgpu@0.19.4`/`wgpu@26.0.1` split.
- [x] `cfd-3d`/`cfd-2d`/`cfd-validation` [patch]: Replace the Apollo-backed
  Fourier wrapper's direct `num_traits` scalar conversion seam with Eunomia
  `FloatElement`/`NumericElement`, and propagate `FloatElement` bounds through
  downstream blood/cavitation consumers that compile the focused 3D spectral
  test target. DES turbulence and the touched validation blood benchmark
  surfaces now use provider-owned scalar conversion where they call migrated
  `cfd-core` blood APIs. Verified with `cargo check -p cfd-3d`, `cargo check
  -p cfd-validation`, `cargo nextest run -p cfd-3d --test
  fourier_validation` (12/12 passed), rustfmt, and focused residue scans.
  Residual provider migration remains broad across other direct `num_traits`,
  `nalgebra`, direct `wgpu`, and execution-provider surfaces.
- [x] `cfd-core` [patch]: Replace
  `physics::fluid::blood::FahraeuasLindqvist` direct
  `num_traits::FromPrimitive` scalar construction and direct generic math with
  Eunomia `FloatElement`/`NumericElement`. The Pries/Secomb relative-viscosity
  formulas, shared `mu_45` fit, relative-viscosity clamp, and tube-hematocrit
  correlation now use provider-owned constants, zero/one identities, `powf`,
  `exp`, `abs`, and max dispatch. Verified with `cargo check -p cfd-core`,
  focused Fåhræus-Lindqvist and blood nextest filters, rustfmt, residue scans,
  and touched-file diff checks. Full cfd-core clippy remains blocked by
  unrelated existing lints in boundary applicator and Rhie-Chow tests. Local
  blood-model scalar `num_traits` construction is closed; the broader fluid
  trait `RealField` boundary remains later provider work.
- [x] `cfd-core` [patch]: Replace
  `physics::fluid::blood::{CassonBlood, CarreauYasudaBlood, BloodModel}`
  direct `num_traits::FromPrimitive` construction and dispatch bounds with
  Eunomia `FloatElement`/`NumericElement`. Casson constructors, hematocrit
  scaling, temperature correction, square-root viscosity formula, validation
  constants, Carreau-Yasuda constants/powers, and the shared blood-model
  selector now use provider-owned scalar construction and math dispatch.
  Verified with `cargo check -p cfd-core`, focused Casson/Carreau-Yasuda/blood
  nextest filters, rustfmt, residue scan, and touched-file diff checks. Full
  cfd-core clippy remains blocked by unrelated existing lints in boundary
  applicator and Rhie-Chow tests. Residual blood-fluid scalar holdout remains
  in Fåhræus-Lindqvist plus the broader fluid trait `RealField` boundary.
- [x] `cfd-core` [patch]: Replace `physics::fluid::blood::CrossBlood` direct
  `num_traits::FromPrimitive` scalar construction with Eunomia
  `FloatElement`/`NumericElement`. The Cross blood constructor and apparent
  viscosity formula now use provider-owned constants, zero/one identities, and
  `powf` dispatch while preserving the broader fluid-trait `RealField`
  boundary for a later slice. Verified with `cargo check -p cfd-core`, focused
  Cross and blood nextest filters, rustfmt, residue scan, and touched-file diff
  checks. Residual blood-fluid scalar holdouts remain in Casson,
  Carreau-Yasuda, Fåhræus-Lindqvist, `BloodModel`, and the fluid trait.
- [x] `cfd-core` [patch]: Close the remaining cavitation scalar-provider
  holdouts. `CavitationModel<T>`/`ZgbParams<T>` now use Eunomia
  `FloatElement`/`NumericElement`, and
  `heterogeneous_inception_threshold_pa` now exposes the `f64` contract that
  its selective cavitation model already computes instead of a fake generic
  `RealField` adapter. Added value-semantic tests for the migrated mass
  transfer models and the heterogeneous legacy adapter. Verified with `cargo
  check -p cfd-core`, `cargo nextest run -p cfd-core cavitation` (40/40
  passed), focused cavitation residue scans, rustfmt, and touched-file diff
  checks. Residual provider migration work remains outside cavitation across
  broader CFDrs nalgebra/ndarray/wgpu/rayon/tokio surfaces.
- [x] `cfd-core` [patch]: Replace
  `cavitation::nuclei_transport` scalar bounds and constants from nalgebra
  `RealField` to Eunomia `FloatElement`/`NumericElement`. Existing value tests
  continue to cover affine nuclei pressure coupling, dissolution rate,
  exponential transit decay, and diffusion-coefficient access after the
  provider migration. Verified with `cargo check -p cfd-core`, `cargo nextest
  run -p cfd-core cavitation` (35/35 passed), focused residue scans, rustfmt,
  and touched-file diff checks. Residual cavitation scalar holdouts were then
  `models.rs` and `heterogeneous_nucleation.rs`; the current cavitation
  closeout removed both holdouts.
- [x] `cfd-core` [patch]: Replace `cavitation::VenturiCavitation` scalar
  bounds and constants from nalgebra `RealField` plus direct
  `num_traits::FromPrimitive` to Eunomia `FloatElement`/`NumericElement`.
  Existing Venturi closed-form tests continue to cover continuity, Bernoulli,
  cavitation number, cavity length, and conical cavity volume after the
  provider migration. Verified with `cargo check -p cfd-core`, `cargo nextest
  run -p cfd-core cavitation` (35/35 passed), focused residue scans, rustfmt,
  and touched-file diff checks. At that point residual cavitation scalar
  holdouts were `models.rs`, `nuclei_transport.rs`, and
  `heterogeneous_nucleation.rs`; the subsequent nuclei-transport slice removed
  `nuclei_transport.rs` from the active residual list.
- [x] `cfd-core` [patch]: Replace the active cavitation Rayleigh/regime,
  biological-damage, cavitation-number, and material-damage scalar surfaces
  from nalgebra `RealField` plus direct `num_traits::FromPrimitive` to Eunomia
  `FloatElement`/`NumericElement`. The migrated files now use provider-owned
  constants and math operations, and added closed-form value tests for the
  newly migrated number/damage APIs. Verified with `cargo check -p cfd-core`,
  `cargo nextest run -p cfd-core cavitation` (35/35 passed), focused residue
  scans, rustfmt, and touched-file diff checks. At that point residual
  cavitation `RealField` surfaces remained in `models.rs`, `venturi.rs`,
  `nuclei_transport.rs`, and `heterogeneous_nucleation.rs`; the subsequent
  Venturi slice removed `venturi.rs` from the active residual list. Full clippy
  is currently blocked by unrelated existing cfd-core lints outside this
  slice.
- [x] `cfd-core` [patch]: Replace hemolysis calculator and platelet activation
  scalar bounds from nalgebra `RealField` plus direct `num_traits`
  construction to Eunomia `FloatElement`/`NumericElement`. The migrated
  hemolysis path now uses provider-owned zero/one/constants and Eunomia
  exponential evaluation while preserving value-semantic NIH, MIH,
  exposure-time, and platelet activation tests. Residual core `RealField`
  surfaces remain outside hemolysis.
- [x] `cfd-2d` [patch]: Replace the structured/unstructured grid center and
  adaptive-resolution boundary with Leto/Eunomia providers. `Grid2D` now
  returns `leto::geometry::Vector2`, `StructuredGrid2D` and
  `AdaptiveGrid2D::effective_resolution` use Eunomia scalar construction, and
  `FvmSolver` no longer carries `nalgebra::RealField` through the structured
  grid contract. Validation benchmark and FDM/LBM consumers now read Leto
  vectors by index. Verified with `cargo check -p cfd-2d`, `cargo check -p
  cfd-validation`, focused grid/FVM nextest runs, residue scans, rustfmt, and
  touched-file diff checks. Residual `cfd-2d` nalgebra/num-traits usage remains
  in momentum, pressure-velocity, SIMPLE/PIMPLE, and other non-grid solver
  surfaces.
- [x] `cfd-2d`/`leto` [patch]: Replace FVM face centers, normals, flux
  velocity arguments, and solver velocity fields from nalgebra `Vector2` to
  `leto::geometry::Vector2`. Extend Leto with a `Vector2` alias plus generic
  fixed-vector norm/normalization operations, and add direct Leto consumption
  to `cfd-2d`. Verified with provider and consumer compile checks, focused
  nextest coverage, residue scans showing no nalgebra Vector2 or num-traits
  holdouts in FVM, and touched-file diff checks. The later grid/FVM provider
  slice removed the inherited `StructuredGrid2D<T>` `RealField` blocker.
- [x] `cfd-2d` [patch]: Replace `solvers::fvm` direct
  `num_traits::FromPrimitive` scalar constants with Eunomia
  `FloatElement`/`NumericElement` across config defaults, solver face-center
  construction, and flux schemes. Add value-semantic default-config coverage
  and non-finite diffusion rejection coverage while preserving the current
  nalgebra `RealField` vector/flux boundary for a later Leto migration. Static
  residue and formatting checks are clean; `cargo check -p cfd-2d` and
  `cargo check -p cfd-3d` pass, and `cargo nextest run -p cfd-2d --lib fvm`
  passes 25/25 focused tests.
- [x] `cfd-2d` [patch]: Replace `stability::CFLCalculator` scalar bounds and
  constants from nalgebra `RealField` plus `num_traits::FromPrimitive` to
  Eunomia `FloatElement`/`NumericElement`, preserving the documented CFL and
  diffusion formulas and validating the existing value-semantic CFL tests.
- [x] `cfd-core` [patch]: Replace material solid/interface trait and value
  contracts from nalgebra `RealField` to Eunomia
  `FloatElement`/`NumericElement`. `MaterialDatabase` now carries Eunomia
  bounds for solid/interface storage and retains `RealField` only through the
  still-unmigrated `Fluid<T>` map. Added value-semantic material tests for
  isotropic shear modulus, constructor constants, wetting angles, and adhesion
  energy.
- [x] `cfd-core`/`leto` [patch]: Replace the vector-backed CFDrs velocity value
  object and physical-parameters gravity vector from nalgebra `Vector3` to
  `leto::geometry::Vector3`, adding direct Leto consumption to `cfd-core`.
  Leto geometry now derives Serde for fixed 3D geometry types so serialized
  CFDrs value objects can use provider-owned vectors without a compatibility
  wrapper. Residual `RealField` in `ProblemAggregate` remains tied to
  `Domain<T>` and fluid contracts.
- [x] `cfd-core` [patch]: Replace scalar-only physics value wrappers
  (`Temperature`, `Pressure`, `ReynoldsNumber`, and `DimensionlessNumber`)
  from direct `nalgebra::RealField`/`ComplexField` contracts to Eunomia
  `FloatElement`/`NumericElement`; propagate the necessary value-wrapper
  bounds through `PhysicalParameters`, `ProblemAggregate`,
  `InitialConditions`, and `SimulationAggregate`. Residual `Velocity` and
  vector/material/hemolysis `RealField` surfaces remain for later Leto/Gaia and
  Eunomia provider slices.
- [x] `cfd-math`/root [patch]: Remove the remaining downstream direct WGPU
  dependency edge outside `cfd-core` by exposing GPU synchronization and
  timestamp-query capability through `cfd_core::compute::gpu::GpuContext`.
  `cfd-math::linear_solver::operators::gpu` now consumes those provider-facing
  context methods instead of importing `wgpu`; `cfd-math/gpu` and the root
  package `gpu` feature no longer activate `dep:wgpu` directly. Residual raw
  WGPU kernel/buffer/pipeline ownership remains in `cfd-core::compute::gpu` for
  Hephaestus migration.
- [x] `cfd-core` [patch]: Replace direct raw WGPU adapter probing in
  `ComputeBackend::detect_gpu_support` with Hephaestus device acquisition via
  optional `hephaestus-wgpu`; verified with `cargo check -p cfd-core --features
  gpu`, `cargo nextest run -p cfd-core --features gpu --lib` (183/183 passed),
  and `cargo tree --workspace -i hephaestus-wgpu` showing the active
  `hephaestus-wgpu -> cfd-core` provider edge. Raw WGPU kernel/buffer/pipeline
  ownership remains for the next Hephaestus migration slice.
- [x] `cfd-1d` [patch]: Remove unused `sprs` and the stale root workspace
  `ndarray` declaration, eliminating the non-Python active `ndarray` graph path
  while preserving cfd-1d's existing `nalgebra-sparse` solver surface for a
  later Leto sparse migration. Propagate `eunomia::FloatElement` into two
  generic cfd-1d test helpers required by migrated cfd-core fluid constructors.
  Verified with `cargo check -p cfd-1d`, `cargo nextest run -p cfd-1d`
  (725/725 passed), and `cargo tree -p cfd-1d -i ndarray`/`cargo tree -p
  cfd-3d -i ndarray` showing no matching package. Remaining workspace
  `ndarray` is `numpy -> cfd-python`.
- [x] `cfd-3d`/`cfd-validation` [patch]: Replace CFDrs' Apollo lock resolution
  with the side-by-side Atlas Apollo provider checkout, removing active
  ndarray-backed Apollo packages from the cfd-3d graph; adapt cfd-3d FEM and
  spectral call sites plus validation generics to Eunomia scalar/complex
  contracts. Verified with `cargo check -p cfd-3d`, `cargo nextest run -p
  cfd-3d` (394/394 passed), Apollo `ndarray` source/manifest/lock scan, and
  active Apollo package trees with no `ndarray` matches. Remaining `ndarray`
  in the cfd-3d graph is the non-Apollo `sprs -> cfd-1d` path.
- [x] `cfd-math` [patch]: Replace geometric multigrid scalar construction and
  transfer weights with Eunomia `FloatElement`/`NumericElement`, remove stale
  AMG `FromPrimitive` bounds after the multigrid provider cleanup, and add
  value-semantic Poisson stencil/full-weighting restriction coverage while
  preserving current nalgebra dense/sparse surfaces for later Leto migration.
- [x] `cfd-math` [patch]: Replace multigrid interpolation direct
  `num_traits::{FromPrimitive, ToPrimitive}` scalar conversion and quality
  metric extraction with Eunomia `FloatElement`/`NumericElement`, preserving
  current nalgebra sparse/vector surfaces; verified with `cargo check -p
  cfd-math` and focused interpolation nextest.
- [x] `cfd-math` [patch]: Replace multigrid coarsening algorithm and
  quality-analysis direct `num_traits::{FromPrimitive, ToPrimitive}` scalar
  conversion with Eunomia `FloatElement`/`NumericElement`, including
  strength-matrix absolute-value dispatch and f64 quality-metric extraction;
  add value-semantic strength-matrix connectivity coverage while preserving
  current nalgebra sparse/vector surfaces and leaving interpolation/GMG for
  later provider slices.
- [x] `cfd-math` [patch]: Replace multigrid smoother scalar thresholds,
  Chebyshev eigenvalue/default constants, and immediate AMG smoother-owner
  constants/complexity filters with Eunomia `FloatElement`/`NumericElement`,
  removing direct scalar-conversion fallbacks from
  `linear_solver::preconditioners::multigrid::smoothers` and the touched AMG
  owner paths; add value-semantic smoother update/eigenvalue tests while
  preserving current nalgebra sparse/vector surfaces and deeper AMG
  `FromPrimitive` contracts for later provider slices.
- [x] `cfd-math` [patch]: Replace linear-solver convergence monitor scalar
  conversions with Eunomia `FloatElement`, removing direct
  `SafeFromF64`/`T::from_f64_or` fallback usage from
  `linear_solver::traits`; add value-semantic convergence-factor,
  CG-bound, and validation-rejection tests while preserving current nalgebra
  vector/operator surfaces for later Leto migration.
- [x] `cfd-math` [patch]: Replace Poisson/Laplacian and momentum/energy
  linear operator scalar constants and provider bounds with Eunomia
  `FloatElement`, removing direct `num_traits::FromPrimitive` and scalar
  conversion fallbacks from `linear_solver::operators::{poisson,momentum}`;
  add value-semantic center-stencil tests for the touched finite-difference
  operators while preserving current nalgebra vector surfaces for later Leto
  migration.
- [x] `cfd-math` [patch]: Replace Schwarz/IncompleteCholesky preconditioner
  provider residue with Eunomia `FloatElement`/`NumericElement`, removing
  Schwarz's stale direct `num_traits::FromPrimitive` requirement and Cholesky's
  direct `T::from_f64(...).unwrap_or_else(...)` symmetry tolerance fallback
  while preserving current nalgebra CSR/vector surfaces for later Leto
  migration.
- [x] `cfd-math` [patch]: Replace SSOR preconditioner default relaxation
  construction and provider bounds with Eunomia `FloatElement`, removing direct
  `num_traits::FromPrimitive` from `linear_solver::preconditioners::ssor`
  while preserving the current nalgebra CSR/vector preconditioner surface for
  later Leto migration.
- [x] `cfd-math` [patch]: Replace basic Jacobi/SOR preconditioner scalar
  constants, diagonal tolerances, omega construction, and absolute-value
  dispatch with Eunomia `FloatElement`/`NumericElement`, removing direct
  `num_traits::FromPrimitive` from `linear_solver::preconditioners::basic`
  while preserving the current nalgebra CSR/vector preconditioner surface for
  later Leto migration.
- [x] `cfd-math` [patch]: Replace GMRES, direct sparse solver,
  linear-solver chain, and block/SIMPLE preconditioner direct
  `num_traits::{Float, FromPrimitive, ToPrimitive}` bounds and scalar
  safeguards with Eunomia `FloatElement`/`NumericElement`, while preserving the
  current nalgebra vector/matrix and rsparse f64-backed direct solver surfaces
  for later Leto/backend migration.
- [x] `cfd-math` [patch]: Replace iterative linear-solver config default
  tolerance construction with Eunomia `FloatElement`, correct stale parallel
  SpMV wording from Rayon to Moirai, and propagate the default-construction
  provider bound through CG, BiCGSTAB, and GMRES while leaving the remaining
  linear-solver `num-traits`/nalgebra migration to later focused slices.
- [x] `cfd-math` [patch]: Replace CFD SIMD central-difference constants and
  field-norm square-root dispatch with Eunomia `FloatElement`/`NumericElement`,
  removing direct `num_traits::FromPrimitive` from `simd` while preserving the
  existing Moirai-backed parallel slice execution and nalgebra scalar API for
  later Leto/eunomia scalar-bound cleanup.
- [x] `cfd-math` [patch]: Replace sparse matrix pattern constants,
  Frobenius norm square-root dispatch, condition-estimate singular-diagonal
  thresholds, and diagonal dominance absolute-value checks with Eunomia
  `FloatElement`/`NumericElement`, removing direct
  `num_traits::{Float, FromPrimitive, Signed}` from `sparse` and correcting
  stale Rayon wording to the actual Moirai parallel slice adapter.
- [x] `cfd-math` [patch]: Replace Anderson/JFNK nonlinear solver default
  constants, perturbation epsilon math, QR/Givens absolute-value and square-root
  dispatch, EW forcing clamp math, and back-substitution diagonal checks with
  Eunomia `FloatElement`/`NumericElement`, removing direct
  `num_traits::{Float, FromPrimitive}` from `nonlinear_solver` while preserving
  the existing nalgebra vector/matrix API for later Leto migration.
- [x] `cfd-math` [patch]: Replace SIMPLE pressure-velocity default tolerance
  and relaxation constants with Eunomia `RealField`/`FloatElement`, remove
  direct nalgebra and `num_traits::FromPrimitive` from `pressure_velocity`,
  and narrow explicit `SIMPLEConfig::new` construction away from
  scalar-conversion bounds.
- [x] `cfd-math` [patch]: Replace iterator stencil coefficients and iterator
  statistics count/math conversions with Eunomia
  `FloatElement`/`NumericElement`, remove direct `num_traits::FromPrimitive`
  from `iterators`, and replace the second-derivative zero placeholder branch
  with real coefficients for every declared 3-point stencil pattern.
- [x] `cfd-math` [patch]: Replace WENO5/WENO7 epsilon defaults, linear
  weights, ENO reconstruction coefficients, and smoothness-indicator constants
  with Eunomia `FloatElement`, remove direct `num_traits::FromPrimitive` bounds
  from the WENO high-order surface, and consolidate denominator squaring through
  a local multiplication helper while preserving the existing WENO formulas.
- [x] `cfd-2d` [patch]: Replace the scheme amplification factor's direct
  `num_complex::Complex<f64>` return value with `eunomia::Complex<f64>`,
  remove direct `num-complex` manifest ownership from cfd-2d, and make the
  existing 2D solver/network scalar-construction paths state their
  `eunomia::FloatElement` requirement explicitly.
- [x] `cfd-3d`/`cfd-validation` [patch]: Repair the cfd-2d nextest dependency
  chain by making the private cfd-3d Leto/Apollo adapter target Apollo's
  reachable ndarray FFT/NUFFT API through Leto conversions, removing the stale
  `MnemosyneStorage` assumption, and propagating `FloatElement` through cfd-3d
  config plus cfd-validation MMS Richardson solver construction.
- [x] `cfd-math`/`cfd-validation` [patch]: Replace the time-stepping stability
  analyzer's public `num_complex::Complex<f64>` von Neumann callback contract
  with `eunomia::Complex<f64>`, remove direct `num-complex` manifest ownership
  from both crates, and evaluate explicit RK stability functions by forward
  substitution over the already-validated lower-triangular Butcher tableau.
- [x] `cfd-math` [patch]: Replace the stability analyzer's direct
  `num_traits::ToPrimitive` formatting/conversion path and
  `T::from_f64(...).unwrap_or_else(num_traits::Zero::zero)` scalar fallbacks
  with Eunomia `FloatElement`/`NumericElement` helpers while preserving the
  existing nalgebra Butcher-tableau surface for the later Leto matrix migration.
- [x] `cfd-math` [patch]: Replace Runge-Kutta scalar constants with Eunomia
  `FloatElement` and correct `LowStorageRK4` to the Carpenter-Kennedy 2N
  residual recurrence, adding value-semantic exponential-decay and zero-RHS
  preservation coverage.
- [x] `cfd-math` [patch]: Replace adaptive time-stepper defaults, PI controller
  constants, and Dormand-Prince tableau coefficients with Eunomia
  `FloatElement`, preserving the existing nalgebra `DVector` stepper API for
  the later Leto vector migration.
- [x] `cfd-math` [patch]: Replace Runge-Kutta-Chebyshev defaults, coefficient
  recurrence constants, adaptive error-control constants, stage/vector length
  conversions, and math dispatch with Eunomia `FloatElement`/`NumericElement`,
  preserving the existing nalgebra `DVector` stepper API for the later Leto
  vector migration.
- [x] `cfd-math` [patch]: Replace IMEX Newton tolerance, ARS343 gamma/delta
  constants, tableau coefficients, and solution weights with Eunomia
  `FloatElement`/`NumericElement`, preserving the existing nalgebra
  `DVector`/`DMatrix` stepper API for the later Leto matrix/vector migration.
- [x] `cfd-math` [patch]: Replace exponential time-differencing ERK4
  constants, phi-function small-argument thresholds, and scaling/squaring
  factorial conversions with Eunomia `FloatElement`, preserving the existing
  nalgebra `DVector`/`DMatrix` exponential-integrator API for the later Leto
  matrix/vector migration.
- [x] `cfd-math` [patch]: Replace integration quadrature constants,
  interval-count conversions, adaptive tolerance/error dispatch, and
  tetrahedral quadrature constants with Eunomia `RealField`/
  `FloatElement`/`NumericElement`; the existing `Quadrature` traits now bind
  on Eunomia scalar contracts and the integration cone has no nalgebra source
  residue.
- [x] `cfd-math` [patch]: Replace finite-difference and gradient scalar
  constants, SIMD helper scalar staging, and 2D Laplacian constants with
  Eunomia `FloatElement`/`NumericElement`; follow-up Leto migration replaced
  nalgebra `DVector`/`Vector3` differentiation result surfaces with Leto
  `Array1`/`Vector3` and removed the type-suffixed SIMD helper name.
- [x] `cfd-math` [patch]: Replace interpolation trait, linear, Lagrange, and
  cubic-spline scalar contracts with Eunomia `RealField`/`FloatElement`,
  removing nalgebra `RealField` source residue from `interpolation` and
  rejecting duplicate Lagrange nodes before denominator division.
- [x] `cfd-core` [patch]: Replace Rhie-Chow interpolation default constants and helper constants in `physics/fluid_dynamics/rhie_chow.rs` with Eunomia `FloatElement`, remove direct `num_traits::FromPrimitive`, and add value-semantic pressure-correction tests for u-face and v-face interpolation.
- [x] `cfd-core` [patch]: Replace boundary geometry measure constants in `physics/boundary/geometry.rs` with Eunomia `FloatElement`, remove silent scalar-conversion fallbacks, preserve narrower conversion-free bounds for `contains_point` and `dimension`, and add value-semantic line/sphere/cylinder/unsupported-measure tests.
- [x] `cfd-core` [patch]: Replace boundary ghost-cell constants and Robin singularity reporting in `physics/boundary/ghost_cells.rs` with Eunomia `FloatElement`/`NumericElement`, remove direct `num_traits::{FromPrimitive, ToPrimitive}`, and reject degenerate Robin coefficients before division.
- [x] `cfd-core` [patch]: Replace staggered-grid coordinate scalar conversions in `geometry/staggered.rs` with Eunomia `FloatElement`, removing direct `num_traits::FromPrimitive` while asserting exact representability for integer grid dimensions and indices.
- [x] `cfd-core` [patch]: Replace boundary time-function and ghost-cell scalar conversions in `physics/boundary/{time_dependent,applicator}.rs` with Eunomia `FloatElement`, route trigonometric/exponential dispatch through Eunomia, and propagate the provider bound through boundary applicator/specification/manager surfaces.
- [x] `cfd-core` [patch]: Replace temperature-dependent fluid `FromPrimitive` bounds and Sutherland exponent conversion in `physics/fluid/temperature.rs` with Eunomia `FloatElement`, and add value-semantic polynomial, Andrade, and Sutherland tests.
- [x] `cfd-core` [patch]: Replace hemolysis calculator and platelet activation scalar constants in `physics/hemolysis/{calculator,trauma}.rs` with Eunomia `FloatElement`, remove direct `FromPrimitive` fallback conversions, and fix platelet activation to use the decaying exponential probability contract.
- [x] `cfd-core` [patch]: Replace mesh quality threshold conversions in `geometry/mesh/quality.rs` with Eunomia `FloatElement`, remove direct `FromPrimitive` fallback thresholds, and add value-semantic quality-level/recommendation tests.
- [x] `cfd-core` [patch]: Replace CPU backend domain-parameter scalar conversions in `compute/cpu.rs` with Eunomia `FloatElement`, delete the local silent fallback conversion helper, and narrow CPU buffer impl bounds away from scalar conversion requirements.
- [x] `cfd-core` [patch]: Replace time integrator scalar constants in `compute/time/integrators.rs` with Eunomia `FloatElement`, remove silent default tolerance conversion fallbacks, and add value-semantic explicit/implicit integrator regression tests.
- [x] `cfd-core` [patch]: Replace time-step controller constants and runtime math helpers in `compute/time/controllers.rs` with Eunomia `FloatElement`, and surface invalid integration order as typed errors instead of silent fallback.
- [x] `cfd-core` [patch]: Replace solver configuration default constants in `compute/solver/config.rs` with Eunomia `FloatElement`, consolidate duplicated builder/default construction, and add value-semantic config default parity tests.
- [x] `cfd-core` [patch]: Replace abstraction default constants in `abstractions/{state,problem}.rs` with Eunomia `FloatElement`, removing direct `FromPrimitive` conversions while keeping non-constructor field-state methods free of scalar-conversion bounds.
- [x] `cfd-core` [patch]: Replace fluid validation threshold conversions with Eunomia `FloatElement`, removing direct `FromPrimitive` fallback constants from `physics/fluid/validation.rs`.
- [x] `cfd-core` [patch]: Replace material/fluid constant constructors and the material database dependency cone with Eunomia `FloatElement`, removing direct `FromPrimitive` scalar conversions and silent zero/one fallback constants from the touched cone.
- [x] `cfd-core` [patch]: Replace physics value-object and dependent management aggregate scalar constant conversions with Eunomia `FloatElement`, remove direct `FromPrimitive` bounds from the touched cone, and stop silently substituting default physical parameter conversion values.
- [x] `cfd-python` [patch]: Replace blood-model PyO3 direct `num-traits` conversions with Eunomia `NumericElement`, remove the direct manifest dependency, and expose actual Rust model field values instead of fallback defaults in Python getters.
- [x] `cfd-io` [patch]: Replace direct `num-traits` scalar bounds and conversions in checkpoint validation, binary helpers, and CSV helpers with Eunomia `RealField`; remove the direct cfd-io `num-traits` dependency and make oversized/zero mesh dimensions fail mass-conservation validation instead of silently falling back to unit spacing.
- [x] `cfd-io` [arch]: Replace checkpoint and binary dense payload ownership with Leto `Array1`/`Array2`, remove direct `nalgebra`/`DMatrix`/`DVector`/`RealField` source usage, and decouple normal cfd-io dependencies from `cfd-core`/`cfd-math` so `cargo tree -p cfd-io -e normal -i nalgebra` has no matching package.
- [x] `cfd-python`: Replace direct `ndarray`/`nalgebra` ownership in 2D PyO3 NumPy-return helper paths with Leto `Array2`; keep NumPy as the external Python boundary representation only.
- [x] `cfd-core`/`cfd-math`: Replace direct `pollster` blocking in GPU context creation, GPU support detection, Poisson residual readback, and cfd-math GPU operator dispatch metrics with Moirai's blocking boundary; remove `pollster` from the workspace graph. Direct Hephaestus GPU-provider replacement remains sequenced behind the `wgpu 0.19` -> `wgpu 26.0` API normalization.
- [x] `cfd-3d`: Replace direct `ndarray` ownership in spectral DNS/forcing/diagnostics/Fourier and IBM NUFFT paths with Leto arrays; keep the residual ndarray conversion isolated in a private Apollo adapter against Apollo's reachable ndarray API until Apollo publishes Leto-native signatures.
- [x] `cfd-schematics`: Consolidate shared components and implement unified access.
- [x] `cfd-schematics`: Enforce Clean Architecture, Dependency Inversion Principle, Single Responsibility Principle, and Single Source of Truth.
- [x] `cfd-schematics`: Replace silent no-op plotters drawer methods with real rendering and output-verified tests.
- [x] `cfd-ui`: Preserve clip-plane slot identity and explicit stale-state errors in clipping commands.
- [x] `cfd-1d`: Preserve transient composition sampling error context in the time-config projection path.
- [x] `cfd-1d`: Preallocate merged timepoint schedules in the transient composition time-config path.
- [x] `cfd-1d`: Remove the silent 1% Quemada viscosity floor and fall back to Secomb when the Quemada domain is invalid.
- [x] `cfd-1d`: Replace simplified margination lift aggregation with separate wall-induced and shear-gradient inertial scaling plus reference-equilibrium tests.
- [x] `cfd-1d`: Remove the margination singular wall cutoff and silent lateral-position clamp by using a bounded inertial-lift envelope with value-semantic regression data.
- [x] `cfd-1d`: Make droplet occupied-channel snapshots a documented projection of finite-length occupancy spans.
- [x] `cfd-1d`: Remove stored point-droplet occupied-channel state so finite-length occupancy spans are the single authoritative droplet representation.
- [x] `cfd-1d`: Restore coupled pressure-event hematocrit flow by finite zero-flow blood viscosity initialization, row-equilibrated pressure solves, and accumulated duplicate sparse entries.
- [x] `cfd-1d`: Route compact plasma-skimming hematocrit through the threshold-aware Pries phase-separation model with Murray-inferred sibling geometry.
- [x] `cfd-1d`: Remove the short-channel Reynolds floor from the Durst entrance correction and use the published low-Re formula directly.
- [x] `cfd-2d`: Reuse a persistent pressure-velocity state workspace and validate repeated SIMPLE iterations.
- [x] `cfd-2d`: Remove the transient PIMPLE outer residual snapshot allocation by reusing the corrected velocity workspace.
- [x] `cfd-2d`: Reset Rhie-Chow coefficient caches on every update and remove the dead SIMPLEC diagonal workspace.
- [x] `cfd-2d`: Wire the pressure-velocity solver boundary-condition and viscosity inputs into the reused state workspace.
- [x] `cfd-2d`: Reject mismatched initial state layouts and non-physical fluid inputs before pressure-velocity stepping.
- [x] `cfd-2d`: Replace the approximate MUSCL3/QUICK right-state reconstruction with exact quadratic face interpolation and value-semantic tests.
- [x] `cfd-2d`: Replace the unimplemented turbulence validation benchmark branch with typed supported-model dispatch and rejection tests.
- [x] `cfd-2d`: Replace serpentine mixing's exponential estimate with the Neumann eigenfunction solution for transverse diffusion and expose nonzero analytical L90/t90 in the discretized solver result.
- [x] `cfd-2d`: Remove silent Pries plasma-skimming clamps and expose checked value-semantic phase-separation evaluation.
- [x] `cfd-2d`: Replace WALE boundary zero-gradient assumptions with second-order one-sided finite-difference gradients.
- [x] `cfd-2d`: Replace Smagorinsky LES zero TKE/dissipation placeholders, zero boundary strain, and default SGS viscosity floor with Yoshizawa SGS diagnostics and second-order strain recovery.
- [x] `cfd-2d`: Remove residual validation-only Smagorinsky SGS viscosity floors and keep LES validation on the physical zero-floor default.
- [x] `cfd-math`: Preserve direct sparse solver conversion failures and reject non-finite fallback output.
- [x] `cfd-2d`/`cfd-3d`/`cfd-validation`/`cfd-python` [patch]: Clear warning-as-error diagnostics that blocked `cargo clippy -p cfd-python --all-targets -- -D warnings` (`dead_code` continuity helpers, modulo-one iteration check, collapsible wall-type matches, assign-op patterns, slice indexing in a fixed forcing spectrum, uninlined format args, and an elidable PyO3 helper lifetime).
- [x] `cfd-optim`: Remove the sonosensitizer activation floor from SDT report metrics and enforce zero-dose propagation at zero cavitation.
- [x] `cfd-optim`: Remove synthetic 1% synergy floors from objective and pool scoring paths.
- [x] `cfd-optim`: Remove synthetic 1% floors from venturi selectivity metric aggregation.
- [x] `cfd-optim`: Make hard-constraint candidate scoring return exact zero for infeasible designs.
- [x] `cfd-optim`: Remove the synthetic non-cavitating floor from the cavitation score.
- [x] `cfd-optim`: Remove the remaining combined-mode, leukapheresis, hydrodynamic-cavitation, and smooth-penalty score floors.
- [x] `cfd-optim`: Rephrase the Milestone 12 report to foreground cancer-cell preferential lysis and healthy-cell protection in the hydrodynamic cavitation narrative.
- [x] `cfd-optim`: Normalize the Milestone 12 report template to use `healthy_cell_protection_index` in the GA scoring tables and summary formulas.
- [x] `cfd-optim`: Refactor objective, search, and venturi scoring paths to consume `healthy_cell_protection_index` directly.
- [x] `cfd-optim`: Normalize Milestone 12 report writers to label the WBC column as recovery and emit the healthy-cell composite consistently.
- [x] `cfd-optim`: Remove heuristic floor values from Milestone 12 SVG report scaling guards.
- [x] `cfd-optim`: Replace residual shielding language in the Milestone 12 narrative with cancer-selective lysis and healthy-cell protection terminology.
- [x] `cfd-optim`: Replace the remaining Milestone 12 process-figure shielding wording and narrative-template exposure wording with healthy-cell protection terminology.
- [x] `cfd-optim`: Clear the Milestone 12 asset-review gate after manual PNG review and rerun the authoritative Option 2 report pipeline.
- [x] `cfd-optim`: Unify GA convergence reporting on the trailing fitness window so the figure and narrative share one trend contract.
- [x] `cfd-optim`: Add validation evidence and artifact traceability to the Milestone 12 report and synchronize the generated narrative/results artifacts.
- [x] `cfd-optim`: Emit the authoritative Milestone 12 narrative in one pass so the release report run completes after manual asset review.
- [x] `cfd-optim`: Correct Milestone 12 cross-mode therapy utility so ultrasound-only Option 1 receives acoustic-cavitation delivery credit instead of the separation-only cap.
- [x] `cfd-3d`: Replace the cavitation-number velocity floor with an exact zero-dynamic-pressure branch that returns infinite cavitation number.
- [x] `cfd-3d`: Fall back to first-order upwind when the WENO5-Z stencil is unavailable on small grids.
- [x] `cfd-3d`: Scale PLIC plane-bisection tolerance with the cell dimensions and reuse precomputed normals for curvature axis selection.
- [x] `cfd-3d`: Remove the 1% void-fraction damage cutoff and accumulate cavitation damage from any nonzero void fraction.
- [x] `cfd-3d`: Classify VOF mixed cells by the exact 0 < α < 1 criterion and drive compression strength from the configured coefficient.
- [x] `cfd-core`: Add Rayleigh-collapse-time-sensitive cavitation bio-damage thresholds for membrane injury and cell death.
- [x] `cfd-core`: Remove the minimum bubble-radius floor and treat collapsed Rayleigh-Plesset bubbles as an absorbing state.
- [x] `cfd-3d`: Route bubble-dynamics updates through the canonical Rayleigh-Plesset absorbing-collapse state instead of a radius floor.
- [x] `cfd-3d`: Reject zero dynamic pressure in the 3D Venturi pressure-coefficient calculation and remove the slice-weight floor.
- [x] `cfd-3d`: Remove the lower hematocrit floor from the cascade viscosity ratio and reuse Picard viscosity buffers.
- [x] `cfd-3d`: Consume the local liquid-density input in the bubble-dynamics Rayleigh-Plesset update and reject nonphysical densities.
- [x] `cfd-core`: Replace the re-equilibrating ambient-pressure seed in Rayleigh-Plesset bubble acceleration with a fixed far-field reference pressure.
- [x] `cfd-3d`: Correct the bubble-dynamics integration test to measure change from the configured initial radius and accept collapsed bubbles as valid absorbing states.
- [x] `cfd-core`: Expose the nuclei diffusion coefficient for solver coupling.
- [x] `cfd-3d`: Couple nuclei diffusion into the cavitation transport update with conservative finite-volume spreading.
- [x] `cfd-3d`: Reuse the cavitation-source workspace across solver steps instead of allocating a fresh source matrix per step.
- [x] `cfd-3d`: Validate flow-field dimensions before cavitation-VOF stepping even when bubble dynamics are disabled.
- [x] `cfd-3d`: Validate pressure and density dimensions in cavitation damage accumulation and walk the damage field via raw column slices.
- [x] `cfd-3d`: Validate cavitation-source dimensions, switch to raw-slice accumulation, and clamp source updates to feasible bounds.
- [x] `cfd-3d`: Validate nuclei transport dimensions and switch cavitation nuclei advection-diffusion to slice-based accumulation.
- [x] `cfd-3d`: Bind the Apollo-backed periodic DNS stepper to a validated reusable `FftPlan3D` instead of per-transform shape dispatch.
- [x] `cfd-3d`: Replace LES turbulent-kinetic-energy aliases with a documented Yoshizawa SGS energy relation shared across eddy-viscosity models.
- [x] `cfd-3d`: Replace the Spalart-Allmaras all-zero turbulent-kinetic-energy path with a Yoshizawa wall-distance diagnostic.
- [x] `cfd-3d`: Reject uninitialized k-epsilon turbulence state instead of synthesizing zero viscosity, TKE, or dissipation fields.
- [x] `cfd-schematics`: Optimize memory layout with indexed node-layout caches and zero-copy parallel-group lookup.

## Rigor & Correctness
- [x] Review all numerical bounds and geometry assumptions in `cfd-schematics`.
