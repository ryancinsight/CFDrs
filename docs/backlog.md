# CFD Suite Backlog

## Sprint 1.96.167: cfd-math native sparse-LU result ownership
**Status**: Completed
**Owner**: Codex `/root`
**Change Class**: [minor] provider API + [patch] consumer allocation removal
**Start Date**: July 22, 2026

### Scope
- Consume `leto_ops::SparseLuSolver::solve_view` from the direct solver.
- Remove the consumer-owned RHS `Vec` staging and solution `Vec` to `Array1`
  copy while preserving the existing fallback and finiteness contracts.
- Synchronize the direct-solver tests, Rustdoc, changelog, and gap evidence.

### Non-Goals
- No change to dense-backed LU arithmetic, sparse storage, pivot policy,
  fallback thresholds, or the legacy provider slice API.
- No broad solver-family migration in this increment.

### Acceptance and Verification
- `DirectSparseSolver::solve` passes `rhs.view()` and returns the provider's
  owned `Array1` result on the primary path.
- Existing small-system, singular-input, dense-fallback, and generic scalar
  value tests pass through the configured native test runner.
- Format, warning-denied check/Clippy, focused Nextest, doctest, Rustdoc, and
  public provider API SemVer checks pass on the exact delivered revisions.

### Claimed Files
`crates/cfd-math/src/linear_solver/direct_solver.rs`, its focused tests, and
the corresponding `docs/{backlog,checklist,gap_audit}.md` plus `CHANGELOG.md`.

### Current Evidence
- Provider `leto-ops` check, warning-denied all-target Clippy, sparse Nextest
  (29/29), doctests (8/8), and Rustdoc pass against the local provider source.
- Consumer `cfd-math` check, lib Clippy, direct-solver Nextest (4/4), doctest,
  and Rustdoc pass. The package fmt check reports six pre-existing import-order
  diffs outside the claimed file; the touched direct-solver file passes
  standalone rustfmt.
- Provider public-surface SemVer classification passes 196/196 checks with
  57 skips, and the consumer pin is updated to merged Leto commit
  `b24fc860864abad84af3118aa2bb27c32bb81265`.

## Sprint 1.96.166: cfd-math IncompleteCholesky Leto CSR
**Status**: Completed
**Start Date**: July 5, 2026

### Sprint Objectives
- Move IncompleteCholesky construction and factor storage to
  `leto_ops::CsrMatrix`.
- Remove direct nalgebra sparse row-offset/get-entry use from the Cholesky
  preconditioner boundary.

### Verification
- cfd-math fmt check
- cfd-math lib check
- cfd-math all-target check
- lib/tests and all-target clippy
- cholesky-filter nextest (5/5 passed)
- preconditioner-filter nextest (76/76 passed)
- targeted Cholesky residue scan

### Outcome
- `cholesky.rs` no longer imports nalgebra sparse or uses direct nalgebra CSR
  row access/construction APIs.
- IncompleteCholesky construction is now Leto CSR-backed; `Preconditioner::
  apply_to` remains Leto `Array1`.
- Residual: Schwarz, direct solver, remaining transitional solver fixtures,
  and the shared solver sparse matrix boundary still expose nalgebra-sparse.

## Sprint 1.96.165: cfd-math ILU Leto CSR
**Status**: Completed
**Start Date**: July 5, 2026

### Sprint Objectives
- Move ILU(0), ILU(k), and triangular solve factor storage to
  `leto_ops::CsrMatrix`.
- Keep remaining shared solver/Schwarz sparse boundaries explicit as
  conversion points until those families move.

### Verification
- cfd-math fmt check
- cfd-math lib check
- cfd-math all-target check
- lib/tests and all-target clippy
- ilu-filter nextest (21/21 passed)
- preconditioner-filter nextest (76/76 passed)
- linear_solver::tests nextest (53/53 passed)
- targeted ILU residue scan

### Outcome
- `preconditioners/ilu` no longer imports nalgebra sparse or uses direct
  nalgebra CSR row access/construction APIs.
- ILU construction is now Leto CSR-backed; `Preconditioner::apply_to` remains
  Leto `Array1`.
- Residual: Schwarz, direct solver, remaining transitional solver fixtures,
  and the shared solver sparse matrix boundary still expose nalgebra-sparse.

## Sprint 1.96.164: cfd-math SSOR Leto CSR
**Status**: Completed
**Start Date**: July 5, 2026

### Sprint Objectives
- Move SSOR preconditioner construction to `leto_ops::CsrMatrix`.
- Keep source preconditioner tests explicit about the remaining legacy solver
  matrix boundary by converting fixtures once before SSOR construction.

### Verification
- cfd-math fmt check
- cfd-math lib check
- cfd-math all-target check
- lib/tests and all-target clippy
- ssor-filter nextest (5/5 passed)
- preconditioner-filter nextest (76/76 passed)
- targeted SSOR residue scan

### Outcome
- `ssor.rs` no longer imports nalgebra sparse or uses direct nalgebra CSR row
  access/construction APIs.
- SSOR construction is now Leto CSR-backed; `Preconditioner::apply_to`
  remains Leto `Array1`.
- Residual: Schwarz, direct solver, integration-test fixtures,
  and the shared solver sparse matrix boundary still expose nalgebra-sparse.

## Sprint 1.96.163: cfd-math Basic Preconditioner Leto CSR
**Status**: Completed
**Start Date**: July 5, 2026

### Sprint Objectives
- Move Jacobi and SOR basic preconditioner construction to
  `leto_ops::CsrMatrix`.
- Keep solver tests explicit about the remaining legacy solver matrix boundary
  by converting fixtures once before basic preconditioner construction.

### Verification
- cfd-math fmt check
- cfd-math lib check
- core solver test check
- cfd-math all-target check
- lib/tests, core solver, and all-target clippy
- linear_solver::tests nextest (53/53 passed)
- core_solver_tests nextest (4/4 passed)
- preconditioner-filter nextest (76/76 passed)
- targeted basic-preconditioner residue scan

### Outcome
- `basic.rs` no longer imports nalgebra sparse or uses direct nalgebra CSR row
  access/construction APIs.
- Jacobi/SOR construction is now Leto CSR-backed; `Preconditioner::apply_to`
  remains Leto `Array1`.
- Residual: Schwarz, direct solver, and the shared solver sparse matrix
  boundary still expose nalgebra-sparse.

## Sprint 1.96.162: cfd-math AMG/Coarsening Leto CSR
**Status**: Completed
**Start Date**: July 5, 2026

### Sprint Objectives
- Move the AMG/coarsening sparse provider boundary to `leto_ops::CsrMatrix`.
- Move coarsening/algebraic-distance benchmarks and AMG tests to Leto CSR
  construction where they exercise AMG.
- Keep the broader Krylov solver sparse boundary explicit as the next
  nalgebra-sparse provider gap.

### Verification
- cfd-math fmt check
- cfd-math lib check
- focused AMG/coarsening test and bench checks
- focused clippy for cfd-math lib, touched tests, and touched benches
- AMG integration nextest (5/5 passed)
- AMG-filter nextest (6/6 passed)
- multigrid::coarsening nextest (10/10 passed)
- targeted AMG/coarsening sparse residue scan

### Outcome
- Multigrid AMG now uses Leto CSR for setup, coarsening, interpolation,
  smoothers, cycles, sparse products, transpose, row access, and SpMV.
- `coarsening_bench.rs`, `algebraic_distance_bench.rs`,
  `amg_coarsening_tests.rs`, and AMG integration preconditioner construction
  now use Leto CSR.
- Residual: broader cfd-math solver/direct/preconditioner sparse matrix
  surfaces still expose nalgebra-sparse; `LinearSolverChain` converts once
  before AMG construction.

## Sprint 1.96.161: cfd-math Leto CSR Benchmarks
**Status**: Completed
**Start Date**: July 5, 2026

### Sprint Objectives
- Move SpMV/CG benchmark CSR construction to `leto_ops::CsrMatrix`.
- Measure direct Leto CSR `LinearOperator::apply` in the SpMV benchmark.
- Leave AMG coarsening benchmark holdouts tied to the production API boundary.

### Verification
- cfd-math fmt check
- focused `spmv_bench`, `cg_bench`, and `math_benchmarks` checks
- focused clippy for all three migrated benches
- cfd-math all-target check
- cfd-math all-target clippy
- sparse-filter nextest (19/19 passed)
- targeted migrated-benchmark residue scan

### Outcome
- `spmv_bench.rs`, `cg_bench.rs`, and the CG section of
  `math_benchmarks.rs` now construct Leto CSR matrices directly.
- SpMV benchmark dispatches through direct Leto CSR `LinearOperator::apply`.
- Residual after Sprint 1.96.162: the coarsening/algebraic-distance AMG
  boundary moved to Leto CSR; the broader solver/direct/preconditioner sparse
  matrix surface still exposes nalgebra-sparse.

## Sprint 1.96.160: cfd-math Leto CSR LinearOperator
**Status**: Completed
**Start Date**: July 5, 2026

### Sprint Objectives
- Let cfd-math iterative solvers consume `leto_ops::CsrMatrix<T>` directly.
- Share one Leto CSR SpMV helper between direct Leto CSR inputs and the
  remaining nalgebra CSR conversion path.
- Migrate the simple GMRES integration test off nalgebra sparse/vector
  construction.

### Verification
- cfd-math fmt check
- simple_gmres test check
- focused simple_gmres nextest (3/3 passed)
- focused simple_gmres clippy
- cfd-math lib check
- cfd-math all-target check
- cfd-math all-target clippy
- sparse-filter nextest (19/19 passed)
- gmres-filter nextest (21/21 passed)
- targeted simple GMRES residue scan

### Outcome
- `leto_ops::CsrMatrix<T>` now implements `LinearOperator<T>`.
- `tests/simple_gmres_tests.rs` now proves GMRES over Leto CSR matrices
  directly.
- Residual: broader cfd-math sparse builders, preconditioners, AMG, direct
  solver, tests, and benches still expose the nalgebra-sparse matrix boundary.

## Sprint 1.96.159: cfd-math Storage-Slice Closure
**Status**: Completed
**Start Date**: July 5, 2026

### Sprint Objectives
- Remove nonlinear mutable `StorageMut` slice helpers.
- Remove multigrid interpolation/smoother storage-slice residue.
- Verify cfd-math source/tests no longer contain the targeted Leto
  storage-slice patterns.

### Verification
- cfd-math fmt check
- cfd-math lib check
- focused cfd-math nonlinear_solver/multigrid nextest (46/46 passed)
- cfd-math lib/tests clippy
- cfd-math all-target check
- cfd-math all-target clippy
- cfd-math source/test storage-slice residue scan

### Outcome
- Nonlinear pivoting and multigrid validation/tests now use direct Leto array
  indexing.
- cfd-math source/tests no longer contain the targeted Leto storage-slice
  residue.
- Residual: remaining cfd-math Atlas-provider work is nalgebra/nalgebra-sparse
  replacement and other provider boundaries.

## Sprint 1.96.158: cfd-math Sparse/Basic Leto Array1
**Status**: Completed
**Start Date**: July 5, 2026

### Sprint Objectives
- Remove `leto::Storage` / `.storage().as_slice()` from sparse operations.
- Remove output `as_slice_mut()` contiguity assumptions from SpMV.
- Remove Jacobi preconditioner diagonal storage-slice iteration.
- Keep sparse operation ownership delegated to Leto CSR operations.

### Verification
- cfd-math fmt check
- cfd-math lib check
- focused cfd-math sparse/preconditioner nextest (95/95 passed)
- cfd-math lib/tests clippy
- cfd-math all-target check
- cfd-math all-target clippy
- targeted sparse/basic residue scans

### Outcome
- `sparse/operations.rs` now stages vector data through direct `Array1`
  indexing before delegating to Leto CSR operations.
- `preconditioners/basic.rs` no longer requires Leto storage slices for Jacobi
  diagonal construction.
- Residual: remaining cfd-math storage-slice owners are nonlinear mutable
  dense-workspace helpers and multigrid interpolation/smoother internals.

## Sprint 1.96.157: cfd-math GPU Operator Leto Array1
**Status**: Completed
**Start Date**: July 5, 2026

### Sprint Objectives
- Remove `leto::Storage` / `.storage().as_slice()` from the GPU linear
  operator.
- Remove output `as_slice_mut()` contiguity assumptions from the GPU operator.
- Preserve the existing Hephaestus-backed `cfd-core` GPU execution path.

### Verification
- cfd-math fmt check
- cfd-math GPU-feature check
- focused cfd-math GPU-feature linear_solver::operators nextest (5/5 passed)
- cfd-math GPU-feature lib clippy
- cfd-math GPU-feature all-target check
- cfd-math GPU-feature all-target clippy
- targeted GPU operator residue scan

### Outcome
- `operators/gpu.rs` now stages GPU upload/readback through direct `Array1`
  indexing and typed dimension checks.
- The GPU operator no longer requires Leto storage slices.
- Residual: remaining cfd-math storage-slice owners are sparse operations and
  multigrid internals; broader raw WGPU provider ownership remains outside
  this operator.

## Sprint 1.96.156: cfd-math Finite-Difference Operators Leto Array1
**Status**: Completed
**Start Date**: July 5, 2026

### Sprint Objectives
- Remove `leto::Storage` / `.storage().as_slice()` from CPU finite-difference
  linear operators.
- Remove output `as_slice_mut()` contiguity assumptions from those operators.
- Keep this slice scoped to CPU finite-difference operators; leave GPU
  provider work to a Hephaestus-specific increment.

### Verification
- cfd-math fmt check
- cfd-math lib check
- focused cfd-math linear_solver::operators nextest (5/5 passed)
- cfd-math lib/tests clippy
- cfd-math all-target check
- cfd-math all-target clippy
- targeted finite-difference operator residue scan

### Outcome
- `operators/poisson.rs` and `operators/momentum.rs` now use direct
  `Array1` indexing for input and output stencil access.
- Migrated CPU operators no longer require Leto storage slices.
- Residual: remaining storage-slice owners are sparse operations, GPU
  operator, and multigrid internals.

## Sprint 1.96.155: cfd-math Nonlinear Linalg Leto Array1
**Status**: Completed
**Start Date**: July 5, 2026

### Sprint Objectives
- Remove immutable `leto::Storage` / `.storage().as_slice()` use from
  `src/nonlinear_solver/linalg.rs`.
- Keep nonlinear solver vector algebra on direct `Array1` indexing.
- Update Anderson acceleration to index the least-squares coefficient vector
  directly.

### Verification
- cfd-math fmt check
- cfd-math lib check
- focused cfd-math nonlinear_solver nextest (9/9 passed)
- cfd-math lib/tests clippy
- cfd-math all-target check
- cfd-math all-target clippy
- targeted nonlinear linalg/Anderson residue scan

### Outcome
- `src/nonlinear_solver/linalg.rs` no longer imports immutable
  `leto::Storage` or exposes an immutable `vector_slice` helper.
- Nonlinear vector helpers read `Array1` values through indexing.
- Residual: mutable dense-workspace helpers still use `StorageMut`; remaining
  immutable storage-slice owners are sparse operations, linear operators, GPU
  operator, and multigrid internals.

## Sprint 1.96.154: cfd-math Production SIMD Vector Leto Array1
**Status**: Completed
**Start Date**: July 5, 2026

### Sprint Objectives
- Remove `leto::Storage` / `.storage().as_slice()` from `src/simd/vector.rs`.
- Keep Moirai `Adaptive` map/reduce dispatch in the SIMD vector operations.
- Verify the production SIMD vector extension and broader SIMD filter.

### Verification
- cfd-math fmt check
- cfd-math lib check
- focused cfd-math simd::vector nextest (1/1 passed)
- cfd-math lib/tests clippy
- cfd-math all-target check
- cfd-math all-target clippy
- broader cfd-math SIMD-filter nextest (26/26 passed)
- targeted source-residue scan

### Outcome
- `src/simd/vector.rs` no longer imports `leto::Storage` or borrows Leto raw
  storage slices.
- SIMD vector operations now stay on the `Array1` indexing surface while
  preserving Moirai `Adaptive` dispatch.
- Residual: source-level storage-slice owners remain in nonlinear linalg,
  sparse operations, linear operators, GPU operator, and multigrid internals.

## Sprint 1.96.153: cfd-math SIMD Integration Test Leto Array1
**Status**: Completed
**Start Date**: July 5, 2026

### Sprint Objectives
- Remove the remaining `leto::Storage`/`.storage().as_slice()` conversion from
  `tests/simd_tests.rs`.
- Keep the existing SIMD slice API assertions value-semantic while sourcing
  `spmv` results from direct Leto array indexing.
- Confirm the cfd-math integration-test vector bridge scan is clean.

### Verification
- cfd-math fmt check
- simd_tests check
- focused cfd-math simd_tests nextest (12/12 passed)
- simd_tests clippy
- cfd-math all-target check
- cfd-math all-target clippy
- broader cfd-math SIMD-filter nextest (26/26 passed)
- targeted cfd-math integration-test provider-residue scans

### Outcome
- `tests/simd_tests.rs` no longer imports `leto::Storage` or converts the Leto
  `spmv` result through `.storage().as_slice()`.
- `crates/cfd-math/tests` now has no `DVector`, nalgebra vector import, local
  preconditioner bridge, `Storage`, or storage-slice conversion residue.
- Residual: remaining provider work is outside the integration-test vector
  bridge layer, including `nalgebra_sparse::CsrMatrix`, dense nalgebra test
  oracles, and source-level Leto storage-slice internals.

## Sprint 1.96.152: cfd-math AMG Integration Test Leto Array1
**Status**: Completed
**Start Date**: July 5, 2026

### Sprint Objectives
- Move AMG exact solution, RHS, solver output, cycle output, and two-grid work
  vectors from nalgebra `DVector` bridge helpers to direct Leto arrays.
- Exercise `Preconditioner::apply_to` directly at the Leto RHS/output
  boundary.
- Keep the current sparse matrix and dense eigenvalue oracle boundaries
  explicit.

### Verification
- cfd-math fmt check
- amg_integration_test check
- focused cfd-math amg_integration_test nextest (5/5 passed)
- amg_integration_test clippy
- cfd-math all-target check
- cfd-math all-target clippy
- broader cfd-math AMG-filter nextest (6/6 passed)
- targeted amg_integration_test provider-residue scan

### Outcome
- `tests/amg_integration_test.rs` no longer imports nalgebra `DVector` or uses
  a local preconditioner bridge helper.
- AMG integration tests now exercise the same Leto RHS/output boundary as the
  production solver and preconditioner traits.
- Residual: this test still uses `nalgebra_sparse::CsrMatrix` for sparse
  storage and nalgebra `DMatrix`/`SymmetricEigen` for the dense energy-norm
  oracle. Remaining cfd-math integration storage-slice residue is in
  `tests/simd_tests.rs`.

## Sprint 1.96.151: cfd-math Preconditioner Edge-Case Tests Leto Array1
**Status**: Completed
**Start Date**: July 5, 2026

### Sprint Objectives
- Move ILU(0), ILU(k), repeated-application, extreme-value, and
  sparsity-preservation preconditioner edge-case tests from nalgebra `DVector`
  bridge helpers to direct Leto arrays.
- Exercise `Preconditioner::apply_to` directly at the Leto RHS/output
  boundary.
- Keep the current sparse matrix storage boundary explicit.

### Verification
- cfd-math fmt check
- preconditioner_edge_cases check
- focused cfd-math preconditioner_edge_cases nextest (6/6 passed)
- preconditioner_edge_cases clippy
- cfd-math all-target check
- cfd-math all-target clippy
- broader cfd-math preconditioner nextest (76/76 passed)
- targeted preconditioner_edge_cases provider-residue scan

### Outcome
- `tests/preconditioner_edge_cases.rs` no longer imports nalgebra `DVector` or
  uses a local preconditioner bridge helper.
- Preconditioner edge-case tests now exercise the same Leto RHS/output
  boundary as the production preconditioner trait.
- Residual: this test still uses `nalgebra_sparse::CsrMatrix` for matrix
  storage; integration-test vector bridge residue was closed by Sprint
  1.96.153.

## Sprint 1.96.150: cfd-math Linear-Solver Test Module Leto Array1
**Status**: Completed
**Start Date**: July 5, 2026

### Sprint Objectives
- Move `src/linear_solver/tests` from nalgebra `DVector` bridge macros to
  direct Leto arrays.
- Verify residuals through the Leto SpMV/helper path.
- Route scalar constants/tolerances in the touched solver/sparse cone through
  Eunomia.

### Verification
- cfd-math fmt check
- cfd-math lib check
- cfd-math all-target check
- focused cfd-math linear_solver::tests nextest (53/53 passed)
- broader cfd-math linear_solver nextest (176/176 passed)
- cfd-math lib/tests clippy
- cfd-math all-target clippy
- targeted source-test provider-residue scan
- targeted old scalar helper residue scan

### Outcome
- `src/linear_solver/tests` no longer imports nalgebra `DVector` or uses local
  solve/preconditioner bridge macros.
- Source linear-solver tests now exercise the same Leto RHS/solution boundary
  as the production solver traits.
- Residual: integration-test vector bridge residue was closed by Sprint
  1.96.153; production sparse/scalar boundaries still include transitional
  nalgebra providers.

## Sprint 1.96.149: cfd-math Core Solver Tests Leto Array1
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move BiCGSTAB, GMRES, preconditioner integration, and condition-number
  robustness tests from nalgebra `DVector` bridge helpers to direct Leto
  arrays.
- Verify residuals through the Leto SpMV path instead of nalgebra
  matrix-vector multiplication.
- Keep the current sparse matrix storage boundary explicit.

### Verification
- cfd-math fmt check
- core_solver_tests check
- focused cfd-math core_solver_tests nextest (4/4 passed)
- core_solver_tests clippy
- cfd-math all-target check
- cfd-math all-target clippy
- targeted core solver test provider-residue scan

### Outcome
- `tests/core_solver_tests.rs` no longer imports nalgebra `DVector` or uses a
  local solve/preconditioner bridge helper.
- Core solver validation tests now exercise the same Leto RHS/solution
  boundary as the production solver traits.
- Residual: this test still uses `nalgebra_sparse::CsrMatrix`/`CooMatrix` for
  matrix storage, and broader cfd-math test diagnostics still contain nalgebra
  `DVector` bridges.

## Sprint 1.96.148: cfd-math Simple GMRES Tests Leto Array1
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move basic, restarted, and preconditioned GMRES integration tests from
  nalgebra `DVector` bridge macros to direct Leto arrays.
- Verify residuals through the Leto SpMV path instead of nalgebra matrix-vector
  multiplication.
- Keep the current sparse matrix storage boundary explicit.

### Verification
- cfd-math fmt check
- simple_gmres test check
- focused cfd-math simple_gmres nextest (3/3 passed)
- simple_gmres clippy
- cfd-math all-target check
- cfd-math all-target clippy
- touched-file `git diff --check`
- targeted simple GMRES test provider-residue scan

### Outcome
- `tests/simple_gmres_tests.rs` no longer imports nalgebra `DVector` or uses a
  local solve bridge macro.
- Simple GMRES tests now exercise the same Leto RHS/solution boundary as the
  production solver traits.
- Residual: this test still uses `nalgebra_sparse::CsrMatrix`/`CooMatrix` for
  matrix storage, and broader cfd-math test diagnostics still contain nalgebra
  `DVector` bridges.

## Sprint 1.96.147: cfd-math Matrix-Free Tests Leto Array1
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move matrix-free CG/GMRES identity tests from nalgebra `DVector` bridge
  macros to direct Leto arrays.
- Move scaled-operator matrix-free integration to direct Leto arrays.
- Assert the exact typed operator-size mismatch error at the Leto boundary.

### Verification
- cfd-math fmt check
- cfd-math all-target check
- cfd-math all-target clippy
- focused cfd-math matrix-free nextest (4/4 passed)
- touched-file `git diff --check`
- targeted matrix-free test provider-residue scan

### Outcome
- `linear_solver/matrix_free/tests.rs` no longer imports nalgebra `DVector` or
  uses a local solve bridge macro.
- Matrix-free solver tests now exercise the same Leto RHS/solution boundary as
  the production solver traits.
- Residual: broader cfd-math linear-solver integration/adversarial/core/
  preconditioner test diagnostics still contain nalgebra `DVector` bridges,
  and production sparse/scalar provider boundaries still include
  `nalgebra_sparse::CsrMatrix` and transitional `nalgebra::RealField`.

## Sprint 1.96.146: cfd-math BiCGSTAB Leto Array1
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move BiCGSTAB residual/search/operator workspaces from nalgebra `DVector` to
  Leto `Array1`.
- Run preconditioned and unpreconditioned BiCGSTAB at the Leto array boundary.
- Remove BiCGSTAB calls to the legacy linear-solver vector bridge helpers.
- Keep `LinearSolverChain` final BiCGSTAB fallback on Leto arrays.
- Consolidate CG/BiCGSTAB Leto vector operations into one helper module.

### Verification
- cfd-math fmt check
- cfd-math all-target check
- cfd-math all-target clippy
- focused cfd-math BiCGSTAB nextest (24/24 passed)
- broader cfd-math linear-solver nextest (176/176 passed)
- cfd-math AMG integration nextest (5/5 passed)
- touched-file `git diff --check`
- targeted BiCGSTAB/chain/traits provider-residue scans

### Outcome
- `linear_solver/bicgstab/mod.rs` no longer owns nalgebra `DVector`
  workspaces or legacy bridge calls.
- `linear_solver/chain.rs` no longer converts the final BiCGSTAB fallback
  through nalgebra vectors.
- `linear_solver/traits.rs` no longer contains obsolete legacy vector bridge
  helpers.
- `linear_solver/array_ops.rs` is the shared Leto vector helper SSOT for CG
  and BiCGSTAB.
- Residual: cfd-math still carries `nalgebra_sparse::CsrMatrix`,
  transitional `nalgebra::RealField` scalar bounds, and nalgebra `DVector` in
  remaining matrix-free/preconditioner/integration test diagnostics.

## Sprint 1.96.145: cfd-math Conjugate Gradient Leto Array1
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move CG residual/direction/operator workspaces from nalgebra `DVector` to
  Leto `Array1`.
- Run preconditioned and unpreconditioned CG at the Leto array boundary.
- Remove CG calls to the legacy linear-solver vector bridge helpers.
- Move CG benchmark call sites to Leto vectors.

### Verification
- cfd-math fmt check
- cfd-math all-target check
- cfd-math all-target clippy
- focused cfd-math conjugate nextest (13/13 passed)
- broader cfd-math linear-solver nextest (176/176 passed)
- targeted CG provider-residue scan

### Outcome
- `linear_solver/conjugate_gradient/mod.rs` no longer owns nalgebra `DVector`
  workspaces or legacy bridge calls.
- `benches/cg_bench.rs` and `benches/math_benchmarks.rs` use Leto vectors at
  the CG API.
- CG value/error behavior is covered by exact tests.
- Residual after the BiCGSTAB follow-up: the shared linear-solver trait family
  still carries the transitional `nalgebra::RealField` scalar bound and sparse
  storage remains on `nalgebra_sparse::CsrMatrix`.

## Sprint 1.96.144: cfd-math Schwarz Preconditioner Leto Array1
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move additive/multiplicative Schwarz local apply methods from nalgebra
  `DVector` to Leto `Array1`.
- Extract local RHS buffers directly into Leto arrays.
- Apply Schwarz at the Leto array boundary without global `DVector` bridges.
- Reject mismatched residual/output lengths with typed configuration errors.

### Verification
- cfd-math fmt check
- cfd-math all-target check
- cfd-math all-target clippy
- focused cfd-math Schwarz nextest (3/3 passed)
- broader cfd-math preconditioner nextest (76/76 passed)
- targeted Schwarz provider-residue scan

### Outcome
- `linear_solver/preconditioners/schwarz.rs` no longer owns `DVector` local
  apply signatures, local RHS workspaces, or default apply conversion bridges.
- Schwarz Leto boundary behavior is covered by exact value and typed-error
  tests.
- Residual: Schwarz still stores and constructs local sparse matrices through
  the shared `nalgebra_sparse::CsrMatrix` boundary and remains constrained by
  the global `Preconditioner<T>` nalgebra scalar bound.

## Sprint 1.96.143: cfd-math ILU Triangular Solve Leto Array1
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move ILU forward/backward substitution workspaces from nalgebra `DVector` to
  Leto `Array1`.
- Apply ILU at the Leto array boundary without residual/intermediate/result
  `DVector` bridges.
- Reject mismatched residual/output lengths with typed configuration errors.
- Route U-solve diagonal identity through Eunomia.

### Verification
- cfd-math fmt check
- cfd-math all-target check
- cfd-math all-target clippy
- focused cfd-math ILU nextest (21/21 passed)
- broader cfd-math preconditioner nextest (74/74 passed)
- targeted ILU provider-residue scan

### Outcome
- `linear_solver/preconditioners/ilu/{types,triangular}.rs` no longer owns
  `DVector` triangular-solve workspaces or conversion bridges.
- ILU mismatch behavior is covered by exact typed-error tests.
- Residual: IncompleteLU still stores the shared
  `nalgebra_sparse::CsrMatrix` LU factor boundary and remains constrained by
  the global `Preconditioner<T>` nalgebra scalar bound.

## Sprint 1.96.142: cfd-math Deflation Preconditioner Leto Array1
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move deflation eigenvector state from nalgebra `DVector` to Leto `Array1`.
- Apply deflation projection corrections without nalgebra work vectors.
- Reject mismatched vector lengths with typed configuration errors.
- Reject zero eigenvalues before projection division.

### Verification
- cfd-math fmt check
- cfd-math all-target check
- cfd-math all-target clippy
- focused cfd-math deflation nextest (3/3 passed)
- broader cfd-math preconditioner nextest (73/73 passed)
- targeted deflation provider-residue scan

### Outcome
- `linear_solver/preconditioners/deflation.rs` no longer owns `DVector`
  eigenvector storage or deflated workspaces.
- Deflation projection, mismatch, and zero-eigenvalue behavior are covered by
  value-semantic tests.
- Residual: Deflation still wraps the base preconditioner behind the existing
  `Box<dyn Preconditioner<T>>`, and the global `Preconditioner<T>` trait still
  uses the transitional nalgebra `RealField` scalar bound.

## Sprint 1.96.141: cfd-math Basic Preconditioners Leto Array1
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move Identity/Jacobi/SOR basic preconditioner residual/output paths from
  nalgebra vector state to Leto `Array1`.
- Store Jacobi inverse diagonal state as `leto::Array1`.
- Reject mismatched residual/output lengths with typed configuration errors.
- Route scalar identities and diagonal tolerance through Eunomia.

### Verification
- cfd-math fmt check
- cfd-math all-target check
- cfd-math all-target clippy
- focused cfd-math basic mismatch nextest (1/1 passed)
- broader cfd-math preconditioner nextest (70/70 passed)
- targeted basic-preconditioner provider-residue scan

### Outcome
- `linear_solver/preconditioners/basic.rs` no longer owns nalgebra `DVector`
  storage or apply bridges for Identity/Jacobi/SOR vector paths.
- Basic preconditioner mismatch behavior is covered by exact typed-error tests.
- Residual: Jacobi and SOR still store the shared
  `nalgebra_sparse::CsrMatrix` boundary and remain constrained by the global
  `Preconditioner<T>` nalgebra scalar bound.

## Sprint 1.96.140: cfd-math IncompleteCholesky Leto Array1
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move IncompleteCholesky residual, intermediate, and solution workspaces from
  nalgebra `DVector` to Leto `Array1`.
- Apply IncompleteCholesky at the Leto array boundary without
  residual/result `DVector` bridges.
- Reject mismatched residual/output lengths with typed configuration errors.
- Route IC(0) square-root dispatch through Eunomia.

### Verification
- cfd-math fmt check
- cfd-math all-target check
- cfd-math all-target clippy
- focused cfd-math Cholesky nextest (5/5 passed)
- targeted Cholesky provider-residue scan

### Outcome
- `linear_solver/preconditioners/cholesky.rs` no longer owns `DVector`
  substitution workspaces or conversion bridges.
- Cholesky direct substitution and trait `apply_to` paths use `leto::Array1`.
- Follow-up: Sprint 1.96.166 moved IncompleteCholesky factor storage to Leto
  CSR. Residual sparse-provider work is now Schwarz, direct solver, and the
  shared solver matrix boundary.

## Sprint 1.96.139: cfd-math SSOR Preconditioner Leto Array1
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move SSOR forward/backward sweep residual and solution workspaces from
  nalgebra `DVector` to Leto `Array1`.
- Apply SSOR at the Leto array boundary without residual/result `DVector`
  bridges.
- Reject mismatched residual/output lengths with typed configuration errors.

### Verification
- cfd-math fmt check
- cfd-math all-target check
- cfd-math all-target clippy
- focused cfd-math SSOR nextest (5/5 passed)
- targeted SSOR provider-residue scan

### Outcome
- `linear_solver/preconditioners/ssor.rs` no longer owns `DVector` sweep
  workspaces or conversion bridges.
- SSOR direct sweeps and trait `apply_to` paths use `leto::Array1`.
- Residual: SSOR still stores the shared `nalgebra_sparse::CsrMatrix`
  boundary and remains constrained by the global `Preconditioner<T>` nalgebra
  scalar bound.

## Sprint 1.96.138: cfd-math Block/SIMPLE Preconditioner Leto Array1
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move block/SIMPLE preconditioner diagonal and Schur vector state from
  nalgebra `DVector` to Leto `Array1`.
- Apply block/SIMPLE preconditioners at the Leto array boundary without
  residual/result `DVector` bridges.
- Replace silent mismatched-vector cloning with typed configuration errors.

### Verification
- cfd-math fmt check
- cfd-math all-target check
- cfd-math all-target clippy
- focused cfd-math block-preconditioner nextest (4/4 passed)
- targeted block-preconditioner provider-residue scan

### Outcome
- `linear_solver/block_preconditioner.rs` no longer owns `DVector` storage or
  conversion bridges.
- Direct `apply` and trait `apply_to` paths use `leto::Array1`.
- Residual: the global `Preconditioner` trait still carries the transitional
  `nalgebra::RealField` scalar bound, and sparse storage still aliases
  `nalgebra_sparse::CsrMatrix`.

## Sprint 1.96.137: cfd-math GMRES Leto Workspace
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move GMRES Arnoldi basis-column extraction, work vectors, preconditioner work,
  and residual checks from nalgebra `DVector` to Leto `Array1`.
- Move GMRES direct solve methods to Leto RHS/solution arrays.
- Keep `LinearSolverChain` GMRES tiers on Leto arrays instead of converting to
  nalgebra before each tier.

### Verification
- cfd-math lib/tests/all-targets check
- cfd-math fmt check
- cfd-math all-target clippy
- focused cfd-math GMRES nextest (21/21 passed)
- cfd-2d no-default lib check
- focused cfd-2d momentum nextest (53/53 passed)
- targeted GMRES `DVector`/legacy bridge residue scan

### Outcome
- `linear_solver/gmres` no longer owns `DVector` workspaces or legacy
  operator/preconditioner bridge calls.
- `LinearSolverChain` uses Leto arrays through GMRES+AMG, GMRES+block,
  unpreconditioned GMRES, and GMRES+ILU tiers.
- Residual: cfd-math still has nalgebra scalar bounds, CG/BiCGSTAB
  workspaces, the final solver-chain BiCGSTAB bridge, and nalgebra-sparse
  storage pending later Leto/Eunomia slices.

## Sprint 1.96.136: cfd-2d Momentum Leto Array1
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move cfd-2d momentum RHS and solution buffers from nalgebra `DVector` to
  Leto `Array1`.
- Remove the obsolete local nalgebra solver bridge.
- Remove direct cfd-2d nalgebra/nalgebra-sparse source and manifest ownership.

### Verification
- cfd-2d fmt check
- cfd-2d no-default lib check
- cfd-2d no-default all-target clippy
- focused cfd-2d momentum nextest (53/53 passed)
- direct cfd-2d source/manifest nalgebra residue scan
- cfd-2d nalgebra/nalgebra-sparse cargo-tree audit

### Residual
- cfd-2d still resolves nalgebra/nalgebra-sparse transitively through upstream
  crates including cfd-1d, cfd-core, cfd-math, cfd-schematics, and Gaia.

---

## Sprint 1.96.135: cfd-2d Pressure-Velocity Leto Array1
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move cfd-2d pressure-velocity pressure-correction RHS and solution caches
  from nalgebra `DVector` to Leto `Array1`.
- Call Leto-native iterative and direct solver APIs directly from
  pressure-velocity dispatch.

### Verification
- cfd-2d fmt check
- cfd-2d no-default lib check
- cfd-2d no-default all-target clippy
- focused cfd-2d pressure-velocity nextest (16/16 passed)
- targeted pressure-velocity DVector/nalgebra/bridge residue scan

### Residual
- direct cfd-2d source/manifest nalgebra ownership is now removed; nalgebra
  remains transitive through upstream owners.

---

## Sprint 1.96.134: cfd-2d SIMPLE Leto Array1
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move cfd-2d SIMPLE pressure-correction RHS and `p_prime` buffers from
  nalgebra `DVector` to Leto `Array1`.
- Call the Leto-native iterative solver boundary directly for SIMPLE pressure
  correction.

### Verification
- cfd-2d fmt check
- cfd-2d no-default lib check
- cfd-2d no-default all-target clippy
- focused cfd-2d SIMPLE nextest (19/19 passed)
- targeted SIMPLE DVector/nalgebra/bridge residue scan

### Residual
- direct cfd-2d source/manifest nalgebra ownership is now removed; nalgebra
  remains transitive through upstream owners.

---

## Sprint 1.96.133: cfd-2d FDM Leto Array1
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move cfd-2d FDM RHS/result vectors from nalgebra `DVector` to Leto
  `Array1`.
- Update Poisson and advection-diffusion stencil assembly to mutate provider
  arrays directly.
- Update shared Gauss-Seidel solve output to return a Leto solution vector.

### Verification
- cfd-2d fmt check
- cfd-2d no-default lib check
- cfd-2d no-default all-target clippy
- focused cfd-2d FDM nextest (2/2 passed)
- targeted FDM DVector/nalgebra residue scan

### Residual
- direct cfd-2d source/manifest nalgebra ownership is now removed; nalgebra
  remains transitive through upstream owners.

---

## Sprint 1.96.132: cfd-2d Time Integration Leto Array1
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move cfd-2d `schemes::time` public state vectors from nalgebra `DVector` to
  Leto `Array1`.
- Update explicit, implicit, multistep, adaptive controller, adaptive
  integrator, and time tests to the same provider boundary.
- Use Leto-owned vector norm calculation for fixed-point convergence checks.

### Verification
- cfd-2d fmt check
- cfd-2d no-default lib check
- cfd-2d no-default all-target clippy
- focused cfd-2d time nextest (29/29 passed)
- targeted time-cone DVector/nalgebra residue scans

### Residual
- direct cfd-2d source/manifest nalgebra ownership is now removed; nalgebra
  remains transitive through upstream owners.

---

## Sprint 1.96.131: cfd-2d DMatrix Leto Array2
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move cfd-2d immersed-boundary force/velocity matrices from nalgebra
  `DMatrix` to `leto::Array2`.
- Move `schemes::Grid2D` storage and dependent scheme callers/tests to Leto
  shape/indexing.
- Update the `blood_venturi` example to use the same provider boundary.

### Verification
- cfd-2d fmt check
- cfd-2d no-default lib check/clippy
- cfd-2d no-default `blood_venturi` example check
- focused immersed-boundary/schemes/upwind/MUSCL nextest (60/60 passed)
- targeted DMatrix and nalgebra-style grid-access residue scans

### Residual
- cfd-2d still retains nalgebra in non-DMatrix `DVector`/sparse linear-system
  boundaries and other provider seams.

---

## Sprint 1.96.130: cfd-1d Vascular Eunomia Complex
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move cfd-1d vascular Bessel recurrence complex values from nalgebra to
  `eunomia::Complex`.
- Move Womersley analytical velocity, wall-shear, and flow-rate complex
  construction to Eunomia.
- Remove the vascular path's `nalgebra::ComplexField` dependency.

### Verification
- cfd-1d lib check/clippy
- focused Bessel/Womersley nextest (26/26 passed)
- cfd-1d fmt check
- targeted nalgebra-complex residue scan

### Residual
- cfd-1d still depends on nalgebra for network sparse/dense linear-system
  storage and the transitional scalar seam.

---

## Sprint 1.96.129: LinearOperator Trait Leto Vectors
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move `LinearOperator::apply` and `apply_transpose` vector surfaces from
  nalgebra `DVector` to Leto `Array1`.
- Keep nalgebra vector conversion confined behind current solver workspaces.
- Update cfd-math operator implementations and the Laplacian CPU benchmark to
  the Leto public API.

### Verification
- cfd-math all-target check/clippy
- focused cfd-math solver/operator nextest (80/80 passed)
- cfd-math fmt check
- targeted DVector operator-signature residue scan

### Residual
- Nalgebra sparse storage and some preconditioner internals still retain local
  `DVector`/`CsrMatrix` conversion bridges after the CG/BiCGSTAB follow-up.
- cfd-validation numerical solver result/error storage still uses `DVector`.
- Broad cfd-validation nextest still fails in the existing venturi
  cross-fidelity convergence tests.

---

## Sprint 1.96.128: Preconditioner Trait Leto Vectors
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move `Preconditioner::apply_to` residual/result surfaces from nalgebra
  `DVector` to Leto `Array1`.
- Keep nalgebra vector conversion confined behind current preconditioner
  matrix internals and solver workspaces.
- Update cfd-1d network preconditioning to the Leto public API.

### Verification
- cfd-math all-target check/clippy
- cfd-1d lib check/clippy
- focused cfd-math solver/preconditioner nextest (131/131 passed)
- cfd-math/cfd-1d fmt check
- targeted `apply_to` DVector-signature residue scan

### Residual
- `LinearOperator::apply`, nalgebra sparse preconditioner internals, and
  iterative solver workspaces still retain `DVector`/`CsrMatrix` conversion
  bridges.
- cfd-validation numerical solver result/error storage still uses `DVector`.
- Broad cfd-validation nextest still fails in the existing venturi
  cross-fidelity convergence tests.

---

## Sprint 1.96.127: Iterative Solver Trait Leto Vectors
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move `IterativeLinearSolver::solve` RHS/result surfaces from nalgebra
  `DVector` to Leto `Array1`.
- Keep nalgebra vector conversion confined behind current solver workspaces
  and local consumer bridges.
- Update cfd-1d, cfd-2d, and cfd-3d iterative solver callers to the Leto
  public API.

### Verification
- cfd-math fmt/check/clippy
- focused cfd-math solver nextest (61/61 passed)
- cfd-1d/cfd-2d/cfd-3d focused checks and clippy
- cfd-validation no-default all-target clippy
- targeted DVector call-site residue scan
- `git diff --check`

### Residual
- `LinearOperator::apply`, `Preconditioner::apply_to`, preconditioners, and
  internal iterative workspaces still expose nalgebra `DVector`/`CsrMatrix`.
- cfd-validation numerical solver result/error storage still uses `DVector`.
- Broad cfd-validation nextest still fails in the existing venturi
  cross-fidelity convergence tests.

---

## Sprint 1.96.126: cfd-math LinearSolver Trait Leto Vectors
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move `LinearSolver::solve_system` RHS/result surfaces from nalgebra
  `DVector` to Leto `Array1`.
- Route CG, BiCGSTAB, and GMRES implementations through the Leto public
  boundary while confining current nalgebra workspaces internally.
- Update cfd-validation numerical solver validation to call the Leto public
  API.

### Verification
- `cargo fmt -p cfd-math -p cfd-validation --check`
- `cargo check -p cfd-math --no-default-features --lib`
- `cargo check -p cfd-validation --no-default-features --lib`
- `cargo clippy -p cfd-math --no-default-features --all-targets -- -D warnings`
- `cargo clippy -p cfd-validation --no-default-features --lib -- -D warnings`
- `cargo clippy -p cfd-validation --no-default-features --all-targets -- -D warnings`
- `cargo clippy -p cfd-2d --no-default-features --all-targets -- -D warnings`
- `cargo nextest run -p cfd-math --no-default-features conjugate_gradient bicgstab gmres --status-level fail`
  (58/58 passed)
- targeted `solve_system` signature and DVector-residue scans
- `git diff --check`

### Residual
- `IterativeLinearSolver`, `LinearOperator`, `Preconditioner`,
  preconditioners, and internal iterative workspaces still expose nalgebra
  `DVector`/`CsrMatrix`.
- cfd-validation numerical solver result/error storage still uses `DVector`.
- Broad cfd-validation nextest still fails in the existing venturi
  cross-fidelity convergence tests.

---

## Sprint 1.96.125: cfd-validation Leto SpMV and Scalar Bounds
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move cfd-validation public SpMV benchmark callers to `leto::Array1`.
- Propagate the crate-local `ValidationScalar` seam through sparse
  linear-solver validation and 1D blood-flow literature validation.
- Clear the cfd-validation blocker in cfd-2d no-default all-target clippy.

### Verification
- `cargo fmt -p cfd-validation --check`
- `cargo check -p cfd-validation --no-default-features --lib`
- `cargo clippy -p cfd-validation --no-default-features --lib -- -D warnings`
- `cargo clippy -p cfd-validation --no-default-features --all-targets -- -D warnings`
- `cargo clippy -p cfd-2d --no-default-features --all-targets -- -D warnings`
- cfd-validation SpMV DVector-residue scan
- `cargo nextest run -p cfd-validation --no-default-features benchmark --status-level fail`
  (40/40 passed)
- `git diff --check`

### Residual
- Broad cfd-validation nextest still fails in two venturi cross-fidelity 2D
  fallback convergence tests unrelated to this SpMV/scalar-bound migration.

---

## Sprint 1.96.124: Solver Chain and FEM Leto Vectors
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move `LinearSolverChain` RHS/result surfaces from nalgebra `DVector` to
  `leto::Array1`.
- Route cfd-2d direct sparse fallback consumers through one Leto bridge.
- Route cfd-3d FEM sparse assembly and chain consumers through one Leto bridge.
- Propagate Leto real-scalar bounds through the cfd-2d/cfd-3d scalar seams and
  the cfd-1d network-solver edge needed by 2D coupling.

### Verification
- `cargo fmt -p cfd-math -p cfd-1d -p cfd-2d -p cfd-3d --check`
- `cargo check -p cfd-math --no-default-features --lib`
- `cargo check -p cfd-1d --no-default-features --lib`
- `cargo check -p cfd-2d --no-default-features --lib`
- `cargo check -p cfd-3d --no-default-features --lib`
- `cargo nextest run -p cfd-math --no-default-features chain direct_solver core_solver simple_gmres --status-level fail`
  (4/4 passed)
- `cargo clippy -p cfd-math --no-default-features --all-targets -- -D warnings`
- `cargo clippy -p cfd-2d --no-default-features --lib -- -D warnings`
- `cargo clippy -p cfd-3d --no-default-features --lib -- -D warnings`

### Residual
- `cfd-validation` still blocks cfd-2d all-target clippy because it has
  nalgebra-vector public SpMV calls and generic Leto scalar-bound gaps exposed
  by the migrated cfd-math/cfd-1d provider contracts.

---

## Sprint 1.96.123: cfd-math Direct Solver Leto Vectors
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move `DirectSparseSolver::solve` RHS and result vectors to `leto::Array1`.
- Remove the obsolete DVector dense-fallback wrapper.
- Keep `LinearSolverChain` as the only direct-tier bridge back to nalgebra
  while that public chain API remains unmigrated.
- Keep the broader iterative solver/preconditioner migration explicit as
  residual work.

### Verification
- `cargo fmt -p cfd-math --check`
- `cargo check -p cfd-math --no-default-features --lib`
- `cargo nextest run -p cfd-math --no-default-features direct_solver chain core_solver simple_gmres --status-level fail`
  (4/4 passed)
- `cargo clippy -p cfd-math --no-default-features --all-targets -- -D warnings`
- `cargo doc -p cfd-math --no-default-features --no-deps`
- `cargo test --doc -p cfd-math --no-default-features` (3 passed, 3 ignored)
- Targeted direct-solver DVector-signature residue scan
- `git diff --check`

---

## Sprint 1.96.122: cfd-math Sparse Builder Leto RHS
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move `SparseMatrixBuilder::build_with_rhs` from nalgebra `DVector` RHS
  mutation to `leto::Array1`.
- Keep Dirichlet column-elimination semantics value-identical.
- Remove the obsolete dummy nalgebra RHS from `SparseMatrixBuilder::build`.
- Update direct-solver and block-preconditioner call sites that build matrices
  through `build_with_rhs`.
- Keep the broader linear-solver trait boundary explicit as residual nalgebra
  `DVector` work.

### Verification
- `cargo fmt -p cfd-math --check`
- `cargo check -p cfd-math --no-default-features --lib`
- `cargo nextest run -p cfd-math --no-default-features sparse direct_solver block_preconditioner --status-level fail`
  (25/25 passed)
- `cargo clippy -p cfd-math --no-default-features --all-targets -- -D warnings`
- Targeted `build_with_rhs` residue scan
- `git diff --check`

---

## Sprint 1.96.121: cfd-math Public SpMV Leto Vectors
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move the then-public SpMV vector arguments to `leto::Array1`; the redundant
  parallel-named wrapper was removed in the 2026-07-10 closure.
- Keep the nalgebra `DVector` SpMV path private for the current
  `LinearOperator` trait implementation.
- Update sparse tests, GMRES/AMG integration tests, interpolation quality
  checks, and the SpMV benchmark to use Leto arrays at public SpMV call sites.
- Keep the public linear-solver/preconditioner trait boundary explicit as
  residual nalgebra `DVector` work.

### Verification
- `cargo fmt -p cfd-math --check`
- `cargo check -p cfd-math --no-default-features --lib`
- `cargo nextest run -p cfd-math --no-default-features sparse spmv interpolation amg simple_gmres core_solver --status-level fail`
  (40/40 passed)
- `cargo clippy -p cfd-math --no-default-features --all-targets -- -D warnings`
- Public SpMV signature scan
- `git diff --check` for the touched sparse/solver test/bench files

---

## Sprint 1.96.120: cfd-math SparseMatrixExt Leto Vectors
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move `SparseMatrixExt::diagonal` to return `leto::Array1`.
- Move `SparseMatrixExt::{set_diagonal,scale_rows,scale_columns}` to accept
  `leto::Array1`.
- Update sparse extension tests to use Leto vectors for row/column scaling.
- Keep the public linear-solver/preconditioner `DVector` boundary explicit as
  residual work.

### Verification
- `cargo fmt -p cfd-math --check`
- `cargo check -p cfd-math --no-default-features --lib`
- `cargo nextest run -p cfd-math --no-default-features sparse basic --status-level fail`
  (21/21 passed)
- `cargo clippy -p cfd-math --no-default-features --all-targets -- -D warnings`
- Targeted `SparseMatrixExt` DVector-signature residue scan
- `git diff --check` for the touched sparse/basic files

---

## Sprint 1.96.119: cfd-math Multigrid Smoother/Cycle Leto Vectors
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move `MultigridLevel`, `AMGHierarchy`, and `MultigridSmoother` vector
  contracts to `leto::Array1`.
- Move Jacobi, Gauss-Seidel, symmetric Gauss-Seidel, SSOR, and Chebyshev
  smoothers to Leto vectors.
- Add a Leto-array sparse SpMV bridge backed by `leto_ops::spmv_into`.
- Move V/W/F multigrid cycles and coarsest-cycle solves to Leto vectors.
- Keep the current public `Preconditioner`/linear-solver boundary explicit as
  residual nalgebra `DVector`/`RealField` work.

### Verification
- `cargo fmt -p cfd-math --check`
- `cargo check -p cfd-math --no-default-features --lib`
- `cargo nextest run -p cfd-math --no-default-features multigrid::cycles smoothers --status-level fail`
  (10/10 passed)
- `cargo clippy -p cfd-math --no-default-features --all-targets -- -D warnings`
- Clean smoother/cycle provider-residue scan
- `git diff --check` for the migrated sparse/multigrid files

---

## Sprint 1.96.118: cfd-math GMG Leto/Eunomia Migration
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move geometric multigrid level matrices to `leto::Array2`.
- Move GMG solve, FAS, transfer, and nonlinear-operator vectors to
  `leto::Array1`.
- Route GMG scalar constants and identities through Eunomia scalar traits.
- Replace nalgebra operator-overload math with explicit Leto-array helper
  kernels for matvec, residual, vector updates, and L2 norm.
- Keep public linear-solver nalgebra boundaries explicit as residual migration
  work.

### Verification
- `cargo fmt -p cfd-math --check`
- `cargo check -p cfd-math --no-default-features --lib`
- `cargo nextest run -p cfd-math --no-default-features gmg --status-level fail`
  (5/5 passed)
- `cargo clippy -p cfd-math --no-default-features --all-targets -- -D warnings`
- Clean GMG provider-residue scan
- `git diff --check` for the GMG files

---

## Sprint 1.96.117: cfd-math GMRES Leto/Eunomia Workspace
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move GMRES Krylov basis and Hessenberg storage to `leto::Array2`.
- Move Givens rotation coefficients, least-squares RHS, and triangular-solve
  output to `leto::Array1`/`Array2`.
- Route GMRES Givens scalar operations through Eunomia
  `RealField`/`NumericElement`.
- Propagate the Eunomia real-field contract through `LinearSolverChain`.
- Keep the remaining public nalgebra `DVector`/`RealField` linear-solver
  boundary explicit for the larger API migration.

### Verification
- `cargo fmt -p cfd-math --check`
- `cargo check -p cfd-math --no-default-features --lib`
- `cargo nextest run -p cfd-math --no-default-features gmres --status-level fail`
  (21/21 passed)
- `cargo clippy -p cfd-math --no-default-features --all-targets -- -D warnings`
- `git diff --check` for the GMRES and chain files

---

## Sprint 1.96.116: cfd-math AMG Restriction Leto Arrays
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Move standalone AMG restriction construction and validation to
  `leto::Array2`.
- Move restriction vector application to `leto::Array1`.
- Route Galerkin restriction projection through `leto_ops::MatrixProduct`.
- Strengthen restriction tests from finite-only checks to concrete `P^T v` and
  `P^T A P` value oracles.
- Keep broader public sparse/linear-solver nalgebra boundaries explicit as
  residual migration work.

### Verification
- `cargo fmt -p cfd-math --check`
- `cargo check -p cfd-math --no-default-features --lib`
- `cargo nextest run -p cfd-math --no-default-features restriction --status-level fail`
  (7/7 passed)
- `cargo clippy -p cfd-math --no-default-features --all-targets -- -D warnings`

---

## Sprint 1.96.115: cfd-math Linear Solver Leto Dense Bridge
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Add one `linear_solver::dense_bridge` path from legacy CSR/vector storage
  into `leto::Array2`/`Array1` and `leto_ops::solve`.
- Route `DirectSparseSolver` small-system dense LU retry through the bridge.
- Route multigrid cycle coarsest small-system dense solve through the bridge.
- Keep the current `nalgebra_sparse::CsrMatrix`/`DVector` public boundary
  explicit until the larger linear-solver API migration.
- Propagate the Leto real-scalar contract through `LinearSolverChain` where it
  calls the direct solver.
- Keep `rsparse` sparse LU recorded as residual primary-path provider work.

### Verification
- `cargo fmt -p cfd-math --check`
- `cargo check -p cfd-math --no-default-features --lib`
- `cargo nextest run -p cfd-math --no-default-features multigrid::cycles direct_solver --status-level fail`
  (9/9 passed)
- `cargo clippy -p cfd-math --no-default-features --all-targets -- -D warnings`

---

## Sprint 1.96.114: cfd-math Sparse Builder Leto Construction Bridge
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Centralize sparse Leto/nalgebra CSR conversion in `sparse::bridge`.
- Route `SparseMatrixBuilder` CSR construction through `leto_ops::CsrMatrix`.
- Route the empty `ParallelAssembly::block_diagonal` construction edge through
  `leto_ops::CsrMatrix`.
- Keep the residual `nalgebra_sparse::CsrMatrix`/`DVector` public storage
  boundary explicit for the next sparse/linear-solver migration slice.

### Verification
- `cargo check -p cfd-math --no-default-features --lib`
- `cargo nextest run -p cfd-math --no-default-features sparse --status-level fail`
  (18/18 passed)
- `cargo clippy -p cfd-math --no-default-features --all-targets -- -D warnings`

---

## Sprint 1.96.113: cfd-math Sparse Extension Leto Provider Consumption
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Route sparse diagonal extraction through `leto_ops::CsrMatrix::diagonal`.
- Route sparse scalar/value, row, and column scaling through Leto CSR provider
  methods.
- Route sparse Frobenius norm, diagonal dominance, and condition-estimate
  utilities through Leto CSR provider methods.
- Keep the residual `nalgebra_sparse::CsrMatrix`/`DVector` storage boundary
  explicit for the next sparse/linear-solver migration slice.

### Verification
- `cargo fmt -p cfd-math --check`
- `cargo check -p cfd-math --no-default-features --lib`
- `cargo nextest run -p cfd-math --no-default-features sparse --status-level fail`
  (18/18 passed)
- `cargo clippy -p cfd-math --no-default-features --all-targets -- -D warnings`

---

## Sprint 1.96.112: cfd-math AMG Leto SpGEMM Consumption
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Route cfd-math CSR×CSR products through `leto_ops::spgemm`.
- Route cfd-math SpMV products through `leto_ops::spmv_into`.
- Use the fallible Leto-backed path for AMG Galerkin operators.
- Route AMG restriction transpose through `leto_ops::CsrMatrix::transpose`.
- Keep the residual `nalgebra_sparse::CsrMatrix`/`DVector` storage boundary
  explicit for the next sparse/linear-solver migration slice.

### Verification
- `cargo fmt -p cfd-math --check`
- `cargo check -p cfd-math --no-default-features --lib`
- `cargo nextest run -p cfd-math --no-default-features sparse --status-level fail`
  (17/17 passed)
- `cargo nextest run -p cfd-math --no-default-features interpolation --status-level fail`
  (15/15 passed)
- `cargo nextest run -p cfd-math --no-default-features amg --status-level fail`
  (6/6 passed)
- `cargo clippy -p cfd-math --no-default-features --all-targets -- -D warnings`

---

## Sprint 1.96.111: Leto-ops CSR Product Provider Gap
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Add a Leto-owned CSR×CSR product for the CFDrs AMG Galerkin path.
- Keep sparse product behavior provider-owned rather than implementing a
  CFDrs-local `nalgebra_sparse` replacement.
- Record the remaining CFDrs consumer migration to `leto_ops::spgemm`,
  `leto_ops::CsrMatrix`, and `leto::Array1`.

### Verification
- `cargo fmt -p leto-ops --check`
- `cargo check -p leto-ops`
- `cargo nextest run -p leto-ops --test ops_tests sparse --status-level fail`
  (14/14 passed)
- `cargo clippy -p leto-ops --all-targets -- -D warnings`
- `cargo doc -p leto-ops --no-deps`

---

## Sprint 1.96.110: cfd-math SIMD Leto/Eunomia Providers
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Replace nalgebra vector ownership in the SIMD helper cone with Leto `Array1`.
- Route SIMD sparse matvec through `leto_ops::CsrMatrix`/`spmv`.
- Route generic SIMD scalar bounds through Eunomia `RealField`/`NumericElement`.
- Record sparse and linear-solver provider residue as separate follow-up work.

### Verification
- `cargo fmt -p cfd-math --check`
- `cargo check -p cfd-math --no-default-features --lib`
- `cargo nextest run -p cfd-math --no-default-features simd --status-level fail`
  (26/26 passed, 318 skipped)
- `cargo clippy -p cfd-math --no-default-features --all-targets -- -D warnings`
- Focused SIMD provider-residue scan.

---

## Sprint 1.96.109: cfd-math Nonlinear Solver Leto Vectors
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Replace nalgebra vector ownership in the nonlinear-solver cone with Leto
  `Array1`/`Array2`.
- Keep Anderson and JFNK vector math on one shared local Leto helper module.
- Route JFNK scalar math through Eunomia `RealField`/`FloatElement`.
- Record sparse and linear-solver provider residue as separate follow-up work.

### Verification
- `cargo fmt -p cfd-math --check`
- `cargo check -p cfd-math --no-default-features --lib`
- `cargo nextest run -p cfd-math --no-default-features nonlinear --status-level fail`
  (9/9 passed, 335 skipped)
- `cargo clippy -p cfd-math --no-default-features --all-targets -- -D warnings`
- Focused nonlinear-solver provider-residue scan.

---

## Sprint 1.96.108: cfd-math DG Leto Dense Arrays
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Replace nalgebra dense vector/matrix ownership in the high-order DG cone
  with Leto `Array1`/`Array2`.
- Route DG projection, derivative, RHS, and implicit Newton correction solves
  through Leto dense solve helpers with typed errors.
- Update DG docs and DG-related benchmarks to use Leto arrays.
- Record sparse and linear-solver provider residue as separate follow-up work.

### Verification
- `cargo fmt -p cfd-math --check`
- `cargo check -p cfd-math --no-default-features --lib`
- `cargo check -p cfd-math --no-default-features --bench dg_benchmarks`
- `cargo check -p cfd-math --no-default-features --bench flux_alloc_bench`
- `cargo nextest run -p cfd-math --no-default-features dg --status-level fail`
  (62/62 passed, 282 skipped)
- `cargo clippy -p cfd-math --no-default-features --all-targets -- -D warnings`
- `cargo test --doc -p cfd-math --no-default-features` (3 passed, 3 ignored)
- Focused high-order/DG-bench provider-residue scan.

---

## Sprint 1.96.107: cfd-math Spectral Leto Dense Arrays
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Replace nalgebra dense vector/matrix ownership in the high-order spectral
  cone with Leto `Array1`/`Array2`.
- Keep derivative, stiffness, dot-product, and matrix-vector operations shared
  through local Leto helpers.
- Surface spectral L2 projection mass-solve failures as typed solver errors
  instead of falling back silently.
- Migrate spectral assembly local dense matrices/RHS values and debug CSR
  materialization to Leto arrays.
- Record sparse and linear-solver provider residue as separate follow-up work.

### Verification
- `cargo fmt -p cfd-math --check`
- `cargo check -p cfd-math --no-default-features --lib`
- `cargo nextest run -p cfd-math --no-default-features spectral --status-level fail`
- `cargo clippy -p cfd-math --no-default-features --all-targets -- -D warnings`
- Focused spectral provider-residue scan.

---

## Sprint 1.96.106: cfd-math WENO Eunomia Scalars
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Remove the remaining nalgebra scalar-trait import from the cfd-math
  high-order WENO cone.
- Keep WENO5/WENO7 constants and nonlinear weights routed through Eunomia
  provider helpers without changing the reconstruction formulas.
- Record DG, sparse, and linear-solver provider residue as separate follow-up
  work.

### Verification
- `cargo fmt -p cfd-math --check`
- `cargo check -p cfd-math --no-default-features --lib`
- `cargo nextest run -p cfd-math --no-default-features weno --status-level fail`
  (6/6 passed, 338 skipped)
- `cargo clippy -p cfd-math --no-default-features --all-targets -- -D warnings`
- Focused WENO direct-provider residue scan.

### Residual Work
- Replace cfd-math sparse and linear-solver nalgebra/Leto storage and solver
  boundaries in separate slices.

---

## Sprint 1.96.105: cfd-1d Solver-Core Eunomia Scalars
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Remove direct `num_traits` bounds from the cfd-1d primary network
  solver-core scalar contract.
- Route Anderson acceleration, convergence checks, linear-system
  equilibration/Jacobi preconditioning, SPD detection, residual norms, and
  diagnostic scalar extraction through Eunomia/cfd-core provider APIs.
- Keep nalgebra/nalgebra-sparse storage replacement as a separate Leto-backed
  solver boundary increment.

### Verification
- `cargo fmt -p cfd-1d --check`
- `cargo check -p cfd-1d`
- `cargo nextest run -p cfd-1d` (725/725 passed, 3 skipped)
- Focused solver-core direct-provider residue scan.

### Residual Work
- Remove remaining cfd-1d direct provider residue in vascular
  Bessel/Womersley, resistance scalar traits, tests, and benches.
- Replace the active nalgebra/nalgebra-sparse linear-system storage boundary
  with Leto-backed dense/sparse provider APIs.
- Continue Hephaestus/Moirai migrations in GPU and execution-policy cones.

---

## Sprint 1.96.104: cfd-1d Network Wrapper Eunomia Scalars
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Remove direct num-traits scalar construction and nalgebra scalar bridge
  conversions from the cfd-1d network wrapper seam.
- Route network characteristic length, Picard resistance refresh,
  hematocrit propagation, Pries phase-separation bridge values, coefficient
  validation, and parallel edge conductance through Atlas scalar providers.
- Keep solver-core scalar compatibility as the next cleanup seam.

### Verification
- `cargo fmt -p cfd-1d --check`
- `cargo check -p cfd-1d`
- `cargo nextest run -p cfd-1d` (725/725 passed, 3 skipped)
- Focused wrapper/matrix/problem direct-provider residue scan.

### Residual Work
- Superseded by Sprint 1.96.105 for the `NetworkSolveScalar` compatibility
  bounds.
- Continue vascular Bessel/Womersley, resistance scalar traits, tests/benches,
  and nalgebra/nalgebra-sparse storage provider replacement.

---

## Sprint 1.96.103: cfd-1d Network Blueprint/Sink Eunomia Scalars
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Remove direct scalar construction from the canonical cfd-1d
  blueprint-to-network entry point.
- Route blueprint resistances, areas, physical constants, boundary values,
  serpentine segment counts, flow probes, and blood defaults through the Atlas
  scalar conversion provider.
- Keep the broader `domain/network/wrapper.rs` direct provider residue as the
  next cleanup seam.

### Verification
- `cargo fmt -p cfd-1d`
- `cargo check -p cfd-1d`
- `cargo nextest run -p cfd-1d` (725/725 passed, 3 skipped)

### Residual Work
- Migrate `domain/network/wrapper.rs` off direct `FromPrimitive`,
  `T::from_f64`, `T::from_usize`, and generic `.abs()` usage.
- Remove the inherited `FromPrimitive` requirement from `Network::new` after
  the wrapper impl is migrated.
- Continue solver-core, vascular Bessel/Womersley, tests/benches, and
  nalgebra/nalgebra-sparse storage provider replacement.

---

## Sprint 1.96.102: cfd-1d Domain Components Eunomia Scalars
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Remove direct `num_traits` bounds from the cfd-1d domain-components seam.
- Route component pressure-drop absolute values and constants/defaults through
  Eunomia/Atlas scalar provider APIs.
- Keep broader solver-core, domain-network, vascular, and nalgebra storage
  residues as separate follow-up work.

### Verification
- `cargo fmt -p cfd-1d --check`
- `cargo check -p cfd-1d`
- `cargo nextest run -p cfd-1d` (725/725 passed, 3 skipped)
- Focused domain-components provider-residue scan.
- `cargo clippy -p cfd-1d --all-targets -- -D warnings` rerun confirms the
  touched `channels.rs` lint is gone; unrelated all-target lint debt remains.

### Residual Work
- Remove remaining direct `num-traits` usage from cfd-1d solver-core,
  domain-network, vascular Bessel/Womersley, tests, and benches.
- Replace remaining cfd-1d nalgebra/nalgebra-sparse storage boundaries with
  Leto where provider functionality is available.
- Clear unrelated all-target clippy debt before using the broad cfd-1d clippy
  gate as closure evidence.

---

## Sprint 1.96.101: cfd-1d Channel/Branching/Analysis Eunomia Scalars
**Status**: Completed
**Start Date**: July 4, 2026

### Sprint Objectives
- Remove direct `num_traits` bounds from the next coherent cfd-1d
  channel/branching/analysis provider seam.
- Route channel solver constants/math and network-analysis scalar conversion
  through `SafeFromF64`, `FloatElement`, and `NumericElement`.
- Keep broader solver-core, network/component, and vascular provider residues
  as separate follow-up work.

### Verification
- `cargo fmt -p cfd-1d`
- `cargo check -p cfd-1d`
- `cargo nextest run -p cfd-1d` (725/725 passed, 3 skipped)
- Focused channel/branching/analyzer and solver-analysis residue scans.

### Residual Work
- Remove remaining direct `num-traits` usage from cfd-1d solver-core,
  domain-network/components, vascular Bessel/Womersley, and package tests.
- Replace remaining cfd-1d nalgebra/nalgebra-sparse storage boundaries with
  Leto where provider functionality is available.
- Clear unrelated all-target clippy debt before using the broad cfd-1d clippy
  gate as closure evidence.

---

## Sprint 1.96.100: cfd-core Fluid Dynamics Operations Eunomia Scalars
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct `num_traits::FromPrimitive` from
  `physics::fluid_dynamics::operations`.
- Route vorticity, divergence, kinetic-energy, and enstrophy scalar constants
  through Eunomia `NumericElement`.
- Preserve current Moirai parallel iteration and record nalgebra vector/storage
  migration separately.

### Verification
- `cargo check -p cfd-core`
- `cargo nextest run -p cfd-core fluid_dynamics::operations` (3/3 passed)
- Touched-file rustfmt and focused residue scans.

### Residual Work
- Migrate `VelocityField<T>` and vector operations away from nalgebra storage.
- Migrate RANS/turbulence traits and Rhie-Chow away from nalgebra scalar/vector
  contracts.

---

## Sprint 1.96.99: cfd-core Fluid Dynamics Service Eunomia Scalars
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct `num_traits::{Float, FromPrimitive}` from
  `physics::fluid_dynamics::service`.
- Route pipe-flow constants, powers, square roots, logarithms, and absolute
  convergence checks through Eunomia `FloatElement`/`NumericElement`.
- Preserve the existing `ConstantPropertyFluid<T>` storage contract and record
  vector/storage migration separately.

### Verification
- `cargo check -p cfd-core`
- `cargo nextest run -p cfd-core fluid_dynamics::service` (2/2 passed)
- Touched-file rustfmt and focused residue scans.

### Residual Work
- Migrate `fluid_dynamics::operations` away from direct `num_traits`.
- Migrate `VelocityField<T>` and vector operations away from nalgebra storage.

---

## Sprint 1.96.98: cfd-core Flow Regime Eunomia Scalars
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct `nalgebra::RealField` and `num_traits::ToPrimitive` from core
  flow-regime classification.
- Route Reynolds and Mach conversion through Eunomia `NumericElement::to_f64`.
- Replace the service wrapper's old direct `num_traits::ToPrimitive` bound
  with the migrated Eunomia real-scalar contract.

### Verification
- `cargo check -p cfd-core`
- `cargo nextest run -p cfd-core flow_regime` (3/3 passed)
- Touched-file rustfmt and focused residue scans.

### Residual Work
- Migrate `FluidDynamicsService` pipe-flow friction-factor formulas away from
  direct `num_traits::{Float, FromPrimitive}`.
- Migrate `fluid_dynamics::operations` and `VelocityField<T>` away from
  nalgebra vector/storage boundaries.

---

## Sprint 1.96.97: cfd-3d Spectral Diagnostics Eunomia Conversion
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct `num_traits::{FromPrimitive, ToPrimitive}` scalar conversion
  from the Apollo/Leto spectral diagnostics helper.
- Route velocity-component staging into Leto arrays through Eunomia
  `NumericElement::to_f64`.
- Record the inherited `cfd-core::VelocityField<T>` nalgebra storage boundary
  as a separate residual.

### Verification
- `cargo check -p cfd-3d`
- `cargo nextest run -p cfd-3d diagnostics` (5/5 passed)
- Touched-file rustfmt and focused residue scans.

### Residual Work
- Replace the upstream `VelocityField<T>` nalgebra `Vector3`/`RealField`
  storage boundary.
- Continue replacing larger `cfd-3d` spectral dense-linalg nalgebra surfaces
  with Leto.

---

## Sprint 1.96.96: cfd-2d MRT/Carreau-Yasuda Eunomia Scalars
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct `nalgebra::RealField` and `num_traits` scalar dispatch from
  the remaining LBM MRT and Carreau-Yasuda rheology/collision paths.
- Route MRT moment constants, relaxation rates, equilibrium moments,
  Carreau-Yasuda powers, square roots, max clamps, and D2Q9 velocity
  conversions through Eunomia/provider helpers.

### Verification
- `cargo check -p cfd-2d`
- `cargo nextest run -p cfd-2d lbm` (31/31 passed)
- `cargo nextest run -p cfd-2d carreau_yasuda` (5/5 passed)
- Touched-file rustfmt and focused residue scans.

### Residual Work
- Continue broader CFDrs replacement of direct `nalgebra`, `num_traits`,
  `ndarray`, raw `wgpu`, `tokio`, and `rayon` surfaces with Atlas providers.

---

## Sprint 1.96.95: cfd-2d LBM Eunomia Scalar Seam
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct `nalgebra::RealField` and `num_traits::FromPrimitive` from the
  authoritative `cfd-2d` LBM macroscopic, D2Q9 equilibrium, collision trait,
  and BGK collision seam.
- Route scalar constants, zero/one identities, weights, pressure/stress
  helpers, and BGK viscosity/omega construction through Eunomia-backed local
  scalar helpers.
- Keep MRT and Carreau-Yasuda broader inherent scalar-provider cleanup as a
  named residual rather than silently broadening this slice. Sprint 1.96.96
  closes that residual.

### Verification
- `cargo check -p cfd-2d`
- `cargo nextest run -p cfd-2d lbm` (31/31 passed)
- Touched-file rustfmt, touched-file `git diff --check`, and focused residue
  scans.

### Residual Work
- MRT and Carreau-Yasuda inherent scalar-provider migration is closed by
  Sprint 1.96.96.
- Continue broader CFDrs replacement of direct `nalgebra`, `num_traits`,
  `ndarray`, raw `wgpu`, `tokio`, and `rayon` surfaces with Atlas providers.

---

## Sprint 1.96.94: cfd-3d Wall Functions Eunomia Scalars
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct `nalgebra::RealField` and `num_traits` scalar math from the
  3D Spalding wall-law helper.
- Route constants, transcendental functions, convergence checks, and clamps
  through Eunomia provider traits.
- Record why direct `GpuContext` Hephaestus replacement needs WGPU
  API-version alignment first.

### Verification
- `cargo check -p cfd-3d`
- `cargo nextest run -p cfd-3d wall_functions` (3/3 passed)
- Touched-file rustfmt, touched-file `git diff --check`, and focused residue
  scans.

### Residual Work
- Continue broader turbulence scalar-provider migration.
- Align CFDrs raw WGPU usage with Hephaestus' WGPU version before replacing
  direct raw-device context fields.

---

## Sprint 1.96.93: cfd-3d Apollo Fourier Eunomia Bounds
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct `num_traits` scalar conversion from the Apollo-backed
  `cfd-3d` Fourier wrapper.
- Propagate migrated `FloatElement` blood/cavitation bounds through the
  downstream 2D/3D/validation consumers needed by the focused 3D spectral gate.

### Verification
- `cargo check -p cfd-3d`
- `cargo check -p cfd-validation`
- `cargo nextest run -p cfd-3d --test fourier_validation` (12/12 passed)
- Touched-file rustfmt and focused residue scans.

### Residual Work
- Continue broader CFDrs migration across remaining direct `num_traits`,
  `nalgebra`, direct `wgpu`, and execution-provider surfaces.

---

## Sprint 1.96.92: cfd-core Fåhræus-Lindqvist Eunomia Scalars
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct `num_traits::FromPrimitive` construction from the
  Fåhræus-Lindqvist blood model.
- Route Pries/Secomb microvascular viscosity math through Eunomia.
- Close the local blood-model scalar `num_traits` construction residue.

### Dependency Audit Findings
- Fåhræus-Lindqvist needs scalar constants, zero/one identities, real powers,
  exponential evaluation, absolute value, and max clamping.
- Cross, Casson, Carreau-Yasuda, and `BloodModel` had already moved to Eunomia;
  this was the last local blood-model scalar holdout.

### Sprint Backlog Items
- [x] Replace Fåhræus-Lindqvist constants with Eunomia scalar construction.
- [x] Replace Pries/Secomb `powf`/`exp`/`abs`/max generic math with Eunomia
  APIs.
- [x] Verify with focused residue scan, compile check, and nextest.

### Residual Work
- Continue migrating the broader fluid trait `RealField` boundary and
  non-blood provider holdouts.
- cfd-core clippy remains blocked by unrelated existing boundary applicator and
  Rhie-Chow lints.

---

## Sprint 1.96.91: cfd-core Casson/Carreau Blood Eunomia Scalars
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct `num_traits::FromPrimitive` construction from Casson blood.
- Route Casson apparent-viscosity and temperature-correction scalar math through
  Eunomia.
- Tighten `BloodModel` dispatch to Eunomia after Casson and Carreau-Yasuda share
  the provider scalar surface.

### Dependency Audit Findings
- Casson needs scalar constants, zero/one identities, square roots, and
  exponential evaluation.
- Carreau-Yasuda had already moved to Eunomia but could not compile through the
  shared dispatch until `BloodModel` stopped requiring the stale bound.

### Sprint Backlog Items
- [x] Replace Casson constructors and hematocrit scaling with Eunomia scalar
  construction.
- [x] Replace Casson `sqrt`/`exp`/zero/one generic math with Eunomia APIs.
- [x] Replace `BloodModel` direct `FromPrimitive` bound with Eunomia
  `FloatElement`.
- [x] Verify with focused residue scan, compile check, and nextest.

### Residual Work
- Continue migrating Fåhræus-Lindqvist and the broader fluid trait `RealField`
  boundary.
- cfd-core clippy remains blocked by unrelated existing boundary applicator and
  Rhie-Chow lints.

---

## Sprint 1.96.90: cfd-core Cross Blood Eunomia Scalars
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct `num_traits::FromPrimitive` construction from the Cross blood
  model.
- Route Cross apparent-viscosity scalar math through Eunomia.

### Dependency Audit Findings
- Cross blood needs scalar constants, comparison against zero, one identity,
  and real-power evaluation.
- The surrounding `Fluid<T>` trait still owns the inherited
  `nalgebra::RealField` boundary.

### Sprint Backlog Items
- [x] Replace Cross normal-blood constants with Eunomia construction.
- [x] Replace Cross zero/one and `powf` calls with Eunomia APIs.
- [x] Verify with focused residue scan, compile check, and nextest.

### Residual Work
- Continue migrating blood-fluid scalar holdouts in Casson, Carreau-Yasuda,
  Fåhræus-Lindqvist, `BloodModel`, and the fluid trait.

---

## Sprint 1.96.89: cfd-core Cavitation Scalar Closeout
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove `nalgebra::RealField` from the remaining cavitation mass-transfer
  models.
- Remove the fake-generic heterogeneous legacy adapter.
- Route mass-transfer constants and math through Eunomia.

### Dependency Audit Findings
- Cavitation mass-transfer models need scalar arithmetic, comparisons,
  constants, powers, roots, absolute value, and min/max.
- The heterogeneous legacy adapter already delegates into a concrete `f64`
  selective-cavitation model, so a generic `T` API was misleading.

### Sprint Backlog Items
- [x] Replace `CavitationModel<T>` and `ZgbParams<T>` bounds with Eunomia.
- [x] Replace direct scalar construction and math with Eunomia APIs.
- [x] Replace `heterogeneous_inception_threshold_pa<T>` with an honest `f64`
  contract.
- [x] Verify with focused residue scan, compile check, and nextest.

### Residual Work
- Continue the broader CFDrs Atlas migration outside cavitation:
  nalgebra/ndarray, raw WGPU, RustFFT, and direct scheduler/memory-provider
  holdouts remain in other modules and crates.

---

## Sprint 1.96.88: cfd-core Nuclei Transport Eunomia Scalars
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove `nalgebra::RealField` from nuclei transport.
- Route nuclei transport constants, zero/one identities, and exponential decay
  through Eunomia.

### Dependency Audit Findings
- Nuclei transport needs scalar arithmetic, comparisons, constants, and
  exponential evaluation.
- Existing tests cover the migrated affine pressure coupling, dissolution
  rate, exponential transit decay, and diffusion accessor.

### Sprint Backlog Items
- [x] Replace `nuclei_adjusted_vapor_pressure<T>` bounds with Eunomia.
- [x] Replace `NucleiTransportConfig<T>` bounds/defaults with Eunomia.
- [x] Replace `NucleiTransport<T>` bounds and source/sink math with Eunomia.
- [x] Verify with focused residue scan, compile check, and nextest.

### Residual Work
- `models.rs` and `heterogeneous_nucleation.rs` were subsequently closed by
  Sprint 1.96.89.

---

## Sprint 1.96.87: cfd-core Venturi Cavitation Eunomia Scalars
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove `nalgebra::RealField` and direct `num_traits::FromPrimitive` from
  `VenturiCavitation<T>`.
- Route Venturi cavitation constants and math through Eunomia.

### Dependency Audit Findings
- Venturi cavitation needs only scalar arithmetic, comparisons, constants,
  integer powers, tangent, and absolute value.
- Existing tests cover the continuity, Bernoulli, cavitation-number, cavity
  length, and conical-volume formulas.

### Sprint Backlog Items
- [x] Replace `VenturiCavitation<T>` bounds with Eunomia.
- [x] Remove direct `num_traits::FromPrimitive` and `T::zero`/`T::one` usage
  from Venturi.
- [x] Verify with focused residue scan, compile check, and nextest.

### Residual Work
- At this point, continue replacing cavitation scalar holdouts in `models.rs`,
  `nuclei_transport.rs`, and `heterogeneous_nucleation.rs`. Sprint 1.96.88
  subsequently removed the `nuclei_transport.rs` holdout.

---

## Sprint 1.96.86: cfd-core Cavitation Eunomia Scalar Cone
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove `nalgebra::RealField` and direct `num_traits::FromPrimitive` from the
  touched cavitation scalar cone.
- Route cavitation scalar constants and math through Eunomia.
- Add value-semantic tests for newly migrated cavitation number and material
  damage APIs.

### Dependency Audit Findings
- Rayleigh-Plesset, biological damage, regime analysis, cavitation number, and
  material damage need only scalar arithmetic, comparisons, constants, powers,
  roots, exponentials, and finite checks.
- `CavitationDamage::incubation_period` keeps the only `f64` conversion at the
  explicit `u64` cycle-count API boundary.

### Sprint Backlog Items
- [x] Replace Rayleigh-Plesset and sonoluminescence scalar bounds with Eunomia.
- [x] Replace biological-damage scalar bounds with Eunomia.
- [x] Replace cavitation-regime scalar bounds with Eunomia.
- [x] Replace cavitation-number scalar bounds with Eunomia.
- [x] Replace material-damage scalar bounds with Eunomia.
- [x] Add closed-form value tests for the newly migrated APIs.
- [x] Verify with focused residue scan, compile check, and nextest.

### Residual Work
- At this point, continue replacing cavitation scalar holdouts in `models.rs`,
  `venturi.rs`, `nuclei_transport.rs`, and `heterogeneous_nucleation.rs`.
  Sprint 1.96.87 subsequently removed the `venturi.rs` holdout.
- Full `cfd-core` clippy remains blocked by unrelated existing lints in
  boundary applicator and Rhie-Chow tests.

---

## Sprint 1.96.85: cfd-core Hemolysis Eunomia Scalars
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove `nalgebra::RealField` from the hemolysis calculator and platelet
  activation scalar boundary.
- Route hemolysis constants, zero/one identities, and exponential evaluation
  through Eunomia.

### Dependency Audit Findings
- Hemolysis calculator logic needs scalar comparisons, arithmetic, constants,
  and geometry-derived area calculations.
- Platelet activation needs scalar comparisons, arithmetic, and `exp`.
- No direct `rayon`/`tokio` source dependency remained in the exact active
  scan; only comments mention those names.

### Sprint Backlog Items
- [x] Replace `HemolysisCalculator<T>` bounds with Eunomia.
- [x] Replace `PlateletActivation<T>` bounds with Eunomia.
- [x] Remove direct `num_traits::FromPrimitive` and `T::zero`/`T::one` usage
  from hemolysis.
- [x] Verify with focused residue scan, compile check, and nextest.

### Residual Work
- Continue replacing `cfd-core` nalgebra scalar/vector contracts in compute,
  boundary, cavitation, fluid, mesh, and fluid-dynamics modules.

---

## Sprint 1.96.84: cfd-2d Grid/FVM Leto-Eunomia Boundary
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove the `StructuredGrid2D<T>` nalgebra `RealField` contract inherited by
  FVM.
- Route grid-center vectors through Leto.
- Route grid scalar construction through Eunomia.

### Dependency Audit Findings
- FVM only needs grid dimensions, spacing, and cell-center access from the
  structured grid.
- `Grid2D::cell_center` was the grid trait vector boundary that still exposed
  nalgebra `Vector2`.
- Several validation/FDM/LBM consumers accessed cell centers through nalgebra
  `.x`/`.y` fields and needed Leto indexing.

### Sprint Backlog Items
- [x] Replace `Grid2D::cell_center` return type with Leto `Vector2`.
- [x] Remove `RealField`/`FromPrimitive` from `StructuredGrid2D`.
- [x] Remove `RealField` from `UnstructuredGrid2D` storage.
- [x] Replace adaptive-grid resolution factors with Eunomia conversion.
- [x] Remove inherited `RealField` from `FvmSolver`.
- [x] Update grid-center consumers.
- [x] Verify with focused residue scans, compile checks, and nextest.

### Residual Work
- Continue replacing remaining `cfd-2d` nalgebra/num-traits surfaces outside
  the grid/FVM cone, especially momentum, pressure-velocity, SIMPLE/PIMPLE, and
  validation-oriented solver paths.
- Wire or remove the currently unlisted `tests_poisson_mms.rs` file in a
  separate FDM test-discovery cleanup slice.

---

## Sprint 1.96.83: cfd-2d FVM Leto Vector Boundary
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Replace FVM-local `Vector2` ownership with Leto geometry.
- Extend Leto only with fixed-vector operations required by the consumer.
- Preserve the remaining grid-level `RealField` dependency as explicit
  residual work.

### Dependency Audit Findings
- FVM face geometry needs 2-D fixed vectors, dot products, norms,
  normalization, zeros, indexing, and construction.
- Leto already owns fixed-vector storage and dot products; it needed a 2-D
  alias and norm/normalization methods for this consumer.
- `FvmSolver` still inherits `RealField` through `StructuredGrid2D<T>`, not
  through FVM vector storage.

### Sprint Backlog Items
- [x] Add Leto `Vector2<T>` alias.
- [x] Add generic Leto fixed-vector norm/normalization operations.
- [x] Replace FVM `Face<T>` vector storage with Leto.
- [x] Replace FVM solver velocity-field vectors with Leto.
- [x] Verify with provider and consumer compile checks.
- [x] Verify with focused nextest coverage.

### Residual Work
- Historical residual closed by Sprint 1.96.84: `StructuredGrid2D<T>` no
  longer requires nalgebra `RealField` for the grid/FVM boundary.
- Continue replacing other `cfd-2d` nalgebra vector/matrix surfaces with Leto
  or Gaia-owned geometry as appropriate.

---

## Sprint 1.96.82: cfd-2d FVM Eunomia Scalar Constants
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Replace FVM configuration, solver face-center, and flux scalar construction
  with Eunomia.
- Propagate the config provider bound through its direct solver owner.
- Preserve existing solver math and flux residuals as explicit follow-up work.
- Add direct value-semantic coverage for configuration defaults.

### Dependency Audit Findings
- `FvmConfig` only needs scalar storage and default literal construction.
- `FvmSolver` face-center construction only needs exact integer-index
  conversion through Eunomia.
- FVM flux schemes only need scalar constants and finite diffusion validation
  beyond their existing `RealField` math contract.

### Sprint Backlog Items
- [x] Replace config `RealField`/`FromPrimitive` imports and bounds.
- [x] Replace default scalar constants with Eunomia construction.
- [x] Add value-semantic default-config coverage.
- [x] Add `FloatElement` to the `FvmSolver` owner bound.
- [x] Replace solver face-center conversions with Eunomia helpers.
- [x] Replace flux scalar constants and diffusion validation with Eunomia.
- [x] Add value-semantic flux invalid-diffusion tests.
- [x] Verify with `cargo check -p cfd-2d`.
- [x] Verify with `cargo check -p cfd-3d`.
- [x] Verify with focused `cargo nextest`.

### Residual Work
- Migrate FVM geometry, solver, and flux contracts away from nalgebra
  `RealField` through the structured-grid contract.

---

## Sprint 1.96.81: cfd-2d CFL Eunomia Scalars
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Replace standalone CFL scalar contracts with Eunomia.
- Preserve the documented advection, diffusion, QUICK, and max-dt stability
  formulas.
- Keep verification focused on the CFL library tests.

### Dependency Audit Findings
- `CFLCalculator` needed scalar constants, absolute value, ordering, and
  literal construction only.
- `NumericElement` owns `ZERO`, `ONE`, and `abs`; `FloatElement` owns
  `from_f64`.

### Sprint Backlog Items
- [x] Replace `RealField`/`FromPrimitive` imports and bounds.
- [x] Replace direct nalgebra/num-traits scalar operations with Eunomia calls.
- [x] Verify with cfd-2d compile and library-filtered CFL nextest.

### Residual Work
- Broader cfd-2d stability/time-scheme, turbulence, validation, and solver
  modules still contain nalgebra/num-traits surfaces.

---

## Sprint 1.96.80: cfd-core Material Eunomia Traits
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove nalgebra scalar trait requirements from material solid/interface
  contracts.
- Keep the fluid database boundary explicit until `Fluid<T>` migrates.
- Add value-semantic tests for migrated formulas and constants.

### Dependency Audit Findings
- Solid/interface material contracts only need scalar constants, arithmetic,
  and trigonometry, which Eunomia owns.
- `MaterialDatabase` still stores `Box<dyn Fluid<T>>`, so its `RealField`
  bound is retained as a fluid-boundary residual instead of silently hidden.

### Sprint Backlog Items
- [x] Replace `SolidProperties` and `InterfaceProperties` `RealField` bounds
  with Eunomia bounds.
- [x] Replace `ElasticSolid`, `WettingProperties`, and
  `FluidSolidInterface` `RealField` bounds with Eunomia bounds.
- [x] Add value-semantic material tests and focused verification.

### Residual Work
- Migrate `Fluid<T>` and `MaterialDatabase` fluid storage away from
  `RealField`.

---

## Sprint 1.96.79: cfd-core Velocity Leto Vector
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Replace the cfd-core velocity value object's vector storage with Leto.
- Move the immediate physical-parameters gravity vector to the same provider.
- Extend Leto geometry serialization instead of adding a CFDrs wrapper.

### Dependency Audit Findings
- Leto already exposes `geometry::Vector3<T>` with `new`, `zeros`, `norm`,
  division, and component access.
- CFDrs `Velocity` derives Serde, so Leto geometry needed provider-owned Serde
  support before it could replace nalgebra at that boundary.

### Sprint Backlog Items
- [x] Add Serde derives to Leto fixed 3D geometry value types.
- [x] Add direct `leto` dependency to `cfd-core`.
- [x] Replace `Velocity` and `PhysicalParameters::gravity` vector storage with
  `leto::geometry::Vector3`.
- [x] Verify with focused source scan, touched-file rustfmt, and
  `cargo check -p cfd-core`, `cargo nextest run -p cfd-core --lib`, and
  downstream `cfd-2d`/`cfd-3d`/`cfd-validation` compile checks.

### Residual Work
- Replace `ProblemAggregate`/`SimulationAggregate` `Domain<T>` and fluid
  `RealField` contracts plus material/hemolysis/fluid-dynamics nalgebra
  surfaces.

---

## Sprint 1.96.78: cfd-core Physics Value Eunomia Scalars
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Replace scalar-only physics value-wrapper `RealField` contracts with Eunomia.
- Keep vector-backed `Velocity` unchanged for the Leto/Gaia vector slice.
- Record remaining nalgebra scalar/vector surfaces explicitly.

### Dependency Audit Findings
- `Temperature`, `Pressure`, `ReynoldsNumber`, and `DimensionlessNumber` only
  needed scalar zero, absolute value, square root, and literal construction.
- `Velocity` still owns `nalgebra::Vector3`, so it is not scalar-only.

### Sprint Backlog Items
- [x] Replace scalar-only value wrapper bounds with Eunomia
  `FloatElement`/`NumericElement`.
- [x] Propagate the resulting bounds through immediate aggregate owners.
- [x] Verify with focused source scan, touched-file rustfmt, and
  `cargo check -p cfd-core` plus `cargo nextest run -p cfd-core --lib`.

### Residual Work
- Replace `Velocity`/vector contracts and material/hemolysis `RealField`
  surfaces with Atlas-owned providers.

---

## Sprint 1.96.77: cfd-math WGPU Boundary Reduction
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove `cfd-math`'s direct optional WGPU dependency.
- Route GPU metric synchronization/capability checks through `cfd-core`.
- Keep existing raw WGPU kernels unchanged in this bounded slice.

### Dependency Audit Findings
- `cfd-math` used WGPU only for `Device::poll` and
  `Features::TIMESTAMP_QUERY` in GPU dispatch metrics.
- Root package `gpu` feature still activated `dep:wgpu` directly.
- `cfd-core::compute::gpu` remains the raw WGPU kernel boundary.

### Sprint Backlog Items
- [x] Add `GpuContext::synchronize`.
- [x] Add `GpuContext::supports_timestamp_queries`.
- [x] Remove `dep:wgpu` from `cfd-math/gpu` and the root package `gpu`
  feature.
- [x] Verify with focused source scans, cfd-math GPU check/nextest, root GPU
  feature check, rustfmt, and diff hygiene.

### Residual Work
- Replace `cfd-core::compute::gpu` raw WGPU buffers, command encoders,
  pipelines, and WGSL kernels with Hephaestus abstractions.

---

## Sprint 1.96.76: cfd-core Hephaestus GPU Probe
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Make Hephaestus the active GPU device-acquisition provider for `cfd-core`.
- Keep existing raw WGPU kernels unchanged in this bounded slice.
- Record the remaining raw buffer/pipeline/shader ownership migration.

### Dependency Audit Findings
- `hephaestus-wgpu`: active through `cfd-core` after this slice.
- `rustfft`: no package exists in the workspace graph.
- `ndarray`: still resolves only through `numpy -> cfd-python`.
- `rayon`: still resolves only through Criterion dev-dependencies.

### Sprint Backlog Items
- [x] Add optional `hephaestus-wgpu` to `cfd-core` under the `gpu` feature.
- [x] Replace direct raw WGPU adapter probing with
  `hephaestus_wgpu::WgpuDevice::try_default`.
- [x] Verify with `cargo check -p cfd-core --features gpu`,
  `cargo nextest run -p cfd-core --features gpu --lib` (183/183 passed), and
  `cargo tree --workspace -i hephaestus-wgpu`.

### Residual Work
- Replace `cfd-core::compute::gpu` raw WGPU buffers, command encoders,
  pipelines, and WGSL kernels with Hephaestus abstractions.
- Replace `cfd-math::linear_solver::operators::gpu` raw WGPU feature checks
  and adapter coupling after the cfd-core GPU surface moves.

---

## Sprint 1.96.75: cfd-1d Non-Python ndarray Path Removal
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove the stale `sprs -> ndarray` path from cfd-1d.
- Remove the unused root workspace `ndarray` declaration.
- Keep cfd-1d test generics aligned with Eunomia scalar provider bounds.

### Dependency Audit Findings
- `sprs`: no source use in cfd-1d; manifest-only residue removed.
- `ndarray`: no active cfd-1d or cfd-3d inverse-tree package remains.
- `numpy`: still pulls workspace `ndarray` for cfd-python.
- `nalgebra-sparse`: remains actively used by cfd-1d and is reserved for a
  later Leto sparse migration slice.

### Sprint Backlog Items

#### cfd-1d ndarray Graph Cleanup
- [x] **CFD1D-DEPS-176 [patch]**: Remove unused `sprs` from cfd-1d.
- [x] **WORKSPACE-DEPS-177 [patch]**: Remove unused root workspace `ndarray`.
- [x] **CFD1D-EUNOMIA-178 [patch]**: Propagate `FloatElement` into affected test helpers.
- [x] **CFD1D-GATE-179 [patch]**: Verify with cfd-1d check, nextest, and ndarray inverse-tree audits.

## Sprint 1.96.74: cfd-3d Apollo ndarray Removal
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Resolve CFDrs' active Apollo packages through the side-by-side Atlas Apollo
  checkout rather than the older ndarray-backed Git revision.
- Keep cfd-3d and cfd-validation compiling against Eunomia scalar/complex
  contracts required by the local provider stack.
- Record the remaining non-Apollo `ndarray` inverse dependency path.

### Dependency Audit Findings
- `apollo-fft` and `apollo-nufft`: active package trees have no `ndarray`
  matches and resolve from `D:\atlas\repos\apollo`.
- `eunomia`: now owns the touched scalar/complex contracts in cfd-3d and
  cfd-validation.
- `ndarray`: remaining active cfd-3d graph path is `sprs -> cfd-1d`, not
  Apollo.
- `leto-ops`: sparse/operator replacement remains for later focused slices.

### Sprint Backlog Items

#### Apollo Provider Boundary
- [x] **CFD3D-APOLLO-172 [patch]**: Patch CFDrs to side-by-side Atlas Apollo/provider checkouts.
- [x] **CFD3D-EUNOMIA-173 [patch]**: Update cfd-3d FEM and spectral call sites for Eunomia provider contracts.
- [x] **VALIDATION-EUNOMIA-174 [patch]**: Propagate `FloatElement` through validation dev-dependency code.
- [x] **CFD3D-GATE-175 [patch]**: Verify with cfd-3d check, cfd-3d nextest, and Apollo ndarray dependency/source scans.

## Sprint 1.96.73: cfd-math Eunomia Geometric Multigrid
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct scalar conversion provider residue from geometric multigrid
  hierarchy construction and transfer weights.
- Remove stale AMG `FromPrimitive` bounds after the multigrid provider cleanup.
- Preserve current nalgebra dense/sparse/vector surfaces for later Leto
  migration.

### Dependency Audit Findings
- `eunomia::FloatElement`: now owns GMG hierarchy scalar construction, Poisson
  matrix constants, and transfer weights.
- `eunomia::NumericElement`: now owns grid-dimension-to-scalar conversion.
- `num-traits`: no direct provider residue appears in
  `linear_solver/preconditioners/multigrid`.
- `nalgebra`: remains for current multigrid dense/sparse/vector surfaces and
  `RealField` bounds pending Leto migration.

### Sprint Backlog Items

#### Geometric Multigrid Provider Boundary
- [x] **MATH-EUNOMIA-167 [patch]**: Replace GMG hierarchy and Poisson scalar constants with Eunomia.
- [x] **MATH-EUNOMIA-168 [patch]**: Replace GMG transfer weights with Eunomia.
- [x] **MATH-EUNOMIA-169 [patch]**: Remove stale AMG `FromPrimitive` import and bounds.
- [x] **MATH-TEST-170 [patch]**: Add value-semantic Poisson stencil and restriction-weight tests.
- [x] **MATH-GATE-171 [patch]**: Verify with rustfmt, cfd-math check, focused GMG/AMG nextest, and multigrid-wide residue scan.

## Sprint 1.96.72: cfd-math Eunomia Multigrid Interpolation
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct scalar conversion provider residue from multigrid interpolation
  operators and quality metrics.
- Preserve current nalgebra sparse/vector surfaces for later Leto migration.
- Verify with focused cfd-math check and interpolation nextest.

### Dependency Audit Findings
- `eunomia::FloatElement`: now owns interpolation scalar constants and
  index-distance conversion.
- `eunomia::NumericElement`: now owns interpolation quality row-sum extraction,
  absolute-value dispatch, and sparsity/average metric conversion.
- `num-traits`: no direct provider residue appears in
  `linear_solver/preconditioners/multigrid/interpolation.rs` after the static
  scan.
- `nalgebra`: remains for current sparse/vector surfaces and `RealField`
  bounds pending Leto migration.

### Sprint Backlog Items

#### Multigrid Interpolation Provider Boundary
- [x] **MATH-EUNOMIA-163 [patch]**: Replace interpolation scalar constants and index-distance conversion with Eunomia.
- [x] **MATH-EUNOMIA-164 [patch]**: Replace interpolation quality metric extraction and absolute-value paths with Eunomia.
- [x] **MATH-TEST-165 [patch]**: Remove touched test debug output and direct index-distance casts.
- [x] **MATH-GATE-166 [patch]**: Verify with `cargo check -p cfd-math` and focused interpolation nextest.

## Sprint 1.96.71: cfd-math Eunomia Multigrid Coarsening
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct scalar conversion provider residue from multigrid coarsening
  algorithms and quality analysis.
- Preserve current nalgebra sparse/vector surfaces for later Leto migration.
- Add value-semantic strength-matrix coverage.

### Dependency Audit Findings
- `eunomia::FloatElement`: now owns coarsening scalar constants, ratio/count
  conversion, distance thresholds, and fallback distance construction.
- `eunomia::NumericElement`: now owns strength-matrix absolute values,
  quality-analysis absolute values, f64 quality-metric extraction, and
  maximum-distance sentinel values.
- `num-traits`: no direct provider residue appears in
  `linear_solver/preconditioners/multigrid/coarsening`.
- `nalgebra`: remains for current sparse/vector surfaces and `RealField`
  bounds pending Leto migration.

### Sprint Backlog Items

#### Multigrid Coarsening Provider Boundary
- [x] **MATH-EUNOMIA-159 [patch]**: Replace coarsening algorithm scalar constants and count conversions with Eunomia.
- [x] **MATH-EUNOMIA-160 [patch]**: Replace coarsening quality sentinels, thresholds, absolute values, and f64 extraction with Eunomia.
- [x] **MATH-TEST-161 [patch]**: Add value-semantic strength-matrix grid-connectivity coverage.
- [x] **MATH-GATE-162 [patch]**: Verify with rustfmt, cfd-math check, focused nextest, and residue scan.

## Sprint 1.96.70: cfd-math Eunomia Multigrid Smoothers
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct scalar conversion fallback usage from multigrid smoother owner
  paths.
- Preserve current nalgebra sparse/vector surfaces for later Leto migration.
- Add value-semantic smoother update and eigenvalue coverage.

### Dependency Audit Findings
- `eunomia::FloatElement`: now owns smoother scalar constants, Chebyshev
  eigenvalue/default constants, and immediate AMG smoother-owner threshold and
  relaxation construction.
- `eunomia::NumericElement`: now owns smoother diagonal and AMG complexity
  absolute-value dispatch in the touched paths.
- `num-traits`: no direct provider residue appears in `smoothers.rs`; `amg.rs`
  still retains `FromPrimitive` for deeper coarsening/interpolation contracts.
- `nalgebra`: remains for the current `DVector` and sparse/vector surfaces.

### Sprint Backlog Items

#### Multigrid Smoother Provider Boundary
- [x] **MATH-EUNOMIA-155 [patch]**: Replace smoother thresholds and Chebyshev constants with Eunomia.
- [x] **MATH-EUNOMIA-156 [patch]**: Replace immediate AMG smoother-owner constants and complexity filters with Eunomia.
- [x] **MATH-TEST-157 [patch]**: Add value-semantic multigrid smoother tests.
- [x] **MATH-GATE-158 [patch]**: Verify with rustfmt, cfd-math check, focused nextest, and residue scans.

## Sprint 1.96.69: cfd-math Eunomia Convergence Monitor
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct scalar conversion fallback usage from linear-solver
  convergence helpers.
- Preserve current nalgebra vector/operator trait surfaces for later Leto
  migration.
- Add value-semantic convergence monitor coverage.

### Dependency Audit Findings
- `eunomia::FloatElement`: now owns convergence-factor exponent, CG
  theoretical-bound factor, validation safety-multiplier construction, and
  convergence-factor exponentiation dispatch in `linear_solver/traits.rs`.
- `cfd_core::conversion::SafeFromF64`: no longer appears in
  `linear_solver/traits.rs`.
- `num-traits`: no direct provider residue appears in
  `linear_solver/traits.rs`.
- `nalgebra`: remains for current `DVector` and `LinearOperator` trait
  surfaces.

### Sprint Backlog Items

#### Convergence Monitor Provider Boundary
- [x] **MATH-EUNOMIA-151 [patch]**: Replace convergence-factor scalar conversion and powf dispatch with Eunomia.
- [x] **MATH-EUNOMIA-152 [patch]**: Replace CG bound and validation safety scalar conversions with Eunomia.
- [x] **MATH-TEST-153 [patch]**: Add value-semantic convergence monitor tests.
- [x] **MATH-GATE-154 [patch]**: Verify with rustfmt, cfd-math check, focused nextest, and residue scan.

## Sprint 1.96.68: cfd-math Eunomia Linear Operators
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct `num-traits` scalar conversion from CPU finite-difference
  linear operators.
- Preserve current nalgebra vector operator surfaces for later Leto migration.
- Add value-semantic operator coverage for the touched stencils.

### Dependency Audit Findings
- `eunomia::FloatElement`: now owns scalar constant construction and provider
  bounds in `linear_solver/operators/{poisson.rs,momentum.rs}`.
- `num-traits`: no longer appears in
  `linear_solver/operators/{poisson.rs,momentum.rs}`.
- `nalgebra`: remains for the current `DVector` operator API.

### Sprint Backlog Items

#### CPU Operator Provider Boundary
- [x] **MATH-EUNOMIA-147 [patch]**: Replace Poisson/Laplacian operator scalar constants with Eunomia.
- [x] **MATH-EUNOMIA-148 [patch]**: Replace momentum/energy operator scalar constants with Eunomia.
- [x] **MATH-TEST-149 [patch]**: Add value-semantic center-stencil tests for touched finite-difference operators.
- [x] **MATH-GATE-150 [patch]**: Verify with rustfmt, cfd-math check, focused nextest, and residue scan.

## Sprint 1.96.67: cfd-math Eunomia Schwarz/Cholesky Preconditioners
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct `num-traits` provider residue from Schwarz and
  IncompleteCholesky preconditioners.
- Preserve current nalgebra CSR/vector preconditioner surfaces for later Leto
  migration.

### Dependency Audit Findings
- `eunomia::FloatElement`: now owns IncompleteCholesky symmetry tolerance
  construction.
- `eunomia::NumericElement`: now owns IncompleteCholesky symmetry residual
  absolute-value dispatch.
- `num-traits`: no longer appears in
  `linear_solver/preconditioners/{schwarz.rs,cholesky.rs}`.
- `nalgebra`/`nalgebra_sparse`: remain for the current preconditioner vector
  and CSR matrix surfaces.

### Sprint Backlog Items

#### Schwarz/Cholesky Provider Boundary
- [x] **MATH-EUNOMIA-144 [patch]**: Remove Schwarz's stale direct `FromPrimitive` provider bound.
- [x] **MATH-EUNOMIA-145 [patch]**: Replace Cholesky symmetry tolerance and absolute-value dispatch with Eunomia.
- [x] **MATH-GATE-146 [patch]**: Verify with rustfmt, cfd-math check, focused nextest, and residue scan.

## Sprint 1.96.66: cfd-math Eunomia SSOR Preconditioner
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct `num-traits` scalar conversion from the SSOR preconditioner.
- Preserve the current nalgebra CSR/vector preconditioner surface for later
  Leto migration.

### Dependency Audit Findings
- `eunomia::FloatElement`: now owns SSOR default relaxation scalar
  construction and preconditioner provider bounds.
- `num-traits`: no longer appears in `linear_solver/preconditioners/ssor.rs`.
- `nalgebra`/`nalgebra_sparse`: remain for the current preconditioner vector
  and CSR matrix surfaces.

### Sprint Backlog Items

#### SSOR Preconditioner Provider Boundary
- [x] **MATH-EUNOMIA-141 [patch]**: Replace SSOR default relaxation construction with Eunomia.
- [x] **MATH-EUNOMIA-142 [patch]**: Replace SSOR provider bounds with Eunomia.
- [x] **MATH-GATE-143 [patch]**: Verify with rustfmt, cfd-math check, focused nextest, and residue scan.

## Sprint 1.96.65: cfd-math Eunomia Basic Preconditioners
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct `num-traits` scalar conversions from basic Jacobi/SOR
  preconditioners.
- Preserve the current nalgebra CSR/vector preconditioner surface for later
  Leto migration.

### Dependency Audit Findings
- `eunomia::FloatElement`: now owns Jacobi tolerance and SOR omega scalar
  construction in `linear_solver/preconditioners/basic.rs`.
- `eunomia::NumericElement`: now owns basic preconditioner absolute-value
  dispatch.
- `num-traits`: no longer appears in `linear_solver/preconditioners/basic.rs`.
- `nalgebra`/`nalgebra_sparse`: remain for the current preconditioner vector
  and CSR matrix surfaces.

### Sprint Backlog Items

#### Basic Preconditioner Provider Boundary
- [x] **MATH-EUNOMIA-138 [patch]**: Replace Jacobi diagonal tolerance construction with Eunomia.
- [x] **MATH-EUNOMIA-139 [patch]**: Replace SOR omega construction and validation constants with Eunomia.
- [x] **MATH-GATE-140 [patch]**: Verify with rustfmt, cfd-math check, focused nextest, and residue scan.

## Sprint 1.96.64: cfd-math Eunomia GMRES Chain
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct `num-traits` scalar bounds from the GMRES-centered
  linear-solver call path.
- Preserve the current nalgebra vector/matrix API for later Leto migration.
- Preserve the current rsparse f64-backed direct solver surface until an
  Atlas-owned direct/backend solver replaces it.

### Dependency Audit Findings
- `eunomia::FloatElement`: now owns GMRES, chain, direct solver, and
  block/SIMPLE preconditioner scalar provider bounds in the touched slice.
- `eunomia::NumericElement`: now owns absolute-value and finite conversion
  dispatch in the touched slice.
- `num-traits`: no longer appears in `gmres`, `chain.rs`,
  `direct_solver.rs`, or `block_preconditioner.rs`.
- `nalgebra`/`rsparse`: remain for current vector/matrix and direct-LU backend
  surfaces pending later Atlas provider migration.

### Sprint Backlog Items

#### GMRES Chain Provider Boundary
- [x] **MATH-EUNOMIA-134 [patch]**: Replace GMRES and Arnoldi scalar bounds with Eunomia.
- [x] **MATH-EUNOMIA-135 [patch]**: Replace linear-solver chain and direct solver scalar conversion bounds with Eunomia.
- [x] **MATH-EUNOMIA-136 [patch]**: Replace block/SIMPLE preconditioner scalar safeguards with Eunomia.
- [x] **MATH-GATE-137 [patch]**: Verify with rustfmt, cfd-math check, focused nextest, and residue scan.

## Sprint 1.96.63: cfd-math Eunomia Linear-Solver Config
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct `num-traits` scalar conversion from iterative linear-solver
  default configuration.
- Preserve the current nalgebra `RealField`/`DVector` solver API for later
  Leto migration.
- Keep parallel SpMV provider wording aligned with Moirai.

### Dependency Audit Findings
- `eunomia::FloatElement`: now owns default solver tolerance construction and
  default-construction bounds for CG, BiCGSTAB, and GMRES.
- `num-traits`: no longer appears in `linear_solver/config.rs`,
  `bicgstab/mod.rs`, or `conjugate_gradient/mod.rs`.
- `rayon`: no longer appears in the iterative solver config doc text.
- `nalgebra`: remains for the current linear-solver vector/matrix surface.

### Sprint Backlog Items

#### Linear-Solver Config Provider Boundary
- [x] **MATH-EUNOMIA-131 [patch]**: Replace iterative solver config defaults with Eunomia.
- [x] **MATH-DOC-132 [patch]**: Correct parallel SpMV config wording from Rayon to Moirai.
- [x] **MATH-GATE-133 [patch]**: Verify with touched-file rustfmt, cfd-math check, focused nextest, and residue scan.

## Sprint 1.96.62: cfd-math Eunomia SIMD Scalars
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct `num-traits` scalar conversions from CFD SIMD helpers.
- Preserve existing Moirai-backed parallel slice execution.
- Preserve the current nalgebra `RealField` SIMD API for later Leto/eunomia
  scalar-bound cleanup.

### Dependency Audit Findings
- `eunomia::FloatElement`: now owns CFD SIMD central-difference constants.
- `eunomia::NumericElement`: now owns CFD SIMD field-norm square-root dispatch.
- `num-traits`: no longer appears in `crates/cfd-math/src/simd`.
- `rayon`: no longer appears in `crates/cfd-math/src/simd`; the code uses
  `moirai::prelude::{ParallelSlice, ParallelSliceMut}`.
- `nalgebra`: remains for the current SIMD scalar/vector surface.

### Sprint Backlog Items

#### SIMD Scalar Provider Boundary
- [x] **MATH-EUNOMIA-128 [patch]**: Replace CFD SIMD central-difference constants with Eunomia.
- [x] **MATH-EUNOMIA-129 [patch]**: Replace CFD field-norm square-root dispatch with Eunomia.
- [x] **MATH-GATE-130 [patch]**: Verify with touched-file rustfmt, cfd-math check, focused SIMD nextest, and SIMD residue scan.

## Sprint 1.96.61: cfd-math Eunomia Sparse Scalars
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct `num-traits` scalar conversions and math dispatch from sparse
  matrix utilities.
- Preserve the current `nalgebra_sparse::CsrMatrix` and nalgebra `DVector`
  sparse API for later Leto sparse/vector migration.
- Keep sparse parallel SpMV documentation aligned with the existing Moirai
  slice adapter.

### Dependency Audit Findings
- `eunomia::FloatElement`: now owns sparse stencil constants and singular
  threshold construction.
- `eunomia::NumericElement`: now owns sparse square-root and absolute-value
  dispatch in the touched sparse math.
- `num-traits`: no longer appears in `crates/cfd-math/src/sparse`.
- `rayon`: no longer appears in `crates/cfd-math/src/sparse`; the code already
  uses `moirai::prelude::ParallelSliceMut`.
- `nalgebra_sparse`/`nalgebra`: remain in these modules for the current sparse
  matrix/vector surface.

### Sprint Backlog Items

#### Sparse Scalar Provider Boundary
- [x] **MATH-EUNOMIA-124 [patch]**: Replace sparse pattern constants with Eunomia.
- [x] **MATH-EUNOMIA-125 [patch]**: Replace sparse norm, condition-estimate, and diagonal-dominance scalar math with Eunomia.
- [x] **MATH-DOC-126 [patch]**: Correct sparse SpMV docs from Rayon wording to Moirai wording.
- [x] **MATH-GATE-127 [patch]**: Verify with touched-file rustfmt, cfd-math check, focused sparse nextest, and sparse residue scan.

## Sprint 1.96.60: cfd-math Eunomia Nonlinear Solver Scalars
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct `num-traits` scalar conversions and math dispatch from
  Anderson/JFNK nonlinear solvers.
- Record that this scalar-only slice preserves the existing nalgebra
  `DVector`/`DMatrix` solver API; superseded by Sprint 1.96.109.
- Route QR/JFNK scalar math through Eunomia provider traits.

### Dependency Audit Findings
- `eunomia::FloatElement`: now owns Anderson/JFNK scalar constant construction.
- `eunomia::NumericElement`: now owns square-root, absolute-value, and
  min/max scalar dispatch in the touched nonlinear-solver math.
- `num-traits`: no longer appears in `crates/cfd-math/src/nonlinear_solver`.
- `nalgebra`: vector/matrix solver ownership noted during this slice is closed
  by Sprint 1.96.109.

### Sprint Backlog Items

#### Nonlinear Solver Scalar Provider Boundary
- [x] **MATH-EUNOMIA-121 [patch]**: Replace Anderson defaults and QR scalar math with Eunomia.
- [x] **MATH-EUNOMIA-122 [patch]**: Replace JFNK defaults, perturbation safeguards, EW forcing clamps, Givens rotations, and back-substitution scalar math with Eunomia.
- [x] **MATH-GATE-123 [patch]**: Verify with touched-file rustfmt, cfd-math check, focused nonlinear-solver nextest, and nonlinear-solver residue scan.

## Sprint 1.96.59: cfd-math Eunomia Pressure-Velocity Scalars
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct `num-traits` scalar conversions from SIMPLE pressure-velocity
  default configuration.
- Preserve the existing SIMPLE config/result API while moving the scalar
  contract from nalgebra `RealField` to Eunomia `RealField`.
- Narrow explicit SIMPLE config construction to the scalar operations it
  actually uses.

### Dependency Audit Findings
- `eunomia::RealField`: now owns the pressure-velocity scalar surface.
- `eunomia::FloatElement`: now owns SIMPLE default scalar constant
  construction.
- `num-traits`: no longer appears in `crates/cfd-math/src/pressure_velocity`.
- `nalgebra`: no longer appears in `crates/cfd-math/src/pressure_velocity`.

### Sprint Backlog Items

#### Pressure-Velocity Scalar Provider Boundary
- [x] **MATH-EUNOMIA-118 [patch]**: Replace SIMPLE default constants with Eunomia.
- [x] **MATH-EUNOMIA-119 [patch]**: Narrow explicit SIMPLE config construction bounds.
- [x] **MATH-EUNOMIA-121 [patch]**: Replace pressure-velocity scalar bound and old identities with Eunomia `RealField`.
- [x] **MATH-GATE-120 [patch]**: Verify with cfd-math fmt/check/clippy, focused pressure-velocity nextest, and pressure-velocity residue scan.

## Sprint 1.96.58: cfd-math Eunomia Iterator Scalars
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct `num-traits` scalar conversions from iterator stencil and
  statistics utilities.
- Preserve the existing iterator trait API for the later Leto/eunomia
  scalar-bound cleanup.
- Replace the second-derivative placeholder zero vector with real coefficients
  for every declared 3-point stencil.

### Dependency Audit Findings
- `eunomia::FloatElement`: now owns iterator scalar constant construction and
  count-to-scalar construction.
- `eunomia::NumericElement`: now owns count staging and standard-deviation
  square-root dispatch.
- `num-traits`: no longer appears in `crates/cfd-math/src/iterators`.
- `nalgebra`: remains in this module for the current `RealField` scalar
  surface.

### Sprint Backlog Items

#### Iterator Scalar Provider Boundary
- [x] **MATH-EUNOMIA-114 [patch]**: Replace stencil coefficient constants with Eunomia.
- [x] **MATH-EUNOMIA-115 [patch]**: Replace iterator statistics scalar conversions and sqrt dispatch with Eunomia.
- [x] **MATH-CORRECTNESS-116 [patch]**: Replace second-derivative zero placeholder with real 3-point coefficients.
- [x] **MATH-GATE-117 [patch]**: Verify with touched-file rustfmt, cfd-math check, focused iterator nextest, and iterator residue scan.

## Sprint 1.96.57: cfd-math Eunomia WENO Scalars
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct `num-traits` scalar conversions from WENO5/WENO7 high-order
  reconstruction.
- Preserve the existing WENO reconstruction formulas while routing scalar
  bounds through Eunomia in the follow-up Sprint 1.96.106.
- Keep denominator squaring explicit through multiplication so provider math
  dispatch does not require a competing `powi` method at call sites.

### Dependency Audit Findings
- `eunomia::FloatElement`: now owns WENO scalar constant construction.
- `num-traits`: no longer appears in `crates/cfd-math/src/high_order`.
- `nalgebra`: the earlier scalar-bound residue in this module is closed by
  Sprint 1.96.106.

### Sprint Backlog Items

#### WENO Scalar Provider Boundary
- [x] **MATH-EUNOMIA-111 [patch]**: Replace WENO5 scalar constants and bounds with Eunomia.
- [x] **MATH-EUNOMIA-112 [patch]**: Replace WENO7 scalar constants and bounds with Eunomia.
- [x] **MATH-GATE-113 [patch]**: Verify with touched-file rustfmt, cfd-math check, focused WENO nextest, and high-order residue scan.

## Sprint 1.96.56: cfd-math Eunomia Interpolation Scalars
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct `num-traits` scalar conversions and nalgebra scalar bounds
  from interpolation.
- Preserve the existing interpolation trait APIs while moving their scalar
  contract to Eunomia `RealField`.
- Keep interpolation algorithms unchanged while moving scalar construction and
  identities to Eunomia.

### Dependency Audit Findings
- `eunomia::RealField`: now owns the interpolation scalar surface.
- `eunomia::FloatElement`: now owns cubic-spline scalar constant construction.
- `num-traits`: no longer appears in `crates/cfd-math/src/interpolation`.
- `nalgebra`: no longer appears in `crates/cfd-math/src/interpolation`.

### Sprint Backlog Items

#### Interpolation Scalar Provider Boundary
- [x] **MATH-EUNOMIA-108 [patch]**: Replace cubic-spline Thomas-algorithm constants with Eunomia.
- [x] **MATH-EUNOMIA-109 [patch]**: Replace cubic-spline `FromPrimitive` bounds with Eunomia.
- [x] **MATH-EUNOMIA-111 [patch]**: Replace interpolation trait, linear, and Lagrange scalar bounds with Eunomia `RealField`.
- [x] **MATH-CORRECTNESS-112 [patch]**: Reject duplicate Lagrange nodes before basis denominator division.
- [x] **MATH-GATE-110 [patch]**: Verify with cfd-math fmt/check/clippy, focused interpolation nextest, and residue scan.

## Sprint 1.96.55: cfd-math Eunomia Differentiation Scalars
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct `num-traits` scalar conversions from finite-difference and
  gradient differentiation operators.
- Preserve derivative formulas while moving the public differentiation result
  surfaces to Leto vectors.
- Keep derivative formulas unchanged while moving scalar construction to
  Eunomia.

### Dependency Audit Findings
- `eunomia::FloatElement`: now owns differentiation scalar constant
  construction and SIMD helper scalar staging.
- `eunomia::NumericElement`: now owns zero/one identities and scalar constants
  used by the Leto-backed differentiation result surfaces.
- `num-traits`: no longer appears in `finite_difference.rs` or `gradient.rs`.
- `leto`: now owns `Array1` derivative outputs and `Vector3`
  gradient/divergence/curl surfaces for this module.

### Sprint Backlog Items

#### Differentiation Scalar Provider Boundary
- [x] **MATH-EUNOMIA-104 [patch]**: Replace finite-difference stencil constants with Eunomia.
- [x] **MATH-EUNOMIA-105 [patch]**: Replace gradient, divergence, and curl constants with Eunomia.
- [x] **MATH-EUNOMIA-106 [patch]**: Replace SIMD helper scalar staging with Eunomia.
- [x] **MATH-LETO-108 [patch]**: Replace nalgebra `DVector`/`Vector3`
  differentiation surfaces with Leto `Array1`/`Vector3` and remove the
  type-suffixed SIMD helper name.
- [x] **MATH-GATE-107 [patch]**: Verify with touched-file rustfmt, cfd-math check, focused differentiation nextest, and residue scan.

## Sprint 1.96.54: cfd-math Eunomia Integration Scalars
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct `num-traits` scalar conversions from integration quadrature.
- Preserve the existing `Quadrature` trait APIs while moving their scalar
  contract from nalgebra `RealField` to Eunomia `RealField`.
- Keep quadrature formulas unchanged while moving scalar construction and math
  dispatch to Eunomia.

### Dependency Audit Findings
- `eunomia::FloatElement`: now owns integration scalar constants and
  floating-point math dispatch.
- `eunomia::NumericElement`: now owns exact count conversion staging and
  adaptive absolute-value dispatch.
- `num-traits`: no longer appears in `crates/cfd-math/src/integration`.
- `nalgebra`: no longer appears in `crates/cfd-math/src/integration`.

### Sprint Backlog Items

#### Integration Scalar Provider Boundary
- [x] **MATH-EUNOMIA-100 [patch]**: Replace 1D quadrature constants and weights with Eunomia.
- [x] **MATH-EUNOMIA-101 [patch]**: Replace composite/adaptive/tensor conversion bounds with Eunomia.
- [x] **MATH-EUNOMIA-102 [patch]**: Replace tetrahedral quadrature constants with Eunomia.
- [x] **MATH-EUNOMIA-104 [patch]**: Replace integration trait bounds and old scalar identities with Eunomia `RealField`.
- [x] **MATH-GATE-103 [patch]**: Verify with cfd-math fmt/check/clippy, focused integration nextest, and residue scan.

## Sprint 1.96.53: cfd-math Eunomia Exponential Scalars
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct `num-traits` scalar conversions from exponential time stepping.
- Preserve the existing `ExponentialTimeDifferencing`/`ExponentialRungeKutta4`
  `DVector`/`DMatrix` API for the later Leto matrix/vector migration.
- Keep ETD/ERK formulas unchanged while moving scalar construction to Eunomia.

### Dependency Audit Findings
- `eunomia::FloatElement`: now owns exponential scalar constant construction
  and integer-to-scalar staging for the power-series factorial.
- `num-traits`: no longer appears in `crates/cfd-math/src/time_stepping/exponential.rs`.
- `nalgebra`: remains in this module for the current `DVector`/`DMatrix`
  exponential integrator API.

### Sprint Backlog Items

#### Exponential Scalar Provider Boundary
- [x] **MATH-EUNOMIA-097 [patch]**: Replace ERK4 scalar coefficients and phi-function thresholds with Eunomia.
- [x] **MATH-EUNOMIA-098 [patch]**: Replace scaling/squaring factorial conversion with Eunomia.
- [x] **MATH-GATE-099 [patch]**: Verify with touched-file rustfmt, cfd-math check, focused exponential nextest, and residue scan.

## Sprint 1.96.52: cfd-math Eunomia IMEX Scalars
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct `num-traits` scalar conversions from IMEX time stepping.
- Preserve the existing `IMEXTimeStepper`/`DVector`/`DMatrix` API for the later Leto matrix/vector migration.
- Keep ARS343 numerical coefficients unchanged while moving scalar construction to Eunomia.

### Dependency Audit Findings
- `eunomia::FloatElement`: now owns IMEX scalar constant construction.
- `eunomia::NumericElement`: now owns ARS343 scalar square-root dispatch.
- `num-traits`: no longer appears in `crates/cfd-math/src/time_stepping/imex.rs`.
- `nalgebra`: remains in this module for the current `DVector`/`DMatrix` time-stepper API.

### Sprint Backlog Items

#### IMEX Scalar Provider Boundary
- [x] **MATH-EUNOMIA-094 [patch]**: Replace IMEX Newton tolerance and ARS343 gamma/delta constants with Eunomia.
- [x] **MATH-EUNOMIA-095 [patch]**: Replace IMEX tableau coefficients and solution weights with Eunomia.
- [x] **MATH-GATE-096 [patch]**: Verify with touched-file rustfmt, cfd-math check, focused IMEX nextest, and residue scan.

## Sprint 1.96.51: cfd-math Eunomia RKC Scalars
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct `num-traits` scalar conversions from Runge-Kutta-Chebyshev time stepping.
- Preserve the existing `RhsFunction`/`DVector` API for the later Leto vector migration.
- Keep RKC recurrence and adaptive error-control behavior unchanged while moving scalar construction and math dispatch to Eunomia.

### Dependency Audit Findings
- `eunomia::FloatElement`: now owns RKC scalar constant construction and transcendental dispatch.
- `eunomia::NumericElement`: now owns RKC absolute value, square root, and exact count-to-scalar conversion staging.
- `num-traits`: no longer appears in `crates/cfd-math/src/time_stepping/rk_chebyshev.rs`.
- `nalgebra`: remains in this module for the current `DVector` time-stepper API.

### Sprint Backlog Items

#### RKC Scalar Provider Boundary
- [x] **MATH-EUNOMIA-091 [patch]**: Replace RKC defaults and Chebyshev recurrence constants with Eunomia `FloatElement`.
- [x] **MATH-EUNOMIA-092 [patch]**: Replace RKC adaptive error-control constants, count conversions, and math dispatch with Eunomia.
- [x] **MATH-GATE-093 [patch]**: Verify with touched-file rustfmt, cfd-math check, focused RKC nextest, and residue scan.

## Sprint 1.96.50: cfd-math Eunomia Adaptive Scalars
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct `num-traits` scalar conversions from adaptive time stepping.
- Preserve the existing `TimeStepper`/`DVector` API for the later Leto vector migration.
- Keep Dormand-Prince numerical coefficients unchanged while moving scalar construction to Eunomia.

### Dependency Audit Findings
- `eunomia::FloatElement`: now owns adaptive scalar constant construction.
- `num-traits`: no longer appears in `crates/cfd-math/src/time_stepping/adaptive.rs`.
- `nalgebra`: remains in this module for the current `DVector` time-stepper API.

### Sprint Backlog Items

#### Adaptive Scalar Provider Boundary
- [x] **MATH-EUNOMIA-088 [patch]**: Replace adaptive-stepper defaults and PI controller constants with Eunomia `FloatElement`.
- [x] **MATH-EUNOMIA-089 [patch]**: Replace Dormand-Prince tableau coefficients with Eunomia `FloatElement`.
- [x] **MATH-GATE-090 [patch]**: Verify with touched-file rustfmt, cfd-math check, focused adaptive nextest, and residue scan.

## Sprint 1.96.49: cfd-math Eunomia Runge-Kutta Scalars
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct `num-traits` scalar conversions from cfd-math Runge-Kutta methods.
- Preserve the existing `TimeStepper`/`DVector` API for the later Leto vector migration.
- Correct low-storage RK4 to the Carpenter-Kennedy 2N recurrence.

### Dependency Audit Findings
- `eunomia::FloatElement`: now owns Runge-Kutta scalar constant construction.
- `num-traits`: no longer appears in `crates/cfd-math/src/time_stepping/runge_kutta.rs`.
- `nalgebra`: remains in this module for the current `DVector` time-stepper API.

### Sprint Backlog Items

#### Runge-Kutta Scalar Provider Boundary
- [x] **MATH-EUNOMIA-085 [patch]**: Replace RK3/RK4/low-storage RK4 scalar constants with Eunomia `FloatElement`.
- [x] **MATH-CORRECTNESS-086 [patch]**: Correct low-storage RK4 to the Carpenter-Kennedy 2N residual recurrence.
- [x] **MATH-GATE-087 [patch]**: Verify with touched-file rustfmt, cfd-math check, focused Runge-Kutta nextest, and residue scan.

## Sprint 1.96.48: cfd-math Eunomia Stability Scalars
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Remove direct `num-traits` scalar conversions from cfd-math stability analysis.
- Preserve the existing RK stability formulas and nalgebra tableau surface.
- Keep the Leto matrix migration as a separate bounded item.

### Dependency Audit Findings
- `eunomia::FloatElement`: now owns stability scalar constant construction.
- `eunomia::NumericElement`: now owns stability diagnostic conversion to f64.
- `num-traits`: no longer appears in `crates/cfd-math/src/time_stepping/stability`.
- `nalgebra`: remains in this module for the current Butcher-tableau API.

### Sprint Backlog Items

#### Stability Scalar Provider Boundary
- [x] **MATH-EUNOMIA-082 [patch]**: Replace stability scalar constants with Eunomia `FloatElement`.
- [x] **MATH-EUNOMIA-083 [patch]**: Replace diagnostic `ToPrimitive` conversions with Eunomia `NumericElement`.
- [x] **MATH-GATE-084 [patch]**: Verify with touched-file rustfmt, cfd-math check, focused stability nextest, residue scan, and diff whitespace check.

## Sprint 1.96.46: cfd-2d Eunomia Scheme Complex
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Move the cfd-2d scheme amplification complex value to Eunomia.
- Remove direct `num-complex` manifest ownership from cfd-2d.
- Close the cfd-2d `FloatElement` compile blockers exposed by cfd-core provider bounds.

### Dependency Audit Findings
- `eunomia`: now owns the touched scheme complex value and scalar-provider bounds needed by cfd-core grid/config construction.
- `num-complex`: no longer appears in cfd-2d source or manifest files.
- `num-traits`: remains in cfd-2d for current conversion contracts outside this slice.
- `cfd-3d`/`cfd-validation`: dependency-chain provider blockers are cleared for the cfd-2d nextest gate.

### Sprint Backlog Items

#### cfd-2d Scheme Complex Provider Boundary
- [x] **2D-EUNOMIA-073 [patch]**: Replace `num_complex::Complex<f64>` amplification factors with `eunomia::Complex<f64>`.
- [x] **2D-EUNOMIA-074 [patch]**: Remove direct `num-complex` manifest ownership from cfd-2d.
- [x] **2D-EUNOMIA-075 [patch]**: Propagate explicit `FloatElement` bounds through 2D network/solver construction paths blocked by cfd-core provider requirements.
- [x] **2D-EUNOMIA-076 [patch]**: Verify with cfd-2d/cfd-3d/cfd-validation format and package checks, cfd-2d nextest, and direct complex source/manifest scan.

## Sprint 1.96.47: cfd-3d Apollo/Leto Adapter Compile Repair
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Make the private cfd-3d Leto/Apollo adapter target Apollo's reachable API.
- Remove stale assumptions about unpublished Apollo Leto-native FFT functions.
- Clear cfd-3d/cfd-validation provider-bound blockers from cfd-2d nextest.

### Dependency Audit Findings
- `leto`: owns cfd-3d algorithm storage in the touched spectral/NUFFT paths.
- `apollo-fft`/`apollo-nufft`: still expose ndarray arrays at the reachable transform boundary.
- `eunomia`: now owns scalar-provider bounds required by cfd-3d config and validation solver construction.

### Sprint Backlog Items

#### cfd-3d Adapter and Validation Provider Bounds
- [x] **3D-ATLAS-077 [patch]**: Replace nonexistent Apollo Leto-native FFT imports with reachable Apollo FFT array APIs behind the private adapter.
- [x] **3D-ATLAS-078 [patch]**: Convert between Leto arrays and Apollo ndarray arrays only inside `atlas_array.rs`.
- [x] **3D-ATLAS-079 [patch]**: Remove the stale direct `MnemosyneStorage` import assumption from cfd-3d.
- [x] **VAL-EUNOMIA-080 [patch]**: Propagate `FloatElement` through cfd-validation MMS Richardson solver construction.
- [x] **ATLAS-GATE-081 [patch]**: Verify with cfd-3d/cfd-validation checks and `cargo nextest run -p cfd-2d`.

## Sprint 1.96.45: cfd-math/cfd-validation Eunomia Stability Complex
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Move time-stepping von Neumann complex callback values to Eunomia.
- Remove direct `num-complex` manifest ownership from cfd-math and cfd-validation.
- Preserve explicit RK stability polynomial semantics.

### Dependency Audit Findings
- `eunomia`: now owns the touched stability-analysis complex values.
- `num-complex`: still resolves transitively through nalgebra/simba and remains
  directly owned by cfd-2d outside this slice.
- `nalgebra`: still owns broader dense/sparse linear algebra surfaces pending
  later Leto migration.

### Sprint Backlog Items

#### cfd-math/cfd-validation Stability Complex Provider Boundary
- [x] **MATH-EUNOMIA-069 [patch]**: Replace `num_complex::Complex<f64>` callback types with `eunomia::Complex<f64>`.
- [x] **MATH-EUNOMIA-070 [patch]**: Replace explicit RK dense complex matrix inversion with forward substitution for `(I - zA)x = 1`.
- [x] **MATH-EUNOMIA-071 [patch]**: Remove direct `num-complex` manifest dependencies from cfd-math and cfd-validation.
- [x] **MATH-EUNOMIA-072 [patch]**: Verify with cfd-math check, focused stability nextest, rustfmt, and source/manifest scan.

## Sprint 1.96.44: cfd-core Eunomia Rhie-Chow Interpolation
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Move Rhie-Chow scalar constants to Eunomia.
- Preserve the existing collocated-grid pressure-correction formulas.
- Add value-semantic tests for both velocity-component face interpolators.

### Dependency Audit Findings
- `eunomia`: now owns scalar conversion for touched Rhie-Chow constants.
- `num-traits`: still remains in cfd-core for modules outside this slice.
- `nalgebra`: still owns current fluid-dynamics scalar trait bounds pending later Leto/Eunomia scalar-bound normalization.

### Sprint Backlog Items

#### cfd-core Rhie-Chow Provider Boundary
- [x] **CORE-EUNOMIA-065 [patch]**: Remove direct `num_traits::FromPrimitive` import and bounds from `physics/fluid_dynamics/rhie_chow.rs`.
- [x] **CORE-EUNOMIA-066 [patch]**: Route default relaxation and face-interpolation `2` constants through Eunomia `FloatElement`.
- [x] **CORE-EUNOMIA-067 [patch]**: Add value-semantic u-face and v-face pressure-correction tests.
- [x] **CORE-EUNOMIA-068 [patch]**: Verify with cfd-core no-default/default checks, focused Rhie-Chow nextest, full cfd-core no-default nextest, rustfmt check, and source scan.

## Sprint 1.96.43: cfd-core Eunomia Boundary Geometry
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Move boundary geometry measure scalar constants to Eunomia.
- Preserve existing boundary geometry formulas and public data types.
- Keep conversion-free geometry methods on their narrower existing scalar bounds.

### Dependency Audit Findings
- `eunomia`: now owns scalar conversion for touched sphere and cylinder measure constants.
- `num-traits`: still remains in cfd-core for modules outside this slice.
- `nalgebra`: still owns current boundary geometry point/vector/scalar types pending later Leto/Gaia geometry normalization.

### Sprint Backlog Items

#### cfd-core Boundary Geometry Provider Boundary
- [x] **CORE-EUNOMIA-060 [patch]**: Remove direct `T::from_f64(...).unwrap_or_else(...)` fallback conversions from `physics/boundary/geometry.rs`.
- [x] **CORE-EUNOMIA-061 [patch]**: Route sphere and cylinder measure constants through Eunomia `FloatElement`.
- [x] **CORE-EUNOMIA-062 [patch]**: Split `BoundaryGeometry` impl bounds so only `measure()` requires scalar conversion capability.
- [x] **CORE-EUNOMIA-063 [patch]**: Add value-semantic line, sphere, cylinder, and unsupported-measure tests.
- [x] **CORE-EUNOMIA-064 [patch]**: Verify with cfd-core no-default/default checks, focused boundary geometry nextest, full cfd-core no-default nextest, rustfmt check, and source scan.

## Sprint 1.96.42: cfd-core Eunomia Boundary Ghost Cells
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Move boundary ghost-cell scalar constants to Eunomia.
- Move Robin singularity value reporting off `num_traits::ToPrimitive`.
- Preserve existing high-order ghost-cell formulas and surface degenerate Robin coefficients explicitly.

### Dependency Audit Findings
- `eunomia`: now owns scalar conversion and finite diagnostic conversion for touched ghost-cell formulas.
- `num-traits`: still remains in cfd-core for modules outside this slice.
- `nalgebra`: still owns current boundary scalar trait bounds pending later Leto/Eunomia scalar-bound normalization.

### Sprint Backlog Items

#### cfd-core Boundary Ghost-Cell Provider Boundary
- [x] **CORE-EUNOMIA-055 [patch]**: Remove direct `FromPrimitive`/`ToPrimitive` import and bounds from `physics/boundary/ghost_cells.rs`.
- [x] **CORE-EUNOMIA-056 [patch]**: Route ghost-cell constants and singularity tolerance through Eunomia `FloatElement`.
- [x] **CORE-EUNOMIA-057 [patch]**: Route Robin singularity diagnostic values through Eunomia `NumericElement`.
- [x] **CORE-EUNOMIA-058 [patch]**: Add value-semantic fourth-order Dirichlet, fourth-order Neumann, and degenerate Robin tests.
- [x] **CORE-EUNOMIA-059 [patch]**: Verify with cfd-core no-default/default checks, focused ghost-cell nextest, full cfd-core no-default nextest, rustfmt check, and source scan.

## Sprint 1.96.41: cfd-core Eunomia Staggered Grid
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Move staggered-grid coordinate scalar conversion to Eunomia.
- Preserve uniform and stretched-grid coordinate formulas.
- Surface non-exact integer grid-index conversion as an invariant violation instead of silently rounding.

### Dependency Audit Findings
- `eunomia`: now owns scalar conversion for touched staggered-grid coordinate constants and indices.
- `num-traits`: still remains in cfd-core for modules outside this slice.
- `nalgebra`: still owns current geometry scalar trait bounds pending later Leto/Eunomia scalar-bound normalization.

### Sprint Backlog Items

#### cfd-core Staggered Grid Provider Boundary
- [x] **CORE-EUNOMIA-051 [patch]**: Remove direct `FromPrimitive` import and bounds from `geometry/staggered.rs`.
- [x] **CORE-EUNOMIA-052 [patch]**: Route staggered-grid half constants and grid-index conversions through Eunomia `FloatElement`.
- [x] **CORE-EUNOMIA-053 [patch]**: Assert exact representability for integer grid dimensions and indices before scalar conversion.
- [x] **CORE-EUNOMIA-054 [patch]**: Verify with cfd-core no-default/default checks, focused staggered nextest, full cfd-core no-default nextest, rustfmt check, diff check, and source scan.

## Sprint 1.96.40: cfd-core Eunomia Boundary Time Functions
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Move boundary time-function scalar conversion and math dispatch to Eunomia.
- Move boundary ghost-cell reflection constants to Eunomia.
- Preserve existing boundary application behavior while propagating the provider bound through boundary surfaces.

### Dependency Audit Findings
- `eunomia`: now owns scalar conversion and explicit `sin`/`exp` dispatch for touched boundary time-function paths.
- `num-traits`: still remains in cfd-core for modules outside this slice.
- `nalgebra`: still owns current boundary scalar trait bounds pending later Leto/Eunomia scalar-bound normalization.

### Sprint Backlog Items

#### cfd-core Boundary Provider Boundary
- [x] **CORE-EUNOMIA-046 [patch]**: Replace boundary time-function conversions and math dispatch with Eunomia `FloatElement`.
- [x] **CORE-EUNOMIA-047 [patch]**: Replace boundary ghost-cell reflection constant conversion with Eunomia `FloatElement`.
- [x] **CORE-EUNOMIA-048 [patch]**: Propagate `FloatElement` through boundary applicator/specification/manager surfaces.
- [x] **CORE-EUNOMIA-049 [patch]**: Add value-semantic sinusoidal, exponential, boundary scaling, and ghost-cell tests.
- [x] **CORE-EUNOMIA-050 [patch]**: Verify with cfd-core no-default/default checks, focused boundary nextest, full cfd-core no-default nextest, rustfmt check, and source scan.

## Sprint 1.96.39: cfd-core Eunomia Temperature-Dependent Fluids
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Move temperature-dependent fluid scalar conversion and math dispatch to Eunomia.
- Remove stale conversion bounds from polynomial fluid models.
- Preserve Arrhenius, Andrade, and Sutherland viscosity formulas.

### Dependency Audit Findings
- `eunomia`: now owns scalar exponent conversion and transcendental dispatch for touched temperature-dependent fluid models.
- `num-traits`: still remains in cfd-core for modules outside this slice.
- `nalgebra`: still owns current fluid scalar trait bounds pending later Leto/Eunomia scalar-bound normalization.

### Sprint Backlog Items

#### cfd-core Temperature Fluid Provider Boundary
- [x] **CORE-EUNOMIA-042 [patch]**: Remove direct `FromPrimitive` import and stale conversion bounds from `physics/fluid/temperature.rs`.
- [x] **CORE-EUNOMIA-043 [patch]**: Route Arrhenius, Andrade, and Sutherland math dispatch through Eunomia `FloatElement`.
- [x] **CORE-EUNOMIA-044 [patch]**: Add value-semantic polynomial, Andrade, and Sutherland tests.
- [x] **CORE-EUNOMIA-045 [patch]**: Verify with cfd-core no-default/default checks, focused temperature nextest, full cfd-core no-default nextest, rustfmt check, and source scan.

## Sprint 1.96.38: cfd-core Eunomia Hemolysis Constants
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Move hemolysis calculator and platelet activation scalar conversion to Eunomia.
- Preserve clinical-index formulas for NIH, MIH, and exposure time.
- Correct platelet activation probability to use the decaying exponential.

### Dependency Audit Findings
- `eunomia`: now owns scalar conversion for touched hemolysis constants.
- `num-traits`: still remains in cfd-core for modules outside this slice.
- `nalgebra`: still owns current hemolysis scalar trait bounds pending later Leto/Eunomia scalar-bound normalization.

### Sprint Backlog Items

#### cfd-core Hemolysis Provider Boundary
- [x] **CORE-EUNOMIA-037 [patch]**: Replace hemolysis calculator constants with Eunomia `FloatElement`.
- [x] **CORE-EUNOMIA-038 [patch]**: Replace platelet activation defaults with Eunomia `FloatElement`.
- [x] **CORE-EUNOMIA-039 [patch]**: Correct platelet activation probability to use `1 - exp(-k * excess_stress * exposure_time)`.
- [x] **CORE-EUNOMIA-040 [patch]**: Add value-semantic NIH, MIH, exposure-time, and platelet activation tests.
- [x] **CORE-EUNOMIA-041 [patch]**: Verify with cfd-core no-default/default checks, focused hemolysis nextest, full cfd-core no-default nextest, rustfmt check, and source scan.

## Sprint 1.96.37: cfd-core Eunomia Mesh Quality Thresholds
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Move mesh quality threshold scalar conversion to Eunomia.
- Preserve existing strict threshold comparison semantics.
- Add value-semantic tests for quality levels and recommendations.

### Dependency Audit Findings
- `eunomia`: now owns scalar conversion for mesh quality thresholds.
- `num-traits`: still remains in cfd-core for modules outside this slice.
- `nalgebra`: still owns current mesh quality scalar trait bounds pending later Leto/Eunomia scalar-bound normalization.

### Sprint Backlog Items

#### cfd-core Mesh Quality Provider Boundary
- [x] **CORE-EUNOMIA-033 [patch]**: Replace mesh quality threshold conversions with Eunomia `FloatElement`.
- [x] **CORE-EUNOMIA-034 [patch]**: Remove direct `FromPrimitive` bounds and fallback threshold conversions from `geometry/mesh/quality.rs`.
- [x] **CORE-EUNOMIA-035 [patch]**: Add value-semantic best/worst/boundary mesh quality tests.
- [x] **CORE-EUNOMIA-036 [patch]**: Verify with cfd-core no-default/default checks, focused mesh quality nextest, full cfd-core no-default nextest, rustfmt check, and source scan.

## Sprint 1.96.36: cfd-core Eunomia CPU Backend
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Move CPU backend domain-parameter scalar conversion to Eunomia.
- Remove the local fallback conversion helper from `compute/cpu.rs`.
- Keep conversion bounds off CPU buffer storage operations.

### Dependency Audit Findings
- `eunomia`: now owns scalar conversion for CPU advection domain parameters.
- `num-traits`: still remains in cfd-core for modules outside this slice.
- `nalgebra`: still owns current CPU backend scalar trait bounds pending later Leto/Eunomia scalar-bound normalization.

### Sprint Backlog Items

#### cfd-core CPU Backend Provider Boundary
- [x] **CORE-EUNOMIA-029 [patch]**: Replace CPU advection domain-parameter conversions with Eunomia `FloatElement`.
- [x] **CORE-EUNOMIA-030 [patch]**: Remove the local silent fallback conversion helper.
- [x] **CORE-EUNOMIA-031 [patch]**: Narrow `CpuBuffer` impl bounds away from scalar conversion traits.
- [x] **CORE-EUNOMIA-032 [patch]**: Verify with cfd-core no-default/default checks, focused CPU advection nextest, full cfd-core no-default nextest, rustfmt check, and source scan.

## Sprint 1.96.35: cfd-core Eunomia Time Integrators
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Move cfd-core time integrator scalar constants to Eunomia.
- Remove direct `FromPrimitive` usage from `compute/time/integrators.rs`.
- Remove silent implicit-solver default tolerance conversion fallbacks.

### Dependency Audit Findings
- `eunomia`: now owns scalar constant conversion for time integration schemes that need literal coefficients.
- `num-traits`: still remains in cfd-core for modules outside this slice.
- `nalgebra`: still owns current time integrator vector/state storage pending later Leto-backed replacement.

### Sprint Backlog Items

#### cfd-core Time Integrator Provider Boundary
- [x] **CORE-EUNOMIA-025 [patch]**: Replace RK2/RK4/Crank-Nicolson scalar constants with Eunomia `FloatElement`.
- [x] **CORE-EUNOMIA-026 [patch]**: Replace implicit default tolerance conversion fallbacks with Eunomia constants.
- [x] **CORE-EUNOMIA-027 [patch]**: Add value-semantic explicit and implicit integrator tests.
- [x] **CORE-EUNOMIA-028 [patch]**: Verify with cfd-core no-default/default checks, focused integrator nextest, full cfd-core no-default nextest, rustfmt check, and source scan.

## Sprint 1.96.34: cfd-core Eunomia Time-Step Controllers
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Move cfd-core time-step controller default constants and runtime math helpers to Eunomia.
- Remove direct `FromPrimitive` and `num_traits::Float` usage from `compute/time/controllers.rs`.
- Replace silent invalid-order fallback with typed errors.

### Dependency Audit Findings
- `eunomia`: now owns scalar constant conversion and power functions for time-step controller calculations.
- `num-traits`: still remains in cfd-core for modules outside this slice, including time integrators.
- `nalgebra`: still owns current time-step scalar trait bounds pending later Leto/Eunomia scalar-bound normalization.

### Sprint Backlog Items

#### cfd-core Time Controller Provider Boundary
- [x] **CORE-EUNOMIA-021 [patch]**: Replace controller default constants with Eunomia `FloatElement`.
- [x] **CORE-EUNOMIA-022 [patch]**: Replace runtime `num_traits::Float` math helpers with Eunomia/local helpers.
- [x] **CORE-EUNOMIA-023 [patch]**: Return typed errors for invalid integration order instead of silently falling back.
- [x] **CORE-EUNOMIA-024 [patch]**: Add value-semantic controller tests and verify with cfd-core gates.

## Sprint 1.96.33: cfd-core Eunomia Solver Config Defaults
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Move cfd-core solver configuration default scalar constants to Eunomia.
- Remove direct `FromPrimitive` bounds and silent zero conversion fallbacks from `compute/solver/config.rs`.
- Consolidate duplicated builder/default construction.

### Dependency Audit Findings
- `eunomia`: now owns scalar default conversion for solver configuration constructors and defaults.
- `num-traits`: still remains in cfd-core for modules outside this slice.
- `nalgebra`: still owns current solver scalar trait boundaries pending later Leto-backed replacement.

### Sprint Backlog Items

#### cfd-core Solver Config Provider Boundary
- [x] **CORE-EUNOMIA-017 [patch]**: Replace `SolverConfig::default` direct conversions with Eunomia `FloatElement`.
- [x] **CORE-EUNOMIA-018 [patch]**: Replace `LinearSolverConfig::default` direct conversions with Eunomia `FloatElement`.
- [x] **CORE-EUNOMIA-019 [patch]**: Consolidate `SolverConfigBuilder::new` through `SolverConfig::default`.
- [x] **CORE-EUNOMIA-020 [patch]**: Add value-semantic default/builder parity tests and verify with cfd-core gates.

## Sprint 1.96.32: cfd-core Eunomia Abstraction Defaults
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Move cfd-core abstraction default scalar constants to Eunomia.
- Remove direct `FromPrimitive` bounds and silent default conversion fallbacks from `abstractions/{state,problem}.rs`.

### Dependency Audit Findings
- `eunomia`: now owns scalar default conversion for abstraction constructors and defaults.
- `num-traits`: still remains in cfd-core for modules outside this slice.
- `nalgebra`: still owns current abstraction vector/scalar types pending later Leto-backed state/vector migration.

### Sprint Backlog Items

#### cfd-core Abstraction Default Provider Boundary
- [x] **CORE-EUNOMIA-013 [patch]**: Replace `FieldState::new` default time-step conversion with Eunomia `FloatElement`.
- [x] **CORE-EUNOMIA-014 [patch]**: Replace `ProblemParameters::default` pressure/temperature conversions with Eunomia `FloatElement`.
- [x] **CORE-EUNOMIA-015 [patch]**: Split field-state constructor/default bounds from field accessor/mutator methods.
- [x] **CORE-EUNOMIA-016 [patch]**: Verify with cfd-core no-default/default checks, focused abstraction nextest, full cfd-core no-default nextest, rustfmt check, and source scan.

## Sprint 1.96.31: cfd-core Eunomia Fluid Validation Thresholds
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Move cfd-core fluid validation threshold constants to Eunomia.
- Remove direct `FromPrimitive` bounds and silent threshold conversion fallbacks from `physics/fluid/validation.rs`.

### Dependency Audit Findings
- `eunomia`: now owns scalar threshold conversion for fluid property validation.
- `num-traits`: still remains in cfd-core for modules outside this slice.
- `nalgebra`: still owns the current physics scalar trait boundary pending later Leto-backed state/vector migration.

### Sprint Backlog Items

#### cfd-core Fluid Validation Provider Boundary
- [x] **CORE-EUNOMIA-010 [patch]**: Replace `PropertyBounds::default` direct conversion thresholds with Eunomia `FloatElement`.
- [x] **CORE-EUNOMIA-011 [patch]**: Replace Reynolds/Prandtl/temperature/pressure validator thresholds with Eunomia constants.
- [x] **CORE-EUNOMIA-012 [patch]**: Verify with cfd-core no-default/default checks, focused validation nextest, full cfd-core no-default nextest, rustfmt/diff checks, and source scan.

## Sprint 1.96.30: cfd-core Eunomia Material/Fluid Constants
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Move cfd-core material and immediate fluid constant constructors to Eunomia.
- Remove direct `FromPrimitive` bounds and silent scalar conversion fallbacks from the touched material/fluid constructor cone.

### Dependency Audit Findings
- `eunomia`: now owns scalar constants for touched material, constant-property fluid, fluid database, and ideal-gas constructors.
- `num-traits`: still remains in cfd-core for modules outside this slice.
- `nalgebra`: still owns cfd-core physics scalar/vector trait boundaries pending later Leto-backed state migration.

### Sprint Backlog Items

#### cfd-core Material/Fluid Constant Provider Boundary
- [x] **CORE-EUNOMIA-006 [patch]**: Replace material solid/interface constructor constants with Eunomia `FloatElement`.
- [x] **CORE-EUNOMIA-007 [patch]**: Replace fluid database constructor constants with Eunomia `FloatElement`.
- [x] **CORE-EUNOMIA-008 [patch]**: Replace constant-property fluid and ideal-gas conversion fallbacks with Eunomia constants/math.
- [x] **CORE-EUNOMIA-009 [patch]**: Verify with `cargo check -p cfd-core --no-default-features`, default `cargo check -p cfd-core`, `cargo nextest run -p cfd-core --no-default-features`, touched-file rustfmt/diff checks, and source scans.

## Sprint 1.96.29: cfd-core Eunomia Physics Value Boundary
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Move cfd-core physics value scalar constants to Eunomia.
- Remove direct `FromPrimitive` bounds and silent scalar conversion fallbacks from the touched value/aggregate cone.

### Dependency Audit Findings
- `eunomia`: now owns scalar constant conversion for touched cfd-core value objects and management aggregates.
- `num-traits`: still remains in cfd-core for modules outside this slice, including `management::conversion`.
- `nalgebra`: still owns vector/scalar math for these value objects pending a later Leto-backed value-state slice.

### Sprint Backlog Items

#### cfd-core Physics Value Scalar Provider Boundary
- [x] **CORE-EUNOMIA-001 [patch]**: Add cfd-core `eunomia.workspace = true`.
- [x] **CORE-EUNOMIA-002 [patch]**: Replace value-object `FromPrimitive` bounds with Eunomia `FloatElement`.
- [x] **CORE-EUNOMIA-003 [patch]**: Replace touched aggregate scalar defaults with Eunomia constants and explicit invariant expectations.
- [x] **CORE-EUNOMIA-004 [patch]**: Make custom Reynolds parameter construction return a typed error instead of silently substituting a default.
- [x] **CORE-EUNOMIA-005 [patch]**: Verify with `cargo check -p cfd-core --no-default-features`, `cargo nextest run -p cfd-core --no-default-features`, touched-file rustfmt, and source scans.

## Sprint 1.96.28: cfd-python Eunomia Blood Binding Boundary
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Move cfd-python blood-model scalar conversions to Eunomia.
- Remove cfd-python's direct `num-traits` dependency and silent getter fallback conversions.

### Dependency Audit Findings
- `eunomia`: now owns the numeric conversion trait boundary for blood-model Python getters.
- `num-traits`: removed from direct cfd-python source and manifest usage.
- `numpy`: remains the external Python ABI array boundary; this sprint does not alter NumPy interop.

### Sprint Backlog Items

#### cfd-python Blood Binding Scalar Provider Boundary
- [x] **PY-EUNOMIA-001 [patch]**: Add cfd-python `eunomia.workspace = true` and remove direct `num-traits = "0.2"`.
- [x] **PY-EUNOMIA-002 [patch]**: Replace `FromPrimitive` shear-rate conversions with direct f64 model calls.
- [x] **PY-EUNOMIA-003 [patch]**: Replace `ToPrimitive` getter fallbacks with Eunomia `NumericElement`.
- [x] **PY-EUNOMIA-004 [patch]**: Verify with `cargo check -p cfd-python`, `cargo nextest run -p cfd-python --no-tests pass`, fmt, and source/manifest scans.
- [x] **PY-EUNOMIA-005 [patch]**: Clear dependency-chain clippy diagnostics in cfd-2d, cfd-3d, cfd-validation, and cfd-python so the full `cargo clippy -p cfd-python --all-targets -- -D warnings` gate passes.

## Sprint 1.96.27: cfd-io Eunomia Scalar Boundary
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Move cfd-io checkpoint, binary, and CSV scalar bounds to Eunomia.
- Remove cfd-io's direct `num-traits` dependency and silent numeric fallback conversions.

### Dependency Audit Findings
- `eunomia`: now owns the scalar trait boundary for cfd-io checkpoint validation, binary helpers, and CSV helper types.
- `num-traits`: removed from direct cfd-io source and manifest usage; it still resolves transitively through provider crates.
- `leto`: remains the dense array provider for checkpoint and binary payloads.

### Sprint Backlog Items

#### cfd-io Scalar Provider Boundary
- [x] **IO-EUNOMIA-001 [patch]**: Replace cfd-io `num_traits::Float` bounds with `eunomia::RealField`.
- [x] **IO-EUNOMIA-002 [patch]**: Replace `FromPrimitive`/`ToPrimitive` checkpoint mass-conservation conversions with Eunomia `from_f64` and exact mesh-dimension validation.
- [x] **IO-EUNOMIA-003 [patch]**: Remove direct cfd-io `num-traits` manifest dependency and add direct `eunomia` dependency.
- [x] **IO-EUNOMIA-004 [patch]**: Verify with `cargo check -p cfd-io`, `cargo check -p cfd-io --all-features`, `cargo nextest run -p cfd-io --no-fail-fast`, `cargo clippy -p cfd-io --all-targets --all-features -- -D warnings`, fmt, and source/manifest scans.

## Sprint 1.96.26: cfd-io Leto Checkpoint/Binary Boundary
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Move cfd-io checkpoint and binary dense payloads to Leto arrays.
- Remove cfd-io's direct and default normal-transitive nalgebra dependency path.

### Dependency Audit Findings
- `leto`: now owns checkpoint `Array2` fields and binary `Array1`/`Array2` payload surfaces.
- `nalgebra`: removed from direct cfd-io source/manifest usage and absent from `cargo tree -p cfd-io -e normal`.
- `cfd-core`/`cfd-math`: removed from cfd-io dependencies; cfd-io now owns local file-format errors instead of importing the domain error crate only for I/O failures.
- `num-traits`: removed from cfd-io in Sprint 1.96.27.

### Sprint Backlog Items

#### cfd-io Dense Payload Provider Boundary
- [x] **IO-LETO-001 [arch]**: Replace checkpoint `DMatrix` fields with Leto `Array2` plus explicit row-major serde payloads.
- [x] **IO-LETO-002 [arch]**: Replace binary `DVector`/`DMatrix` helpers with Leto `Array1`/`Array2`.
- [x] **IO-LETO-003 [arch]**: Remove direct cfd-io `nalgebra` dependency and default normal-transitive `cfd-core`/`cfd-math` nalgebra path.
- [x] **IO-LETO-004 [patch]**: Update checkpoint roundtrip tests to assert Leto shape/value/checksum semantics.
- [x] **IO-LETO-005 [patch]**: Verify with `cargo check -p cfd-io`, `cargo nextest run -p cfd-io --no-fail-fast`, `cargo check -p cfd-io --all-features`, fmt, tree, and source scans.

## Sprint 1.96.25: cfd-python Leto PyO3 Boundary
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Move the bounded cfd-python 2D NumPy-return helpers to Leto-owned Rust arrays.
- Keep NumPy as the external Python boundary representation, not the Rust domain storage.

### Dependency Audit Findings
- `leto`: now owns the Rust dense 2D arrays for Poiseuille analytical velocity and Ghia benchmark helper data.
- `numpy`: remains required for Python-visible `numpy.ndarray` results.
- `ndarray`/`nalgebra`: removed as direct `cfd-python` dependencies.
- `num-traits`: still used by cfd-python blood-model bindings and remains a later Eunomia migration item.

### Sprint Backlog Items

#### cfd-python PyO3 Dense Array Boundary
- [x] **PY-LETO-001 [patch]**: Add a private Leto-to-NumPy PyO3 boundary helper for 2D arrays.
- [x] **PY-LETO-002 [patch]**: Replace Poiseuille analytical velocity `DMatrix`/`ndarray` construction with Leto `Array2`.
- [x] **PY-LETO-003 [patch]**: Replace cavity Ghia benchmark `ndarray` construction with Leto `Array2`.
- [x] **PY-LETO-004 [patch]**: Remove direct `ndarray` and `nalgebra` cfd-python dependencies.
- [x] **PY-LETO-005 [patch]**: Verify with `cargo check -p cfd-python`, `cargo nextest run -p cfd-python --no-tests pass`, `cargo fmt -p cfd-python --check`, and source/dependency scans.

## Sprint 1.96.24: CFDrs Moirai GPU Boundary
**Status**: Completed
**Start Date**: July 3, 2026

### Sprint Objectives
- Move direct `pollster` GPU wait/probe boundaries in `cfd-core` and `cfd-math` to the Atlas Moirai runtime.
- Preserve the current raw WGPU kernel surface until the Hephaestus WGPU version boundary can be normalized.

### Dependency Audit Findings
- `moirai`: now owns the synchronous blocking boundary for `cfd-core` GPU context creation, support detection, Poisson residual readback, and `cfd-math` GPU operator dispatch.
- `pollster`: removed from the workspace dependency graph.
- `hephaestus-wgpu`: is the target GPU provider, but direct replacement is blocked in this slice because CFDrs uses `wgpu 0.19` and Hephaestus uses `wgpu 26.0`.

### Sprint Backlog Items

#### cfd-core GPU Blocking Boundary
- [x] **GPU-MOIRAI-001 [patch]**: Replace direct `pollster::block_on` calls in `cfd-core` with `moirai::block_on`.
- [x] **GPU-MOIRAI-002 [patch]**: Replace direct `pollster::block_on` calls in `cfd-math` with `moirai::block_on`.
- [x] **GPU-MOIRAI-003 [patch]**: Remove `pollster` from `cfd-core`, `cfd-math`, and the workspace dependency graph.
- [x] **GPU-MOIRAI-004 [patch]**: Record the Hephaestus WGPU version-normalization blocker.
- [x] **GPU-MOIRAI-005 [patch]**: Verify with `cargo check -p cfd-core --features gpu`, `cargo nextest run -p cfd-core --features gpu --lib`, `cargo check -p cfd-math --features gpu`, and `cargo nextest run -p cfd-math --features gpu --lib`.

## Sprint 1.96.23: cfd-3d Atlas Array Migration
**Status**: Completed
**Start Date**: July 2, 2026

### Sprint Objectives
- Move the next bounded `cfd-3d` spectral/NUFFT dense-array path to Atlas-owned Leto arrays.
- Keep current Apollo FFT/NUFFT usage provider-owned without relying on unpublished local Apollo commits.

### Dependency Audit Findings
- `leto`: now owns `cfd-3d` spectral DNS, forcing, diagnostics, Fourier wrapper, and IBM NUFFT dense 1D/3D array storage in this slice.
- `apollo-fft` and `apollo-nufft`: remain the FFT/NUFFT providers; their reachable published APIs still use `ndarray`, so conversion is isolated in `crates/cfd-3d/src/atlas_array.rs`.
- `ndarray`: remains a boundary-only compatibility dependency for Apollo in this slice, not algorithm-owned dense storage in the migrated paths.

### Sprint Backlog Items

#### cfd-3d Leto Dense-Array Boundary
- [x] **ATLAS-ARRAY-001 [patch]**: Add a private Leto/Apollo adapter for current Apollo FFT and NUFFT APIs.
- [x] **ATLAS-ARRAY-002 [patch]**: Migrate spectral DNS, forcing, diagnostics, and Fourier wrapper dense arrays to Leto.
- [x] **ATLAS-ARRAY-003 [patch]**: Migrate IBM NUFFT scalar/velocity sampling and force spreading dense arrays to Leto.
- [x] **ATLAS-ARRAY-004 [patch]**: Verify with `cargo check -p cfd-3d --no-default-features` and `cargo nextest run -p cfd-3d --no-default-features --lib`.

## Sprint 1.96.22: cfd-1d Dependency-Aware Physics Audit
**Status**: Completed
**Start Date**: May 12, 2026

### Sprint Objectives
- Audit `cfd-1d` and the crates it uses for network flow analysis, hydraulic diagnostics, topology conversion, and numerical solving.
- Correct the next highest-risk local physics defect found during the audit without changing solver state semantics.

### Dependency Audit Findings
- `cfd-core`: appropriate for fluid properties used in Reynolds and shear diagnostics; the analysis layer must preserve their physical definitions.
- `cfd-math`: appropriate for reusable solver infrastructure; post-solve flow-regime classification remains local to `cfd-1d` because it depends on network edge state and reduced-order hydraulic properties.
- `cfd-schematics`: appropriate as topology and channel-geometry input authority; `cfd-1d` consumes converted edge area and hydraulic diameter for diagnostics.
- External crates: `petgraph` remains appropriate for graph topology traversal, `nalgebra` for scalar algebra bounds, `rayon` for independent analysis paths, and `serde`/`serde_json` for persisted scenarios.

### Sprint Backlog Items

#### Flow-Analysis Reynolds Orientation
- [x] **FLOW-RE-001 [patch]**: Classify flow regimes with Reynolds number computed from velocity magnitude, not signed velocity.
- [x] **FLOW-RE-002 [patch]**: Add real-network tests for reverse-flow transitional classification.
- [x] **FLOW-RE-003 [patch]**: Add forward/reverse reciprocity tests for scalar diagnostics: Reynolds number, shear rate, and flow regime.
- [x] **FLOW-RE-004 [patch]**: Verify the touched `cfd-1d` flow-analysis path with bounded Cargo check, unit test, clippy, and nextest runs.

## Sprint 1.96.21: cfd-1d Dependency-Aware Physics Audit
**Status**: Completed
**Start Date**: May 12, 2026

### Sprint Objectives
- Audit `cfd-1d` and the crates it uses for one-dimensional hydraulics, two-phase regime classification, rheology, topology conversion, and numerical solving.
- Correct the next highest-risk local physics defect found during the audit without leaving compatibility wrappers or fake finite states.

### Dependency Audit Findings
- `cfd-core`: appropriate for typed physics errors and shared fluid contracts; droplet-regime dimensionless groups should use these errors when dimensional inputs violate physical domains.
- `cfd-math`: appropriate for reusable numerical kernels; capillary, Weber, and Ohnesorge group definitions remain local to `cfd-1d` two-phase millifluidic physics.
- `cfd-schematics`: appropriate as topology and geometry-input authority; `cfd-1d` consumes characteristic length and velocity values for regime classification.
- External crates: `petgraph` remains appropriate for network topology, `nalgebra`/`nalgebra-sparse`/`sprs` for algebra and sparse storage, `rayon` for independent branch analysis, and `serde`/`serde_json` for scenario persistence.

### Sprint Backlog Items

#### Droplet-Regime Physical Domains
- [x] **DROP-DOMAIN-001 [patch]**: Reject nonpositive or nonfinite surface tension instead of clamping denominators in capillary and Weber numbers.
- [x] **DROP-DOMAIN-002 [patch]**: Reject nonpositive/nonfinite characteristic length and density where required by Weber and Ohnesorge numbers.
- [x] **DROP-DOMAIN-003 [patch]**: Reject nonfinite or negative capillary numbers before flow-regime classification.
- [x] **DROP-DOMAIN-004 [patch]**: Add value-semantic tests for invalid surface tension, invalid length, and invalid capillary-number classification.
- [x] **DROP-DOMAIN-005 [patch]**: Verify the touched `cfd-1d` droplet-regime path with bounded Cargo check, unit test, clippy, and nextest runs.

## Sprint 1.96.20: cfd-1d Dependency-Aware Physics Audit
**Status**: Completed
**Start Date**: May 12, 2026

### Sprint Objectives
- Audit `cfd-1d` and the crates it uses for hydraulic-network physics, non-Newtonian rheology, schematic topology conversion, and numerical solving.
- Correct the next highest-risk local physics defect found during the audit without changing public APIs.

### Dependency Audit Findings
- `cfd-core`: appropriate for shared fluid properties, Casson blood rheology, physical constants, and typed physics errors; junction-loss viscosity must use these rheology contracts with a nonnegative shear-rate magnitude.
- `cfd-math`: appropriate for reusable solver and algebra kernels; K-factor junction constitutive physics remains local to `cfd-1d` because it depends on 1D branch area, hydraulic diameter, and minor-loss coefficients.
- `cfd-schematics`: appropriate as the source of channel topology and geometry inputs; `cfd-1d` owns conversion into branch hydraulic properties and junction-loss state.
- External crates: `petgraph` remains appropriate for graph topology, `nalgebra`/`nalgebra-sparse`/`sprs` for algebra and sparse storage, `rayon` for independent branch analysis, and `serde`/`serde_json` for scenario persistence.

### Sprint Backlog Items

#### Junction-Loss Shear-Rate Physics
- [x] **JUNC-SHEAR-001 [patch]**: Derive junction rheology shear rate from velocity magnitude instead of signed velocity.
- [x] **JUNC-SHEAR-002 [patch]**: Use explicit nonnegative `FlowConditions::shear_rate` for junction non-Newtonian viscosity when supplied.
- [x] **JUNC-SHEAR-003 [patch]**: Reject negative explicit junction wall shear rates before rheology lookup.
- [x] **JUNC-SHEAR-004 [patch]**: Add Casson blood tests for explicit shear-rate dependence and reverse-flow resistance reciprocity.
- [x] **JUNC-SHEAR-005 [patch]**: Verify the touched `cfd-1d` junction-loss path with bounded Cargo check, unit test, clippy, and nextest runs.

## Sprint 1.96.19: cfd-1d Dependency-Aware Physics Audit
**Status**: Completed
**Start Date**: May 11, 2026

### Sprint Objectives
- Audit `cfd-1d` and the crates it uses for one-dimensional hydraulic-network physics, rheology, topology conversion, and numerical solving.
- Correct the next highest-risk local physics defect found during the audit without changing public APIs.

### Dependency Audit Findings
- `cfd-core`: appropriate for shared Newtonian and non-Newtonian fluid contracts, Casson blood rheology, physical constants, and typed errors; Darcy-Weisbach resistance must propagate these rheology results instead of substituting baseline viscosity.
- `cfd-math`: appropriate for reusable numerical kernels and solver infrastructure; friction-factor and resistance-model equations remain local because they are one-dimensional hydraulic constitutive laws.
- `cfd-schematics`: appropriate as the source of topology and channel geometry inputs; `cfd-1d` owns conversion into hydraulic diameter, area, and resistance state.
- External crates: `petgraph` remains appropriate for network topology, `nalgebra`/`nalgebra-sparse`/`sprs` for algebra and sparse storage, `rayon` for independent branch analysis, and `serde`/`serde_json` for persisted scenarios.

### Sprint Backlog Items

#### Darcy-Weisbach Shear-Rate Physics
- [x] **DW-SHEAR-001 [patch]**: Use explicit nonnegative `FlowConditions::shear_rate` for Darcy-Weisbach non-Newtonian viscosity when supplied.
- [x] **DW-SHEAR-002 [patch]**: Reject negative explicit Darcy-Weisbach wall shear rates before rheology lookup.
- [x] **DW-SHEAR-003 [patch]**: Propagate `cfd-core` rheology errors instead of replacing them with baseline fluid viscosity.
- [x] **DW-SHEAR-004 [patch]**: Add Casson blood tests for explicit shear-rate dependence and reverse-flow auto-Reynolds reciprocity.
- [x] **DW-SHEAR-005 [patch]**: Verify the touched `cfd-1d` Darcy-Weisbach path with bounded Cargo check, unit test, clippy, and nextest runs.

## Sprint 1.96.18: cfd-1d Dependency-Aware Physics Audit
**Status**: Completed
**Start Date**: May 11, 2026

### Sprint Objectives
- Audit `cfd-1d` and the crates it uses for hydraulic-network physics, rheology, schematic topology conversion, and numerical solving.
- Correct the next highest-risk local physics defect found during the audit without changing public APIs.

### Dependency Audit Findings
- `cfd-core`: appropriate for shared fluid properties, Casson/non-Newtonian blood rheology, physical constants, and error contracts; `cfd-1d` should call these rheology contracts rather than duplicating viscosity laws.
- `cfd-math`: appropriate for reusable linear algebra and solver support; network assembly and hydraulic resistance model selection remain local because they encode one-dimensional circuit physics.
- `cfd-schematics`: appropriate as the topology and cross-section source of truth for millifluidic channels; `cfd-1d` converts schematics into network elements rather than owning schematic geometry.
- External crates: `petgraph` is appropriate for graph topology, `nalgebra`/`nalgebra-sparse`/`sprs` for matrix assembly and sparse storage, `rayon` for independent branch calculations, and `serde` for persisted network state.

### Sprint Backlog Items

#### Rectangular-Channel Shear-Rate Physics
- [x] **RECT-SHEAR-001 [patch]**: Use explicit nonnegative `FlowConditions::shear_rate` for rectangular-channel non-Newtonian viscosity when supplied.
- [x] **RECT-SHEAR-002 [patch]**: Reject negative explicit rectangular-channel wall shear rates before rheology lookup.
- [x] **RECT-SHEAR-003 [patch]**: Add Casson blood tests for explicit shear-rate dependence and reverse-flow resistance reciprocity.
- [x] **RECT-SHEAR-004 [patch]**: Verify the touched `cfd-1d` rectangular resistance path with bounded Cargo check, unit test, clippy, and nextest runs.

## Sprint 1.96.17: cfd-3d Dependency-Aware Physics Audit
**Status**: Completed
**Start Date**: May 11, 2026

### Sprint Objectives
- Audit `cfd-3d` and its direct dependencies (`cfd-core`, `cfd-math`, `cfd-mesh`, `cfd-io`, `cfd-1d`, `cfd-2d`, `cfd-schematics`, and numerical support crates) for role alignment in 3D FEM, IBM, level-set, VOF, spectral, and domain solvers.
- Correct the next highest-risk local physics defect found during the audit without changing public APIs.

### Dependency Audit Findings
- `cfd-core`: appropriate for shared fluid property types, boundary-condition contracts, and errors; `cfd-3d` must validate these inputs before FEM assembly consumes them.
- `cfd-math`: appropriate for reusable sparse and iterative solvers; Stokes problem validation remains local because it owns FEM pressure-space metadata and element viscosity fields.
- `cfd-mesh`: appropriate as mesh topology and geometry authority; `cfd-3d` correctly consumes `IndexedMesh` while validating solver-specific pressure DOFs and per-cell fields.
- `cfd-io`: appropriate for output and checkpoints, outside FEM assembly physics.
- `cfd-1d` and `cfd-2d`: appropriate as lower-fidelity references and validation baselines, not replacements for 3D FEM well-posedness checks.
- `cfd-schematics`: appropriate as topology input authority for blueprint-driven preprocessing.
- External crates: `apollofft`, `ndarray`, `nalgebra`, `nalgebra-sparse`, `rayon`, `crossbeam`, and `serde` remain appropriate for spectral transforms, dense fields, algebra, sparse systems, parallel execution, and persisted configuration.

### Sprint Backlog Items

#### FEM Stokes Problem Physical Invariants
- [x] **FEM-INV-001 [patch]**: Reject nonpositive or nonfinite fluid density and dynamic viscosity before Stokes assembly.
- [x] **FEM-INV-002 [patch]**: Reject invalid pressure corner-node counts.
- [x] **FEM-INV-003 [patch]**: Reject nonfinite body-force and boundary-condition scalar/vector values.
- [x] **FEM-INV-004 [patch]**: Reject per-element viscosity fields with incorrect length or nonpositive/nonfinite values.
- [x] **FEM-INV-005 [patch]**: Verify the touched `cfd-3d` FEM problem module with bounded Cargo check, unit test, clippy, and nextest runs.

## Sprint 1.96.16: cfd-3d Dependency-Aware Physics Audit
**Status**: Completed
**Start Date**: May 11, 2026

### Sprint Objectives
- Audit `cfd-3d` and its direct dependencies (`cfd-core`, `cfd-math`, `cfd-mesh`, `cfd-io`, `cfd-1d`, `cfd-2d`, `cfd-schematics`, and numerical support crates) for role alignment in 3D FEM, IBM, level-set, VOF, spectral, and domain solvers.
- Correct the next highest-risk local physics defect found during the audit without changing public APIs.

### Dependency Audit Findings
- `cfd-core`: appropriate for shared error contracts, fluid/cavitation physics, and boundary concepts; `cfd-3d` should reuse these contracts for transport validation.
- `cfd-math`: appropriate for reusable sparse, spectral, and high-order numerical kernels; level-set Hamilton-Jacobi transport remains local because it depends on 3D grid layout and WENO stencil ownership.
- `cfd-mesh`: appropriate for mesh/CSG authority and scalar trait bounds; it is not a substitute for 3D level-set field validation.
- `cfd-io`: appropriate for visualization and checkpoint output outside the transport update path.
- `cfd-1d` and `cfd-2d`: appropriate as cross-fidelity references and validation baselines, not replacements for 3D level-set interface transport.
- `cfd-schematics`: appropriate as topology input authority for blueprint-driven preprocessing.
- External crates: `apollofft`, `ndarray`, `nalgebra`, `nalgebra-sparse`, `rayon`, `crossbeam`, and `serde` are appropriate for spectral transforms, dense fields, algebra, sparse systems, parallel execution, and persisted configuration.

### Sprint Backlog Items

#### Level-Set Hamilton-Jacobi Transport Preconditions
- [x] **LS-PRE-001 [patch]**: Reject nonpositive and nonfinite level-set time steps before WENO/SSPRK3 or upwind transport.
- [x] **LS-PRE-002 [patch]**: Reject nonpositive/nonfinite grid spacing before Hamilton-Jacobi transport.
- [x] **LS-PRE-003 [patch]**: Reject nonfinite level-set velocity components before derivative reconstruction.
- [x] **LS-PRE-004 [patch]**: Add value-semantic tests for NaN velocity, zero time step, and negative time step rejection.
- [x] **LS-PRE-005 [patch]**: Verify the touched `cfd-3d` level-set module with bounded Cargo check, integration test, clippy, and nextest runs.

## Sprint 1.96.15: cfd-3d Dependency-Aware Physics Audit
**Status**: Completed
**Start Date**: May 11, 2026

### Sprint Objectives
- Audit `cfd-3d` and its direct dependencies (`cfd-core`, `cfd-math`, `cfd-mesh`, `cfd-io`, `cfd-1d`, `cfd-2d`, `cfd-schematics`, and numerical support crates) for role alignment in 3D FEM, IBM, level-set, VOF, spectral, and domain solvers.
- Correct the highest-risk local physics defect found during the audit without changing public APIs.

### Dependency Audit Findings
- `cfd-core`: appropriate for shared fluid properties, cavitation physics, boundary conditions, errors, and physical constants; `cfd-3d` should consume these contracts instead of duplicating material models.
- `cfd-math`: appropriate for sparse solvers, spectral operators, and reusable numerical kernels; 3D solver assembly remains local where element topology or grid storage is 3D-specific.
- `cfd-mesh`: appropriate for unstructured geometry, CSG, and manufacturing mesh inputs; FEM and domain solvers correctly treat it as mesh authority rather than duplicating mesh topology.
- `cfd-io`: appropriate for visualization and checkpoint output; it is not part of the time-integration physics path.
- `cfd-1d` and `cfd-2d`: appropriate for cross-fidelity references, reduced-order seeding, and shared turbulence validation, not substitutes for 3D field solves.
- `cfd-schematics`: appropriate as millifluidic topology/layout authority for blueprint-driven preprocessing.
- External crates: `apollofft` is appropriate for spectral transforms, `ndarray` for dense field data, `nalgebra` and `nalgebra-sparse` for vector/matrix algebra, `rayon` and `crossbeam` for independent field/element work, and `serde` for persisted configurations.

### Sprint Backlog Items

#### VOF Explicit Transport Preconditions
- [x] **VOF-PRE-001 [patch]**: Reject nonpositive and nonfinite VOF time steps before explicit interface transport.
- [x] **VOF-PRE-002 [patch]**: Reject nonfinite VOF velocity components before CFL evaluation and geometric flux calculation.
- [x] **VOF-PRE-003 [patch]**: Add value-semantic tests for NaN velocity, infinite velocity, zero time step, and negative time step rejection.
- [x] **VOF-PRE-004 [patch]**: Verify the touched `cfd-3d` VOF module with bounded Cargo check, integration test, clippy, and nextest runs.

## Sprint 1.96.14: cfd-1d Dependency-Aware Physics Audit
**Status**: Completed
**Start Date**: May 11, 2026

### Sprint Objectives
- Audit `cfd-1d` and its direct internal dependencies (`cfd-core`, `cfd-math`, `cfd-schematics`) for role alignment in 1D millifluidic network solves.
- Correct the highest-risk local physics defect found during the audit without changing public APIs.

### Dependency Audit Findings
- `cfd-core`: appropriate for fluid properties, non-Newtonian blood rheology, boundary-condition concepts, shared errors, and physical constants; `cfd-1d` should consume these contracts instead of duplicating rheology.
- `cfd-math`: appropriate for sparse and iterative numerical kernels in network solves; local graph-to-matrix assembly remains in `cfd-1d` because it owns Kirchhoff topology and boundary elimination.
- `cfd-schematics`: appropriate as topology and geometry input authority; `cfd-1d` converts schematics into hydraulic networks and does not own schematic rendering or layout semantics.
- External crates: `petgraph` is appropriate for directed network topology, `nalgebra` and `nalgebra-sparse` for scalar/vector and sparse matrix operations, `sprs` for sparse storage compatibility, `serde` for persisted design state, and `rayon` only for independent network analyses.

### Sprint Backlog Items

#### Hagen-Poiseuille Non-Newtonian Shear Magnitude
- [x] **HP-SHEAR-001 [patch]**: Derive wall shear rate from velocity or flow-rate magnitude for straight circular channels.
- [x] **HP-SHEAR-002 [patch]**: Reject explicitly negative wall shear-rate inputs because shear rate is a scalar magnitude.
- [x] **HP-SHEAR-003 [patch]**: Add value-semantic Casson blood regression tests proving reverse-flow resistance reciprocity.
- [x] **HP-SHEAR-004 [patch]**: Verify the touched `cfd-1d` Hagen-Poiseuille module with bounded Cargo check, unit test, clippy, and nextest runs.

## Sprint 1.96.13: cfd-2d Dependency-Aware Physics Audit
**Status**: Completed
**Start Date**: May 4, 2026

### Sprint Objectives
- Audit `cfd-2d` and its direct internal dependencies (`cfd-core`, `cfd-math`, `cfd-mesh`, `cfd-io`, `cfd-1d`, `cfd-schematics`) for role alignment in 2D Navier-Stokes, LBM, and schematics-projected channel solves.
- Correct the highest-risk local physics defect found during the audit without changing the public solver API.

### Dependency Audit Findings
- `cfd-core`: appropriate for shared fluid physics, boundary condition types, errors, constants, and GPU feature plumbing; `cfd-2d` should not duplicate those foundation contracts.
- `cfd-math`: appropriate for sparse/iterative numerical kernels and high-order schemes; local discretization remains in `cfd-2d` where it depends on 2D stencil geometry.
- `cfd-mesh`: appropriate for mesh-backed geometry inputs and unstructured contexts; structured solver kernels correctly keep 2D grid storage local.
- `cfd-io`: appropriate for field output and visualization export; it is not part of the time-integration physics path.
- `cfd-1d`: appropriate as a reduced-order pressure/flow reference for channel projection and cross-fidelity seeding, not as a replacement for 2D PDE solves.
- `cfd-schematics`: appropriate as topology and layout authority; `cfd-2d` consumes blueprints through projection metadata instead of owning schematic generation.

### Sprint Backlog Items

#### 2D LBM Low-Mach Physics
- [x] **LBM-MACH-001 [patch]**: Enforce `Ma <= 0.1` during LBM initialization.
- [x] **LBM-MACH-002 [patch]**: Enforce `Ma <= 0.1` for Zou-He velocity and pressure boundary reconstruction.
- [x] **LBM-MACH-003 [patch]**: Propagate boundary-condition Mach violations through `LbmSolver::step`.
- [x] **LBM-MACH-004 [patch]**: Add value-semantic tests for high-Mach initialization and inlet rejection.
- [x] **LBM-MACH-005 [patch]**: Verify the touched `cfd-2d` LBM module with bounded Cargo check, unit test, clippy, and nextest runs.

## Sprint 1.96.12: cfd-1d Dependency-Aware Physics Audit
**Status**: Completed
**Start Date**: May 4, 2026

### Sprint Objectives
- Audit `cfd-1d` and its direct internal dependencies (`cfd-core`, `cfd-math`, `cfd-schematics`) for role alignment in the 1D millifluidic solver.
- Correct the highest-risk local physics defect found during the audit without changing public APIs.

### Dependency Audit Findings
- `cfd-core`: appropriate as the source for fluid properties, constants, boundary conditions, shared errors, and solver traits; `cfd-1d` uses it at physics and solver boundaries rather than duplicating fluid databases.
- `cfd-math`: appropriate for sparse and iterative linear algebra in network solves; `cfd-1d` also uses `nalgebra-sparse` directly for assembly, which is acceptable because it owns the graph-to-matrix mapping.
- `cfd-schematics`: appropriate as topology and cross-section authority; `cfd-1d` consumes `NodeSpec`, `ChannelSpec`, and render metadata but keeps hydraulic resistance and flow-state logic local.
- External crates: `petgraph` is appropriate for directed network topology, `nalgebra` for scalar/vector algebra, `serde` for persisted design state, `rayon` only where independent network paths can be parallelized.

### Sprint Backlog Items

#### 1D Serpentine Reverse-Flow Physics
- [x] **SERP-REV-001 [patch]**: Use flow-speed magnitude for serpentine shear rate, Reynolds number, Dean number, friction losses, and bend minor losses.
- [x] **SERP-REV-002 [patch]**: Preserve orientation outside scalar resistance coefficients by keeping resistance magnitudes invariant under `Q -> -Q`.
- [x] **SERP-REV-003 [patch]**: Add value-semantic tests for reverse-flow coefficient and analysis invariance.
- [x] **SERP-REV-004 [patch]**: Verify the touched `cfd-1d` serpentine module with bounded Cargo check, unit test, clippy, and nextest runs.

## Sprint 1.96.11: cfd-3d Sigma SGS Energy Physics
**Status**: Completed
**Start Date**: May 4, 2026

### Sprint Objectives
- Correct `cfd-3d` Sigma LES SGS kinetic-energy diagnostics so they use the same Yoshizawa eddy-viscosity inversion as the other algebraic LES closures.
- Consolidate SGS energy computation through the shared turbulence helper instead of retaining a Sigma-specific strain-rate expression with incorrect dimensions.

### Sprint Backlog Items

#### 3D LES Closure Physics
- [x] **SIGMA-SGS-001 [patch]**: Replace `nu_t |S| / Delta` Sigma SGS energy with `(nu_t / (C_k Delta))^2`.
- [x] **SIGMA-SGS-002 [patch]**: Route Sigma SGS energy through the shared `kinetic_energy_from_eddy_viscosity` SSOT.
- [x] **SIGMA-SGS-003 [patch]**: Add value-semantic regression coverage for the Yoshizawa relation and the former dimensional mismatch.
- [x] **SIGMA-SGS-004 [patch]**: Verify the touched `cfd-3d` Sigma module with bounded Cargo check, unit test, clippy, and nextest runs.

## Sprint 1.96.10: cfd-3d VOF Directional CFL Physics
**Status**: Completed
**Start Date**: May 4, 2026

### Sprint Objectives
- Correct `cfd-3d` VOF timestep selection so it targets the same summed directional CFL invariant enforced by geometric advection.
- Add value-semantic regression coverage for diagonal flow that distinguishes the corrected reciprocal summed-rate bound from the former norm/min-spacing estimate.

### Sprint Backlog Items

#### 3D VOF Stability Physics
- [x] **VOF-CFL-001 [patch]**: Replace speed/min-spacing timestep selection with `C / max(|u_x|/dx + |u_y|/dy + |u_z|/dz)`.
- [x] **VOF-CFL-002 [patch]**: Document the explicit VOF CFL theorem at the solver timestep API.
- [x] **VOF-CFL-003 [patch]**: Add regression coverage that accepts the corrected diagonal-flow timestep and rejects the former estimate.
- [x] **VOF-CFL-004 [patch]**: Verify the touched `cfd-3d` VOF module with bounded Cargo check, integration test, and nextest runs.

## Sprint 1.96.9: cfd-2d Upwind Coefficient Orientation
**Status**: Completed
**Start Date**: May 4, 2026

### Sprint Objectives
- Correct `cfd-2d` finite-volume first-order upwind coefficients so west-face convection uses the Patankar face-orientation sign convention.
- Add value-semantic tests that distinguish positive and negative west-face flux orientation.

### Sprint Backlog Items

#### 2D Convection Discretization Physics
- [x] **CONV-001 [patch]**: Replace `a_W = D_W + max(-F_W, 0)` with `a_W = D_W + max(F_W, 0)` in `FirstOrderUpwind`.
- [x] **CONV-002 [patch]**: Document the finite-volume upwind coefficient invariant for east and west faces.
- [x] **CONV-003 [patch]**: Add regression coverage for west-face positive and negative flux orientation.
- [x] **CONV-004 [patch]**: Verify the touched `cfd-2d` convection module with bounded Cargo check, unit test, and nextest runs.

## Sprint 1.96.8: cfd-1d Branch Reverse-Flow Physics
**Status**: Completed
**Start Date**: May 4, 2026

### Sprint Objectives
- Correct `cfd-1d` branch-junction solvers so reversed parent flow preserves daughter orientation instead of being rejected.
- Keep wall-shear and non-Newtonian viscosity inputs magnitude-based while pressure drops remain signed.

### Sprint Backlog Items

#### 1D Branch-Junction Orientation Physics
- [x] **BRC-001 [patch]**: Remove nonnegative parent-flow rejection from two-way and three-way branch solvers.
- [x] **BRC-002 [patch]**: Solve branch splits on `|Q_parent|` and reapply flow orientation to daughter flows.
- [x] **BRC-003 [patch]**: Compute wall-shear and apparent viscosity from flow magnitude while preserving signed pressure drops.
- [x] **BRC-004 [patch]**: Add value-semantic reverse-flow tests for pressure-balanced and prescribed two-way/three-way branches.
- [x] **BRC-005 [patch]**: Verify the touched `cfd-1d` branch module with bounded Cargo check, unit test, and nextest runs.

## Sprint 1.96.7: cfd-3d Venturi Pressure-Coefficient Physics
**Status**: Completed
**Start Date**: May 4, 2026

### Sprint Objectives
- Correct `cfd-3d` Venturi pressure coefficients so the solver follows its documented throat dynamic-pressure definition.
- Add value-semantic tests for coefficient scaling and undefined zero-flux rejection.

### Sprint Backlog Items

#### 3D Venturi Coefficient Physics
- [x] **V3D-001 [patch]**: Replace inlet dynamic-pressure scaling with throat dynamic pressure derived from face-integrated flow and throat area.
- [x] **V3D-002 [patch]**: Preserve typed rejection when coefficient scaling has no positive throat flow or area.
- [x] **V3D-003 [patch]**: Add direct regression tests for throat coefficient values and zero-flux rejection.
- [x] **V3D-004 [patch]**: Verify the touched `cfd-3d` Venturi module with bounded Cargo check, unit test, and nextest runs.

## Sprint 1.96.6: cfd-2d Explicit Stability Physics
**Status**: Completed
**Start Date**: May 4, 2026

### Sprint Objectives
- Correct `cfd-2d` explicit time-step physics so 2D advection and diffusion bounds follow summed-rate CFL and von Neumann stability conditions.
- Add value-semantic regression tests that verify `max_stable_dt` satisfies the same stability invariants reported by `advection_cfl` and `diffusion_number`.

### Sprint Backlog Items

#### 2D Stability Physics
- [x] **CFL-001 [patch]**: Replace componentwise minimum advection time-step bound with `1 / (|u|/dx + |v|/dy)`.
- [x] **CFL-002 [patch]**: Replace `0.5*min(dx²,dy²)/ν` diffusion bound with `0.5 / (ν(1/dx² + 1/dy²))`.
- [x] **CFL-003 [patch]**: Add regression tests for summed 2D advection and diffusion stability limits.
- [x] **CFL-004 [patch]**: Verify the touched `cfd-2d` stability module with bounded Cargo check, unit test, and nextest runs.

## Sprint 1.96.5: cfd-1d Venturi Reverse-Flow Physics
**Status**: Completed
**Start Date**: May 4, 2026

### Sprint Objectives
- Correct `cfd-1d` Venturi resistance physics so reverse-flow inputs use the same scalar loss magnitudes as forward flow.
- Preserve signed reported velocities while making Reynolds number, viscosity, friction factor, and pressure-loss coefficients depend on velocity magnitude.

### Sprint Backlog Items

#### 1D Venturi Resistance Physics
- [x] **VNT-001 [patch]**: Make `VenturiModel::calculate_coefficients` decompose losses with `|V|` and `|Q|` instead of treating negative inlet velocity as zero flow.
- [x] **VNT-002 [patch]**: Make detailed Venturi analysis use shear-rate and Reynolds magnitudes for reverse-flow cases.
- [x] **VNT-003 [patch]**: Add regression tests for symmetric-geometry coefficient and pressure-loss orientation invariance.
- [x] **VNT-004 [patch]**: Verify the touched `cfd-1d` venturi module with bounded Cargo check, unit test, and nextest runs.

## Sprint 1.96.4: Geometric Conservation Residual Verification
**Status**: Completed
**Start Date**: May 4, 2026

### Sprint Objectives
- Replace copy-through GCL evolution in `cfd-validation` with an executable conservative finite-volume residual.
- Verify Euler and Runge-Kutta constant-state preservation while proving non-constant fields are evaluated by the residual.

### Sprint Backlog Items

#### Conservation Validation Physics
- [x] **GCL-001 [patch]**: Implement a conservative second-order face-flux residual for the geometric conservation checker.
- [x] **GCL-002 [patch]**: Run Euler, midpoint, SSPRK3, and RK4 checks through real residual-based stage updates.
- [x] **GCL-003 [patch]**: Add value-semantic tests for unsupported RK stages and quadratic-field residual sensitivity.
- [x] **GCL-004 [patch]**: Verify the touched library target with bounded Cargo check, unit test, and nextest runs.

## Sprint 1.96.3: Womersley Analytical SSOT
**Status**: Completed
**Start Date**: May 3, 2026

### Sprint Objectives
- Replace the validation-local approximate Womersley profile with the canonical exact Bessel implementation from `cfd-1d`.
- Add value-semantic no-slip verification for the wall velocity invariant.

### Sprint Backlog Items

#### Analytical Validation SSOT
- [x] **WOM-001 [patch]**: Delegate `cfd-validation` Womersley velocity, wall shear stress, and flow rate to the `cfd-1d` exact `WomersleyProfile`.
- [x] **WOM-002 [patch]**: Add wall no-slip regression coverage for multiple phases.
- [x] **WOM-003 [patch]**: Verify the touched analytical module with bounded Cargo checks and tests.

## Sprint 1.96.2: Optimization Terminology Contract Cleanup
**Status**: Completed
**Start Date**: May 3, 2026

### Sprint Objectives
- Remove misleading unsupported/stub-like wording from executable boundary and optimization contracts without changing numerical behavior.
- Keep the `cfd-schematics` serpentine optimization path generator named by role rather than by reduced implementation detail.

### Sprint Backlog Items

#### Contract Terminology Cleanup
- [x] **TERM-001 [patch]**: Rename internal `generate_simplified_serpentine_path` to `generate_optimization_serpentine_path` and update all optimization call sites.
- [x] **TERM-002 [patch]**: Replace boundary stencil `"not implemented"` error text with explicit unsupported-order wording.
- [x] **TERM-003 [patch]**: Verify touched crates with bounded Cargo check and record nextest compile-bound status.

## Sprint 1.96.1: Workspace SSOT Cleanup
**Status**: Completed
**Start Date**: April 29, 2026

### Sprint Objectives
- Remove obsolete tracked root-source artifacts that duplicate canonical `gaia` CSG and mesh assembly code.
- Preserve the Cargo workspace authority: only files referenced by package manifests, module trees, examples, tests, docs, or report-generation assets remain as source artifacts.

### Sprint Backlog Items

#### Source Authority Cleanup
- [x] **SSOT-001 [patch]**: Delete unreferenced root historical Rust artifacts `old_assemble.rs`, `old_arrangement.rs`, `old_phase2.rs`, `old_operations.rs`, `old_indexed.rs`, `old_gwn_bvh.rs`, `old_seam.rs`, `old_phase4.rs`, and empty `csg_bi.rs`.
- [x] **SSOT-002 [patch]**: Verify no Cargo manifest, crate module, example, test, documentation, or report artifact references the deleted files.
- [x] **SSOT-003 [patch]**: Configure Windows GNU builds to use MSYS2 clang/lld and LLVM archive tools for Rust linking and C/C++ build scripts.
- [x] **SSOT-004 [patch]**: Remove misleading audit-trigger terminology from explicit unsupported-operation paths and already-bounded model comments.

## Sprint 1.96.0: Hydrodynamic Cavitation On-a-Chip (HCOC) Modeling
**Status**: Completed
**Start Date**: Next

### Sprint Objectives
- Enable mathematical modeling of distinct cellular injury zones (lysis, necrosis, permeabilization) from hydrodynamic cavitation in micro-orifices.
- Implement stiffness-coupled heterogeneous nucleation to predict altered cavitation inception times for Circulating Tumor Cells (CTCs).
- Provide unified mathematical specifications and zero-cost `GhostCell` topological integration for the new physics models.

### Sprint Backlog Items

#### Cellular Injury Modeling
- [x] **HCOC-001**: Implement `bio_damage.rs` with scalar fractional damage models representing cell membrane strain, mapping cavitation intensity to lysis/necrosis ratios.
- [x] **HCOC-002**: Derive and integrate property tests validating conservation of cellular mass across injury states.

#### CTC Nucleation Transport
- [x] **HCOC-003**: Implement `heterogeneous_nucleation.rs` replacing the linear $k_n$ scalar with a tensorial or struct-based multi-population interfacial tension approach.
- [x] **HCOC-004**: Validate modified effective vapor pressure inception offsets against analytical expected bounds.
- [x] **HCOC-005**: Bind Apollo-backed periodic DNS transforms to a validated reusable `FftPlan3D`.
- [x] **HCOC-006**: Correct `cfd-2d` MUSCL3/QUICK face reconstruction and document the polynomial reproduction theorem.
- [x] **HCOC-007**: Separate `cfd-1d` margination wall-induced and shear-gradient inertial lift scaling and document the reference-equilibrium proof.
- [x] **HCOC-008**: Replace `cfd-2d` turbulence benchmark placeholder dispatch with a closed supported-model registry.
- [x] **HCOC-009**: Restore coupled cfd-1d blood pressure-event solves with finite startup viscosity and row-equilibrated hydraulic linear systems.
- [x] **HCOC-010**: Replace compact `cfd-1d` plasma-skimming screening with the threshold-aware Pries phase-separation model and Murray-inferred sibling geometry.
- [x] **HCOC-011**: Replace `cfd-2d` serpentine mixing's exponential estimate with the Neumann transverse-diffusion eigenfunction solution.
- [x] **HCOC-012**: Replace `cfd-3d` LES turbulent-kinetic-energy viscosity aliases with a shared Yoshizawa SGS energy relation.
- [x] **HCOC-013**: Remove silent clamps from `cfd-2d` Pries plasma-skimming and add checked physical-envelope regression tests.
- [x] **HCOC-014**: Replace `cfd-2d` WALE boundary zero-gradient reductions with second-order one-sided gradient stencils.
- [x] **HCOC-015**: Replace `cfd-3d` Spalart-Allmaras all-zero TKE output with a wall-distance Yoshizawa diagnostic.
- [x] **HCOC-016**: Reject uninitialized `cfd-3d` k-epsilon state instead of synthesizing zero turbulence fields.
- [x] **HCOC-017**: Replace `cfd-2d` Smagorinsky LES zero SGS energy/dissipation placeholders with Yoshizawa diagnostics and second-order boundary strain.
- [x] **HCOC-018**: Remove the `cfd-1d` margination singular wall cutoff and silent lateral-position clamp while preserving inertial-lift scaling tests.
- [x] **HCOC-019**: Remove stored `cfd-1d` droplet occupied-channel state so finite-length occupancy spans are the only authoritative representation.
- [x] **HCOC-020**: Remove residual `cfd-2d` Smagorinsky validation SGS floors and keep LES validation on physical zero-floor defaults.
- [x] **HCOC-021**: Correct Milestone 12 cross-track therapy utility so ultrasound-only SDT receives acoustic-cavitation delivery credit from channel resonance, cancer-center enrichment, and active therapy-zone fraction.

## Sprint 1.95.1: CFD-MESH 3D Performance & Memory Optimization
**Status**: In Progress
**Start Date**: March 29, 2026

### Sprint Objectives
- Eliminate WASM memory exhaustion (OOM) during high-resolution 3D mesh generation.
- Refactor `BowyerWatson3D` kernel for strictly zero-allocation execution inside the point insertion loop.
- Flatten memory structures (O(N) contiguous buffers) for WASM/heap efficiency.
- Prove invariants mathematically via property tests, ensuring exact numerical determinism.

### Sprint Backlog Items

#### Core Memory Flattening
- [x] **MESH-PERF-001**: Lift `HashMap` and `Vec` allocations out of `BowyerWatson3D::insert_point`.
- [x] **MESH-PERF-002**: Pre-allocate `IndexedMesh` internal arrays in `build_volume` using AABB heuristic.
- [x] **MESH-PERF-003**: Optimize `SdfMesher` BCC lattice generation loops to prevent over-allocation.

## Sprint 1.91.0: Advanced Validation Framework Expansion
**Status**: Archived
**Start Date**: November 2, 2025

### Sprint Objectives
- Extend Method of Manufactured Solutions (MMS) test coverage for complex geometries
- Implement Richardson extrapolation automation for error analysis
- Add performance benchmarking framework for production scaling validation
- Enhance validation suite with comprehensive error metrics and convergence studies

### Sprint Backlog Items

#### MMS Test Coverage Expansion
- [x] **MMS-001**: Implement MMS for non-rectangular domains (circular, annular geometries) - Tasks 1.1, 1.2 & 1.3 COMPLETED
- [ ] **MMS-002**: Add MMS validation for complex boundary conditions (mixed BC types)
- [ ] **MMS-003**: Extend MMS to turbulent flow regimes with manufactured turbulence quantities
- [ ] **MMS-004**: Implement MMS for multi-physics coupling (heat transfer, species transport)

#### Richardson Extrapolation Automation
- [x] **RE-001**: Create Richardson extrapolation core library for grid convergence studies ✅ COMPLETED
- [x] **RE-002**: Implement automated grid refinement strategies (geometric progression) ✅ COMPLETED
- [x] **RE-003**: Add error estimation and convergence rate calculation ✅ COMPLETED
- [x] **RE-004**: Integrate Richardson extrapolation with existing MMS framework ✅ COMPLETED

#### Performance Benchmarking Framework
- [ ] **PB-001**: Design comprehensive benchmarking suite for CFD operations
- [ ] **PB-002**: Implement memory usage profiling and optimization tracking
- [ ] **PB-003**: Add parallel scaling analysis for multi-core/multi-GPU configurations
- [ ] **PB-004**: Create performance regression detection and alerting

#### Validation Infrastructure
- [ ] **VI-001**: Enhance validation test organization with clear taxonomy
- [ ] **VI-002**: Implement automated validation report generation (HTML/PDF)
- [ ] **VI-003**: Add validation metrics dashboard and visualization
- [ ] **VI-004**: Create validation configuration management system

### Completed Items
- [x] **PRE-001**: Resolve all build, test, and example errors (except known SIMPLEC convergence issues)
- [x] **PRE-002**: Clean up compiler warnings and code quality issues
- [x] **PRE-003**: Fix boundary condition handling in momentum solver

### Known Issues
- **SIMPLEC Convergence**: SIMPLEC/PIMPLE solver has matrix construction issues resolved, but pressure Poisson system remains singular due to incorrect reference pressure handling. Requires deeper investigation of boundary condition implementation and pressure-velocity coupling algorithm.
- **GPU Feature Gaps**: Some GPU-related code is unused due to feature flag configuration

### Sprint Metrics
- **Test Coverage**: Maintain >80% coverage across all crates
- **Performance**: No regression in benchmark performance
- **Code Quality**: Zero warnings, comprehensive documentation
- **Validation**: All new MMS tests pass with expected convergence rates
